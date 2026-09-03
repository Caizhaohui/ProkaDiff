//! CIGAR pileup into in-memory columns (no FASTQ). BAM records map onto `AlignedRead`.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::fasta::FastaRecord;
use crate::jc::SubAlignment;
use crate::jc_seq::{MIN_CLIP_FOR_PLACE, MIN_CLIP_FOR_SEED};
use crate::ra::{BaseObs, PileupColumn, Strand};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarKind {
    Match,
    Ins,
    Del,
    SoftClip,
    HardClip,
    Skip,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CigarOp {
    pub kind: CigarKind,
    pub len: usize,
}

#[derive(Clone, Debug)]
pub struct AlignedRead {
    pub contig_idx: usize,
    /// 0-based leftmost reference position.
    pub ref_start_0: i64,
    pub minus: bool,
    pub seq: Vec<u8>,
    pub cigar: Vec<CigarOp>,
    pub mapq: u8,
}

pub const UNIQUE_MAPQ: u8 = 10;
pub const END_TRIM: usize = 1;
/// Softclips shorter than this are not recorded as JC seed sides.
/// Second-pass accept still requires ≥14 bp into each side (`jc::accept_junction`).
pub const MIN_CLIP_FOR_JC: usize = 14;

#[derive(Clone, Debug)]
pub struct SplitCandidate {
    pub contig_idx: usize,
    pub side2_contig_idx: usize,
    pub minus: bool,
    pub side1_minus: bool,
    pub side2_minus: bool,
    pub left: SubAlignment,
    pub right: SubAlignment,
    pub side1_pos_1: u64,
    pub side2_pos_1: u64,
    pub overlap: i64,
    pub origin: SplitOrigin,
}

/// Where a split candidate came from. Used by the engine's redundancy
/// folding to tell "one read, several placement hits" apart from genuinely
/// independent support.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SplitOrigin {
    /// Long CIGAR I/D: the candidate is its own identity and must never be
    /// folded away by the multi-copy placement rule.
    Cigar,
    /// Softclip placement: index of the originating `SoftclipHint`. A read
    /// whose clip has several exact placement hits (repeat copies)
    /// contributes one candidate per hit, all sharing this identity.
    Softclip(usize),
    /// Candidate-junction second pass: hash of the read name. Kept in its
    /// own variant so a qname hash can never collide with a first-pass
    /// `SoftclipHint` index in the supporting-read identity set.
    SecondPass(u64),
}

#[derive(Clone, Debug)]
pub struct SoftclipHint {
    pub contig_idx: usize,
    /// Read strand (drives plus/minus read counts in `JunctionSupport`).
    pub minus: bool,
    /// 1-based reference position of the aligned base next to the clip.
    /// SAM/BAM stores SEQ reference-forward (minus-strand reads are
    /// reverse-complemented), so this is strand-independent: last match
    /// for a trailing clip, first match for a leading clip.
    pub aligned_pos_1: u64,
    /// True if the clip lies at LOWER reference coordinates than the
    /// aligned block (leading clip of the stored sequence). Strand-
    /// independent for the same reason.
    pub clip_is_left: bool,
    /// Clipped bases in stored (reference-forward) orientation;
    /// `place_softclip` searches both strands.
    pub clip_seq: Vec<u8>,
    pub aligned_q_len: usize,
}

/// Increment per-base depth for every CIGAR match (any MAPQ). Used to tell
/// unique-depth-0 repeats (multi-mapping still covers) from true deletions.
impl SplitCandidate {
    fn local(
        contig_idx: usize,
        minus: bool,
        left: SubAlignment,
        right: SubAlignment,
        side1_pos_1: u64,
        side2_pos_1: u64,
        overlap: i64,
    ) -> Self {
        Self {
            contig_idx,
            side2_contig_idx: contig_idx,
            minus,
            // breseq side-strand convention: each side's strand is the
            // direction its sequence extends away from the junction.
            // side1 (left flank end) extends leftward, side2 (right flank
            // start) extends rightward — matches the measured oracle
            // deletion-flank JC (e.g. Clonal 3289961 -1 -> 3289978 1).
            side1_minus: true,
            side2_minus: false,
            left,
            right,
            side1_pos_1,
            side2_pos_1,
            overlap,
            origin: SplitOrigin::Cigar,
        }
    }
}

/// Increment per-base depth for every CIGAR match (any MAPQ). Used to tell
/// unique-depth-0 repeats (multi-mapping still covers) from true deletions.
fn add_match_depth(read: &AlignedRead, depth: &mut [u32]) {
    if read.ref_start_0 < 0 {
        return;
    }
    let mut ref_pos = read.ref_start_0 as usize;
    for op in &read.cigar {
        match op.kind {
            CigarKind::Match => {
                for _ in 0..op.len {
                    if ref_pos < depth.len() {
                        depth[ref_pos] = depth[ref_pos].saturating_add(1);
                    }
                    ref_pos += 1;
                }
            }
            CigarKind::Del | CigarKind::Skip => {
                ref_pos += op.len;
            }
            CigarKind::Ins | CigarKind::SoftClip | CigarKind::HardClip => {}
        }
    }
}

pub fn apply_read(
    read: &AlignedRead,
    columns: &mut [PileupColumn],
    unique_depth: &mut [u32],
    total_depth: &mut [u32],
    splits: &mut Vec<SplitCandidate>,
    clips: &mut Vec<SoftclipHint>,
) {
    if read.ref_start_0 < 0 {
        return;
    }
    add_match_depth(read, total_depth);
    if read.mapq < UNIQUE_MAPQ {
        return;
    }
    let strand = if read.minus {
        Strand::Minus
    } else {
        Strand::Plus
    };
    let read_len = read.seq.len();
    let mut ref_pos = read.ref_start_0 as usize;
    let mut qpos = 0usize;
    let mut hits: Vec<(usize, usize, u8)> = Vec::new(); // ref, q, base-or-dash
    let mut ins_at: Vec<(usize, Vec<u8>)> = Vec::new(); // after this ref index (0-based)

    let mut last_match_ref: Option<usize> = None;
    let mut last_match_qend = 0usize;
    let mut first_match_ref: Option<usize> = None;
    let mut pending_left_clip: Option<Vec<u8>> = None;

    for op in &read.cigar {
        match op.kind {
            CigarKind::Match => {
                if first_match_ref.is_none() {
                    first_match_ref = Some(ref_pos);
                }
                for _ in 0..op.len {
                    if ref_pos < columns.len() && qpos < read.seq.len() {
                        hits.push((ref_pos, qpos, read.seq[qpos].to_ascii_uppercase()));
                    }
                    last_match_ref = Some(ref_pos);
                    last_match_qend = qpos + 1;
                    ref_pos += 1;
                    qpos += 1;
                }
            }
            CigarKind::Ins => {
                if op.len > 2 {
                    if let Some(left_ref) = last_match_ref {
                        let left = SubAlignment {
                            read_start: 0,
                            read_end: last_match_qend,
                        };
                        let right = SubAlignment {
                            read_start: qpos,
                            read_end: read_len,
                        };
                        splits.push(SplitCandidate::local(
                            read.contig_idx,
                            read.minus,
                            left,
                            right,
                            left_ref as u64 + 1,
                            left_ref as u64 + 1,
                            0,
                        ));
                    }
                } else if op.len >= 1 && ref_pos > 0 {
                    let mut oligo = Vec::new();
                    for _ in 0..op.len {
                        if qpos < read.seq.len() {
                            oligo.push(read.seq[qpos].to_ascii_uppercase());
                        }
                        qpos += 1;
                    }
                    ins_at.push((ref_pos - 1, oligo));
                    continue;
                }
                qpos += op.len;
            }
            CigarKind::Del => {
                if op.len > 2 {
                    if let Some(left_ref) = last_match_ref {
                        let left = SubAlignment {
                            read_start: 0,
                            read_end: last_match_qend,
                        };
                        let right_q = qpos;
                        let right = SubAlignment {
                            read_start: right_q,
                            read_end: read_len,
                        };
                        splits.push(SplitCandidate::local(
                            read.contig_idx,
                            read.minus,
                            left,
                            right,
                            left_ref as u64 + 1,
                            (ref_pos + op.len) as u64 + 1,
                            -(op.len as i64),
                        ));
                    }
                } else {
                    for _ in 0..op.len {
                        if ref_pos < columns.len() {
                            hits.push((ref_pos, qpos, b'-'));
                        }
                        ref_pos += 1;
                    }
                    continue;
                }
                ref_pos += op.len;
            }
            CigarKind::SoftClip => {
                if op.len >= MIN_CLIP_FOR_SEED && qpos + op.len <= read.seq.len() {
                    let clip = read.seq[qpos..qpos + op.len].to_ascii_uppercase();
                    if last_match_ref.is_none() {
                        pending_left_clip = Some(clip);
                    } else if let Some(pos) = last_match_ref {
                        // Trailing clip of the stored sequence. SAM/BAM stores SEQ
                        // reference-forward (minus-strand reads are reverse-
                        // complemented), so a trailing clip lies at HIGHER
                        // reference coordinates on BOTH strands and the
                        // junction-adjacent aligned base is the last match.
                        // (Measured, synth_is_mob 2371412: minus M+S reads ending
                        // at 602 carry motif-prefix clips; keying them at the
                        // first match scattered minus candidates over 517-580 and
                        // stranded every junction single-strand.)
                        clips.push(SoftclipHint {
                            contig_idx: read.contig_idx,
                            minus: read.minus,
                            aligned_pos_1: pos as u64 + 1,
                            clip_is_left: false,
                            clip_seq: clip,
                            aligned_q_len: last_match_qend,
                        });
                    }
                }
                qpos += op.len;
            }
            CigarKind::HardClip => {}
            CigarKind::Skip => {
                ref_pos += op.len;
            }
        }
    }

    if let (Some(clip), Some(first)) = (pending_left_clip, first_match_ref) {
        let aligned_q_len = last_match_qend.saturating_sub(clip.len());
        if aligned_q_len >= MIN_CLIP_FOR_SEED {
            // Leading clip of the stored (reference-forward) sequence: lies
            // at LOWER reference coordinates on both strands; the
            // junction-adjacent aligned base is the first match.
            clips.push(SoftclipHint {
                contig_idx: read.contig_idx,
                minus: read.minus,
                aligned_pos_1: first as u64 + 1,
                clip_is_left: true,
                clip_seq: clip,
                aligned_q_len,
            });
        }
    }

    if hits.len() > END_TRIM * 2 {
        hits = hits[END_TRIM..hits.len() - END_TRIM].to_vec();
    } else {
        hits.clear();
    }

    for (r, _, base) in &hits {
        if *r < columns.len() {
            columns[*r].observations.push(BaseObs {
                base: *base,
                strand,
            });
            unique_depth[*r] = unique_depth[*r].saturating_add(1);
        }
    }
    for (r, oligo) in ins_at {
        if r < columns.len() {
            columns[r].insertions.push((oligo, strand));
        }
    }
}

fn rc_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            x => x,
        })
        .collect()
}

fn find_exact(hay: &[u8], needle: &[u8]) -> Vec<usize> {
    if needle.is_empty() || needle.len() > hay.len() {
        return Vec::new();
    }
    hay.windows(needle.len())
        .enumerate()
        .filter(|(_, w)| *w == needle)
        .map(|(i, _)| i)
        .collect()
}

/// Seed length of the placement index. Every searched clip is at least
/// [`MIN_CLIP_FOR_PLACE`] long, so its first `PLACE_K` bases always exist and
/// no needle needs a shorter seed. Capped at 32 by the 2-bit-per-base `u64`
/// packing. Tying it to the clip floor means lowering that floor
/// automatically makes the seed shorter (less selective, still correct)
/// rather than silently breaking the lookup.
const PLACE_K: usize = if MIN_CLIP_FOR_PLACE < 32 {
    MIN_CLIP_FOR_PLACE
} else {
    32
};

fn base_2bit(b: u8) -> Option<u64> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Pack the first `PLACE_K` bases into 2 bits each. `None` if the seed is
/// short or holds a non-ACGT base (ambiguity codes are not representable).
fn pack_seed(seq: &[u8]) -> Option<u64> {
    if seq.len() < PLACE_K {
        return None;
    }
    let mut key = 0u64;
    for &b in &seq[..PLACE_K] {
        key = (key << 2) | base_2bit(b)?;
    }
    Some(key)
}

struct PlaceIndex {
    hays: Vec<Vec<u8>>,
    /// `(packed first-PLACE_K-mer, contig, 0-based start)`, sorted by key.
    ///
    /// Replaces a genome-wide `hay.windows(needle.len())` scan per clip per
    /// strand, which cost O(genome x clips) byte comparisons — the measured
    /// cause of Clonal 2372539 running 19 h on 43k short clips before being
    /// cancelled (`benchmark/results/clonal_leftover.md`). A seed lookup is
    /// a binary search plus verification of the few positions sharing the
    /// seed.
    ///
    /// Flat and sorted rather than `HashMap<u64, Vec<_>>`: there is one entry
    /// per reference position, so on a 4.6 Mbp chromosome this is ~74 MB
    /// against several hundred for per-key heap allocations (peak RSS is
    /// already ~2.9 GB).
    ///
    /// Forward strand only. A reverse-complement query is answered by looking
    /// up the RC of the needle, which lands on exactly the coordinates a
    /// forward scan for the RC sequence would report.
    seeds: Vec<(u64, u32, u32)>,
}

pub const MAX_CLIP_HITS: usize = 20;

impl PlaceIndex {
    fn from_fasta(fasta: &[FastaRecord]) -> Self {
        let hays: Vec<Vec<u8>> = fasta
            .iter()
            .map(|r| {
                r.seq
                    .iter()
                    .copied()
                    .map(|b| b.to_ascii_uppercase())
                    .collect()
            })
            .collect();
        let total: usize = hays.iter().map(|h| h.len()).sum();
        let mut seeds: Vec<(u64, u32, u32)> = Vec::with_capacity(total);
        let mask = if PLACE_K >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * PLACE_K)) - 1
        };
        for (ci, hay) in hays.iter().enumerate() {
            let mut key = 0u64;
            // Bases since the last non-ACGT: a seed is only emitted once
            // PLACE_K consecutive representable bases are in the window.
            let mut run = 0usize;
            for (i, &b) in hay.iter().enumerate() {
                match base_2bit(b) {
                    Some(v) => {
                        key = ((key << 2) | v) & mask;
                        run += 1;
                    }
                    None => {
                        run = 0;
                        key = 0;
                        continue;
                    }
                }
                if run >= PLACE_K {
                    let start = i + 1 - PLACE_K;
                    seeds.push((key, ci as u32, start as u32));
                }
            }
        }
        seeds.sort_unstable();
        Self { hays, seeds }
    }

    /// Exact hits of `needle_upper` on either strand, as
    /// `(contig, 0-based start, matched_reverse_complement)`.
    fn search(&self, needle_upper: &[u8]) -> Vec<(usize, usize, bool)> {
        if needle_upper.len() < MIN_CLIP_FOR_PLACE {
            return Vec::new();
        }
        let rc = rc_dna(needle_upper);
        let mut hits = Vec::new();
        self.collect(needle_upper, false, &mut hits);
        self.collect(&rc, true, &mut hits);
        if hits.len() > MAX_CLIP_HITS {
            // Highly repetitive across genome (e.g. simple sequence / high copy);
            // skip to avoid combinatorial explosion.
            return Vec::new();
        }
        // Per contig: forward hits ascending, then RC hits ascending.
        hits.sort_unstable_by_key(|&(ci, pos, is_rc)| (ci, is_rc, pos));
        hits
    }

    fn collect(&self, needle: &[u8], is_rc: bool, out: &mut Vec<(usize, usize, bool)>) {
        let Some(key) = pack_seed(needle) else {
            // Seed holds a non-ACGT base, so it is absent from the index.
            // Rare; scan so results stay identical to a full search.
            for (ci, hay) in self.hays.iter().enumerate() {
                for h in find_exact(hay, needle) {
                    out.push((ci, h, is_rc));
                }
            }
            return;
        };
        let lo = self.seeds.partition_point(|e| e.0 < key);
        for &(_, ci, pos) in self.seeds[lo..].iter().take_while(|e| e.0 == key) {
            let (ci, pos) = (ci as usize, pos as usize);
            let hay = &self.hays[ci];
            if hay.len() - pos >= needle.len() && &hay[pos..pos + needle.len()] == needle {
                out.push((ci, pos, is_rc));
            }
        }
    }
}

fn clip_key(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .copied()
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

fn is_continuation(hint: &SoftclipHint, hit_contig: usize, hit_0: usize, clip_len: usize) -> bool {
    if hit_contig != hint.contig_idx {
        return false;
    }
    if hint.clip_is_left {
        let hit_end_1 = hit_0 as u64 + clip_len as u64;
        hit_end_1 == hint.aligned_pos_1 || hit_end_1 + 1 == hint.aligned_pos_1
    } else {
        let hit_start_1 = hit_0 as u64 + 1;
        hit_start_1 == hint.aligned_pos_1 + 1 || hit_start_1 == hint.aligned_pos_1
    }
}

fn hits_to_splits(
    hint: &SoftclipHint,
    hint_id: usize,
    clip_len: usize,
    hits: &[(usize, usize, bool)],
) -> Vec<SplitCandidate> {
    let mut out = Vec::new();
    for &(idx, hit_0, is_rc) in hits {
        if is_continuation(hint, idx, hit_0, clip_len) {
            continue;
        }
        // content_rc_hit: the stored clip is reference-forward (SAM/BAM
        // reverse-complements minus-strand reads), so the clip content
        // reads along the REVERSE strand of the hit region iff the
        // reverse-complement search matched — independent of read strand.
        // The junction-adjacent base is the hit start iff the clip side
        // and the content orientation agree, otherwise the hit end.
        let content_rc_hit = is_rc;
        let side2_pos_1 = if hint.clip_is_left == content_rc_hit {
            hit_0 as u64 + 1
        } else {
            hit_0 as u64 + clip_len as u64
        };
        // Emitted side strands follow the breseq convention: a side's
        // strand is the direction its sequence extends AWAY from the
        // junction in reference coordinates. side1's aligned block
        // extends opposite the clip side; side2 extends toward lower
        // coordinates iff the junction-adjacent base is the hit end.
        // (Reproduces the measured oracle strands of both synth_is_mob
        // 2371412 junctions; inverted-hit placements are unmeasured.)
        out.push(SplitCandidate {
            contig_idx: hint.contig_idx,
            side2_contig_idx: idx,
            minus: hint.minus,
            side1_minus: !hint.clip_is_left,
            side2_minus: hint.clip_is_left ^ content_rc_hit,
            left: SubAlignment {
                read_start: 0,
                read_end: hint.aligned_q_len,
            },
            right: SubAlignment {
                read_start: 0,
                read_end: clip_len,
            },
            side1_pos_1: hint.aligned_pos_1,
            side2_pos_1,
            overlap: 0,
            origin: SplitOrigin::Softclip(hint_id),
        });
    }
    out
}

/// Place softclips onto the reference (exact match, both strands).
///
/// Clips shorter than [`MIN_CLIP_FOR_PLACE`] are not searched (recorded
/// seeds still exist for unique-side evidence). Identical clip sequences
/// are searched once against an uppercase-once index.
pub fn place_softclips(clips: &[SoftclipHint], fasta: &[FastaRecord]) -> Vec<SplitCandidate> {
    let index = PlaceIndex::from_fasta(fasta);
    let mut keys: Vec<Vec<u8>> = Vec::new();
    let mut seen: HashMap<Vec<u8>, ()> = HashMap::new();
    for hint in clips {
        if hint.clip_seq.len() < MIN_CLIP_FOR_PLACE || hint.aligned_q_len < MIN_CLIP_FOR_SEED {
            continue;
        }
        let key = clip_key(&hint.clip_seq);
        if seen.insert(key.clone(), ()).is_none() {
            keys.push(key);
        }
    }
    keys.sort();
    let searched: HashMap<Vec<u8>, Vec<(usize, usize, bool)>> = keys
        .into_par_iter()
        .map(|k| {
            let hits = index.search(&k);
            (k, hits)
        })
        .collect();
    let mut out = Vec::new();
    for (hint_id, hint) in clips.iter().enumerate() {
        if hint.clip_seq.len() < MIN_CLIP_FOR_PLACE || hint.aligned_q_len < MIN_CLIP_FOR_SEED {
            continue;
        }
        let key = clip_key(&hint.clip_seq);
        if let Some(hits) = searched.get(&key) {
            out.extend(hits_to_splits(hint, hint_id, key.len(), hits));
        }
    }
    out
}

/// Place a single softclip. Repeat copies are allowed. Skip the trivial
/// continuation of the same alignment. `hint_id` tags every produced
/// candidate so the engine can fold multi-copy placements.
pub fn place_softclip(
    hint: &SoftclipHint,
    fasta: &[FastaRecord],
    hint_id: usize,
) -> Vec<SplitCandidate> {
    let mut out = place_softclips(std::slice::from_ref(hint), fasta);
    for c in &mut out {
        c.origin = SplitOrigin::Softclip(hint_id);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ra::{call_consensus, RaOptions};

    fn cols(seq: &[u8]) -> Vec<PileupColumn> {
        seq.iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: Vec::new(),
                insertions: Vec::new(),
            })
            .collect()
    }

    fn match_read(seq: &[u8], start: i64, minus: bool) -> AlignedRead {
        AlignedRead {
            contig_idx: 0,
            ref_start_0: start,
            minus,
            seq: seq.to_vec(),
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: seq.len(),
            }],
            mapq: 40,
        }
    }

    #[test]
    fn pileup_snp_from_several_reads() {
        let ref_seq = b"ACGTACGTAC";
        let mut columns = cols(ref_seq);
        let mut depth = vec![0u32; ref_seq.len()];
        let mut splits = Vec::new();
        let mut clips = Vec::new();
        let mut alt = ref_seq.to_vec();
        alt[4] = b'T';
        for minus in [false, false, false, true, true, true] {
            apply_read(
                &match_read(&alt, 0, minus),
                &mut columns,
                &mut depth,
                &mut vec![0u32; ref_seq.len()],
                &mut splits,
                &mut clips,
            );
        }
        assert_eq!(
            call_consensus(&columns[4], &RaOptions::default()),
            crate::ra::ConsensusCall::Snp { alt: b'T' }
        );
        assert!(splits.is_empty());
    }

    #[test]
    fn long_deletion_emits_split_candidate() {
        let ref_seq = vec![b'A'; 40];
        let mut columns = cols(&ref_seq);
        let mut depth = vec![0u32; 40];
        let mut splits = Vec::new();
        let read = AlignedRead {
            contig_idx: 0,
            ref_start_0: 0,
            minus: false,
            seq: vec![b'A'; 30],
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: 15,
                },
                CigarOp {
                    kind: CigarKind::Del,
                    len: 8,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 15,
                },
            ],
            mapq: 40,
        };
        apply_read(
            &read,
            &mut columns,
            &mut depth,
            &mut [0u32; 40],
            &mut splits,
            &mut Vec::new(),
        );
        assert_eq!(splits.len(), 1);
        assert!(splits[0].overlap < 0);
    }

    #[test]
    fn three_prime_softclip_places_onto_repeat_copy() {
        let mut ref_seq = vec![b'C'; 160];
        let motif = *b"ACGTACGTACGTACGTACGT";
        ref_seq[40..60].copy_from_slice(&motif);
        ref_seq[100..120].copy_from_slice(&motif);
        let fasta = [FastaRecord {
            name: "chr".into(),
            seq: ref_seq.clone(),
        }];
        let mut clip_seq = vec![b'C'; 20];
        clip_seq.extend_from_slice(&motif);
        let read = AlignedRead {
            contig_idx: 0,
            ref_start_0: 70,
            minus: false,
            seq: clip_seq,
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: 20,
                },
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: 20,
                },
            ],
            mapq: 40,
        };
        let mut columns = cols(&ref_seq);
        let mut depth = vec![0u32; 160];
        let mut splits = Vec::new();
        let mut clips = Vec::new();
        apply_read(
            &read,
            &mut columns,
            &mut depth,
            &mut [0u32; 160],
            &mut splits,
            &mut clips,
        );
        assert_eq!(clips.len(), 1);
        assert_eq!(clips[0].aligned_pos_1, 90);
        assert!(!clips[0].clip_is_left);
        let placed = place_softclip(&clips[0], &fasta, 0);
        assert!(
            placed
                .iter()
                .any(|s| (41..=60).contains(&s.side2_pos_1) || (101..=120).contains(&s.side2_pos_1)),
            "expected a hit on an IS copy, got {placed:?}"
        );
    }

    #[test]
    fn minus_strand_leading_clip_keys_at_first_match() {
        // SAM/BAM stores SEQ reference-forward (minus-strand reads are
        // reverse-complemented), so a leading (5′-of-SEQ) softclip lies at
        // LOWER reference coordinates on both strands and the junction
        // base is the FIRST aligned reference base — the same rule as for
        // plus reads. Models a minus read crossing an insert RIGHT
        // boundary: stored SEQ = [motif suffix][right flank].
        // (Measured: synth_is_mob 2371412 minus S+M reads start at 599.)
        let ref_seq = vec![b'C'; 160];
        let mut seq = b"TGCAAGTCGATCGTTAGCCA".to_vec();
        seq.extend(std::iter::repeat_n(b'C', 20));
        let read = AlignedRead {
            contig_idx: 0,
            ref_start_0: 10, // match spans ref 11..30 (1-based)
            minus: true,
            seq,
            cigar: vec![
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: 20,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 20,
                },
            ],
            mapq: 40,
        };
        let mut columns = cols(&ref_seq);
        let mut clips = Vec::new();
        apply_read(
            &read,
            &mut columns,
            &mut vec![0u32; 160],
            &mut [0u32; 160],
            &mut Vec::new(),
            &mut clips,
        );
        assert_eq!(clips.len(), 1);
        assert_eq!(
            clips[0].aligned_pos_1, 11,
            "leading clip keys at first_match_ref+1 on both strands"
        );
        assert!(
            clips[0].clip_is_left,
            "leading clip lies at lower reference coordinates"
        );
        assert_eq!(
            clips[0].clip_seq, b"TGCAAGTCGATCGTTAGCCA",
            "stored SEQ is reference-forward, not reverse-complemented"
        );
    }

    #[test]
    fn minus_strand_trailing_clip_keys_at_last_match() {
        // Minus read, trailing (3′-of-SEQ) clip (CIGAR M+S): the stored SEQ
        // is reference-forward, so the clip lies at HIGHER reference
        // coordinates and the junction base is the LAST aligned reference
        // base — the same rule as for plus reads. Models a minus read
        // crossing an insert LEFT boundary: stored SEQ =
        // [left flank][motif prefix]. (Measured: synth_is_mob 2371412 minus
        // M+S reads end at 602 carrying motif-prefix clips; keying them at
        // the first match was the zero-JC root cause.)
        let ref_seq = vec![b'C'; 160];
        let mut seq = vec![b'C'; 20];
        seq.extend_from_slice(b"TGCAAGTCGATCGTTAGCCA");
        let read = AlignedRead {
            contig_idx: 0,
            ref_start_0: 10, // match spans ref 11..30 (1-based)
            minus: true,
            seq,
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: 20,
                },
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: 20,
                },
            ],
            mapq: 40,
        };
        let mut columns = cols(&ref_seq);
        let mut clips = Vec::new();
        apply_read(
            &read,
            &mut columns,
            &mut vec![0u32; 160],
            &mut [0u32; 160],
            &mut Vec::new(),
            &mut clips,
        );
        assert_eq!(clips.len(), 1);
        assert_eq!(
            clips[0].aligned_pos_1, 30,
            "trailing clip keys at last_match_ref+1 on both strands"
        );
        assert!(
            !clips[0].clip_is_left,
            "trailing clip lies at higher reference coordinates"
        );
    }

    fn motif_clip_read(motif: &[u8], match_len: usize) -> (Vec<u8>, AlignedRead) {
        let mut ref_seq = vec![b'C'; 160];
        let m = motif.len();
        ref_seq[40..40 + m].copy_from_slice(motif);
        ref_seq[100..100 + m].copy_from_slice(motif);
        let mut seq = vec![b'C'; match_len];
        seq.extend_from_slice(motif);
        let read = AlignedRead {
            contig_idx: 0,
            ref_start_0: 70,
            minus: false,
            seq,
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: match_len,
                },
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: m,
                },
            ],
            mapq: 40,
        };
        (ref_seq, read)
    }

    fn recorded_clip(ref_seq: &[u8], read: &AlignedRead) -> SoftclipHint {
        let mut columns = cols(ref_seq);
        let mut clips = Vec::new();
        apply_read(
            read,
            &mut columns,
            &mut vec![0u32; 160],
            &mut [0u32; 160],
            &mut Vec::new(),
            &mut clips,
        );
        assert_eq!(clips.len(), 1);
        clips.pop().unwrap()
    }

    #[test]
    fn five_bp_clip_is_below_the_seed_floor_and_not_recorded() {
        // MIN_CLIP_FOR_SEED == MIN_CLIP_FOR_PLACE == 6 (the published floor):
        // below it, a clip is not even recorded as a seed hint.
        let motif = *b"ACGTA";
        let (ref_seq, read) = motif_clip_read(&motif, 20);
        let mut columns = cols(&ref_seq);
        let mut clips = Vec::new();
        apply_read(
            &read,
            &mut columns,
            &mut vec![0u32; ref_seq.len()],
            &mut vec![0u32; ref_seq.len()],
            &mut Vec::new(),
            &mut clips,
        );
        assert!(
            clips.is_empty(),
            "S=5 must not be a seed hint; got {clips:?}"
        );
    }

    #[test]
    fn six_bp_clip_is_recorded_but_not_placed_genome_wide() {
        // MIN_CLIP_FOR_SEED is 6 (recorded as seed hint for unique side),
        // but MIN_CLIP_FOR_PLACE is 12 (documented in docs/parity.md):
        // 6 bp clips must NOT be placed genome-wide due to 4^6 = 4096 non-uniqueness.
        let motif = *b"ACGTAC";
        let (ref_seq, read) = motif_clip_read(&motif, 20);
        let fasta = [FastaRecord {
            name: "chr".into(),
            seq: ref_seq.clone(),
        }];
        let hint = recorded_clip(&ref_seq, &read);
        assert_eq!(hint.clip_seq, motif);
        let placed = place_softclip(&hint, &fasta, 0);
        assert!(
            placed.is_empty(),
            "6 bp clip is below the placement floor of 12 and must not place; got {placed:?}"
        );
    }

    #[test]
    fn ten_bp_clip_is_recorded_but_not_placed_genome_wide() {
        let motif = *b"ACGTACGTAC";
        let (ref_seq, read) = motif_clip_read(&motif, 20);
        let fasta = [FastaRecord {
            name: "chr".into(),
            seq: ref_seq.clone(),
        }];
        let hint = recorded_clip(&ref_seq, &read);
        assert_eq!(hint.clip_seq.len(), 10);
        let placed = place_softclip(&hint, &fasta, 0);
        assert!(
            placed.is_empty(),
            "10 bp clip is below the placement floor of 12 and must not place; got {placed:?}"
        );
    }

    #[test]
    fn highly_repetitive_clip_is_skipped() {
        // A 12 bp clip that matches > MAX_CLIP_HITS (20) times across the reference
        // must be skipped to protect against combinatorial explosion.
        let motif = *b"ACGTACGTACGT";
        let mut ref_seq = vec![b'A'; 500];
        for i in 0..25 {
            ref_seq[i * 15..i * 15 + 12].copy_from_slice(&motif);
        }
        let fasta = [FastaRecord {
            name: "chr".into(),
            seq: ref_seq.clone(),
        }];
        let hint = SoftclipHint {
            contig_idx: 0,
            minus: false,
            aligned_pos_1: 450,
            clip_is_left: false,
            clip_seq: motif.to_vec(),
            aligned_q_len: 12,
        };
        let placed = place_softclip(&hint, &fasta, 0);
        assert!(
            placed.is_empty(),
            "clip matching >20 times must be skipped, got {} hits",
            placed.len()
        );
    }

    #[test]
    fn twelve_bp_softclip_places_onto_repeat_copy() {
        // Clonal 36 bp: maxS=12, S≥14=0. Placement starts at 12 bp.
        let motif = *b"ACGTACGTACGT";
        let (ref_seq, read) = motif_clip_read(&motif, 20);
        let fasta = [FastaRecord {
            name: "chr".into(),
            seq: ref_seq.clone(),
        }];
        let hint = recorded_clip(&ref_seq, &read);
        let placed = place_softclip(&hint, &fasta, 0);
        assert!(
            placed
                .iter()
                .any(|s| (41..=52).contains(&s.side2_pos_1) || (101..=112).contains(&s.side2_pos_1)),
            "12 bp clip must place onto a motif copy, got {placed:?}"
        );
    }
}
