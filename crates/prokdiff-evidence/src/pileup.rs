//! CIGAR pileup into in-memory columns (no FASTQ). BAM records map onto `AlignedRead`.

use crate::fasta::FastaRecord;
use crate::jc::SubAlignment;
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
/// Softclips shorter than the published JC overlap floor are not mosaic sides.
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
}

#[derive(Clone, Debug)]
pub struct SoftclipHint {
    pub contig_idx: usize,
    pub minus: bool,
    /// 1-based reference position of the aligned base next to the clip
    /// (reference-oriented: for a minus-strand read the junction-adjacent
    /// base is on the opposite read end).
    pub aligned_pos_1: u64,
    /// Reference-oriented: true if the clip lies at LOWER reference
    /// coordinates than the aligned block (`clip_on_read_left XOR minus`).
    pub clip_is_left: bool,
    /// Clipped bases in READ orientation (`place_softclip` searches both strands).
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
            side1_minus: minus,
            side2_minus: minus,
            left,
            right,
            side1_pos_1,
            side2_pos_1,
            overlap,
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
                if op.len >= MIN_CLIP_FOR_JC && qpos + op.len <= read.seq.len() {
                    let clip = read.seq[qpos..qpos + op.len].to_ascii_uppercase();
                    if last_match_ref.is_none() {
                        pending_left_clip = Some(clip);
                    } else if let Some(pos) = last_match_ref {
                        // 3′ (read-right) clip. On a minus-strand read it lies at
                        // LOWER reference coordinates; the junction-adjacent
                        // aligned base is then the first match, not the last.
                        let ref_clip_is_left = read.minus;
                        let junction_ref = if ref_clip_is_left {
                            first_match_ref.expect("last match implies a first match")
                        } else {
                            pos
                        };
                        clips.push(SoftclipHint {
                            contig_idx: read.contig_idx,
                            minus: read.minus,
                            aligned_pos_1: junction_ref as u64 + 1,
                            clip_is_left: ref_clip_is_left,
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
        if aligned_q_len >= MIN_CLIP_FOR_JC {
            // 5′ (read-left) clip. On a minus-strand read it lies at HIGHER
            // reference coordinates; the junction-adjacent aligned base is
            // then the last match, not the first.
            let ref_clip_is_left = !read.minus;
            let junction_ref = if ref_clip_is_left {
                first
            } else {
                last_match_ref.expect("left clip implies a match")
            };
            clips.push(SoftclipHint {
                contig_idx: read.contig_idx,
                minus: read.minus,
                aligned_pos_1: junction_ref as u64 + 1,
                clip_is_left: ref_clip_is_left,
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
        .filter(|(_, w)| {
            w.iter()
                .zip(needle.iter())
                .all(|(a, b)| a.eq_ignore_ascii_case(b))
        })
        .map(|(i, _)| i)
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

/// Place a softclip onto the reference (exact match, both strands). Repeat
/// copies are allowed (unassigned JC). Skip the trivial continuation of the
/// same alignment.
pub fn place_softclip(hint: &SoftclipHint, fasta: &[FastaRecord]) -> Vec<SplitCandidate> {
    if hint.clip_seq.len() < MIN_CLIP_FOR_JC || hint.aligned_q_len < MIN_CLIP_FOR_JC {
        return Vec::new();
    }
    let clip: Vec<u8> = hint
        .clip_seq
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let clip_rc = rc_dna(&clip);
    let mut out = Vec::new();
    for (idx, rec) in fasta.iter().enumerate() {
        let hay: Vec<u8> = rec.seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        let hits = find_exact(&hay, &clip)
            .into_iter()
            .map(|h| (h, false))
            .chain(find_exact(&hay, &clip_rc).into_iter().map(|h| (h, true)));
        for (hit_0, is_rc) in hits {
            if is_continuation(hint, idx, hit_0, clip.len()) {
                continue;
            }
            // content_rc_hit: the reference-oriented clip content reads along
            // the REVERSE strand of the hit region. The junction-adjacent base
            // is the hit start iff the clip side and the content orientation
            // agree, otherwise the hit end.
            let content_rc_hit = is_rc ^ hint.minus;
            let side2_pos_1 = if hint.clip_is_left == content_rc_hit {
                hit_0 as u64 + 1
            } else {
                hit_0 as u64 + clip.len() as u64
            };
            out.push(SplitCandidate {
                contig_idx: hint.contig_idx,
                side2_contig_idx: idx,
                minus: hint.minus,
                side1_minus: hint.minus,
                side2_minus: content_rc_hit,
                left: SubAlignment {
                    read_start: 0,
                    read_end: hint.aligned_q_len,
                },
                right: SubAlignment {
                    read_start: 0,
                    read_end: clip.len(),
                },
                side1_pos_1: hint.aligned_pos_1,
                side2_pos_1,
                overlap: 0,
            });
        }
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
        let placed = place_softclip(&clips[0], &fasta);
        assert!(
            placed
                .iter()
                .any(|s| (41..=60).contains(&s.side2_pos_1) || (101..=120).contains(&s.side2_pos_1)),
            "expected a hit on an IS copy, got {placed:?}"
        );
    }

    #[test]
    fn minus_strand_left_clip_keys_at_last_match() {
        // Minus read, 5′ clip on the read (CIGAR S+M): in reference
        // orientation the clip lies at HIGHER coordinates, so the junction
        // base is the LAST aligned reference base and clip_is_left = false.
        let ref_seq = vec![b'C'; 160];
        let mut seq = rc_dna(b"TGCAAGTCGATCGTTAGCCA");
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
            clips[0].aligned_pos_1, 30,
            "minus 5′ clip must key at last_match_ref+1 (reference orientation)"
        );
        assert!(
            !clips[0].clip_is_left,
            "minus 5′ clip lies at higher reference coordinates"
        );
        assert_eq!(clips[0].clip_seq, rc_dna(b"TGCAAGTCGATCGTTAGCCA"));
    }

    #[test]
    fn minus_strand_right_clip_keys_at_first_match() {
        // Minus read, 3′ clip on the read (CIGAR M+S): in reference
        // orientation the clip lies at LOWER coordinates, so the junction
        // base is the FIRST aligned reference base and clip_is_left = true.
        // Models an insert right boundary: SEQ = rc([motif][ref 11..30]).
        let ref_seq = vec![b'C'; 160];
        let mut seq = vec![b'C'; 20];
        seq.extend_from_slice(&rc_dna(b"TGCAAGTCGATCGTTAGCCA"));
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
            clips[0].aligned_pos_1, 11,
            "minus 3′ clip must key at first_match_ref+1 (reference orientation)"
        );
        assert!(
            clips[0].clip_is_left,
            "minus 3′ clip lies at lower reference coordinates"
        );
    }
}
