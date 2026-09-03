//! Single-sample consensus engine: Bowtie2 → BAM (noodles) → RA / MC / JC → Genome Diff.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind as NoodlesKind;
use prokadiff_gd::{GdEntry, GdKind, GenomeDiff};
use rayon::prelude::*;

use crate::align::{align_to_bam, AlignKind, FastqInput};
use crate::error::{EvidenceError, Result};
use crate::fasta::{read_reference, FastaRecord};
use crate::jc::{accept_junction, JunctionSupport};
use crate::mc::{
    call_missing_coverage, promote_mc_to_del_with_jc, true_missing_core, MC_DEL_MIN_LEN,
};
use crate::normalize::{right_align_del, right_align_ins};
use crate::pileup::{
    apply_read, place_softclips, AlignedRead, CigarKind, CigarOp, SplitCandidate, SplitOrigin,
};
use crate::ra::{call_consensus, ConsensusCall, PileupColumn, RaOptions};

#[derive(Clone, Debug)]
pub struct EngineOptions {
    pub threads: usize,
    pub keep_bam: bool,
    pub ra: RaOptions,
    pub mc_min_len: usize,
    /// Unique-mapping + total-depth-0 gap must be at least this long to become
    /// a consensus DEL without JC support (see `docs/parity.md`).
    pub mc_del_min_len: usize,
}

type ContigPileup = (Vec<PileupColumn>, Vec<u32>, Vec<u32>, Vec<SplitCandidate>);

/// Two split candidates whose side1 AND side2 both land within this many bp
/// are clustered into one junction before `accept_junction`. IS terminal /
/// TSD register ambiguity keys reads of one event a few bp apart (measured
/// in the synth_is_mob 2371412 BAM: plus reads of the same junction end at
/// 601 and 602; breseq's own two JC for the event land at 599 and 602,
/// 3 bp apart), so exact-key aggregation can strand reads of one junction
/// in different buckets. Tied to `jc_tol_bp=5` in
/// `scripts/layer2_gd_compare.py` (oracle-side matching tolerance); see
/// `docs/parity.md`.
pub const JC_CLUSTER_TOL_BP: u64 = 5;

/// One junction cluster: representative coordinates are those of the FIRST
/// member (smallest side1/side2 after the deterministic sort; the
/// breseq-side comparison matches greedily within ±5 bp, so any member
/// coordinate would do — first-member is the simplest deterministic choice).
/// `JunctionSupport` is aggregated over ALL members, merging strands across
/// the few-bp register split.
struct JcCluster {
    rep_side1: u64,
    rep_side2: u64,
    rep_overlap: i64,
    rep_side1_minus: bool,
    rep_side2_minus: bool,
    support: JunctionSupport,
    /// (hint contig, origin) of each supporting read. Members from the same
    /// underlying read share their id even when the clip placed onto several
    /// repeat copies. The origin variant is part of the key so a
    /// second-pass qname hash can never alias a first-pass hint index.
    ids: HashSet<(usize, SplitOrigin)>,
    /// True if any member came from a long CIGAR I/D. Such clusters have
    /// unique identity and are never folded by the multi-copy rule.
    has_cigar: bool,
}

impl JcCluster {
    fn new(sp: &SplitCandidate) -> Self {
        let mut c = Self {
            rep_side1: sp.side1_pos_1,
            rep_side2: sp.side2_pos_1,
            rep_overlap: sp.overlap,
            rep_side1_minus: sp.side1_minus,
            rep_side2_minus: sp.side2_minus,
            support: JunctionSupport::default(),
            ids: HashSet::new(),
            has_cigar: false,
        };
        c.add(sp);
        c
    }

    fn add(&mut self, sp: &SplitCandidate) {
        match sp.origin {
            SplitOrigin::Cigar => self.has_cigar = true,
            SplitOrigin::Softclip(_) | SplitOrigin::SecondPass(_) => {
                self.ids.insert((sp.contig_idx, sp.origin));
            }
        }
        // Real per-side extensions. Do NOT clamp: `jc::accept_junction`
        // rule 4 rejects a junction whose smallest per-side extension is
        // below 3 bp, and a floor here would make that rule vacuous.
        let ov1 = sp.left.len();
        let ov2 = sp.right.len();
        // Best-read semantics: one unbalanced read must not kill the
        // junction (published accept rules are max-of-min per read).
        let m = ov1.min(ov2);
        let e = &mut self.support;
        e.best_min_overlap = e.best_min_overlap.max(m);
        if sp.minus {
            e.minus_reads += 1;
            e.minus_best_min = e.minus_best_min.max(m);
        } else {
            e.plus_reads += 1;
            e.plus_best_min = e.plus_best_min.max(m);
        }
        e.min_overlap_side1 = if e.min_overlap_side1 == 0 {
            ov1
        } else {
            e.min_overlap_side1.min(ov1)
        };
        e.min_overlap_side2 = if e.min_overlap_side2 == 0 {
            ov2
        } else {
            e.min_overlap_side2.min(ov2)
        };
    }
}

/// An accepted junction cluster plus the identity of its supporting
/// evidence, carried through the redundancy folding below.
struct AcceptedJc {
    c1: usize,
    p1: u64,
    m1: bool,
    c2: usize,
    p2: u64,
    m2: bool,
    overlap: i64,
    support: JunctionSupport,
    /// (hint contig, origin) of each supporting read; CIGAR-split members
    /// contribute no id and set `has_cigar` instead.
    ids: HashSet<(usize, SplitOrigin)>,
    has_cigar: bool,
}

/// (contig, 1-based pos, minus) of one junction side.
type SideKey = (usize, u64, bool);
/// Junction direction key: (side1, side2).
type DirKey = (SideKey, SideKey);

fn dir_key(j: &AcceptedJc) -> DirKey {
    ((j.c1, j.p1, j.m1), (j.c2, j.p2, j.m2))
}

/// Canonical direction: lexicographically smaller of the forward key and
/// the swapped-sides key. A side's strand is the direction its sequence
/// extends away from the junction — intrinsic to the junction geometry,
/// not to the reporting read — so swapping sides keeps each strand
/// attached to its position. (Measured, synth_is_mob 2372015: the reverse
/// report `2560(-1)->601(+1)` duplicates `599(+1)->2558(-1)` with strands
/// preserved per position, shifted only by the ±2 bp register.)
fn canonical_dir_key(j: &AcceptedJc) -> DirKey {
    let f = dir_key(j);
    let r = (f.1, f.0);
    f.min(r)
}

/// Same physical event seen from both sides: equal contigs and strands,
/// both coordinates within JC_CLUSTER_TOL_BP (register shifts).
fn dir_key_within_tol(a: &DirKey, b: &DirKey) -> bool {
    a.0 .0 == b.0 .0
        && a.0 .2 == b.0 .2
        && a.0 .1.abs_diff(b.0 .1) <= JC_CLUSTER_TOL_BP
        && a.1 .0 == b.1 .0
        && a.1 .2 == b.1 .2
        && a.1 .1.abs_diff(b.1 .1) <= JC_CLUSTER_TOL_BP
}

fn absorb_support(a: &mut JunctionSupport, b: &JunctionSupport) {
    a.plus_reads += b.plus_reads;
    a.minus_reads += b.minus_reads;
    a.best_min_overlap = a.best_min_overlap.max(b.best_min_overlap);
    a.plus_best_min = a.plus_best_min.max(b.plus_best_min);
    a.minus_best_min = a.minus_best_min.max(b.minus_best_min);
    // 0 is the unset sentinel (see JcCluster::add).
    a.min_overlap_side1 = match (a.min_overlap_side1, b.min_overlap_side1) {
        (0, y) => y,
        (x, 0) => x,
        (x, y) => x.min(y),
    };
    a.min_overlap_side2 = match (a.min_overlap_side2, b.min_overlap_side2) {
        (0, y) => y,
        (x, 0) => x,
        (x, y) => x.min(y),
    };
}

fn support_total(s: &JunctionSupport) -> usize {
    s.plus_reads + s.minus_reads
}

/// Merge `j` into group representative `g`, keeping the better-supported
/// report's own emitted coordinates (ties keep the existing
/// representative, i.e. the smallest canonical key since the input is
/// sorted). Reads and ids of the folded-away report are retained.
fn absorb_into_rep(g: &mut AcceptedJc, j: AcceptedJc) {
    if support_total(&j.support) > support_total(&g.support) {
        let mut rep = j;
        absorb_support(&mut rep.support, &g.support);
        rep.ids.extend(g.ids.drain());
        rep.has_cigar |= g.has_cigar;
        *g = rep;
    } else {
        absorb_support(&mut g.support, &j.support);
        g.ids.extend(j.ids);
        g.has_cigar |= j.has_cigar;
    }
}

/// Fold reverse-direction duplicates: one physical junction reported once
/// per side (a read aligned on side1 with its clip placed on side2, and
/// another read aligned on side2 with its clip placed on side1). Reports
/// of the same event match after swapping sides within JC_CLUSTER_TOL_BP.
/// One report is kept per event: the better-supported one (ties: smallest
/// canonical key), emitted in its own orientation — the oracle-side
/// comparison canonicalizes sides itself and matches within ±5 bp, so any
/// register/orientation of the same event compares equal.
/// Deterministic: the input is sorted before the greedy single pass.
fn fold_reverse_duplicates(mut accepted: Vec<AcceptedJc>) -> Vec<AcceptedJc> {
    accepted.sort_by_key(canonical_dir_key);
    let mut folded: Vec<AcceptedJc> = Vec::new();
    for j in accepted {
        let ck = canonical_dir_key(&j);
        let mut hit = None;
        for (i, g) in folded.iter().enumerate() {
            if dir_key_within_tol(&canonical_dir_key(g), &ck) {
                hit = Some(i);
                break;
            }
        }
        match hit {
            Some(i) => absorb_into_rep(&mut folded[i], j),
            None => folded.push(j),
        }
    }
    folded
}

/// One shared junction side (same contig, coordinate within 10 bp tolerance),
/// in either emitted orientation.
fn shares_one_side(a: &AcceptedJc, b: &AcceptedJc) -> bool {
    let sa = [(a.c1, a.p1), (a.c2, a.p2)];
    let sb = [(b.c1, b.p1), (b.c2, b.p2)];
    sa.iter()
        .any(|x| sb.iter().any(|y| x.0 == y.0 && x.1.abs_diff(y.1) <= 10))
}

/// Fold multi-copy placement duplicates: a read whose clipped bases have
/// several exact placement hits (repeat copies) contributes one candidate
/// per hit, so the same junction appears once per copy. Two accepted
/// junctions sharing one side and >50% of the smaller supporting-read
/// identity set are one event. (The identity set is the real
/// discriminator: two genuinely distinct junctions can share a repeat
/// copy as one side, but never the same underlying reads.) The emitted
/// placement is the one with the most supporting reads — measured,
/// synth_is_mob 2372015: the copy breseq kept is the one carrying extra
/// copy-anchored reads — ties broken by canonical-key order.
/// CIGAR-split junctions have unique identity and are never folded by
/// this rule. Deterministic: sorted before the greedy single pass.
fn fold_multicopy_placements(mut junctions: Vec<AcceptedJc>) -> Vec<AcceptedJc> {
    junctions.sort_by_key(canonical_dir_key);
    let mut groups: Vec<AcceptedJc> = Vec::new();
    for j in junctions {
        let mut hit = None;
        for (i, g) in groups.iter().enumerate() {
            if g.has_cigar || j.has_cigar {
                continue;
            }
            if !shares_one_side(g, &j) {
                continue;
            }
            let shared = g.ids.intersection(&j.ids).count();
            let smaller = g.ids.len().min(j.ids.len());
            if smaller == 0 || shared * 2 <= smaller {
                continue;
            }
            hit = Some(i);
            break;
        }
        match hit {
            Some(i) => absorb_into_rep(&mut groups[i], j),
            None => groups.push(j),
        }
    }
    groups
}

impl Default for EngineOptions {
    fn default() -> Self {
        Self {
            threads: 8,
            keep_bam: false,
            ra: RaOptions::default(),
            mc_min_len: 3,
            mc_del_min_len: MC_DEL_MIN_LEN,
        }
    }
}

/// Align `reads` to `ref_fa` and write a Genome Diff to `outdir/output.gd`.
pub fn run_sample(
    ref_fa: &Path,
    reads: &FastqInput,
    outdir: &Path,
    opts: &EngineOptions,
) -> Result<PathBuf> {
    std::fs::create_dir_all(outdir)?;
    let work = outdir.join("work");
    std::fs::create_dir_all(&work)?;
    let bam_path = outdir.join("aligned.bam");
    eprintln!("prokadiff: primary alignment");
    align_to_bam(
        ref_fa,
        reads,
        &bam_path,
        opts.threads,
        &work,
        AlignKind::Primary,
    )?;
    let fasta = read_reference(ref_fa)?;
    eprintln!(
        "prokadiff: pileup + junction seeds (place S>={})",
        crate::jc_seq::MIN_CLIP_FOR_PLACE
    );
    let aligned = read_aligned_bam(&bam_path, &fasta)?;
    let contig_results = pileup_contigs(&fasta, &aligned, opts);
    eprintln!("prokadiff: candidate-junction second pass");
    let extra = second_pass_splits(
        &fasta,
        reads,
        &bam_path,
        &contig_results,
        &aligned,
        opts,
        &work,
    )?;
    let gd = emit_from_pileup(&fasta, contig_results, opts, &extra);
    let gd_path = outdir.join("output.gd");
    gd.write_path(&gd_path)?;
    if !opts.keep_bam {
        let _ = std::fs::remove_file(&bam_path);
        let _ = std::fs::remove_dir_all(&work);
    }
    Ok(gd_path)
}

pub fn call_from_bam(
    bam_path: &Path,
    fasta: &[FastaRecord],
    opts: &EngineOptions,
) -> Result<GenomeDiff> {
    let aligned = read_aligned_bam(bam_path, fasta)?;
    Ok(call_from_aligned(fasta, &aligned, opts))
}

fn read_aligned_bam(bam_path: &Path, fasta: &[FastaRecord]) -> Result<Vec<AlignedRead>> {
    let name_to_idx: HashMap<&str, usize> = fasta
        .iter()
        .enumerate()
        .map(|(i, r)| (r.name.as_str(), i))
        .collect();

    let mut reader = File::open(bam_path).map(bam::io::Reader::new)?;
    let header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(n, _)| n.to_string())
        .collect();

    let mut aligned = Vec::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        if rec.flags().is_unmapped() || rec.flags().is_secondary() {
            continue;
        }
        let Some(Ok(rid)) = rec.reference_sequence_id() else {
            continue;
        };
        let contig_name = ref_names.get(rid).map(String::as_str).unwrap_or("");
        let Some(&contig_idx) = name_to_idx.get(contig_name) else {
            continue;
        };
        let Some(Ok(start)) = rec.alignment_start() else {
            continue;
        };
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(0);
        let minus = rec.flags().is_reverse_complemented();
        let seq_len = rec.sequence().len();
        let mut seq = Vec::with_capacity(seq_len);
        for i in 0..seq_len {
            if let Some(b) = rec.sequence().get(i) {
                seq.push(b);
            }
        }
        let cigar = noodles_cigar(&rec)?;
        aligned.push(AlignedRead {
            contig_idx,
            ref_start_0: start.get() as i64 - 1,
            minus,
            seq,
            cigar,
            mapq,
        });
    }
    Ok(aligned)
}

fn second_pass_read_names(bam_path: &Path) -> Result<(HashSet<String>, HashSet<String>)> {
    use crate::align::{keep_for_second_pass, normalize_qname};

    let mut reader = File::open(bam_path).map(bam::io::Reader::new)?;
    let _header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    let mut keep = HashSet::new();
    let mut seen = HashSet::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        if rec.flags().is_secondary() {
            continue;
        }
        let name = rec.name().map(|n| n.to_string()).unwrap_or_default();
        if name.is_empty() {
            continue;
        }
        let q = normalize_qname(&name).to_string();
        seen.insert(q.clone());
        let cigar = noodles_cigar(&rec)?;
        let max_s = cigar
            .iter()
            .filter(|op| op.kind == CigarKind::SoftClip)
            .map(|op| op.len)
            .max()
            .unwrap_or(0);
        if keep_for_second_pass(
            rec.flags().is_unmapped(),
            rec.flags().is_mate_unmapped(),
            max_s,
        ) {
            keep.insert(q);
        }
    }
    Ok((keep, seen))
}

fn cigar_counts(cigar: &[CigarOp]) -> (usize, usize, usize, usize) {
    let mut n_match = 0usize;
    let mut n_ins = 0usize;
    let mut n_del = 0usize;
    for op in cigar {
        match op.kind {
            CigarKind::Match => n_match += op.len,
            CigarKind::Ins => n_ins += op.len,
            CigarKind::Del => n_del += op.len,
            _ => {}
        }
    }
    let q_cover = n_match + n_ins;
    (n_match, n_ins, n_del, q_cover)
}

fn construct_span(ref_start_0: i64, cigar: &[CigarOp]) -> (usize, usize) {
    let mut pos = ref_start_0.max(0) as usize;
    let start = pos;
    for op in cigar {
        match op.kind {
            CigarKind::Match | CigarKind::Del | CigarKind::Skip => pos += op.len,
            CigarKind::Ins | CigarKind::SoftClip | CigarKind::HardClip => {}
        }
    }
    (start, pos)
}

fn noodles_cigar(rec: &noodles::bam::Record) -> Result<Vec<CigarOp>> {
    let mut cigar = Vec::new();
    for op in rec.cigar().iter() {
        let op = op.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        let kind = match op.kind() {
            NoodlesKind::Match | NoodlesKind::SequenceMatch | NoodlesKind::SequenceMismatch => {
                CigarKind::Match
            }
            NoodlesKind::Insertion => CigarKind::Ins,
            NoodlesKind::Deletion => CigarKind::Del,
            NoodlesKind::SoftClip => CigarKind::SoftClip,
            NoodlesKind::HardClip => CigarKind::HardClip,
            NoodlesKind::Skip => CigarKind::Skip,
            _ => continue,
        };
        cigar.push(CigarOp {
            kind,
            len: op.len(),
        });
    }
    Ok(cigar)
}

fn primary_alignment_scores(bam_path: &Path) -> Result<HashMap<String, i32>> {
    let mut reader = File::open(bam_path).map(bam::io::Reader::new)?;
    let _header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    let mut scores = HashMap::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        if rec.flags().is_unmapped() || rec.flags().is_secondary() {
            continue;
        }
        let name = rec.name().map(|n| n.to_string()).unwrap_or_default();
        if name.is_empty() {
            continue;
        }
        let norm_name = crate::align::normalize_qname(&name).to_string();
        let cigar = noodles_cigar(&rec)?;
        let (m, i, d, _) = cigar_counts(&cigar);
        let score = crate::jc_seq::cigar_match_indel_score(m, i, d);
        scores
            .entry(norm_name)
            .and_modify(|s: &mut i32| *s = (*s).max(score))
            .or_insert(score);
    }
    Ok(scores)
}

/// Minimum construct bases a second-pass alignment must place on EACH side
/// of the breakpoint to count as junction evidence at all.
///
/// Matched to `jc::accept_junction` rule 4 (smallest per-side extension
/// ≥3 bp), which is a min over ALL supporting reads while rules 2 and 3 are
/// max-of-min (best read / best per strand). A read placing 1–2 bp past the
/// breakpoint is noise, so it is excluded here rather than admitted and then
/// allowed to drag the cluster minimum under 3 and reject a junction that a
/// dozen well-balanced reads support. The engine used to admit such reads and
/// clamp their extensions up to 3, which fabricated evidence and made rule 4
/// unreachable; excluding them keeps the rule live without inventing bases.
const JC_MIN_SPAN_EACH_SIDE: usize = 3;

/// Why a second-pass alignment did not become junction support. Counted so a
/// single Clonal job can tell "no spanning reads" from "the score gate ate
/// them" without guessing (JC=0 root-cause work, `benchmark/results/`).
#[derive(Clone, Copy, Debug, Default)]
struct SecondPassReject {
    considered: usize,
    short_cover: usize,
    not_spanning: usize,
    worse_than_primary: usize,
    kept: usize,
}

/// Hash of a read name, used as second-pass supporting-read identity.
fn qname_id(qname: &str) -> u64 {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    let mut hasher = DefaultHasher::new();
    qname.hash(&mut hasher);
    hasher.finish()
}

/// Promote spanning second-pass alignments to `SplitCandidate`s.
///
/// `candidates` pairs each retained candidate junction with the breakpoint
/// offset of ITS construct (`JunctionConstruct::breakpoint`), which is below
/// `flank` whenever a flank was clamped at a contig boundary.
fn parse_junction_bam(
    bam_path: &Path,
    candidates: &[(crate::jc_seq::CandidateJunction, usize)],
    ref_scores: &HashMap<String, i32>,
) -> Result<(Vec<SplitCandidate>, SecondPassReject)> {
    use crate::jc::SubAlignment;
    use crate::jc_seq::{spans_breakpoint, JC_MIN_COVER_BASES};

    let mut stats = SecondPassReject::default();
    let mut reader = File::open(bam_path).map(bam::io::Reader::new)?;
    let header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(n, _)| n.to_string())
        .collect();
    let mut extra = Vec::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        if rec.flags().is_unmapped() {
            continue;
        }
        let Some(Ok(rid)) = rec.reference_sequence_id() else {
            continue;
        };
        let Some(name) = ref_names.get(rid) else {
            continue;
        };
        let Some(idx_str) = name.strip_prefix("jc") else {
            continue;
        };
        let Ok(ci) = idx_str.parse::<usize>() else {
            continue;
        };
        let Some((cand, breakpoint)) = candidates.get(ci) else {
            continue;
        };
        let breakpoint = *breakpoint;
        let Some(Ok(start)) = rec.alignment_start() else {
            continue;
        };
        stats.considered += 1;
        let cigar = noodles_cigar(&rec)?;
        let (n_match, n_ins, n_del, q_cover) = cigar_counts(&cigar);
        if q_cover < JC_MIN_COVER_BASES {
            stats.short_cover += 1;
            continue;
        }
        let (cstart, cend) = construct_span(start.get() as i64 - 1, &cigar);
        if !spans_breakpoint(cstart, cend, breakpoint, JC_MIN_SPAN_EACH_SIDE) {
            stats.not_spanning += 1;
            continue;
        }
        let jc_score = crate::jc_seq::cigar_match_indel_score(n_match, n_ins, n_del);
        let qname = rec.name().map(|n| n.to_string()).unwrap_or_default();
        let norm_qname = crate::align::normalize_qname(&qname);
        if let Some(&rs) = ref_scores.get(norm_qname) {
            if jc_score < rs {
                stats.worse_than_primary += 1;
                continue;
            }
        }
        // Real per-side extension. Do NOT clamp: `jc::accept_junction` rule 4
        // requires ≥3 bp on the smallest side, and a floor here would make
        // that rule vacuous.
        let ov1 = breakpoint.saturating_sub(cstart);
        let ov2 = cend.saturating_sub(breakpoint);
        stats.kept += 1;
        extra.push(SplitCandidate {
            contig_idx: cand.side1_contig,
            side2_contig_idx: cand.side2_contig,
            minus: rec.flags().is_reverse_complemented(),
            side1_minus: cand.side1_minus,
            side2_minus: cand.side2_minus,
            left: SubAlignment {
                read_start: 0,
                read_end: ov1,
            },
            right: SubAlignment {
                read_start: ov1,
                read_end: ov1 + ov2,
            },
            side1_pos_1: cand.side1_pos_1,
            side2_pos_1: cand.side2_pos_1,
            overlap: cand.overlap,
            origin: SplitOrigin::SecondPass(qname_id(&qname)),
        });
    }
    Ok((extra, stats))
}

fn cluster_key_splits(splits: &[SplitCandidate]) -> Vec<((usize, usize, i8), JcCluster)> {
    let mut groups: HashMap<(usize, usize, i8), Vec<&SplitCandidate>> = HashMap::new();
    for sp in splits {
        groups
            .entry((
                sp.contig_idx,
                sp.side2_contig_idx,
                sp.overlap.signum() as i8,
            ))
            .or_default()
            .push(sp);
    }
    let mut group_keys: Vec<(usize, usize, i8)> = groups.keys().copied().collect();
    group_keys.sort_unstable();
    let mut out = Vec::new();
    for gk in group_keys {
        let mut members: Vec<&SplitCandidate> = groups[&gk].clone();
        members.sort_by_key(|sp| (sp.side1_pos_1, sp.side2_pos_1));
        let mut clusters: Vec<JcCluster> = Vec::new();
        for sp in members {
            // Clusters are appended in non-decreasing order of rep_side1.
            // Any cluster with c.rep_side1 + JC_CLUSTER_TOL_BP < sp.side1_pos_1 can never match
            // this or any future member; prune them from the search.
            let start_idx = clusters.partition_point(|c| {
                c.rep_side1.saturating_add(JC_CLUSTER_TOL_BP) < sp.side1_pos_1
            });
            let mut hit = None;
            for (rel_i, c) in clusters[start_idx..].iter().enumerate() {
                if sp.side1_pos_1.abs_diff(c.rep_side1) <= JC_CLUSTER_TOL_BP
                    && sp.side2_pos_1.abs_diff(c.rep_side2) <= JC_CLUSTER_TOL_BP
                {
                    hit = Some(start_idx + rel_i);
                    break;
                }
            }
            match hit {
                Some(i) => clusters[i].add(sp),
                None => clusters.push(JcCluster::new(sp)),
            }
        }
        for c in clusters {
            out.push((gk, c));
        }
    }
    out
}

fn seed_candidates_from_pileup(
    contig_results: &[ContigPileup],
) -> Vec<crate::jc_seq::CandidateJunction> {
    let mut cands = Vec::new();
    for (_, _, _, splits) in contig_results {
        for (gk, c) in cluster_key_splits(splits) {
            let n = support_total(&c.support);
            if n == 0 {
                continue;
            }
            cands.push(crate::jc_seq::CandidateJunction {
                side1_contig: gk.0,
                side1_pos_1: c.rep_side1,
                side1_minus: c.rep_side1_minus,
                side2_contig: gk.1,
                side2_pos_1: c.rep_side2,
                side2_minus: c.rep_side2_minus,
                overlap: c.rep_overlap,
                seed_reads: n,
            });
        }
    }
    cands
}

fn second_pass_splits(
    fasta: &[FastaRecord],
    reads: &FastqInput,
    primary_bam: &Path,
    contig_results: &[ContigPileup],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
    work: &Path,
) -> Result<Vec<SplitCandidate>> {
    use crate::fasta::write_records;
    use crate::jc_seq::{junction_construct, rank_and_cap, CandidateJunction, JC_MIN_COVER_BASES};

    let flank = aligned
        .iter()
        .map(|r| r.seq.len())
        .max()
        .unwrap_or(36)
        .max(28);
    let mut all_cands = seed_candidates_from_pileup(contig_results);

    // Extract split-read candidates from Stage 2 sensitive alignments if present
    let contig_names: Vec<String> = fasta.iter().map(|r| r.name.clone()).collect();
    for (sam_name, default_mate) in [
        ("stage2_r1.sam", 1),
        ("stage2_r2.sam", 2),
        ("stage2_pe.sam", 0),
        ("stage2_se.sam", 0),
        ("stage2.sam", 0),
    ] {
        let sam_path = work.join(sam_name);
        if sam_path.is_file() {
            if let Ok(split_cands) = crate::split_seed::extract_candidate_junctions_from_sam(
                &sam_path,
                &contig_names,
                default_mate,
            ) {
                eprintln!(
                    "prokadiff: extracted {} split-read candidate junctions from {}",
                    split_cands.len(),
                    sam_name
                );
                all_cands.extend(split_cands);
            }
        }
    }

    // Merge duplicate candidate junctions across sources (R1, R2, pileup clips)
    // Canonicalize key so (A, B) and (B, A) aggregate their seed_reads.
    type JcMergeKey = (usize, u64, bool, usize, u64, bool, i64);
    let mut merged_cands: HashMap<JcMergeKey, CandidateJunction> = HashMap::new();
    for cand in all_cands {
        let (c1, p1, m1, c2, p2, m2) =
            if (cand.side1_contig, cand.side1_pos_1) <= (cand.side2_contig, cand.side2_pos_1) {
                (
                    cand.side1_contig,
                    cand.side1_pos_1,
                    cand.side1_minus,
                    cand.side2_contig,
                    cand.side2_pos_1,
                    cand.side2_minus,
                )
            } else {
                (
                    cand.side2_contig,
                    cand.side2_pos_1,
                    !cand.side2_minus,
                    cand.side1_contig,
                    cand.side1_pos_1,
                    !cand.side1_minus,
                )
            };
        let key = (c1, p1, m1, c2, p2, m2, cand.overlap);
        merged_cands
            .entry(key)
            .and_modify(|e| e.seed_reads += cand.seed_reads)
            .or_insert(CandidateJunction {
                side1_contig: c1,
                side1_pos_1: p1,
                side1_minus: m1,
                side2_contig: c2,
                side2_pos_1: p2,
                side2_minus: m2,
                overlap: cand.overlap,
                seed_reads: cand.seed_reads,
            });
    }
    let all_cands: Vec<CandidateJunction> = merged_cands.into_values().collect();

    let ref_len: usize = fasta.iter().map(|r| r.seq.len()).sum();
    let cands = rank_and_cap(all_cands, ref_len, flank);
    if cands.is_empty() {
        eprintln!("prokadiff: second-pass skipped (no seed junctions)");
        return Ok(Vec::new());
    }
    let mut recs = Vec::new();
    // Each construct carries its OWN breakpoint: `junction_flank` clamps at
    // contig boundaries, so a construct seeded near a contig end has its
    // breakpoint below `flank`. Using `flank` there mislocates the junction
    // and corrupts the spanning / per-side-overlap checks downstream.
    let mut kept: Vec<(CandidateJunction, usize)> = Vec::new();
    for cand in &cands {
        let construct = junction_construct(fasta, cand, flank);
        if construct.len() < JC_MIN_COVER_BASES {
            continue;
        }
        recs.push(FastaRecord {
            name: format!("jc{}", kept.len()),
            seq: construct.sequence(),
        });
        kept.push((cand.clone(), construct.breakpoint()));
    }
    if recs.is_empty() {
        return Ok(Vec::new());
    }
    let jc_dir = work.join("jc");
    std::fs::create_dir_all(&jc_dir)?;
    let fa = jc_dir.join("junctions.fa");
    write_records(&recs, &fa)?;
    eprintln!(
        "prokadiff: second-pass {} junction constructs; filtering clip/unmapped reads",
        recs.len()
    );
    let gzip = reads.files.iter().any(|p| {
        p.file_name()
            .and_then(|n| n.to_str())
            .is_some_and(|n| n.ends_with(".gz"))
    });
    let filtered_owner;
    let pass_reads: &FastqInput = if gzip {
        eprintln!("prokadiff: gzip FASTQ; second-pass uses all reads (no clip/unmapped filter)");
        reads
    } else {
        let (keep, seen) = second_pass_read_names(primary_bam)?;
        let pass_dir = jc_dir.join("reads");
        let (filtered, nkeep) = crate::align::filter_fastq_input(reads, &keep, &seen, &pass_dir)?;
        eprintln!("prokadiff: second-pass {nkeep} FASTQ records after clip/unmapped filter");
        if nkeep == 0 {
            return Ok(Vec::new());
        }
        filtered_owner = filtered;
        &filtered_owner
    };
    let jc_bam = jc_dir.join("aligned.bam");
    align_to_bam(
        &fa,
        pass_reads,
        &jc_bam,
        opts.threads,
        &jc_dir,
        AlignKind::Junction,
    )?;
    let ref_scores = primary_alignment_scores(primary_bam)?;
    let (extra, stats) = parse_junction_bam(&jc_bam, &kept, &ref_scores)?;
    eprintln!(
        "prokadiff: second-pass records considered={} kept={} \
         rejected(short_cover={} not_spanning={} worse_than_primary={})",
        stats.considered,
        stats.kept,
        stats.short_cover,
        stats.not_spanning,
        stats.worse_than_primary
    );
    Ok(extra)
}

/// In-memory consensus (no BAM). Used by `call_from_bam` and layer-0 tests.
pub fn call_from_aligned(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
) -> GenomeDiff {
    call_from_aligned_extra(fasta, aligned, opts, &[])
}

/// Like `call_from_aligned`, then merge extra split candidates (second-pass
/// junction hits) before clustering and `accept_junction`.
pub fn call_from_aligned_extra(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
    extra_splits: &[SplitCandidate],
) -> GenomeDiff {
    let contig_results = pileup_contigs(fasta, aligned, opts);
    emit_from_pileup(fasta, contig_results, opts, extra_splits)
}

fn pileup_contigs(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
) -> Vec<ContigPileup> {
    let pool = if opts.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(opts.threads)
            .build()
            .ok()
    } else {
        None
    };

    let call_contigs = || {
        fasta
            .par_iter()
            .enumerate()
            .map(|(idx, rec)| {
                let n = rec.seq.len();
                let mut columns: Vec<PileupColumn> = rec
                    .seq
                    .iter()
                    .map(|&b| PileupColumn {
                        ref_base: b.to_ascii_uppercase(),
                        observations: Vec::new(),
                        insertions: Vec::new(),
                    })
                    .collect();
                let mut unique_depth = vec![0u32; n];
                let mut total_depth = vec![0u32; n];
                let mut splits = Vec::new();
                let mut clips = Vec::new();
                for read in aligned.iter().filter(|r| r.contig_idx == idx) {
                    apply_read(
                        read,
                        &mut columns,
                        &mut unique_depth,
                        &mut total_depth,
                        &mut splits,
                        &mut clips,
                    );
                }
                splits.extend(place_softclips(&clips, fasta));
                (columns, unique_depth, total_depth, splits)
            })
            .collect()
    };
    match &pool {
        Some(p) => p.install(call_contigs),
        None => call_contigs(),
    }
}

fn emit_from_pileup(
    fasta: &[FastaRecord],
    mut contig_results: Vec<ContigPileup>,
    opts: &EngineOptions,
    extra_splits: &[SplitCandidate],
) -> GenomeDiff {
    for sp in extra_splits {
        if sp.contig_idx < contig_results.len() {
            contig_results[sp.contig_idx].3.push(sp.clone());
        }
    }

    let mut gd = GenomeDiff::new();
    gd.metadata
        .push(("PROGRAM".into(), "prokadiff-evidence".into()));
    let mut next_id = 1u32;
    // Accepted junctions are collected across ALL contigs and emitted only
    // after redundancy folding: a reverse-direction duplicate keyed on
    // contig B can belong to a physical event first reported on contig A.
    let mut accepted_all: Vec<AcceptedJc> = Vec::new();

    for (idx, rec) in fasta.iter().enumerate() {
        let (columns, unique_depth, total_depth, splits) = &contig_results[idx];
        let mut pending_del: Option<(u64, u8)> = None;
        let flush_del = |gd: &mut GenomeDiff, id: &mut u32, start: u64, size: u8| {
            // Match breseq mutation coordinates (RA evidence may stay 5′).
            let start = right_align_del(&rec.seq, start, u64::from(size));
            gd.entries
                .push(GdEntry::del(*id, rec.name.clone(), start, u64::from(size)));
            *id += 1;
        };

        for (i, col) in columns.iter().enumerate() {
            let pos = i as u64 + 1;
            match call_consensus(col, &opts.ra) {
                ConsensusCall::Snp { alt } => {
                    if let Some((s, sz)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz);
                    }
                    gd.entries.push(GdEntry::snp(
                        next_id,
                        rec.name.clone(),
                        pos,
                        (alt as char).to_string(),
                    ));
                    next_id += 1;
                }
                ConsensusCall::Ins { seq } => {
                    if let Some((s, sz)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz);
                    }
                    let (pos, seq) = right_align_ins(&rec.seq, pos, &seq);
                    let s = String::from_utf8_lossy(&seq).into_owned();
                    gd.entries
                        .push(GdEntry::ins(next_id, rec.name.clone(), pos, s));
                    next_id += 1;
                }
                ConsensusCall::Del { size } => match pending_del {
                    Some((s, sz)) if s + u64::from(sz) == pos && sz + size <= 2 => {
                        pending_del = Some((s, sz + size));
                    }
                    Some((s, sz)) => {
                        flush_del(&mut gd, &mut next_id, s, sz);
                        pending_del = Some((pos, size));
                    }
                    None => pending_del = Some((pos, size)),
                },
                _ => {
                    if let Some((s, sz)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz);
                    }
                }
            }
        }
        if let Some((s, sz)) = pending_del.take() {
            flush_del(&mut gd, &mut next_id, s, sz);
        }

        // Junction clustering (replaces exact-key aggregation). Do not key
        // on read strand: plus and minus must aggregate (accept_junction).
        // Group by (contig, side2 contig, overlap sign) — the exact overlap
        // value is register-dependent, but its sign is biologically distinct
        // (deletion-like < 0, blunt = 0, insertion-like > 0). Within a group,
        // sort by (side1, side2) and greedily cluster: a candidate joins the
        // first cluster whose FIRST member is within JC_CLUSTER_TOL_BP on
        // both sides (single-pass against the first member is deterministic;
        // documented in docs/parity.md).
        let mut jc_endpoints_1 = Vec::new();
        for (gk, c) in cluster_key_splits(splits) {
            if !accept_junction(&c.support) {
                continue;
            }
            if gk.0 == idx {
                jc_endpoints_1.push(c.rep_side1);
            }
            if gk.1 == idx {
                jc_endpoints_1.push(c.rep_side2);
            }
            accepted_all.push(AcceptedJc {
                c1: gk.0,
                p1: c.rep_side1,
                m1: c.rep_side1_minus,
                c2: gk.1,
                p2: c.rep_side2,
                m2: c.rep_side2_minus,
                overlap: c.rep_overlap,
                support: c.support,
                ids: c.ids,
                has_cigar: c.has_cigar,
            });
        }

        for span in call_missing_coverage(unique_depth, opts.mc_min_len) {
            if span.len() < 3 {
                continue;
            }
            // Contig ends look like unique-depth 0; breseq reports UN, not DEL.
            if span.start_0 == 0 || span.end_0 == unique_depth.len() {
                continue;
            }
            let start = span.start_1();
            let end = span.end_1();
            let mc_id = next_id;
            next_id += 1;
            gd.entries
                .push(GdEntry::mc(mc_id, rec.name.clone(), start, end, 0, 0));
            if !promote_mc_to_del_with_jc(
                span,
                unique_depth,
                total_depth,
                opts.mc_del_min_len,
                opts.mc_min_len,
                &jc_endpoints_1,
            ) {
                continue;
            }
            let core = true_missing_core(
                span,
                unique_depth,
                total_depth,
                crate::mc::coverage_cutoff(unique_depth),
            )
            .unwrap_or(span);
            let mut del_start = core.start_1();
            let mut del_end = core.end_1();
            // If an accepted JC flanks this missing coverage span,
            // the JC breakpoint defines the exact deletion boundary:
            for j in &accepted_all {
                if j.c1 == idx && j.c2 == idx {
                    let (jp1, jp2) = if j.p1 <= j.p2 {
                        (j.p1, j.p2)
                    } else {
                        (j.p2, j.p1)
                    };
                    // Shift jp1/jp2 if overlap > 0 (assign exclusively to one side)
                    let (norm_jp1, norm_jp2) = if j.overlap > 0 {
                        let k = j.overlap as u64;
                        if j.m1 {
                            (jp1.saturating_sub(k), jp2)
                        } else if !j.m2 {
                            (jp1, jp2.saturating_add(k))
                        } else {
                            (jp1, jp2)
                        }
                    } else {
                        (jp1, jp2)
                    };
                    if norm_jp1.abs_diff(span.start_1().saturating_sub(1)) <= 50
                        && norm_jp2.abs_diff(span.end_1().saturating_add(1)) <= 50
                    {
                        del_start = norm_jp1 + 1;
                        del_end = norm_jp2.saturating_sub(1);
                        break;
                    }
                }
            }
            let size = del_end.saturating_sub(del_start) + 1;
            let del_start = right_align_del(&rec.seq, del_start, size);
            gd.entries
                .push(GdEntry::del(next_id, rec.name.clone(), del_start, size));
            if let Some(last) = gd.entries.last_mut() {
                last.parent_ids = vec![mc_id];
            }
            next_id += 1;
        }
    }

    // Redundancy folding (synth_is_mob 2372015: 6 JC emitted for 2
    // physical junctions — 2 multi-copy placements + 2 reverse-direction
    // duplicates). Reverse first: folding copy placements before direction
    // duplicates would strand the copy-side reports. Then emit one JC per
    // physical event, in deterministic coordinate order.
    let mut folded = fold_multicopy_placements(fold_reverse_duplicates(accepted_all));
    folded.sort_by_key(|j| (j.c1, j.p1, j.m1, j.c2, j.p2, j.m2));
    for mut j in folded {
        // breseq GD convention: normalize overlap into endpoints when microhomology exists
        if j.overlap > 0 {
            let k = j.overlap as u64;
            if j.m1 {
                j.p1 = j.p1.saturating_sub(k);
                j.overlap = 0;
            } else if !j.m2 {
                j.p2 = j.p2.saturating_add(k);
                j.overlap = 0;
            }
        }
        let s1 = if j.m1 { "-1" } else { "1" };
        let s2 = if j.m2 { "-1" } else { "1" };
        let n1 = fasta[j.c1].name.clone();
        let n2 = fasta[j.c2].name.clone();
        gd.entries
            .push(GdEntry::jc(next_id, n1, j.p1, s1, n2, j.p2, s2, j.overlap));
        next_id += 1;
    }

    // Subsume/mask any SNP that falls inside a called DEL
    let dels: Vec<(String, u64, u64)> = gd
        .entries
        .iter()
        .filter(|e| e.kind == GdKind::Del)
        .filter_map(|e| {
            let seq_id = e.fields.first()?.clone();
            let pos: u64 = e.fields.get(1)?.parse().ok()?;
            let sz: u64 = e.fields.get(2)?.parse().ok()?;
            Some((seq_id, pos, pos + sz.saturating_sub(1)))
        })
        .collect();
    if !dels.is_empty() {
        gd.entries.retain(|e| {
            if e.kind == GdKind::Snp {
                let Some(seq_id) = e.fields.first() else {
                    return true;
                };
                let Some(pos) = e.fields.get(1).and_then(|p| p.parse::<u64>().ok()) else {
                    return true;
                };
                !dels
                    .iter()
                    .any(|(d_seq, start, end)| d_seq == seq_id && pos >= *start && pos <= *end)
            } else {
                true
            }
        });
    }

    gd
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jc::SubAlignment;
    use crate::pileup::{apply_read, UNIQUE_MAPQ};

    #[test]
    fn ra_calls_from_in_memory_pileup_without_bam() {
        let fasta = [FastaRecord {
            name: "chr".into(),
            seq: b"ACGTACGTAC".to_vec(),
        }];
        let n = fasta[0].seq.len();
        let mut columns: Vec<PileupColumn> = fasta[0]
            .seq
            .iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: Vec::new(),
                insertions: Vec::new(),
            })
            .collect();
        let mut depth = vec![0u32; n];
        let mut splits = Vec::new();
        let mut alt = fasta[0].seq.clone();
        alt[2] = b'T';
        for minus in [false, false, false, true, true, true] {
            let read = AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus,
                seq: alt.clone(),
                cigar: vec![CigarOp {
                    kind: CigarKind::Match,
                    len: alt.len(),
                }],
                mapq: UNIQUE_MAPQ,
            };
            apply_read(
                &read,
                &mut columns,
                &mut depth,
                &mut vec![0u32; n],
                &mut splits,
                &mut Vec::new(),
            );
        }
        assert_eq!(
            call_consensus(&columns[2], &RaOptions::default()),
            ConsensusCall::Snp { alt: b'T' }
        );
    }

    fn fasta_chr(seq: &[u8]) -> [FastaRecord; 1] {
        [FastaRecord {
            name: "chr".into(),
            seq: seq.to_vec(),
        }]
    }

    fn covering_reads(
        seq: &[u8],
        start: i64,
        len: usize,
        n_plus: usize,
        n_minus: usize,
        mapq: u8,
    ) -> Vec<AlignedRead> {
        let slice = seq[start as usize..start as usize + len].to_vec();
        let mut out = Vec::new();
        for minus in std::iter::repeat_n(false, n_plus).chain(std::iter::repeat_n(true, n_minus)) {
            out.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: start,
                minus,
                seq: slice.clone(),
                cigar: vec![CigarOp {
                    kind: CigarKind::Match,
                    len,
                }],
                mapq,
            });
        }
        out
    }

    fn gd_dels(gd: &GenomeDiff) -> Vec<&prokadiff_gd::GdEntry> {
        gd.entries
            .iter()
            .filter(|e| e.kind == prokadiff_gd::GdKind::Del)
            .collect()
    }

    fn opts_single_thread() -> EngineOptions {
        EngineOptions {
            threads: 1,
            ..EngineOptions::default()
        }
    }

    #[test]
    fn short_internal_coverage_gap_does_not_emit_del() {
        let seq = vec![b'A'; 40];
        let mut reads = covering_reads(&seq, 0, 15, 6, 6, UNIQUE_MAPQ);
        reads.extend(covering_reads(&seq, 20, 20, 6, 6, UNIQUE_MAPQ));
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(
            gd_dels(&gd).is_empty(),
            "short unique+total gap must stay MC, not consensus DEL; got {:?}",
            gd_dels(&gd).iter().map(|e| &e.fields).collect::<Vec<_>>()
        );
    }

    #[test]
    fn contig_terminal_zero_does_not_emit_del() {
        let seq = vec![b'A'; 40];
        let reads = covering_reads(&seq, 10, 30, 6, 6, UNIQUE_MAPQ);
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(gd_dels(&gd).is_empty());
    }

    #[test]
    fn repeat_like_low_mapq_coverage_is_not_del() {
        let seq = vec![b'A'; 40];
        let mut reads = covering_reads(&seq, 0, 10, 6, 6, UNIQUE_MAPQ);
        reads.extend(covering_reads(&seq, 10, 16, 6, 6, 3)); // MAPQ < unique
        reads.extend(covering_reads(&seq, 26, 14, 6, 6, UNIQUE_MAPQ));
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(
            gd_dels(&gd).is_empty(),
            "unique-depth 0 with remaining total coverage must not be DEL; got {:?}",
            gd_dels(&gd).iter().map(|e| &e.fields).collect::<Vec<_>>()
        );
    }

    #[test]
    fn long_internal_true_gap_emits_del() {
        let seq = vec![b'A'; 80];
        let mut reads = covering_reads(&seq, 0, 10, 6, 6, UNIQUE_MAPQ);
        reads.extend(covering_reads(&seq, 70, 10, 6, 6, UNIQUE_MAPQ));
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let dels = gd_dels(&gd);
        assert_eq!(
            dels.len(),
            1,
            "expected one unique-mapping DEL, got {dels:?}"
        );
        let start: u64 = dels[0].fields[1].parse().unwrap();
        let size: u64 = dels[0].fields[2].parse().unwrap();
        assert!(start > 1, "DEL must not be contig-terminal, start={start}");
        assert!(
            size >= 50,
            "unique+total-zero core must be long, size={size}"
        );
    }

    fn gd_jcs(gd: &GenomeDiff) -> Vec<&prokadiff_gd::GdEntry> {
        gd.entries
            .iter()
            .filter(|e| e.kind == prokadiff_gd::GdKind::Jc)
            .collect()
    }

    #[test]
    fn both_strand_long_deletion_emits_jc() {
        let seq = vec![b'A'; 50];
        let mut reads = Vec::new();
        for minus in [false, false, false, true, true, true] {
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus,
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
                mapq: UNIQUE_MAPQ,
            });
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(
            !gd_jcs(&gd).is_empty(),
            "plus+minus long-D splits must aggregate into an accepted JC; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn softclip_onto_repeat_copy_emits_jc() {
        // Unique flanks + two identical 20 bp IS-like copies (first-period mosaic).
        let mut seq = vec![b'C'; 160];
        let motif = *b"ACGTACGTACGTACGTACGT";
        seq[40..60].copy_from_slice(&motif);
        seq[100..120].copy_from_slice(&motif);
        // Reads crossing the insert left boundary: match ref 71..90, 3′ clip.
        // SAM/BAM stores SEQ reference-forward (minus reads are reverse-
        // complemented), so plus and minus reads crossing the same boundary
        // share the SAME stored SEQ and CIGAR; only the flag differs.
        let mut read_seq = vec![b'C'; 20];
        read_seq.extend_from_slice(&motif);
        let mut reads = Vec::new();
        for minus in [false, false, false, true, true, true] {
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 70,
                minus,
                seq: read_seq.clone(),
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
                mapq: UNIQUE_MAPQ,
            });
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert!(
            !jcs.is_empty(),
            "softclipped IS motif must place a JC onto a repeat copy; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let side2: u64 = jcs[0].fields[4].parse().unwrap();
        assert!(
            (41..=60).contains(&side2) || (101..=120).contains(&side2),
            "side2 should land in an IS copy, got {side2}"
        );
    }

    #[test]
    fn same_reads_placed_on_two_copies_fold_to_one_jc() {
        // Multi-copy placement folding (job 2372015). Regression base for
        // the synth_is_mob empty-GD failures (jobs 2371258, 2371412): plus
        // and minus reads with 3′ softclips cross the SAME insert left
        // boundary. SAM/BAM stores SEQ reference-forward, so both strands
        // share the same stored SEQ/CIGAR and must key side1 at the
        // junction base (30). The SAME reads' clips hit both motif copies,
        // so after folding exactly ONE JC remains (not one per copy); with
        // no copy-anchored reads the two placements tie on support and the
        // deterministic tie-break keeps the smallest side2 (first copy).
        let motif = *b"TGCAAGTCGATCGTTAGCCA";
        let mut seq = vec![b'C'; 160];
        seq[40..60].copy_from_slice(&motif);
        seq[100..120].copy_from_slice(&motif);
        // Junction: unique flank up to ref 30 (1-based), then motif content.
        let mut read_seq = vec![b'C'; 20];
        read_seq.extend_from_slice(&motif);
        let mut reads = Vec::new();
        for _ in 0..3 {
            for minus in [false, true] {
                reads.push(AlignedRead {
                    contig_idx: 0,
                    ref_start_0: 10, // 20M covers ref 11..30 (1-based)
                    minus,
                    seq: read_seq.clone(),
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
                    mapq: UNIQUE_MAPQ,
                });
            }
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            1,
            "same reads placed on two copies must fold to one JC; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let side1: u64 = jcs[0].fields[1].parse().unwrap();
        let side2: u64 = jcs[0].fields[4].parse().unwrap();
        assert_eq!(
            side1, 30,
            "both strands must key side1 at the junction base"
        );
        assert_eq!(
            side2, 41,
            "support tie keeps the lexicographically smallest side2"
        );
    }

    #[test]
    fn few_bp_strand_split_junctions_cluster_into_one_jc() {
        // Regression for the synth_is_mob job 2371412 zero-JC failure. IS
        // terminal / TSD register ambiguity keys reads of one insertion
        // junction a few bp apart (measured in the 2371412 BAM: plus M+S
        // reads end at 601/602, and breseq's own two JC for the event land
        // at 599/602, 3 bp apart). Plus reads: M+S ending at 599. Minus
        // reads: same reference-forward stored SEQ, M+S ending at 602.
        // Exact-key aggregation put the registers in different buckets;
        // ±5 bp clustering must merge them into exactly ONE accepted JC.
        let motif = *b"TGCAAGTCGATCGTTAGCCA";
        let mut seq = vec![b'C'; 700];
        seq[100..120].copy_from_slice(&motif); // one existing copy -> side2 = 101
                                               // Stored SEQ (reference-forward on both strands): match = unique
                                               // flank, 3′ clip = motif content.
        let mut read_seq = vec![b'C'; 20];
        read_seq.extend_from_slice(&motif);
        let mut reads = Vec::new();
        for _ in 0..3 {
            // Plus reads: match ref 580..599 (1-based) -> side1 = 599.
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 579,
                minus: false,
                seq: read_seq.clone(),
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
                mapq: UNIQUE_MAPQ,
            });
            // Minus reads: match ref 583..602 (1-based) -> side1 = 602.
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 582,
                minus: true,
                seq: read_seq.clone(),
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
                mapq: UNIQUE_MAPQ,
            });
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            1,
            "599/602 strand split must cluster into exactly one JC; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let side1: u64 = jcs[0].fields[1].parse().unwrap();
        let side2: u64 = jcs[0].fields[4].parse().unwrap();
        assert_eq!(
            side1, 599,
            "representative = first member after the (side1, side2) sort"
        );
        assert_eq!(side2, 101, "clip must place onto the motif copy start");
    }

    #[test]
    fn junctions_beyond_cluster_tolerance_do_not_merge() {
        // Two junctions 10 bp apart (> JC_CLUSTER_TOL_BP), each with both
        // strands, must stay two clusters -> two JC.
        let motif = *b"TGCAAGTCGATCGTTAGCCA";
        let mut seq = vec![b'C'; 700];
        seq[100..120].copy_from_slice(&motif);
        let mut read_seq = vec![b'C'; 20];
        read_seq.extend_from_slice(&motif);
        let mut reads = Vec::new();
        // Junction A: plus ends 599 / minus ends 602.
        // Junction B: plus ends 609 / minus ends 612 (10 bp from A).
        for (plus_start0, minus_start0) in [(579i64, 582i64), (589, 592)] {
            for _ in 0..3 {
                for (start0, minus) in [(plus_start0, false), (minus_start0, true)] {
                    reads.push(AlignedRead {
                        contig_idx: 0,
                        ref_start_0: start0,
                        minus,
                        seq: read_seq.clone(),
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
                        mapq: UNIQUE_MAPQ,
                    });
                }
            }
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            2,
            "junctions 10 bp apart must not merge; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let mut side1s: Vec<u64> = jcs.iter().map(|jc| jc.fields[1].parse().unwrap()).collect();
        side1s.sort_unstable();
        assert_eq!(side1s, vec![599, 609]);
    }

    /// 80 bp IS-like motif with 2 bp register coincidences against the
    /// flank markers used in the folding tests: the first 2 bp equal the
    /// right-flank start ("TT") and the last 2 bp equal the left-flank end
    /// ("AA"), mirroring the measured synth_is_mob motif/ref coincidences
    /// (job 2371412 diagnostic: motif prefix GT == ref 601-602, motif
    /// suffix GC == ref 599-600).
    fn is_motif_80() -> Vec<u8> {
        let mut m = Vec::new();
        m.extend_from_slice(b"TTGCAAGTCGATCGTTAGCC"); // [0..20]
        m.extend_from_slice(b"GATCCGATGCTAGCTTACGA"); // [20..40]
        m.extend_from_slice(b"ACGTGGATCCTTAGCGATCA"); // [40..60]
        m.extend_from_slice(b"GCTAGGCCATGATCGATAAA"); // [60..80]
        m
    }

    fn softclip_read(
        ref_start_0: i64,
        minus: bool,
        seq: &[u8],
        cigar: &[(CigarKind, usize)],
    ) -> AlignedRead {
        AlignedRead {
            contig_idx: 0,
            ref_start_0,
            minus,
            seq: seq.to_vec(),
            cigar: cigar
                .iter()
                .map(|&(kind, len)| CigarOp { kind, len })
                .collect(),
            mapq: UNIQUE_MAPQ,
        }
    }

    #[test]
    fn synth_is_mob_2372015_multicopy_and_reverse_duplicates_fold() {
        // End-to-end geometry of job 2372015 (4000 bp genome, 80 bp motif
        // copies at 1201-1280 / 2481-2560, third copy inserted at 601).
        // Each insert boundary is reported (a) once per existing motif copy
        // (multi-copy placement of the same clipped reads) and (b) once in
        // the reverse direction by reads aligned on a copy with the unique
        // flank clipped. Before folding this emitted 6 JC for 2 physical
        // junctions (4 redundant over-JC vs breseq: 599->1278, 602->1203,
        // 2480->599, 2560->601); after folding exactly 2 canonical JC must
        // remain.
        let motif = is_motif_80();
        // Unique flank markers (NOT homopolymers: copy-anchored reads clip
        // the flank, and place_softclip searches both strands — an A/T
        // homopolymer pair would reverse-complement onto each other).
        // lmark ends in "AA" and rmark starts with "TT", mirroring the
        // measured 2 bp register coincidences with the motif ends.
        let lmark = *b"ACGATCGGTTAACTGCGGAA"; // genome 581..600 (1-based)
        let rmark = *b"TTGCGATCGATCCGGAATCG"; // genome 601..620 (1-based)
        let mut seq = vec![b'C'; 4000];
        seq[580..600].copy_from_slice(&lmark);
        seq[600..620].copy_from_slice(&rmark);
        seq[1200..1280].copy_from_slice(&motif); // copy1, 1-based 1201..1280
        seq[2480..2560].copy_from_slice(&motif); // copy2, 1-based 2481..2560

        // Stored SEQ (reference-forward on both strands).
        let mut left_junction_seq = lmark.to_vec();
        left_junction_seq.extend_from_slice(&motif[0..20]);
        let mut right_junction_seq = Vec::new();
        right_junction_seq.extend_from_slice(&motif[60..80]);
        right_junction_seq.extend_from_slice(&rmark);

        let mut reads = Vec::new();
        // Left junction, flank-anchored: M on the left flank, trailing clip
        // = motif prefix. Register 600 (20M+20S) and register 602 (22M+18S;
        // the M block extends 2 bp into the motif because motif[0..2]=="TT"
        // equals the right-flank start).
        for minus in [false, false, true, true] {
            reads.push(softclip_read(
                580,
                minus,
                &left_junction_seq,
                &[(CigarKind::Match, 20), (CigarKind::SoftClip, 20)],
            ));
        }
        for minus in [false, true] {
            reads.push(softclip_read(
                580,
                minus,
                &left_junction_seq,
                &[(CigarKind::Match, 22), (CigarKind::SoftClip, 18)],
            ));
        }
        // Left junction, copy-anchored (bowtie2 placed the motif block on
        // copy2): leading clip = left-flank suffix.
        for minus in [false, true] {
            reads.push(softclip_read(
                2480,
                minus,
                &left_junction_seq,
                &[(CigarKind::SoftClip, 20), (CigarKind::Match, 20)],
            ));
        }
        // Right junction, flank-anchored: leading clip = motif suffix, M on
        // the right flank. Register 601 (20S+20M) and register 599
        // (18S+22M; motif[78..80]=="AA" equals the left-flank end).
        for minus in [false, false, true, true] {
            reads.push(softclip_read(
                600,
                minus,
                &right_junction_seq,
                &[(CigarKind::SoftClip, 20), (CigarKind::Match, 20)],
            ));
        }
        for minus in [false, true] {
            reads.push(softclip_read(
                598,
                minus,
                &right_junction_seq,
                &[(CigarKind::SoftClip, 18), (CigarKind::Match, 22)],
            ));
        }
        // Right junction, copy-anchored: M on copy2 end, trailing clip =
        // right-flank prefix.
        for minus in [false, true] {
            reads.push(softclip_read(
                2540,
                minus,
                &right_junction_seq,
                &[(CigarKind::Match, 20), (CigarKind::SoftClip, 20)],
            ));
        }

        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            2,
            "6 redundant reports must fold to exactly 2 canonical JC; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        // Canonical coordinates: the copy2 placements carry the extra
        // copy-anchored reads, so they win the multi-copy fold.
        assert_eq!(
            jcs[0].fields,
            vec!["chr", "599", "1", "chr", "2558", "-1", "0"],
            "left junction canonical JC"
        );
        assert_eq!(
            jcs[1].fields,
            vec!["chr", "600", "-1", "chr", "2481", "1", "0"],
            "right junction canonical JC"
        );
    }

    #[test]
    fn distinct_junctions_apart_beyond_tol_are_not_folded() {
        // Two genuinely distinct junctions 40 bp apart (> tol), each
        // supported by its own read set whose clips place on both motif
        // copies. Each junction's two copy placements fold (same reads);
        // the two junctions must NOT fold into each other.
        let motif = is_motif_80();
        let mut seq = vec![b'C'; 4000];
        seq[580..600].fill(b'A');
        seq[620..640].fill(b'G');
        seq[1200..1280].copy_from_slice(&motif);
        seq[2480..2560].copy_from_slice(&motif);
        let mut reads = Vec::new();
        for (flank_base, start0) in [(b'A', 580i64), (b'G', 620)] {
            let mut read_seq = vec![flank_base; 20];
            read_seq.extend_from_slice(&motif[0..20]);
            for minus in [false, false, false, true, true, true] {
                reads.push(softclip_read(
                    start0,
                    minus,
                    &read_seq,
                    &[(CigarKind::Match, 20), (CigarKind::SoftClip, 20)],
                ));
            }
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            2,
            "junctions 40 bp apart must not fold; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let mut side1s: Vec<u64> = jcs.iter().map(|jc| jc.fields[1].parse().unwrap()).collect();
        side1s.sort_unstable();
        assert_eq!(side1s, vec![600, 640]);
        for jc in &jcs {
            // Equal support on both copies -> smallest side2 (copy1).
            assert_eq!(jc.fields[4], "1201");
        }
    }

    #[test]
    fn distinct_read_sets_same_flank_are_not_folded() {
        // Same side1 (within tol) but DISJOINT supporting read sets placing
        // on two different copies: not the same underlying reads, so the
        // multi-copy rule must not fold them (identity gate).
        let motif_p = is_motif_80();
        let mut motif_q = Vec::new();
        motif_q.extend_from_slice(b"GCGTTAACCGGATCGTACGT");
        motif_q.extend_from_slice(b"CCGGAATTAACCGGTTAACC");
        motif_q.extend_from_slice(b"TTGGCCAATTGGCCAATTGG");
        motif_q.extend_from_slice(b"AACCGGTTCCTTAAGGCCAA");
        let mut seq = vec![b'C'; 4000];
        seq[580..600].fill(b'A');
        seq[1200..1280].copy_from_slice(&motif_p);
        seq[2480..2560].copy_from_slice(&motif_q);
        let mut reads = Vec::new();
        for motif in [&motif_p, &motif_q] {
            let mut read_seq = vec![b'A'; 20];
            read_seq.extend_from_slice(&motif[0..20]);
            for minus in [false, false, false, true, true, true] {
                reads.push(softclip_read(
                    580,
                    minus,
                    &read_seq,
                    &[(CigarKind::Match, 20), (CigarKind::SoftClip, 20)],
                ));
            }
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            2,
            "disjoint read sets must not fold even at the same side1; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let mut side2s: Vec<u64> = jcs.iter().map(|jc| jc.fields[4].parse().unwrap()).collect();
        side2s.sort_unstable();
        assert_eq!(side2s, vec![1201, 2481]);
    }

    #[test]
    fn cigar_split_junctions_are_never_multicopy_folded() {
        // Two long-D deletions sharing side1 (within tol) but reaching
        // different side2 positions: distinct events. CIGAR splits have
        // unique identity and must survive as 2 JC.
        let seq = vec![b'A'; 400];
        let mut reads = Vec::new();
        for del_len in [8usize, 100] {
            for minus in [false, false, false, true, true, true] {
                reads.push(softclip_read(
                    0,
                    minus,
                    &[b'A'; 30],
                    &[
                        (CigarKind::Match, 15),
                        (CigarKind::Del, del_len),
                        (CigarKind::Match, 15),
                    ],
                ));
            }
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            2,
            "CIGAR-split junctions must never be multi-copy folded; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let mut side2s: Vec<u64> = jcs.iter().map(|jc| jc.fields[4].parse().unwrap()).collect();
        side2s.sort_unstable();
        assert_eq!(side2s, vec![24, 116]);
    }

    #[test]
    fn jc_supported_cigar_gap_promotes_short_del() {
        let seq = vec![b'A'; 50];
        let mut reads = Vec::new();
        for minus in [false, false, false, true, true, true] {
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus,
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
                mapq: UNIQUE_MAPQ,
            });
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(
            !gd_dels(&gd).is_empty(),
            "short true gap with both-flank JC must become DEL; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
    }

    fn short_clip_is_seed_reads() -> (Vec<u8>, Vec<AlignedRead>) {
        let mut seq = vec![b'C'; 160];
        let motif = *b"ACGTACGTAC";
        seq[40..50].copy_from_slice(&motif);
        seq[100..110].copy_from_slice(&motif);
        let mut read_seq = vec![b'C'; 20];
        read_seq.extend_from_slice(&motif);
        let mut reads = Vec::new();
        for minus in [false, false, false, true, true, true] {
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 70,
                minus,
                seq: read_seq.clone(),
                cigar: vec![
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 20,
                    },
                    CigarOp {
                        kind: CigarKind::SoftClip,
                        len: 10,
                    },
                ],
                mapq: UNIQUE_MAPQ,
            });
        }
        (seq, reads)
    }

    #[test]
    fn ten_bp_clips_seed_but_do_not_accept_without_second_pass() {
        let (seq, reads) = short_clip_is_seed_reads();
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(
            gd_jcs(&gd).is_empty(),
            "10 bp clips are below the 14 bp accept floor; got {:?}",
            gd_jcs(&gd).iter().map(|e| &e.fields).collect::<Vec<_>>()
        );
    }

    #[test]
    fn second_pass_hits_promote_short_clip_seed_to_jc() {
        let (seq, reads) = short_clip_is_seed_reads();
        // Seed keys side1 at last match (90). Extra hits model spanning
        // second-pass alignments (18 bp into each side, both strands).
        let extra: Vec<SplitCandidate> = [false, false, false, true, true, true]
            .into_iter()
            .map(|minus| SplitCandidate {
                contig_idx: 0,
                side2_contig_idx: 0,
                minus,
                side1_minus: true,
                side2_minus: false,
                left: SubAlignment {
                    read_start: 0,
                    read_end: 18,
                },
                right: SubAlignment {
                    read_start: 18,
                    read_end: 36,
                },
                side1_pos_1: 90,
                side2_pos_1: 41,
                overlap: 0,
                origin: SplitOrigin::Softclip(10_000),
            })
            .collect();
        let gd = call_from_aligned_extra(&fasta_chr(&seq), &reads, &opts_single_thread(), &extra);
        assert!(
            !gd_jcs(&gd).is_empty(),
            "second-pass spanning hits must lift the seed past accept_junction; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
    }

    // --- candidate-junction second pass (`parse_junction_bam`) ---------------
    //
    // Fixture BAMs are authored as SAM text and converted with the production
    // `sam_to_sorted_bam`, so these cover the real decode path.

    fn write_fixture_bam(dir: &str, sq: &[(&str, usize)], recs: &[String]) -> PathBuf {
        let d = std::env::temp_dir().join(dir);
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        let mut text = String::from("@HD\tVN:1.6\n");
        for (name, len) in sq {
            text.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
        }
        for r in recs {
            text.push_str(r);
            text.push('\n');
        }
        let sam = d.join("in.sam");
        std::fs::write(&sam, text).unwrap();
        let bam = d.join("in.bam");
        crate::align::sam_to_sorted_bam(&sam, &bam).unwrap();
        bam
    }

    fn sam_rec(
        qname: &str,
        flag: u16,
        rname: &str,
        pos: usize,
        cigar: &str,
        qlen: usize,
    ) -> String {
        format!(
            "{qname}\t{flag}\t{rname}\t{pos}\t44\t{cigar}\t*\t0\t0\t{}\t{}",
            "A".repeat(qlen),
            "I".repeat(qlen)
        )
    }

    /// Candidate paired with the breakpoint of ITS construct.
    fn jc_cand(breakpoint: usize) -> (crate::jc_seq::CandidateJunction, usize) {
        (
            crate::jc_seq::CandidateJunction {
                side1_contig: 0,
                side1_pos_1: 90,
                side1_minus: true,
                side2_contig: 0,
                side2_pos_1: 41,
                side2_minus: false,
                overlap: 0,
                seed_reads: 3,
            },
            breakpoint,
        )
    }

    #[test]
    fn second_pass_uses_per_construct_breakpoint_not_flank() {
        // Construct seeded at a contig start: left flank clamped to 10 bp, so
        // the breakpoint is 10 even though `flank` was 36. A 30M hit at
        // construct pos 1 puts 10 bp left and 20 bp right of the true
        // breakpoint — genuine spanning evidence.
        let bam = write_fixture_bam(
            "prokdiff-jcbam-breakpoint",
            &[("jc0", 30)],
            &[sam_rec("r1", 0, "jc0", 1, "30M", 30)],
        );
        let scores = HashMap::new();

        let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(10)], &scores).unwrap();
        assert_eq!(stats.considered, 1);
        assert_eq!(stats.kept, 1);
        assert_eq!(extra.len(), 1);
        assert_eq!(
            extra[0].left.len(),
            10,
            "side 1 extension is breakpoint-cstart"
        );
        assert_eq!(extra[0].right.len(), 20);

        // Same record against the pre-fix constant breakpoint (`flank` = 36):
        // the whole hit looks left-of-breakpoint and the read is discarded.
        // This is the JC=0 mechanism for every construct near a contig edge.
        let (extra_flank, stats_flank) = parse_junction_bam(&bam, &[jc_cand(36)], &scores).unwrap();
        assert!(extra_flank.is_empty());
        assert_eq!(stats_flank.not_spanning, 1);
    }

    #[test]
    fn second_pass_rejects_two_bp_overhang_instead_of_clamping_it_up() {
        // 36M at construct pos 35 lands 2 bp past a breakpoint of 36. The
        // engine used to admit this and clamp both extensions to >=3, which
        // fabricated evidence and made `accept_junction` rule 4 unreachable.
        let bam = write_fixture_bam(
            "prokdiff-jcbam-overhang",
            &[("jc0", 72)],
            &[sam_rec("r1", 0, "jc0", 35, "36M", 36)],
        );
        let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
        assert!(extra.is_empty(), "2 bp past the breakpoint is not evidence");
        assert_eq!(stats.not_spanning, 1);
        assert_eq!(stats.kept, 0);
    }

    #[test]
    fn second_pass_reports_real_overlaps_with_no_floor() {
        // Exactly 3 bp on the short side: admitted, and reported as 3 — not
        // rounded, not clamped.
        let bam = write_fixture_bam(
            "prokdiff-jcbam-exact3",
            &[("jc0", 72)],
            &[sam_rec("r1", 0, "jc0", 34, "36M", 36)],
        );
        let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
        assert_eq!(stats.kept, 1);
        assert_eq!(extra[0].left.len(), 3);
        assert_eq!(extra[0].right.len(), 33);
    }

    #[test]
    fn second_pass_counts_short_cover_rejects() {
        // 20 aligned query bases is below JC_MIN_COVER_BASES (28).
        let bam = write_fixture_bam(
            "prokdiff-jcbam-shortcover",
            &[("jc0", 72)],
            &[sam_rec("r1", 0, "jc0", 25, "20M16S", 36)],
        );
        let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
        assert!(extra.is_empty());
        assert_eq!(stats.short_cover, 1);
        assert_eq!(stats.not_spanning, 0, "cover is checked before spanning");
    }

    #[test]
    fn second_pass_counts_primary_score_gate_separately_from_spanning() {
        // A spanning 36M hit scores 36. The gate keeps it against a weaker
        // primary alignment and drops it against a stronger one; the counter
        // distinguishes the two so a Clonal job can tell which gate ate the
        // reads without re-running.
        let bam = write_fixture_bam(
            "prokdiff-jcbam-scoregate",
            &[("jc0", 72)],
            &[sam_rec("r1", 0, "jc0", 19, "36M", 36)],
        );
        let cands = [jc_cand(36)];

        let weaker = HashMap::from([("r1".to_string(), 30)]);
        let (extra, stats) = parse_junction_bam(&bam, &cands, &weaker).unwrap();
        assert_eq!(stats.kept, 1);
        assert_eq!(stats.worse_than_primary, 0);
        assert_eq!(extra.len(), 1);

        let stronger = HashMap::from([("r1".to_string(), 40)]);
        let (extra, stats) = parse_junction_bam(&bam, &cands, &stronger).unwrap();
        assert!(extra.is_empty());
        assert_eq!(stats.worse_than_primary, 1);
        assert_eq!(stats.not_spanning, 0);
    }

    #[test]
    fn second_pass_origins_are_per_read_and_never_alias_softclip_ids() {
        let bam = write_fixture_bam(
            "prokdiff-jcbam-origins",
            &[("jc0", 72)],
            &[
                sam_rec("r1", 0, "jc0", 19, "36M", 36),
                sam_rec("r2", 16, "jc0", 20, "36M", 36),
            ],
        );
        let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
        assert_eq!(stats.kept, 2);
        // Two reads, two identities: the multi-copy fold must be able to tell
        // them apart.
        let ids: HashSet<SplitOrigin> = extra.iter().map(|c| c.origin).collect();
        assert_eq!(ids.len(), 2);
        // Variant-tagged, so a qname hash can never collide with a first-pass
        // `SoftclipHint` index in the supporting-read identity set.
        assert!(extra
            .iter()
            .all(|c| matches!(c.origin, SplitOrigin::SecondPass(_))));
        // Read strand is carried through; side strands come from the candidate.
        assert!(extra.iter().any(|c| c.minus));
        assert!(extra.iter().any(|c| !c.minus));
    }

    #[test]
    fn second_pass_ignores_records_on_unknown_constructs() {
        // A reference name that is not `jc<idx>`, and an index past the end of
        // the retained candidate list: neither may be counted or promoted.
        let bam = write_fixture_bam(
            "prokdiff-jcbam-unknown",
            &[("chr", 72), ("jc7", 72)],
            &[
                sam_rec("r1", 0, "chr", 19, "36M", 36),
                sam_rec("r2", 0, "jc7", 19, "36M", 36),
            ],
        );
        let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
        assert!(extra.is_empty());
        assert_eq!(stats.considered, 0);
    }
}
