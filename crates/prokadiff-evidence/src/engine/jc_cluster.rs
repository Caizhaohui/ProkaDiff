//! Junction clustering, redundancy folding, and candidate construction.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use crate::align::{align_to_bam, AlignKind, FastqInput};
use crate::engine::bam_io::{parse_junction_bam, PrimaryBamData};
use crate::engine::EngineOptions;
use crate::error::Result;
use crate::fasta::FastaRecord;
use crate::jc::JunctionSupport;
use crate::pileup::{SplitCandidate, SplitOrigin};

/// Two split candidates whose side1 AND side2 both land within this many bp
/// are clustered into one junction before `accept_junction`.
pub const JC_CLUSTER_TOL_BP: u64 = 5;

/// One junction cluster: representative coordinates are those of the FIRST
/// member. `JunctionSupport` is aggregated over ALL members.
pub(crate) struct JcCluster {
    pub(crate) rep_side1: u64,
    pub(crate) rep_side2: u64,
    pub(crate) rep_overlap: i64,
    pub(crate) rep_side1_minus: bool,
    pub(crate) rep_side2_minus: bool,
    pub(crate) support: JunctionSupport,
    pub(crate) ids: HashSet<(usize, SplitOrigin)>,
    pub(crate) has_cigar: bool,
}

impl JcCluster {
    pub(crate) fn new(sp: &SplitCandidate) -> Self {
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

    pub(crate) fn add(&mut self, sp: &SplitCandidate) {
        match sp.origin {
            SplitOrigin::Cigar => self.has_cigar = true,
            SplitOrigin::Softclip(_) | SplitOrigin::SecondPass(_) => {
                self.ids.insert((sp.contig_idx, sp.origin));
            }
        }
        let ov1 = sp.left.len();
        let ov2 = sp.right.len();
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
/// evidence, carried through redundancy folding.
pub(crate) struct AcceptedJc {
    pub(crate) c1: usize,
    pub(crate) p1: u64,
    pub(crate) m1: bool,
    pub(crate) c2: usize,
    pub(crate) p2: u64,
    pub(crate) m2: bool,
    pub(crate) overlap: i64,
    pub(crate) support: JunctionSupport,
    pub(crate) ids: HashSet<(usize, SplitOrigin)>,
    pub(crate) has_cigar: bool,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct JcSideRef<'a> {
    pub(crate) pos: u64,
    pub(crate) minus: bool,
    pub(crate) contig: &'a str,
}

impl AcceptedJc {
    pub(crate) fn side1<'a>(&self, fasta: &'a [FastaRecord]) -> JcSideRef<'a> {
        JcSideRef {
            pos: self.p1,
            minus: self.m1,
            contig: &fasta[self.c1].name,
        }
    }
    pub(crate) fn side2<'a>(&self, fasta: &'a [FastaRecord]) -> JcSideRef<'a> {
        JcSideRef {
            pos: self.p2,
            minus: self.m2,
            contig: &fasta[self.c2].name,
        }
    }
}

pub(crate) type SideKey = (usize, u64, bool);
pub(crate) type DirKey = (SideKey, SideKey);

pub(crate) fn dir_key(j: &AcceptedJc) -> DirKey {
    ((j.c1, j.p1, j.m1), (j.c2, j.p2, j.m2))
}

pub(crate) fn canonical_dir_key(j: &AcceptedJc) -> DirKey {
    let f = dir_key(j);
    let r = (f.1, f.0);
    f.min(r)
}

pub(crate) fn dir_key_within_tol(a: &DirKey, b: &DirKey) -> bool {
    a.0 .0 == b.0 .0
        && a.0 .2 == b.0 .2
        && a.0 .1.abs_diff(b.0 .1) <= JC_CLUSTER_TOL_BP
        && a.1 .0 == b.1 .0
        && a.1 .2 == b.1 .2
        && a.1 .1.abs_diff(b.1 .1) <= JC_CLUSTER_TOL_BP
}

pub(crate) fn absorb_support(a: &mut JunctionSupport, b: &JunctionSupport) {
    a.plus_reads += b.plus_reads;
    a.minus_reads += b.minus_reads;
    a.best_min_overlap = a.best_min_overlap.max(b.best_min_overlap);
    a.plus_best_min = a.plus_best_min.max(b.plus_best_min);
    a.minus_best_min = a.minus_best_min.max(b.minus_best_min);
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

pub(crate) fn support_total(s: &JunctionSupport) -> usize {
    s.plus_reads + s.minus_reads
}

pub(crate) fn absorb_into_rep(g: &mut AcceptedJc, j: AcceptedJc) {
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

pub(crate) fn fold_reverse_duplicates(mut accepted: Vec<AcceptedJc>) -> Vec<AcceptedJc> {
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

pub(crate) fn shares_one_side(a: &AcceptedJc, b: &AcceptedJc) -> bool {
    let sa = [(a.c1, a.p1, a.m1), (a.c2, a.p2, a.m2)];
    let sb = [(b.c1, b.p1, b.m1), (b.c2, b.p2, b.m2)];
    sa.iter().any(|x| {
        sb.iter()
            .any(|y| x.0 == y.0 && x.1.abs_diff(y.1) <= 2 && x.2 == y.2)
    })
}

pub(crate) fn find_repeat_at<'a>(
    repeats: &'a [crate::fasta::RepeatRegion],
    contig: &str,
    pos: u64,
) -> Option<&'a crate::fasta::RepeatRegion> {
    repeats
        .iter()
        .find(|r| r.seq_id == contig && (pos.abs_diff(r.start) <= 20 || pos.abs_diff(r.end) <= 20))
}

pub(crate) fn fold_multicopy_placements(
    mut junctions: Vec<AcceptedJc>,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) -> Vec<AcceptedJc> {
    junctions.sort_by_key(canonical_dir_key);
    let mut groups: Vec<AcceptedJc> = Vec::new();
    for j in junctions {
        let mut hit = None;
        for (i, g) in groups.iter().enumerate() {
            if g.has_cigar && j.has_cigar {
                continue;
            }
            if !shares_one_side(g, &j) {
                continue;
            }
            let shared = g.ids.intersection(&j.ids).count();
            let smaller = g.ids.len().min(j.ids.len());
            let shared_reads = smaller > 0 && shared * 2 > smaller;

            let same_repeat_family = if !repeats.is_empty() {
                let mut shared_unique = None;
                if g.c1 == j.c1 && g.p1.abs_diff(j.p1) <= 2 && g.m1 == j.m1 {
                    shared_unique = Some(((g.c2, g.p2), (j.c2, j.p2)));
                } else if g.c1 == j.c2 && g.p1.abs_diff(j.p2) <= 2 && g.m1 == j.m2 {
                    shared_unique = Some(((g.c2, g.p2), (j.c1, j.p1)));
                } else if g.c2 == j.c1 && g.p2.abs_diff(j.p1) <= 2 && g.m2 == j.m1 {
                    shared_unique = Some(((g.c1, g.p1), (j.c2, j.p2)));
                } else if g.c2 == j.c2 && g.p2.abs_diff(j.p2) <= 2 && g.m2 == j.m2 {
                    shared_unique = Some(((g.c1, g.p1), (j.c1, j.p1)));
                }

                if let Some(((g_rep_c, g_rep_p), (j_rep_c, j_rep_p))) = shared_unique {
                    let g_rep = find_repeat_at(repeats, &fasta[g_rep_c].name, g_rep_p);
                    let j_rep = find_repeat_at(repeats, &fasta[j_rep_c].name, j_rep_p);
                    match (g_rep, j_rep) {
                        (Some(r1), Some(r2)) => r1.name == r2.name,
                        _ => false,
                    }
                } else {
                    false
                }
            } else {
                false
            };

            if shared_reads || same_repeat_family {
                hit = Some(i);
                break;
            }
        }
        match hit {
            Some(i) => absorb_into_rep(&mut groups[i], j),
            None => groups.push(j),
        }
    }
    groups
}

pub(crate) fn cluster_key_splits(
    splits: &[SplitCandidate],
) -> Vec<((usize, usize, i8), JcCluster)> {
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

pub(crate) fn seed_candidates_from_pileup(
    contig_results: &[crate::engine::ContigPileup],
) -> Vec<crate::jc_seq::CandidateJunction> {
    let mut cands = Vec::new();
    for cp in contig_results {
        for (gk, c) in cluster_key_splits(&cp.splits) {
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

pub(crate) fn second_pass_splits(
    fasta: &[FastaRecord],
    reads: &FastqInput,
    primary_data: &PrimaryBamData,
    contig_results: &[crate::engine::ContigPileup],
    opts: &EngineOptions,
    work: &Path,
) -> Result<Vec<SplitCandidate>> {
    use crate::fasta::write_records;
    use crate::jc_seq::{junction_construct, rank_and_cap, CandidateJunction, JC_MIN_COVER_BASES};

    let flank = primary_data
        .aligned
        .iter()
        .map(|r| r.seq.len())
        .max()
        .unwrap_or(36)
        .max(28);
    let mut all_cands = seed_candidates_from_pileup(contig_results);

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
        let pass_dir = jc_dir.join("reads");
        let (filtered, nkeep) = crate::align::filter_fastq_input(
            reads,
            &primary_data.second_pass_keep,
            &primary_data.second_pass_seen,
            &pass_dir,
        )?;
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
    let (extra, stats) = parse_junction_bam(&jc_bam, &kept, &primary_data.primary_scores)?;
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

pub(crate) fn normalize_jc_coords(
    s1: JcSideRef<'_>,
    s2: JcSideRef<'_>,
    overlap: i64,
    repeats: &[crate::fasta::RepeatRegion],
) -> (u64, u64) {
    let mut np1 = s1.pos;
    let mut np2 = s2.pos;
    if overlap > 0 && overlap <= 4 {
        let k = overlap as u64;
        let rep1 = find_repeat_at(repeats, s1.contig, s1.pos);
        let rep2 = find_repeat_at(repeats, s2.contig, s2.pos);
        match (rep1, rep2) {
            (Some(_), None) => {
                if s2.minus {
                    np2 = np2.saturating_sub(k);
                } else {
                    np2 = np2.saturating_add(k);
                }
            }
            (None, Some(_)) => {
                if s1.minus {
                    np1 = np1.saturating_sub(k);
                } else {
                    np1 = np1.saturating_add(k);
                }
            }
            _ => {
                if s1.minus {
                    np1 = np1.saturating_sub(k);
                } else if !s2.minus {
                    np2 = np2.saturating_add(k);
                } else {
                    np1 = np1.saturating_add(k);
                }
            }
        }
    }
    (np1, np2)
}
