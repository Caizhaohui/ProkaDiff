//! Junction clustering, redundancy folding, and candidate construction.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use tracing::{debug, info};

use crate::align::{align_to_bam, AlignKind, FastqInput};
use crate::engine::bam_io::{parse_junction_bam, PrimaryBamData};
use crate::engine::EngineOptions;
use crate::error::Result;
use crate::fasta::FastaRecord;
use crate::jc::JunctionSupport;
use crate::pileup::{SplitCandidate, SplitOrigin};
use crate::repeat_ambiguous::RepeatAmbiguousSeed;

/// Two split candidates whose side1 AND side2 both land within this many bp
/// are clustered into one junction before `accept_junction`.
pub const JC_CLUSTER_TOL_BP: u64 = 5;

/// Identical element ends share reads. A spanning-count gap no larger than
/// this is alignment enumeration, not copy-discriminating evidence.
pub(crate) const COPY_SUPPORT_TIE_GAP: usize = 2;

/// Repeat region fully containing `pos` (any interior position). Unlike
/// [`find_repeat_at`], this is not limited to copy termini.
pub(crate) fn repeat_containing<'a>(
    repeats: &'a [crate::fasta::RepeatRegion],
    contig: &str,
    pos: u64,
) -> Option<&'a crate::fasta::RepeatRegion> {
    repeats
        .iter()
        .find(|r| r.seq_id == contig && pos >= r.start && pos <= r.end)
}

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
    /// Strongest Stage-2 seed support among members (see
    /// `SplitCandidate::seed_reads`).
    pub(crate) seed_reads: usize,
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
            seed_reads: 0,
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
        self.seed_reads = self.seed_reads.max(sp.seed_reads);
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

/// One observed placement of a junction. Scores are frozen when the
/// placement is accepted and are never rewritten by family aggregation.
#[derive(Clone, Debug)]
pub(crate) struct JcPlacement {
    pub(crate) c1: usize,
    pub(crate) p1: u64,
    pub(crate) m1: bool,
    pub(crate) c2: usize,
    pub(crate) p2: u64,
    pub(crate) m2: bool,
    pub(crate) overlap: i64,
    pub(crate) own_support_reads: usize,
    pub(crate) own_seed_reads: usize,
}

/// An accepted junction cluster plus the identity of its supporting
/// evidence, carried through redundancy folding.
///
/// `support` / `seed_reads` are the deduplicated family totals.
/// `own_support_reads` / `own_seed_reads` belong to the displayed
/// placement and are the only scores used to rank representatives.
#[derive(Clone)]
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
    /// Strongest Stage-2 seed support merged into this family (see
    /// `SplitCandidate::seed_reads`). Not a ranking key.
    pub(crate) seed_reads: usize,
    /// Spanning-read count of the displayed placement, frozen at accept.
    pub(crate) own_support_reads: usize,
    /// Stage-2 seed count of the displayed placement, frozen at accept.
    pub(crate) own_seed_reads: usize,
    /// Demonstrated-equivalent placements, including the displayed one.
    pub(crate) alternatives: Vec<JcPlacement>,
}

#[derive(Default)]
pub(crate) struct SecondPassResult {
    pub(crate) splits: Vec<SplitCandidate>,
    pub(crate) repeat_seeds: Vec<RepeatAmbiguousSeed>,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct JcSideRef<'a> {
    pub(crate) pos: u64,
    pub(crate) minus: bool,
    pub(crate) contig: &'a str,
}

impl AcceptedJc {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_parts(
        c1: usize,
        p1: u64,
        m1: bool,
        c2: usize,
        p2: u64,
        m2: bool,
        overlap: i64,
        support: JunctionSupport,
        ids: HashSet<(usize, SplitOrigin)>,
        has_cigar: bool,
        seed_reads: usize,
    ) -> Self {
        let own_support_reads = support_total(&support);
        let own_seed_reads = seed_reads;
        let mut j = Self {
            c1,
            p1,
            m1,
            c2,
            p2,
            m2,
            overlap,
            support,
            ids,
            has_cigar,
            seed_reads,
            own_support_reads,
            own_seed_reads,
            alternatives: Vec::new(),
        };
        j.alternatives.push(j.placement_snapshot());
        j
    }

    fn placement_snapshot(&self) -> JcPlacement {
        JcPlacement {
            c1: self.c1,
            p1: self.p1,
            m1: self.m1,
            c2: self.c2,
            p2: self.p2,
            m2: self.m2,
            overlap: self.overlap,
            own_support_reads: self.own_support_reads,
            own_seed_reads: self.own_seed_reads,
        }
    }

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

fn own_rank(j: &AcceptedJc) -> (usize, usize) {
    (j.own_support_reads, j.own_seed_reads)
}

/// Circular distance on one contig. `len == 0` is ordinary linear distance.
pub(crate) fn circular_bp_distance(a: u64, b: u64, len: u64) -> u64 {
    let d = a.abs_diff(b);
    if len == 0 {
        d
    } else {
        d.min(len.saturating_sub(d))
    }
}

/// Distance from a unique flank to the repeat-terminus side, when the
/// junction has exactly one of each. `None` when the geometry is not a
/// unique↔terminus pair (both sides repeat, both unique, or cross-contig).
pub(crate) fn unique_repeat_terminus_distance(
    j: &AcceptedJc,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) -> Option<u64> {
    if j.c1 != j.c2 || j.c1 >= fasta.len() {
        return None;
    }
    let name = fasta[j.c1].name.as_str();
    let term = |pos: u64| find_repeat_at(repeats, name, pos).is_some();
    let inside = |pos: u64| repeat_containing(repeats, name, pos).is_some();
    let unique = |pos: u64| !inside(pos) && !term(pos);
    let (t1, t2) = (term(j.p1), term(j.p2));
    let (u1, u2) = (unique(j.p1), unique(j.p2));
    if (t1 && u2) || (t2 && u1) {
        Some(circular_bp_distance(
            j.p1,
            j.p2,
            fasta[j.c1].seq.len() as u64,
        ))
    } else {
        None
    }
}

fn position_in_repeat(
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
    contig: usize,
    pos: u64,
) -> bool {
    let Some(rec) = fasta.get(contig) else {
        return false;
    };
    find_repeat_at(repeats, &rec.name, pos).is_some()
        || repeat_containing(repeats, &rec.name, pos).is_some()
}

/// The side these two junctions share is a unique flank, not a repeat.
fn shares_unique_flank(
    a: &AcceptedJc,
    b: &AcceptedJc,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) -> bool {
    let a_sides = [(a.c1, a.p1, a.m1), (a.c2, a.p2, a.m2)];
    let b_sides = [(b.c1, b.p1, b.m1), (b.c2, b.p2, b.m2)];
    a_sides.iter().any(|&side| {
        b_sides.iter().any(|&other| {
            side_matches(side, other) && !position_in_repeat(repeats, fasta, side.0, side.1)
        })
    })
}

/// Repeat side sits 1–5 bp outside the terminus that MOB promotion shifts.
fn outside_terminus_shortfall(
    j: &AcceptedJc,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) -> bool {
    for (contig, pos) in [(j.c1, j.p1), (j.c2, j.p2)] {
        let Some(rec) = fasta.get(contig) else {
            continue;
        };
        let Some(rep) = find_repeat_at(repeats, &rec.name, pos) else {
            continue;
        };
        let outside = if rep.strand >= 0 {
            pos < rep.start && rep.start.saturating_sub(pos) <= 5
        } else {
            pos > rep.end && pos.saturating_sub(rep.end) <= 5
        };
        if outside {
            return true;
        }
    }
    false
}

/// Display `challenger` instead of `existing` for an already-equivalent pair.
///
/// A support gap above [`COPY_SUPPORT_TIE_GAP`] wins outright. Inside that
/// window the counts are not copy-discriminating: a 1–5 bp pre-terminus
/// shortfall stays displayed so the MOB duplication is measured from that
/// register, and the remaining unique-flank tie displays the nearer
/// annotated end. Stage-2 seed does not move an equal-support flank onto a
/// farther identical terminus.
fn challenger_replaces_displayed(
    challenger: &AcceptedJc,
    existing: &AcceptedJc,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) -> bool {
    let max_own = challenger.own_support_reads.max(existing.own_support_reads);
    let challenger_eligible = challenger
        .own_support_reads
        .saturating_add(COPY_SUPPORT_TIE_GAP)
        >= max_own;
    let existing_eligible = existing
        .own_support_reads
        .saturating_add(COPY_SUPPORT_TIE_GAP)
        >= max_own;
    if challenger_eligible != existing_eligible {
        return challenger_eligible;
    }
    let challenger_short = outside_terminus_shortfall(challenger, repeats, fasta);
    let existing_short = outside_terminus_shortfall(existing, repeats, fasta);
    if challenger_short != existing_short {
        return challenger_short;
    }
    match (
        unique_repeat_terminus_distance(challenger, repeats, fasta),
        unique_repeat_terminus_distance(existing, repeats, fasta),
    ) {
        (Some(cd), Some(ed)) if cd != ed => return cd < ed,
        _ => {}
    }
    if challenger.own_support_reads != existing.own_support_reads {
        return challenger.own_support_reads > existing.own_support_reads;
    }
    challenger.own_seed_reads > existing.own_seed_reads
}

fn absorb_support_quality(a: &mut JunctionSupport, b: &JunctionSupport) {
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

fn placement_key(p: &JcPlacement) -> (usize, u64, bool, usize, u64, bool, i64) {
    (p.c1, p.p1, p.m1, p.c2, p.p2, p.m2, p.overlap)
}

fn push_alternative(alts: &mut Vec<JcPlacement>, alt: JcPlacement) {
    if let Some(existing) = alts
        .iter_mut()
        .find(|a| placement_key(a) == placement_key(&alt))
    {
        if (alt.own_support_reads, alt.own_seed_reads)
            > (existing.own_support_reads, existing.own_seed_reads)
        {
            *existing = alt;
        }
    } else {
        alts.push(alt);
    }
}

/// Add `src` into the family of `dst` without double-counting observations.
/// Displayed coordinates and `own_*` scores of `dst` stay put.
fn merge_family(dst: &mut AcceptedJc, src: AcceptedJc) {
    let fresh_n = src.ids.difference(&dst.ids).count();
    let src_ids_n = src.ids.len();
    if src_ids_n == 0 || fresh_n == src_ids_n {
        absorb_support(&mut dst.support, &src.support);
    } else if fresh_n == 0 {
        absorb_support_quality(&mut dst.support, &src.support);
    } else {
        dst.support.plus_reads += src.support.plus_reads * fresh_n / src_ids_n;
        dst.support.minus_reads += src.support.minus_reads * fresh_n / src_ids_n;
        absorb_support_quality(&mut dst.support, &src.support);
    }
    dst.seed_reads = dst.seed_reads.max(src.seed_reads);
    dst.ids.extend(src.ids);
    dst.has_cigar |= src.has_cigar;
    for alt in src.alternatives {
        push_alternative(&mut dst.alternatives, alt);
    }
}

fn refresh_displayed_own(g: &mut AcceptedJc) {
    g.own_support_reads = support_total(&g.support);
    g.own_seed_reads = g.seed_reads;
    let displayed = (g.c1, g.p1, g.m1, g.c2, g.p2, g.m2, g.overlap);
    if let Some(alt) = g
        .alternatives
        .iter_mut()
        .find(|alt| (alt.c1, alt.p1, alt.m1, alt.c2, alt.p2, alt.m2, alt.overlap) == displayed)
    {
        alt.own_support_reads = g.own_support_reads;
        alt.own_seed_reads = g.own_seed_reads;
    }
}

pub(crate) fn absorb_into_rep(g: &mut AcceptedJc, mut j: AcceptedJc) {
    // Rank the displayed placement by its own frozen scores. Family
    // totals are applied only after the winner is chosen, so absorbing
    // a support-10 sibling cannot bury a later placement whose own
    // support is 15.
    if own_rank(&j) > own_rank(g) {
        std::mem::swap(g, &mut j);
    }
    merge_family(g, j);
}

/// Merge two reports of one coordinate (a reverse duplicate). Different
/// reads of that same placement become its own score. Cross-copy folding
/// does not use this: there the frozen per-copy scores stay separate.
fn absorb_same_coordinate(g: &mut AcceptedJc, j: AcceptedJc) {
    absorb_into_rep(g, j);
    refresh_displayed_own(g);
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
            Some(i) => absorb_same_coordinate(&mut folded[i], j),
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
            .any(|y| x.0 == y.0 && x.1.abs_diff(y.1) <= JC_CLUSTER_TOL_BP && x.2 == y.2)
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

fn overlap_compatible(a: &AcceptedJc, b: &AcceptedJc) -> bool {
    a.overlap.signum() == b.overlap.signum()
        && (a.overlap - b.overlap).unsigned_abs() <= JC_CLUSTER_TOL_BP
}

fn shared_evidence(a: &AcceptedJc, b: &AcceptedJc) -> bool {
    let shared = a.ids.intersection(&b.ids).count();
    let smaller = a.ids.len().min(b.ids.len());
    smaller > 0 && shared * 2 > smaller
}

fn side_matches(a: (usize, u64, bool), b: (usize, u64, bool)) -> bool {
    a.0 == b.0 && a.2 == b.2 && a.1.abs_diff(b.1) <= JC_CLUSTER_TOL_BP
}

/// The side of `j` that is not the flank shared with `other`.
fn unmatched_side(j: &AcceptedJc, other: &AcceptedJc) -> Option<(usize, u64, bool)> {
    let sides = [(j.c1, j.p1, j.m1), (j.c2, j.p2, j.m2)];
    let other_sides = [
        (other.c1, other.p1, other.m1),
        (other.c2, other.p2, other.m2),
    ];
    let unmatched: Vec<_> = sides
        .into_iter()
        .filter(|side| !other_sides.iter().any(|other| side_matches(*side, *other)))
        .collect();
    if unmatched.len() == 1 {
        Some(unmatched[0])
    } else {
        None
    }
}

/// Element end a repeat-side coordinate is entering: `5` or `3`.
///
/// Forward copy: 5' is `start` on `+`, 3' is `end` on `-`.
/// Reverse copy: 5' is `end` on `-`, 3' is `start` on `+`.
pub(crate) fn element_end_kind(
    rep: &crate::fasta::RepeatRegion,
    pos: u64,
    minus: bool,
) -> Option<u8> {
    let d_start = pos.abs_diff(rep.start);
    let d_end = pos.abs_diff(rep.end);
    let at_start = d_start <= JC_CLUSTER_TOL_BP;
    let at_end = d_end <= JC_CLUSTER_TOL_BP;
    let nearer_start = at_start && (!at_end || d_start <= d_end);
    let nearer_end = at_end && !nearer_start;
    if rep.strand >= 0 {
        if nearer_start && !minus {
            Some(5)
        } else if nearer_end && minus {
            Some(3)
        } else {
            None
        }
    } else if nearer_end && minus {
        Some(5)
    } else if nearer_start && !minus {
        Some(3)
    } else {
        None
    }
}

pub(crate) fn place_element_end(rep: &crate::fasta::RepeatRegion, kind: u8) -> (u64, bool) {
    let five = kind == 5;
    if rep.strand >= 0 {
        if five {
            (rep.start, false)
        } else {
            (rep.end, true)
        }
    } else if five {
        (rep.end, true)
    } else {
        (rep.start, false)
    }
}

struct ElementLocus {
    name: String,
    offset: i64,
    end_kind: Option<u8>,
}

fn element_locus(
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
    contig: usize,
    pos: u64,
    minus: bool,
) -> Option<ElementLocus> {
    let rec = fasta.get(contig)?;
    let rep = find_repeat_at(repeats, &rec.name, pos)
        .or_else(|| repeat_containing(repeats, &rec.name, pos))?;
    let offset = if rep.strand >= 0 {
        pos as i64 - rep.start as i64
    } else {
        rep.end as i64 - pos as i64
    };
    Some(ElementLocus {
        name: rep.name.clone(),
        offset,
        end_kind: element_end_kind(rep, pos, minus),
    })
}

fn loci_compatible(a: &ElementLocus, b: &ElementLocus) -> bool {
    if a.name != b.name {
        return false;
    }
    match (a.end_kind, b.end_kind) {
        (Some(left), Some(right)) => left == right,
        (None, None) => (a.offset - b.offset).unsigned_abs() <= JC_CLUSTER_TOL_BP,
        _ => false,
    }
}

fn ends_match(
    a: &AcceptedJc,
    b: &AcceptedJc,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) -> bool {
    match (unmatched_side(a, b), unmatched_side(b, a)) {
        (Some((ca, pa, ma)), Some((cb, pb, mb))) => {
            match (
                element_locus(repeats, fasta, ca, pa, ma),
                element_locus(repeats, fasta, cb, pb, mb),
            ) {
                (Some(left), Some(right)) => loci_compatible(&left, &right),
                (None, None) => true,
                _ => false,
            }
        }
        (None, None) => true,
        _ => false,
    }
}

/// Same unique flank placed on one element end, with reads partitioned by
/// copy. Different unique positions that only share a repeat side do not
/// qualify: the shared side has to be the unique flank, and both unmatched
/// sides have to be that same annotated end.
fn partitioned_unique_flank_copies(
    a: &AcceptedJc,
    b: &AcceptedJc,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) -> bool {
    let (Some((ca, pa, ma)), Some((cb, pb, mb))) = (unmatched_side(a, b), unmatched_side(b, a))
    else {
        return false;
    };
    let (Some(left), Some(right)) = (
        element_locus(repeats, fasta, ca, pa, ma),
        element_locus(repeats, fasta, cb, pb, mb),
    ) else {
        return false;
    };
    left.end_kind.is_some()
        && loci_compatible(&left, &right)
        && shares_unique_flank(a, b, repeats, fasta)
}

/// Two junctions are the same event only when the overlap register and the
/// oriented element end agree, and either the reads overlap or the pair is
/// one unique flank split across identical ends. The repeat name by itself
/// is not an equivalence.
fn demonstrated_equivalent(
    a: &AcceptedJc,
    b: &AcceptedJc,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) -> bool {
    if a.has_cigar && b.has_cigar {
        return false;
    }
    if !shares_one_side(a, b) || !overlap_compatible(a, b) || !ends_match(a, b, repeats, fasta) {
        return false;
    }
    shared_evidence(a, b) || partitioned_unique_flank_copies(a, b, repeats, fasta)
}

/// Among alternatives inside the support tie window, keep a pre-terminus
/// shortfall, otherwise the nearer unique↔repeat terminus.
pub(crate) fn refresh_displayed_nearest_supported_copy(
    j: &mut AcceptedJc,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) {
    if unique_repeat_terminus_distance(j, repeats, fasta).is_none() {
        return;
    }
    let max_own = j
        .alternatives
        .iter()
        .map(|a| a.own_support_reads)
        .chain(std::iter::once(j.own_support_reads))
        .max()
        .unwrap_or(j.own_support_reads);
    let displayed_key = copy_display_key(j, repeats, fasta);
    let mut pick: Option<JcPlacement> = None;
    let mut best_key = displayed_key;
    for alt in &j.alternatives {
        if alt.own_support_reads.saturating_add(COPY_SUPPORT_TIE_GAP) < max_own {
            continue;
        }
        let mut trial = j.clone();
        trial.c1 = alt.c1;
        trial.p1 = alt.p1;
        trial.m1 = alt.m1;
        trial.c2 = alt.c2;
        trial.p2 = alt.p2;
        trial.m2 = alt.m2;
        trial.overlap = alt.overlap;
        let key = copy_display_key(&trial, repeats, fasta);
        if key > best_key {
            best_key = key;
            pick = Some(alt.clone());
        }
    }
    if let Some(alt) = pick {
        j.c1 = alt.c1;
        j.p1 = alt.p1;
        j.m1 = alt.m1;
        j.c2 = alt.c2;
        j.p2 = alt.p2;
        j.m2 = alt.m2;
        j.overlap = alt.overlap;
        j.own_support_reads = alt.own_support_reads;
        j.own_seed_reads = alt.own_seed_reads;
    }
}

/// Higher is better: a shortfall outranks an exact terminus, then the
/// nearer circular distance. `None` distance sorts last.
fn copy_display_key(
    j: &AcceptedJc,
    repeats: &[crate::fasta::RepeatRegion],
    fasta: &[FastaRecord],
) -> (u8, i64) {
    let short = u8::from(outside_terminus_shortfall(j, repeats, fasta));
    let near = unique_repeat_terminus_distance(j, repeats, fasta)
        .map(|d| -(d as i64))
        .unwrap_or(i64::MIN);
    (short, near)
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
            if demonstrated_equivalent(g, &j, repeats, fasta) {
                hit = Some(i);
                break;
            }
        }
        match hit {
            Some(i) => {
                // Equivalence is already established. Support, then a
                // pre-terminus shortfall, then distance to the unique flank
                // choose the displayed copy. Seed is only the last tie-break.
                if challenger_replaces_displayed(&j, &groups[i], repeats, fasta) {
                    let existing = std::mem::replace(&mut groups[i], j);
                    merge_family(&mut groups[i], existing);
                } else {
                    merge_family(&mut groups[i], j);
                }
            }
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

struct UniqueRepeatClass {
    unique: (usize, u64, bool),
    repeat_name: String,
    /// Bounds of the specific repeat COPY the repeat side joins. Two
    /// registers of one unique flank may join different copies of the
    /// same family; those are distinct constructs and must not be
    /// clustered together (Wave 9, Target B fallback).
    repeat_span: (u64, u64),
    end_kind: u8,
    dist: u64,
    repeat_pos: u64,
}

/// One unique flank joined to one annotated element end. Both-repeat
/// junctions and interiors return `None` so they are not collapsed.
fn classify_unique_repeat_candidate(
    cand: &crate::jc_seq::CandidateJunction,
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> Option<UniqueRepeatClass> {
    let sides = [
        (cand.side1_contig, cand.side1_pos_1, cand.side1_minus),
        (cand.side2_contig, cand.side2_pos_1, cand.side2_minus),
    ];
    let mut unique = None;
    let mut repeat_side = None;
    for (contig, pos, minus) in sides {
        let name = fasta.get(contig)?.name.as_str();
        if let Some(rep) = find_repeat_at(repeats, name, pos) {
            if repeat_side.is_some() {
                return None;
            }
            repeat_side = Some((contig, pos, minus, rep));
        } else if repeat_containing(repeats, name, pos).is_some() || unique.is_some() {
            return None;
        } else {
            unique = Some((contig, pos, minus));
        }
    }
    let (u_contig, u_pos, u_minus) = unique?;
    let (r_contig, r_pos, r_minus, rep) = repeat_side?;
    if u_contig != r_contig {
        return None;
    }
    let end_kind = element_end_kind(rep, r_pos, r_minus)?;
    let len = fasta.get(u_contig)?.seq.len() as u64;
    Some(UniqueRepeatClass {
        unique: (u_contig, u_pos, u_minus),
        repeat_name: rep.name.clone(),
        repeat_span: (rep.start, rep.end),
        end_kind,
        dist: circular_bp_distance(u_pos, r_pos, len),
        repeat_pos: r_pos,
    })
}

fn nearer_unique_repeat_representative(
    group: &[crate::jc_seq::CandidateJunction],
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> Option<usize> {
    if group.len() < 2 {
        return None;
    }
    let classes = group
        .iter()
        .map(|cand| classify_unique_repeat_candidate(cand, fasta, repeats))
        .collect::<Option<Vec<_>>>()?;
    let head = &classes[0];
    let same_flank = classes.iter().all(|class| {
        class.repeat_name == head.repeat_name
            && class.end_kind == head.end_kind
            && class.unique.0 == head.unique.0
            && class.unique.2 == head.unique.2
            && class.unique.1.abs_diff(head.unique.1) <= JC_CLUSTER_TOL_BP
    });
    if !same_flank {
        return None;
    }
    let mut best = 0usize;
    for (i, class) in classes.iter().enumerate().skip(1) {
        let closer = class.dist < classes[best].dist
            || (class.dist == classes[best].dist && class.repeat_pos < classes[best].repeat_pos);
        if closer {
            best = i;
        }
    }
    Some(best)
}

/// Whether any side of `cand` sits at (or inside) an annotated repeat.
fn candidate_touches_repeat(
    cand: &crate::jc_seq::CandidateJunction,
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> bool {
    [
        (cand.side1_contig, cand.side1_pos_1),
        (cand.side2_contig, cand.side2_pos_1),
    ]
    .into_iter()
    .any(|(contig, pos)| {
        fasta
            .get(contig)
            .map(|rec| rec.name.as_str())
            .is_some_and(|name| {
                find_repeat_at(repeats, name, pos).is_some()
                    || repeat_containing(repeats, name, pos).is_some()
            })
    })
}

/// The side of `cand` that is NOT at or inside a repeat, when exactly one
/// side qualifies. Used to check that unclassifiable members of an
/// identical-sequence group still describe the same unique flank.
fn unique_side_outside_repeats(
    cand: &crate::jc_seq::CandidateJunction,
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> Option<(usize, u64, bool)> {
    let sides = [
        (cand.side1_contig, cand.side1_pos_1, cand.side1_minus),
        (cand.side2_contig, cand.side2_pos_1, cand.side2_minus),
    ];
    let mut unique = None;
    for (contig, pos, minus) in sides {
        let name = fasta.get(contig)?.name.as_str();
        if find_repeat_at(repeats, name, pos).is_some()
            || repeat_containing(repeats, name, pos).is_some()
        {
            if unique.is_some() {
                return None;
            }
            continue;
        }
        unique = Some((contig, pos, minus));
    }
    unique
}

/// An identical-sequence group describes ONE junction when every member
/// shares one unique flank (within the cluster tolerance) or every member
/// is repeat-anchored on both sides. Groups pairing two distinct, far
/// apart unique flanks are distinct junctions and stay intact.
fn identical_group_has_one_flank(
    group: &[crate::jc_seq::CandidateJunction],
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> bool {
    let classes: Vec<Option<UniqueRepeatClass>> = group
        .iter()
        .map(|cand| classify_unique_repeat_candidate(cand, fasta, repeats))
        .collect();
    let classifiable: Vec<&UniqueRepeatClass> = classes.iter().flatten().collect();
    if classifiable.is_empty() {
        // Repeat↔repeat (or repeat-interior) groups collapse: identical
        // copies of one element end are indistinguishable to the aligner.
        // Pure unique↔unique identical groups are separate junctions.
        return group
            .iter()
            .all(|cand| candidate_touches_repeat(cand, fasta, repeats));
    }
    let head = classifiable[0];
    let same_flank = classifiable.iter().all(|class| {
        class.unique.0 == head.unique.0
            && class.unique.2 == head.unique.2
            && class.unique.1.abs_diff(head.unique.1) <= JC_CLUSTER_TOL_BP
    });
    if !same_flank {
        return false;
    }
    for (cand, class) in group.iter().zip(&classes) {
        if class.is_none() {
            match unique_side_outside_repeats(cand, fasta, repeats) {
                Some((_, pos, _)) if pos.abs_diff(head.unique.1) <= JC_CLUSTER_TOL_BP => {}
                _ => return false,
            }
        }
    }
    true
}

/// Representative of an identical-sequence group: prefer a member whose
/// repeat side anchors on an annotated terminus, then the strongest seed,
/// then the nearer unique↔terminus distance, then deterministic
/// coordinates.
fn pick_identical_group_representative(
    group: &[crate::jc_seq::CandidateJunction],
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> usize {
    let rank = |idx: usize| {
        let cand = &group[idx];
        let class = classify_unique_repeat_candidate(cand, fasta, repeats);
        let anchored = class.is_some() as usize;
        let dist = class.as_ref().map(|c| c.dist).unwrap_or(u64::MAX);
        (
            anchored,
            cand.seed_reads,
            u64::MAX - dist,
            (u64::MAX - cand.side1_pos_1),
            (u64::MAX - cand.side2_pos_1),
        )
    };
    let mut best = 0usize;
    let mut best_key = rank(0);
    for idx in 1..group.len() {
        let key = rank(idx);
        if key > best_key {
            best = idx;
            best_key = key;
        }
    }
    best
}

/// Cluster unique↔repeat representatives whose unique flanks sit on the
/// same contig and strand within [`JC_CLUSTER_TOL_BP`] and join the SAME
/// repeat copy (`(rep.start, rep.end)`) of one family at the same element
/// end (5' or 3'). One junction whose flank register is enumerated at
/// several 1–5 bp offsets otherwise becomes one construct per offset; the
/// second-pass aligner splits the reads across those near-identical
/// constructs and each diluted copy can be accepted as a separate
/// (spurious) junction. Registers joining DIFFERENT copies of the family
/// stay separate constructs: each copy's junction is independent evidence
/// and collapsing across copies deletes a copy's true junction while its
/// byte-identical siblings capture the reads (Wave 9, Wave 8 risk-register
/// fallback). Each cluster keeps one representative: strongest seed, then
/// the nearer unique↔terminus copy, then deterministic coordinates. The
/// representative carries the cluster's strongest seed. Registers joining
/// one copy are made CONTIGUOUS before the linear scan (copy identity
/// sorts BEFORE the flank position, Wave 10 / Codex Sol [P2]): a purely
/// positional sort lets another copy's register interleave between two
/// ≤5 bp registers of one copy, and the scan — which only compares the
/// cluster's last member — flushes there, re-emitting one construct per
/// offset.
fn cluster_unique_flank_registers(
    reps: Vec<crate::jc_seq::CandidateJunction>,
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> Vec<crate::jc_seq::CandidateJunction> {
    if reps.len() < 2 {
        return reps;
    }
    let mut keyed: Vec<(UniqueRepeatClass, crate::jc_seq::CandidateJunction)> = Vec::new();
    let mut passthrough: Vec<crate::jc_seq::CandidateJunction> = Vec::new();
    for cand in reps {
        match classify_unique_repeat_candidate(&cand, fasta, repeats) {
            Some(class) => keyed.push((class, cand)),
            None => passthrough.push(cand),
        }
    }
    // Wave 10 (Codex Sol [P2]): the repeat-copy identity must be the
    // PRIMARY sort key, not the unique-flank position. Sorting by position
    // first lets a register of copy B sit between two ≤5 bp registers of
    // copy A; the linear scan below only compares against the cluster's
    // LAST member, so it flushed at the copy-B interloper and emitted one
    // construct per offset again (the exact defect Target B of Wave 8 was
    // meant to remove). Grouping by (contig, strand, family, copy span,
    // element end) first makes every register joining one copy contiguous,
    // and the position sort inside that run lets all ≤5 bp registers of
    // that copy merge regardless of what other copies do at nearby
    // coordinates. Cross-copy merging remains impossible: copy span is
    // still a required-equal extend key below (Wave 9 invariant).
    keyed.sort_by(|a, b| {
        (
            a.0.unique.0,
            a.0.unique.2,
            &a.0.repeat_name,
            a.0.repeat_span,
            a.0.end_kind,
            a.0.unique.1,
        )
            .cmp(&(
                b.0.unique.0,
                b.0.unique.2,
                &b.0.repeat_name,
                b.0.repeat_span,
                b.0.end_kind,
                b.0.unique.1,
            ))
    });
    let mut out: Vec<crate::jc_seq::CandidateJunction> = Vec::with_capacity(keyed.len());
    let flush = |cluster: &mut Vec<(UniqueRepeatClass, crate::jc_seq::CandidateJunction)>,
                 out: &mut Vec<crate::jc_seq::CandidateJunction>| {
        if cluster.is_empty() {
            return;
        }
        let mut best = 0usize;
        for idx in 1..cluster.len() {
            let better = (cluster[idx].1.seed_reads, u64::MAX - cluster[idx].0.dist)
                > (cluster[best].1.seed_reads, u64::MAX - cluster[best].0.dist);
            if better {
                best = idx;
            }
        }
        let max_seed = cluster
            .iter()
            .map(|(_, cand)| cand.seed_reads)
            .max()
            .unwrap_or(0);
        let mut kept = cluster.swap_remove(best).1;
        kept.seed_reads = kept.seed_reads.max(max_seed);
        cluster.clear();
        out.push(kept);
    };
    let mut cluster: Vec<(UniqueRepeatClass, crate::jc_seq::CandidateJunction)> = Vec::new();
    for item in keyed {
        let extend = match cluster.last() {
            Some((head, _)) => {
                head.unique.0 == item.0.unique.0
                    && head.unique.2 == item.0.unique.2
                    && head.unique.1.abs_diff(item.0.unique.1) <= JC_CLUSTER_TOL_BP
                    && head.repeat_name == item.0.repeat_name
                    && head.repeat_span == item.0.repeat_span
                    && head.end_kind == item.0.end_kind
            }
            None => false,
        };
        if extend {
            cluster.push(item);
        } else {
            flush(&mut cluster, &mut out);
            cluster.push(item);
        }
    }
    flush(&mut cluster, &mut out);
    out.extend(passthrough);
    out.sort_by_key(|cand| {
        (
            cand.side1_contig,
            cand.side1_pos_1,
            cand.side2_contig,
            cand.side2_pos_1,
        )
    });
    out
}

/// Identical IS ends make one junction sequence per copy. A second-pass
/// aligner reports one of those references. Keep the nearer copy of that
/// sequence and give it the group's strongest seed, so a higher-seed far
/// copy cannot capture the reads. Groups that fail the unique↔repeat
/// classification (repeat↔repeat constructs, or members whose repeat side
/// landed in a copy interior) also collapse to ONE representative when
/// they describe one flank or one repeat↔repeat event: keeping N
/// byte-identical constructs in `junctions.fa` only dilutes the
/// second-pass reads across indistinguishable copies (Wave 8, Target A).
/// The grouping key is STRAND-CANONICAL: a construct and its reverse
/// complement are one double-stranded reference (bowtie2 indexes both
/// strands), so opposite-orientation duplicates must land in the same
/// group and collapse too, even when their unique flanks are distinct —
/// the aligner cannot tell them apart and emits both, splitting the
/// reads (Wave 11: the ldr↔IS150 reverse-complement duplicates).
/// Afterwards, unique-flank registers enumerated at 1–5 bp offsets are
/// clustered so near-identical constructs cannot split the reads
/// (Wave 8, Target B); the clustering is restricted to registers joining
/// the SAME repeat copy, so distinct copies of one family are never
/// merged (Wave 9).
pub(crate) fn collapse_identical_unique_repeat_copies(
    cands: Vec<crate::jc_seq::CandidateJunction>,
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
    flank: usize,
) -> Vec<crate::jc_seq::CandidateJunction> {
    // Key = lexicographically smaller of (sequence, reverse complement);
    // `fwd` records whether the member's own construct is the canonical
    // orientation. A group holding BOTH orientations of one canonical
    // sequence is an aligner-indistinguishable reverse-complement pair.
    let mut groups: BTreeMap<Vec<u8>, Vec<(crate::jc_seq::CandidateJunction, bool)>> =
        BTreeMap::new();
    for cand in cands {
        let seq = crate::jc_seq::junction_construct(fasta, &cand, flank).sequence();
        let rc_seq = crate::jc_seq::rc_dna(&seq);
        let (key_seq, fwd) = if seq <= rc_seq {
            (seq, true)
        } else {
            (rc_seq, false)
        };
        groups.entry(key_seq).or_default().push((cand, fwd));
    }
    let mut reps = Vec::new();
    for (_seq, group) in groups {
        let rc_paired = group.iter().any(|(_, fwd)| *fwd) && group.iter().any(|(_, fwd)| !*fwd);
        let group: Vec<crate::jc_seq::CandidateJunction> =
            group.into_iter().map(|(cand, _)| cand).collect();
        if let Some(winner) = nearer_unique_repeat_representative(&group, fasta, repeats) {
            let mut kept = group[winner].clone();
            kept.seed_reads = group
                .iter()
                .map(|cand| cand.seed_reads)
                .max()
                .unwrap_or(kept.seed_reads);
            reps.push(kept);
        } else if group.len() > 1
            && (rc_paired || identical_group_has_one_flank(&group, fasta, repeats))
        {
            let best = pick_identical_group_representative(&group, fasta, repeats);
            let mut kept = group[best].clone();
            kept.seed_reads = group
                .iter()
                .map(|cand| cand.seed_reads)
                .max()
                .unwrap_or(kept.seed_reads);
            reps.push(kept);
        } else {
            reps.extend(group);
        }
    }
    cluster_unique_flank_registers(reps, fasta, repeats)
}

pub(crate) fn second_pass_splits(
    fasta: &[FastaRecord],
    reads: &FastqInput,
    primary_data: &PrimaryBamData,
    contig_results: &[crate::engine::ContigPileup],
    opts: &EngineOptions,
    work: &Path,
) -> Result<SecondPassResult> {
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
    let mut repeat_seeds = Vec::new();

    let contig_names: Vec<String> = fasta.iter().map(|r| r.name.clone()).collect();
    for sam_name in [
        "stage2_r1.sam",
        "stage2_r2.sam",
        "stage2_pe.sam",
        "stage2_se.sam",
        "stage2.sam",
    ] {
        let sam_path = work.join(sam_name);
        if sam_path.is_file() {
            let contig_lens: Vec<u64> = fasta.iter().map(|r| r.seq.len() as u64).collect();
            if let Ok(split_result) =
                crate::split_seed::extract_candidate_junctions_with_repeat_evidence_from_sam_anchored(
                    &sam_path,
                    &contig_names,
                    &opts.repeats,
                    &contig_lens,
                )
            {
                debug!(
                    "prokadiff: extracted {} split-read candidate junctions from {}",
                    split_result.candidates.len(),
                    sam_name
                );
                all_cands.extend(split_result.candidates);
                repeat_seeds.extend(split_result.repeat_seeds);
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
    let all_cands = collapse_identical_unique_repeat_copies(all_cands, fasta, &opts.repeats, flank);

    let ref_len: usize = fasta.iter().map(|r| r.seq.len()).sum();
    let cands = rank_and_cap(all_cands, ref_len, flank);
    if cands.is_empty() {
        info!("prokadiff: second-pass skipped (no seed junctions)");
        return Ok(SecondPassResult {
            splits: Vec::new(),
            repeat_seeds,
        });
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
        return Ok(SecondPassResult {
            splits: Vec::new(),
            repeat_seeds,
        });
    }
    let jc_dir = work.join("jc");
    std::fs::create_dir_all(&jc_dir)?;
    let fa = jc_dir.join("junctions.fa");
    write_records(&recs, &fa)?;
    debug!(
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
        debug!("prokadiff: gzip FASTQ; second-pass uses all reads (no clip/unmapped filter)");
        reads
    } else {
        let pass_dir = jc_dir.join("reads");
        let (filtered, nkeep) = crate::align::filter_fastq_input(
            reads,
            &primary_data.second_pass_keep,
            &primary_data.second_pass_seen,
            &pass_dir,
        )?;
        debug!("prokadiff: second-pass {nkeep} FASTQ records after clip/unmapped filter");
        if nkeep == 0 {
            return Ok(SecondPassResult {
                splits: Vec::new(),
                repeat_seeds,
            });
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
    let (extra, stats) = parse_junction_bam(
        &jc_bam,
        &kept,
        &primary_data.primary_scores,
        &opts.repeats,
        &contig_names,
    )?;
    info!(
        "prokadiff: second-pass records considered={} kept={} \
         rejected(short_cover={} not_spanning={} worse_than_primary={})",
        stats.considered,
        stats.kept,
        stats.short_cover,
        stats.not_spanning,
        stats.worse_than_primary
    );
    Ok(SecondPassResult {
        splits: extra,
        repeat_seeds,
    })
}

/// Normalize JC coordinates against the repeat table:
///
/// 1. Microhomology slide (existing, `0 < overlap <= 4`): slide one side by
///    the overlap so the emitted record carries `overlap = 0`, matching
///    breseq 0.40.x absolute-0 normalization.
/// 2. Terminus snap: a side within [`JC_CLUSTER_TOL_BP`] of an annotated
///    repeat-copy terminus is snapped onto the terminus. Repeat-side
///    breakpoints anchor at copy termini (the junction read continues into
///    the element boundary); Stage-2 local alignments may land a few bp
///    inside or outside the terminus.
/// 3. Large-overlap absorption (`overlap > 4`): when exactly one side lies
///    inside an annotated repeat, sliding that side by the full overlap is
///    accepted only if the slid coordinate lands on that repeat's terminus
///    (within tolerance). A register that over-counted microhomology at a
///    copy boundary (e.g. pos = terminus + overlap) is thereby folded back
///    to the terminus with `overlap = 0`. Otherwise the original overlap is
///    preserved.
///
/// Returns the normalized `(pos1, pos2, overlap)`.
///
/// Terminus snap is included. Callers that turn a 1–5 bp pre-terminus
/// shortfall into a target-site duplication must run [`slide_jc_overlap`]
/// first, measure that shortfall, then [`snap_jc_termini`]. Snapping first
/// erases the shortfall and the duplication collapses.
pub(crate) fn normalize_jc_coords(
    s1: JcSideRef<'_>,
    s2: JcSideRef<'_>,
    overlap: i64,
    repeats: &[crate::fasta::RepeatRegion],
) -> (u64, u64, i64) {
    let (np1, np2, out_overlap) = slide_jc_overlap(s1, s2, overlap, repeats);
    let (np1, np2) = snap_jc_termini(
        JcSideRef {
            pos: np1,
            minus: s1.minus,
            contig: s1.contig,
        },
        JcSideRef {
            pos: np2,
            minus: s2.minus,
            contig: s2.contig,
        },
        repeats,
    );
    (np1, np2, out_overlap)
}

/// Overlap slide and large-overlap absorption, without terminus snap.
pub(crate) fn slide_jc_overlap(
    s1: JcSideRef<'_>,
    s2: JcSideRef<'_>,
    overlap: i64,
    repeats: &[crate::fasta::RepeatRegion],
) -> (u64, u64, i64) {
    let mut np1 = s1.pos;
    let mut np2 = s2.pos;
    let mut out_overlap = overlap;
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
        // The slide assigns the shared bases to exactly one side, so the
        // emitted overlap is an absolute 0 (docs/parity.md).
        out_overlap = 0;
    } else if overlap > 4 {
        // Slide the repeat side only, and only onto that repeat's terminus.
        // NB: the register can sit just OUTSIDE the copy (e.g. terminus +
        // overlap), so the termini-adjacent `find_repeat_at` window is the
        // right matcher here, not containment.
        let k = overlap as u64;
        let rep1 = find_repeat_at(repeats, s1.contig, s1.pos);
        let rep2 = find_repeat_at(repeats, s2.contig, s2.pos);
        match (rep1, rep2) {
            (Some(r1), None) => {
                let slid = if s1.minus {
                    np1.saturating_sub(k)
                } else {
                    np1.saturating_add(k)
                };
                if slid.abs_diff(r1.start) <= JC_CLUSTER_TOL_BP
                    || slid.abs_diff(r1.end) <= JC_CLUSTER_TOL_BP
                {
                    np1 = if slid.abs_diff(r1.start) <= slid.abs_diff(r1.end) {
                        r1.start
                    } else {
                        r1.end
                    };
                    out_overlap = 0;
                }
            }
            (None, Some(r2)) => {
                let slid = if s2.minus {
                    np2.saturating_sub(k)
                } else {
                    np2.saturating_add(k)
                };
                if slid.abs_diff(r2.start) <= JC_CLUSTER_TOL_BP
                    || slid.abs_diff(r2.end) <= JC_CLUSTER_TOL_BP
                {
                    np2 = if slid.abs_diff(r2.start) <= slid.abs_diff(r2.end) {
                        r2.start
                    } else {
                        r2.end
                    };
                    out_overlap = 0;
                }
            }
            _ => {}
        }
    }

    (np1, np2, out_overlap)
}

/// Land sides that sit within [`JC_CLUSTER_TOL_BP`] of an annotated copy
/// terminus exactly on that terminus. Independent of overlap.
pub(crate) fn snap_jc_termini(
    s1: JcSideRef<'_>,
    s2: JcSideRef<'_>,
    repeats: &[crate::fasta::RepeatRegion],
) -> (u64, u64) {
    let snap = |contig: &str, pos: u64| -> u64 {
        let Some(r) = find_repeat_at(repeats, contig, pos) else {
            return pos;
        };
        let d_start = pos.abs_diff(r.start);
        let d_end = pos.abs_diff(r.end);
        if d_start <= JC_CLUSTER_TOL_BP && d_start <= d_end {
            r.start
        } else if d_end <= JC_CLUSTER_TOL_BP {
            r.end
        } else {
            pos
        }
    };
    (snap(s1.contig, s1.pos), snap(s2.contig, s2.pos))
}
