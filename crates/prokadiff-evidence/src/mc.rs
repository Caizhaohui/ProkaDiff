//! Missing-coverage (MC) seed-and-extend on unique-mapping depth.
//!
//! Methods: deletions seed at unique coverage 0 and extend until unique coverage
//! exceeds a cutoff derived from the coverage distribution. First-period MVP uses
//! a simple cutoff: max(2, median(nonzero unique depth) / 4).
//!
//! MC evidence is not automatically a consensus DEL. Unique-depth 0 also occurs
//! in multi-copy repeats / IS flanks (reads still map, but MAPQ is low). See
//! [`promote_mc_to_del`].

/// Default minimum unique-mapping gap length to promote MC → DEL without JC.
pub const MC_DEL_MIN_LEN: usize = 50;

/// JC endpoint must fall within this many bp of a gap flank (1-based).
pub const JC_DEL_WINDOW: u64 = 20;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MissingSpan {
    /// 0-based inclusive start.
    pub start_0: usize,
    /// 0-based exclusive end.
    pub end_0: usize,
}

impl MissingSpan {
    pub fn len(self) -> usize {
        self.end_0.saturating_sub(self.start_0)
    }

    pub fn is_empty(self) -> bool {
        self.len() == 0
    }

    /// 1-based closed interval start.
    pub fn start_1(self) -> u64 {
        self.start_0 as u64 + 1
    }

    /// 1-based closed interval end.
    pub fn end_1(self) -> u64 {
        self.end_0 as u64
    }
}

pub fn coverage_cutoff(unique_depth: &[u32]) -> u32 {
    let mut nz: Vec<u32> = unique_depth.iter().copied().filter(|&d| d > 0).collect();
    if nz.is_empty() {
        return 2;
    }
    nz.sort_unstable();
    let median = nz[nz.len() / 2];
    (median / 4).max(2)
}

/// Seed at unique depth 0; extend while depth < cutoff; keep spans of length ≥ `min_len`.
pub fn call_missing_coverage(unique_depth: &[u32], min_len: usize) -> Vec<MissingSpan> {
    if unique_depth.is_empty() {
        return Vec::new();
    }
    let cutoff = coverage_cutoff(unique_depth);
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < unique_depth.len() {
        if unique_depth[i] == 0 {
            let mut j = i + 1;
            while j < unique_depth.len() && unique_depth[j] < cutoff {
                j += 1;
            }
            if j - i >= min_len {
                out.push(MissingSpan {
                    start_0: i,
                    end_0: j,
                });
            }
            i = j;
        } else {
            i += 1;
        }
    }
    out
}

/// Shrink an MC span to the unique-depth-0 core (drop low-but-nonzero unique edges).
pub fn unique_zero_core(span: MissingSpan, unique_depth: &[u32]) -> Option<MissingSpan> {
    if span.end_0 > unique_depth.len() || span.start_0 >= span.end_0 {
        return None;
    }
    let mut a = span.start_0;
    let mut b = span.end_0;
    while a < b && unique_depth[a] != 0 {
        a += 1;
    }
    while b > a && unique_depth[b - 1] != 0 {
        b -= 1;
    }
    if b > a {
        Some(MissingSpan {
            start_0: a,
            end_0: b,
        })
    } else {
        None
    }
}

/// Unique-depth 0 **and** total-depth below cutoff (true missing, not a repeat).
pub fn true_missing_core(
    span: MissingSpan,
    unique_depth: &[u32],
    total_depth: &[u32],
    cutoff: u32,
) -> Option<MissingSpan> {
    if span.end_0 > unique_depth.len()
        || span.end_0 > total_depth.len()
        || span.start_0 >= span.end_0
    {
        return None;
    }
    let missing = |i: usize| unique_depth[i] == 0 && total_depth[i] < cutoff;
    let mut a = span.start_0;
    let mut b = span.end_0;
    while a < b && !missing(a) {
        a += 1;
    }
    while b > a && !missing(b - 1) {
        b -= 1;
    }
    if b > a {
        Some(MissingSpan {
            start_0: a,
            end_0: b,
        })
    } else {
        None
    }
}

fn flanks_have_unique_coverage(span: MissingSpan, unique_depth: &[u32], cutoff: u32) -> bool {
    if span.start_0 == 0 || span.end_0 >= unique_depth.len() {
        return false;
    }
    // Pileup unique depth is end-trimmed by 1 bp; allow that margin when
    // looking for unique flanks of a true-missing core.
    const TRIM: usize = 1;
    let left_ok = (0..span.start_0)
        .rev()
        .take(1 + TRIM)
        .any(|i| unique_depth[i] >= cutoff);
    let right_ok = (span.end_0..unique_depth.len())
        .take(1 + TRIM)
        .any(|i| unique_depth[i] >= cutoff);
    left_ok && right_ok
}

fn interior_total_below_cutoff(span: MissingSpan, total_depth: &[u32], cutoff: u32) -> bool {
    if span.end_0 > total_depth.len() || span.start_0 >= span.end_0 {
        return false;
    }
    total_depth[span.start_0..span.end_0]
        .iter()
        .all(|&d| d < cutoff)
}

fn jc_supports_both_flanks(span: MissingSpan, jc_endpoints_1: &[u64]) -> bool {
    if jc_endpoints_1.is_empty() {
        return false;
    }
    let left_anchor = span.start_1().saturating_sub(1);
    let right_anchor = span.end_1().saturating_add(1);
    let left_ok = jc_endpoints_1
        .iter()
        .any(|&p| p.abs_diff(left_anchor) <= JC_DEL_WINDOW && p <= span.start_1());
    let right_ok = jc_endpoints_1
        .iter()
        .any(|&p| p.abs_diff(right_anchor) <= JC_DEL_WINDOW && p >= span.end_1());
    left_ok && right_ok
}

/// Promote an MC span to consensus DEL (no JC). See [`promote_mc_to_del_with_jc`].
pub fn promote_mc_to_del(
    span: MissingSpan,
    unique_depth: &[u32],
    total_depth: &[u32],
    del_min_len: usize,
    min_span: usize,
) -> bool {
    promote_mc_to_del_with_jc(span, unique_depth, total_depth, del_min_len, min_span, &[])
}

/// First-period MC→DEL rule (clean-room vs published methods, not breseq source):
///
/// 1. Contig-terminal unique-depth 0 is not a DEL (breseq reports UN).
/// 2. Both unique flanks must have unique depth ≥ coverage cutoff.
/// 3. Interior **total** depth (any MAPQ) must stay below cutoff — unique-only
///    holes in repeats / IS copies stay MC.
/// 4. Without JC: unique-0 core length ≥ `del_min_len` (default 50).
/// 5. With JC endpoints near **both** flanks: length ≥ `min_span` (default 3).
pub fn promote_mc_to_del_with_jc(
    span: MissingSpan,
    unique_depth: &[u32],
    total_depth: &[u32],
    del_min_len: usize,
    min_span: usize,
    jc_endpoints_1: &[u64],
) -> bool {
    if unique_depth.is_empty() || total_depth.len() != unique_depth.len() {
        return false;
    }
    if span.len() < min_span {
        return false;
    }
    if span.start_0 == 0 || span.end_0 == unique_depth.len() {
        return false;
    }
    let cutoff = coverage_cutoff(unique_depth);
    let Some(core) = true_missing_core(span, unique_depth, total_depth, cutoff) else {
        return false;
    };
    if core.start_0 == 0 || core.end_0 == unique_depth.len() {
        return false;
    }
    if !flanks_have_unique_coverage(core, unique_depth, cutoff) {
        return false;
    }
    if !interior_total_below_cutoff(core, total_depth, cutoff) {
        return false;
    }
    if jc_supports_both_flanks(core, jc_endpoints_1) {
        return core.len() >= min_span;
    }
    core.len() >= del_min_len
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seeds_at_zero_and_extends_through_low_coverage() {
        let depth = [12u32, 12, 0, 0, 0, 1, 12, 12];
        let spans = call_missing_coverage(&depth, 3);
        assert_eq!(spans.len(), 1);
        assert_eq!(spans[0].start_0, 2);
        assert_eq!(spans[0].end_0, 6);
        assert_eq!(spans[0].start_1(), 3);
        assert_eq!(spans[0].end_1(), 6);
    }

    #[test]
    fn ignores_short_dips() {
        let depth = [10u32, 0, 10];
        assert!(call_missing_coverage(&depth, 3).is_empty());
    }

    /// Repeat-like: unique depth 0 but multi-mapping (total) coverage remains.
    /// Must not become a consensus DEL (Clonal leftover: 134 MC→DEL overcalls).
    #[test]
    fn repeat_like_unique_gap_is_not_promoted_to_del() {
        // 12 unique, then 8 bp unique-0 with total>0, then unique again.
        let unique = [12u32, 12, 12, 0, 0, 0, 0, 0, 0, 0, 0, 12, 12, 12];
        let total = [12u32, 12, 12, 8, 8, 9, 8, 7, 8, 8, 8, 12, 12, 12];
        let spans = call_missing_coverage(&unique, 3);
        assert_eq!(spans.len(), 1);
        assert!(!promote_mc_to_del(spans[0], &unique, &total, 50, 3));
    }

    #[test]
    fn short_internal_true_gap_is_not_promoted_to_del() {
        let unique = [12u32, 12, 0, 0, 0, 0, 12, 12];
        let total = unique;
        let spans = call_missing_coverage(&unique, 3);
        assert_eq!(spans.len(), 1);
        assert_eq!(spans[0].len(), 4);
        assert!(!promote_mc_to_del(spans[0], &unique, &total, 50, 3));
    }

    #[test]
    fn contig_terminal_zero_depth_is_not_promoted_to_del() {
        let unique = [0u32, 0, 0, 0, 12, 12, 12];
        let total = unique;
        let spans = call_missing_coverage(&unique, 3);
        assert_eq!(spans[0].start_0, 0);
        assert!(!promote_mc_to_del(spans[0], &unique, &total, 50, 3));
        let unique_r = [12u32, 12, 12, 0, 0, 0, 0];
        let spans_r = call_missing_coverage(&unique_r, 3);
        assert_eq!(spans_r[0].end_0, unique_r.len());
        assert!(!promote_mc_to_del(spans_r[0], &unique_r, &unique_r, 50, 3));
    }

    #[test]
    fn long_internal_zero_total_gap_promotes_to_del() {
        let mut unique = vec![12u32; 80];
        for d in unique.iter_mut().take(70).skip(10) {
            *d = 0;
        }
        let total = unique.clone();
        let spans = call_missing_coverage(&unique, 3);
        assert_eq!(spans.len(), 1);
        assert!(spans[0].len() >= 50);
        assert!(promote_mc_to_del(spans[0], &unique, &total, 50, 3));
    }

    #[test]
    fn jc_supported_short_gap_promotes_to_del() {
        let unique = [12u32, 12, 0, 0, 0, 0, 12, 12];
        let total = unique;
        let spans = call_missing_coverage(&unique, 3);
        // With JC at both unique flanks, even a short true gap may become DEL.
        assert!(promote_mc_to_del_with_jc(
            spans[0],
            &unique,
            &total,
            50,
            3,
            &[2, 7], // 1-based unique flanks just outside the gap
        ));
    }
}
