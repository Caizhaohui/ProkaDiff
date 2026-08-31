//! Missing-coverage (MC) seed-and-extend on unique-mapping depth.
//!
//! Methods: deletions seed at unique coverage 0 and extend until unique coverage
//! exceeds a cutoff derived from the coverage distribution. First-period MVP uses
//! a simple cutoff: max(2, median(nonzero unique depth) / 4).

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
}
