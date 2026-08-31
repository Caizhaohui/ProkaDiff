//! JC candidate and accept/reject rules (published methods, no FASTQ).
//!
//! Identifying candidate junctions (read-level mosaic pairs):
//! 1. One alignment begins with the first base of the read.
//! 2. Both alignments together cover >2 bases more than any single alignment.
//! 3. Both have ≥5 read bases that do not overlap the other.
//! 4. One has ≥10 unique read bases.
//! 5. ≤20 bp unique to the read between the matches.
//!
//! Accept (consensus):
//! 1. Reads on both strands of the predicted junction.
//! 2. ≥14 bp into each side of the reference.
//! 3. Each strand extends ≥9 bp into each side.
//! 4. Smallest reference overlap ≥3 bp on each side.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SubAlignment {
    /// 0-based inclusive start on the read.
    pub read_start: usize,
    /// 0-based exclusive end on the read.
    pub read_end: usize,
}

impl SubAlignment {
    pub fn len(self) -> usize {
        self.read_end.saturating_sub(self.read_start)
    }

    pub fn is_empty(self) -> bool {
        self.len() == 0
    }

    pub fn begins_at_read_start(self) -> bool {
        self.read_start == 0
    }

    /// Bases of `self` not overlapping `other` on the read.
    pub fn unique_vs(self, other: Self) -> usize {
        let overlap = overlap_len(self, other);
        self.len().saturating_sub(overlap)
    }
}

fn overlap_len(a: SubAlignment, b: SubAlignment) -> usize {
    let start = a.read_start.max(b.read_start);
    let end = a.read_end.min(b.read_end);
    end.saturating_sub(start)
}

fn gap_between(a: SubAlignment, b: SubAlignment) -> usize {
    if overlap_len(a, b) > 0 {
        return 0;
    }
    if a.read_end <= b.read_start {
        b.read_start - a.read_end
    } else {
        a.read_start - b.read_end
    }
}

/// Whether a pair of sub-alignments is a candidate mosaic junction.
pub fn is_candidate_junction(read_len: usize, a: SubAlignment, b: SubAlignment) -> bool {
    if read_len == 0 || a.is_empty() || b.is_empty() {
        return false;
    }
    if !(a.begins_at_read_start() || b.begins_at_read_start()) {
        return false;
    }
    let union = a.len() + b.len() - overlap_len(a, b);
    let longest_single = a.len().max(b.len());
    if union <= longest_single + 2 {
        return false;
    }
    let u_a = a.unique_vs(b);
    let u_b = b.unique_vs(a);
    if u_a < 5 || u_b < 5 {
        return false;
    }
    if u_a.max(u_b) < 10 {
        return false;
    }
    if gap_between(a, b) > 20 {
        return false;
    }
    true
}

#[derive(Clone, Debug, Default)]
pub struct JunctionSupport {
    pub plus_reads: usize,
    pub minus_reads: usize,
    /// Minimum overlap (bp) of supporting reads into side 1 of the reference.
    pub min_overlap_side1: usize,
    pub min_overlap_side2: usize,
    pub plus_side1: usize,
    pub plus_side2: usize,
    pub minus_side1: usize,
    pub minus_side2: usize,
}

pub fn accept_junction(s: &JunctionSupport) -> bool {
    if s.plus_reads == 0 || s.minus_reads == 0 {
        return false;
    }
    if s.min_overlap_side1 < 14 || s.min_overlap_side2 < 14 {
        return false;
    }
    if s.plus_side1 < 9 || s.plus_side2 < 9 || s.minus_side1 < 9 || s.minus_side2 < 9 {
        return false;
    }
    if s.min_overlap_side1.min(s.min_overlap_side2) < 3 {
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_textbook_mosaic_pair() {
        // 80 bp read: first 40 bp match left, last 45 bp match right (5 bp overlap on read).
        let a = SubAlignment {
            read_start: 0,
            read_end: 40,
        };
        let b = SubAlignment {
            read_start: 35,
            read_end: 80,
        };
        assert!(is_candidate_junction(80, a, b));
    }

    #[test]
    fn rejects_pair_that_does_not_start_at_read_begin() {
        let a = SubAlignment {
            read_start: 5,
            read_end: 40,
        };
        let b = SubAlignment {
            read_start: 35,
            read_end: 80,
        };
        assert!(!is_candidate_junction(80, a, b));
    }

    #[test]
    fn rejects_too_little_unique_sequence() {
        let a = SubAlignment {
            read_start: 0,
            read_end: 6,
        };
        let b = SubAlignment {
            read_start: 4,
            read_end: 10,
        };
        assert!(!is_candidate_junction(10, a, b));
    }

    #[test]
    fn accept_junction_requires_both_strands_and_overlap() {
        let good = JunctionSupport {
            plus_reads: 4,
            minus_reads: 3,
            min_overlap_side1: 20,
            min_overlap_side2: 18,
            plus_side1: 20,
            plus_side2: 18,
            minus_side1: 16,
            minus_side2: 15,
        };
        assert!(accept_junction(&good));

        let one_strand = JunctionSupport {
            minus_reads: 0,
            ..good
        };
        assert!(!accept_junction(&one_strand));

        let short = JunctionSupport {
            min_overlap_side1: 10,
            ..good
        };
        assert!(!accept_junction(&short));
    }
}
