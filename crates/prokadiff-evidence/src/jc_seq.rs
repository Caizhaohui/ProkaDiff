//! Candidate-junction sequence construction (published JC methods, no FASTQ).
//!
//! A candidate junction sequence is two reference flanks of length `read_len`
//! concatenated after assigning overlap bases to one side:
//! `2 * read_len - max(overlap, 0)` (36 bp reads overlapping 5 bp → 67 bp).
//!
//! Side strand is the direction the flank extends AWAY from the junction
//! (same convention as `SplitCandidate.side1_minus` / `side2_minus`).

use crate::fasta::FastaRecord;

/// Published: alignments covering fewer than this many read bases are discarded
/// when scoring second-pass hits.
pub const JC_MIN_COVER_BASES: usize = 28;

/// Softclips this long may seed a candidate (mosaic unique-bases floor).
/// Accept still requires ≥14 bp after second-pass support is aggregated.
pub const MIN_CLIP_FOR_SEED: usize = 6;

/// Genome-wide exact placement of the clipped sequence.
/// 6–11 mer in bacterial genomes is non-unique (millions of combinations / combinatorial explosion),
/// so genome-wide exact placement is restricted to S >= 12 (documented in docs/parity.md).
pub const MIN_CLIP_FOR_PLACE: usize = 12;

/// Cap on retained candidates (published: 5000 or 0.1× reference length).
pub const JC_MAX_CANDIDATES: usize = 5000;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CandidateJunction {
    pub side1_contig: usize,
    pub side1_pos_1: u64,
    pub side1_minus: bool,
    pub side2_contig: usize,
    pub side2_pos_1: u64,
    pub side2_minus: bool,
    pub overlap: i64,
    pub seed_reads: usize,
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

/// Reference bases of length `flank` entering the junction from Side 1.
/// In breseq GenomeDiff convention:
/// - `minus == true` (strand -1): connects to bases <= pos_1. The flank extends from lower
///   coordinates up to pos_1 in the forward orientation.
/// - `minus == false` (strand 1): connects to bases >= pos_1. The flank extends from higher
///   coordinates down to pos_1 (reverse-complemented).
pub fn side1_flank(seq: &[u8], pos_1: u64, minus: bool, flank: usize) -> Vec<u8> {
    if flank == 0 || seq.is_empty() || pos_1 == 0 {
        return Vec::new();
    }
    let pos_0 = pos_1 as usize - 1;
    if pos_0 >= seq.len() {
        return Vec::new();
    }
    if minus {
        let start = pos_0.saturating_add(1).saturating_sub(flank);
        seq[start..=pos_0]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect()
    } else {
        let end = (pos_0 + flank).min(seq.len());
        rc_dna(&seq[pos_0..end])
    }
}

/// Reference bases of length `flank` exiting the junction on Side 2.
/// In breseq GenomeDiff convention:
/// - `minus == false` (strand 1): connects to bases >= pos_2. The flank extends from pos_2
///   to higher coordinates in the forward orientation.
/// - `minus == true` (strand -1): connects to bases <= pos_2. The flank extends from pos_2
///   to lower coordinates (reverse-complemented).
pub fn side2_flank(seq: &[u8], pos_1: u64, minus: bool, flank: usize) -> Vec<u8> {
    if flank == 0 || seq.is_empty() || pos_1 == 0 {
        return Vec::new();
    }
    let pos_0 = pos_1 as usize - 1;
    if pos_0 >= seq.len() {
        return Vec::new();
    }
    if minus {
        let start = pos_0.saturating_add(1).saturating_sub(flank);
        rc_dna(&seq[start..=pos_0])
    } else {
        let end = (pos_0 + flank).min(seq.len());
        seq[pos_0..end]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect()
    }
}

/// A candidate-junction construct: the side 1 flank followed by the side 2
/// flank with `max(overlap, 0)` bases trimmed.
///
/// The breakpoint offset is [`Self::breakpoint`] = `left.len()`, NOT `flank`.
/// `side1_flank` clamps at contig boundaries and returns FEWER than
/// `flank` bases there, so a construct seeded near a contig end has its
/// breakpoint at a lower offset than `flank`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct JunctionConstruct {
    pub left: Vec<u8>,
    pub right: Vec<u8>,
}

impl JunctionConstruct {
    /// Offset of the junction breakpoint in [`Self::sequence`] (end of side 1).
    pub fn breakpoint(&self) -> usize {
        self.left.len()
    }

    pub fn len(&self) -> usize {
        self.left.len() + self.right.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn sequence(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.len());
        out.extend_from_slice(&self.left);
        out.extend_from_slice(&self.right);
        out
    }
}

/// Build both flanks of a candidate junction, dropping `max(overlap, 0)`
/// bases from side 2. Use [`JunctionConstruct::breakpoint`] rather than
/// `flank` to locate the junction in the emitted sequence.
pub fn junction_construct(
    fasta: &[FastaRecord],
    cand: &CandidateJunction,
    flank: usize,
) -> JunctionConstruct {
    let empty = JunctionConstruct {
        left: Vec::new(),
        right: Vec::new(),
    };
    let Some(s1) = fasta.get(cand.side1_contig) else {
        return empty;
    };
    let Some(s2) = fasta.get(cand.side2_contig) else {
        return empty;
    };
    let left = side1_flank(&s1.seq, cand.side1_pos_1, cand.side1_minus, flank);
    let right = side2_flank(&s2.seq, cand.side2_pos_1, cand.side2_minus, flank);
    let trim = cand.overlap.max(0) as usize;
    let right = if trim >= right.len() {
        Vec::new()
    } else {
        right[trim..].to_vec()
    };
    JunctionConstruct { left, right }
}

/// Concatenate two flanks, dropping `max(overlap, 0)` bases from side 2.
pub fn junction_sequence(fasta: &[FastaRecord], cand: &CandidateJunction, flank: usize) -> Vec<u8> {
    junction_construct(fasta, cand, flank).sequence()
}

/// Matched bases minus indel operations (published second-pass score).
pub fn cigar_match_indel_score(n_match: usize, n_ins: usize, n_del: usize) -> i32 {
    n_match as i32 - n_ins as i32 - n_del as i32
}

/// Whether a query alignment covering `construct_start_0..construct_end_0`
/// (0-based, exclusive end on the construct) spans the breakpoint with at
/// least `min_each` bases on each side.
pub fn spans_breakpoint(
    construct_start_0: usize,
    construct_end_0: usize,
    breakpoint: usize,
    min_each: usize,
) -> bool {
    if construct_end_0 <= construct_start_0 {
        return false;
    }
    let left = breakpoint.saturating_sub(construct_start_0);
    let right = construct_end_0.saturating_sub(breakpoint);
    left >= min_each && right >= min_each
}

/// Keep top-scoring seeds until 5000 candidates or 0.1× reference length
/// (published candidate-junction cap).
pub fn rank_and_cap(
    mut cands: Vec<CandidateJunction>,
    ref_len: usize,
    flank: usize,
) -> Vec<CandidateJunction> {
    if cands.iter().any(|c| c.seed_reads >= 2) {
        cands.retain(|c| c.seed_reads >= 2);
    }
    cands.sort_by(|a, b| {
        b.seed_reads
            .cmp(&a.seed_reads)
            .then(a.side1_contig.cmp(&b.side1_contig))
            .then(a.side1_pos_1.cmp(&b.side1_pos_1))
            .then(a.side2_contig.cmp(&b.side2_contig))
            .then(a.side2_pos_1.cmp(&b.side2_pos_1))
    });
    let max_cum = (ref_len / 10).max(flank.saturating_mul(2));
    let mut out = Vec::new();
    let mut acc = 0usize;
    for cand in cands {
        if out.len() >= JC_MAX_CANDIDATES {
            break;
        }
        let slen = (2 * flank).saturating_sub(cand.overlap.max(0) as usize);
        if !out.is_empty() && acc.saturating_add(slen) > max_cum {
            break;
        }
        acc = acc.saturating_add(slen);
        out.push(cand);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::FastaRecord;

    fn rec(seq: &[u8]) -> [FastaRecord; 1] {
        [FastaRecord {
            name: "chr".into(),
            seq: seq.to_vec(),
        }]
    }

    #[test]
    fn plus_flank_starts_at_junction_and_extends_higher() {
        // 1-based: 1=A 2=C 3=G 4=T 5=A 6=C 7=G 8=T 9=A 10=C. Plus @5 len 5 = ACGTA.
        let seq = b"ACGTACGTAC";
        assert_eq!(side2_flank(seq, 5, false, 5), b"ACGTA");
        assert_eq!(side1_flank(seq, 5, false, 5), b"TACGT");
    }

    #[test]
    fn minus_flank_is_rc_of_bases_ending_at_junction() {
        // Minus @5: bases 1–5 ACGTA, RC = TACGT.
        let seq = b"ACGTACGTAC";
        assert_eq!(side2_flank(seq, 5, true, 5), b"TACGT");
        assert_eq!(side1_flank(seq, 5, true, 5), b"ACGTA");
    }

    #[test]
    fn sequence_length_is_two_flanks_minus_positive_overlap() {
        // Published example: 36 bp reads, 5 bp overlap → 67 bp.
        let seq = vec![b'A'; 80];
        let fasta = rec(&seq);
        let cand = CandidateJunction {
            side1_contig: 0,
            side1_pos_1: 40,
            side1_minus: true,
            side2_contig: 0,
            side2_pos_1: 41,
            side2_minus: false,
            overlap: 5,
            seed_reads: 1,
        };
        let j = junction_sequence(&fasta, &cand, 36);
        assert_eq!(j.len(), 36 * 2 - 5);
        // Interior candidate: both flanks are full, so breakpoint == flank.
        assert_eq!(junction_construct(&fasta, &cand, 36).breakpoint(), 36);
    }

    #[test]
    fn breakpoint_follows_truncated_left_flank_at_contig_start() {
        // side1 is minus @10 with flank 36: only bases 1..=10 exist, so the
        // left arm is 10 bp, not 36. The breakpoint must be 10.
        let seq = vec![b'A'; 80];
        let fasta = rec(&seq);
        let cand = CandidateJunction {
            side1_contig: 0,
            side1_pos_1: 10,
            side1_minus: true,
            side2_contig: 0,
            side2_pos_1: 11,
            side2_minus: false,
            overlap: 0,
            seed_reads: 1,
        };
        let c = junction_construct(&fasta, &cand, 36);
        assert_eq!(c.left.len(), 10);
        assert_eq!(c.breakpoint(), 10);
        assert_eq!(c.len(), c.sequence().len());
        // A hit covering 5..20 straddles the true breakpoint (10) but would
        // look entirely left-of-breakpoint under the old `flank` constant.
        assert!(spans_breakpoint(5, 20, c.breakpoint(), 3));
        assert!(!spans_breakpoint(5, 20, 36, 3));
    }

    #[test]
    fn breakpoint_follows_truncated_right_flank_at_contig_end() {
        // side2 is plus @75 with flank 36: only 6 bases remain to the end.
        let seq = vec![b'C'; 80];
        let fasta = rec(&seq);
        let cand = CandidateJunction {
            side1_contig: 0,
            side1_pos_1: 40,
            side1_minus: true,
            side2_contig: 0,
            side2_pos_1: 75,
            side2_minus: false,
            overlap: 0,
            seed_reads: 1,
        };
        let c = junction_construct(&fasta, &cand, 36);
        assert_eq!(c.breakpoint(), 36);
        assert_eq!(c.right.len(), 6);
        assert_eq!(c.len(), 42);
    }

    #[test]
    fn rank_and_cap_survives_overlap_larger_than_two_flanks() {
        // `overlap` is an arbitrary i64 on a pub fn; 2*flank - overlap must
        // not underflow.
        let mut cand = dummy_cand(1, 1);
        cand.overlap = 10_000;
        let kept = rank_and_cap(vec![cand], 1_000, 20);
        assert_eq!(kept.len(), 1);
    }

    #[test]
    fn zero_overlap_keeps_both_full_flanks() {
        let seq = vec![b'T'; 80];
        let fasta = rec(&seq);
        let cand = CandidateJunction {
            side1_contig: 0,
            side1_pos_1: 40,
            side1_minus: true,
            side2_contig: 0,
            side2_pos_1: 41,
            side2_minus: false,
            overlap: 0,
            seed_reads: 1,
        };
        assert_eq!(junction_sequence(&fasta, &cand, 20).len(), 40);
    }

    #[test]
    fn negative_overlap_does_not_inflate_or_trim() {
        let seq = vec![b'G'; 80];
        let fasta = rec(&seq);
        let cand = CandidateJunction {
            side1_contig: 0,
            side1_pos_1: 30,
            side1_minus: true,
            side2_contig: 0,
            side2_pos_1: 40,
            side2_minus: false,
            overlap: -8,
            seed_reads: 1,
        };
        assert_eq!(junction_sequence(&fasta, &cand, 20).len(), 40);
    }

    #[test]
    fn match_indel_score_subtracts_gaps() {
        assert_eq!(cigar_match_indel_score(36, 0, 0), 36);
        assert_eq!(cigar_match_indel_score(30, 2, 1), 27);
    }

    #[test]
    fn spanning_requires_bases_on_both_sides_of_breakpoint() {
        assert!(spans_breakpoint(20, 50, 36, 3));
        assert!(!spans_breakpoint(0, 36, 36, 3)); // ends at breakpoint
        assert!(!spans_breakpoint(36, 72, 36, 3)); // starts at breakpoint
        assert!(!spans_breakpoint(30, 40, 36, 10)); // only 6 bp each side
    }

    fn dummy_cand(seed: usize, pos: u64) -> CandidateJunction {
        CandidateJunction {
            side1_contig: 0,
            side1_pos_1: pos,
            side1_minus: false,
            side2_contig: 0,
            side2_pos_1: pos + 100,
            side2_minus: true,
            overlap: 0,
            seed_reads: seed,
        }
    }

    #[test]
    fn rank_and_cap_keeps_highest_seed_then_stops_at_tenth_of_ref() {
        let cands = vec![dummy_cand(8, 3), dummy_cand(10, 1), dummy_cand(9, 2)];
        let kept = rank_and_cap(cands, 100, 20); // 0.1× = 10 bp; each seq = 40
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].seed_reads, 10);
        assert_eq!(kept[0].side1_pos_1, 1);
    }
}
