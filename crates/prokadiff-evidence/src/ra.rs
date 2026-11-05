//! Clonal RA consensus on a synthetic pileup column (no BAM).
//!
//! Methods (https://breseq.barricklab.org/0.40.1/methods/): haploid consensus;
//! alignments with I/D > 2 bp are treated as JC, so RA only considers ≤2 bp
//! indels. Minimum end-trim of 1 is assumed already applied to observations.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BaseObs {
    /// Reference-aligned base: A/C/G/T or `b'-'` for a 1 bp deletion in this column.
    pub base: u8,
    pub strand: Strand,
    /// Identity of the aligned record that contributed this observation.
    /// Unique per read within one contig pileup (an ordinal assigned by the
    /// pileup builder, not the molecule id: supplementary alignments of one
    /// molecule share a molecule id but are distinct aligned records).
    /// Additive M11 subclonal field; consensus callers must not use it.
    pub read_id: u64,
}

/// One read's inserted sequence *after* a reference position.
/// `read_id` follows the same identity contract as [`BaseObs::read_id`]
/// and survives right-normalization and repeat masking because it is
/// attached to the entry before those steps.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InsertionObs {
    /// Inserted oligo (uppercase, reference-forward, right-normalized).
    pub seq: Vec<u8>,
    pub strand: Strand,
    pub read_id: u64,
}

/// One read's deletion EVENT spanning `span` reference bases starting at this
/// position (M11 Round-4 P1).
///
/// Structural record created by the pileup builder: exactly one entry per
/// read per deletion run, attached at the run's START column after the final
/// right-normalization position is fixed. The span is deliberately
/// independent of the per-base `b'-'` `BaseObs` entries: end-trimming
/// legitimately removes edge dash observations (e.g. a `70M2D1M` read whose
/// deletion right-normalizes onto the read's 3′ end loses its last dash),
/// so consumers must take deletion spans from these records and never
/// re-infer them from surviving dash runs.
/// `read_id` follows the same identity contract as [`BaseObs::read_id`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DeletionObs {
    /// Deleted reference bases. Recorded spans are 1–2: CIGAR deletions
    /// longer than 2 bp are diverted to the JC split path and never pushed
    /// as dash observations, so no longer run can be recorded here.
    pub span: u32,
    pub strand: Strand,
    pub read_id: u64,
}

#[derive(Clone, Debug, Default)]
pub struct PileupColumn {
    pub ref_base: u8,
    pub observations: Vec<BaseObs>,
    /// Inserted sequence *after* this reference position, one per supporting read.
    pub insertions: Vec<InsertionObs>,
    /// Deletion events STARTING at this reference position, one per
    /// supporting read (M11 subclonal; default-mode consensus callers must
    /// not use it).
    pub deletions: Vec<DeletionObs>,
}

#[derive(Clone, Debug)]
pub struct RaOptions {
    pub min_coverage: usize,
    pub min_frequency: f64,
}

impl Default for RaOptions {
    fn default() -> Self {
        Self {
            min_coverage: 5,
            min_frequency: 0.75,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ConsensusCall {
    MatchRef,
    Snp {
        alt: u8,
    },
    /// 1 or 2 bp deletion relative to the reference (RA subset).
    Del {
        size: u8,
    },
    /// 1 or 2 bp insertion after this position.
    Ins {
        seq: Vec<u8>,
    },
    RejectedStrandBias,
    InsufficientCoverage,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RaMetrics {
    pub depth: usize,
    pub support_reads: usize,
    pub frequency: f64,
}

pub fn call_consensus(col: &PileupColumn, opts: &RaOptions) -> ConsensusCall {
    call_consensus_with_metrics(col, opts).0
}

pub fn call_consensus_with_metrics(
    col: &PileupColumn,
    opts: &RaOptions,
) -> (ConsensusCall, Option<RaMetrics>) {
    if !col.insertions.is_empty() {
        return call_insertion_with_metrics(col, opts);
    }
    let cov = col.observations.len();
    if cov < opts.min_coverage {
        return (ConsensusCall::InsufficientCoverage, None);
    }
    let mut counts = [0usize; 5]; // A C G T -
    let mut plus = [0usize; 5];
    let mut minus = [0usize; 5];
    for obs in &col.observations {
        if let Some(i) = base_index(obs.base) {
            counts[i] += 1;
            match obs.strand {
                Strand::Plus => plus[i] += 1,
                Strand::Minus => minus[i] += 1,
            }
        }
    }
    let (best_i, best_n) = counts.iter().enumerate().max_by_key(|(_, n)| *n).unwrap();
    if *best_n == 0 {
        return (ConsensusCall::MatchRef, None);
    }
    let freq = *best_n as f64 / cov as f64;
    if freq < opts.min_frequency {
        return (ConsensusCall::MatchRef, None);
    }
    let consensus = index_base(best_i);
    if consensus == col.ref_base {
        return (ConsensusCall::MatchRef, None);
    }
    let plus_cov: usize = plus.iter().sum();
    let minus_cov: usize = minus.iter().sum();
    // Strand bias only when both strands are observed and the consensus is one-sided.
    if plus_cov > 0 && minus_cov > 0 && (plus[best_i] == 0 || minus[best_i] == 0) {
        return (ConsensusCall::RejectedStrandBias, None);
    }
    let metrics = RaMetrics {
        depth: cov,
        support_reads: *best_n,
        frequency: freq,
    };
    if consensus == b'-' {
        (ConsensusCall::Del { size: 1 }, Some(metrics))
    } else {
        (ConsensusCall::Snp { alt: consensus }, Some(metrics))
    }
}

fn call_insertion_with_metrics(
    col: &PileupColumn,
    opts: &RaOptions,
) -> (ConsensusCall, Option<RaMetrics>) {
    let cov = col.observations.len().max(col.insertions.len());
    if cov < opts.min_coverage {
        return (ConsensusCall::InsufficientCoverage, None);
    }
    // Majority inserted oligo of length 1–2.
    let mut best: Option<Vec<u8>> = None;
    let mut best_n = 0usize;
    let mut plus = 0usize;
    let mut minus = 0usize;
    for ins in &col.insertions {
        if ins.seq.is_empty() || ins.seq.len() > 2 {
            continue;
        }
        let n = col
            .insertions
            .iter()
            .filter(|other| other.seq == ins.seq)
            .count();
        if n > best_n {
            best_n = n;
            best = Some(ins.seq.clone());
            plus = col
                .insertions
                .iter()
                .filter(|other| other.seq == ins.seq && other.strand == Strand::Plus)
                .count();
            minus = col
                .insertions
                .iter()
                .filter(|other| other.seq == ins.seq && other.strand == Strand::Minus)
                .count();
        }
    }
    let Some(seq) = best else {
        return (ConsensusCall::MatchRef, None);
    };
    let freq = best_n as f64 / cov as f64;
    if freq < opts.min_frequency {
        return (ConsensusCall::MatchRef, None);
    }
    if plus == 0 || minus == 0 {
        return (ConsensusCall::RejectedStrandBias, None);
    }
    let metrics = RaMetrics {
        depth: cov,
        support_reads: best_n,
        frequency: freq,
    };
    (ConsensusCall::Ins { seq }, Some(metrics))
}

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        b'-' => Some(4),
        _ => None,
    }
}

fn index_base(i: usize) -> u8 {
    b"ACGT-"[i]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn col(ref_base: u8, bases: &[(u8, Strand)]) -> PileupColumn {
        PileupColumn {
            ref_base,
            observations: bases
                .iter()
                .enumerate()
                .map(|(i, (b, s))| BaseObs {
                    base: *b,
                    strand: *s,
                    read_id: i as u64,
                })
                .collect(),
            insertions: Vec::new(),
            deletions: Vec::new(),
        }
    }

    #[test]
    fn consensus_snp_when_majority_differs_on_both_strands() {
        let c = col(
            b'A',
            &[
                (b'G', Strand::Plus),
                (b'G', Strand::Plus),
                (b'G', Strand::Plus),
                (b'G', Strand::Minus),
                (b'G', Strand::Minus),
                (b'G', Strand::Minus),
            ],
        );
        assert_eq!(
            call_consensus(&c, &RaOptions::default()),
            ConsensusCall::Snp { alt: b'G' }
        );
    }

    #[test]
    fn one_bp_deletion_on_both_strands() {
        let c = col(
            b'A',
            &[
                (b'-', Strand::Plus),
                (b'-', Strand::Plus),
                (b'-', Strand::Plus),
                (b'-', Strand::Minus),
                (b'-', Strand::Minus),
                (b'-', Strand::Minus),
            ],
        );
        assert_eq!(
            call_consensus(&c, &RaOptions::default()),
            ConsensusCall::Del { size: 1 }
        );
    }

    #[test]
    fn two_bp_insertion_on_both_strands() {
        let mut c = col(b'T', &[(b'T', Strand::Plus); 6]);
        c.insertions = vec![
            InsertionObs {
                seq: b"AT".to_vec(),
                strand: Strand::Plus,
                read_id: 0,
            },
            InsertionObs {
                seq: b"AT".to_vec(),
                strand: Strand::Plus,
                read_id: 1,
            },
            InsertionObs {
                seq: b"AT".to_vec(),
                strand: Strand::Plus,
                read_id: 2,
            },
            InsertionObs {
                seq: b"AT".to_vec(),
                strand: Strand::Minus,
                read_id: 3,
            },
            InsertionObs {
                seq: b"AT".to_vec(),
                strand: Strand::Minus,
                read_id: 4,
            },
            InsertionObs {
                seq: b"AT".to_vec(),
                strand: Strand::Minus,
                read_id: 5,
            },
        ];
        assert_eq!(
            call_consensus(&c, &RaOptions::default()),
            ConsensusCall::Ins {
                seq: b"AT".to_vec()
            }
        );
    }

    #[test]
    fn snp_called_when_the_other_strand_has_no_coverage() {
        let c = col(b'T', &[(b'C', Strand::Plus); 6]);
        assert_eq!(
            call_consensus(&c, &RaOptions::default()),
            ConsensusCall::Snp { alt: b'C' }
        );
    }

    #[test]
    fn strand_bias_rejects_one_sided_variant() {
        let c = col(
            b'A',
            &[
                (b'G', Strand::Plus),
                (b'G', Strand::Plus),
                (b'G', Strand::Plus),
                (b'G', Strand::Plus),
                (b'G', Strand::Plus),
                (b'A', Strand::Minus),
            ],
        );
        assert_eq!(
            call_consensus(&c, &RaOptions::default()),
            ConsensusCall::RejectedStrandBias
        );
    }

    #[test]
    fn ambiguous_n_does_not_convert_to_a_or_call_fake_snp() {
        // Ref is G. We have 2 G's and 10 N's.
        // N must not be counted as A (which would call a G->A SNP).
        let mut obs = vec![(b'G', Strand::Plus), (b'G', Strand::Minus)];
        for _ in 0..5 {
            obs.push((b'N', Strand::Plus));
            obs.push((b'N', Strand::Minus));
        }
        let c = col(b'G', &obs);
        assert_eq!(
            call_consensus(&c, &RaOptions::default()),
            ConsensusCall::MatchRef
        );
    }

    #[test]
    fn ambiguous_n_only_column_matches_ref() {
        let obs = vec![
            (b'N', Strand::Plus),
            (b'N', Strand::Plus),
            (b'N', Strand::Minus),
            (b'N', Strand::Minus),
            (b'N', Strand::Plus),
            (b'N', Strand::Minus),
        ];
        let c = col(b'C', &obs);
        assert_eq!(
            call_consensus(&c, &RaOptions::default()),
            ConsensusCall::MatchRef
        );
    }

    #[test]
    fn majority_valid_base_with_some_n_calls_majority_correctly() {
        // Ref is A. We have 9 G's (both strands) and 1 N. Total cov = 10.
        // G frequency = 9/10 = 0.9 >= 0.8 min_frequency. Calls SNP A -> G.
        let obs = vec![
            (b'G', Strand::Plus),
            (b'G', Strand::Plus),
            (b'G', Strand::Plus),
            (b'G', Strand::Plus),
            (b'G', Strand::Plus),
            (b'G', Strand::Minus),
            (b'G', Strand::Minus),
            (b'G', Strand::Minus),
            (b'G', Strand::Minus),
            (b'N', Strand::Plus),
        ];
        let c = col(b'A', &obs);
        assert_eq!(
            call_consensus(&c, &RaOptions::default()),
            ConsensusCall::Snp { alt: b'G' }
        );
    }
}
