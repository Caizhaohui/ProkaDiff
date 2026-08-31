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
}

#[derive(Clone, Debug, Default)]
pub struct PileupColumn {
    pub ref_base: u8,
    pub observations: Vec<BaseObs>,
    /// Inserted sequence *after* this reference position, one per supporting read.
    pub insertions: Vec<(Vec<u8>, Strand)>,
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
            min_frequency: 0.8,
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

pub fn call_consensus(col: &PileupColumn, opts: &RaOptions) -> ConsensusCall {
    if !col.insertions.is_empty() {
        return call_insertion(col, opts);
    }
    let cov = col.observations.len();
    if cov < opts.min_coverage {
        return ConsensusCall::InsufficientCoverage;
    }
    let mut counts = [0usize; 5]; // A C G T -
    let mut plus = [0usize; 5];
    let mut minus = [0usize; 5];
    for obs in &col.observations {
        let i = base_index(obs.base);
        counts[i] += 1;
        match obs.strand {
            Strand::Plus => plus[i] += 1,
            Strand::Minus => minus[i] += 1,
        }
    }
    let (best_i, best_n) = counts.iter().enumerate().max_by_key(|(_, n)| *n).unwrap();
    let freq = *best_n as f64 / cov as f64;
    if freq < opts.min_frequency {
        return ConsensusCall::MatchRef;
    }
    let consensus = index_base(best_i);
    if consensus == col.ref_base {
        return ConsensusCall::MatchRef;
    }
    if plus[best_i] == 0 || minus[best_i] == 0 {
        return ConsensusCall::RejectedStrandBias;
    }
    if consensus == b'-' {
        return ConsensusCall::Del { size: 1 };
    }
    ConsensusCall::Snp { alt: consensus }
}

fn call_insertion(col: &PileupColumn, opts: &RaOptions) -> ConsensusCall {
    let cov = col.observations.len().max(col.insertions.len());
    if cov < opts.min_coverage {
        return ConsensusCall::InsufficientCoverage;
    }
    // Majority inserted oligo of length 1–2.
    let mut best: Option<Vec<u8>> = None;
    let mut best_n = 0usize;
    let mut plus = 0usize;
    let mut minus = 0usize;
    for (seq, _strand) in &col.insertions {
        if seq.is_empty() || seq.len() > 2 {
            continue;
        }
        let n = col.insertions.iter().filter(|(s, _)| s == seq).count();
        if n > best_n {
            best_n = n;
            best = Some(seq.clone());
            plus = col
                .insertions
                .iter()
                .filter(|(s, st)| s == seq && *st == Strand::Plus)
                .count();
            minus = col
                .insertions
                .iter()
                .filter(|(s, st)| s == seq && *st == Strand::Minus)
                .count();
        }
    }
    let Some(seq) = best else {
        return ConsensusCall::MatchRef;
    };
    let freq = best_n as f64 / cov as f64;
    if freq < opts.min_frequency {
        return ConsensusCall::MatchRef;
    }
    if plus == 0 || minus == 0 {
        return ConsensusCall::RejectedStrandBias;
    }
    ConsensusCall::Ins { seq }
}

fn base_index(b: u8) -> usize {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        b'-' => 4,
        _ => 0,
    }
}

fn index_base(i: usize) -> u8 {
    [b'A', b'C', b'G', b'T', b'-'][i]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn col(ref_base: u8, bases: &[(u8, Strand)]) -> PileupColumn {
        PileupColumn {
            ref_base,
            observations: bases
                .iter()
                .map(|(b, s)| BaseObs {
                    base: *b,
                    strand: *s,
                })
                .collect(),
            insertions: Vec::new(),
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
            (b"AT".to_vec(), Strand::Plus),
            (b"AT".to_vec(), Strand::Plus),
            (b"AT".to_vec(), Strand::Plus),
            (b"AT".to_vec(), Strand::Minus),
            (b"AT".to_vec(), Strand::Minus),
            (b"AT".to_vec(), Strand::Minus),
        ];
        assert_eq!(
            call_consensus(&c, &RaOptions::default()),
            ConsensusCall::Ins {
                seq: b"AT".to_vec()
            }
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
}
