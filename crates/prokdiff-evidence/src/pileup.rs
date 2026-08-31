//! CIGAR pileup into in-memory columns (no FASTQ). BAM records map onto `AlignedRead`.

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

#[derive(Clone, Debug)]
pub struct SplitCandidate {
    pub contig_idx: usize,
    pub minus: bool,
    pub left: SubAlignment,
    pub right: SubAlignment,
    pub side1_pos_1: u64,
    pub side2_pos_1: u64,
    pub overlap: i64,
}

pub fn apply_read(
    read: &AlignedRead,
    columns: &mut [PileupColumn],
    unique_depth: &mut [u32],
    splits: &mut Vec<SplitCandidate>,
) {
    if read.ref_start_0 < 0 || read.mapq < UNIQUE_MAPQ {
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

    for op in &read.cigar {
        match op.kind {
            CigarKind::Match => {
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
                        splits.push(SplitCandidate {
                            contig_idx: read.contig_idx,
                            minus: read.minus,
                            left,
                            right,
                            side1_pos_1: left_ref as u64 + 1,
                            side2_pos_1: left_ref as u64 + 1,
                            overlap: 0,
                        });
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
                        splits.push(SplitCandidate {
                            contig_idx: read.contig_idx,
                            minus: read.minus,
                            left,
                            right,
                            side1_pos_1: left_ref as u64 + 1,
                            side2_pos_1: (ref_pos + op.len) as u64 + 1,
                            overlap: -(op.len as i64),
                        });
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
                qpos += op.len;
            }
            CigarKind::HardClip => {}
            CigarKind::Skip => {
                ref_pos += op.len;
            }
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
        let mut alt = ref_seq.to_vec();
        alt[4] = b'T';
        for minus in [false, false, false, true, true, true] {
            apply_read(
                &match_read(&alt, 0, minus),
                &mut columns,
                &mut depth,
                &mut splits,
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
        apply_read(&read, &mut columns, &mut depth, &mut splits);
        assert_eq!(splits.len(), 1);
        assert!(splits[0].overlap < 0);
    }
}
