//! Split-read candidate junction discovery from Stage 2 SAM alignments (published breseq methods).
//!
//! Evaluates pairs of sub-alignments from the same sequencing read using the published
//! mosaic pair criteria (`jc::is_candidate_junction`), computes breakpoint genomic coordinates
//! and strand orientations, and emits deduplicated `CandidateJunction` records.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind as NoodlesKind;

use crate::align::normalize_qname;
use crate::error::{EvidenceError, Result};
use crate::jc::{is_candidate_junction, SubAlignment};
use crate::jc_seq::CandidateJunction;

/// A local sub-alignment of a read against the reference genome.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SplitReadAlignment {
    pub qname: String,
    pub mate: u8,
    pub contig_idx: usize,
    pub ref_start_1: u64,
    pub ref_span: usize,
    pub read_start: usize,
    pub read_end: usize,
    pub read_len: usize,
    pub is_rc: bool,
}

/// Convert a SAM record into a `SplitReadAlignment` with canonical read coordinates.
pub fn parse_sam_split_record(
    rec: &sam::Record,
    name_to_idx: &HashMap<&str, usize>,
    default_mate: u8,
) -> Option<SplitReadAlignment> {
    let flags = rec.flags().ok()?;
    if flags.is_unmapped() {
        return None;
    }
    let name = rec.name().map(|n| n.to_string()).unwrap_or_default();
    if name.is_empty() {
        return None;
    }
    let mate = if flags.is_segmented() {
        if flags.is_first_segment() {
            1
        } else if flags.is_last_segment() {
            2
        } else {
            default_mate
        }
    } else {
        default_mate
    };
    let ref_name_bytes = rec.reference_sequence_name()?;
    let contig_name = std::str::from_utf8(ref_name_bytes.as_ref()).ok()?;
    let &contig_idx = name_to_idx.get(contig_name)?;
    let ref_start_1 = rec.alignment_start()?.ok()?.get() as u64;
    let is_rc = flags.is_reverse_complemented();

    let cigar = rec.cigar();
    let mut leading_s = 0usize;
    let mut trailing_s = 0usize;
    let mut query_match = 0usize;
    let mut ref_span = 0usize;
    let mut saw_match = false;

    for op in cigar.iter() {
        let op = op.ok()?;
        let len = op.len();
        match op.kind() {
            NoodlesKind::SoftClip => {
                if !saw_match {
                    leading_s += len;
                } else {
                    trailing_s += len;
                }
            }
            NoodlesKind::Match | NoodlesKind::SequenceMatch | NoodlesKind::SequenceMismatch => {
                saw_match = true;
                query_match += len;
                ref_span += len;
            }
            NoodlesKind::Insertion => {
                saw_match = true;
                query_match += len;
            }
            NoodlesKind::Deletion => {
                saw_match = true;
                ref_span += len;
            }
            _ => {}
        }
    }

    if query_match < 5 {
        return None;
    }

    let read_len = leading_s + query_match + trailing_s;
    if read_len == 0 {
        return None;
    }

    let (read_start, read_end) = if is_rc {
        (trailing_s, trailing_s + query_match)
    } else {
        (leading_s, leading_s + query_match)
    };

    Some(SplitReadAlignment {
        qname: normalize_qname(&name).to_string(),
        mate,
        contig_idx,
        ref_start_1,
        ref_span,
        read_start,
        read_end,
        read_len,
        is_rc,
    })
}

type JunctionKey = (usize, u64, bool, usize, u64, bool, i64);

/// Pair up sub-alignments of the same read, test mosaic pair criteria, and emit candidate junctions.
pub fn find_candidate_junctions(alignments: &[SplitReadAlignment]) -> Vec<CandidateJunction> {
    // Group by (qname, mate)
    let mut groups: HashMap<(&str, u8), Vec<&SplitReadAlignment>> = HashMap::new();
    for aln in alignments {
        groups.entry((&aln.qname, aln.mate)).or_default().push(aln);
    }

    let mut jc_counts: HashMap<JunctionKey, usize> = HashMap::new();

    for (_key, members) in groups {
        if members.len() < 2 || members.len() > 20 {
            continue;
        }
        let mut read_jcs = HashSet::new();
        for i in 0..members.len() {
            for j in (i + 1)..members.len() {
                let (a, b) = if members[i].read_start <= members[j].read_start {
                    (members[i], members[j])
                } else {
                    (members[j], members[i])
                };

                let read_len = a.read_len.max(b.read_len);
                let sub_a = SubAlignment {
                    read_start: a.read_start,
                    read_end: a.read_end,
                };
                let sub_b = SubAlignment {
                    read_start: b.read_start,
                    read_end: b.read_end,
                };

                if !is_candidate_junction(read_len, sub_a, sub_b) {
                    continue;
                }

                // Compute reference breakpoint and side strands.
                // Side 1 (from piece A, 5' on read):
                let (side1_pos_1, side1_minus) = if !a.is_rc {
                    (a.ref_start_1 + a.ref_span as u64 - 1, true)
                } else {
                    (a.ref_start_1, false)
                };

                // Side 2 (from piece B, 3' on read):
                let (side2_pos_1, side2_minus) = if !b.is_rc {
                    (b.ref_start_1, false)
                } else {
                    (b.ref_start_1 + b.ref_span as u64 - 1, true)
                };

                let overlap = a.read_end as i64 - b.read_start as i64;

                // Skip trivial collinear continuations along the same contig and strand
                if a.contig_idx == b.contig_idx
                    && side1_minus
                    && !side2_minus
                    && side1_pos_1 + 1 == side2_pos_1
                    && overlap == 0
                {
                    continue;
                }

                // Canonicalize: side1 <= side2. When sides are swapped, strands invert.
                let (s1_c, s1_p, s1_m, s2_c, s2_p, s2_m) =
                    if (a.contig_idx, side1_pos_1) <= (b.contig_idx, side2_pos_1) {
                        (
                            a.contig_idx,
                            side1_pos_1,
                            side1_minus,
                            b.contig_idx,
                            side2_pos_1,
                            side2_minus,
                        )
                    } else {
                        (
                            b.contig_idx,
                            side2_pos_1,
                            !side2_minus,
                            a.contig_idx,
                            side1_pos_1,
                            !side1_minus,
                        )
                    };

                read_jcs.insert((s1_c, s1_p, s1_m, s2_c, s2_p, s2_m, overlap));
            }
        }
        for k in read_jcs {
            *jc_counts.entry(k).or_insert(0) += 1;
        }
    }

    let mut out: Vec<CandidateJunction> = jc_counts
        .into_iter()
        .map(
            |((s1_c, s1_p, s1_m, s2_c, s2_p, s2_m, overlap), seed_reads)| CandidateJunction {
                side1_contig: s1_c,
                side1_pos_1: s1_p,
                side1_minus: s1_m,
                side2_contig: s2_c,
                side2_pos_1: s2_p,
                side2_minus: s2_m,
                overlap,
                seed_reads,
            },
        )
        .collect();

    out.sort_by(|a, b| {
        b.seed_reads
            .cmp(&a.seed_reads)
            .then_with(|| a.side1_contig.cmp(&b.side1_contig))
            .then_with(|| a.side1_pos_1.cmp(&b.side1_pos_1))
            .then_with(|| a.side2_contig.cmp(&b.side2_contig))
            .then_with(|| a.side2_pos_1.cmp(&b.side2_pos_1))
    });

    out
}

/// Extract candidate junctions from a SAM file (e.g. Stage 2 sensitive alignments).
pub fn extract_candidate_junctions_from_sam(
    sam_path: &Path,
    contig_names: &[String],
    default_mate: u8,
) -> Result<Vec<CandidateJunction>> {
    let name_to_idx: HashMap<&str, usize> = contig_names
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();

    let mut reader = File::open(sam_path)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let _header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;

    let mut alignments = Vec::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        if let Some(aln) = parse_sam_split_record(&rec, &name_to_idx, default_mate) {
            alignments.push(aln);
        }
    }

    Ok(find_candidate_junctions(&alignments))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_alignments_pair_generates_candidate_junction() {
        // 36 bp read:
        // Sub-alignment A: bases 0..20 map forward to contig 0 at 100..119 (ref_span=20).
        // Sub-alignment B: bases 20..36 map reverse to contig 1 at 200..215 (ref_span=16).
        let aln_a = SplitReadAlignment {
            qname: "read1".into(),
            mate: 0,
            contig_idx: 0,
            ref_start_1: 100,
            ref_span: 20,
            read_start: 0,
            read_end: 20,
            read_len: 36,
            is_rc: false,
        };
        let aln_b = SplitReadAlignment {
            qname: "read1".into(),
            mate: 0,
            contig_idx: 1,
            ref_start_1: 200,
            ref_span: 16,
            read_start: 20,
            read_end: 36,
            read_len: 36,
            is_rc: true,
        };

        let cands = find_candidate_junctions(&[aln_a, aln_b]);
        assert_eq!(cands.len(), 1);
        let c = &cands[0];
        // Side 1: contig 0, pos 100+20-1 = 119, minus=true
        assert_eq!(c.side1_contig, 0);
        assert_eq!(c.side1_pos_1, 119);
        assert!(c.side1_minus);
        // Side 2: contig 1, pos 200+16-1 = 215, minus=true (since is_rc=true)
        assert_eq!(c.side2_contig, 1);
        assert_eq!(c.side2_pos_1, 215);
        assert!(c.side2_minus);
        assert_eq!(c.overlap, 0);
        assert_eq!(c.seed_reads, 1);
    }

    #[test]
    fn rejects_pair_with_large_insertion_gap() {
        let aln_a = SplitReadAlignment {
            qname: "read1".into(),
            mate: 0,
            contig_idx: 0,
            ref_start_1: 100,
            ref_span: 10,
            read_start: 0,
            read_end: 10,
            read_len: 50,
            is_rc: false,
        };
        let aln_b = SplitReadAlignment {
            qname: "read1".into(),
            mate: 0,
            contig_idx: 0,
            ref_start_1: 200,
            ref_span: 15,
            read_start: 35, // gap = 35 - 10 = 25 > 20 bp limit
            read_end: 50,
            read_len: 50,
            is_rc: false,
        };

        let cands = find_candidate_junctions(&[aln_a, aln_b]);
        assert!(cands.is_empty(), "gap > 20 bp must be rejected");
    }

    #[test]
    fn rejects_trivial_collinear_continuation() {
        // Two consecutive matches along the same strand with 0 overlap
        let aln_a = SplitReadAlignment {
            qname: "read1".into(),
            mate: 0,
            contig_idx: 0,
            ref_start_1: 100,
            ref_span: 18,
            read_start: 0,
            read_end: 18,
            read_len: 36,
            is_rc: false,
        };
        let aln_b = SplitReadAlignment {
            qname: "read1".into(),
            mate: 0,
            contig_idx: 0,
            ref_start_1: 118, // 100 + 18 = 118, exact next base!
            ref_span: 18,
            read_start: 18,
            read_end: 36,
            read_len: 36,
            is_rc: false,
        };

        let cands = find_candidate_junctions(&[aln_a, aln_b]);
        assert!(
            cands.is_empty(),
            "trivial continuation along same strand must be skipped"
        );
    }
}
