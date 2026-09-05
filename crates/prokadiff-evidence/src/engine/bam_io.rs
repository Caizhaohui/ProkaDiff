//! BAM I/O and second-pass junction alignment parsing (noodles).

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind as NoodlesKind;

use crate::align::{keep_for_second_pass, normalize_qname};
use crate::error::{EvidenceError, Result};
use crate::fasta::FastaRecord;
use crate::jc::SubAlignment;
use crate::jc_seq::{spans_breakpoint, JC_MIN_COVER_BASES};
use crate::pileup::{AlignedRead, CigarKind, CigarOp, SplitCandidate, SplitOrigin};

pub(crate) struct PrimaryBamData {
    pub(crate) aligned: Vec<AlignedRead>,
    pub(crate) second_pass_keep: HashSet<String>,
    pub(crate) second_pass_seen: HashSet<String>,
    pub(crate) primary_scores: HashMap<String, i32>,
}

pub(crate) fn read_aligned_bam(bam_path: &Path, fasta: &[FastaRecord]) -> Result<Vec<AlignedRead>> {
    Ok(read_primary_bam(bam_path, fasta)?.aligned)
}

pub(crate) fn read_primary_bam(bam_path: &Path, fasta: &[FastaRecord]) -> Result<PrimaryBamData> {
    let name_to_idx: HashMap<&str, usize> = fasta
        .iter()
        .enumerate()
        .map(|(i, r)| (r.name.as_str(), i))
        .collect();

    let mut reader = File::open(bam_path).map(bam::io::Reader::new)?;
    let header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(n, _)| n.to_string())
        .collect();

    let mut aligned = Vec::new();
    let mut second_pass_keep = HashSet::new();
    let mut second_pass_seen = HashSet::new();
    let mut primary_scores = HashMap::new();

    for rec in reader.records() {
        let rec = rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        if rec.flags().is_secondary() {
            continue;
        }

        let name = rec.name().map(|n| n.to_string()).unwrap_or_default();
        let norm_name = if !name.is_empty() {
            let q = normalize_qname(&name).to_string();
            second_pass_seen.insert(q.clone());
            Some(q)
        } else {
            None
        };

        let cigar = noodles_cigar(&rec)?;
        let max_s = cigar
            .iter()
            .filter(|op| op.kind == CigarKind::SoftClip)
            .map(|op| op.len)
            .max()
            .unwrap_or(0);

        if let Some(ref q) = norm_name {
            if keep_for_second_pass(
                rec.flags().is_unmapped(),
                rec.flags().is_mate_unmapped(),
                max_s,
            ) {
                second_pass_keep.insert(q.clone());
            }
        }

        if rec.flags().is_unmapped() {
            continue;
        }

        if let Some(ref q) = norm_name {
            let (m, i, d, _) = cigar_counts(&cigar);
            let score = crate::jc_seq::cigar_match_indel_score(m, i, d);
            primary_scores
                .entry(q.clone())
                .and_modify(|s: &mut i32| *s = (*s).max(score))
                .or_insert(score);
        }

        let Some(Ok(rid)) = rec.reference_sequence_id() else {
            continue;
        };
        let contig_name = ref_names.get(rid).map(String::as_str).unwrap_or("");
        let Some(&contig_idx) = name_to_idx.get(contig_name) else {
            continue;
        };
        let Some(Ok(start)) = rec.alignment_start() else {
            continue;
        };
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(0);
        let as_score = rec
            .data()
            .get(&noodles::sam::alignment::record::data::field::Tag::ALIGNMENT_SCORE)
            .and_then(|r| r.ok())
            .and_then(|v| v.as_int());
        let xs_tag = noodles::sam::alignment::record::data::field::Tag::from([b'X', b'S']);
        let xs_score = rec
            .data()
            .get(&xs_tag)
            .and_then(|r| r.ok())
            .and_then(|v| v.as_int());
        let mapq = if let (Some(as_s), Some(xs_s)) = (as_score, xs_score) {
            if xs_s >= as_s {
                0
            } else {
                mapq
            }
        } else {
            mapq
        };
        let minus = rec.flags().is_reverse_complemented();
        let seq_len = rec.sequence().len();
        let mut seq = Vec::with_capacity(seq_len);
        for i in 0..seq_len {
            if let Some(b) = rec.sequence().get(i) {
                seq.push(b);
            }
        }
        aligned.push(AlignedRead {
            contig_idx,
            ref_start_0: start.get() as i64 - 1,
            minus,
            seq,
            cigar,
            mapq,
        });
    }
    Ok(PrimaryBamData {
        aligned,
        second_pass_keep,
        second_pass_seen,
        primary_scores,
    })
}

pub(crate) fn cigar_counts(cigar: &[CigarOp]) -> (usize, usize, usize, usize) {
    let mut n_match = 0usize;
    let mut n_ins = 0usize;
    let mut n_del = 0usize;
    for op in cigar {
        match op.kind {
            CigarKind::Match => n_match += op.len,
            CigarKind::Ins => n_ins += op.len,
            CigarKind::Del => n_del += op.len,
            _ => {}
        }
    }
    let q_cover = n_match + n_ins;
    (n_match, n_ins, n_del, q_cover)
}

pub(crate) fn construct_span(ref_start_0: i64, cigar: &[CigarOp]) -> (usize, usize) {
    let mut pos = ref_start_0.max(0) as usize;
    let start = pos;
    for op in cigar {
        match op.kind {
            CigarKind::Match | CigarKind::Del | CigarKind::Skip => pos += op.len,
            CigarKind::Ins | CigarKind::SoftClip | CigarKind::HardClip => {}
        }
    }
    (start, pos)
}

pub(crate) fn noodles_cigar(rec: &noodles::bam::Record) -> Result<Vec<CigarOp>> {
    let mut cigar = Vec::new();
    for op in rec.cigar().iter() {
        let op = op.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        let kind = match op.kind() {
            NoodlesKind::Match | NoodlesKind::SequenceMatch | NoodlesKind::SequenceMismatch => {
                CigarKind::Match
            }
            NoodlesKind::Insertion => CigarKind::Ins,
            NoodlesKind::Deletion => CigarKind::Del,
            NoodlesKind::SoftClip => CigarKind::SoftClip,
            NoodlesKind::HardClip => CigarKind::HardClip,
            NoodlesKind::Skip => CigarKind::Skip,
            _ => continue,
        };
        cigar.push(CigarOp {
            kind,
            len: op.len(),
        });
    }
    Ok(cigar)
}

/// Minimum construct bases a second-pass alignment must place on EACH side
/// of the breakpoint to count as junction evidence at all.
pub(crate) const JC_MIN_SPAN_EACH_SIDE: usize = 3;

/// Why a second-pass alignment did not become junction support.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct SecondPassReject {
    pub(crate) considered: usize,
    pub(crate) short_cover: usize,
    pub(crate) not_spanning: usize,
    pub(crate) worse_than_primary: usize,
    pub(crate) kept: usize,
}

/// Hash of a read name, used as second-pass supporting-read identity.
pub(crate) fn qname_id(qname: &str) -> u64 {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    let mut hasher = DefaultHasher::new();
    qname.hash(&mut hasher);
    hasher.finish()
}

/// Promote spanning second-pass alignments to `SplitCandidate`s.
pub(crate) fn parse_junction_bam(
    bam_path: &Path,
    candidates: &[(crate::jc_seq::CandidateJunction, usize)],
    ref_scores: &HashMap<String, i32>,
) -> Result<(Vec<SplitCandidate>, SecondPassReject)> {
    let mut stats = SecondPassReject::default();
    let mut reader = File::open(bam_path).map(bam::io::Reader::new)?;
    let header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(n, _)| n.to_string())
        .collect();
    let mut extra = Vec::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        if rec.flags().is_unmapped() {
            continue;
        }
        let Some(Ok(rid)) = rec.reference_sequence_id() else {
            continue;
        };
        let Some(name) = ref_names.get(rid) else {
            continue;
        };
        let Some(idx_str) = name.strip_prefix("jc") else {
            continue;
        };
        let Ok(ci) = idx_str.parse::<usize>() else {
            continue;
        };
        let Some((cand, breakpoint)) = candidates.get(ci) else {
            continue;
        };
        let breakpoint = *breakpoint;
        let Some(Ok(start)) = rec.alignment_start() else {
            continue;
        };
        stats.considered += 1;
        let cigar = noodles_cigar(&rec)?;
        let (n_match, n_ins, n_del, q_cover) = cigar_counts(&cigar);
        if q_cover < JC_MIN_COVER_BASES {
            stats.short_cover += 1;
            continue;
        }
        let (cstart, cend) = construct_span(start.get() as i64 - 1, &cigar);
        if !spans_breakpoint(cstart, cend, breakpoint, JC_MIN_SPAN_EACH_SIDE) {
            stats.not_spanning += 1;
            continue;
        }
        let jc_score = crate::jc_seq::cigar_match_indel_score(n_match, n_ins, n_del);
        let qname = rec.name().map(|n| n.to_string()).unwrap_or_default();
        let norm_qname = crate::align::normalize_qname(&qname);
        if let Some(&rs) = ref_scores.get(norm_qname) {
            if jc_score < rs {
                stats.worse_than_primary += 1;
                continue;
            }
        }
        let ov1 = breakpoint.saturating_sub(cstart);
        let ov2 = cend.saturating_sub(breakpoint);
        stats.kept += 1;
        extra.push(SplitCandidate {
            contig_idx: cand.side1_contig,
            side2_contig_idx: cand.side2_contig,
            minus: rec.flags().is_reverse_complemented(),
            side1_minus: cand.side1_minus,
            side2_minus: cand.side2_minus,
            left: SubAlignment {
                read_start: 0,
                read_end: ov1,
            },
            right: SubAlignment {
                read_start: ov1,
                read_end: ov1 + ov2,
            },
            side1_pos_1: cand.side1_pos_1,
            side2_pos_1: cand.side2_pos_1,
            overlap: cand.overlap,
            origin: SplitOrigin::SecondPass(qname_id(&qname)),
        });
    }
    Ok((extra, stats))
}
