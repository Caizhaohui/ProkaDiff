//! Split-read candidate junction discovery from Stage 2 SAM alignments (published breseq methods).
//!
//! Evaluates pairs of sub-alignments from the same sequencing read using the published
//! mosaic pair criteria (`jc::is_candidate_junction`), computes breakpoint genomic coordinates
//! and strand orientations, and emits deduplicated `CandidateJunction` records.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind as NoodlesKind;

use crate::error::{EvidenceError, Result};
use crate::jc::{is_candidate_junction, SubAlignment};
use crate::jc_seq::CandidateJunction;
use crate::repeat_ambiguous::{
    canonical_molecule_identity, normalized_molecule_id, RepeatAmbiguousSeed, RepeatCopyPlacement,
    RepeatEvidenceSource, MAX_RETAINED_REPEAT_PLACEMENTS,
};

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
) -> Option<SplitReadAlignment> {
    let flags = rec.flags().ok()?;
    if flags.is_unmapped() {
        return None;
    }
    let name = rec.name().map(|n| n.to_string()).unwrap_or_default();
    if name.is_empty() {
        return None;
    }
    let (norm_qname, mate, _) = canonical_molecule_identity(&name, Some(&flags));
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
        qname: norm_qname,
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

/// Computational safety bound pending empirical benchmarking; not a claim of biological optimality.
pub const MAX_STAGE2_FAMILY_GEOMETRIES: usize = 5_000;

#[derive(Clone, Debug, Default)]
pub struct SplitSeedResult {
    pub candidates: Vec<CandidateJunction>,
    pub repeat_seeds: Vec<RepeatAmbiguousSeed>,
    pub rejected_trivial_continuation: usize,
}

/// Pair up sub-alignments of the same read, test mosaic pair criteria, and emit candidate junctions.
fn find_candidate_junctions_concrete(alignments: &[SplitReadAlignment]) -> Vec<CandidateJunction> {
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

pub fn find_candidate_junctions_with_repeat_evidence(
    alignments: &[SplitReadAlignment],
) -> SplitSeedResult {
    find_candidate_junctions_with_repeat_evidence_limit(alignments, MAX_STAGE2_FAMILY_GEOMETRIES)
}

pub fn find_candidate_junctions_with_repeat_evidence_limit(
    alignments: &[SplitReadAlignment],
    max_family_geometries: usize,
) -> SplitSeedResult {
    let candidates = find_candidate_junctions_concrete(alignments);
    let mut by_read: HashMap<(&str, u8), Vec<&SplitReadAlignment>> = HashMap::new();
    for alignment in alignments {
        by_read
            .entry((&alignment.qname, alignment.mate))
            .or_default()
            .push(alignment);
    }
    let mut repeat_seeds = Vec::new();
    let mut rejected_trivial_continuation = 0usize;
    for ((qname, mate), members) in by_read {
        if members.len() <= 20 {
            continue;
        }
        let mut families: BTreeMap<(usize, usize, usize, usize, bool), Vec<&SplitReadAlignment>> =
            BTreeMap::new();
        for member in members {
            families
                .entry((
                    member.read_start,
                    member.read_end,
                    member.read_len,
                    member.ref_span,
                    member.is_rc,
                ))
                .or_default()
                .push(member);
        }
        let mut family_values: Vec<Vec<&SplitReadAlignment>> = families
            .into_values()
            .map(|mut family| {
                family.sort_by_key(|alignment| {
                    (
                        alignment.contig_idx,
                        alignment.ref_start_1,
                        alignment.ref_span,
                        alignment.is_rc,
                    )
                });
                family
            })
            .collect();
        family_values.sort_by_key(|family| stage2_family_sort_key(family));
        let mut family_geometries = 0usize;
        let mut overflowed = false;
        let mut read_seeds = Vec::new();
        'family_pairs: for (i, left_family) in family_values.iter().enumerate() {
            for right_family in family_values.iter().skip(i + 1) {
                family_geometries = family_geometries.saturating_add(1);
                if family_geometries > max_family_geometries {
                    overflowed = true;
                    break 'family_pairs;
                }
                let (a, b) = if left_family[0].read_start <= right_family[0].read_start {
                    (left_family, right_family)
                } else {
                    (right_family, left_family)
                };
                let sub_a = SubAlignment {
                    read_start: a[0].read_start,
                    read_end: a[0].read_end,
                };
                let sub_b = SubAlignment {
                    read_start: b[0].read_start,
                    read_end: b[0].read_end,
                };
                if !is_candidate_junction(a[0].read_len.max(b[0].read_len), sub_a, sub_b) {
                    continue;
                }
                if a.len() > 1 && b.len() > 1 {
                    // Both sides are multicopy: neither side supplies a qualified unique anchor.
                    // Do not invent a unique anchor. Do not produce copy-specific resolution.
                    // Preserve observation as diagnostic evidence (CANDIDATE state).
                    let mut observed_placements = Vec::new();
                    for aln in a {
                        let (pos, minus) = side1_endpoint(aln);
                        observed_placements.push(RepeatCopyPlacement {
                            contig_idx: aln.contig_idx,
                            position_1: pos,
                            minus,
                        });
                    }
                    for aln in b {
                        let (pos, minus) = side2_endpoint(aln);
                        observed_placements.push(RepeatCopyPlacement {
                            contig_idx: aln.contig_idx,
                            position_1: pos,
                            minus,
                        });
                    }
                    observed_placements.sort();
                    observed_placements.dedup();
                    observed_placements.truncate(MAX_RETAINED_REPEAT_PLACEMENTS);
                    let placement_count = a.len() + b.len();
                    let placement_set_complete = placement_count <= MAX_RETAINED_REPEAT_PLACEMENTS;
                    let resource_complete = placement_set_complete;
                    let family = canonical_both_sides_family(a, b);

                    read_seeds.push(RepeatAmbiguousSeed {
                        source: RepeatEvidenceSource::Stage2Mosaic,
                        anchor_contig_idx: 0,
                        anchor_position_1: 0,
                        anchor_minus: false,
                        observed_placements,
                        eligible_copies: Vec::new(),
                        placement_count,
                        placement_family: family,
                        placement_sequence: Vec::new(),
                        placement_set_complete,
                        resource_complete,
                        resolution_complete: false,
                        overlap: a[0].read_end as i64 - b[0].read_start as i64,
                        molecule_id: normalized_molecule_id(qname, mate),
                        molecule_minus: a[0].is_rc,
                        clip_length: sub_b.len(),
                        aligned_length: sub_a.len(),
                        effective_mapq: 0,
                        alignment_score: None,
                        pair_geometry: "STAGE2_MOSAIC",
                        unique_anchor_qualified: false,
                        copy_resolved: false,
                        resolution_molecule_id: None,
                        reciprocal_evidence: false,
                        rejection_reason: Some("BOTH_SIDES_MULTICOPY"),
                        resource_limit_reason: None,
                    });
                    continue;
                }

                let (anchor, copies, anchor_is_unique, copied_from_a) =
                    if a.len() == 1 && b.len() > 1 {
                        let (a_pos, a_minus) = side1_endpoint(a[0]);
                        ((a[0].contig_idx, a_pos, a_minus), b, true, false)
                    } else if b.len() == 1 && a.len() > 1 {
                        let (b_pos, b_minus) = side2_endpoint(b[0]);
                        ((b[0].contig_idx, b_pos, b_minus), a, true, true)
                    } else {
                        continue;
                    };

                let mut copies = copies.to_vec();
                copies.sort_by_key(|copy| {
                    (copy.contig_idx, copy.ref_start_1, copy.ref_span, copy.is_rc)
                });
                let placement_count = copies.len();
                let placement_set_complete = placement_count <= MAX_RETAINED_REPEAT_PLACEMENTS;
                let resource_complete = placement_set_complete;

                let observed_placements: Vec<RepeatCopyPlacement> = copies
                    .iter()
                    .take(MAX_RETAINED_REPEAT_PLACEMENTS)
                    .map(|copy| {
                        let (position_1, minus) = if copied_from_a {
                            side1_endpoint(copy)
                        } else {
                            side2_endpoint(copy)
                        };
                        RepeatCopyPlacement {
                            contig_idx: copy.contig_idx,
                            position_1,
                            minus,
                        }
                    })
                    .collect();

                let eligible_copies: Vec<RepeatCopyPlacement> = copies
                    .iter()
                    .filter(|copy| {
                        let is_trivial = if copied_from_a {
                            is_trivial_continuation(copy, b[0])
                        } else {
                            is_trivial_continuation(a[0], copy)
                        };
                        !is_trivial
                    })
                    .take(MAX_RETAINED_REPEAT_PLACEMENTS)
                    .map(|copy| {
                        let (position_1, minus) = if copied_from_a {
                            side1_endpoint(copy)
                        } else {
                            side2_endpoint(copy)
                        };
                        RepeatCopyPlacement {
                            contig_idx: copy.contig_idx,
                            position_1,
                            minus,
                        }
                    })
                    .collect();

                if eligible_copies.is_empty() {
                    rejected_trivial_continuation += 1;
                    read_seeds.push(RepeatAmbiguousSeed {
                        source: RepeatEvidenceSource::Stage2Mosaic,
                        anchor_contig_idx: anchor.0,
                        anchor_position_1: anchor.1,
                        anchor_minus: anchor.2,
                        observed_placements,
                        eligible_copies: Vec::new(),
                        placement_count,
                        placement_family: canonical_stage2_family(&copies),
                        placement_sequence: Vec::new(),
                        placement_set_complete,
                        resource_complete,
                        resolution_complete: false,
                        overlap: a[0].read_end as i64 - b[0].read_start as i64,
                        molecule_id: normalized_molecule_id(qname, mate),
                        molecule_minus: a[0].is_rc,
                        clip_length: sub_b.len(),
                        aligned_length: sub_a.len(),
                        effective_mapq: 0,
                        alignment_score: None,
                        pair_geometry: "STAGE2_MOSAIC",
                        unique_anchor_qualified: false,
                        copy_resolved: false,
                        resolution_molecule_id: None,
                        reciprocal_evidence: false,
                        rejection_reason: Some("TRIVIAL_CONTINUATION_ONLY"),
                        resource_limit_reason: None,
                    });
                    continue;
                }

                let family = canonical_stage2_family(&copies);
                read_seeds.push(RepeatAmbiguousSeed {
                    source: RepeatEvidenceSource::Stage2Mosaic,
                    anchor_contig_idx: anchor.0,
                    anchor_position_1: anchor.1,
                    anchor_minus: anchor.2,
                    observed_placements,
                    eligible_copies,
                    placement_count,
                    placement_family: family,
                    placement_sequence: Vec::new(),
                    placement_set_complete,
                    resource_complete,
                    resolution_complete: false,
                    overlap: a[0].read_end as i64 - b[0].read_start as i64,
                    molecule_id: normalized_molecule_id(qname, mate),
                    molecule_minus: a[0].is_rc,
                    clip_length: sub_b.len(),
                    aligned_length: sub_a.len(),
                    effective_mapq: 0,
                    alignment_score: None,
                    pair_geometry: "STAGE2_MOSAIC",
                    unique_anchor_qualified: anchor_is_unique,
                    copy_resolved: false,
                    resolution_molecule_id: None,
                    reciprocal_evidence: false,
                    rejection_reason: None,
                    resource_limit_reason: None,
                });
            }
        }
        if overflowed {
            for seed in &mut read_seeds {
                seed.resource_limit_reason = Some("RESOURCE_LIMIT_STAGE2_FAMILY_GEOMETRY");
                seed.placement_set_complete = false;
                seed.resource_complete = false;
                seed.resolution_complete = false;
                seed.eligible_copies.clear();
            }
            let limit_seed = stage2_resource_limit_seed(qname, mate, &family_values);
            repeat_seeds.extend(read_seeds);
            repeat_seeds.push(limit_seed);
        } else {
            repeat_seeds.extend(read_seeds);
        }
    }
    repeat_seeds.sort_by_key(|seed| {
        (
            seed.anchor_contig_idx,
            seed.anchor_position_1,
            seed.anchor_minus,
            seed.placement_family.clone(),
            seed.molecule_id,
        )
    });
    SplitSeedResult {
        candidates,
        repeat_seeds,
        rejected_trivial_continuation,
    }
}

fn stage2_resource_limit_seed(
    qname: &str,
    mate: u8,
    families: &[Vec<&SplitReadAlignment>],
) -> RepeatAmbiguousSeed {
    let first = families.first().and_then(|family| family.first()).copied();
    let (anchor_contig_idx, anchor_position_1, anchor_minus, aligned_length) = match first {
        Some(alignment) => {
            let (position_1, minus) = side1_endpoint(alignment);
            (
                alignment.contig_idx,
                position_1,
                minus,
                alignment.read_end.saturating_sub(alignment.read_start),
            )
        }
        None => (0, 0, false, 0),
    };
    RepeatAmbiguousSeed {
        source: RepeatEvidenceSource::Stage2Mosaic,
        anchor_contig_idx,
        anchor_position_1,
        anchor_minus,
        observed_placements: Vec::new(),
        eligible_copies: Vec::new(),
        placement_count: families.len(),
        placement_family: format!("STAGE2_RESOURCE:{qname}:{mate}").into_bytes(),
        placement_sequence: Vec::new(),
        placement_set_complete: false,
        resource_complete: false,
        resolution_complete: false,
        overlap: 0,
        molecule_id: normalized_molecule_id(qname, mate),
        molecule_minus: false,
        clip_length: 0,
        aligned_length,
        effective_mapq: 0,
        alignment_score: None,
        pair_geometry: "STAGE2_RESOURCE_LIMIT",
        unique_anchor_qualified: false,
        copy_resolved: false,
        resolution_molecule_id: None,
        reciprocal_evidence: false,
        rejection_reason: Some("RESOURCE_LIMIT_STAGE2_FAMILY_GEOMETRY"),
        resource_limit_reason: Some("RESOURCE_LIMIT_STAGE2_FAMILY_GEOMETRY"),
    }
}

pub fn find_candidate_junctions(alignments: &[SplitReadAlignment]) -> Vec<CandidateJunction> {
    find_candidate_junctions_with_repeat_evidence(alignments).candidates
}

fn side1_endpoint(alignment: &SplitReadAlignment) -> (u64, bool) {
    if alignment.is_rc {
        (alignment.ref_start_1, false)
    } else {
        (alignment.ref_start_1 + alignment.ref_span as u64 - 1, true)
    }
}

fn side2_endpoint(alignment: &SplitReadAlignment) -> (u64, bool) {
    if alignment.is_rc {
        (alignment.ref_start_1 + alignment.ref_span as u64 - 1, true)
    } else {
        (alignment.ref_start_1, false)
    }
}

fn is_trivial_continuation(a: &SplitReadAlignment, b: &SplitReadAlignment) -> bool {
    let (a_pos, a_minus) = side1_endpoint(a);
    let (b_pos, b_minus) = side2_endpoint(b);
    a.contig_idx == b.contig_idx
        && a_minus
        && !b_minus
        && a_pos + 1 == b_pos
        && a.read_end == b.read_start
}

fn stage2_family_sort_key(family: &[&SplitReadAlignment]) -> (usize, usize, usize, usize, bool) {
    let first = family[0];
    (
        first.read_start,
        first.read_end,
        first.read_len,
        first.ref_span,
        first.is_rc,
    )
}

fn canonical_stage2_family(copies: &[&SplitReadAlignment]) -> Vec<u8> {
    let mut family = String::from("STAGE2");
    for copy in copies {
        use std::fmt::Write as _;
        let _ = write!(
            family,
            ":{}:{}:{}:{}:{}:{}:{}",
            copy.contig_idx,
            copy.ref_start_1,
            copy.ref_span,
            copy.read_start,
            copy.read_end,
            copy.read_len,
            copy.is_rc
        );
    }
    family.into_bytes()
}

fn canonical_both_sides_family(a: &[&SplitReadAlignment], b: &[&SplitReadAlignment]) -> Vec<u8> {
    let fam_a = canonical_stage2_family(a);
    let fam_b = canonical_stage2_family(b);
    let (first, second) = if fam_a <= fam_b {
        (fam_a, fam_b)
    } else {
        (fam_b, fam_a)
    };
    let mut out = b"BOTH_SIDES:".to_vec();
    out.extend_from_slice(&first);
    out.push(b'|');
    out.extend_from_slice(&second);
    out
}

/// Extract candidate junctions from a SAM file (e.g. Stage 2 sensitive alignments).
pub fn extract_candidate_junctions_from_sam(
    sam_path: &Path,
    contig_names: &[String],
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
        if let Some(aln) = parse_sam_split_record(&rec, &name_to_idx) {
            alignments.push(aln);
        }
    }

    Ok(find_candidate_junctions(&alignments))
}

pub fn extract_candidate_junctions_with_repeat_evidence_from_sam(
    sam_path: &Path,
    contig_names: &[String],
) -> Result<SplitSeedResult> {
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
    for record in reader.records() {
        let record = record.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        if let Some(alignment) = parse_sam_split_record(&record, &name_to_idx) {
            alignments.push(alignment);
        }
    }
    Ok(find_candidate_junctions_with_repeat_evidence(&alignments))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::{apply_read, AlignedRead, CigarKind, CigarOp};
    use crate::ra::PileupColumn;
    use crate::FastaRecord;

    fn repeat_group(copy_starts: &[u64]) -> Vec<SplitReadAlignment> {
        let anchor = SplitReadAlignment {
            qname: "repeat-read".into(),
            mate: 1,
            contig_idx: 0,
            ref_start_1: 100,
            ref_span: 20,
            read_start: 0,
            read_end: 20,
            read_len: 40,
            is_rc: false,
        };
        let mut alignments = vec![anchor];
        alignments.extend(
            copy_starts
                .iter()
                .copied()
                .map(|ref_start_1| SplitReadAlignment {
                    qname: "repeat-read".into(),
                    mate: 1,
                    contig_idx: 0,
                    ref_start_1,
                    ref_span: 20,
                    read_start: 20,
                    read_end: 40,
                    read_len: 40,
                    is_rc: false,
                }),
        );
        alignments
    }

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
    fn stage_two_repeat_groups_over_twenty_are_preserved_without_concrete_candidates() {
        let anchor = SplitReadAlignment {
            qname: "repeat-read".into(),
            mate: 1,
            contig_idx: 0,
            ref_start_1: 100,
            ref_span: 20,
            read_start: 0,
            read_end: 20,
            read_len: 40,
            is_rc: false,
        };
        let mut alignments = vec![anchor];
        for copy in 0..21 {
            alignments.push(SplitReadAlignment {
                qname: "repeat-read".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 1_000 + copy * 100,
                ref_span: 20,
                read_start: 20,
                read_end: 40,
                read_len: 40,
                is_rc: false,
            });
        }
        let result = find_candidate_junctions_with_repeat_evidence(&alignments);
        assert!(result.candidates.is_empty());
        assert_eq!(result.repeat_seeds.len(), 1);
        assert_eq!(result.repeat_seeds[0].placement_count, 21);
        assert!(result.repeat_seeds[0].unique_anchor_qualified);
    }

    #[test]
    fn stage_two_family_keeps_later_nontrivial_copy_after_trivial_first_copy() {
        let mut copies = vec![120];
        copies.extend((0..20).map(|index| 1_000 + index * 100));
        let result = find_candidate_junctions_with_repeat_evidence(&repeat_group(&copies));
        assert!(result.candidates.is_empty());
        assert_eq!(result.repeat_seeds.len(), 1);
        assert_eq!(result.repeat_seeds[0].placement_count, 21);
        // Trivial placement (120) remains observable in diagnostic provenance
        assert!(result.repeat_seeds[0]
            .observed_placements
            .iter()
            .any(|placement| placement.position_1 == 120));
        // Trivial placement is NOT in eligible_copies (not resolution-eligible)
        assert!(!result.repeat_seeds[0]
            .eligible_copies
            .iter()
            .any(|placement| placement.position_1 == 120));
        assert_eq!(result.repeat_seeds[0].eligible_copies.len(), 20);
    }

    #[test]
    fn stage_two_family_is_independent_of_raw_alignment_order() {
        let mut copies: Vec<u64> = (0..20).map(|index| 1_000 + index * 100).collect();
        copies.push(120);
        let forward = repeat_group(&copies);
        let mut reordered = forward.clone();
        reordered[1..].reverse();
        let forward_result = find_candidate_junctions_with_repeat_evidence(&forward);
        let reordered_result = find_candidate_junctions_with_repeat_evidence(&reordered);
        let forward_seed = &forward_result.repeat_seeds[0];
        let reordered_seed = &reordered_result.repeat_seeds[0];
        assert_eq!(
            forward_seed.placement_family,
            reordered_seed.placement_family
        );
        assert_eq!(
            forward_seed.observed_placements,
            reordered_seed.observed_placements
        );
        assert_eq!(forward_seed.eligible_copies, reordered_seed.eligible_copies);
        assert_eq!(forward_seed.placement_count, 21);
        assert_eq!(reordered_seed.placement_count, 21);
    }

    #[test]
    fn blocker_1_trivial_continuation_not_resolution_eligible_and_reciprocal_resolves_valid_only() {
        use crate::repeat_ambiguous::{
            resolve_repeat_copy, CopyResolutionEvidence, CopyResolutionSource,
            RepeatAmbiguousDiagnostics, RepeatAmbiguousState,
        };

        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 2000],
        }];

        // Anchor: 100..119 (endpoint 119, minus=true)
        // Copies:
        // - 120 (trivial continuation: endpoint 120, minus=false, 119+1==120)
        // - 1_000 (valid nontrivial copy 1)
        // - 1_100 (valid nontrivial copy 2)
        // Plus 19 extra copies to exceed the 20-member repeat threshold
        let mut copies = vec![120, 1_000, 1_100];
        copies.extend((0..18).map(|i| 1_200 + i * 50));

        let forward_result = find_candidate_junctions_with_repeat_evidence(&repeat_group(&copies));
        assert_eq!(forward_result.repeat_seeds.len(), 1);
        let seed = &forward_result.repeat_seeds[0];

        // 1. Trivial placement remains in diagnostic provenance
        assert!(seed.observed_placements.iter().any(|p| p.position_1 == 120));
        // 2. Trivial placement is NOT resolution-eligible
        assert!(!seed.eligible_copies.iter().any(|p| p.position_1 == 120));
        assert!(seed.eligible_copies.iter().any(|p| p.position_1 == 1_000));
        assert!(seed.eligible_copies.iter().any(|p| p.position_1 == 1_100));

        // 3. Reciprocal evidence targeting the trivial placement CANNOT resolve
        let trivial_resolution = resolve_repeat_copy(
            seed,
            &[CopyResolutionEvidence {
                source: CopyResolutionSource::ReciprocalBreakpoint,
                molecule_id: 999,
                copy: RepeatCopyPlacement {
                    contig_idx: 0,
                    position_1: 120,
                    minus: false,
                },
            }],
        );
        assert!(!trivial_resolution.copy_resolved);
        assert_eq!(
            trivial_resolution.eligible_copies.len(),
            seed.eligible_copies.len()
        );

        // 4. Reciprocal evidence targeting a valid nontrivial placement DOES resolve
        let valid_resolution = resolve_repeat_copy(
            seed,
            &[CopyResolutionEvidence {
                source: CopyResolutionSource::ReciprocalBreakpoint,
                molecule_id: 999,
                copy: RepeatCopyPlacement {
                    contig_idx: 0,
                    position_1: 1_000,
                    minus: false,
                },
            }],
        );
        assert!(valid_resolution.copy_resolved);
        assert_eq!(
            valid_resolution.eligible_copies,
            vec![RepeatCopyPlacement {
                contig_idx: 0,
                position_1: 1_000,
                minus: false,
            }]
        );

        // 5. Test with reordered input: reverse copies
        let mut reordered_copies = copies.clone();
        reordered_copies.reverse();
        let reordered_result =
            find_candidate_junctions_with_repeat_evidence(&repeat_group(&reordered_copies));
        let reordered_seed = &reordered_result.repeat_seeds[0];
        assert_eq!(seed.observed_placements, reordered_seed.observed_placements);
        assert_eq!(seed.eligible_copies, reordered_seed.eligible_copies);
        assert_eq!(seed.placement_family, reordered_seed.placement_family);

        // In aggregate diagnostics with sufficient support across both strands, state is SUPPORTED_AMBIGUOUS (never exact JC)
        let mut seed_minus = seed.clone();
        seed_minus.molecule_id = 2;
        seed_minus.molecule_minus = true;
        let mut seed_plus2 = seed.clone();
        seed_plus2.molecule_id = 3;
        seed_plus2.molecule_minus = false;
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(
            &fasta,
            vec![seed.clone(), seed_minus, seed_plus2],
        );
        assert_eq!(
            diagnostics.records[0].state,
            RepeatAmbiguousState::SupportedAmbiguous
        );
        assert_eq!(
            diagnostics.records[0].feasible_copies.len(),
            seed.eligible_copies.len()
        );
        // Diagnostic provenance retains observed placements
        assert!(diagnostics.records[0]
            .observed_placements
            .iter()
            .any(|p| p.position_1 == 120));
        assert!(!diagnostics.records[0]
            .feasible_copies
            .iter()
            .any(|p| p.position_1 == 120));
    }

    #[test]
    fn blocker_2_stage_two_overflow_leaves_no_usable_partial_seeds() {
        use crate::repeat_ambiguous::{RepeatAmbiguousDiagnostics, RepeatAmbiguousState};

        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 2000],
        }];

        // Construct 102 distinct alignment blocks for one read so that family_geometries
        // easily exceeds MAX_STAGE2_FAMILY_GEOMETRIES (102 * 101 / 2 = 5151 > 5000).
        let mut alignments = Vec::new();
        for index in 0..102 {
            alignments.push(SplitReadAlignment {
                qname: "overflow-read".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 100 + index as u64 * 100,
                ref_span: 20 + index,
                read_start: 0,
                read_end: 20,
                read_len: 40,
                is_rc: false,
            });
        }
        let result = find_candidate_junctions_with_repeat_evidence(&alignments);

        // 1. RESOURCE_LIMIT exists
        assert!(result
            .repeat_seeds
            .iter()
            .any(|s| { s.resource_limit_reason == Some("RESOURCE_LIMIT_STAGE2_FAMILY_GEOMETRY") }));

        // 2. Earlier members of that SAME family cannot resolve (eligible_copies empty, set incomplete)
        for s in &result.repeat_seeds {
            assert!(!s.placement_set_complete);
            assert!(!s.resource_complete);
            assert!(s.eligible_copies.is_empty());
            assert!(!s.copy_resolved);
        }

        // 3. Diagnostics from these seeds must be entirely RESOURCE_LIMIT, never RESOLVED
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(&fasta, result.repeat_seeds);
        assert!(!diagnostics.records.is_empty());
        for rec in &diagnostics.records {
            assert_eq!(rec.state, RepeatAmbiguousState::ResourceLimit);
            assert!(!rec.exact_breakpoint_supported);
            assert!(rec.feasible_copies.is_empty());
        }
    }

    #[test]
    fn blocker_3_both_sides_multicopy_rejects_without_inventing_anchor() {
        use crate::repeat_ambiguous::{RepeatAmbiguousDiagnostics, RepeatAmbiguousState};

        let fasta = vec![FastaRecord {
            name: "chr1".into(),
            seq: vec![b'A'; 5000],
        }];

        // Construct alignments where:
        // A. unique side (read_start 0..20, 1 placement) × multicopy side (read_start 20..40, 21 placements)
        let unique_a = SplitReadAlignment {
            qname: "read-unique".into(),
            mate: 1,
            contig_idx: 0,
            ref_start_1: 100,
            ref_span: 20,
            read_start: 0,
            read_end: 20,
            read_len: 40,
            is_rc: false,
        };
        let mut alns_a = vec![unique_a];
        for i in 0..21 {
            alns_a.push(SplitReadAlignment {
                qname: "read-unique".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 1_000 + i * 100,
                ref_span: 20,
                read_start: 20,
                read_end: 40,
                read_len: 40,
                is_rc: false,
            });
        }
        let res_a = find_candidate_junctions_with_repeat_evidence(&alns_a);
        assert_eq!(res_a.repeat_seeds.len(), 1);
        assert!(res_a.repeat_seeds[0].unique_anchor_qualified);
        assert_eq!(res_a.repeat_seeds[0].anchor_position_1, 119); // 100 + 20 - 1

        // B. multicopy side (2 placements) × multicopy side (21 placements)
        let multi_b1 = SplitReadAlignment {
            qname: "read-both-multi".into(),
            mate: 1,
            contig_idx: 0,
            ref_start_1: 100,
            ref_span: 20,
            read_start: 0,
            read_end: 20,
            read_len: 40,
            is_rc: false,
        };
        let multi_b2 = SplitReadAlignment {
            qname: "read-both-multi".into(),
            mate: 1,
            contig_idx: 0,
            ref_start_1: 500,
            ref_span: 20,
            read_start: 0,
            read_end: 20,
            read_len: 40,
            is_rc: false,
        };
        let mut alns_b = vec![multi_b1.clone(), multi_b2.clone()];
        for i in 0..21 {
            alns_b.push(SplitReadAlignment {
                qname: "read-both-multi".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 1_000 + i * 100,
                ref_span: 20,
                read_start: 20,
                read_end: 40,
                read_len: 40,
                is_rc: false,
            });
        }
        let res_b = find_candidate_junctions_with_repeat_evidence(&alns_b);
        // Requirement 2: Both sides multicopy MUST NOT disappear!
        // Diagnostic CANDIDATE survives
        assert_eq!(res_b.repeat_seeds.len(), 1);
        let seed_b = &res_b.repeat_seeds[0];
        assert!(!seed_b.unique_anchor_qualified);
        assert_eq!(seed_b.anchor_position_1, 0); // D: no arbitrary coordinate becomes anchor identity
        assert_eq!(seed_b.anchor_contig_idx, 0);
        assert_eq!(seed_b.eligible_copies.len(), 0);
        assert_eq!(seed_b.observed_placements.len(), 23); // contains BOTH sides (2 + 21)
        assert_eq!(seed_b.rejection_reason, Some("BOTH_SIDES_MULTICOPY"));

        // A & B: Diagnostics remain CANDIDATE and produce NO exact structural product
        let diag_b = RepeatAmbiguousDiagnostics::from_seeds(&fasta, res_b.repeat_seeds.clone());
        assert_eq!(diag_b.records.len(), 1);
        assert_eq!(diag_b.records[0].state, RepeatAmbiguousState::Candidate);
        assert!(!diag_b.records[0].exact_breakpoint_supported);
        assert!(res_b.candidates.is_empty(), "B: zero exact JC product");

        // C. Reordered placements of B -> logically identical result
        let mut alns_c = alns_b.clone();
        alns_c.reverse();
        let res_c = find_candidate_junctions_with_repeat_evidence(&alns_c);
        assert_eq!(res_c.repeat_seeds.len(), 1);
        assert_eq!(res_c.repeat_seeds[0], *seed_b);

        // D. Coordinate swapping on side A -> anchor position remains 0 (no coordinate becomes anchor)
        let multi_d1 = SplitReadAlignment {
            ref_start_1: 600,
            ..multi_b1
        };
        let multi_d2 = SplitReadAlignment {
            ref_start_1: 50,
            ..multi_b2
        };
        let mut alns_d = vec![multi_d1, multi_d2];
        for i in 0..21 {
            alns_d.push(SplitReadAlignment {
                qname: "read-both-multi".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 1_000 + i * 100,
                ref_span: 20,
                read_start: 20,
                read_end: 40,
                read_len: 40,
                is_rc: false,
            });
        }
        let res_d = find_candidate_junctions_with_repeat_evidence(&alns_d);
        assert_eq!(res_d.repeat_seeds.len(), 1);
        assert_eq!(res_d.repeat_seeds[0].anchor_position_1, 0);
        assert_eq!(res_d.repeat_seeds[0].anchor_contig_idx, 0);
        assert!(!res_d.repeat_seeds[0].unique_anchor_qualified);
    }

    #[test]
    fn blocker_4_sam_split_record_cross_path_molecule_identity() {
        use noodles::sam;
        use std::io::BufReader;

        let sam_data = b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
readA/1\t0\tchr1\t100\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
readA\t65\tchr1\t100\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
readA/2\t0\tchr1\t200\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
readA\t129\tchr1\t200\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
readPlain\t0\tchr1\t300\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n";

        let mut name_to_idx = HashMap::new();
        name_to_idx.insert("chr1", 0);

        let mut reader = sam::io::Reader::new(BufReader::new(&sam_data[..]));
        let _header = reader.read_header().unwrap();

        let mut records = Vec::new();
        for rec in reader.records() {
            let rec = rec.unwrap();
            let aln = parse_sam_split_record(&rec, &name_to_idx).expect("must parse");
            records.push(aln);
        }
        assert_eq!(records.len(), 5);

        // 1. readA/1 with missing mate flag vs readA with first-mate flag (65 = SEGMENTED | FIRST_SEGMENT)
        assert_eq!(records[0].qname, "readA");
        assert_eq!(records[0].mate, 1);
        assert_eq!(records[1].qname, "readA");
        assert_eq!(records[1].mate, 1);
        let id_rec0 = normalized_molecule_id(&records[0].qname, records[0].mate);
        let id_rec1 = normalized_molecule_id(&records[1].qname, records[1].mate);
        assert_eq!(id_rec0, id_rec1);

        // 2. readA/2 with missing mate flag vs readA with second-mate flag (129 = SEGMENTED | LAST_SEGMENT)
        assert_eq!(records[2].qname, "readA");
        assert_eq!(records[2].mate, 2);
        assert_eq!(records[3].qname, "readA");
        assert_eq!(records[3].mate, 2);
        let id_rec2 = normalized_molecule_id(&records[2].qname, records[2].mate);
        let id_rec3 = normalized_molecule_id(&records[3].qname, records[3].mate);
        assert_eq!(id_rec2, id_rec3);

        // 3. Different mates from same pair remain distinguishable
        assert_ne!(id_rec0, id_rec2);

        // 4. Primary BAM read "readPlain" (flags empty) vs Stage-2 SAM read "readPlain" (flags empty)
        // -> both get MATE_UNKNOWN (0) and identical molecule identity!
        let (prim_name, prim_mate, prim_id) = canonical_molecule_identity("readPlain", None);
        assert_eq!(prim_mate, 0);
        assert_eq!(prim_name, "readPlain");
        let stage2_id = normalized_molecule_id(&records[4].qname, records[4].mate);
        assert_eq!(records[4].mate, 0);
        assert_eq!(prim_id, stage2_id);
    }

    #[test]
    fn primary_and_stage_two_paths_share_one_molecule_identity() {
        let molecule_id = normalized_molecule_id("repeat-read", 1);
        let read = AlignedRead {
            contig_idx: 0,
            ref_start_0: 70,
            minus: false,
            seq: [vec![b'C'; 20], b"ACGTACGTACGTACGTACGT".to_vec()].concat(),
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: 20,
                },
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: 20,
                },
            ],
            mapq: 42,
            molecule_id,
        };
        let mut columns = vec![
            PileupColumn {
                ref_base: b'C',
                observations: Vec::new(),
                insertions: Vec::new(),
            };
            200
        ];
        let mut clips = Vec::new();
        apply_read(
            &read,
            &mut columns,
            &mut vec![0; 200],
            &mut vec![0; 200],
            &mut Vec::new(),
            &mut clips,
            None,
        );
        assert_eq!(clips[0].molecule_id, molecule_id);
        let copies: Vec<u64> = (0..21).map(|index| 1_000 + index * 100).collect();
        let result = find_candidate_junctions_with_repeat_evidence(&repeat_group(&copies));
        assert_eq!(result.repeat_seeds[0].molecule_id, molecule_id);
    }

    #[test]
    fn stage_two_family_geometry_overflow_is_a_resource_diagnostic() {
        let mut alignments = Vec::new();
        for index in 0..102 {
            alignments.push(SplitReadAlignment {
                qname: "overflow-read".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 100 + index as u64 * 100,
                ref_span: 20 + index,
                read_start: 0,
                read_end: 20,
                read_len: 40,
                is_rc: false,
            });
        }
        let result = find_candidate_junctions_with_repeat_evidence(&alignments);
        assert!(result.repeat_seeds.iter().any(|seed| {
            !seed.placement_set_complete
                && seed.rejection_reason == Some("RESOURCE_LIMIT_STAGE2_FAMILY_GEOMETRY")
        }));
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

    #[test]
    fn max_stage2_family_geometries_boundary_n_minus_1_n_n_plus_1() {
        // Boundary tests around MAX_STAGE2_FAMILY_GEOMETRIES = 5,000.
        // For a read with K distinct families, total pairs evaluated is K * (K - 1) / 2.
        // With K = 100: 100 * 99 / 2 = 4,950 pairs <= 5,000 (N-1 / N regime). Does NOT overflow.
        let mut alns_100 = Vec::new();
        for index in 0..100 {
            alns_100.push(SplitReadAlignment {
                qname: "test-read".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 100 + index as u64 * 100,
                ref_span: 20 + index,
                read_start: 0,
                read_end: 20,
                read_len: 40,
                is_rc: false,
            });
        }
        let res_100 = find_candidate_junctions_with_repeat_evidence(&alns_100);
        assert!(!res_100
            .repeat_seeds
            .iter()
            .any(|s| s.resource_limit_reason.is_some()));

        // With K = 101: 101 * 100 / 2 = 5,050 pairs.
        // Reaches 5,001 (N + 1) during evaluation, triggering RESOURCE_LIMIT_STAGE2_FAMILY_GEOMETRY!
        let mut alns_101 = alns_100.clone();
        alns_101.push(SplitReadAlignment {
            qname: "test-read".into(),
            mate: 1,
            contig_idx: 0,
            ref_start_1: 100 + 100_u64 * 100,
            ref_span: 20 + 100,
            read_start: 0,
            read_end: 20,
            read_len: 40,
            is_rc: false,
        });
        let res_101 = find_candidate_junctions_with_repeat_evidence(&alns_101);
        assert!(res_101
            .repeat_seeds
            .iter()
            .any(|s| { s.resource_limit_reason == Some("RESOURCE_LIMIT_STAGE2_FAMILY_GEOMETRY") }));
        for s in &res_101.repeat_seeds {
            assert!(!s.resource_complete);
            assert!(!s.resolution_complete);
            assert!(s.eligible_copies.is_empty());
        }

        // Test parameterized limit at exact G - 1, G, G + 1 boundary:
        // For a fixture with G = 10 geometries (5 families -> 10 pairs):
        // Ensure members.len() > 20 so Stage-2 repeat family analysis activates.
        let mut alns_5 = Vec::new();
        for i in 0..21 {
            alns_5.push(SplitReadAlignment {
                qname: "bound-read".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 100 + i as u64 * 10,
                ref_span: 20,
                read_start: 0,
                read_end: 20,
                read_len: 40,
                is_rc: false,
            });
        }
        for index in 1..5 {
            alns_5.push(SplitReadAlignment {
                qname: "bound-read".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 100 + index as u64 * 100,
                ref_span: 20 + index,
                read_start: 0,
                read_end: 20,
                read_len: 40,
                is_rc: false,
            });
        }
        // G - 1 = 9 -> overflows at 10th geometry
        let res_under = find_candidate_junctions_with_repeat_evidence_limit(&alns_5, 9);
        assert!(res_under
            .repeat_seeds
            .iter()
            .any(|s| s.resource_limit_reason.is_some()));

        // G = 10 -> exactly accommodates all 10 geometries without overflow
        let res_exact = find_candidate_junctions_with_repeat_evidence_limit(&alns_5, 10);
        assert!(!res_exact
            .repeat_seeds
            .iter()
            .any(|s| s.resource_limit_reason.is_some()));

        // G + 1 = 11 -> normal complete behavior
        let res_over = find_candidate_junctions_with_repeat_evidence_limit(&alns_5, 11);
        assert!(!res_over
            .repeat_seeds
            .iter()
            .any(|s| s.resource_limit_reason.is_some()));
    }

    #[test]
    fn trivial_only_family_diagnostics_observable_and_rejected() {
        use crate::repeat_ambiguous::RepeatAmbiguousDiagnostics;

        let fasta = vec![FastaRecord {
            name: "chr1".into(),
            seq: vec![b'A'; 2000],
        }];
        // Anchor piece A: read 0..20, ref 100..119, is_rc: false. (endpoint 119, -)
        let anchor = SplitReadAlignment {
            qname: "read-trivial".into(),
            mate: 1,
            contig_idx: 0,
            ref_start_1: 100,
            ref_span: 20,
            read_start: 0,
            read_end: 20,
            read_len: 40,
            is_rc: false,
        };
        // 21 placements of piece B, ALL of which are trivial collinear continuations:
        // piece B: read 20..40, ref 120..139, is_rc: false (endpoint 120, +).
        // Since 119 + 1 == 120 and read_end == read_start == 20: trivial continuation!
        let mut alns = vec![anchor];
        for _ in 0..21 {
            alns.push(SplitReadAlignment {
                qname: "read-trivial".into(),
                mate: 1,
                contig_idx: 0,
                ref_start_1: 120,
                ref_span: 20,
                read_start: 20,
                read_end: 40,
                read_len: 40,
                is_rc: false,
            });
        }
        let result = find_candidate_junctions_with_repeat_evidence(&alns);
        // Trivial continuation is recorded in SplitSeedResult
        assert_eq!(result.rejected_trivial_continuation, 1);
        assert_eq!(result.repeat_seeds.len(), 1);
        assert_eq!(
            result.repeat_seeds[0].rejection_reason,
            Some("TRIVIAL_CONTINUATION_ONLY")
        );

        // When processed through RepeatAmbiguousDiagnostics:
        // Observable in diagnostics count, but NOT promoted to CANDIDATE in records
        let diag = RepeatAmbiguousDiagnostics::from_seeds(&fasta, result.repeat_seeds);
        assert_eq!(diag.trivial_continuation_count, 1);
        assert_eq!(
            diag.records.len(),
            0,
            "must NOT promote trivial continuation to CANDIDATE"
        );
    }

    #[test]
    fn real_parser_reciprocal_resolution_is_reachable() {
        use crate::repeat_ambiguous::{
            RepeatAmbiguousDiagnostics, RepeatAmbiguousState, RepeatCopyPlacement,
        };
        use noodles::sam;
        use std::io::BufReader;

        let fasta = vec![FastaRecord {
            name: "chr1".into(),
            seq: vec![b'A'; 10_000],
        }];

        // Construct in-memory SAM records for two distinct molecules:
        // Molecule 1 ("read1"):
        //   - piece 1: anchor chr1:1000..1049 (len 50, + strand) -> endpoint (1049, -)
        //   - piece 2: 21 repeat copies starting at chr1:2000 (copy B1), 2100, 2200...
        // Molecule 2 ("read2"):
        //   - piece 1: anchor chr1:2000..2049 (len 50, - strand) -> endpoint (2000, +)
        //   - piece 2: 21 repeat copies starting at chr1:1000 (copy A1), 1100, 1200...
        let mut sam_str = String::from("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n");
        // read1 anchor: 50M50S at 1000
        sam_str.push_str("read1\t0\tchr1\t1000\t60\t50M50S\t*\t0\t0\t");
        sam_str.push_str(&"A".repeat(100));
        sam_str.push_str("\t*\n");
        // read1 21 copies of piece 2: 50S50M at 2000 + i * 100
        for i in 0..21 {
            sam_str.push_str(&format!(
                "read1\t0\tchr1\t{}\t0\t50S50M\t*\t0\t0\t{}\t*\n",
                2000 + i * 100,
                "A".repeat(100)
            ));
        }

        // read2 anchor: 50S50M at 2000, is_rc (flag 16)
        sam_str.push_str("read2\t16\tchr1\t2000\t60\t50S50M\t*\t0\t0\t");
        sam_str.push_str(&"A".repeat(100));
        sam_str.push_str("\t*\n");
        // read2 21 copies of piece 2: 50M50S at 1000 + i * 100, is_rc (flag 16)
        for i in 0..21 {
            sam_str.push_str(&format!(
                "read2\t16\tchr1\t{}\t0\t50M50S\t*\t0\t0\t{}\t*\n",
                1000 + i * 100,
                "A".repeat(100)
            ));
        }

        // Molecule 3 ("read3"): minus-strand read also anchored at 1049
        // Providing the minus-strand read required by the pre-SR1 accept_junction consensus gate
        // piece b anchor: 50M50S at 1000, flag 16 -> trailing_s = 50, read_start = 50, read_end = 100
        // side2_endpoint for is_rc: (1000 + 50 - 1, true) = (1049, true)
        sam_str.push_str("read3\t16\tchr1\t1000\t60\t50M50S\t*\t0\t0\t");
        sam_str.push_str(&"A".repeat(100));
        sam_str.push_str("\t*\n");
        // piece a: 21 copies 50S50M at 2000 + i * 100, flag 16
        // with flag 16: trailing_s = 0, read_start = 0, read_end = 50
        // side1_endpoint for is_rc: (2000 + i * 100, false)
        for i in 0..21 {
            sam_str.push_str(&format!(
                "read3\t16\tchr1\t{}\t0\t50S50M\t*\t0\t0\t{}\t*\n",
                2000 + i * 100,
                "A".repeat(100)
            ));
        }

        let mut name_to_idx = HashMap::new();
        name_to_idx.insert("chr1", 0);

        let mut reader = sam::io::Reader::new(BufReader::new(sam_str.as_bytes()));
        let _header = reader.read_header().unwrap();

        let mut alignments = Vec::new();
        for rec in reader.records() {
            let rec = rec.unwrap();
            let aln = parse_sam_split_record(&rec, &name_to_idx).expect("must parse");
            alignments.push(aln);
        }

        let split_res = find_candidate_junctions_with_repeat_evidence(&alignments);
        assert_eq!(split_res.repeat_seeds.len(), 3);
        let s1 = &split_res.repeat_seeds[0];
        let s2 = &split_res.repeat_seeds[1];
        let s3 = &split_res.repeat_seeds[2];
        assert_ne!(s1.molecule_id, s2.molecule_id);
        assert_ne!(s1.molecule_id, s3.molecule_id);
        assert_ne!(s2.molecule_id, s3.molecule_id);
        assert!(s1.unique_anchor_qualified);
        assert!(s2.unique_anchor_qualified);
        assert!(s3.unique_anchor_qualified);
        assert_eq!(s1.overlap, 0);
        assert_eq!(s2.overlap, 0);
        assert_eq!(s3.overlap, 0);
        assert_eq!(s1.pair_geometry, "STAGE2_MOSAIC");
        assert_eq!(s2.pair_geometry, "STAGE2_MOSAIC");
        assert_eq!(s3.pair_geometry, "STAGE2_MOSAIC");

        // Now resolve through RepeatAmbiguousDiagnostics
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(&fasta, split_res.repeat_seeds);
        assert_eq!(diagnostics.records.len(), 2);
        // Find record anchored at 1049 (from read1 + read3)
        let rec1 = diagnostics
            .records
            .iter()
            .find(|r| r.anchor_position_1 == 1049)
            .expect("record anchored at 1049 must exist");
        assert_eq!(rec1.state, RepeatAmbiguousState::Resolved);
        assert!(rec1.exact_breakpoint_supported);
        assert!(rec1.reciprocal_evidence);
        assert_eq!(rec1.plus_molecules, 1);
        assert_eq!(rec1.minus_molecules, 1);
        assert_eq!(rec1.molecule_ids.len(), 2);
        assert_eq!(
            rec1.feasible_copies,
            vec![RepeatCopyPlacement {
                contig_idx: 0,
                position_1: 2000,
                minus: false,
            }]
        );
    }

    #[test]
    fn duplicate_alignment_same_mate_deduplicates_molecule() {
        let sam_data = b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
readDup\t0\tchr1\t100\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
readDup\t0\tchr1\t200\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n";

        let mut name_to_idx = HashMap::new();
        name_to_idx.insert("chr1", 0);

        let mut reader = noodles::sam::io::Reader::new(std::io::BufReader::new(&sam_data[..]));
        let _header = reader.read_header().unwrap();

        let mut records = Vec::new();
        for rec in reader.records() {
            let rec = rec.unwrap();
            let aln = parse_sam_split_record(&rec, &name_to_idx).expect("must parse");
            records.push(aln);
        }
        assert_eq!(records.len(), 2);
        let id0 = normalized_molecule_id(&records[0].qname, records[0].mate);
        let id1 = normalized_molecule_id(&records[1].qname, records[1].mate);
        assert_eq!(
            id0, id1,
            "duplicate alignments of same mate share one molecule id"
        );
    }
}
