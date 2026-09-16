use std::collections::HashSet;

use prokadiff_gd::{GdEntry, GdKind};

use crate::classify::{ClassifiedMutation, MutationClass};
use crate::intended::IntendedEditAssessment;

/// Classification of a variant's origin with respect to editing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OriginStatus {
    /// Variant present in starter strain (parent background).
    StarterBackground,
    /// Variant confirmed to be part of declared intended edit.
    Intended,
    /// Post-edit differential variant (observed in edited clone, absent in starter).
    PostEditDifferential,
}

impl OriginStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::StarterBackground => "STARTER_BACKGROUND",
            Self::Intended => "INTENDED",
            Self::PostEditDifferential => "POST_EDIT_DIFF",
        }
    }
}

/// Physical scale of a genomic variant.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SizeClass {
    /// Small point mutation or short indel (<= 2 bp).
    Small,
    /// Large structural variant (> 2 bp DEL, MOB, JC, AMP, CON).
    Structural,
}

impl SizeClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Small => "SMALL",
            Self::Structural => "STRUCTURAL",
        }
    }
}

/// Relationship between an observed variant and declared intended edits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntendedRelation {
    /// Not associated with any declared intended locus.
    None,
    /// Matches declared on-target edit completely.
    Expected,
    /// Constituent event of an intended edit (e.g. one of two cassette junctions).
    PartOfExpectedEdit,
    /// Aberrant / unexpected event occurring at or near an intended locus.
    UnexpectedAtOnTargetLocus,
}

impl IntendedRelation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "NONE",
            Self::Expected => "EXPECTED",
            Self::PartOfExpectedEdit => "PART_OF_EXPECTED",
            Self::UnexpectedAtOnTargetLocus => "UNEXPECTED_AT_TARGET",
        }
    }
}

/// Spatial and sequence association with a predicted guide-homologous site.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GuideRelation {
    None,
    OnTarget,
    CandidateOffTarget {
        site_id: String,
        spacer_mismatches: u32,
        pam: String,
        distance_to_site: u64,
    },
}

impl GuideRelation {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::None => "NONE",
            Self::OnTarget => "ON_TARGET",
            Self::CandidateOffTarget { .. } => "CANDIDATE_OFF_TARGET",
        }
    }
}

/// Priority recommendation for manual review by experimentalists.
/// (Expresses human review priority, NOT biological safety).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ReviewPriority {
    Info,          // Expected complete edit or benign event
    Review,        // Distal small variant or coding missense
    HighAttention, // Structural rearrangement, MOB insertion, or unexpected on-target structure
}

impl ReviewPriority {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Info => "INFO",
            Self::Review => "REVIEW",
            Self::HighAttention => "HIGH_ATTENTION",
        }
    }
}

/// Summary of sequencing evidence supporting the variant call.
#[derive(Clone, Debug, PartialEq)]
pub struct EvidenceSummary {
    pub ra: bool,
    pub mc: bool,
    pub jc: bool,
    pub supporting_reads: Option<u64>,
    pub coverage: Option<f64>,
}

impl EvidenceSummary {
    pub fn format_brief(&self) -> String {
        format!(
            "RA={};MC={};JC={}",
            if self.ra { 1 } else { 0 },
            if self.mc { 1 } else { 0 },
            if self.jc { 1 } else { 0 }
        )
    }
}

/// Annotation for mobile genetic element (IS transposon) insertions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MobileElementAnnotation {
    pub family: Option<String>,
    pub element_name: Option<String>,
    pub insertion_site: Option<u64>,
    pub target_site_duplication: Option<String>,
    pub source_copy: Option<String>,
}

/// Annotation for repetitive genomic regions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepeatAnnotation {
    pub repeat_name: String,
    pub copy_count: Option<u32>,
}

/// Annotation for genomic features / genes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneAnnotation {
    pub locus_tag: Option<String>,
    pub gene_name: Option<String>,
    pub feature_type: String, // "CDS", "tRNA", "rRNA", "intergenic"
    pub product: Option<String>,
    pub consequence: Option<String>,
}

/// Unified multi-dimensionally annotated variant.
#[derive(Clone, Debug)]
pub struct AnnotatedVariant {
    pub variant_id: String, // E.g. "VAR_0001"
    pub entry: GdEntry,

    pub origin_status: OriginStatus,
    pub size_class: SizeClass,
    pub intended_relation: IntendedRelation,
    pub guide_relation: GuideRelation,

    pub mobile_element_relation: Option<MobileElementAnnotation>,
    pub repeat_relation: Option<RepeatAnnotation>,
    pub gene_annotation: Option<GeneAnnotation>,

    pub evidence: EvidenceSummary,
    pub review_priority: ReviewPriority,

    pub legacy_class: Option<MutationClass>,
}

/// Metadata about the analyzed sample and run context.
#[derive(Clone, Debug)]
pub struct SampleMetadata {
    pub starter_names: Vec<String>,
    pub edited_names: Vec<String>,
    pub reference_names: Vec<String>,
    pub editor: String,
    pub spacer: Option<String>,
    pub pam: Option<String>,
    pub threads: usize,
}

/// Verifiable technical provenance and oracle gating status.
#[derive(Clone, Debug)]
pub struct AnalysisProvenance {
    pub prokadiff_version: String,
    pub git_commit: String,
    pub reference_sha256: Option<String>,
    pub bowtie2_version: Option<String>,
    pub offtarget_search_status: String,
    pub cfd_scoring_status: String,
    pub hsu_scoring_status: String,
    pub bulge_search_status: String,
    pub run_timestamp: String,
}

/// Unified single aggregate root for post-edit genome audit results.
#[derive(Clone, Debug)]
pub struct AuditResult {
    pub sample: SampleMetadata,
    pub intended_edits: Vec<IntendedEditAssessment>,
    pub variants: Vec<AnnotatedVariant>,
    pub provenance: AnalysisProvenance,
}

/// Build a unified `AuditResult` from the classified mutations and intended assessments.
pub fn build_audit_result(
    sample: SampleMetadata,
    intended_assessments: Vec<IntendedEditAssessment>,
    unintended_mutations: &[ClassifiedMutation],
    intended_observed: &[GdEntry],
    provenance: AnalysisProvenance,
) -> AuditResult {
    let mut variants = Vec::new();
    let mut id_counter = 1usize;

    // Collect IDs belonging to intended edits
    let intended_matched_ids: HashSet<u32> = intended_assessments
        .iter()
        .flat_map(|a| a.matched_event_ids.iter().copied())
        .collect();
    let unexpected_ids: HashSet<u32> = intended_assessments
        .iter()
        .flat_map(|a| a.unexpected_event_ids.iter().copied())
        .collect();

    // 1. Process intended observed variants
    for entry in intended_observed {
        let is_part_of_multi = intended_assessments
            .iter()
            .any(|a| a.matched_event_ids.contains(&entry.id) && a.matched_event_ids.len() > 1);
        let intended_relation = if is_part_of_multi {
            IntendedRelation::PartOfExpectedEdit
        } else {
            IntendedRelation::Expected
        };

        let size_class = classify_size(&entry.kind, entry);
        let evidence = determine_evidence(entry);

        variants.push(AnnotatedVariant {
            variant_id: format!("VAR_{id_counter:04}"),
            entry: entry.clone(),
            origin_status: OriginStatus::Intended,
            size_class,
            intended_relation,
            guide_relation: GuideRelation::OnTarget,
            mobile_element_relation: extract_mobile_element(entry),
            repeat_relation: None,
            gene_annotation: None,
            evidence,
            review_priority: ReviewPriority::Info,
            legacy_class: None,
        });
        id_counter += 1;
    }

    // 2. Process post-edit differential unintended variants
    for cm in unintended_mutations {
        let entry = &cm.entry;
        let is_unexpected_on_target = unexpected_ids.contains(&entry.id);
        let is_intended_matched = intended_matched_ids.contains(&entry.id);

        let origin_status = if is_intended_matched {
            OriginStatus::Intended
        } else {
            OriginStatus::PostEditDifferential
        };

        let intended_relation = if is_unexpected_on_target {
            IntendedRelation::UnexpectedAtOnTargetLocus
        } else if is_intended_matched {
            IntendedRelation::Expected
        } else {
            IntendedRelation::None
        };

        let size_class = classify_size(&entry.kind, entry);
        let evidence = determine_evidence(entry);

        let guide_relation = match cm.class {
            MutationClass::NearHomolog => GuideRelation::CandidateOffTarget {
                site_id: format!("SITE_{id_counter:04}"),
                spacer_mismatches: cm.offtarget_mismatch.unwrap_or(0),
                pam: cm.pam_profile.clone().unwrap_or_default(),
                distance_to_site: cm.distance_to_site.unwrap_or(0),
            },
            _ => GuideRelation::None,
        };

        let review_priority = if is_unexpected_on_target {
            ReviewPriority::HighAttention
        } else {
            match size_class {
                SizeClass::Structural => ReviewPriority::HighAttention,
                SizeClass::Small => match guide_relation {
                    GuideRelation::CandidateOffTarget {
                        spacer_mismatches, ..
                    } if spacer_mismatches <= 1 => ReviewPriority::HighAttention,
                    _ => ReviewPriority::Review,
                },
            }
        };

        variants.push(AnnotatedVariant {
            variant_id: format!("VAR_{id_counter:04}"),
            entry: entry.clone(),
            origin_status,
            size_class,
            intended_relation,
            guide_relation,
            mobile_element_relation: extract_mobile_element(entry),
            repeat_relation: None,
            gene_annotation: None,
            evidence,
            review_priority,
            legacy_class: Some(cm.class),
        });
        id_counter += 1;
    }

    AuditResult {
        sample,
        intended_edits: intended_assessments,
        variants,
        provenance,
    }
}

fn classify_size(kind: &GdKind, entry: &GdEntry) -> SizeClass {
    match kind {
        GdKind::Mob | GdKind::Jc | GdKind::Amp | GdKind::Con => SizeClass::Structural,
        GdKind::Del => {
            let size = entry
                .fields
                .get(2)
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(1);
            if size > 2 {
                SizeClass::Structural
            } else {
                SizeClass::Small
            }
        }
        _ => SizeClass::Small,
    }
}

fn determine_evidence(entry: &GdEntry) -> EvidenceSummary {
    let mut ra = false;
    let mut mc = false;
    let mut jc = false;
    match entry.kind {
        GdKind::Snp | GdKind::Ins | GdKind::Sub => ra = true,
        GdKind::Del => {
            let size = entry
                .fields
                .get(2)
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(1);
            if size <= 2 {
                ra = true;
            } else {
                mc = true;
            }
        }
        GdKind::Jc => jc = true,
        GdKind::Mob => {
            jc = true;
        }
        GdKind::Amp | GdKind::Con => {
            mc = true;
        }
        GdKind::Ra => ra = true,
        GdKind::Mc => mc = true,
        GdKind::Un => {}
    }
    EvidenceSummary {
        ra,
        mc,
        jc,
        supporting_reads: None,
        coverage: None,
    }
}

fn extract_mobile_element(entry: &GdEntry) -> Option<MobileElementAnnotation> {
    if entry.kind != GdKind::Mob {
        return None;
    }
    // Fields for MOB: [seq_id, position, repeat_name, strand, ...]
    let element_name = entry.fields.get(2).cloned();
    let tsd = entry.attrs.get("repeat_seq").cloned();
    Some(MobileElementAnnotation {
        family: None,
        element_name,
        insertion_site: entry.position(),
        target_site_duplication: tsd,
        source_copy: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IntendedEditStatus;

    #[test]
    fn test_build_audit_result_combines_intended_and_unintended() {
        let sample = SampleMetadata {
            starter_names: vec!["s1.fq".into()],
            edited_names: vec!["e1.fq".into()],
            reference_names: vec!["ref.fa".into()],
            editor: "cas9".into(),
            spacer: Some("ACGT".into()),
            pam: Some("NGG".into()),
            threads: 4,
        };
        let prov = AnalysisProvenance {
            prokadiff_version: "0.2.0".into(),
            git_commit: "abc1234".into(),
            reference_sha256: None,
            bowtie2_version: Some("2.5.4".into()),
            offtarget_search_status: "EXACT".into(),
            cfd_scoring_status: "DISABLED".into(),
            hsu_scoring_status: "DISABLED".into(),
            bulge_search_status: "DISABLED".into(),
            run_timestamp: "2026-09-16T12:00:00Z".into(),
        };

        let intended_edit = IntendedEditAssessment {
            edit_id: "edit_1".into(),
            kind: "del".into(),
            seq_id: "chr".into(),
            expected_start: 1000,
            expected_end: 1200,
            status: IntendedEditStatus::Complete,
            matched_event_ids: vec![1],
            left_boundary: None,
            right_boundary: None,
            expected_size: Some(201),
            observed_size: Some(201),
            unexpected_event_ids: vec![],
            notes: vec![],
        };

        let intended_entry = GdEntry::del(1, "chr", 1000, 201);
        let unintended_cm = ClassifiedMutation {
            entry: GdEntry::snp(2, "chr", 5000, "T"),
            class: MutationClass::ScatteredSnv,
            pam_profile: None,
            offtarget_mismatch: None,
            distance_to_site: None,
            hypothesis: None,
        };

        let audit = build_audit_result(
            sample,
            vec![intended_edit],
            &[unintended_cm],
            &[intended_entry],
            prov,
        );

        assert_eq!(audit.variants.len(), 2);
        assert_eq!(audit.variants[0].origin_status, OriginStatus::Intended);
        assert_eq!(audit.variants[0].review_priority, ReviewPriority::Info);

        assert_eq!(
            audit.variants[1].origin_status,
            OriginStatus::PostEditDifferential
        );
        assert_eq!(audit.variants[1].review_priority, ReviewPriority::Review);
    }
}
