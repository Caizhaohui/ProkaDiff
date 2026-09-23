use prokadiff_gd::{GdEntry, GenomeDiff};

use crate::{
    build_audit_result, classify, AnalysisProvenance, ClassifyOptions, ClassifyResult, EditorKind,
    EvidenceObservation, IntendedEdit, IntendedEditStatus, IntendedEventRole, OriginStatus,
    RefContig, SampleMetadata,
};

fn genome(entries: Vec<GdEntry>) -> GenomeDiff {
    GenomeDiff {
        metadata: Vec::new(),
        entries,
    }
}

fn reference() -> Vec<RefContig> {
    vec![RefContig {
        name: "chr".into(),
        seq: vec![b'A'; 10_000],
    }]
}

fn options() -> ClassifyOptions {
    ClassifyOptions {
        editor: EditorKind::Dsb,
        spacer: None,
        pam: None,
        near_distance: 50,
        max_mismatches: 4,
        hypothesis: false,
    }
}

fn edit() -> IntendedEdit {
    IntendedEdit {
        edit_id: "edit_1".into(),
        seq_id: "chr".into(),
        start: 100,
        end: 100,
        ref_allele: "A".into(),
        alt: "T".into(),
        kind: "snp".into(),
    }
}

fn run(edited: Vec<GdEntry>, starter: Vec<GdEntry>, edits: Vec<IntendedEdit>) -> ClassifyResult {
    classify(
        &genome(edited),
        &genome(starter),
        &edits,
        &reference(),
        &options(),
    )
    .expect("classification succeeds")
}

#[test]
fn exact_intended_event_uses_stable_event_id() {
    let result = classify(
        &genome(vec![GdEntry::snp(42, "chr", 100, "T")]),
        &genome(vec![]),
        &[edit()],
        &reference(),
        &options(),
    )
    .expect("classification succeeds");
    let assessment = &result.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(
        assessment.matched_event_ids,
        vec![result.differential_events[0].event_id.clone()]
    );
    assert_eq!(
        assessment.event_relationships[0].role,
        IntendedEventRole::ExpectedConstituent
    );
    assert!(result.unintended.is_empty());
}

#[test]
fn absent_and_inherited_events_cannot_verify_edit() {
    let absent = run(vec![], vec![], vec![edit()]);
    let inherited = run(
        vec![GdEntry::snp(42, "chr", 100, "T")],
        vec![GdEntry::snp(8, "chr", 100, "T")],
        vec![edit()],
    );
    for result in [absent, inherited] {
        assert!(result.differential_events.is_empty());
        assert_eq!(
            result.intended_edit_assessments.as_ref().expect("declared")[0].status,
            IntendedEditStatus::Missing
        );
    }
}

#[test]
fn discordant_and_expected_events_at_same_locus_remain_distinct() {
    let result = run(
        vec![
            GdEntry::snp(11, "chr", 100, "T"),
            GdEntry::snp(12, "chr", 100, "C"),
        ],
        vec![],
        vec![edit()],
    );
    let assessment = &result.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(assessment.status, IntendedEditStatus::Partial);
    assert_eq!(result.intended_event_ids.len(), 1);
    assert_eq!(result.unintended.len(), 1);
    assert_eq!(result.unintended[0].entry.fields[2], "C");
    assert_ne!(
        result.intended_event_ids[0],
        result.unintended[0]
            .event_id
            .clone()
            .expect("differential id")
    );
}

#[test]
fn different_sequence_event_at_intended_locus_is_related_and_unmasked() {
    let result = run(
        vec![
            GdEntry::snp(11, "chr", 100, "T"),
            GdEntry::ins(12, "chr", 100, "G"),
        ],
        vec![],
        vec![edit()],
    );
    let assessment = &result.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(assessment.status, IntendedEditStatus::Partial);
    assert_eq!(assessment.unexpected_event_ids.len(), 1);
    assert_eq!(result.unintended.len(), 1);
    assert_eq!(result.unintended[0].entry.kind, prokadiff_gd::GdKind::Ins);
}

#[test]
fn overlapping_aberrant_structure_is_not_masked() {
    let result = run(
        vec![
            GdEntry::snp(1, "chr", 100, "T"),
            GdEntry::del(2, "chr", 99, 4),
        ],
        vec![],
        vec![edit()],
    );
    let assessment = &result.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(assessment.status, IntendedEditStatus::UnexpectedStructure);
    assert_eq!(assessment.unexpected_event_ids.len(), 1);
    assert_eq!(result.unintended.len(), 1);
    assert_eq!(
        result.unintended[0].event_id.as_ref(),
        assessment.unexpected_event_ids.first()
    );
}

#[test]
fn one_event_can_verify_multiple_declarations() {
    let mut second = edit();
    second.edit_id = "edit_2".into();
    let result = run(
        vec![GdEntry::snp(5, "chr", 100, "T")],
        vec![],
        vec![edit(), second],
    );
    let assessments = result.intended_edit_assessments.as_ref().expect("declared");
    assert_eq!(assessments.len(), 2);
    assert_eq!(
        assessments[0].matched_event_ids,
        assessments[1].matched_event_ids
    );
    assert_eq!(
        assessments[0].matched_event_ids[0],
        result.differential_events[0].event_id
    );
    assert_eq!(result.intended_event_ids.len(), 1);
}

#[test]
fn intended_identity_does_not_depend_on_raw_gd_id() {
    let first = run(vec![GdEntry::snp(1, "chr", 100, "T")], vec![], vec![edit()]);
    let second = run(
        vec![GdEntry::snp(900, "chr", 100, "T")],
        vec![],
        vec![edit()],
    );
    assert_eq!(
        first.intended_edit_assessments.as_ref().expect("declared")[0].matched_event_ids,
        second.intended_edit_assessments.as_ref().expect("declared")[0].matched_event_ids
    );
}

#[test]
fn deletion_candidate_selection_is_independent_of_input_order() {
    let intended = IntendedEdit {
        edit_id: "deletion".into(),
        seq_id: "chr".into(),
        start: 100,
        end: 109,
        ref_allele: ".".into(),
        alt: ".".into(),
        kind: "del".into(),
    };
    let nearby = GdEntry::del(20, "chr", 102, 8);
    let exact = GdEntry::del(10, "chr", 100, 10);
    let first = run(
        vec![nearby.clone(), exact.clone()],
        vec![],
        vec![intended.clone()],
    );
    let second = run(vec![exact, nearby], vec![], vec![intended]);
    let a = &first.intended_edit_assessments.as_ref().expect("declared")[0];
    let b = &second.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(a.matched_event_ids, b.matched_event_ids);
    assert_eq!(a.unexpected_event_ids, b.unexpected_event_ids);
    assert_eq!(a.left_boundary, b.left_boundary);
    assert_eq!(a.status, IntendedEditStatus::UnexpectedStructure);
    assert_eq!(a.matched_event_ids.len(), 1);
    assert_eq!(a.unexpected_event_ids.len(), 1);
}

#[test]
fn cassette_two_junctions_keep_distinct_event_ids() {
    let intended = IntendedEdit {
        edit_id: "cassette".into(),
        seq_id: "chr".into(),
        start: 500,
        end: 500,
        ref_allele: ".".into(),
        alt: ".".into(),
        kind: "cassette".into(),
    };
    let result = run(
        vec![
            GdEntry::jc(1, "chr", 500, "+", "chr", 800, "-", 0),
            GdEntry::jc(2, "chr", 500, "-", "chr", 900, "+", 0),
        ],
        vec![],
        vec![intended],
    );
    let assessment = &result.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(assessment.status, IntendedEditStatus::Complete);
    assert_eq!(assessment.matched_event_ids.len(), 2);
    assert!(assessment.matched_event_ids.iter().all(|id| result
        .differential_events
        .iter()
        .any(|event| &event.event_id == id)));
}

#[test]
fn cassette_pair_is_independent_of_gd_line_order() {
    let intended = IntendedEdit {
        edit_id: "cassette".into(),
        seq_id: "chr".into(),
        start: 500,
        end: 500,
        ref_allele: ".".into(),
        alt: ".".into(),
        kind: "cassette".into(),
    };
    let left = GdEntry::jc(1, "chr", 500, "+", "chr", 800, "-", 0);
    let right = GdEntry::jc(2, "chr", 500, "-", "chr", 900, "+", 0);
    let first = run(
        vec![left.clone(), right.clone()],
        vec![],
        vec![intended.clone()],
    );
    let second = run(vec![right, left], vec![], vec![intended]);
    assert_eq!(
        first.intended_edit_assessments,
        second.intended_edit_assessments
    );
    assert_eq!(first.intended_event_ids, second.intended_event_ids);
}

#[test]
fn cassette_failed_boundary_remains_unintended() {
    let intended = IntendedEdit {
        edit_id: "cassette".into(),
        seq_id: "chr".into(),
        start: 500,
        end: 500,
        ref_allele: ".".into(),
        alt: ".".into(),
        kind: "cassette".into(),
    };
    let result = run(
        vec![
            GdEntry::jc(1, "chr", 506, "+", "chr", 800, "-", 0),
            GdEntry::jc(2, "chr", 500, "-", "chr", 900, "+", 0),
        ],
        vec![],
        vec![intended],
    );
    let assessment = &result.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(assessment.status, IntendedEditStatus::Partial);
    assert!(assessment
        .event_relationships
        .iter()
        .any(|relation| relation.role == IntendedEventRole::PartialObservation));
    assert_eq!(result.unintended.len(), 1);
}

#[test]
fn deletion_junction_with_extra_variant_is_unexpected_structure() {
    let intended = IntendedEdit {
        edit_id: "deletion".into(),
        seq_id: "chr".into(),
        start: 100,
        end: 109,
        ref_allele: ".".into(),
        alt: ".".into(),
        kind: "del".into(),
    };
    let result = run(
        vec![
            GdEntry::jc(1, "chr", 100, "+", "chr", 109, "-", 0),
            GdEntry::snp(2, "chr", 105, "T"),
        ],
        vec![],
        vec![intended],
    );
    let assessment = &result.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(assessment.status, IntendedEditStatus::UnexpectedStructure);
    assert_eq!(assessment.unexpected_event_ids.len(), 1);
    assert!(result
        .unintended
        .iter()
        .any(|mutation| { mutation.event_id.as_ref() == assessment.unexpected_event_ids.first() }));
}

#[test]
fn mc_is_diagnostic_only_and_missing_evidence_stays_unknown() {
    let result = run(
        vec![GdEntry::mc(90, "chr", 90, 110, 0, 0)],
        vec![],
        vec![edit()],
    );
    let assessment = &result.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(assessment.status, IntendedEditStatus::Missing);
    assert!(assessment.matched_event_ids.is_empty());
    assert!(assessment.unexpected_event_ids.is_empty());
    assert_eq!(assessment.mc_observation, EvidenceObservation::Observed);
    assert_eq!(assessment.mc_diagnostics[0].gd_id, 90);
    assert!(assessment.notes[0].contains("differential mutation unconfirmed"));
    assert!(result.differential_events.is_empty());

    let no_mc = run(vec![], vec![], vec![edit()]);
    assert_eq!(
        no_mc.intended_edit_assessments.as_ref().expect("declared")[0].mc_observation,
        EvidenceObservation::Unknown
    );
}

#[test]
fn inherited_mc_does_not_create_edited_diagnostic() {
    let edited_mc = GdEntry::mc(90, "chr", 90, 110, 0, 0);
    let starter_mc = GdEntry::mc(3, "chr", 90, 110, 0, 0);
    let result = run(vec![edited_mc], vec![starter_mc], vec![edit()]);
    let assessment = &result.intended_edit_assessments.as_ref().expect("declared")[0];
    assert_eq!(assessment.mc_observation, EvidenceObservation::Unknown);
    assert!(assessment.mc_diagnostics.is_empty());
}

#[test]
fn event_id_and_evidence_reach_audit_without_offtarget_rewrite() {
    let mut snp = GdEntry::snp(7, "chr", 100, "T");
    snp.parent_ids.push(50);
    let ra = GdEntry::ra(50, "chr", 100, 0, "A", "T");
    let result = run(vec![ra, snp], vec![], vec![edit()]);
    let event_id = result.differential_events[0].event_id.clone();
    assert_eq!(
        result.differential_events[0]
            .evidence
            .ra
            .as_ref()
            .expect("RA reference")[0]
            .gd_id,
        50
    );
    let audit = build_audit_result(
        SampleMetadata {
            starter_names: vec![],
            edited_names: vec![],
            reference_names: vec![],
            editor: "dsb".into(),
            spacer: None,
            pam: None,
            threads: 1,
        },
        result.intended_edit_assessments.clone().expect("declared"),
        &result.unintended,
        &result.intended_event_ids,
        &result.differential_events,
        &[],
        &[],
        AnalysisProvenance {
            prokadiff_version: "test".into(),
            git_commit: "test".into(),
            reference_sha256: None,
            bowtie2_version: None,
            offtarget_search_status: "NOT_REQUESTED".into(),
            cfd_scoring_status: "DISABLED".into(),
            hsu_scoring_status: "DISABLED".into(),
            bulge_search_status: "DISABLED".into(),
            run_timestamp: "test".into(),
        },
        &[],
    );
    assert_eq!(audit.intended_edits[0].matched_event_ids, vec![event_id]);
    assert_eq!(audit.variants[0].origin_status, OriginStatus::Intended);
    assert_eq!(audit.variants[0].evidence.ra, Some(true));
    assert_eq!(audit.variants[0].evidence.mc, None);
    assert_eq!(
        audit.variants[0].evidence.format_brief(),
        "RA=1;MC=NA;JC=NA"
    );
}
