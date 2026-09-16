use prokadiff_classify::{
    AnalysisProvenance, AnnotatedVariant, AuditResult, BoundaryAssessment, EvidenceSummary,
    GeneAnnotation, GuideRelation, IntendedEditAssessment, IntendedEditStatus, IntendedRelation,
    MobileElementAnnotation, OriginStatus, RefContig, ReviewPriority, SampleMetadata, SizeClass,
};
use prokadiff_gd::GdEntry;
use prokadiff_report::{
    write_edit_outcomes_tsv, write_markdown_report, write_post_edit_variants_tsv,
    write_provenance_tsv,
};

fn mock_sample() -> SampleMetadata {
    SampleMetadata {
        starter_names: vec!["starter_R1.fq.gz".into(), "starter_R2.fq.gz".into()],
        edited_names: vec!["edited_R1.fq.gz".into(), "edited_R2.fq.gz".into()],
        reference_names: vec!["NC_000913.3.gbk".into()],
        editor: "cas9".into(),
        spacer: Some("GAGTTCATCTACGCCGTGAA".into()),
        pam: Some("NGG".into()),
        threads: 8,
    }
}

fn mock_provenance(cfd_status: &str) -> AnalysisProvenance {
    AnalysisProvenance {
        prokadiff_version: "0.2.0".into(),
        git_commit: "9f8e7d6c5b4a".into(),
        reference_sha256: Some(
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855".into(),
        ),
        bowtie2_version: Some("2.5.4".into()),
        offtarget_search_status: "VALIDATED_EXACT_MATCH".into(),
        cfd_scoring_status: cfd_status.into(),
        hsu_scoring_status: "DISABLED_UNVALIDATED_ORACLE".into(),
        bulge_search_status: "EXACT_UNGAPPED".into(),
        run_timestamp: "2026-09-16T12:00:00Z".into(),
    }
}

fn mock_refs() -> Vec<RefContig> {
    vec![RefContig {
        name: "NC_000913.3".into(),
        seq: b"ATCGATCGATCG".to_vec(),
    }]
}

#[test]
fn test_fixture_1_clean_edit() {
    let tmp_report = std::env::temp_dir().join("fixture1_report.md");

    let intended = vec![IntendedEditAssessment {
        edit_id: "edit_del".into(),
        kind: "del".into(),
        seq_id: "NC_000913.3".into(),
        expected_start: 100000,
        expected_end: 101200,
        status: IntendedEditStatus::Complete,
        matched_event_ids: vec![1],
        left_boundary: Some(BoundaryAssessment {
            expected_pos: 100000,
            observed_pos: Some(100000),
            diff_bp: 0,
            passed: true,
        }),
        right_boundary: Some(BoundaryAssessment {
            expected_pos: 101200,
            observed_pos: Some(101200),
            diff_bp: 0,
            passed: true,
        }),
        expected_size: Some(1201),
        observed_size: Some(1201),
        unexpected_event_ids: vec![],
        notes: vec!["Both junctions confirmed".into()],
    }];

    // 0 additional unintended variants, 1 intended
    let intended_var = AnnotatedVariant {
        variant_id: "VAR_0001".into(),
        entry: GdEntry::del(1, "NC_000913.3", 100000, 1201),
        origin_status: OriginStatus::Intended,
        size_class: SizeClass::Structural,
        intended_relation: IntendedRelation::Expected,
        guide_relation: GuideRelation::OnTarget,
        mobile_element_relation: None,
        repeat_relation: None,
        gene_annotation: None,
        evidence: EvidenceSummary {
            ra: false,
            mc: true,
            jc: true,
            supporting_reads: Some(45),
            coverage: Some(60.0),
        },
        review_priority: ReviewPriority::Info,
        legacy_class: None,
    };

    let audit = AuditResult {
        sample: mock_sample(),
        intended_edits: intended,
        variants: vec![intended_var],
        guide_sites: vec![],
        variant_site_links: vec![],
        provenance: mock_provenance("DISABLED_UNVALIDATED_ORACLE"),
    };

    write_markdown_report(&tmp_report, &audit, &mock_refs()).unwrap();
    let content = std::fs::read_to_string(&tmp_report).unwrap();

    assert!(content.contains("**Intended Edit Outcome: PASS**"));
    assert!(content.contains(
        "No additional post-edit differential variants were detected outside declared edits"
    ));
    assert!(content.contains("## 1. Executive Summary"));
    assert!(content.contains("## 8. Analysis Limitations & Cautions"));

    let tmp_outcomes = std::env::temp_dir().join("fixture1_outcomes.tsv");
    write_edit_outcomes_tsv(&tmp_outcomes, &audit.intended_edits).unwrap();
    let outcomes_content = std::fs::read_to_string(&tmp_outcomes).unwrap();
    assert!(outcomes_content.contains("edit_del\tdel\tNC_000913.3\t100000\t101200\tCOMPLETE"));

    let _ = std::fs::remove_file(tmp_report);
    let _ = std::fs::remove_file(tmp_outcomes);
}

#[test]
fn test_fixture_2_partial_cassette() {
    let tmp_report = std::env::temp_dir().join("fixture2_report.md");

    let intended = vec![IntendedEditAssessment {
        edit_id: "cassette_int".into(),
        kind: "cassette".into(),
        seq_id: "NC_000913.3".into(),
        expected_start: 500000,
        expected_end: 500000,
        status: IntendedEditStatus::Partial,
        matched_event_ids: vec![10],
        left_boundary: Some(BoundaryAssessment {
            expected_pos: 500000,
            observed_pos: Some(500000),
            diff_bp: 0,
            passed: true,
        }),
        right_boundary: None,
        expected_size: Some(3500),
        observed_size: None,
        unexpected_event_ids: vec![],
        notes: vec!["single junction detected (partial integration)".into()],
    }];

    let audit = AuditResult {
        sample: mock_sample(),
        intended_edits: intended,
        variants: vec![],
        guide_sites: vec![],
        variant_site_links: vec![],
        provenance: mock_provenance("DISABLED_UNVALIDATED_ORACLE"),
    };

    write_markdown_report(&tmp_report, &audit, &mock_refs()).unwrap();
    let content = std::fs::read_to_string(&tmp_report).unwrap();

    assert!(content.contains("**Intended Edit Outcome: PARTIAL**"));
    assert!(!content.contains("Intended Edit Outcome: PASS"));
    assert!(content.contains("single junction detected (partial integration)"));

    let _ = std::fs::remove_file(tmp_report);
}

#[test]
fn test_fixture_3_mobile_element_insertion() {
    let tmp_report = std::env::temp_dir().join("fixture3_report.md");

    let is_var = AnnotatedVariant {
        variant_id: "VAR_0001".into(),
        entry: GdEntry::mob(5, "NC_000913.3", 3445210, "IS1", "+", 9),
        origin_status: OriginStatus::PostEditDifferential,
        size_class: SizeClass::Structural,
        intended_relation: IntendedRelation::None,
        guide_relation: GuideRelation::None,
        mobile_element_relation: Some(MobileElementAnnotation {
            family: Some("IS1".into()),
            element_name: Some("IS1".into()),
            insertion_site: Some(3445210),
            target_site_duplication: Some("9bp".into()),
            source_copy: None,
        }),
        repeat_relation: None,
        gene_annotation: None,
        evidence: EvidenceSummary {
            ra: false,
            mc: false,
            jc: true,
            supporting_reads: Some(30),
            coverage: Some(50.0),
        },
        review_priority: ReviewPriority::HighAttention,
        legacy_class: None,
    };

    let audit = AuditResult {
        sample: mock_sample(),
        intended_edits: vec![],
        variants: vec![is_var],
        guide_sites: vec![],
        variant_site_links: vec![],
        provenance: mock_provenance("DISABLED_UNVALIDATED_ORACLE"),
    };

    write_markdown_report(&tmp_report, &audit, &mock_refs()).unwrap();
    let content = std::fs::read_to_string(&tmp_report).unwrap();

    assert!(content.contains("Mobile element IS1 insertion (TSD=9bp)"));
    assert!(content.contains("**HIGH_ATTENTION**"));
    assert!(content.contains("3445210"));

    let _ = std::fs::remove_file(tmp_report);
}

#[test]
fn test_fixture_4_candidate_offtarget() {
    let tmp_report = std::env::temp_dir().join("fixture4_report.md");

    let ot_var = AnnotatedVariant {
        variant_id: "VAR_0002".into(),
        entry: GdEntry::snp(8, "NC_000913.3", 2135820, "G"),
        origin_status: OriginStatus::PostEditDifferential,
        size_class: SizeClass::Small,
        intended_relation: IntendedRelation::None,
        guide_relation: GuideRelation::CandidateOffTarget {
            site_id: "SITE_0001".into(),
            spacer_mismatches: 2,
            pam: "AGG".into(),
            distance_to_site: 3,
        },
        mobile_element_relation: None,
        repeat_relation: None,
        gene_annotation: None,
        evidence: EvidenceSummary {
            ra: true,
            mc: false,
            jc: false,
            supporting_reads: Some(40),
            coverage: Some(55.0),
        },
        review_priority: ReviewPriority::Review,
        legacy_class: None,
    };

    let audit = AuditResult {
        sample: mock_sample(),
        intended_edits: vec![],
        variants: vec![ot_var],
        guide_sites: vec![],
        variant_site_links: vec![],
        provenance: mock_provenance("DISABLED_UNVALIDATED_ORACLE"),
    };

    write_markdown_report(&tmp_report, &audit, &mock_refs()).unwrap();
    let content = std::fs::read_to_string(&tmp_report).unwrap();

    assert!(content.contains("## 6. Candidate Guide-dependent Off-target Events"));
    assert!(content.contains("VAR_0002"));
    assert!(content.contains("AGG"));
    assert!(content.contains("3 bp"));
    assert!(content.contains("Spatial proximity and sequence homology to predicted guide sites indicate candidate association, but do not establish nuclease cleavage causality without experimental controls."));

    let _ = std::fs::remove_file(tmp_report);
}

#[test]
fn test_fixture_5_unvalidated_cfd_gate() {
    let tmp_report = std::env::temp_dir().join("fixture5_report.md");

    let audit = AuditResult {
        sample: mock_sample(),
        intended_edits: vec![],
        variants: vec![],
        guide_sites: vec![],
        variant_site_links: vec![],
        provenance: mock_provenance("DISABLED_UNVALIDATED_ORACLE"),
    };

    write_markdown_report(&tmp_report, &audit, &mock_refs()).unwrap();
    let content = std::fs::read_to_string(&tmp_report).unwrap();

    assert!(content.contains("CFD Scoring Unavailable"));
    assert!(content.contains("external-oracle parity validation was not enabled"));

    let _ = std::fs::remove_file(tmp_report);
}

#[test]
fn test_post_edit_variants_and_provenance_tsv() {
    let tmp_var_tsv = std::env::temp_dir().join("test_post_edit_variants.tsv");
    let tmp_prov_tsv = std::env::temp_dir().join("test_provenance.tsv");

    let ot_var = AnnotatedVariant {
        variant_id: "VAR_0001".into(),
        entry: GdEntry::snp(1, "NC_000913.3", 100, "T"),
        origin_status: OriginStatus::PostEditDifferential,
        size_class: SizeClass::Small,
        intended_relation: IntendedRelation::None,
        guide_relation: GuideRelation::CandidateOffTarget {
            site_id: "SITE_0001".into(),
            spacer_mismatches: 1,
            pam: "NGG".into(),
            distance_to_site: 5,
        },
        mobile_element_relation: None,
        repeat_relation: None,
        gene_annotation: None,
        evidence: EvidenceSummary {
            ra: true,
            mc: false,
            jc: false,
            supporting_reads: None,
            coverage: None,
        },
        review_priority: ReviewPriority::HighAttention,
        legacy_class: None,
    };

    write_post_edit_variants_tsv(&tmp_var_tsv, &[ot_var], &mock_refs()).unwrap();
    let content = std::fs::read_to_string(&tmp_var_tsv).unwrap();
    assert!(content.contains("variant_id\tseq_id\tposition"));
    assert!(content.contains("VAR_0001\tNC_000913.3\t100\t100\tSNP"));
    assert!(content
        .contains("POST_EDIT_DIFF\tNONE\tSMALL\tHIGH_ATTENTION\tCANDIDATE_OFF_TARGET\t1\tNGG\t5"));

    let prov = mock_provenance("DISABLED_UNVALIDATED_ORACLE");
    write_provenance_tsv(&tmp_prov_tsv, &prov).unwrap();
    let prov_content = std::fs::read_to_string(&tmp_prov_tsv).unwrap();
    assert!(prov_content.contains("prokadiff_version\t0.2.0"));
    assert!(prov_content.contains("cfd_scoring_status\tDISABLED_UNVALIDATED_ORACLE"));

    let _ = std::fs::remove_file(tmp_var_tsv);
    let _ = std::fs::remove_file(tmp_prov_tsv);
}

#[test]
fn test_gene_annotation_rendering_in_report_and_tsv() {
    let tmp_report = std::env::temp_dir().join("test_gene_report.md");
    let tmp_tsv = std::env::temp_dir().join("test_gene_post_edit.tsv");

    let gene_var = AnnotatedVariant {
        variant_id: "VAR_0001".into(),
        entry: GdEntry::snp(1, "NC_000913.3", 150, "C"),
        origin_status: OriginStatus::PostEditDifferential,
        size_class: SizeClass::Small,
        intended_relation: IntendedRelation::None,
        guide_relation: GuideRelation::None,
        mobile_element_relation: None,
        repeat_relation: None,
        gene_annotation: Some(GeneAnnotation {
            locus_tag: Some("b0001".into()),
            gene_name: Some("dnaA".into()),
            feature_type: "CDS".into(),
            product: Some("replication initiator".into()),
            consequence: Some("within CDS".into()),
        }),
        evidence: EvidenceSummary {
            ra: true,
            mc: false,
            jc: false,
            supporting_reads: Some(25),
            coverage: Some(40.0),
        },
        review_priority: ReviewPriority::Review,
        legacy_class: None,
    };

    let audit = AuditResult {
        sample: mock_sample(),
        intended_edits: vec![],
        variants: vec![gene_var.clone()],
        guide_sites: vec![],
        variant_site_links: vec![],
        provenance: mock_provenance("DISABLED_UNVALIDATED_ORACLE"),
    };

    write_markdown_report(&tmp_report, &audit, &mock_refs()).unwrap();
    let rep_content = std::fs::read_to_string(&tmp_report).unwrap();
    assert!(rep_content.contains("dnaA (b0001) [CDS]"));

    write_post_edit_variants_tsv(&tmp_tsv, &[gene_var], &mock_refs()).unwrap();
    let tsv_content = std::fs::read_to_string(&tmp_tsv).unwrap();
    assert!(tsv_content.contains("dnaA\tb0001\tCDS"));

    let _ = std::fs::remove_file(tmp_report);
    let _ = std::fs::remove_file(tmp_tsv);
}
