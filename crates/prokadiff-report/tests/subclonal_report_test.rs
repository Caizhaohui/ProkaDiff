use prokadiff_classify::{
    subclonal::{AlleleCounts, StrandCounts, SubclonalClass, SubclonalVariant, WilsonCi},
    AnalysisProvenance, AuditResult, AuditSummary, CandidateSearchStatus, SampleMetadata,
};
use prokadiff_report::{write_markdown_report, write_subclonal_audit_tsv, SchemaVersion, SUBCLONAL_CAVEAT_TEXT};

fn mock_empty_audit() -> AuditResult {
    AuditResult {
        sample: SampleMetadata {
            starter_names: vec!["s.fq".into()],
            edited_names: vec!["e.fq".into()],
            reference_names: vec!["ref.fa".into()],
            editor: "cas9".into(),
            spacer: None,
            pam: None,
            threads: 1,
        },
        intended_provided: false,
        intended_edits: vec![],
        intended_observed_event_ids: vec![],
        starter_vs_reference: 0,
        hypothesis_enabled: false,
        event_index: std::collections::BTreeMap::new(),
        event_display: std::collections::BTreeMap::new(),
        unintended: vec![],
        variants: vec![],
        guide_sites: vec![],
        variant_site_links: vec![],
        associations: vec![],
        candidate_search_status: CandidateSearchStatus::NotPerformed,
        summary: AuditSummary {
            intended_provided: false,
            intended_declared: None,
            intended_edits_complete: None,
            intended_edits_partial: None,
            intended_edits_missing: None,
            intended_edits_unexpected: None,
            intended_events_observed: None,
            intended_status: None,
            intended_missing: None,
            starter_vs_reference: 0,
            unintended_count: 0,
            post_edit_variant_count: 0,
            structural_count: 0,
            near_homolog_count: 0,
            scattered_snv_count: 0,
            mobile_element_count: 0,
            non_mobile_structural_count: 0,
            candidate_offtarget_count: 0,
            unassociated_small_count: 0,
            high_attention_count: 0,
            review_count: 0,
            info_count: 0,
            association_count: 0,
            candidate_search_status: CandidateSearchStatus::NotPerformed,
        },
        provenance: AnalysisProvenance {
            prokadiff_version: "0.2.0".into(),
            git_commit: "test_sha".into(),
            reference_sha256: None,
            bowtie2_version: None,
            offtarget_search_status: "OK".into(),
            cfd_scoring_status: "OK".into(),
            hsu_scoring_status: "OK".into(),
            bulge_search_status: "OK".into(),
            run_timestamp: "2026-11-05".into(),
        },
        subclonal_enabled: false,
        subclonal_variants: vec![],
    }
}

#[test]
fn test_default_mode_report_does_not_contain_subclonal_section() {
    let audit = mock_empty_audit();
    let report_path = std::env::temp_dir().join("test_default_subclonal_report.md");
    write_markdown_report(&report_path, &audit, SchemaVersion::V2).unwrap();

    let content = std::fs::read_to_string(&report_path).unwrap();
    let _ = std::fs::remove_file(report_path);
    assert!(!content.contains("Subclonal / Pool-seq Variant Audit"));
    assert!(content.contains("## 8. Analysis Limitations & Cautions"));
}

#[test]
fn test_subclonal_mode_report_contains_verbatim_caveat_and_table() {
    let mut audit = mock_empty_audit();
    let subclonal_var = SubclonalVariant {
        contig: "NC_000913.3".to_string(),
        pos: 12345,
        ref_seq: "C".to_string(),
        alt_seq: "T".to_string(),
        variant_type: "SNP".to_string(),
        starter_counts: AlleleCounts::new(StrandCounts::new(50, 50), StrandCounts::new(0, 0), 0),
        edited_counts: AlleleCounts::new(StrandCounts::new(45, 45), StrandCounts::new(5, 5), 0),
        starter_af: Some(0.0),
        edited_af: Some(0.10),
        starter_ci: Some(WilsonCi { low: 0.0, high: 0.036 }),
        edited_ci: Some(WilsonCi { low: 0.055, high: 0.174 }),
        starter_strand_bias_p: 1.0,
        edited_strand_bias_p: 0.5,
        starter_error_rate: 0.005,
        edited_error_rate: 0.005,
        starter_p_err: 1.0,
        starter_q_err: 1.0,
        edited_p_err: 1e-4,
        edited_q_err: 1e-4,
        delta_af: 0.10,
        delta_ci_low: 0.02,
        delta_ci_high: 0.18,
        p_diff: 1e-4,
        q_diff: 1e-4,
        starter_exclusion_bound: Some(0.036),
        classification: SubclonalClass::EditedSpecific,
        flags: vec!["NONE".to_string()],
    };
    audit.subclonal_variants.push(subclonal_var);

    let report_path = std::env::temp_dir().join("test_subclonal_report.md");
    let tsv_path = std::env::temp_dir().join("test_subclonal_audit.tsv");

    write_markdown_report(&report_path, &audit, SchemaVersion::V2).unwrap();
    write_subclonal_audit_tsv(&tsv_path, &audit.subclonal_variants).unwrap();

    let content = std::fs::read_to_string(&report_path).unwrap();
    let _ = std::fs::remove_file(report_path);
    assert!(content.contains("## Subclonal / Pool-seq Variant Audit (Target AF 0.05–0.70)"));
    assert!(content.contains(SUBCLONAL_CAVEAT_TEXT));
    assert!(content.contains("NC_000913.3"));
    assert!(content.contains("`EDITED_SPECIFIC`"));
    assert!(content.contains("## 8. Analysis Limitations & Cautions"));

    let tsv_content = std::fs::read_to_string(&tsv_path).unwrap();
    let _ = std::fs::remove_file(tsv_path);
    let header = tsv_content.lines().next().unwrap();
    assert!(header.starts_with("contig\tposition\tvariant_type\tref\talt\tstarter_depth"));
    assert!(tsv_content.contains("EDITED_SPECIFIC"));
}

#[test]
fn test_mode_switch_cleans_up_stale_subclonal_tsv_with_rollback() {
    let outdir = std::env::temp_dir().join(format!("test_mode_switch_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&outdir);
    std::fs::create_dir_all(&outdir).unwrap();

    // 1. First run: Subclonal mode enabled with a variant -> writes subclonal_audit.tsv
    let mut subclonal_audit = mock_empty_audit();
    subclonal_audit.subclonal_enabled = true;
    subclonal_audit.subclonal_variants.push(SubclonalVariant {
        contig: "NC_000913.3".into(),
        pos: 100,
        ref_seq: "A".into(),
        alt_seq: "C".into(),
        variant_type: "SNP".into(),
        starter_counts: AlleleCounts::new(StrandCounts::new(50, 50), StrandCounts::new(0, 0), 0),
        edited_counts: AlleleCounts::new(StrandCounts::new(45, 45), StrandCounts::new(5, 5), 0),
        starter_af: Some(0.0),
        edited_af: Some(0.10),
        starter_ci: None,
        edited_ci: None,
        starter_strand_bias_p: 1.0,
        edited_strand_bias_p: 1.0,
        starter_error_rate: 0.005,
        edited_error_rate: 0.005,
        starter_p_err: 1.0,
        starter_q_err: 1.0,
        edited_p_err: 1e-4,
        edited_q_err: 1e-4,
        delta_af: 0.10,
        delta_ci_low: 0.02,
        delta_ci_high: 0.18,
        p_diff: 1e-4,
        q_diff: 1e-4,
        starter_exclusion_bound: Some(0.036),
        classification: SubclonalClass::EditedSpecific,
        flags: vec![],
    });

    prokadiff_report::write_product_outputs(&subclonal_audit, SchemaVersion::V2, &outdir).unwrap();
    assert!(outdir.join("subclonal_audit.tsv").is_file());

    // 2. Second run: Default mode (subclonal disabled) into SAME directory -> subclonal_audit.tsv must be cleaned up
    let default_audit = mock_empty_audit();
    prokadiff_report::write_product_outputs(&default_audit, SchemaVersion::V2, &outdir).unwrap();
    assert!(
        !outdir.join("subclonal_audit.tsv").exists(),
        "stale subclonal_audit.tsv must be removed on switch back to default mode"
    );

    let _ = std::fs::remove_dir_all(&outdir);
}
