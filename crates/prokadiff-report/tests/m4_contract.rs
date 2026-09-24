use prokadiff_classify::{
    build_complete_audit_result, AnalysisProvenance, AuditResult, CandidateSearchStatus,
    CanonicalEvent, ClassifiedMutation, DifferentialEvent, EventId, EvidenceObservation,
    EvidenceReferences, IntendedEditAssessment, IntendedEditStatus, IntendedEventRelationship,
    IntendedEventRole, MutationClass, MutationSiteAssociation, RefContig, SampleMetadata,
};
use prokadiff_gd::GdEntry;
use prokadiff_offtarget::{BulgeType, OffTargetSite, Strand};
use prokadiff_report::{write_product_outputs, SchemaVersion};

fn refs() -> Vec<RefContig> {
    vec![RefContig {
        name: "chr".into(),
        seq: vec![b'A'; 1000],
    }]
}

fn sample() -> SampleMetadata {
    SampleMetadata {
        starter_names: vec!["starter.fq".into()],
        edited_names: vec!["edited.fq".into()],
        reference_names: vec!["ref.fa".into()],
        editor: "cas9".into(),
        spacer: Some("ACGT".into()),
        pam: Some("NGG".into()),
        threads: 1,
    }
}

fn provenance() -> AnalysisProvenance {
    AnalysisProvenance {
        prokadiff_version: "test".into(),
        git_commit: "test".into(),
        reference_sha256: None,
        bowtie2_version: None,
        offtarget_search_status: "FIXTURE_VALIDATED".into(),
        cfd_scoring_status: "DISABLED".into(),
        hsu_scoring_status: "DISABLED".into(),
        bulge_search_status: "EXACT_UNGAPPED".into(),
        run_timestamp: "2026-09-24T00:00:00Z".into(),
    }
}

fn event(entry: GdEntry, refs: &[RefContig]) -> DifferentialEvent {
    let canonical = CanonicalEvent::from_gd_entry(&entry, refs).expect("valid test event");
    let (event_id, _) = canonical.compute_event_id();
    DifferentialEvent {
        event_id,
        canonical,
        representative: entry,
        merged_source_ids: vec![],
        evidence: EvidenceReferences::default(),
    }
}

fn mutation(
    event: &DifferentialEvent,
    class: MutationClass,
    hypothesis: Option<&str>,
) -> ClassifiedMutation {
    ClassifiedMutation {
        entry: event.representative.clone(),
        class,
        pam_profile: None,
        offtarget_mismatch: None,
        distance_to_site: None,
        hypothesis: hypothesis.map(str::to_string),
        event_id: Some(event.event_id.clone()),
    }
}

fn site(id: &str, start: u64, end: u64) -> OffTargetSite {
    OffTargetSite {
        site_id: id.into(),
        seq_id: "chr".into(),
        start,
        end,
        strand: Strand::Plus,
        guide: "ACGT".into(),
        target_seq: "ACGTAGG".into(),
        pam: "AGG".into(),
        mismatches: 1,
        bulge_type: BulgeType::Dna,
        bulge_size: 1,
        search_backend: "rust_bulge".into(),
        cfd_score: Some(0.5),
        hsu_score: None,
    }
}

fn association(
    event: &DifferentialEvent,
    site: &OffTargetSite,
    side: Option<prokadiff_classify::JunctionSideTag>,
) -> MutationSiteAssociation {
    MutationSiteAssociation {
        event_id: event.event_id.clone(),
        mutation_type: event.representative.kind.as_str().into(),
        mutation_position: event.representative.position().expect("position"),
        junction_side: side,
        site_id: site.site_id.clone(),
        site_seq_id: site.seq_id.clone(),
        site_start: site.start,
        site_end: site.end,
        site_strand: site.strand,
        distance_to_site: 3,
        mismatches: site.mismatches,
        pam: site.pam.clone(),
        cfd_score: site.cfd_score,
        hsu_score: site.hsu_score,
        search_backend: site.search_backend.clone(),
        bulge_type: site.bulge_type,
        bulge_size: site.bulge_size,
        target_seq: site.target_seq.clone(),
        association_window: 50,
    }
}

fn assessment(
    edit_id: &str,
    expected: &EventId,
    unexpected: Option<&EventId>,
) -> IntendedEditAssessment {
    let mut relationships = vec![IntendedEventRelationship {
        event_id: expected.clone(),
        role: IntendedEventRole::ExpectedConstituent,
    }];
    if let Some(unexpected) = unexpected {
        relationships.push(IntendedEventRelationship {
            event_id: unexpected.clone(),
            role: IntendedEventRole::UnexpectedAtLocus,
        });
    }
    IntendedEditAssessment {
        edit_id: edit_id.into(),
        kind: "snp".into(),
        seq_id: "chr".into(),
        expected_start: 10,
        expected_end: 10,
        status: if unexpected.is_some() {
            IntendedEditStatus::UnexpectedStructure
        } else {
            IntendedEditStatus::Complete
        },
        matched_event_ids: vec![expected.clone()],
        event_relationships: relationships,
        left_boundary: None,
        right_boundary: None,
        expected_size: None,
        observed_size: None,
        unexpected_event_ids: unexpected.into_iter().cloned().collect(),
        mc_diagnostics: vec![],
        mc_observation: EvidenceObservation::Unknown,
        notes: vec![],
    }
}

#[allow(clippy::too_many_arguments)]
fn audit(
    intended_provided: bool,
    assessments: Vec<IntendedEditAssessment>,
    intended_events: Vec<EventId>,
    unintended: Vec<ClassifiedMutation>,
    events: Vec<DifferentialEvent>,
    sites: Vec<OffTargetSite>,
    associations: Vec<MutationSiteAssociation>,
    hypothesis_enabled: bool,
    status: CandidateSearchStatus,
) -> AuditResult {
    build_complete_audit_result(
        sample(),
        intended_provided,
        assessments,
        intended_events,
        7,
        hypothesis_enabled,
        &unintended,
        &events,
        &sites,
        &associations,
        status,
        provenance(),
        &[],
        &refs(),
    )
    .expect("valid aggregate")
}

fn outdir(name: &str) -> std::path::PathBuf {
    let path = std::env::temp_dir().join(format!("prokadiff_m4_{name}_{}", std::process::id()));
    std::fs::create_dir_all(&path).expect("create test output directory");
    path
}

const PRODUCT_FILES: [&str; 7] = [
    "summary.txt",
    "unintended.tsv",
    "edit_outcomes.tsv",
    "post_edit_variants.tsv",
    "offtarget_sites.tsv",
    "mutation_offtarget_links.tsv",
    "report.md",
];

fn m3_v1_golden(case: &str, file: &str) -> &'static str {
    match (case, file) {
        ("unexpected_structure", "summary.txt") => {
            include_str!("golden/m3_v1/unexpected_structure/summary.txt")
        }
        ("unexpected_structure", "unintended.tsv") => {
            include_str!("golden/m3_v1/unexpected_structure/unintended.tsv")
        }
        ("unexpected_structure", "edit_outcomes.tsv") => {
            include_str!("golden/m3_v1/unexpected_structure/edit_outcomes.tsv")
        }
        ("unexpected_structure", "post_edit_variants.tsv") => {
            include_str!("golden/m3_v1/unexpected_structure/post_edit_variants.tsv")
        }
        ("unexpected_structure", "offtarget_sites.tsv") => {
            include_str!("golden/m3_v1/unexpected_structure/offtarget_sites.tsv")
        }
        ("unexpected_structure", "mutation_offtarget_links.tsv") => {
            include_str!("golden/m3_v1/unexpected_structure/mutation_offtarget_links.tsv")
        }
        ("unexpected_structure", "report.md") => {
            include_str!("golden/m3_v1/unexpected_structure/report.md")
        }
        ("zero_row_intended", "summary.txt") => {
            include_str!("golden/m3_v1/zero_row_intended/summary.txt")
        }
        ("zero_row_intended", "unintended.tsv") => {
            include_str!("golden/m3_v1/zero_row_intended/unintended.tsv")
        }
        ("zero_row_intended", "edit_outcomes.tsv") => {
            include_str!("golden/m3_v1/zero_row_intended/edit_outcomes.tsv")
        }
        ("zero_row_intended", "post_edit_variants.tsv") => {
            include_str!("golden/m3_v1/zero_row_intended/post_edit_variants.tsv")
        }
        ("zero_row_intended", "offtarget_sites.tsv") => {
            include_str!("golden/m3_v1/zero_row_intended/offtarget_sites.tsv")
        }
        ("zero_row_intended", "mutation_offtarget_links.tsv") => {
            include_str!("golden/m3_v1/zero_row_intended/mutation_offtarget_links.tsv")
        }
        ("zero_row_intended", "report.md") => {
            include_str!("golden/m3_v1/zero_row_intended/report.md")
        }
        _ => panic!("unknown M3 golden fixture {case}/{file}"),
    }
}

fn compare_v1_to_m3_golden(case: &str, aggregate: &AuditResult) {
    let path = outdir(&format!("m3_golden_{case}"));
    write_product_outputs(aggregate, SchemaVersion::V1, &path).expect("render explicit V1 set");
    for file in PRODUCT_FILES {
        let actual = std::fs::read_to_string(path.join(file)).expect("read rendered product");
        assert_eq!(
            actual,
            m3_v1_golden(case, file),
            "V1 differs from M3 {file}"
        );
    }
}

fn event_ids_in(text: &str) -> std::collections::BTreeSet<String> {
    text.split(|character: char| !character.is_ascii_alphanumeric() && character != '_')
        .filter(|token| token.starts_with("DV1_") && token.len() == 36)
        .map(str::to_string)
        .collect()
}

fn tsv_event_ids(path: &std::path::Path) -> std::collections::BTreeSet<String> {
    let text = std::fs::read_to_string(path).expect("read event-id TSV");
    let mut lines = text.lines();
    let header: Vec<_> = lines.next().expect("TSV header").split('\t').collect();
    let Some(index) = header.iter().position(|column| *column == "event_id") else {
        return std::collections::BTreeSet::new();
    };
    lines
        .filter_map(|line| line.split('\t').nth(index))
        .filter(|event_id| event_id.starts_with("DV1_"))
        .map(str::to_string)
        .collect()
}

#[test]
fn schema_versions_preserve_v1_and_append_v2_with_hypothesis_modes() {
    let refs = refs();
    let differential_event = event(GdEntry::snp(1, "chr", 10, "T"), &refs);
    for (version, hypothesis_enabled) in [
        (SchemaVersion::V1, false),
        (SchemaVersion::V1, true),
        (SchemaVersion::V2, false),
        (SchemaVersion::V2, true),
    ] {
        let audit = audit(
            false,
            vec![],
            vec![],
            vec![mutation(
                &differential_event,
                MutationClass::ScatteredSnv,
                Some("sos_widney2014"),
            )],
            vec![differential_event.clone()],
            vec![],
            vec![],
            hypothesis_enabled,
            CandidateSearchStatus::NotPerformed,
        );
        let path = outdir(&format!("schema_{version:?}_{hypothesis_enabled}"));
        write_product_outputs(&audit, version, &path).expect("render products");
        let header = std::fs::read_to_string(path.join("unintended.tsv"))
            .expect("unintended output")
            .lines()
            .next()
            .expect("header")
            .to_string();
        let expected_v1 = if hypothesis_enabled {
            "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\tside2_seq_id\tside2_position\thypothesis"
        } else {
            "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\tside2_seq_id\tside2_position"
        };
        match version {
            SchemaVersion::V1 => assert_eq!(header, expected_v1),
            SchemaVersion::V2 => assert_eq!(header, format!("{expected_v1}\tevent_id")),
        }
        let post_header = std::fs::read_to_string(path.join("post_edit_variants.tsv"))
            .expect("post output")
            .lines()
            .next()
            .expect("header")
            .to_string();
        let expected_post = "variant_id\tseq_id\tposition\tend\tgd_type\tref\talt\torigin_status\tintended_relation\tsize_class\treview_priority\tguide_relation\tguide_mismatches\tpam\tdist_to_site\tmobile_element\tevidence\tgene\tlocus_tag\tfeature_type\tlegacy_class";
        match version {
            SchemaVersion::V1 => assert_eq!(post_header, expected_post),
            SchemaVersion::V2 => assert_eq!(post_header, format!("{expected_post}\tevent_id")),
        }
        let outcomes_header = std::fs::read_to_string(path.join("edit_outcomes.tsv"))
            .expect("outcomes output")
            .lines()
            .next()
            .expect("header")
            .to_string();
        let links_header = std::fs::read_to_string(path.join("mutation_offtarget_links.tsv"))
            .expect("links output")
            .lines()
            .next()
            .expect("header")
            .to_string();
        match version {
            SchemaVersion::V1 => {
                assert_eq!(
                    outcomes_header,
                    "edit_id\tkind\tseq_id\texpected_start\texpected_end\tstatus\tmatched_event_ids\tleft_boundary_status\tright_boundary_status\texpected_size\tobserved_size\tunexpected_events\tnotes"
                );
                assert_eq!(
                    links_header,
                    "mutation_id\tsite_id\tmutation_type\tmutation_position\tsite_start\tsite_end\tdistance_to_site\tmismatches\tpam\tcfd_score\tassociation_window"
                );
            }
            SchemaVersion::V2 => {
                assert_eq!(
                    outcomes_header,
                    "edit_id\tkind\tseq_id\texpected_start\texpected_end\tstatus\tmatched_event_ids\tleft_boundary_status\tright_boundary_status\texpected_size\tobserved_size\tunexpected_events\tnotes\tmatched_event_ids_dv1\tunexpected_event_ids_dv1"
                );
                assert_eq!(
                    links_header,
                    "mutation_id\tsite_id\tmutation_type\tmutation_position\tsite_start\tsite_end\tdistance_to_site\tmismatches\tpam\tcfd_score\tassociation_window\tevent_id\tjunction_side\tsearch_backend\tbulge_type\tbulge_size"
                );
            }
        }
        let report = std::fs::read_to_string(path.join("report.md")).expect("report output");
        if version.is_v2() {
            assert!(report.contains(differential_event.event_id.as_str()));
        } else {
            assert!(!report.contains(differential_event.event_id.as_str()));
        }
    }
}

#[test]
fn aggregate_distinguishes_intended_modes_and_retains_same_locus_relationships() {
    let refs = refs();
    let expected = event(GdEntry::snp(1, "chr", 10, "T"), &refs);
    let unexpected = event(GdEntry::snp(2, "chr", 10, "C"), &refs);
    let assessment_one = assessment("edit_one", &expected.event_id, Some(&unexpected.event_id));
    let assessment_two = assessment("edit_two", &expected.event_id, None);
    let aggregate = audit(
        true,
        vec![assessment_one, assessment_two],
        vec![expected.event_id.clone()],
        vec![mutation(&unexpected, MutationClass::ScatteredSnv, None)],
        vec![expected.clone(), unexpected.clone()],
        vec![],
        vec![],
        false,
        CandidateSearchStatus::PerformedNoCandidates,
    );
    let path = outdir("intended");
    write_product_outputs(&aggregate, SchemaVersion::V2, &path).expect("render products");
    let summary = std::fs::read_to_string(path.join("summary.txt")).expect("summary");
    let outcomes = std::fs::read_to_string(path.join("edit_outcomes.tsv")).expect("outcomes");
    let unintended = std::fs::read_to_string(path.join("unintended.tsv")).expect("unintended");
    assert!(summary.contains("intended_provided\tyes"));
    assert!(summary.contains("intended_edits_declared\t2"));
    assert!(outcomes.contains(expected.event_id.as_str()));
    assert!(outcomes.contains(unexpected.event_id.as_str()));
    assert!(unintended.contains(unexpected.event_id.as_str()));
    assert!(std::fs::read_to_string(path.join("report.md"))
        .expect("report")
        .contains(expected.event_id.as_str()));
    assert_eq!(aggregate.variants.len(), 2);

    let no_table = audit(
        false,
        vec![],
        vec![],
        vec![mutation(&expected, MutationClass::ScatteredSnv, None)],
        vec![expected.clone()],
        vec![],
        vec![],
        false,
        CandidateSearchStatus::PerformedNoCandidates,
    );
    let empty_table = audit(
        true,
        vec![],
        vec![],
        vec![mutation(&expected, MutationClass::ScatteredSnv, None)],
        vec![expected],
        vec![],
        vec![],
        false,
        CandidateSearchStatus::PerformedWithCandidates,
    );
    let no_table_path = outdir("no_table");
    let empty_table_path = outdir("empty_table");
    write_product_outputs(&no_table, SchemaVersion::V2, &no_table_path).expect("render no table");
    write_product_outputs(&empty_table, SchemaVersion::V2, &empty_table_path)
        .expect("render empty table");
    assert!(std::fs::read_to_string(no_table_path.join("summary.txt"))
        .expect("no-table summary")
        .contains("intended_provided\tno"));
    assert!(
        std::fs::read_to_string(empty_table_path.join("summary.txt"))
            .expect("empty-table summary")
            .contains("intended_edits_declared\t0")
    );
    assert!(std::fs::read_to_string(empty_table_path.join("report.md"))
        .expect("empty-table report")
        .contains("the supplied intended-edit table contains zero rows"));
    assert!(std::fs::read_to_string(no_table_path.join("report.md"))
        .expect("no-candidate report")
        .contains("performed search found no candidate guide sites"));
    assert!(std::fs::read_to_string(empty_table_path.join("report.md"))
        .expect("unlinked report")
        .contains("without nearby predicted guide-homologous sites"));
}

#[test]
fn v2_links_preserve_each_junction_side_and_site_evidence() {
    let refs = refs();
    let junction = event(GdEntry::jc(3, "chr", 20, "+", "chr", 30, "-", 0), &refs);
    let site_one = site("SITE_000001", 15, 18);
    let site_two = site("SITE_000002", 32, 35);
    let site_three = site("SITE_000003", 16, 19);
    let aggregate = audit(
        false,
        vec![],
        vec![],
        vec![mutation(&junction, MutationClass::Structural, None)],
        vec![junction.clone()],
        vec![site_one.clone(), site_two.clone(), site_three.clone()],
        vec![
            association(
                &junction,
                &site_one,
                Some(prokadiff_classify::JunctionSideTag::Side1),
            ),
            association(
                &junction,
                &site_two,
                Some(prokadiff_classify::JunctionSideTag::Side2),
            ),
            association(
                &junction,
                &site_three,
                Some(prokadiff_classify::JunctionSideTag::Side1),
            ),
        ],
        false,
        CandidateSearchStatus::PerformedWithCandidates,
    );
    let path = outdir("junction");
    write_product_outputs(&aggregate, SchemaVersion::V2, &path).expect("render products");
    let links = std::fs::read_to_string(path.join("mutation_offtarget_links.tsv")).expect("links");
    assert_eq!(links.lines().count(), 4);
    assert!(links.contains("side_1"));
    assert!(links.contains("side_2"));
    assert!(links.contains("rust_bulge\tDNA\t1"));
    assert!(links.contains(junction.event_id.as_str()));
    assert!(std::fs::read_to_string(path.join("report.md"))
        .expect("report")
        .contains(junction.event_id.as_str()));
}

#[test]
fn mob_companion_junction_resolves_without_creating_a_variant_row() {
    let refs = refs();
    let mobile_element = event(GdEntry::mob(4, "chr", 100, "IS1", "+", 3), &refs);
    let companion = event(GdEntry::jc(5, "chr", 100, "+", "chr", 200, "-", 0), &refs);
    let candidate_site = site("SITE_000010", 196, 199);
    let aggregate = audit(
        false,
        vec![],
        vec![],
        vec![mutation(&mobile_element, MutationClass::Structural, None)],
        vec![mobile_element.clone(), companion.clone()],
        vec![candidate_site.clone()],
        vec![association(
            &companion,
            &candidate_site,
            Some(prokadiff_classify::JunctionSideTag::Side2),
        )],
        false,
        CandidateSearchStatus::PerformedWithCandidates,
    );
    assert_eq!(aggregate.event_index.len(), 2);
    assert!(aggregate.event_index.contains_key(&companion.event_id));
    assert_eq!(aggregate.variants.len(), 1);
    let path = outdir("mob_companion");
    write_product_outputs(&aggregate, SchemaVersion::V2, &path).expect("render products");
    assert_eq!(
        std::fs::read_to_string(path.join("post_edit_variants.tsv"))
            .expect("post variants")
            .lines()
            .count(),
        2
    );
    assert_eq!(
        std::fs::read_to_string(path.join("unintended.tsv"))
            .expect("unintended")
            .lines()
            .count(),
        2
    );
    assert!(
        std::fs::read_to_string(path.join("mutation_offtarget_links.tsv"))
            .expect("links")
            .contains(companion.event_id.as_str())
    );
}

#[test]
fn one_audit_result_has_consistent_public_ids_counts_and_intended_status() {
    let refs = refs();
    let expected = event(GdEntry::snp(41, "chr", 10, "T"), &refs);
    let unexpected = event(GdEntry::snp(42, "chr", 20, "C"), &refs);
    let mobile_element = event(GdEntry::mob(43, "chr", 100, "IS1", "+", 3), &refs);
    let companion = event(GdEntry::jc(44, "chr", 100, "+", "chr", 200, "-", 0), &refs);
    let candidate_site = site("SITE_000041", 18, 22);
    let aggregate = audit(
        true,
        vec![assessment(
            "edit",
            &expected.event_id,
            Some(&unexpected.event_id),
        )],
        vec![expected.event_id.clone()],
        vec![
            mutation(&unexpected, MutationClass::ScatteredSnv, None),
            mutation(&mobile_element, MutationClass::Structural, None),
        ],
        vec![
            expected.clone(),
            unexpected.clone(),
            mobile_element.clone(),
            companion.clone(),
        ],
        vec![candidate_site.clone()],
        vec![
            association(&unexpected, &candidate_site, None),
            association(
                &companion,
                &candidate_site,
                Some(prokadiff_classify::JunctionSideTag::Side2),
            ),
        ],
        false,
        CandidateSearchStatus::PerformedWithCandidates,
    );
    let path = outdir("cross_output_consistency");
    write_product_outputs(&aggregate, SchemaVersion::V2, &path).expect("render V2 products");

    let post = std::fs::read_to_string(path.join("post_edit_variants.tsv")).expect("variants");
    let post_ids = tsv_event_ids(&path.join("post_edit_variants.tsv"));
    let unintended = std::fs::read_to_string(path.join("unintended.tsv")).expect("unintended");
    let unintended_ids = tsv_event_ids(&path.join("unintended.tsv"));
    let outcomes = std::fs::read_to_string(path.join("edit_outcomes.tsv")).expect("outcomes");
    let outcome_rows: Vec<Vec<&str>> = outcomes
        .lines()
        .map(|line| line.split('\t').collect())
        .collect();
    let header = &outcome_rows[0];
    let matched = header
        .iter()
        .position(|column| *column == "matched_event_ids_dv1")
        .expect("matched EventId column");
    let unexpected_column = header
        .iter()
        .position(|column| *column == "unexpected_event_ids_dv1")
        .expect("unexpected EventId column");
    let status_column = header
        .iter()
        .position(|column| *column == "status")
        .expect("status column");
    let matched_ids = event_ids_in(outcome_rows[1][matched]);
    let unexpected_ids_from_edit = event_ids_in(outcome_rows[1][unexpected_column]);
    assert_eq!(matched_ids, [expected.event_id.to_string()].into());
    assert_eq!(
        unexpected_ids_from_edit,
        [unexpected.event_id.to_string()].into()
    );
    assert_eq!(outcome_rows[1][status_column], "UNEXPECTED_STRUCTURE");

    let links = std::fs::read_to_string(path.join("mutation_offtarget_links.tsv")).expect("links");
    let link_ids = tsv_event_ids(&path.join("mutation_offtarget_links.tsv"));
    let report = std::fs::read_to_string(path.join("report.md")).expect("report");
    let report_ids = event_ids_in(&report);
    let summary = std::fs::read_to_string(path.join("summary.txt")).expect("summary");

    assert_eq!(
        post_ids,
        [
            expected.event_id.to_string(),
            unexpected.event_id.to_string(),
            mobile_element.event_id.to_string()
        ]
        .into()
    );
    assert_eq!(
        unintended_ids,
        [
            unexpected.event_id.to_string(),
            mobile_element.event_id.to_string()
        ]
        .into()
    );
    assert_eq!(
        link_ids,
        [
            unexpected.event_id.to_string(),
            companion.event_id.to_string()
        ]
        .into()
    );
    assert!(report_ids.is_superset(&matched_ids));
    assert!(report_ids.is_superset(&unexpected_ids_from_edit));
    assert!(report_ids.contains(mobile_element.event_id.as_str()));
    assert!(aggregate
        .associations
        .iter()
        .all(|association| aggregate.event_index.contains_key(&association.event_id)));
    assert!(aggregate.event_index.contains_key(&companion.event_id));
    assert!(!post_ids.contains(companion.event_id.as_str()));
    assert!(!unintended_ids.contains(companion.event_id.as_str()));

    let class_counts = unintended
        .lines()
        .skip(1)
        .filter_map(|line| line.split('\t').nth(6))
        .fold(
            std::collections::BTreeMap::<&str, usize>::new(),
            |mut counts, class| {
                *counts.entry(class).or_default() += 1;
                counts
            },
        );
    assert_eq!(
        class_counts.get("structural"),
        Some(&aggregate.summary.structural_count)
    );
    assert_eq!(
        class_counts.get("scattered_snv"),
        Some(&aggregate.summary.scattered_snv_count)
    );
    assert_eq!(
        aggregate.summary.unintended_count,
        unintended.lines().count() - 1
    );
    assert_eq!(
        aggregate.summary.post_edit_variant_count,
        post.lines().count() - 2
    );
    assert_eq!(
        aggregate.summary.association_count,
        links.lines().count() - 1
    );
    assert!(summary.contains(&format!(
        "structural\t{}",
        aggregate.summary.structural_count
    )));
    assert!(summary.contains(&format!(
        "scattered_snv\t{}",
        aggregate.summary.scattered_snv_count
    )));
    assert!(summary.contains("intended_status\tpartial"));
    assert!(report.contains("Intended Edit Outcome: ALERT"));
    assert!(report.contains(&format!(
        "Total post-edit differential variants:** {}",
        aggregate.summary.post_edit_variant_count
    )));
    assert!(report.contains("**UNEXPECTED_STRUCTURE**"));
}

#[test]
fn event_id_is_stable_across_raw_gd_renumbering_and_unknown_ids_fail() {
    let refs = refs();
    let first = event(GdEntry::snp(1, "chr", 44, "T"), &refs);
    let renumbered = event(GdEntry::snp(900, "chr", 44, "T"), &refs);
    let second = event(GdEntry::snp(2, "chr", 55, "G"), &refs);
    let second_renumbered = event(GdEntry::snp(800, "chr", 55, "G"), &refs);
    assert_eq!(first.event_id, renumbered.event_id);
    assert_eq!(second.event_id, second_renumbered.event_id);
    let third = event(GdEntry::snp(3, "chr", 66, "C"), &refs);
    let third_renumbered = event(GdEntry::snp(700, "chr", 66, "C"), &refs);
    let guide_site = site("SITE_000020", 43, 46);
    let original = audit(
        true,
        vec![assessment("edit", &first.event_id, Some(&second.event_id))],
        vec![first.event_id.clone()],
        vec![
            mutation(&second, MutationClass::ScatteredSnv, None),
            mutation(&third, MutationClass::ScatteredSnv, None),
        ],
        vec![first.clone(), second.clone(), third.clone()],
        vec![guide_site.clone()],
        vec![association(&second, &guide_site, None)],
        false,
        CandidateSearchStatus::PerformedWithCandidates,
    );
    let reordered_and_renumbered = audit(
        true,
        vec![assessment(
            "edit",
            &renumbered.event_id,
            Some(&second_renumbered.event_id),
        )],
        vec![renumbered.event_id.clone()],
        vec![
            mutation(&third_renumbered, MutationClass::ScatteredSnv, None),
            mutation(&second_renumbered, MutationClass::ScatteredSnv, None),
        ],
        vec![
            third_renumbered.clone(),
            second_renumbered.clone(),
            renumbered.clone(),
        ],
        vec![guide_site.clone()],
        vec![association(&second_renumbered, &guide_site, None)],
        false,
        CandidateSearchStatus::PerformedWithCandidates,
    );
    assert_eq!(
        original.event_index.keys().collect::<Vec<_>>(),
        reordered_and_renumbered
            .event_index
            .keys()
            .collect::<Vec<_>>()
    );
    let original_path = outdir("reorder_public_original");
    let reordered_path = outdir("reorder_public_renumbered");
    write_product_outputs(&original, SchemaVersion::V2, &original_path).expect("original V2");
    write_product_outputs(
        &reordered_and_renumbered,
        SchemaVersion::V2,
        &reordered_path,
    )
    .expect("reordered V2");
    for file in [
        "post_edit_variants.tsv",
        "unintended.tsv",
        "mutation_offtarget_links.tsv",
    ] {
        assert_eq!(
            tsv_event_ids(&original_path.join(file)),
            tsv_event_ids(&reordered_path.join(file)),
            "public EventId set changed in {file}"
        );
    }
    let outcome_event_ids = |path: &std::path::Path| {
        let text = std::fs::read_to_string(path.join("edit_outcomes.tsv")).expect("outcomes");
        let mut lines = text.lines();
        let header: Vec<_> = lines.next().expect("header").split('\t').collect();
        let matched = header
            .iter()
            .position(|column| *column == "matched_event_ids_dv1");
        let unexpected = header
            .iter()
            .position(|column| *column == "unexpected_event_ids_dv1");
        let row: Vec<_> = lines.next().expect("edit row").split('\t').collect();
        (
            row[matched.expect("matched DV1 column")].to_string(),
            row[unexpected.expect("unexpected DV1 column")].to_string(),
        )
    };
    assert_eq!(
        outcome_event_ids(&original_path),
        outcome_event_ids(&reordered_path)
    );
    let variant_labels = |path: &std::path::Path| {
        let text = std::fs::read_to_string(path.join("post_edit_variants.tsv")).expect("variants");
        text.lines()
            .skip(1)
            .filter_map(|line| {
                let columns: Vec<_> = line.split('\t').collect();
                columns
                    .last()
                    .filter(|event_id| event_id.starts_with("DV1_"))
                    .map(|event_id| ((*event_id).to_string(), columns[0].to_string()))
            })
            .collect::<std::collections::BTreeMap<_, _>>()
    };
    let original_labels = variant_labels(&original_path);
    let reordered_labels = variant_labels(&reordered_path);
    assert_eq!(original_labels.len(), 3);
    assert_eq!(reordered_labels.len(), 3);
    assert_ne!(
        original_labels[second.event_id.as_str()],
        reordered_labels[second.event_id.as_str()]
    );
    assert_ne!(
        original_labels[third.event_id.as_str()],
        reordered_labels[third.event_id.as_str()]
    );
    let legacy_outcome_cells = |path: &std::path::Path| {
        let text = std::fs::read_to_string(path.join("edit_outcomes.tsv")).expect("outcomes");
        let row: Vec<_> = text.lines().nth(1).expect("edit row").split('\t').collect();
        (row[6].to_string(), row[11].to_string())
    };
    assert_ne!(
        legacy_outcome_cells(&original_path),
        legacy_outcome_cells(&reordered_path)
    );
    let original_report = std::fs::read_to_string(original_path.join("report.md")).expect("report");
    let reordered_report =
        std::fs::read_to_string(reordered_path.join("report.md")).expect("report");
    assert_eq!(
        event_ids_in(&original_report),
        event_ids_in(&reordered_report)
    );
    assert!(event_ids_in(&original_report).contains(first.event_id.as_str()));
    assert!(event_ids_in(&original_report).contains(second.event_id.as_str()));
    assert!(event_ids_in(&original_report).contains(third.event_id.as_str()));
    let unresolved = EventId::parse("DV1_0123456789abcdef0123456789abcdef").expect("valid id");
    let result = build_complete_audit_result(
        sample(),
        false,
        vec![],
        vec![],
        0,
        false,
        &[ClassifiedMutation {
            entry: GdEntry::snp(7, "chr", 50, "G"),
            class: MutationClass::ScatteredSnv,
            pam_profile: None,
            offtarget_mismatch: None,
            distance_to_site: None,
            hypothesis: None,
            event_id: Some(unresolved.clone()),
        }],
        std::slice::from_ref(&first),
        &[],
        &[],
        CandidateSearchStatus::PerformedNoCandidates,
        provenance(),
        &[],
        &refs,
    );
    assert!(matches!(
        result,
        Err(prokadiff_classify::AssociationProjectionError::UnknownEventId(event_id)) if event_id == unresolved
    ));

    let unknown_site = site("SITE_999999", 10, 12);
    let result = build_complete_audit_result(
        sample(),
        false,
        vec![],
        vec![],
        0,
        false,
        &[],
        std::slice::from_ref(&renumbered),
        &[],
        &[association(&renumbered, &unknown_site, None)],
        CandidateSearchStatus::PerformedWithCandidates,
        provenance(),
        &[],
        &refs,
    );
    assert!(matches!(
        result,
        Err(prokadiff_classify::AssociationProjectionError::UnknownSiteId(site_id)) if site_id == "SITE_999999"
    ));
}

#[test]
fn explicit_v1_product_outputs_match_m3_golden_sets() {
    let refs = refs();
    let expected = event(GdEntry::snp(11, "chr", 10, "T"), &refs);
    let unexpected = event(GdEntry::jc(12, "chr", 20, "+", "chr", 30, "-", 0), &refs);
    let candidate_site = site("SITE_000001", 18, 22);
    let unexpected_structure = audit(
        true,
        vec![assessment(
            "edit_one",
            &expected.event_id,
            Some(&unexpected.event_id),
        )],
        vec![expected.event_id.clone()],
        vec![mutation(&unexpected, MutationClass::Structural, None)],
        vec![expected, unexpected.clone()],
        vec![candidate_site.clone()],
        vec![association(&unexpected, &candidate_site, None)],
        false,
        CandidateSearchStatus::PerformedWithCandidates,
    );
    compare_v1_to_m3_golden("unexpected_structure", &unexpected_structure);
    let summary =
        std::fs::read_to_string(outdir("m3_golden_unexpected_structure").join("summary.txt"))
            .expect("summary");
    assert!(summary.contains("intended_status\tpartial"));
    assert!(summary.contains("structural\t1"));

    let zero_rows = audit(
        true,
        vec![],
        vec![],
        vec![],
        vec![],
        vec![],
        vec![],
        false,
        CandidateSearchStatus::NotPerformed,
    );
    compare_v1_to_m3_golden("zero_row_intended", &zero_rows);
    let zero_summary =
        std::fs::read_to_string(outdir("m3_golden_zero_row_intended").join("summary.txt"))
            .expect("zero-row summary");
    assert!(zero_summary.contains("intended_declared\t0"));
    assert!(zero_summary.contains("intended_missing\tNA"));
    let report =
        std::fs::read_to_string(outdir("m3_golden_unexpected_structure").join("report.md"))
            .expect("legacy report");
    assert!(report.contains("**Intended Edit Outcome: ALERT**"));
    assert!(!report.contains("DV1_"));
}

#[test]
fn aggregate_rendering_is_pure_and_v1_rows_are_v2_prefixes() {
    let refs = refs();
    let expected = event(GdEntry::snp(11, "chr", 10, "T"), &refs);
    let unexpected = event(GdEntry::snp(12, "chr", 20, "G"), &refs);
    let candidate_site = site("SITE_000001", 18, 22);
    let aggregate = audit(
        true,
        vec![assessment(
            "edit",
            &expected.event_id,
            Some(&unexpected.event_id),
        )],
        vec![expected.event_id.clone()],
        vec![mutation(
            &unexpected,
            MutationClass::ScatteredSnv,
            Some("sos_widney2014"),
        )],
        vec![expected.clone(), unexpected.clone()],
        vec![candidate_site.clone()],
        vec![association(&unexpected, &candidate_site, None)],
        true,
        CandidateSearchStatus::PerformedWithCandidates,
    );
    let v1 = outdir("pure_v1");
    let v2a = outdir("pure_v2a");
    let v2b = outdir("pure_v2b");
    write_product_outputs(&aggregate, SchemaVersion::V1, &v1).expect("v1 render");
    write_product_outputs(&aggregate, SchemaVersion::V2, &v2a).expect("v2 render");
    write_product_outputs(&aggregate, SchemaVersion::V2, &v2b).expect("second v2 render");
    for name in [
        "summary.txt",
        "unintended.tsv",
        "edit_outcomes.tsv",
        "post_edit_variants.tsv",
        "offtarget_sites.tsv",
        "mutation_offtarget_links.tsv",
        "report.md",
    ] {
        assert_eq!(
            std::fs::read_to_string(v2a.join(name)).expect("v2a"),
            std::fs::read_to_string(v2b.join(name)).expect("v2b"),
            "{name} must be a pure aggregate rendering"
        );
    }
    for name in [
        "unintended.tsv",
        "edit_outcomes.tsv",
        "post_edit_variants.tsv",
        "offtarget_sites.tsv",
        "mutation_offtarget_links.tsv",
    ] {
        let v1_text = std::fs::read_to_string(v1.join(name)).expect("v1");
        let v1_rows: Vec<Vec<&str>> = v1_text
            .lines()
            .map(|line| line.split('\t').collect())
            .collect();
        let v2_text = std::fs::read_to_string(v2a.join(name)).expect("v2");
        let v2_rows: Vec<Vec<&str>> = v2_text
            .lines()
            .map(|line| line.split('\t').collect())
            .collect();
        assert_eq!(v1_rows.len(), v2_rows.len());
        for (v1_row, v2_row) in v1_rows.iter().zip(v2_rows) {
            assert_eq!(*v1_row, v2_row[..v1_row.len()]);
        }
    }
    assert!(std::fs::read_to_string(v2a.join("report.md"))
        .expect("v2 report")
        .contains(expected.event_id.as_str()));
    assert!(!std::fs::read_to_string(v1.join("report.md"))
        .expect("v1 report")
        .contains(expected.event_id.as_str()));
    let summary = std::fs::read_to_string(v2a.join("summary.txt")).expect("summary");
    let report = std::fs::read_to_string(v2a.join("report.md")).expect("report");
    let links = std::fs::read_to_string(v2a.join("mutation_offtarget_links.tsv")).expect("links");
    assert!(summary.contains("intended_status\tpartial"));
    assert!(summary.contains("scattered_snv\t1"));
    assert!(report.contains("Intended Edit Outcome: ALERT"));
    assert!(report.contains("Total post-edit differential variants:** 1"));
    assert!(report.contains(expected.event_id.as_str()));
    assert!(report.contains(unexpected.event_id.as_str()));
    assert_eq!(aggregate.summary.association_count, 1);
    assert_eq!(links.lines().count(), 2);
}

#[test]
fn summary_preserves_legacy_zero_row_and_unexpected_structure_semantics() {
    let refs = refs();
    let expected = event(GdEntry::snp(21, "chr", 10, "T"), &refs);
    let unexpected = event(GdEntry::jc(22, "chr", 10, "+", "chr", 30, "-", 0), &refs);
    let zero_rows = audit(
        true,
        vec![],
        vec![],
        vec![],
        vec![],
        vec![],
        vec![],
        false,
        CandidateSearchStatus::NotPerformed,
    );
    let unexpected_structure = audit(
        true,
        vec![assessment(
            "edit",
            &expected.event_id,
            Some(&unexpected.event_id),
        )],
        vec![],
        vec![mutation(&unexpected, MutationClass::Structural, None)],
        vec![expected, unexpected],
        vec![],
        vec![],
        false,
        CandidateSearchStatus::NotPerformed,
    );
    let zero_path = outdir("summary_zero");
    let unexpected_path = outdir("summary_unexpected");
    write_product_outputs(&zero_rows, SchemaVersion::V1, &zero_path).expect("zero render");
    write_product_outputs(&unexpected_structure, SchemaVersion::V1, &unexpected_path)
        .expect("unexpected render");
    assert!(std::fs::read_to_string(zero_path.join("summary.txt"))
        .expect("zero summary")
        .contains("intended_missing\tNA"));
    let summary =
        std::fs::read_to_string(unexpected_path.join("summary.txt")).expect("unexpected summary");
    assert!(summary.contains("intended_status\tpartial"));
    assert!(std::fs::read_to_string(unexpected_path.join("report.md"))
        .expect("unexpected report")
        .contains("Intended Edit Outcome: ALERT"));
}

#[test]
fn failed_publish_leaves_no_partial_product_set() {
    let refs = refs();
    let differential_event = event(GdEntry::snp(31, "chr", 10, "T"), &refs);
    let mut aggregate = audit(
        false,
        vec![],
        vec![],
        vec![mutation(
            &differential_event,
            MutationClass::ScatteredSnv,
            None,
        )],
        vec![differential_event],
        vec![],
        vec![],
        false,
        CandidateSearchStatus::NotPerformed,
    );
    aggregate.unintended[0].event_id = None;
    let path = outdir("render_failure");
    assert!(write_product_outputs(&aggregate, SchemaVersion::V2, &path).is_err());
    assert!(!path.join("unintended.tsv").exists());
    assert!(!path.join("summary.txt").exists());
    assert!(!path.join("report.md").exists());
}
