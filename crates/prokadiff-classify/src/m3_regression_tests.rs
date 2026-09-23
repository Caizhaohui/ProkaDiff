use prokadiff_gd::{GdEntry, GenomeDiff};
use prokadiff_offtarget::{BulgeType, OffTargetSite, Strand};

use crate::association::{
    associate_events_to_sites, canonical_sort_sites, linear_associate_events_to_sites,
    primary_association, project_associations_to_links, CandidateSearchStatus,
};
use crate::audit::{build_audit_result, AnalysisProvenance, GuideRelation, SampleMetadata};
use crate::classify::{classify, ClassifyOptions, MutationClass};
use crate::differential::{DifferentialEvent, EventId};
use crate::{EditorKind, RefContig};

fn make_event(entry: GdEntry) -> DifferentialEvent {
    let refs = vec![
        RefContig {
            name: "chr".into(),
            seq: vec![b'A'; 20_000],
        },
        RefContig {
            name: "chr1".into(),
            seq: vec![b'A'; 20_000],
        },
        RefContig {
            name: "chr2".into(),
            seq: vec![b'A'; 20_000],
        },
    ];
    let canonical =
        crate::CanonicalEvent::from_gd_entry(&entry, &refs).expect("valid canonical event");
    let (event_id, _) = canonical.compute_event_id();
    DifferentialEvent {
        event_id,
        canonical,
        representative: entry,
        merged_source_ids: vec![1],
        evidence: crate::EvidenceReferences::default(),
    }
}

fn sample_metadata() -> SampleMetadata {
    SampleMetadata {
        starter_names: vec!["starter.fastq".into()],
        edited_names: vec!["edited.fastq".into()],
        reference_names: vec!["ref.fa".into()],
        editor: "cas9".into(),
        spacer: Some("GACTGACTGACTGACTGACT".into()),
        pam: Some("NGG".into()),
        threads: 1,
    }
}

fn sample_provenance() -> AnalysisProvenance {
    AnalysisProvenance {
        prokadiff_version: "0.1.0".into(),
        git_commit: "test".into(),
        reference_sha256: None,
        bowtie2_version: None,
        offtarget_search_status: "PERFORMED".into(),
        cfd_scoring_status: "DISABLED".into(),
        hsu_scoring_status: "DISABLED".into(),
        bulge_search_status: "DISABLED".into(),
        run_timestamp: "2026-09-23T00:00:00Z".into(),
    }
}

fn mock_site(id: &str, seq_id: &str, start: u64, end: u64, pam: &str, mm: u32) -> OffTargetSite {
    OffTargetSite {
        site_id: id.into(),
        seq_id: seq_id.into(),
        start,
        end,
        strand: Strand::Plus,
        guide: "GACTGACTGACTGACTGACT".into(),
        target_seq: format!("GACTGACTGACTGACTGACT{pam}"),
        pam: pam.into(),
        mismatches: mm,
        bulge_type: BulgeType::None,
        bulge_size: 0,
        search_backend: "rust_exact".into(),
        cfd_score: Some(1.0),
        hsu_score: Some(100.0),
    }
}

// -----------------------------------------------------------------------------
// Gate 1: No SITE_UNKNOWN anywhere in emitted output or AuditResult
// -----------------------------------------------------------------------------
#[test]
fn test_no_site_unknown_anywhere_in_emitted_output() {
    let snp = GdEntry::snp(1, "chr", 500, "C");
    let event = make_event(snp.clone());

    let classified = crate::classify::ClassifiedMutation {
        entry: snp,
        class: MutationClass::ScatteredSnv,
        pam_profile: None,
        offtarget_mismatch: None,
        distance_to_site: None,
        hypothesis: None,
        event_id: Some(event.event_id.clone()),
    };

    let audit = build_audit_result(
        sample_metadata(),
        vec![],
        &[classified],
        &[],
        &[event],
        &[],
        &[], // no associations
        sample_provenance(),
        &[],
    );

    assert_eq!(audit.variants.len(), 1);
    assert_eq!(audit.variants[0].guide_relation, GuideRelation::None);

    // Verify debug/display strings do not contain SITE_UNKNOWN
    let rendered = format!("{audit:?}");
    assert!(
        !rendered.contains("SITE_UNKNOWN"),
        "output must not contain SITE_UNKNOWN"
    );
}

// -----------------------------------------------------------------------------
// Gate 2: Zero-site mutation produces zero association records
// -----------------------------------------------------------------------------
#[test]
fn test_zero_site_mutation_produces_zero_association_records() {
    let snp = GdEntry::snp(2, "chr", 1000, "G");
    let event = make_event(snp);

    // Sites are far away (>50 bp)
    let sites = vec![
        mock_site("SITE_000001", "chr", 100, 123, "TGG", 0),
        mock_site("SITE_000002", "chr", 200, 223, "CGG", 1),
    ];

    let assocs = associate_events_to_sites(&[event], &sites, 50);
    assert!(
        assocs.is_empty(),
        "isolated mutation must produce 0 associations"
    );
}

// -----------------------------------------------------------------------------
// Gate 3: Search status: NotPerformed vs PerformedNoCandidates vs PerformedWithCandidates
// -----------------------------------------------------------------------------
#[test]
fn test_search_status_three_states() {
    let refs = vec![RefContig {
        name: "chr".into(),
        seq: b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
    }];
    let starter = GenomeDiff {
        metadata: vec![],
        entries: vec![],
    };
    let edited = GenomeDiff {
        metadata: vec![],
        entries: vec![GdEntry::snp(1, "chr", 10, "C")],
    };

    // State 1: NotPerformed (DSB or missing spacer)
    let mut opts_dsb = ClassifyOptions::new(EditorKind::Dsb);
    opts_dsb.spacer = None;
    let res_dsb = classify(&edited, &starter, &[], &refs, &opts_dsb).unwrap();
    assert_eq!(
        res_dsb.candidate_search.status,
        CandidateSearchStatus::NotPerformed
    );
    assert!(res_dsb.candidate_search.sites.is_empty());
    assert!(res_dsb.associations.is_empty());

    // State 2: PerformedNoCandidates (Cas9 with spacer, but poly-A ref has no NGG PAM)
    let mut opts_cas9 = ClassifyOptions::new(EditorKind::Cas9);
    opts_cas9.spacer = Some("GACTGACTGACTGACTGACT".into());
    opts_cas9.pam = Some("NGG".into());
    let res_no_cand = classify(&edited, &starter, &[], &refs, &opts_cas9).unwrap();
    assert_eq!(
        res_no_cand.candidate_search.status,
        CandidateSearchStatus::PerformedNoCandidates
    );
    assert!(res_no_cand.candidate_search.sites.is_empty());
    assert!(res_no_cand.associations.is_empty());

    // State 3: PerformedWithCandidates, but zero associations (site is at 200..=223, mutation at 10)
    let mut ref_seq = b"AAAAAAAAAAAAAAAAAAAA".to_vec();
    ref_seq.resize(200, b'A');
    ref_seq.extend_from_slice(b"GACTGACTGACTGACTGACTTGG");
    ref_seq.extend(std::iter::repeat_n(b'A', 100));
    let refs_with_site = vec![RefContig {
        name: "chr".into(),
        seq: ref_seq,
    }];
    let res_with_cand = classify(&edited, &starter, &[], &refs_with_site, &opts_cas9).unwrap();
    assert_eq!(
        res_with_cand.candidate_search.status,
        CandidateSearchStatus::PerformedWithCandidates
    );
    assert_eq!(res_with_cand.candidate_search.sites.len(), 1);
    assert!(
        res_with_cand.associations.is_empty(),
        "mutation at 10 is >50 bp from site at 200"
    );
}

// -----------------------------------------------------------------------------
// Gate 4: Multiple qualifying sites all survive
// -----------------------------------------------------------------------------
#[test]
fn test_multiple_qualifying_sites_all_survive() {
    let snp = GdEntry::snp(4, "chr", 150, "C");
    let event = make_event(snp);

    // 3 sites within window 50: distances are 10 bp, 20 bp, 30 bp
    let sites = vec![
        mock_site("SITE_000001", "chr", 100, 120, "TGG", 1), // dist = 150 - 120 = 30
        mock_site("SITE_000002", "chr", 110, 130, "CGG", 2), // dist = 150 - 130 = 20
        mock_site("SITE_000003", "chr", 120, 140, "AGG", 0), // dist = 150 - 140 = 10
    ];

    let assocs = associate_events_to_sites(&[event], &sites, 50);
    assert_eq!(
        assocs.len(),
        3,
        "all 3 qualifying candidate associations must be preserved"
    );

    let primary = primary_association(&assocs).expect("primary exists");
    assert_eq!(primary.site_id, "SITE_000003");
    assert_eq!(primary.distance_to_site, 10);
    assert_eq!(primary.pam, "AGG");
}

// -----------------------------------------------------------------------------
// Gate 5: JC breakpoint with multiple qualifying sites produces >2 rows
// -----------------------------------------------------------------------------
#[test]
fn test_jc_breakpoint_multiple_qualifying_sites_produces_more_than_two_rows() {
    // Side1 at chr1:100, Side2 at chr2:500
    let jc = GdEntry::jc(5, "chr1", 100, "+", "chr2", 500, "-", 0);
    let event = make_event(jc);

    let sites = vec![
        // 2 sites near Side 1 (chr1:100)
        mock_site("SITE_000001", "chr1", 70, 90, "TGG", 0), // dist = 10
        mock_site("SITE_000002", "chr1", 110, 130, "CGG", 1), // dist = 10
        // 2 sites near Side 2 (chr2:500)
        mock_site("SITE_000003", "chr2", 470, 490, "TGG", 2), // dist = 10
        mock_site("SITE_000004", "chr2", 510, 530, "GGG", 0), // dist = 10
    ];

    let assocs = associate_events_to_sites(std::slice::from_ref(&event), &sites, 50);
    assert_eq!(
        assocs.len(),
        4,
        "JC with 2 qualifying sites per breakpoint must yield 4 rows"
    );

    let links = project_associations_to_links(&assocs, &[event]);
    assert_eq!(
        links.len(),
        4,
        "projected links must also contain all 4 rows"
    );
    for link in &links {
        assert_eq!(link.mutation_id, "mut_5");
    }
}

// -----------------------------------------------------------------------------
// Gate 6: MC-only evidence produces zero MutationSiteAssociation
// -----------------------------------------------------------------------------
#[test]
fn test_mc_only_evidence_produces_no_mutation_site_association() {
    let mc_entry = GdEntry::mc(6, "chr1", 100, 200, 0, 0);
    let refs = [RefContig {
        name: "chr1".into(),
        seq: vec![b'A'; 1000],
    }];
    let dummy_contig = crate::CanonicalContig::from_ref_contig(&refs[0]).unwrap();
    let mc_event = DifferentialEvent {
        event_id: EventId::parse("DV1_00000000000000000000000000000006").unwrap(),
        canonical: crate::CanonicalEvent::Snp {
            contig: dummy_contig,
            position: 100,
            new_base: b'T',
        },
        representative: mc_entry,
        merged_source_ids: vec![6],
        evidence: Default::default(),
    };

    let sites = vec![mock_site("SITE_000001", "chr1", 110, 130, "TGG", 0)];

    let assocs = associate_events_to_sites(std::slice::from_ref(&mc_event), &sites, 50);
    assert!(
        assocs.is_empty(),
        "MC-only record must produce zero associations"
    );

    let linear_assocs = linear_associate_events_to_sites(&[mc_event], &sites, 50);
    assert!(
        linear_assocs.is_empty(),
        "MC-only record must produce zero linear associations"
    );
}

// -----------------------------------------------------------------------------
// Gate 7: Reported site_id, PAM, mismatch, and distance all come from the same OffTargetSite
// -----------------------------------------------------------------------------
#[test]
fn test_reported_site_id_pam_mismatch_distance_from_same_site() {
    let snp = GdEntry::snp(7, "chr1", 100, "C");
    let event = make_event(snp);

    // Site 1: dist = 15, mm = 3, pam = "TGG"
    // Site 2: dist = 5, mm = 1, pam = "CGG" (Primary because closer)
    let sites = vec![
        mock_site("SITE_000001", "chr1", 80, 85, "TGG", 3),
        mock_site("SITE_000002", "chr1", 90, 95, "CGG", 1),
    ];

    let assocs = associate_events_to_sites(&[event], &sites, 50);
    let primary = primary_association(&assocs).expect("has primary");

    assert_eq!(primary.site_id, "SITE_000002");
    assert_eq!(primary.pam, "CGG");
    assert_eq!(primary.mismatches, 1);
    assert_eq!(primary.distance_to_site, 5);
    assert_eq!(primary.site_start, 90);
    assert_eq!(primary.site_end, 95);
}

// -----------------------------------------------------------------------------
// Gate 8: Query PAM NGG with observed PAM CGG outputs CGG
// -----------------------------------------------------------------------------
#[test]
fn test_query_pam_ngg_with_observed_pam_cgg_outputs_cgg() {
    // Spacer: 20 bp, PAM in genome is CGG
    let spacer = "GACTGACTGACTGACTGACT";
    let mut ref_seq = b"AAAAAAAAAA".to_vec();
    ref_seq.extend_from_slice(spacer.as_bytes());
    ref_seq.extend_from_slice(b"CGG"); // genomic PAM is CGG, query is NGG
    ref_seq.extend(std::iter::repeat_n(b'A', 100));

    let refs = vec![RefContig {
        name: "chr1".into(),
        seq: ref_seq,
    }];
    let starter = GenomeDiff {
        metadata: vec![],
        entries: vec![],
    };
    // Mutation at position 40 (7 bp past site 11..=33)
    let edited = GenomeDiff {
        metadata: vec![],
        entries: vec![GdEntry::snp(1, "chr1", 40, "C")],
    };

    let mut opts = ClassifyOptions::new(EditorKind::Cas9);
    opts.spacer = Some(spacer.into());
    opts.pam = Some("NGG".into()); // Query PAM pattern is NGG

    let res = classify(&edited, &starter, &[], &refs, &opts).unwrap();
    assert_eq!(res.unintended.len(), 1);
    assert_eq!(res.unintended[0].class, MutationClass::NearHomolog);
    // Observed PAM must be CGG, NOT NGG!
    assert_eq!(res.unintended[0].pam_profile.as_deref(), Some("CGG"));
    assert_eq!(res.associations.len(), 1);
    assert_eq!(res.associations[0].pam, "CGG");
}

// -----------------------------------------------------------------------------
// Gate 9: Equal-distance/equal-mismatch tie remains deterministic without dropping alternatives
// -----------------------------------------------------------------------------
#[test]
fn test_equal_distance_equal_mismatch_tie_remains_deterministic_without_dropping() {
    let snp = GdEntry::snp(9, "chr1", 100, "C");
    let event = make_event(snp);

    // Two sites equidistant from 100:
    // Left site: 85..=90 -> distance = 10
    // Right site: 110..=115 -> distance = 10
    let site_left = mock_site("SITE_000002", "chr1", 85, 90, "TGG", 1);
    let site_right = mock_site("SITE_000001", "chr1", 110, 115, "CGG", 1);
    let mut sites = vec![site_left, site_right];
    canonical_sort_sites(&mut sites);

    let assocs = associate_events_to_sites(&[event], &sites, 50);
    assert_eq!(assocs.len(), 2, "both alternatives must be preserved");

    let primary = primary_association(&assocs).expect("primary exists");
    // Canonical tie-breaker picks smaller site_id
    assert_eq!(primary.site_id, "SITE_000001");
}

// -----------------------------------------------------------------------------
// Gate 10: Exact vs bulge evidence metadata remains distinguishable internally
// -----------------------------------------------------------------------------
#[test]
fn test_exact_vs_bulge_evidence_metadata_distinguishable_internally() {
    let snp = GdEntry::snp(10, "chr1", 100, "C");
    let event = make_event(snp);

    let mut site_exact = mock_site("SITE_000001", "chr1", 90, 95, "TGG", 0);
    site_exact.search_backend = "rust_exact".into();
    site_exact.bulge_type = BulgeType::None;
    site_exact.bulge_size = 0;

    let mut site_bulge = mock_site("SITE_000002", "chr1", 105, 110, "TGG", 0);
    site_bulge.search_backend = "rust_bulge".into();
    site_bulge.bulge_type = BulgeType::Dna;
    site_bulge.bulge_size = 1;

    let assocs = associate_events_to_sites(&[event], &[site_exact, site_bulge], 50);
    assert_eq!(assocs.len(), 2);

    assert_eq!(assocs[0].search_backend, "rust_exact");
    assert_eq!(assocs[0].bulge_type, BulgeType::None);
    assert_eq!(assocs[0].bulge_size, 0);

    assert_eq!(assocs[1].search_backend, "rust_bulge");
    assert_eq!(assocs[1].bulge_type, BulgeType::Dna);
    assert_eq!(assocs[1].bulge_size, 1);
}

// -----------------------------------------------------------------------------
// Gate 11: Shuffled GD input gives identical association output
// -----------------------------------------------------------------------------
#[test]
fn test_shuffled_gd_input_gives_identical_association_output() {
    let e1 = make_event(GdEntry::snp(1, "chr1", 100, "C"));
    let e2 = make_event(GdEntry::del(2, "chr1", 200, 10));
    let e3 = make_event(GdEntry::jc(3, "chr1", 300, "+", "chr1", 400, "-", 0));

    let sites = vec![
        mock_site("SITE_000001", "chr1", 95, 115, "TGG", 0),
        mock_site("SITE_000002", "chr1", 205, 225, "CGG", 1),
        mock_site("SITE_000003", "chr1", 305, 325, "AGG", 0),
    ];

    let order1 = vec![e1.clone(), e2.clone(), e3.clone()];
    let order2 = vec![e3.clone(), e1.clone(), e2.clone()];
    let order3 = vec![e2.clone(), e3.clone(), e1.clone()];

    let assocs1 = associate_events_to_sites(&order1, &sites, 50);
    let assocs2 = associate_events_to_sites(&order2, &sites, 50);
    let assocs3 = associate_events_to_sites(&order3, &sites, 50);

    assert_eq!(
        assocs1, assocs2,
        "permuting input events must yield identical associations"
    );
    assert_eq!(
        assocs1, assocs3,
        "permuting input events must yield identical associations"
    );
}

// -----------------------------------------------------------------------------
// Gate 12: Indexed search equals a linear reference implementation on fixtures
// -----------------------------------------------------------------------------
#[test]
fn test_indexed_search_equals_linear_reference_implementation() {
    // Generate 50 sites across 2 contigs
    let mut sites = Vec::new();
    for i in 1..=25 {
        sites.push(mock_site(
            &format!("SITE_{i:06}"),
            "chr1",
            i * 100,
            i * 100 + 22,
            "TGG",
            (i % 4) as u32,
        ));
    }
    for i in 26..=50 {
        sites.push(mock_site(
            &format!("SITE_{i:06}"),
            "chr2",
            (i - 25) * 100,
            (i - 25) * 100 + 22,
            "CGG",
            (i % 3) as u32,
        ));
    }
    canonical_sort_sites(&mut sites);

    // Generate diverse events (Points, Spans, Junctions)
    let mut events = Vec::new();
    // Points
    events.push(make_event(GdEntry::snp(1, "chr1", 105, "C")));
    events.push(make_event(GdEntry::snp(2, "chr1", 130, "G")));
    events.push(make_event(GdEntry::snp(3, "chr1", 5000, "T"))); // isolated
                                                                 // Spans
    events.push(make_event(GdEntry::del(4, "chr2", 210, 15)));
    events.push(make_event(GdEntry::sub(5, "chr2", 400, 20, "ACTG")));
    // Junctions
    events.push(make_event(GdEntry::jc(
        6, "chr1", 305, "+", "chr2", 305, "-", 0,
    )));
    events.push(make_event(GdEntry::jc(
        7, "chr1", 900, "+", "chr1", 950, "-", 0,
    )));

    for window in [10, 25, 50, 100] {
        let linear = linear_associate_events_to_sites(&events, &sites, window);
        let indexed = associate_events_to_sites(&events, &sites, window);
        assert_eq!(
            indexed, linear,
            "indexed search must be bitwise identical to linear reference for window {window}"
        );
    }
}

// -----------------------------------------------------------------------------
// Gate 13: Regression test: No positional association fallback when event_id is None
// Two unrelated events sharing the same coordinate cannot inherit each other's association
// -----------------------------------------------------------------------------
#[test]
fn test_positional_association_fallback_removed_two_unrelated_events_same_coordinate() {
    let snp1 = GdEntry::snp(1, "chr1", 100, "C");
    let ev1 = make_event(snp1.clone());

    let snp2 = GdEntry::snp(2, "chr1", 100, "G"); // Same coordinate 100, different mutation

    let site = mock_site("SITE_000001", "chr1", 90, 95, "CGG", 0);
    let assocs =
        associate_events_to_sites(std::slice::from_ref(&ev1), std::slice::from_ref(&site), 50);
    assert_eq!(assocs.len(), 1);
    assert_eq!(assocs[0].event_id, ev1.event_id);
    assert_eq!(assocs[0].mutation_position, 100);

    // M1 has real event_id, M2 has event_id: None (unbacked)
    let cm1 = crate::classify::ClassifiedMutation {
        entry: snp1,
        class: MutationClass::NearHomolog,
        pam_profile: Some("CGG".into()),
        offtarget_mismatch: Some(0),
        distance_to_site: Some(5),
        hypothesis: None,
        event_id: Some(ev1.event_id.clone()),
    };
    let cm2 = crate::classify::ClassifiedMutation {
        entry: snp2,
        class: MutationClass::ScatteredSnv,
        pam_profile: None,
        offtarget_mismatch: None,
        distance_to_site: None,
        hypothesis: None,
        event_id: None,
    };

    let audit = build_audit_result(
        sample_metadata(),
        vec![],
        &[cm1, cm2],
        &[],
        &[ev1],
        &[site],
        &assocs,
        sample_provenance(),
        &[],
    );

    assert_eq!(audit.variants.len(), 2);
    // Variant 1 with real EventId gets the association
    assert!(matches!(
        audit.variants[0].guide_relation,
        GuideRelation::CandidateOffTarget { ref site_id, .. } if site_id == "SITE_000001"
    ));
    // Variant 2 sharing coordinate 100 but having event_id: None MUST NOT inherit association
    assert_eq!(
        audit.variants[1].guide_relation,
        GuideRelation::None,
        "mutation with event_id=None must have GuideRelation::None and cannot inherit association by position"
    );
}

// -----------------------------------------------------------------------------
// Gate 14: Regression test: None / Some(0.0) / Some(0.5) do not affect M3 association ranking
// -----------------------------------------------------------------------------
#[test]
fn test_unvalidated_scores_do_not_affect_association_ranking() {
    let dummy_id = EventId::parse("DV1_00000000000000000000000000000001").unwrap();

    let make_assoc = |site_id: &str, cfd: Option<f64>, hsu: Option<f64>| {
        crate::association::MutationSiteAssociation {
            event_id: dummy_id.clone(),
            mutation_type: "SNP".into(),
            mutation_position: 100,
            junction_side: None,
            site_id: site_id.into(),
            site_seq_id: "chr1".into(),
            site_start: 90,
            site_end: 95,
            site_strand: Strand::Plus,
            distance_to_site: 5,
            mismatches: 1,
            pam: "CGG".into(),
            cfd_score: cfd,
            hsu_score: hsu,
            search_backend: "rust_exact".into(),
            bulge_type: BulgeType::None,
            bulge_size: 0,
            target_seq: "GACTGACTGACTGACTGACTCGG".into(),
            association_window: 50,
        }
    };

    let a1 = make_assoc("SITE_000001", None, None);
    let a2 = make_assoc("SITE_000002", Some(0.5), Some(50.0));
    let a3 = make_assoc("SITE_000003", Some(0.0), Some(0.0));

    for perm in [
        vec![a1.clone(), a2.clone(), a3.clone()],
        vec![a2.clone(), a3.clone(), a1.clone()],
        vec![a3.clone(), a1.clone(), a2.clone()],
        vec![a2.clone(), a1.clone(), a3.clone()],
    ] {
        let primary = primary_association(&perm).expect("primary exists");
        assert_eq!(
            primary.site_id, "SITE_000001",
            "scores None/Some(0.0)/Some(0.5) must not affect primary ranking; SITE_000001 must win"
        );

        let mut sorted = perm;
        crate::association::canonical_sort_associations(&mut sorted);
        assert_eq!(sorted[0].site_id, "SITE_000001");
        assert_eq!(sorted[1].site_id, "SITE_000002");
        assert_eq!(sorted[2].site_id, "SITE_000003");
    }
}

// -----------------------------------------------------------------------------
// Gate 15: Real-path MC exclusion test through build_differential_events()
// -----------------------------------------------------------------------------
#[test]
fn test_real_path_mc_exclusion_through_build_differential_events() {
    use crate::differential::{build_differential_events, EvidenceKind};

    let refs = vec![RefContig {
        name: "chr1".into(),
        seq: vec![b'A'; 5000],
    }];

    let starter = GenomeDiff {
        metadata: vec![],
        entries: vec![],
    };

    // Edited GD contains:
    // - MC evidence entry (id: 1)
    // - DEL product mutation (id: 2) with parent_ids: [1]
    let mc_entry = GdEntry::mc(1, "chr1", 500, 600, 0, 0);
    let mut del_entry = GdEntry::del(2, "chr1", 520, 20);
    del_entry.parent_ids = vec![1];

    let edited = GenomeDiff {
        metadata: vec![],
        entries: vec![mc_entry, del_entry],
    };

    let diff_res = build_differential_events(&starter, &edited, &refs)
        .expect("build_differential_events succeeds");

    // 1. Verify MC does NOT become DifferentialEvent
    assert_eq!(diff_res.events.len(), 1);
    let event = &diff_res.events[0];
    assert_eq!(event.representative.kind, prokadiff_gd::GdKind::Del);
    assert_ne!(event.representative.kind, prokadiff_gd::GdKind::Mc);

    // 2. Verify MC receives no EventId
    for e in &diff_res.events {
        assert_ne!(
            e.representative.id, 1,
            "MC entry id 1 must not be an event representative"
        );
    }

    // 3. Verify MC produces no MutationSiteAssociation
    let site = mock_site("SITE_000001", "chr1", 510, 530, "CGG", 0);
    let assocs = associate_events_to_sites(&diff_res.events, &[site], 50);
    assert_eq!(assocs.len(), 1);
    assert_eq!(assocs[0].mutation_type, "DEL");
    assert_eq!(assocs[0].event_id, event.event_id);
    for assoc in &assocs {
        assert_ne!(assoc.mutation_type, "MC");
    }

    // 4. Verify evidence-only behavior remains intact:
    // DEL's evidence references contain MC reference with id 1
    assert!(event.evidence.mc.is_some());
    let mc_refs = event.evidence.mc.as_ref().unwrap();
    assert_eq!(mc_refs.len(), 1);
    assert_eq!(mc_refs[0].gd_id, 1);
    assert_eq!(mc_refs[0].kind, EvidenceKind::Mc);

    // Additional check: An edited GD with ONLY MC produces exactly 0 DifferentialEvents
    let edited_only_mc = GenomeDiff {
        metadata: vec![],
        entries: vec![GdEntry::mc(10, "chr1", 100, 200, 0, 0)],
    };
    let diff_only_mc = build_differential_events(&starter, &edited_only_mc, &refs).unwrap();
    assert!(
        diff_only_mc.events.is_empty(),
        "MC-only GD must produce zero differential events"
    );
}

// -----------------------------------------------------------------------------
// Gate 16: Index boundary oracle test (indexed vs linear byte-identical)
// -----------------------------------------------------------------------------
#[test]
fn test_index_boundary_oracle_coverage() {
    let contig_len = 10_000u64;

    // Events:
    // 1. Query SNP at position 100
    let e_100 = make_event(GdEntry::snp(1, "chr1", 100, "C"));
    // 2. Query at contig left-edge: position 1
    let e_left = make_event(GdEntry::snp(2, "chr1", 1, "G"));
    // 3. Query at contig right-edge: position contig_len
    let e_right = make_event(GdEntry::snp(3, "chr1", contig_len, "T"));

    let events = vec![e_100, e_left, e_right];

    let mut sites = vec![
        // 1. Site exactly at +window boundary for query at 100:
        // dist = s_start - q_end = 150 - 100 = 50 == window (when window=50)
        mock_site("SITE_000001", "chr1", 150, 172, "CGG", 0),
        // Site just outside +window boundary for query at 100: dist = 151 - 100 = 51 > 50
        mock_site("SITE_000002", "chr1", 151, 173, "CGG", 0),
        // 2. Site ending exactly at query_start - window for query at 100:
        // end = 100 - 50 = 50. Site start 28, end 50 (len 23). dist = 100 - 50 = 50 == window
        mock_site("SITE_000003", "chr1", 28, 50, "CGG", 1),
        // Site ending at query_start - window - 1: end = 49. dist = 100 - 49 = 51 > 50
        mock_site("SITE_000004", "chr1", 27, 49, "CGG", 1),
        // 3. Long site whose start is outside naive range (start < q_start - window = 50),
        // but whose interval still overlaps the allowed window:
        // start = 10, end = 60 (length 51). dist = 100 - 60 = 40 <= 50.
        mock_site("SITE_000005", "chr1", 10, 60, "TGG", 2),
        // Another long site: start = 5, end = 49. dist = 100 - 49 = 51 > 50 (outside)
        mock_site("SITE_000006", "chr1", 5, 49, "TGG", 2),
        // 4. Contig-edge coordinates:
        // Left contig edge: site at 1..23 (overlaps query at 1, dist = 0)
        mock_site("SITE_000007", "chr1", 1, 23, "CGG", 0),
        // Left contig edge: site at 51..73 (dist to 1 is 51 - 1 = 50 == window)
        mock_site("SITE_000008", "chr1", 51, 73, "CGG", 0),
        // Right contig edge: site at contig_len - 22 ..= contig_len (overlaps query at contig_len, dist = 0)
        mock_site("SITE_000009", "chr1", contig_len - 22, contig_len, "CGG", 0),
        // Right contig edge: site ending at contig_len - 50 (dist to contig_len is 50 == window)
        mock_site(
            "SITE_000010",
            "chr1",
            contig_len - 72,
            contig_len - 50,
            "CGG",
            0,
        ),
    ];
    canonical_sort_sites(&mut sites);

    for w in [10, 25, 50, 100] {
        let linear = linear_associate_events_to_sites(&events, &sites, w);
        let indexed = associate_events_to_sites(&events, &sites, w);
        assert_eq!(
            indexed, linear,
            "indexed search must be bitwise identical to linear reference for window {w}"
        );
    }
}

// -----------------------------------------------------------------------------
// Gate 17: Primary association final tie test with JC event
// -----------------------------------------------------------------------------
#[test]
fn test_primary_association_final_tie_deterministic_jc() {
    let dummy_id = EventId::parse("DV1_00000000000000000000000000000009").unwrap();

    let make_jc_assoc =
        |site_id: &str, side: crate::geometry::JunctionSideTag, cfd: Option<f64>| {
            crate::association::MutationSiteAssociation {
                event_id: dummy_id.clone(),
                mutation_type: "JC".into(),
                mutation_position: 100,
                junction_side: Some(side),
                site_id: site_id.into(),
                site_seq_id: "chr1".into(),
                site_start: 90,
                site_end: 95,
                site_strand: Strand::Plus,
                distance_to_site: 5,
                mismatches: 1,
                pam: "CGG".into(),
                cfd_score: cfd,
                hsu_score: None,
                search_backend: "rust_exact".into(),
                bulge_type: BulgeType::None,
                bulge_size: 0,
                target_seq: "GACTGACTGACTGACTGACTCGG".into(),
                association_window: 50,
            }
        };

    // Sub-case A: Side tie-break (Side1 vs Side2)
    // Candidate A: Side1, site_id "SITE_000002"
    // Candidate B: Side2, site_id "SITE_000001"
    // They tie on distance (5), mismatch (1), and score (Some(0.5)).
    // Because junction_side is checked before site_id, Side1 (Candidate A) MUST win,
    // even though Candidate B has lexicographically smaller site_id!
    let assoc_side1 = make_jc_assoc(
        "SITE_000002",
        crate::geometry::JunctionSideTag::Side1,
        Some(0.5),
    );
    let assoc_side2 = make_jc_assoc(
        "SITE_000001",
        crate::geometry::JunctionSideTag::Side2,
        Some(0.5),
    );

    let pair_a_b = [assoc_side1.clone(), assoc_side2.clone()];
    let pair_b_a = [assoc_side2, assoc_side1];
    let p_a_b = primary_association(&pair_a_b).unwrap();
    let p_b_a = primary_association(&pair_b_a).unwrap();
    assert_eq!(
        p_a_b.junction_side,
        Some(crate::geometry::JunctionSideTag::Side1)
    );
    assert_eq!(p_a_b.site_id, "SITE_000002");
    assert_eq!(
        p_b_a.junction_side,
        Some(crate::geometry::JunctionSideTag::Side1)
    );
    assert_eq!(p_b_a.site_id, "SITE_000002");

    // Sub-case B: site_id tie-break when junction_side also ties (both Side1)
    // Candidate C: Side1, site_id "SITE_000001"
    // Candidate D: Side1, site_id "SITE_000002"
    let assoc_c = make_jc_assoc("SITE_000001", crate::geometry::JunctionSideTag::Side1, None);
    let assoc_d = make_jc_assoc("SITE_000002", crate::geometry::JunctionSideTag::Side1, None);

    let pair_c_d = [assoc_c.clone(), assoc_d.clone()];
    let pair_d_c = [assoc_d, assoc_c];
    let p_c_d = primary_association(&pair_c_d).unwrap();
    let p_d_c = primary_association(&pair_d_c).unwrap();
    assert_eq!(p_c_d.site_id, "SITE_000001");
    assert_eq!(p_d_c.site_id, "SITE_000001");
}
