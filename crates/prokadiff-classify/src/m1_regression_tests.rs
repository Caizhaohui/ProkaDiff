use prokadiff_gd::{GdEntry, GenomeDiff};

use crate::{
    build_audit_result, build_differential_events, classify, AnalysisProvenance, CanonicalEvent,
    ClassifyOptions, DifferentialError, EditorKind, RefContig, SampleMetadata,
};

fn gd(entries: Vec<GdEntry>) -> GenomeDiff {
    GenomeDiff {
        metadata: vec![("GENOME_DIFF".into(), "1.0".into())],
        entries,
    }
}

fn refs() -> Vec<RefContig> {
    let sequence: Vec<u8> = b"ACGT".iter().copied().cycle().take(1_000).collect();
    vec![
        RefContig {
            name: "chr1".into(),
            seq: sequence.clone(),
        },
        RefContig {
            name: "chr2".into(),
            seq: sequence,
        },
    ]
}

fn event_ids(starter: Vec<GdEntry>, edited: Vec<GdEntry>) -> Vec<String> {
    build_differential_events(&gd(starter), &gd(edited), &refs())
        .expect("valid differential inputs")
        .events
        .into_iter()
        .map(|event| event.event_id.to_string())
        .collect()
}

fn dsb_options() -> ClassifyOptions {
    ClassifyOptions {
        editor: EditorKind::Dsb,
        spacer: None,
        pam: None,
        near_distance: 50,
        max_mismatches: 3,
        hypothesis: false,
    }
}

#[test]
fn classify_requires_real_reference_contigs() {
    // Given a product mutation whose contig is absent from the reference skeleton.
    let edited = gd(vec![GdEntry::snp(1, "chr1", 100, "T")]);

    // When classification attempts to build a differential event.
    let result = classify(&edited, &gd(vec![]), &[], &[], &dsb_options());

    // Then the typed differential error reaches the caller.
    assert!(matches!(result, Err(DifferentialError::ContigNotFound(name)) if name == "chr1"));
}

#[test]
fn unrelated_event_insertion_and_removal_does_not_change_event_id() {
    // Given the same target event with and without an unrelated mutation.
    let target = GdEntry::snp(1, "chr1", 100, "T");
    let without_unrelated = classify(
        &gd(vec![target.clone()]),
        &gd(vec![]),
        &[],
        &refs(),
        &dsb_options(),
    )
    .expect("real reference is present");
    let with_unrelated = classify(
        &gd(vec![GdEntry::del(2, "chr1", 300, 7), target]),
        &gd(vec![]),
        &[],
        &refs(),
        &dsb_options(),
    )
    .expect("real reference is present");

    // When the target event is located in each result.
    let target_id = |result: &crate::ClassifyResult| {
        result
            .differential_events
            .iter()
            .find(|event| event.representative.id == 1)
            .expect("target event remains present")
            .event_id
            .clone()
    };

    // Then unrelated coordinates do not alter its content-derived identity.
    assert_eq!(target_id(&without_unrelated), target_id(&with_unrelated));
}

#[test]
fn ra_mc_and_jc_evidence_reaches_audit_variants() {
    // Given three product mutations backed by actual RA, MC, and JC records.
    let mut ra = GdEntry::ra(10, "chr1", 100, 0, "A", "T");
    ra.attrs.insert("pd_ra_support_reads".into(), "21".into());
    let mut snp = GdEntry::snp(1, "chr1", 100, "T");
    snp.parent_ids = vec![10];

    let mc = GdEntry::mc(20, "chr1", 200, 206, 0, 0);
    let mut deletion = GdEntry::del(2, "chr1", 200, 7);
    deletion.parent_ids = vec![20];

    let mut junction = GdEntry::jc(3, "chr1", 400, "+", "chr2", 500, "-", 0);
    junction
        .attrs
        .insert("pd_support_reads".into(), "17".into());
    let edited = gd(vec![ra, snp, mc, deletion, junction]);

    // When classification is converted into the unified audit model.
    let classified = classify(&edited, &gd(vec![]), &[], &refs(), &dsb_options())
        .expect("evidence-backed mutations are valid");
    let audit = build_audit_result(
        SampleMetadata {
            starter_names: vec!["starter.fastq".into()],
            edited_names: vec!["edited.fastq".into()],
            reference_names: vec!["ref.fa".into()],
            editor: "dsb".into(),
            spacer: None,
            pam: None,
            threads: 1,
        },
        vec![],
        &classified.unintended,
        &classified.intended_event_ids,
        &classified.differential_events,
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
            run_timestamp: "2026-09-22T00:00:00Z".into(),
        },
        &[],
    );

    // Then audit/report evidence summaries come from DifferentialEvent evidence.
    let snp_audit = audit
        .variants
        .iter()
        .find(|variant| variant.entry.id == 1)
        .expect("SNP audit variant");
    assert_eq!(snp_audit.evidence.ra, Some(true));
    assert_eq!(snp_audit.evidence.supporting_reads, Some(21));

    let del_audit = audit
        .variants
        .iter()
        .find(|variant| variant.entry.id == 2)
        .expect("DEL audit variant");
    assert_eq!(del_audit.evidence.mc, Some(true));

    let jc_audit = audit
        .variants
        .iter()
        .find(|variant| variant.entry.id == 3)
        .expect("JC audit variant");
    assert_eq!(jc_audit.evidence.jc, Some(true));
    assert_eq!(jc_audit.evidence.supporting_reads, Some(17));
}

#[test]
fn ins_position_zero_is_rejected() {
    // Given an insertion whose 1-based anchor is zero.
    let entry = GdEntry::ins(1, "chr1", 0, "A");

    // When it is canonicalized.
    let result = CanonicalEvent::from_gd_entry(&entry, &refs());

    // Then the coordinate is rejected at the boundary.
    assert!(matches!(
        result,
        Err(DifferentialError::InvalidPosition { position: 0, .. })
    ));
}

#[test]
fn snp_sub_and_ins_reject_non_iupac_alleles() {
    // Given mutation alleles outside the documented DNA/IUPAC alphabet.
    let invalid = [
        GdEntry::snp(1, "chr1", 10, "Z"),
        GdEntry::sub(2, "chr1", 10, 2, "AZ"),
        GdEntry::ins(3, "chr1", 10, "A?"),
    ];

    // When each mutation is canonicalized, then each invalid allele is rejected.
    for entry in invalid {
        assert!(matches!(
            CanonicalEvent::from_gd_entry(&entry, &refs()),
            Err(DifferentialError::InvalidAllele { .. })
        ));
    }
}

#[test]
fn overflowing_coordinate_span_is_rejected_without_wrapping() {
    let entry = GdEntry::del(1, "chr1", 2, u64::MAX);
    let result = CanonicalEvent::from_gd_entry(&entry, &refs());
    assert!(matches!(result, Err(DifferentialError::InvalidSize { .. })));
}

#[test]
fn structural_del_tolerance_accepts_four_and_five_but_rejects_six_bp() {
    // Given one inherited seven-base deletion.
    let starter = vec![GdEntry::del(1, "chr1", 100, 7)];

    // When edited coordinates differ by four, five, and six bases.
    let within_four = event_ids(starter.clone(), vec![GdEntry::del(2, "chr1", 104, 7)]);
    let within_five = event_ids(starter.clone(), vec![GdEntry::del(3, "chr1", 105, 7)]);
    let beyond_six = event_ids(starter.clone(), vec![GdEntry::del(4, "chr1", 106, 7)]);
    let within_negative_four = event_ids(starter.clone(), vec![GdEntry::del(5, "chr1", 96, 7)]);
    let within_negative_five = event_ids(starter.clone(), vec![GdEntry::del(6, "chr1", 95, 7)]);
    let beyond_negative_six = event_ids(starter, vec![GdEntry::del(7, "chr1", 94, 7)]);

    // Then only the event beyond the fixed five-base boundary survives subtraction.
    assert!(within_four.is_empty());
    assert!(within_five.is_empty());
    assert_eq!(beyond_six.len(), 1);
    assert!(within_negative_four.is_empty());
    assert!(within_negative_five.is_empty());
    assert_eq!(beyond_negative_six.len(), 1);
}

#[test]
fn mixed_short_and_structural_deletions_do_not_tolerant_match() {
    // Given a short starter deletion and a structural edited deletion with nearby boundaries.
    let starter = vec![GdEntry::del(1, "chr1", 100, 2)];
    let edited = vec![GdEntry::del(2, "chr1", 99, 3)];

    // When differential subtraction is performed.
    let surviving = event_ids(starter, edited);

    // Then tolerance is not applied unless both deletions are structural.
    assert_eq!(surviving.len(), 1);
}

#[test]
fn ambiguous_two_by_two_matching_maximizes_valid_matches() {
    // Given a graph where the first edited JC can match both starter JCs, while the
    // second edited JC can match only the first starter JC.
    let starter = vec![
        GdEntry::jc(1, "chr1", 10, "+", "chr2", 10, "-", 0),
        GdEntry::jc(2, "chr1", 15, "+", "chr2", 15, "-", 0),
    ];
    let edited = vec![
        GdEntry::jc(3, "chr1", 11, "+", "chr2", 11, "-", 0),
        GdEntry::jc(4, "chr1", 12, "+", "chr2", 5, "-", 0),
    ];

    // When tolerant subtraction chooses a one-to-one matching.
    let surviving = event_ids(starter, edited);

    // Then it finds the cardinality-two matching instead of the greedy cardinality-one match.
    assert!(surviving.is_empty());
}

#[test]
fn jc_signed_overlap_must_be_compatible_for_tolerant_matching() {
    // Given coordinate-compatible junctions with different signed overlap semantics.
    let starter = vec![GdEntry::jc(1, "chr1", 100, "+", "chr2", 200, "-", -2)];
    let edited = vec![GdEntry::jc(2, "chr1", 103, "+", "chr2", 204, "-", 2)];

    // When tolerant subtraction is performed, then the edited JC survives.
    assert_eq!(event_ids(starter, edited).len(), 1);
}

#[test]
fn mob_duplication_size_must_be_compatible_for_tolerant_matching() {
    // Given coordinate-compatible MOB calls with different target-site duplications.
    let starter = vec![GdEntry::mob(1, "chr1", 100, "IS1", "+", 3)];
    let edited = vec![GdEntry::mob(2, "chr1", 104, "IS1", "+", 4)];

    // When tolerant subtraction is performed, then the edited MOB survives.
    assert_eq!(event_ids(starter, edited).len(), 1);
}

#[test]
fn starter_mob_absorbs_its_constituent_jc() {
    // Given each sample reports a MOB plus its explicitly linked constituent JC.
    let mut starter_mob = GdEntry::mob(2, "chr1", 100, "IS1", "+", 3);
    starter_mob.parent_ids = vec![1];
    let starter = gd(vec![
        GdEntry::jc(1, "chr1", 95, "+", "chr2", 205, "-", 0),
        starter_mob,
    ]);
    let mut edited_mob = GdEntry::mob(12, "chr1", 103, "IS1", "+", 3);
    edited_mob.parent_ids = vec![11];
    let edited = gd(vec![
        GdEntry::jc(11, "chr1", 98, "+", "chr2", 202, "-", 0),
        edited_mob,
    ]);

    // When the canonical differential set is built.
    let result = build_differential_events(&starter, &edited, &refs()).expect("valid MOB lineage");

    // Then both samples absorb constituent JCs before subtraction.
    assert_eq!(result.starter_mutation_count, 1);
    assert!(result.events.is_empty());
}

#[test]
fn dangling_parent_ids_remain_attached_to_the_differential_event() {
    // Given a mutation referencing evidence records that are absent from its GD document.
    let mut snp = GdEntry::snp(1, "chr1", 100, "T");
    snp.parent_ids = vec![900, 800, 900];

    // When the differential event is built.
    let result = build_differential_events(&gd(vec![]), &gd(vec![snp]), &refs())
        .expect("valid mutation with dangling lineage");

    // Then unresolved lineage is preserved deterministically without inferred evidence.
    assert_eq!(
        result.events[0].evidence.unresolved_parent_ids,
        vec![800, 900]
    );
    assert_eq!(result.events[0].evidence.ra, None);
    assert_eq!(result.events[0].evidence.mc, None);
    assert_eq!(result.events[0].evidence.jc, None);
}

#[test]
fn full_input_permutations_produce_the_same_differential_events() {
    // Given an ambiguous structural graph plus an edited-only SNP.
    let starter = vec![
        GdEntry::jc(1, "chr1", 10, "+", "chr2", 10, "-", 0),
        GdEntry::jc(2, "chr1", 15, "+", "chr2", 15, "-", 0),
    ];
    let edited = vec![
        GdEntry::jc(3, "chr1", 11, "+", "chr2", 11, "-", 0),
        GdEntry::jc(4, "chr1", 12, "+", "chr2", 5, "-", 0),
        GdEntry::snp(5, "chr1", 300, "T"),
    ];
    let starter_orders = [
        starter.clone(),
        vec![starter[1].clone(), starter[0].clone()],
    ];
    let edited_orders = [
        edited.clone(),
        vec![edited[0].clone(), edited[2].clone(), edited[1].clone()],
        vec![edited[1].clone(), edited[0].clone(), edited[2].clone()],
        vec![edited[1].clone(), edited[2].clone(), edited[0].clone()],
        vec![edited[2].clone(), edited[0].clone(), edited[1].clone()],
        vec![edited[2].clone(), edited[1].clone(), edited[0].clone()],
    ];

    // When every input ordering is evaluated.
    let baseline = event_ids(starter_orders[0].clone(), edited_orders[0].clone());
    for starter_order in starter_orders {
        for edited_order in &edited_orders {
            // Then the canonical differential set is invariant to both input orders.
            assert_eq!(
                event_ids(starter_order.clone(), edited_order.clone()),
                baseline
            );
        }
    }
}

#[test]
fn con_is_exact_only_and_rejects_tolerance() {
    // CON events with 1 bp coordinate jitter must NOT be tolerantly subtracted.
    let starter = vec![GdEntry::con(1, "chr1", 100, 10, "chr1:200-209")];
    let edited = vec![GdEntry::con(2, "chr1", 101, 10, "chr1:200-209")];
    let surviving = event_ids(starter, edited);
    assert_eq!(surviving.len(), 1);
}

#[test]
fn jc_tolerance_accepts_four_and_five_but_rejects_six_bp() {
    let starter = vec![GdEntry::jc(1, "chr1", 100, "+", "chr2", 200, "-", 0)];

    // Positive and negative coordinate jitter at 4, 5, and 6 bp
    assert!(event_ids(
        starter.clone(),
        vec![GdEntry::jc(2, "chr1", 104, "+", "chr2", 204, "-", 0)]
    )
    .is_empty());
    assert!(event_ids(
        starter.clone(),
        vec![GdEntry::jc(3, "chr1", 105, "+", "chr2", 205, "-", 0)]
    )
    .is_empty());
    assert_eq!(
        event_ids(
            starter.clone(),
            vec![GdEntry::jc(4, "chr1", 106, "+", "chr2", 200, "-", 0)]
        )
        .len(),
        1
    );
    assert!(event_ids(
        starter.clone(),
        vec![GdEntry::jc(5, "chr1", 96, "+", "chr2", 196, "-", 0)]
    )
    .is_empty());
    assert!(event_ids(
        starter.clone(),
        vec![GdEntry::jc(6, "chr1", 95, "+", "chr2", 195, "-", 0)]
    )
    .is_empty());
    assert_eq!(
        event_ids(
            starter,
            vec![GdEntry::jc(7, "chr1", 94, "+", "chr2", 200, "-", 0)]
        )
        .len(),
        1
    );
}

#[test]
fn mob_tolerance_accepts_four_and_five_but_rejects_six_bp() {
    let starter = vec![GdEntry::mob(1, "chr1", 100, "IS1", "+", 3)];

    assert!(event_ids(
        starter.clone(),
        vec![GdEntry::mob(2, "chr1", 104, "IS1", "+", 3)]
    )
    .is_empty());
    assert!(event_ids(
        starter.clone(),
        vec![GdEntry::mob(3, "chr1", 105, "IS1", "+", 3)]
    )
    .is_empty());
    assert_eq!(
        event_ids(
            starter.clone(),
            vec![GdEntry::mob(4, "chr1", 106, "IS1", "+", 3)]
        )
        .len(),
        1
    );
    assert!(event_ids(
        starter.clone(),
        vec![GdEntry::mob(5, "chr1", 96, "IS1", "+", 3)]
    )
    .is_empty());
    assert!(event_ids(
        starter.clone(),
        vec![GdEntry::mob(6, "chr1", 95, "IS1", "+", 3)]
    )
    .is_empty());
    assert_eq!(
        event_ids(starter, vec![GdEntry::mob(7, "chr1", 94, "IS1", "+", 3)]).len(),
        1
    );
}

#[test]
fn duplicate_edited_events_merge_source_ids_and_evidence() {
    let mut ra1 = GdEntry::ra(10, "chr1", 100, 0, "A", "T");
    ra1.attrs.insert("pd_ra_support_reads".into(), "10".into());
    let mut ra2 = GdEntry::ra(20, "chr1", 100, 0, "A", "T");
    ra2.attrs.insert("pd_ra_support_reads".into(), "15".into());

    let mut snp1 = GdEntry::snp(1, "chr1", 100, "T");
    snp1.parent_ids = vec![10];
    let mut snp2 = GdEntry::snp(2, "chr1", 100, "T");
    snp2.parent_ids = vec![20];

    let edited = gd(vec![ra1, ra2, snp1, snp2]);
    let result = build_differential_events(&gd(vec![]), &edited, &refs()).expect("valid inputs");

    assert_eq!(result.events.len(), 1);
    let ev = &result.events[0];
    assert_eq!(ev.merged_source_ids, vec![1, 2]);
    assert_eq!(ev.evidence.ra.as_ref().map(|v| v.len()), Some(2));
}

#[test]
fn deterministic_tie_breaking_favors_smaller_canonical_bytes() {
    // Two starter candidates with identical coordinate cost to an edited candidate.
    // Starter A is at 99 (delta 1, max_delta 1, sum_delta 1)
    // Starter B is at 101 (delta 1, max_delta 1, sum_delta 1)
    // Both have identical cost (1, 1).
    let starter_a = GdEntry::mob(1, "chr1", 99, "IS1", "+", 3);
    let starter_b = GdEntry::mob(2, "chr1", 101, "IS1", "+", 3);
    let edited = GdEntry::mob(3, "chr1", 100, "IS1", "+", 3);

    // Run with starter in [A, B] order and [B, A] order.
    let res1 = build_differential_events(
        &gd(vec![starter_a.clone(), starter_b.clone()]),
        &gd(vec![edited.clone()]),
        &refs(),
    )
    .expect("valid");
    let res2 =
        build_differential_events(&gd(vec![starter_b, starter_a]), &gd(vec![edited]), &refs())
            .expect("valid");

    assert!(res1.events.is_empty());
    assert!(res2.events.is_empty());
}

#[test]
fn non_jc_record_with_pd_support_reads_does_not_set_jc_evidence() {
    use crate::audit::EvidenceSummary;

    let mut snp = GdEntry::snp(1, "chr1", 100, "T");
    snp.attrs.insert("pd_support_reads".into(), "25".into());

    let result = build_differential_events(&gd(vec![]), &gd(vec![snp]), &refs()).unwrap();
    let ev = &result.events[0];
    let summary = EvidenceSummary::from_evidence_references(&ev.evidence);
    // Sol finding: pd_support_reads on non-JC must NOT set jc=Some(true)
    assert_eq!(summary.jc, None);
}

#[test]
fn contig_boundary_checks_ins_and_del() {
    // refs() contigs have length 1_000
    // INS at position 1000 is valid (after last base)
    let ins_at_end = GdEntry::ins(1, "chr1", 1000, "A");
    assert!(CanonicalEvent::from_gd_entry(&ins_at_end, &refs()).is_ok());

    // INS at position 1001 is invalid (beyond contig length)
    let ins_beyond_end = GdEntry::ins(2, "chr1", 1001, "A");
    assert!(matches!(
        CanonicalEvent::from_gd_entry(&ins_beyond_end, &refs()),
        Err(DifferentialError::InvalidPosition { .. })
    ));

    // DEL ending exactly at 1000 (position 995, size 6 -> 995..=1000) is valid
    let del_at_end = GdEntry::del(3, "chr1", 995, 6);
    assert!(CanonicalEvent::from_gd_entry(&del_at_end, &refs()).is_ok());

    // DEL ending at 1001 (position 995, size 7 -> 995..=1001) is invalid
    let del_beyond_end = GdEntry::del(4, "chr1", 995, 7);
    assert!(matches!(
        CanonicalEvent::from_gd_entry(&del_beyond_end, &refs()),
        Err(DifferentialError::InvalidSize { .. })
    ));
}

#[test]
fn snp_and_one_base_sub_produce_distinct_canonical_event_ids_under_m1_contract() {
    // Under the documented M1 contract (docs/schema.md: "SUB 不拆分为 SNP/INS/DEL，不去除公共前缀后缀"),
    // a 1-base SUB (e.g. SUB chr1 100 1 T) and a SNP (SNP chr1 100 T) have separate physical
    // preimage encodings (kind tag 0x02 vs 0x01) and therefore produce distinct DV1 EventIds.
    let snp = GdEntry::snp(1, "chr1", 100, "T");
    let sub = GdEntry::sub(2, "chr1", 100, 1, "T");

    let snp_canonical = CanonicalEvent::from_gd_entry(&snp, &refs()).expect("valid SNP");
    let sub_canonical = CanonicalEvent::from_gd_entry(&sub, &refs()).expect("valid 1-base SUB");

    let (snp_id, snp_preimage) = snp_canonical.compute_event_id();
    let (sub_id, sub_preimage) = sub_canonical.compute_event_id();

    assert_ne!(snp_id, sub_id);
    assert_ne!(snp_preimage, sub_preimage);
    assert_eq!(snp_preimage[14], 0x01); // kind SNP
    assert_eq!(sub_preimage[14], 0x02); // kind SUB
}

#[test]
fn sub_does_not_trim_common_prefix_or_suffix() {
    // Under schema.md, SUB preserves its exact coordinates and replacement sequence
    // without stripping common prefix or suffix bases.
    let sub = GdEntry::sub(1, "chr1", 100, 3, "ATG");
    let canonical = CanonicalEvent::from_gd_entry(&sub, &refs()).expect("valid SUB");
    match canonical {
        CanonicalEvent::Sub {
            position,
            size,
            new_seq,
            ..
        } => {
            assert_eq!(position, 100);
            assert_eq!(size, 3);
            assert_eq!(new_seq, b"ATG");
        }
        _ => panic!("expected Sub canonical event"),
    }
}
