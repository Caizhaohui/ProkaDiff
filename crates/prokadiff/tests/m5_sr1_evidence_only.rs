use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use noodles::sam;
use prokadiff_classify::{associate_events_to_sites, build_differential_events, RefContig};
use prokadiff_evidence::align::sam_to_sorted_bam;
use prokadiff_evidence::engine::bam_io::read_primary_bam;
use prokadiff_evidence::pileup::{AlignedRead, CigarKind, CigarOp};
use prokadiff_evidence::repeat_ambiguous::{
    normalized_molecule_id, RepeatAmbiguousDiagnostics, RepeatAmbiguousState,
};
use prokadiff_evidence::split_seed::{
    find_candidate_junctions_with_repeat_evidence, parse_sam_split_record, SplitReadAlignment,
};
use prokadiff_evidence::{
    call_from_aligned_with_diagnostics, call_from_aligned_with_extra_and_diagnostics,
    EngineOptions, FastaRecord,
};
use prokadiff_gd::{GdKind, GenomeDiff};

fn repeated_clip_input(copies: usize) -> (Vec<FastaRecord>, Vec<AlignedRead>) {
    let motif = b"ACGTTGCATGACCTGATCGA";
    let mut sequence = vec![b'T'; copies * 30 + 40];
    for index in 0..copies {
        let start = index * 30;
        sequence[start..start + motif.len()].copy_from_slice(motif);
    }
    let anchor_start = (copies * 30) as i64;
    let read_sequence = [vec![b'T'; 20], motif.to_vec()].concat();
    let reads = [false, false, false, true, true, true]
        .into_iter()
        .enumerate()
        .map(|(molecule_id, minus)| AlignedRead {
            contig_idx: 0,
            ref_start_0: anchor_start,
            minus,
            seq: read_sequence.clone(),
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: 20,
                },
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: motif.len(),
                },
            ],
            mapq: 42,
            molecule_id: molecule_id as u64,
        })
        .collect();
    (
        vec![FastaRecord {
            name: "chr".into(),
            seq: sequence,
        }],
        reads,
    )
}

fn assert_zero_product_isolation(gd: &GenomeDiff, fasta: &[FastaRecord]) {
    assert_eq!(
        gd.entries
            .iter()
            .filter(|entry| matches!(entry.kind, GdKind::Jc))
            .count(),
        0,
        "zero exact JC required"
    );
    assert_eq!(
        gd.entries
            .iter()
            .filter(|entry| matches!(entry.kind, GdKind::Del))
            .count(),
        0,
        "zero DEL required"
    );
    assert_eq!(
        gd.entries
            .iter()
            .filter(|entry| matches!(entry.kind, GdKind::Mob))
            .count(),
        0,
        "zero MOB required"
    );

    let references = vec![RefContig {
        name: fasta[0].name.clone(),
        seq: fasta[0].seq.clone(),
    }];
    let result = build_differential_events(&GenomeDiff::default(), gd, &references)
        .expect("evidence-only records are valid differential input");

    assert_eq!(result.events.len(), 0, "zero DifferentialEvent required");
    assert!(result.events.is_empty(), "zero EventId required");

    let associations = associate_events_to_sites(&result.events, &[], 50);
    assert_eq!(associations.len(), 0, "zero association required");
}

fn assert_evidence_only(copies: usize, expected_state: RepeatAmbiguousState) {
    let (fasta, reads) = repeated_clip_input(copies);
    let (gd, diagnostics) =
        call_from_aligned_with_diagnostics(&fasta, &reads, &EngineOptions::default());

    // 1. Confirm state matches expected_state
    assert_eq!(diagnostics.records.len(), 1);
    assert_eq!(diagnostics.records[0].state, expected_state);

    // 2. From the SAME fixture, assert zero direct downstream product
    assert_zero_product_isolation(&gd, &fasta);
}

#[test]
fn supported_ambiguous_repeat_evidence_has_no_product_event() {
    assert_evidence_only(21, RepeatAmbiguousState::SupportedAmbiguous);
}

#[test]
fn resource_limited_repeat_evidence_has_no_product_event() {
    assert_evidence_only(4_097, RepeatAmbiguousState::ResourceLimit);
}

#[test]
fn stage_two_overflow_repeat_evidence_has_no_product_event() {
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'A'; 200_000],
    }];
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
    let split_res = find_candidate_junctions_with_repeat_evidence(&alignments);
    assert!(split_res.candidates.is_empty());
    assert!(!split_res.repeat_seeds.is_empty());

    // Exercise the SAME overflow repeat evidence through the downstream-capable path
    // with actual aligned reads:
    let aligned_reads = vec![
        AlignedRead {
            contig_idx: 0,
            ref_start_0: 100,
            minus: false,
            seq: vec![b'A'; 40],
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: 40,
            }],
            mapq: 42,
            molecule_id: 1,
        },
        AlignedRead {
            contig_idx: 0,
            ref_start_0: 100,
            minus: true,
            seq: vec![b'A'; 40],
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: 40,
            }],
            mapq: 42,
            molecule_id: 2,
        },
    ];

    let (gd, diagnostics) = call_from_aligned_with_extra_and_diagnostics(
        &fasta,
        &aligned_reads,
        &EngineOptions::default(),
        &[],
        split_res.repeat_seeds,
    );

    // 1. Diagnostics from the SAME fixture must confirm RESOURCE_LIMIT
    assert_eq!(diagnostics.records.len(), 1);
    assert_eq!(
        diagnostics.records[0].state,
        RepeatAmbiguousState::ResourceLimit
    );
    assert!(!diagnostics.records[0].exact_breakpoint_supported);

    // 2. From the SAME fixture, assert zero direct downstream product
    assert_zero_product_isolation(&gd, &fasta);
}

fn resolved_reciprocal_clip_input() -> (Vec<FastaRecord>, Vec<AlignedRead>) {
    let motif_a = b"ACGTTGCATGACCTGATCGA";
    let motif_b = b"TGCAGTGCTACGATGACGTA";
    let copies = 21;
    let mut sequence = vec![b'T'; 3000];
    for i in 0..copies {
        let start = i * 40;
        sequence[start..start + 20].copy_from_slice(motif_a);
    }
    for j in 0..copies {
        let start = 1000 + j * 40;
        sequence[start..start + 20].copy_from_slice(motif_b);
    }

    let forward_seq = [motif_a.to_vec(), motif_b.to_vec()].concat();
    let reciprocal_seq = [motif_a.to_vec(), motif_b.to_vec()].concat();

    let mut reads = Vec::new();
    // 6 forward reads (3 plus strand, 3 minus strand) anchored at copy 5 of motif_a (pos 200)
    for (idx, minus) in [false, false, false, true, true, true]
        .into_iter()
        .enumerate()
    {
        reads.push(AlignedRead {
            contig_idx: 0,
            ref_start_0: 200,
            minus,
            seq: forward_seq.clone(),
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
            molecule_id: (idx + 1) as u64,
        });
    }

    // 1 reciprocal read (molecule 99) anchored at copy 3 of motif_b (pos 1120)
    reads.push(AlignedRead {
        contig_idx: 0,
        ref_start_0: 1120,
        minus: false,
        seq: reciprocal_seq,
        cigar: vec![
            CigarOp {
                kind: CigarKind::SoftClip,
                len: 20,
            },
            CigarOp {
                kind: CigarKind::Match,
                len: 20,
            },
        ],
        mapq: 42,
        molecule_id: 99,
    });

    (
        vec![FastaRecord {
            name: "chr".into(),
            seq: sequence,
        }],
        reads,
    )
}

#[test]
fn resolved_repeat_evidence_has_no_product_event() {
    let (fasta, reads) = resolved_reciprocal_clip_input();
    let (gd, diagnostics) =
        call_from_aligned_with_diagnostics(&fasta, &reads, &EngineOptions::default());

    // 1. Genuinely reaches RepeatAmbiguousState::Resolved through the real parser/evidence path
    let resolved_record = diagnostics
        .records
        .iter()
        .find(|record| record.state == RepeatAmbiguousState::Resolved)
        .expect("genuine RepeatAmbiguousState::Resolved record required");

    assert!(resolved_record.exact_breakpoint_supported);
    assert!(resolved_record.reciprocal_evidence);
    assert_eq!(resolved_record.feasible_copies.len(), 1);

    // 2. From that SAME logical fixture, assert zero direct downstream products:
    assert_zero_product_isolation(&gd, &fasta);
}

fn truncated_reciprocal_clip_input() -> (Vec<FastaRecord>, Vec<AlignedRead>) {
    let motif_a = b"ACGTTGCATGACCTGATCGA";
    let motif_b = b"TGCAGTGCTACGATGACGTA";
    let copies_a = 4_097; // exceeds MAX_RETAINED_REPEAT_PLACEMENTS (4096)
    let copies_b = 21;
    let mut sequence = vec![b'T'; copies_a * 25 + 2000];
    for i in 0..copies_a {
        let start = i * 25;
        sequence[start..start + 20].copy_from_slice(motif_a);
    }
    let offset_b = copies_a * 25 + 100;
    for j in 0..copies_b {
        let start = offset_b + j * 40;
        sequence[start..start + 20].copy_from_slice(motif_b);
    }

    let forward_seq = [motif_b.to_vec(), motif_a.to_vec()].concat();
    let reciprocal_seq = [motif_b.to_vec(), motif_a.to_vec()].concat();

    let mut reads = Vec::new();
    // Forward reads anchored at motif_b copy 0, with softclip motif_a (which has 4,097 copies)
    for (idx, minus) in [false, false, false, true, true, true]
        .into_iter()
        .enumerate()
    {
        reads.push(AlignedRead {
            contig_idx: 0,
            ref_start_0: offset_b as i64,
            minus,
            seq: forward_seq.clone(),
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
            molecule_id: (idx + 1) as u64,
        });
    }

    // Reciprocal read anchored at copy 5 of motif_a (pos 5 * 25 = 125, within the retained 4096 subset),
    // matching motif_a and softclipped with motif_b.
    // If placement_set_complete were true, this reciprocal read would select copy 5.
    reads.push(AlignedRead {
        contig_idx: 0,
        ref_start_0: 125,
        minus: false,
        seq: reciprocal_seq,
        cigar: vec![
            CigarOp {
                kind: CigarKind::SoftClip,
                len: 20,
            },
            CigarOp {
                kind: CigarKind::Match,
                len: 20,
            },
        ],
        mapq: 42,
        molecule_id: 99,
    });

    (
        vec![FastaRecord {
            name: "chr".into(),
            seq: sequence,
        }],
        reads,
    )
}

#[test]
fn truncation_exceeding_bound_with_reciprocal_evidence_cannot_resolve() {
    let (fasta, reads) = truncated_reciprocal_clip_input();
    let (gd, diagnostics) =
        call_from_aligned_with_diagnostics(&fasta, &reads, &EngineOptions::default());

    // Record for the truncated family (placement_count 4097 > 4096)
    let truncated_record = diagnostics
        .records
        .iter()
        .find(|record| record.placement_count > 4096)
        .expect("truncated family record required");

    // Assert:
    // state != RESOLVED (must be ResourceLimit)
    assert_ne!(truncated_record.state, RepeatAmbiguousState::Resolved);
    assert_eq!(truncated_record.state, RepeatAmbiguousState::ResourceLimit);
    // exact_breakpoint_supported = false
    assert!(!truncated_record.exact_breakpoint_supported);
    // placement_set_complete is false
    assert!(!truncated_record.placement_set_complete);
    // eligible partial copies cannot become a singleton exact call
    assert!(truncated_record.feasible_copies.is_empty());

    // and zero direct downstream product from the repeat diagnostic path
    assert_zero_product_isolation(&gd, &fasta);
}

#[test]
fn real_parser_reciprocal_index_resource_limit_has_no_product_event() {
    let motif = b"ACGTTGCATGACCTGATCGA";
    let copies = 4_001; // 4001 <= 4096 so placement_set_complete = true and resource_complete = true
    let match_anchor = b"CCCCCCCCCCCCCCCCCCCC";
    let mut sequence = vec![b'T'; copies * 25 + 2000];
    sequence[500..520].copy_from_slice(match_anchor);
    for index in 0..copies {
        let start = 1000 + index * 25;
        sequence[start..start + motif.len()].copy_from_slice(motif);
    }
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: sequence,
    }];

    // Generate SAM input with 51 reads
    let read_seq = [match_anchor.to_vec(), motif.to_vec()].concat();
    let read_seq_str = std::str::from_utf8(&read_seq).unwrap();

    let mut sam_content = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr\tLN:{}\n",
        copies * 25 + 2000
    );
    for mol in 0..51 {
        sam_content.push_str(&format!(
            "read_{mol}\t0\tchr\t501\t42\t20M20S\t*\t0\t0\t{read_seq_str}\t*\n"
        ));
    }

    let temp_dir = std::env::temp_dir();
    let nonce = format!("{}_{}", std::process::id(), 1);
    let sam_path = temp_dir.join(format!("prokadiff_recip_index_{nonce}.sam"));
    let bam_path = temp_dir.join(format!("prokadiff_recip_index_{nonce}.bam"));

    std::fs::write(&sam_path, &sam_content).unwrap();
    sam_to_sorted_bam(&sam_path, &bam_path).unwrap();

    let bam_data = read_primary_bam(&bam_path, &fasta).unwrap();
    let _ = std::fs::remove_file(&sam_path);
    let _ = std::fs::remove_file(&bam_path);

    let (gd, diagnostics) =
        call_from_aligned_with_diagnostics(&fasta, &bam_data.aligned, &EngineOptions::default());

    let limited_record = diagnostics
        .records
        .iter()
        .find(|record| {
            record
                .rejection_reasons
                .contains(&"RESOURCE_LIMIT_RECIPROCAL_INDEX")
        })
        .expect("must contain record reaching RESOURCE_LIMIT_RECIPROCAL_INDEX");

    assert_eq!(limited_record.state, RepeatAmbiguousState::ResourceLimit);
    assert!(!limited_record.exact_breakpoint_supported);
    assert!(limited_record.feasible_copies.is_empty());

    assert_zero_product_isolation(&gd, &fasta);
}

#[test]
fn real_parser_reciprocal_candidates_resource_limit_has_no_product_event() {
    let motif_a = b"ACGTTGCATGACCTGATCGA";
    let motif_b = b"TGCAGTGCTACGATGACGTA";
    let cand_match = b"CCCCCCCCCCCCCCCCCCCC";
    let copies = 21;
    let mut sequence = vec![b'T'; 35_000];
    for i in 0..copies {
        let start = i * 40;
        sequence[start..start + 20].copy_from_slice(motif_a);
    }
    for j in 0..copies {
        let start = 1000 + j * 40;
        sequence[start..start + 20].copy_from_slice(motif_b);
    }
    sequence[5000..5020].copy_from_slice(cand_match);
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: sequence,
    }];

    // Record 0: forward read anchored at copy 5 of motif_a (pos 200..220), softclip motif_b
    // Anchor position: 220, minus: true
    let fwd_seq = [motif_a.to_vec(), motif_b.to_vec()].concat();
    let fwd_seq_str = std::str::from_utf8(&fwd_seq).unwrap();

    // 5,001 candidate reads: anchored at 5000 (matching cand_match), softclip motif_a on left (20S20M)
    // Left clip motif_a placed at hit 200 has position 200 + 20 = 220, minus: true (exact match to forward anchor!)
    let cand_seq = [motif_a.to_vec(), cand_match.to_vec()].concat();
    let cand_seq_str = std::str::from_utf8(&cand_seq).unwrap();

    let mut sam_content = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr\tLN:35000\nread_fwd\t0\tchr\t201\t42\t20M20S\t*\t0\t0\t{fwd_seq_str}\t*\n"
    );
    for i in 1..=5001 {
        sam_content.push_str(&format!(
            "read_{i}\t0\tchr\t5001\t42\t20S20M\t*\t0\t0\t{cand_seq_str}\t*\n"
        ));
    }

    let temp_dir = std::env::temp_dir();
    let nonce = format!("{}_{}", std::process::id(), 2);
    let sam_path = temp_dir.join(format!("prokadiff_recip_cand_{nonce}.sam"));
    let bam_path = temp_dir.join(format!("prokadiff_recip_cand_{nonce}.bam"));

    std::fs::write(&sam_path, &sam_content).unwrap();
    sam_to_sorted_bam(&sam_path, &bam_path).unwrap();

    let bam_data = read_primary_bam(&bam_path, &fasta).unwrap();
    let _ = std::fs::remove_file(&sam_path);
    let _ = std::fs::remove_file(&bam_path);

    let (gd, diagnostics) =
        call_from_aligned_with_diagnostics(&fasta, &bam_data.aligned, &EngineOptions::default());

    let limited_record = diagnostics
        .records
        .iter()
        .find(|record| {
            record
                .rejection_reasons
                .contains(&"RESOURCE_LIMIT_RECIPROCAL_CANDIDATES")
        })
        .expect("must contain record reaching RESOURCE_LIMIT_RECIPROCAL_CANDIDATES");

    assert_eq!(limited_record.state, RepeatAmbiguousState::ResourceLimit);
    assert!(!limited_record.exact_breakpoint_supported);
    assert!(limited_record.feasible_copies.is_empty());

    assert_zero_product_isolation(&gd, &fasta);
}

#[test]
fn primary_bam_and_stage2_sam_cross_path_molecule_consistency() {
    let fasta = vec![FastaRecord {
        name: "chr1".into(),
        seq: vec![b'A'; 20_000],
    }];

    // SAM records covering:
    // 1. FIRST_SEGMENT flag (65) vs /1 suffix (0)
    // 2. LAST_SEGMENT flag (129) vs /2 suffix (0)
    // 3. no flags and no suffix -> UNKNOWN mate (0)
    // 4. duplicate alignments of same molecule
    let sam_text = "\
@HD\tVN:1.6\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:20000\n\
mol_a\t65\tchr1\t100\t42\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
mol_a/1\t0\tchr1\t100\t42\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
mol_b\t129\tchr1\t200\t42\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
mol_b/2\t0\tchr1\t200\t42\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
mol_plain\t0\tchr1\t300\t42\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
mol_dup\t65\tchr1\t400\t42\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n\
mol_dup\t65\tchr1\t450\t42\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n";

    let temp_dir = std::env::temp_dir();
    let nonce = format!("{}_{}", std::process::id(), 3);
    let sam_path = temp_dir.join(format!("prokadiff_cross_path_{nonce}.sam"));
    let bam_path = temp_dir.join(format!("prokadiff_cross_path_{nonce}.bam"));

    std::fs::write(&sam_path, sam_text).unwrap();
    sam_to_sorted_bam(&sam_path, &bam_path).unwrap();

    // 1. Primary BAM parser
    let bam_data = read_primary_bam(&bam_path, &fasta).unwrap();

    // 2. Stage-2 SAM parser
    let mut name_to_idx = HashMap::new();
    name_to_idx.insert("chr1", 0);
    let sam_file = File::open(&sam_path).unwrap();
    let mut sam_reader = sam::io::Reader::new(BufReader::new(sam_file));
    let _ = sam_reader.read_header().unwrap();
    let mut stage2_alignments = Vec::new();
    for rec in sam_reader.records() {
        let rec = rec.unwrap();
        if let Some(aln) = parse_sam_split_record(&rec, &name_to_idx) {
            stage2_alignments.push(aln);
        }
    }

    let _ = std::fs::remove_file(&sam_path);
    let _ = std::fs::remove_file(&bam_path);

    // Assert coverage of all cases:
    // First segment (flag 65) vs /1 suffix (flag 0)
    let prim_mol_a_flag = bam_data
        .aligned
        .iter()
        .find(|r| r.ref_start_0 == 99 && r.molecule_id == normalized_molecule_id("mol_a", 1))
        .expect("mol_a flag");
    let stage2_mol_a_flag = stage2_alignments
        .iter()
        .find(|a| a.qname == "mol_a" && a.mate == 1)
        .expect("stage2 mol_a flag");
    assert_eq!(
        prim_mol_a_flag.molecule_id,
        normalized_molecule_id(&stage2_mol_a_flag.qname, stage2_mol_a_flag.mate)
    );

    let prim_mol_a_suf = bam_data
        .aligned
        .iter()
        .find(|r| r.molecule_id == normalized_molecule_id("mol_a", 1))
        .expect("mol_a suf");
    let stage2_mol_a_suf = stage2_alignments
        .iter()
        .find(|a| a.qname == "mol_a" && a.mate == 1)
        .expect("stage2 mol_a suf");
    assert_eq!(
        prim_mol_a_suf.molecule_id,
        normalized_molecule_id(&stage2_mol_a_suf.qname, stage2_mol_a_suf.mate)
    );
    assert_eq!(prim_mol_a_flag.molecule_id, prim_mol_a_suf.molecule_id);

    // Last segment (flag 129) vs /2 suffix (flag 0)
    let prim_mol_b_flag = bam_data
        .aligned
        .iter()
        .find(|r| r.molecule_id == normalized_molecule_id("mol_b", 2))
        .expect("mol_b flag");
    let stage2_mol_b_flag = stage2_alignments
        .iter()
        .find(|a| a.qname == "mol_b" && a.mate == 2)
        .expect("stage2 mol_b flag");
    assert_eq!(
        prim_mol_b_flag.molecule_id,
        normalized_molecule_id(&stage2_mol_b_flag.qname, stage2_mol_b_flag.mate)
    );

    let prim_mol_b_suf = bam_data
        .aligned
        .iter()
        .find(|r| r.molecule_id == normalized_molecule_id("mol_b", 2))
        .expect("mol_b suf");
    let stage2_mol_b_suf = stage2_alignments
        .iter()
        .find(|a| a.qname == "mol_b" && a.mate == 2)
        .expect("stage2 mol_b suf");
    assert_eq!(
        prim_mol_b_suf.molecule_id,
        normalized_molecule_id(&stage2_mol_b_suf.qname, stage2_mol_b_suf.mate)
    );
    assert_eq!(prim_mol_b_flag.molecule_id, prim_mol_b_suf.molecule_id);

    // Plain read (no flags, no suffix -> UNKNOWN mate)
    let prim_mol_plain = bam_data
        .aligned
        .iter()
        .find(|r| r.molecule_id == normalized_molecule_id("mol_plain", 0))
        .expect("mol_plain");
    let stage2_mol_plain = stage2_alignments
        .iter()
        .find(|a| a.qname == "mol_plain" && a.mate == 0)
        .expect("stage2 mol_plain");
    assert_eq!(
        prim_mol_plain.molecule_id,
        normalized_molecule_id(&stage2_mol_plain.qname, stage2_mol_plain.mate)
    );

    // Duplicate alignments: two alignments share identical molecule identity
    let prim_dups: Vec<_> = bam_data
        .aligned
        .iter()
        .filter(|r| r.molecule_id == normalized_molecule_id("mol_dup", 1))
        .collect();
    assert_eq!(prim_dups.len(), 2);
    assert_eq!(prim_dups[0].molecule_id, prim_dups[1].molecule_id);

    let stage2_dups: Vec<_> = stage2_alignments
        .iter()
        .filter(|a| a.qname == "mol_dup" && a.mate == 1)
        .collect();
    assert_eq!(stage2_dups.len(), 2);
    assert_eq!(
        normalized_molecule_id(&stage2_dups[0].qname, stage2_dups[0].mate),
        normalized_molecule_id(&stage2_dups[1].qname, stage2_dups[1].mate)
    );
    assert_eq!(
        prim_dups[0].molecule_id,
        normalized_molecule_id(&stage2_dups[0].qname, stage2_dups[0].mate)
    );

    // No inflation of independent molecule support:
    // Feed duplicate alignments into RepeatAmbiguousDiagnostics
    let seed1 = prokadiff_evidence::repeat_ambiguous::RepeatAmbiguousSeed {
        source: prokadiff_evidence::repeat_ambiguous::RepeatEvidenceSource::PrimaryClip,
        anchor_contig_idx: 0,
        anchor_position_1: 400,
        anchor_minus: false,
        observed_placements: vec![prokadiff_evidence::repeat_ambiguous::RepeatCopyPlacement {
            contig_idx: 0,
            position_1: 800,
            minus: false,
        }],
        eligible_copies: vec![prokadiff_evidence::repeat_ambiguous::RepeatCopyPlacement {
            contig_idx: 0,
            position_1: 800,
            minus: false,
        }],
        placement_count: 21,
        placement_family: b"FAMILY_DUP".to_vec(),
        placement_sequence: b"FAMILY_DUP".to_vec(),
        placement_set_complete: true,
        resource_complete: true,
        resolution_complete: false,
        overlap: 0,
        molecule_id: prim_dups[0].molecule_id,
        molecule_minus: false,
        clip_length: 20,
        aligned_length: 20,
        effective_mapq: 42,
        alignment_score: Some(40),
        pair_geometry: "PRIMARY_SOFTCLIP",
        unique_anchor_qualified: true,
        copy_resolved: false,
        resolution_molecule_id: None,
        reciprocal_evidence: false,
        rejection_reason: None,
        resource_limit_reason: None,
    };
    let mut seed2 = seed1.clone();
    seed2.anchor_position_1 = 450;
    let diag = RepeatAmbiguousDiagnostics::from_seeds(&fasta, vec![seed1, seed2]);
    for rec in &diag.records {
        assert_eq!(
            rec.molecule_ids.len(),
            1,
            "duplicate alignments must not inflate molecule count"
        );
        assert_eq!(rec.plus_molecules, 1, "support must not inflate");
    }

    // Reordered input: reverse order produces identical canonical molecule identity
    let sam_lines: Vec<&str> = sam_text.lines().collect();
    let header_lines = &sam_lines[..2];
    let body_lines: Vec<&str> = sam_lines[2..].iter().rev().copied().collect();
    let mut rev_sam = header_lines.join("\n");
    rev_sam.push('\n');
    rev_sam.push_str(&body_lines.join("\n"));
    rev_sam.push('\n');

    let nonce_rev = format!("{}_{}", std::process::id(), 4);
    let sam_rev_path = temp_dir.join(format!("prokadiff_cross_path_rev_{nonce_rev}.sam"));
    let bam_rev_path = temp_dir.join(format!("prokadiff_cross_path_rev_{nonce_rev}.bam"));

    std::fs::write(&sam_rev_path, rev_sam).unwrap();
    sam_to_sorted_bam(&sam_rev_path, &bam_rev_path).unwrap();

    let bam_rev_data = read_primary_bam(&bam_rev_path, &fasta).unwrap();
    let _ = std::fs::remove_file(&sam_rev_path);
    let _ = std::fs::remove_file(&bam_rev_path);

    let mut ids_orig: Vec<u64> = bam_data.aligned.iter().map(|r| r.molecule_id).collect();
    let mut ids_rev: Vec<u64> = bam_rev_data.aligned.iter().map(|r| r.molecule_id).collect();
    ids_orig.sort();
    ids_rev.sort();
    assert_eq!(
        ids_orig, ids_rev,
        "reordered input must yield identical canonical molecule identities"
    );
}
