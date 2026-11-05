//! Synthetic Titration Test Suite for Subclonal / Pool-seq Variant Audit Mode.
//!
//! Milestone M11 Task M11-T05:
//! Verifies:
//! 1. In-silico titration gradient: 50%, 25%, 10%, 5%, and 0% (negative control).
//! 2. Precision and false-positive control: 0% negative controls never emit subclonal variants.
//! 3. Strand bias flagging on extreme asymmetric support.
//! 4. Starter subtraction logic: pre-existing background vs. edited enrichment.
//! 5. Depth gating: starter depth < 20 yields INDETERMINATE_STARTER_DEPTH.
//! 6. CLI integration: argument parsing, help text, and parameter validation.

use std::path::PathBuf;
use std::process::Command;

use prokadiff_classify::subclonal::{
    evaluate_subclonal_candidates, RawSubclonalObservation,
};
use prokadiff_evidence::subclonal::{AlleleCounts, StrandCounts, SubclonalClass, SubclonalParams};

fn prokadiff_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_prokadiff"))
}

#[test]
fn test_synthetic_titration_gradient_50_25_10_05_00() {
    let params = SubclonalParams {
        min_depth: 20,
        min_alt_reads: 3,
        min_af: 0.05,
        max_af: 0.7499,
        min_delta_af: 0.05,
        max_strand_bias_pvalue: 0.001,
        fdr_alpha: 0.01,
        error_floor_snp: 0.005,
        error_floor_indel: 0.010,
    };

    // 200x depth simulation across titration gradient:
    // 50% (100 reads), 25% (50 reads), 10% (20 reads), 5% (10 reads), 0% (0 reads).
    let candidates = vec![
        // Site 1: 50% AF in edited (100/200), 0% in starter
        RawSubclonalObservation {
            contig: "chr1".to_string(),
            pos: 1000,
            ref_seq: "A".to_string(),
            alt_seq: "G".to_string(),
            variant_type: "SNP".to_string(),
            starter_counts: AlleleCounts::new(
                StrandCounts::new(100, 100),
                StrandCounts::new(0, 0),
                0,
            ),
            edited_counts: AlleleCounts::new(
                StrandCounts::new(50, 50),
                StrandCounts::new(50, 50),
                0,
            ),
            starter_error_rate: 0.005,
            edited_error_rate: 0.005,
        },
        // Site 2: 25% AF in edited (50/200), 0% in starter
        RawSubclonalObservation {
            contig: "chr1".to_string(),
            pos: 2000,
            ref_seq: "C".to_string(),
            alt_seq: "T".to_string(),
            variant_type: "SNP".to_string(),
            starter_counts: AlleleCounts::new(
                StrandCounts::new(100, 100),
                StrandCounts::new(0, 0),
                0,
            ),
            edited_counts: AlleleCounts::new(
                StrandCounts::new(75, 75),
                StrandCounts::new(25, 25),
                0,
            ),
            starter_error_rate: 0.005,
            edited_error_rate: 0.005,
        },
        // Site 3: 10% AF in edited (20/200), 0% in starter
        RawSubclonalObservation {
            contig: "chr1".to_string(),
            pos: 3000,
            ref_seq: "G".to_string(),
            alt_seq: "A".to_string(),
            variant_type: "SNP".to_string(),
            starter_counts: AlleleCounts::new(
                StrandCounts::new(100, 100),
                StrandCounts::new(0, 0),
                0,
            ),
            edited_counts: AlleleCounts::new(
                StrandCounts::new(90, 90),
                StrandCounts::new(10, 10),
                0,
            ),
            starter_error_rate: 0.005,
            edited_error_rate: 0.005,
        },
        // Site 4: 5% AF in edited (10/200), 0% in starter (at the 5% boundary)
        RawSubclonalObservation {
            contig: "chr1".to_string(),
            pos: 4000,
            ref_seq: "T".to_string(),
            alt_seq: "C".to_string(),
            variant_type: "SNP".to_string(),
            starter_counts: AlleleCounts::new(
                StrandCounts::new(100, 100),
                StrandCounts::new(0, 0),
                0,
            ),
            edited_counts: AlleleCounts::new(
                StrandCounts::new(95, 95),
                StrandCounts::new(5, 5),
                0,
            ),
            starter_error_rate: 0.005,
            edited_error_rate: 0.005,
        },
        // Site 5: 0% negative control (pure wild-type)
        RawSubclonalObservation {
            contig: "chr1".to_string(),
            pos: 5000,
            ref_seq: "A".to_string(),
            alt_seq: "T".to_string(),
            variant_type: "SNP".to_string(),
            starter_counts: AlleleCounts::new(
                StrandCounts::new(100, 100),
                StrandCounts::new(0, 0),
                0,
            ),
            edited_counts: AlleleCounts::new(
                StrandCounts::new(100, 100),
                StrandCounts::new(0, 0),
                0,
            ),
            starter_error_rate: 0.005,
            edited_error_rate: 0.005,
        },
    ];

    let results = evaluate_subclonal_candidates(candidates, &params);
    // 4 true subclonal sites pass; 0% negative control (Site 5) is NotSubclonal and strictly excluded
    assert_eq!(results.len(), 4);
    assert!(!results.iter().any(|r| r.pos == 5000));

    // Site 1 (50%)
    let v50 = &results[0];
    assert_eq!(v50.pos, 1000);
    assert_eq!(v50.classification, SubclonalClass::EditedSpecific);
    assert!((v50.edited_af.unwrap() - 0.50).abs() < 1e-4);
    assert!(v50.delta_ci_low > 0.35);

    // Site 2 (25%)
    let v25 = &results[1];
    assert_eq!(v25.pos, 2000);
    assert_eq!(v25.classification, SubclonalClass::EditedSpecific);
    assert!((v25.edited_af.unwrap() - 0.25).abs() < 1e-4);
    assert!(v25.delta_ci_low > 0.15);

    // Site 3 (10%)
    let v10 = &results[2];
    assert_eq!(v10.pos, 3000);
    assert_eq!(v10.classification, SubclonalClass::EditedSpecific);
    assert!((v10.edited_af.unwrap() - 0.10).abs() < 1e-4);
    assert!(v10.delta_ci_low > 0.03);

    // Site 4 (5%)
    let v05 = &results[3];
    assert_eq!(v05.pos, 4000);
    assert_eq!(v05.classification, SubclonalClass::EditedSpecific);
    assert!((v05.edited_af.unwrap() - 0.05).abs() < 1e-4);
    assert!(v05.delta_ci_low > 0.0);
}

#[test]
fn test_strand_bias_detection_and_flagging() {
    let params = SubclonalParams::default();

    // 12% AF with 12 reads all on the Forward strand (0 on reverse strand), total depth 100
    // Fisher exact test p-value: 0.000973 < 0.001
    let obs = RawSubclonalObservation {
        contig: "chr1".to_string(),
        pos: 1234,
        ref_seq: "G".to_string(),
        alt_seq: "C".to_string(),
        variant_type: "SNP".to_string(),
        starter_counts: AlleleCounts::new(
            StrandCounts::new(50, 50),
            StrandCounts::new(0, 0),
            0,
        ),
        edited_counts: AlleleCounts::new(
            StrandCounts::new(44, 44),
            StrandCounts::new(12, 0), // 12 on fwd, 0 on rev -> severe strand bias!
            0,
        ),
        starter_error_rate: 0.005,
        edited_error_rate: 0.005,
    };

    let results = evaluate_subclonal_candidates(vec![obs], &params);
    assert_eq!(results.len(), 1);
    let v = &results[0];
    assert!(v.edited_strand_bias_p < params.max_strand_bias_pvalue);
    assert!(v.flags.contains(&"STRAND_BIAS".to_string()));
    assert!(v.flags.contains(&"SINGLE_STRAND_ALT".to_string()));
}

#[test]
fn test_starter_inheritance_and_enrichment() {
    let params = SubclonalParams::default();

    // Case A: Pre-existing at 10% in both starter and edited
    let obs_inherited = RawSubclonalObservation {
        contig: "chr1".to_string(),
        pos: 100,
        ref_seq: "A".to_string(),
        alt_seq: "T".to_string(),
        variant_type: "SNP".to_string(),
        starter_counts: AlleleCounts::new(
            StrandCounts::new(45, 45),
            StrandCounts::new(5, 5),
            0,
        ),
        edited_counts: AlleleCounts::new(
            StrandCounts::new(45, 45),
            StrandCounts::new(5, 5),
            0,
        ),
        starter_error_rate: 0.005,
        edited_error_rate: 0.005,
    };

    // Case B: Pre-existing at 10% in starter, enriched to 50% in edited
    let obs_enriched = RawSubclonalObservation {
        contig: "chr1".to_string(),
        pos: 200,
        ref_seq: "C".to_string(),
        alt_seq: "G".to_string(),
        variant_type: "SNP".to_string(),
        starter_counts: AlleleCounts::new(
            StrandCounts::new(45, 45),
            StrandCounts::new(5, 5),
            0,
        ),
        edited_counts: AlleleCounts::new(
            StrandCounts::new(25, 25),
            StrandCounts::new(25, 25),
            0,
        ),
        starter_error_rate: 0.005,
        edited_error_rate: 0.005,
    };

    let results = evaluate_subclonal_candidates(vec![obs_inherited, obs_enriched], &params);
    // Un-enriched starter background (pos 100) is subtracted and excluded; enriched allele (pos 200) is reported
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].pos, 200);
    assert_eq!(results[0].classification, SubclonalClass::EditedEnriched);
}

#[test]
fn test_indeterminate_starter_depth_gating() {
    let params = SubclonalParams::default();

    // Edited has 20% AF with depth 100, but starter depth is only 15 (< 20)
    let obs = RawSubclonalObservation {
        contig: "chr1".to_string(),
        pos: 100,
        ref_seq: "G".to_string(),
        alt_seq: "A".to_string(),
        variant_type: "SNP".to_string(),
        starter_counts: AlleleCounts::new(
            StrandCounts::new(8, 7), // depth 15
            StrandCounts::new(0, 0),
            0,
        ),
        edited_counts: AlleleCounts::new(
            StrandCounts::new(40, 40),
            StrandCounts::new(10, 10),
            0,
        ),
        starter_error_rate: 0.005,
        edited_error_rate: 0.005,
    };

    let results = evaluate_subclonal_candidates(vec![obs], &params);
    assert_eq!(results[0].classification, SubclonalClass::IndeterminateStarterDepth);
}

#[test]
fn test_cli_help_documents_subclonal_flags() {
    let output = Command::new(prokadiff_bin())
        .arg("--help")
        .output()
        .expect("run prokadiff --help");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--subclonal"));
    assert!(stdout.contains("--min-subclonal-af"));
    assert!(stdout.contains("--min-subclonal-alt-reads"));
    assert!(stdout.contains("--max-strand-bias-pvalue"));
}

#[test]
fn test_subclonal_cli_e2e_run() {
    let repo_root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf();
    let fixture_dir = repo_root.join("testdata/generated/synth_parent_child");

    let outdir = std::env::temp_dir().join(format!("prokadiff_subclonal_e2e_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&outdir);

    let mut path_env = std::env::var("PATH").unwrap_or_default();
    path_env = format!("/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:{path_env}");

    let output = Command::new(prokadiff_bin())
        .env("PATH", path_env)
        .arg("--starter")
        .arg(fixture_dir.join("starter_R1.fastq"))
        .arg(fixture_dir.join("starter_R2.fastq"))
        .arg("--edited")
        .arg(fixture_dir.join("edited_R1.fastq"))
        .arg(fixture_dir.join("edited_R2.fastq"))
        .arg("--ref")
        .arg(fixture_dir.join("ref.fa"))
        .arg("--editor")
        .arg("dsb")
        .arg("--outdir")
        .arg(&outdir)
        .arg("--subclonal")
        .output()
        .expect("run prokadiff with --subclonal");

    assert!(output.status.success(), "STDERR: {}", String::from_utf8_lossy(&output.stderr));

    // Verify subclonal_audit.tsv exists
    let tsv_path = outdir.join("subclonal_audit.tsv");
    assert!(tsv_path.is_file(), "subclonal_audit.tsv must be generated");
    let tsv_content = std::fs::read_to_string(&tsv_path).expect("read subclonal_audit.tsv");
    assert!(tsv_content.starts_with("contig\tposition\tvariant_type\tref\talt\t"));

    // Verify report.md contains mandatory caveat block
    let report_path = outdir.join("report.md");
    assert!(report_path.is_file(), "report.md must be generated");
    let report_content = std::fs::read_to_string(&report_path).expect("read report.md");
    assert!(report_content.contains("Subclonal / Pool-seq Variant Audit (Target AF 0.05–0.70)"));
    assert!(report_content.contains("Subclonal audit mode reports read-level ALT fractions at sites outside the consensus threshold."));

    let _ = std::fs::remove_dir_all(&outdir);
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
        .collect()
}

#[test]
fn test_production_read_titration_gradient_with_negative_control() {
    let test_dir = std::env::temp_dir().join(format!("prokadiff_titration_e2e_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&test_dir);
    std::fs::create_dir_all(&test_dir).unwrap();

    // 1. Construct 5000 bp synthetic reference contig (pseudorandom non-repeating sequence)
    let mut rng = 123456789u64;
    let mut ref_bytes = Vec::with_capacity(5000);
    let bases = b"ACGT";
    for _ in 0..5000 {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
        ref_bytes.push(bases[((rng >> 33) & 3) as usize]);
    }

    // Set known reference bases at target loci (1-based positions)
    ref_bytes[999] = b'A';  // pos 1000: REF 'A', ALT 'G' (50% in edited)
    ref_bytes[1999] = b'C'; // pos 2000: REF 'C', ALT 'T' (25% in edited)
    ref_bytes[2999] = b'G'; // pos 3000: REF 'G', ALT 'A' (10% in edited)
    ref_bytes[3999] = b'T'; // pos 4000: REF 'T', ALT 'C' (5% in edited)
    ref_bytes[4499] = b'A'; // pos 4500: REF 'A' (0% negative control)

    // Pin indel/INS locus contexts (0-based indices): distinct neighbours so
    // right-alignment cannot slide the emitted allele, and no homopolymer at
    // the anchor so repeat masking cannot hide the base observation.
    ref_bytes[1198..1202].copy_from_slice(b"ACTG"); // INS anchor 1199='C'; ref[1200]='T' != 'G' -> no INS slide
    ref_bytes[2498..2502].copy_from_slice(b"TACG"); // 1 bp DEL anchor 2499='A'; ref[2500]='C' -> no slide
    ref_bytes[3498..3502].copy_from_slice(b"GTAC"); // 2 bp DEL anchor 3499..3501="TA"; ref[3501]='C' -> no slide

    let ref_path = test_dir.join("ref.fa");
    std::fs::write(&ref_path, format!(">chr1\n{}\n", String::from_utf8_lossy(&ref_bytes))).unwrap();

    // 2. Generate titrated reads covering target loci (the titrated INS and
    //    DEL loci exercise the production INS/DEL counting paths end to end)
    #[derive(Clone, Copy)]
    enum TitrationVariant {
        Snp { alt_b: u8 },
        Del { len: usize },
        Ins { oligo: &'static [u8] },
    }

    let loci: Vec<(usize, TitrationVariant, usize, usize)> = vec![
        (1000, TitrationVariant::Snp { alt_b: b'G' }, 50, 100),  // 50% AF
        (1200, TitrationVariant::Ins { oligo: b"GT" }, 15, 100), // 15% AF
        (2000, TitrationVariant::Snp { alt_b: b'T' }, 25, 100),  // 25% AF
        (2500, TitrationVariant::Del { len: 1 }, 25, 100),       // 25% AF
        (3000, TitrationVariant::Snp { alt_b: b'A' }, 10, 100),  // 10% AF
        (3500, TitrationVariant::Del { len: 2 }, 20, 100),       // 20% AF
        (4000, TitrationVariant::Snp { alt_b: b'C' }, 10, 200),  // 5% AF (200 fragments)
        (4500, TitrationVariant::Snp { alt_b: b'T' }, 0, 100),   // 0% negative control
    ];

    let mut starter_r1 = Vec::new();
    let mut starter_r2 = Vec::new();
    let mut edited_r1 = Vec::new();
    let mut edited_r2 = Vec::new();

    let mut read_id = 0u64;

    // Write one PE fragment (R1 + revcomp R2) into the given FASTQ buffers.
    // R1 = fragment[0..70]; R2 = revcomp(fragment[10..min(80, len)]).
    let push_pair = |buf1: &mut Vec<u8>, buf2: &mut Vec<u8>, id: u64, fragment: &[u8]| {
        let r1_end = 70.min(fragment.len());
        let r2_start = 10.min(fragment.len().saturating_sub(1));
        let r2_end = 80.min(fragment.len());
        let r1 = &fragment[0..r1_end];
        let r2 = revcomp(&fragment[r2_start..r2_end]);
        let qual1 = "I".repeat(r1.len());
        let qual2 = "I".repeat(r2.len());
        buf1.extend_from_slice(
            format!("@read_{id}/1\n{}\n+\n{qual1}\n", String::from_utf8_lossy(r1)).as_bytes(),
        );
        buf2.extend_from_slice(
            format!("@read_{id}/2\n{}\n+\n{qual2}\n", String::from_utf8_lossy(&r2)).as_bytes(),
        );
    };

    for &(pos, variant, n_alt_edited, n_total) in &loci {
        let start_idx = pos - 40;
        let pos0 = pos - 1;
        let ref_fragment = ref_bytes[start_idx..start_idx + 80].to_vec();

        // ALT fragment for this locus (REF reads reuse ref_fragment).
        let alt_fragment = match variant {
            TitrationVariant::Snp { alt_b, .. } => {
                let mut fragment = ref_fragment.clone();
                fragment[39] = alt_b; // mutate 1-based pos (offset 39 from start_idx)
                fragment
            }
            TitrationVariant::Del { len } => {
                // Omit the deleted reference base(s) from the fragment center.
                let mut fragment = ref_bytes[start_idx..pos0].to_vec();
                fragment.extend_from_slice(&ref_bytes[pos0 + len..start_idx + 80]);
                fragment
            }
            TitrationVariant::Ins { oligo } => {
                // Insert the extra base(s) after the anchor reference base.
                let mut fragment = ref_bytes[start_idx..=pos0].to_vec();
                fragment.extend_from_slice(oligo);
                fragment.extend_from_slice(&ref_bytes[pos0 + 1..start_idx + 80]);
                fragment
            }
        };

        // Starter reads (all REF)
        for _ in 0..n_total {
            read_id += 1;
            push_pair(&mut starter_r1, &mut starter_r2, read_id, &ref_fragment);
        }

        // Edited reads (titrated ALT reads + REF reads)
        for i in 0..n_total {
            read_id += 1;
            let fragment = if i < n_alt_edited { &alt_fragment } else { &ref_fragment };
            push_pair(&mut edited_r1, &mut edited_r2, read_id, fragment);
        }
    }

    let starter_r1_path = test_dir.join("starter_R1.fq");
    let starter_r2_path = test_dir.join("starter_R2.fq");
    let edited_r1_path = test_dir.join("edited_R1.fq");
    let edited_r2_path = test_dir.join("edited_R2.fq");

    std::fs::write(&starter_r1_path, starter_r1).unwrap();
    std::fs::write(&starter_r2_path, starter_r2).unwrap();
    std::fs::write(&edited_r1_path, edited_r1).unwrap();
    std::fs::write(&edited_r2_path, edited_r2).unwrap();

    let outdir = test_dir.join("out");

    let mut path_env = std::env::var("PATH").unwrap_or_default();
    path_env = format!("/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:{path_env}");

    // 3. Execute ProkaDiff with --subclonal
    let output = Command::new(prokadiff_bin())
        .env("PATH", path_env)
        .arg("--starter")
        .arg(&starter_r1_path)
        .arg(&starter_r2_path)
        .arg("--edited")
        .arg(&edited_r1_path)
        .arg(&edited_r2_path)
        .arg("--ref")
        .arg(&ref_path)
        .arg("--editor")
        .arg("dsb")
        .arg("--outdir")
        .arg(&outdir)
        .arg("--subclonal")
        .output()
        .expect("execute prokadiff --subclonal titration");

    assert!(output.status.success(), "STDERR: {}", String::from_utf8_lossy(&output.stderr));

    // 4. Verify subclonal_audit.tsv
    let tsv_path = outdir.join("subclonal_audit.tsv");
    assert!(tsv_path.is_file(), "subclonal_audit.tsv must be present");
    let tsv_content = std::fs::read_to_string(&tsv_path).unwrap();
    let data_lines: Vec<&str> = tsv_content
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .collect();

    // Exactly 7 titration variants must be reported; the 0% negative control
    // (pos 4500) must be absent. Depths are deterministic: every locus emits
    // n_total PE fragments and BOTH mates observe the locus, so column depth
    // is 2*n_total and ALT observations are 2*n_alt.
    assert_eq!(
        data_lines.len(),
        7,
        "expected exactly 7 reported subclonal variants; found {}: {:?}",
        data_lines.len(),
        data_lines
    );

    // Row lookup by exact (position, ref, alt) identity.
    let find_row = |target_pos: u64, target_ref: &str, target_alt: &str| -> Vec<&str> {
        let line = data_lines
            .iter()
            .find(|l| {
                let parts: Vec<&str> = l.split('\t').collect();
                parts.get(1).map(|p| p.parse::<u64>().unwrap_or(0)) == Some(target_pos)
                    && parts.get(3) == Some(&target_ref)
                    && parts.get(4) == Some(&target_alt)
            })
            .unwrap_or_else(|| {
                panic!("missing subclonal variant row for pos {target_pos} ref {target_ref} alt {target_alt}; rows: {data_lines:?}")
            });
        line.split('\t').collect()
    };

    // Locus assertion helper: exact depths/counts, AF within +/-0.02, class.
    let assert_locus = |row: &[&str],
                        expect_type: &str,
                        expect_starter_depth: u64,
                        expect_edited_depth: u64,
                        expect_edited_alt: u64,
                        expect_af: f64| {
        assert_eq!(row[2], expect_type, "variant_type at pos {}: {}", row[1], row[2]);
        assert_eq!(
            row[5].parse::<u64>().unwrap(),
            expect_starter_depth,
            "starter_depth at pos {}",
            row[1]
        );
        assert_eq!(row[6], "0", "starter_alt at pos {} must be 0", row[1]);
        assert_eq!(
            row[10].parse::<u64>().unwrap(),
            expect_edited_depth,
            "edited_depth at pos {}",
            row[1]
        );
        assert_eq!(
            row[11].parse::<u64>().unwrap(),
            expect_edited_alt,
            "edited_alt at pos {}",
            row[1]
        );
        let af: f64 = row[12].parse().unwrap();
        assert!(
            (af - expect_af).abs() <= 0.02,
            "AF at pos {} = {af}, expected ~{expect_af} (+/-0.02)",
            row[1]
        );
        assert_eq!(row[24], "EDITED_SPECIFIC", "classification at pos {}", row[1]);
    };

    // Site 1: 50% AF SNP (pos 1000, A -> G): 100 fragments x 2 mates.
    let row1 = find_row(1000, "A", "G");
    assert_locus(&row1, "SNP", 200, 200, 100, 0.50);

    // Site 2: 15% AF INS (pos 1200, insert "GT"): depth 200, ALT entries 30.
    let row2 = find_row(1200, "-", "GT");
    assert_locus(&row2, "INS", 200, 200, 30, 0.15);

    // Site 3: 25% AF SNP (pos 2000, C -> T).
    let row3 = find_row(2000, "C", "T");
    assert_locus(&row3, "SNP", 200, 200, 50, 0.25);

    // Site 4: 25% AF 1 bp DEL (pos 2500, ref 'A' -> '-'): span-1 runs only.
    let row4 = find_row(2500, "A", "-");
    assert_locus(&row4, "DEL", 200, 200, 50, 0.25);

    // Site 5: 10% AF SNP (pos 3000, G -> A).
    let row5 = find_row(3000, "G", "A");
    assert_locus(&row5, "SNP", 200, 200, 20, 0.10);

    // Site 6: 20% AF 2 bp DEL (pos 3500, ref "TA" -> '-'): span-2 runs only.
    let row6 = find_row(3500, "TA", "-");
    assert_locus(&row6, "DEL", 200, 200, 40, 0.20);

    // Site 7: 5% AF SNP (pos 4000, T -> C, 200 fragments => depth 400).
    let row7 = find_row(4000, "T", "C");
    assert_locus(&row7, "SNP", 400, 400, 20, 0.05);

    // Negative control (pos 4500): MUST be strictly absent
    assert!(
        !data_lines.iter().any(|l| l.split('\t').nth(1) == Some("4500")),
        "0% negative control at pos 4500 must not be emitted"
    );

    let _ = std::fs::remove_dir_all(&test_dir);
}
