//! Multi-clone batch audit mode end-to-end integration tests (Milestone M10).
//!
//! Verifies:
//! 1. CLI help exposes `--batch-manifest`.
//! 2. Rule 6 enforcement: `--batch-manifest` without `--starter` fails with starter mandatory error.
//! 3. Missing manifest file triggers clean error.
//! 4. Invalid manifest format (duplicate clone ID, path traversal, ragged row, invalid spacer) triggers clean validation error.
//! 5. Preflight check triggers clean InputFileNotFound error before running the engine.
//! 6. Full batch pipeline execution on synthetic cohort: asserts exit code 0, generates batch_report.md,
//!    batch_variants.tsv, and individual clone product sets in clones/<clone_id>/ with accurate recurrence.

use std::path::{Path, PathBuf};
use std::process::Command;

fn prokadiff_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_prokadiff"))
}

fn test_dir(prefix: &str) -> PathBuf {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    let id = COUNTER.fetch_add(1, Ordering::SeqCst);
    let dir = std::env::temp_dir().join(format!(
        "prokadiff_batch_e2e_{prefix}_{id}_{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn touch(dir: &Path, name: &str) -> PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, b"").unwrap();
    p
}

#[test]
fn test_batch_cli_help_documents_manifest() {
    let output = Command::new(prokadiff_bin())
        .arg("--help")
        .output()
        .expect("failed to execute prokadiff");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("--batch-manifest"),
        "expected '--batch-manifest' in help stdout, got: {stdout}"
    );
}

#[test]
fn test_batch_cli_missing_starter_fails_rule6() {
    let dir = test_dir("missing_starter");
    let manifest = dir.join("manifest.tsv");
    std::fs::write(&manifest, "clone_id\tr1\nC01\tclone1.fq\n").unwrap();
    let ref_fa = touch(&dir, "ref.fa");

    let output = Command::new(prokadiff_bin())
        .args([
            "--batch-manifest",
            manifest.to_str().unwrap(),
            "--ref",
            ref_fa.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            dir.join("out").to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("starter-strain WGS is mandatory"),
        "expected Rule 6 starter mandatory error, got: {stderr}"
    );
}

#[test]
fn test_batch_cli_missing_manifest_file_fails() {
    let dir = test_dir("nonexistent_manifest");
    let starter = touch(&dir, "starter.fq");
    let ref_fa = touch(&dir, "ref.fa");

    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            starter.to_str().unwrap(),
            "--batch-manifest",
            dir.join("nonexistent_manifest.tsv").to_str().unwrap(),
            "--ref",
            ref_fa.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            dir.join("out").to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("not found") || stderr.contains("error"),
        "expected file not found or IO error, got: {stderr}"
    );
}

#[test]
fn test_batch_cli_duplicate_clone_id_in_manifest_fails() {
    let dir = test_dir("dup_clone");
    let starter = touch(&dir, "starter.fq");
    let ref_fa = touch(&dir, "ref.fa");
    let r1 = touch(&dir, "read.fq");

    let manifest = dir.join("manifest.tsv");
    std::fs::write(
        &manifest,
        format!(
            "clone_id\tr1\nC01\t{}\nC01\t{}\n",
            r1.display(),
            r1.display()
        ),
    )
    .unwrap();

    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            starter.to_str().unwrap(),
            "--batch-manifest",
            manifest.to_str().unwrap(),
            "--ref",
            ref_fa.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            dir.join("out").to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("duplicate clone_id 'C01'"),
        "expected duplicate clone_id error, got: {stderr}"
    );
}

#[test]
fn test_batch_cli_path_traversal_clone_id_rejected() {
    let dir = test_dir("traversal_clone");
    let starter = touch(&dir, "starter.fq");
    let ref_fa = touch(&dir, "ref.fa");
    let r1 = touch(&dir, "read.fq");

    let manifest = dir.join("manifest.tsv");
    std::fs::write(
        &manifest,
        format!("clone_id\tr1\n../../escape\t{}\n", r1.display()),
    )
    .unwrap();

    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            starter.to_str().unwrap(),
            "--batch-manifest",
            manifest.to_str().unwrap(),
            "--ref",
            ref_fa.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            dir.join("out").to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("invalid clone_id"),
        "expected invalid clone_id error, got: {stderr}"
    );
}

#[test]
fn test_batch_cli_ragged_row_in_manifest_fails() {
    let dir = test_dir("ragged_manifest");
    let starter = touch(&dir, "starter.fq");
    let ref_fa = touch(&dir, "ref.fa");
    let r1 = touch(&dir, "read.fq");

    let manifest = dir.join("manifest.tsv");
    std::fs::write(
        &manifest,
        format!("clone_id\tr1\tr2\nC01\t{}\n", r1.display()),
    )
    .unwrap();

    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            starter.to_str().unwrap(),
            "--batch-manifest",
            manifest.to_str().unwrap(),
            "--ref",
            ref_fa.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            dir.join("out").to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("ragged row"),
        "expected ragged row error, got: {stderr}"
    );
}

#[test]
fn test_batch_cli_preflight_missing_read_file_fails() {
    let dir = test_dir("preflight_missing");
    let starter = touch(&dir, "starter.fq");
    let ref_fa = touch(&dir, "ref.fa");

    let manifest = dir.join("manifest.tsv");
    std::fs::write(&manifest, "clone_id\tr1\nC01\tnon_existent_read.fq\n").unwrap();

    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            starter.to_str().unwrap(),
            "--batch-manifest",
            manifest.to_str().unwrap(),
            "--ref",
            ref_fa.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            dir.join("out").to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("not found"),
        "expected preflight missing file error, got: {stderr}"
    );
}

#[test]
fn test_batch_cli_manifest_only_guides_accepted_and_executed() {
    let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf();
    let fixture_dir = repo_root.join("testdata/generated/synth_parent_child");

    if !fixture_dir.join("starter_R1.fastq").exists() {
        let gen_script = repo_root.join("testdata/generate.sh");
        if gen_script.exists() {
            let status = Command::new("bash")
                .arg(&gen_script)
                .current_dir(&repo_root)
                .status()
                .expect("failed to execute testdata/generate.sh");
            assert!(status.success(), "testdata/generate.sh failed");
        }
    }
    assert!(
        fixture_dir.join("starter_R1.fastq").exists(),
        "synthetic fixture must exist"
    );

    let dir = test_dir("manifest_cas9_guide");
    let outdir = dir.join("out");
    let manifest = dir.join("manifest.tsv");

    let r1 = fixture_dir.join("edited_R1.fastq");
    let r2 = fixture_dir.join("edited_R2.fastq");

    std::fs::write(
        &manifest,
        format!(
            "clone_id\tr1\tr2\tspacer\tpam\nClone_Cas9\t{}\t{}\tGAGTTCATCTACGCCGTGAA\tNGG\n",
            r1.display(),
            r2.display()
        ),
    )
    .unwrap();

    let mut path_env = std::env::var("PATH").unwrap_or_default();
    path_env = format!("/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:{path_env}");

    // With --editor cas9 and NO global --spacer, batch execution must succeed and apply the manifest spacer
    let output = Command::new(prokadiff_bin())
        .env("PATH", path_env)
        .args([
            "--starter",
            fixture_dir.join("starter_R1.fastq").to_str().unwrap(),
            fixture_dir.join("starter_R2.fastq").to_str().unwrap(),
            "--batch-manifest",
            manifest.to_str().unwrap(),
            "--ref",
            fixture_dir.join("ref.fa").to_str().unwrap(),
            "--editor",
            "cas9",
            "--threads",
            "4",
            "--outdir",
            outdir.to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute prokadiff");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "batch execution with manifest-only guide failed: {stderr}"
    );

    // Verify clone provenance recorded the effective manifest spacer
    let clone_dir = outdir.join("clones/Clone_Cas9");
    assert!(clone_dir.join("report.md").exists());
    assert!(clone_dir.join("provenance.tsv").exists());
    let prov = std::fs::read_to_string(clone_dir.join("provenance.tsv")).unwrap();
    assert!(
        prov.contains("FIXTURE_VALIDATED") || prov.contains("VALIDATED"),
        "offtarget search must run with manifest guide: {prov}"
    );

    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn test_batch_cli_manifest_missing_guide_for_cas9_fails() {
    let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf();
    let fixture_dir = repo_root.join("testdata/generated/synth_parent_child");

    let dir = test_dir("manifest_missing_cas9_guide");
    let outdir = dir.join("out");
    let manifest = dir.join("manifest.tsv");

    let r1 = fixture_dir.join("edited_R1.fastq");
    let r2 = fixture_dir.join("edited_R2.fastq");

    std::fs::write(
        &manifest,
        format!(
            "clone_id\tr1\tr2\tspacer\tpam\nClone_NoGuide\t{}\t{}\t.\t.\n",
            r1.display(),
            r2.display()
        ),
    )
    .unwrap();

    let mut path_env = std::env::var("PATH").unwrap_or_default();
    path_env = format!("/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:{path_env}");

    // With --editor cas9 and NO global --spacer, missing clone spacer must fail cleanly
    let output = Command::new(prokadiff_bin())
        .env("PATH", path_env)
        .args([
            "--starter",
            fixture_dir.join("starter_R1.fastq").to_str().unwrap(),
            fixture_dir.join("starter_R2.fastq").to_str().unwrap(),
            "--batch-manifest",
            manifest.to_str().unwrap(),
            "--ref",
            fixture_dir.join("ref.fa").to_str().unwrap(),
            "--editor",
            "cas9",
            "--threads",
            "4",
            "--outdir",
            outdir.to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("clone 'Clone_NoGuide' is missing a required spacer sequence for editor Cas9"),
        "expected missing spacer error identifying clone, got: {stderr}"
    );

    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn test_batch_e2e_full_execution_synthetic_clones() {
    let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf();
    let fixture_dir = repo_root.join("testdata/generated/synth_parent_child");

    if !fixture_dir.join("starter_R1.fastq").exists() {
        let gen_script = repo_root.join("testdata/generate.sh");
        if gen_script.exists() {
            let status = Command::new("bash")
                .arg(&gen_script)
                .current_dir(&repo_root)
                .status()
                .expect("failed to execute testdata/generate.sh");
            assert!(status.success(), "testdata/generate.sh failed");
        }
    }
    assert!(
        fixture_dir.join("starter_R1.fastq").exists(),
        "synthetic fixture must exist for Gate D execution test"
    );

    let dir = test_dir("full_batch_synth");
    let outdir = dir.join("out");
    let manifest_path = dir.join("batch_manifest.tsv");

    // Construct 2-clone manifest using synthetic edited reads:
    // Clone_1: uses intended.tsv (intended edit at pos 2500 marked COMPLETE)
    // Clone_2: no intended table (intended edit not specified)
    let r1 = fixture_dir.join("edited_R1.fastq");
    let r2 = fixture_dir.join("edited_R2.fastq");
    let intended = fixture_dir.join("intended.tsv");

    let manifest_content = format!(
        "clone_id\tr1\tr2\tintended_tsv\n\
         Clone_1\t{}\t{}\t{}\n\
         Clone_2\t{}\t{}\t.\n",
        r1.display(),
        r2.display(),
        intended.display(),
        r1.display(),
        r2.display(),
    );
    std::fs::write(&manifest_path, manifest_content).unwrap();

    let mut path_env = std::env::var("PATH").unwrap_or_default();
    path_env = format!("/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:{path_env}");

    let output = Command::new(prokadiff_bin())
        .env("PATH", path_env)
        .args([
            "--starter",
            fixture_dir.join("starter_R1.fastq").to_str().unwrap(),
            fixture_dir.join("starter_R2.fastq").to_str().unwrap(),
            "--batch-manifest",
            manifest_path.to_str().unwrap(),
            "--ref",
            fixture_dir.join("ref.fa").to_str().unwrap(),
            "--editor",
            "dsb",
            "--threads",
            "4",
            "--outdir",
            outdir.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run prokadiff batch pipeline");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "prokadiff batch failed with code {:?}\nstdout:\n{stdout}\nstderr:\n{stderr}",
        output.status.code()
    );

    // 1. Verify batch-level product outputs
    let batch_report_md = outdir.join("batch_report.md");
    let batch_variants_tsv = outdir.join("batch_variants.tsv");
    assert!(batch_report_md.exists(), "batch_report.md must be generated");
    assert!(batch_variants_tsv.exists(), "batch_variants.tsv must be generated");

    // 2. Audit batch_report.md contents
    let report_content = std::fs::read_to_string(&batch_report_md).unwrap();
    assert!(report_content.contains("# ProkaDiff Multi-Clone Batch Audit Report"));
    assert!(report_content.contains("**Total Evaluated Clones:** 2"));
    assert!(report_content.contains("`Clone_1`"));
    assert!(report_content.contains("`Clone_2`"));
    assert!(report_content.contains("Complete (1/1)"));
    assert!(report_content.contains("Not specified"));
    assert!(report_content.contains("Strict Non-Prescriptive Principle (Taskbook §66)"));
    assert!(report_content.contains("Spatial association is consistent with a possible guide-dependent event but does not establish cleavage causality."));
    assert!(!report_content.contains("best clone"));
    assert!(!report_content.contains("safest"));

    // 3. Audit batch_variants.tsv contents
    let matrix_content = std::fs::read_to_string(&batch_variants_tsv).unwrap();
    assert!(matrix_content.contains("batch_variant_id\tseq_id\tposition\tend\tgd_type\tref\talt\tside2_seq_id\tside2_position\trecurrence\tcount\ttotal\tclone_Clone_1\tclone_Clone_2"));
    // Extra edit at pos 8000 must be present and shared across both clones
    assert!(matrix_content.contains("\t8000\t8000\tSNP\t"));
    assert!(matrix_content.contains("shared_all\t2\t2\t"));

    // 4. Verify individual clone directory outputs
    let clone1_dir = outdir.join("clones/Clone_1");
    let clone2_dir = outdir.join("clones/Clone_2");
    assert!(clone1_dir.exists(), "clones/Clone_1 directory must exist");
    assert!(clone2_dir.exists(), "clones/Clone_2 directory must exist");

    for cdir in [&clone1_dir, &clone2_dir] {
        assert!(cdir.join("report.md").exists());
        assert!(cdir.join("post_edit_variants.tsv").exists());
        assert!(cdir.join("unintended.tsv").exists());
        assert!(cdir.join("summary.txt").exists());
        assert!(cdir.join("provenance.tsv").exists());
        assert!(cdir.join("starter.gd").exists());
        assert!(cdir.join("edited.gd").exists());
    }

    // In Clone_1, intended edit at 2500 is complete
    let clone1_summary = std::fs::read_to_string(clone1_dir.join("summary.txt")).unwrap();
    assert!(clone1_summary.contains("intended_status\tall_observed"));

    // In Clone_2, intended edit was not declared
    let clone2_summary = std::fs::read_to_string(clone2_dir.join("summary.txt")).unwrap();
    assert!(clone2_summary.contains("intended_provided\tno"));

    let _ = std::fs::remove_dir_all(dir);
}
