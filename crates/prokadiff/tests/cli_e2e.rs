//! End-to-end (E2E) CLI automated integration tests.
//!
//! Verifies:
//! 1. `--version` and `--help` text and exit code 0.
//! 2. Mandatory arguments: missing `--starter`, `--edited`, `--ref`, `--outdir`.
//! 3. Editor validation: `cast`/`is110` roadmap errors, unknown editor errors.
//! 4. Spacer and PAM validation: missing spacer for cas9/cas12a, DNA alphabet enforcement.
//! 5. FASTQ pairing validation: odd number of FASTQ files rejected with friendly message.
//! 6. Input file existence: non-existent files reported with clean `FileNotFound` error.
//! 7. Single-sample `evidence` subcommand validation.
//! 8. Full E2E pipeline execution on `testdata/generated/synth_parent_child` fixture:
//!    asserting exit code 0, subtract correctness, intended masking, and output formats.

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
        "prokadiff_e2e_{prefix}_{id}_{}",
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
fn e2e_cli_version() {
    let output = Command::new(prokadiff_bin())
        .arg("--version")
        .output()
        .expect("failed to execute prokadiff");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("prokadiff 0.1.0"),
        "expected 'prokadiff 0.1.0' in stdout, got: {stdout}"
    );
}

#[test]
fn e2e_cli_help() {
    let output = Command::new(prokadiff_bin())
        .arg("--help")
        .output()
        .expect("failed to execute prokadiff");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("starter vs edited WGS"));
    assert!(stdout.contains("--starter"));
    assert!(stdout.contains("--edited"));
    assert!(stdout.contains("--ref"));
    assert!(stdout.contains("--intended"));
    assert!(stdout.contains("--editor"));
    assert!(stdout.contains("--spacer"));
    assert!(stdout.contains("--pam"));
    assert!(stdout.contains("--outdir"));
    assert!(stdout.contains("evidence"));
}

#[test]
fn e2e_cli_missing_starter_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "--edited", "e.fq", "--ref", "r.fa", "--editor", "dsb", "--outdir", "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("starter-strain WGS is mandatory"),
        "expected starter mandatory message, got: {stderr}"
    );
}

#[test]
fn e2e_cli_missing_edited_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            "s.fq",
            "--ref",
            "r.fa",
            "--editor",
            "dsb",
            "--outdir",
            "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("starter-strain WGS is mandatory"),
        "expected starter mandatory message, got: {stderr}"
    );
}

#[test]
fn e2e_cli_missing_ref_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            "s.fq",
            "--edited",
            "e.fq",
            "--editor",
            "dsb",
            "--outdir",
            "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--ref: at least one reference file"),
        "expected missing ref message, got: {stderr}"
    );
}

#[test]
fn e2e_cli_missing_outdir_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            "s.fq",
            "--edited",
            "e.fq",
            "--ref",
            "r.fa",
            "--editor",
            "dsb",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--outdir: output directory is required"),
        "expected missing outdir message, got: {stderr}"
    );
}

#[test]
fn e2e_cli_missing_spacer_for_cas9_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            "s.fq",
            "--edited",
            "e.fq",
            "--ref",
            "r.fa",
            "--editor",
            "cas9",
            "--outdir",
            "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--spacer is required for --editor cas9 and cas12a"),
        "expected missing spacer message, got: {stderr}"
    );
}

#[test]
fn e2e_cli_roadmap_editor_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            "s.fq",
            "--edited",
            "e.fq",
            "--ref",
            "r.fa",
            "--editor",
            "cast",
            "--outdir",
            "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("unsupported --editor 'cast'"),
        "expected unsupported cast editor, got: {stderr}"
    );
    assert!(
        stderr.contains("roadmap"),
        "expected roadmap reference, got: {stderr}"
    );
}

#[test]
fn e2e_cli_unknown_editor_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            "s.fq",
            "--edited",
            "e.fq",
            "--ref",
            "r.fa",
            "--editor",
            "alien_crispr",
            "--outdir",
            "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("unknown --editor 'alien_crispr'"),
        "expected unknown editor message, got: {stderr}"
    );
}

#[test]
fn e2e_cli_invalid_spacer_uracil_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            "s.fq",
            "--edited",
            "e.fq",
            "--ref",
            "r.fa",
            "--editor",
            "cas9",
            "--spacer",
            "AUGCAUGCAUGCAUGCAUGC",
            "--outdir",
            "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("use 'T' instead of 'U'"),
        "expected uracil warning, got: {stderr}"
    );
}

#[test]
fn e2e_cli_odd_fastq_count_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            "s1.fq",
            "--starter",
            "s2.fq",
            "--starter",
            "s3.fq",
            "--edited",
            "e.fq",
            "--ref",
            "r.fa",
            "--editor",
            "dsb",
            "--outdir",
            "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("expected 1 file (SE) or an even number of files"),
        "expected odd fastq count error, got: {stderr}"
    );
    assert!(
        stderr.contains("got 3 file(s)"),
        "expected file count in error, got: {stderr}"
    );
}

#[test]
fn e2e_cli_file_not_found_fails() {
    let d = test_dir("fnf_e2e");
    let s = touch(&d, "s.fq");
    let e = touch(&d, "e.fq");
    let missing_ref = d.join("non_existent_ref.fa");

    let output = Command::new(prokadiff_bin())
        .args([
            "--starter",
            s.to_str().unwrap(),
            "--edited",
            e.to_str().unwrap(),
            "--ref",
            missing_ref.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--ref: file") && stderr.contains("not found"),
        "expected file not found error, got: {stderr}"
    );
    let _ = std::fs::remove_dir_all(d);
}

#[test]
fn e2e_cli_evidence_help() {
    let output = Command::new(prokadiff_bin())
        .args(["evidence", "--help"])
        .output()
        .expect("failed to execute prokadiff");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Single-sample RA/MC/JC engine"));
    assert!(stdout.contains("--fastq"));
}

#[test]
fn e2e_cli_evidence_missing_file_fails() {
    let output = Command::new(prokadiff_bin())
        .args([
            "evidence",
            "--ref",
            "non_existent_ref_123.fa",
            "--fastq",
            "test.fq",
            "--outdir",
            "out",
        ])
        .output()
        .expect("failed to execute prokadiff");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--ref: file") && stderr.contains("not found"),
        "expected file not found error, got: {stderr}"
    );
}

#[test]
fn e2e_pipeline_synth_parent_child() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fixture_dir = manifest_dir.join("../../testdata/generated/synth_parent_child");
    if !fixture_dir.join("starter_R1.fastq").exists() {
        eprintln!("Skipping e2e_pipeline_synth_parent_child: fixture not generated yet");
        return;
    }

    // Ensure bowtie2 is in PATH (check conda env if needed)
    let conda_bin = PathBuf::from("/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin");
    let mut current_path = std::env::var("PATH").unwrap_or_default();
    if conda_bin.join("bowtie2").exists() {
        current_path = format!("{}:{}", conda_bin.display(), current_path);
    }

    // Check if bowtie2 is available
    let has_bowtie2 = Command::new("bowtie2")
        .arg("--version")
        .env("PATH", &current_path)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);

    if !has_bowtie2 {
        eprintln!("Skipping e2e_pipeline_synth_parent_child: bowtie2 not found on PATH");
        return;
    }

    let outdir = test_dir("synth_parent_child_pipeline");
    let output = Command::new(prokadiff_bin())
        .env("PATH", &current_path)
        .args([
            "--starter",
            fixture_dir.join("starter_R1.fastq").to_str().unwrap(),
            "--starter",
            fixture_dir.join("starter_R2.fastq").to_str().unwrap(),
            "--edited",
            fixture_dir.join("edited_R1.fastq").to_str().unwrap(),
            "--edited",
            fixture_dir.join("edited_R2.fastq").to_str().unwrap(),
            "--ref",
            fixture_dir.join("ref.fa").to_str().unwrap(),
            "--intended",
            fixture_dir.join("intended.tsv").to_str().unwrap(),
            "--editor",
            "dsb",
            "--threads",
            "4",
            "--outdir",
            outdir.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run prokadiff pipeline");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "prokadiff failed with code {:?}\nstdout:\n{stdout}\nstderr:\n{stderr}",
        output.status.code()
    );

    let unintended_tsv = outdir.join("unintended.tsv");
    let summary_txt = outdir.join("summary.txt");
    let starter_gd = outdir.join("starter.gd");
    let edited_gd = outdir.join("edited.gd");

    assert!(unintended_tsv.exists(), "unintended.tsv must be generated");
    assert!(summary_txt.exists(), "summary.txt must be generated");
    assert!(starter_gd.exists(), "starter.gd must be generated");
    assert!(edited_gd.exists(), "edited.gd must be generated");

    let tsv_content = std::fs::read_to_string(&unintended_tsv).unwrap();
    let lines: Vec<&str> = tsv_content.lines().collect();
    assert!(!lines.is_empty(), "unintended.tsv must have header");
    assert_eq!(
        lines[0],
        "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\tside2_seq_id\tside2_position\thypothesis"
    );

    // Historical SNPs (from historical_snps.txt) should NOT be in unintended.tsv
    let hist_content = std::fs::read_to_string(fixture_dir.join("historical_snps.txt")).unwrap();
    for line in hist_content.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let pos = line.split('\t').next().unwrap();
        assert!(
            !tsv_content.contains(&format!("\t{pos}\t")),
            "historical SNP at pos {pos} must not be in unintended.tsv"
        );
    }

    // Intended edit (pos 2500) must be masked out
    assert!(
        !tsv_content.contains("\t2500\t"),
        "intended edit at pos 2500 must be masked out"
    );

    // Extra edit (pos 8000) must be in unintended.tsv as scattered_snv
    assert!(
        tsv_content.contains("\t8000\t8000\tSNP\t"),
        "extra edit at pos 8000 must be present in unintended.tsv"
    );
    assert!(
        tsv_content.contains("scattered_snv"),
        "extra edit must be classified as scattered_snv"
    );

    // Check summary.txt
    let summary_content = std::fs::read_to_string(&summary_txt).unwrap();
    assert!(summary_content.contains("intended_provided\tyes"));
    assert!(summary_content.contains("intended_declared\t1"));
    assert!(summary_content.contains("intended_observed\t1"));
    assert!(summary_content.contains("intended_status\tall_observed"));
    assert!(summary_content.contains("intended_missing\t0"));
    assert!(summary_content.contains("scattered_snv\t1"));

    let _ = std::fs::remove_dir_all(outdir);
}
