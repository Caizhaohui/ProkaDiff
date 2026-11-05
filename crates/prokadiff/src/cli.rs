//! CLI parsing and validation (no Bowtie2).

use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand};
use prokadiff_report::SchemaVersion;

pub const STARTER_MANDATORY: &str =
    "starter-strain WGS is mandatory; provide both --starter and --edited. \
     NCBI/reference-only comparison is not supported.";

pub const EDITOR_ROADMAP: &str =
    "CAST / bRNA-IS110 editors are on the post-v0.1 roadmap; first-period --editor is \
     cas9 | cas12a | dsb only (see docs/schema.md).";

#[derive(Parser, Debug)]
#[command(
    name = "prokadiff",
    version,
    about = "Post-edit genome audit and edit verification for engineered prokaryotes — starter vs edited WGS",
    subcommand_negates_reqs = true
)]
pub struct Cli {
    /// Starter-strain FASTQ. Repeatable: two consecutive files = one PE pair (R1 then R2).
    #[arg(long, action = ArgAction::Append, num_args = 1..)]
    pub starter: Vec<PathBuf>,
    /// Edited-strain FASTQ. Same pairing rule as --starter.
    #[arg(long, action = ArgAction::Append, num_args = 1..)]
    pub edited: Vec<PathBuf>,
    /// Coordinate-skeleton reference (FASTA or GenBank). Repeatable.
    #[arg(long = "ref", action = ArgAction::Append, num_args = 1..)]
    pub refs: Vec<PathBuf>,
    /// Optional path to multi-clone batch manifest TSV (enables multi-clone engineering batch audit mode).
    #[arg(long = "batch-manifest")]
    pub batch_manifest: Option<PathBuf>,
    /// Optional intended-edit table. Omitted → every post-subtract call is unintended.
    #[arg(long)]
    pub intended: Option<PathBuf>,
    /// First period: cas9 | cas12a | dsb. cast / is110 error.
    #[arg(long)]
    pub editor: Option<String>,
    /// Target guide sequence (DNA alphabet, e.g. 20 nt; required for cas9/cas12a).
    #[arg(long)]
    pub spacer: Option<String>,
    /// PAM motif (IUPAC DNA, e.g. NGG for cas9, TTTV for cas12a; defaults provided).
    #[arg(long)]
    pub pam: Option<String>,
    /// Number of parallel threads for alignment and pileup.
    #[arg(long, default_value_t = 8)]
    pub threads: usize,
    /// Explicitly opt-in to experimental mechanistic hypothesis annotations (e.g. SOS-induced distal mutations).
    /// Default: false (production path remains strictly observational).
    #[arg(
        long = "experimental-hypothesis-annotation",
        visible_alias = "hypothesis",
        default_value_t = false
    )]
    pub experimental_hypothesis_annotation: bool,
    /// Deprecated legacy flag: omit hypothesis (hypothesis is now disabled by default).
    #[arg(long = "no-hypothesis", hide = true, default_value_t = false)]
    pub no_hypothesis: bool,
    #[arg(long, default_value = "v2")]
    pub schema_version: SchemaVersion,
    /// Directory where diff outputs (unintended.tsv, summary.txt, .gd) will be written.
    #[arg(long)]
    pub outdir: Option<PathBuf>,
    /// Keep intermediate BAM / Bowtie2 index under outdir.
    #[arg(long, default_value_t = false)]
    pub keep_bam: bool,
    /// Suppress informational progress messages (-q / --quiet).
    #[arg(short, long, default_value_t = false, global = true)]
    pub quiet: bool,
    /// Distance window in bp to associate mutations with off-target candidate sites (default: 50).
    #[arg(long = "offtarget-association-window", default_value_t = 50)]
    pub offtarget_association_window: u64,
    /// Maximum DNA bulge size (insertion in genomic DNA relative to guide) for off-target search (default: 0).
    #[arg(long = "max-dna-bulge", default_value_t = 0)]
    pub max_dna_bulge: u32,
    /// Maximum RNA bulge size (deletion in genomic DNA relative to guide) for off-target search (default: 0).
    #[arg(long = "max-rna-bulge", default_value_t = 0)]
    pub max_rna_bulge: u32,
    /// Minimum total supporting reads required to emit a JC record (default: 3).
    #[arg(long = "jc-min-support-reads", default_value_t = 3)]
    pub jc_min_support_reads: usize,
    /// Minimum allele frequency relative to local coverage for JC records (default: 0.05).
    #[arg(long = "jc-min-frequency", default_value_t = 0.05)]
    pub jc_min_frequency: f64,
    /// Enable subclonal (low-frequency) small variant audit mode (SNP, SUB, INS/DEL <= 2bp).
    #[arg(long = "subclonal", default_value_t = false)]
    pub subclonal: bool,
    /// Minimum alternate allele frequency threshold for subclonal audit mode (default: 0.05).
    #[arg(long = "min-subclonal-af", default_value_t = 0.05)]
    pub min_subclonal_af: f64,
    /// Minimum supporting ALT reads for subclonal audit mode (default: 3).
    #[arg(long = "min-subclonal-alt-reads", default_value_t = 3)]
    pub min_subclonal_alt_reads: u64,
    /// Maximum strand bias p-value for subclonal warning flag (default: 0.001).
    #[arg(long = "max-strand-bias-pvalue", default_value_t = 0.001)]
    pub max_strand_bias_pvalue: f64,
    /// Increase diagnostic logging verbosity (-v for DEBUG, -vv for TRACE).
    #[arg(short, long, action = ArgAction::Count, global = true)]
    pub verbose: u8,
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Single-sample RA/MC/JC engine (Bowtie2 + noodles). For oracle / parity jobs.
    Evidence(EvidenceArgs),
}

#[derive(Parser, Debug)]
pub struct EvidenceArgs {
    #[arg(long = "ref", action = ArgAction::Append, num_args = 1.., required = true)]
    pub refs: Vec<PathBuf>,
    /// FASTQ: one file = SE; two consecutive files = one PE pair.
    #[arg(long = "fastq", visible_alias = "reads", action = ArgAction::Append, num_args = 1.., required = true)]
    pub fastq: Vec<PathBuf>,
    #[arg(long, default_value_t = 8)]
    pub threads: usize,
    #[arg(long, required = true)]
    pub outdir: PathBuf,
    #[arg(long, default_value_t = false)]
    pub keep_bam: bool,
    /// Minimum total supporting reads required to emit a JC record (default: 3).
    #[arg(long = "jc-min-support-reads", default_value_t = 3)]
    pub jc_min_support_reads: usize,
    /// Minimum allele frequency relative to local coverage for JC records (default: 0.05).
    #[arg(long = "jc-min-frequency", default_value_t = 0.05)]
    pub jc_min_frequency: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Editor {
    Cas9,
    Cas12a,
    Dsb,
}

#[derive(Debug, PartialEq, Eq)]
pub enum CliError {
    StarterMandatory,
    EditorRoadmap { value: String },
    BadEditor { value: String },
    MissingRef,
    MissingOutdir,
    EmptyFastq,
    MissingSpacer,
    InvalidSpacer { reason: String },
    InvalidPam { reason: String },
    OddFastqCount { flag: &'static str, count: usize },
    FileNotFound { flag: &'static str, path: PathBuf },
    InvalidParameter { flag: &'static str, reason: String },
}

impl std::fmt::Display for CliError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::StarterMandatory => f.write_str(STARTER_MANDATORY),
            Self::EditorRoadmap { value } => {
                write!(f, "unsupported --editor '{value}': {EDITOR_ROADMAP}")
            }
            Self::BadEditor { value } => {
                write!(
                    f,
                    "unknown --editor '{value}'; supported first-period editors are cas9 | cas12a | dsb"
                )
            }
            Self::MissingRef => {
                f.write_str("--ref: at least one reference file (FASTA or GenBank) is required")
            }
            Self::MissingOutdir => f.write_str("--outdir: output directory is required"),
            Self::EmptyFastq => f.write_str("--fastq: at least one FASTQ file is required"),
            Self::MissingSpacer => f.write_str(
                "--spacer is required for --editor cas9 and cas12a (omit for --editor dsb)",
            ),
            Self::InvalidSpacer { reason } => write!(f, "invalid --spacer: {reason}"),
            Self::InvalidPam { reason } => write!(f, "invalid --pam: {reason}"),
            Self::OddFastqCount { flag, count } => write!(
                f,
                "{flag}: expected 1 file (SE) or an even number of files (PE pairs, R1 then R2), but got {count} file(s)"
            ),
            Self::FileNotFound { flag, path } => write!(
                f,
                "{flag}: file '{}' not found or not accessible",
                path.display()
            ),
            Self::InvalidParameter { flag, reason } => {
                write!(f, "invalid {flag}: {reason}")
            }
        }
    }
}

impl std::error::Error for CliError {}

pub fn parse_editor(raw: &str) -> Result<Editor, CliError> {
    match raw.to_ascii_lowercase().as_str() {
        "cas9" => Ok(Editor::Cas9),
        "cas12a" => Ok(Editor::Cas12a),
        "dsb" => Ok(Editor::Dsb),
        "cast" | "is110" => Err(CliError::EditorRoadmap {
            value: raw.to_string(),
        }),
        other => Err(CliError::BadEditor {
            value: other.to_string(),
        }),
    }
}

pub fn validate_spacer(raw: &str) -> Result<String, CliError> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(CliError::InvalidSpacer {
            reason: "spacer sequence cannot be empty".into(),
        });
    }
    if trimmed.contains('U') || trimmed.contains('u') {
        return Err(CliError::InvalidSpacer {
            reason: "spacer must use DNA alphabet (use 'T' instead of 'U')".into(),
        });
    }
    let upper = trimmed.to_ascii_uppercase();
    if let Some(ch) = upper
        .chars()
        .find(|c| !matches!(c, 'A' | 'C' | 'G' | 'T' | 'N'))
    {
        return Err(CliError::InvalidSpacer {
            reason: format!("spacer contains invalid character '{ch}'; expected A, C, G, T"),
        });
    }
    Ok(upper)
}

pub fn validate_pam(raw: &str) -> Result<String, CliError> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(CliError::InvalidPam {
            reason: "PAM sequence cannot be empty".into(),
        });
    }
    if trimmed.contains('U') || trimmed.contains('u') {
        return Err(CliError::InvalidPam {
            reason: "PAM must use DNA alphabet (use 'T' instead of 'U')".into(),
        });
    }
    let upper = trimmed.to_ascii_uppercase();
    if let Some(ch) = upper.chars().find(|c| {
        !matches!(
            c,
            'A' | 'C' | 'G' | 'T' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | 'N'
        )
    }) {
        return Err(CliError::InvalidPam {
            reason: format!("PAM contains invalid character '{ch}'; expected IUPAC DNA code"),
        });
    }
    Ok(upper)
}

#[derive(Debug)]
pub struct ProductJob {
    pub starter: Vec<PathBuf>,
    pub edited: Vec<PathBuf>,
    #[allow(dead_code)]
    pub batch_manifest: Option<PathBuf>,
    pub refs: Vec<PathBuf>,
    pub editor: Editor,
    pub spacer: Option<String>,
    pub pam: Option<String>,
    pub threads: usize,
    pub outdir: PathBuf,
    pub keep_bam: bool,
    pub intended: Option<PathBuf>,
    pub hypothesis: bool,
    pub schema_version: SchemaVersion,
    pub offtarget_association_window: u64,
    pub max_dna_bulge: u32,
    pub max_rna_bulge: u32,
    pub jc_min_support_reads: usize,
    pub jc_min_frequency: f64,
    #[allow(dead_code)]
    pub subclonal: bool,
    #[allow(dead_code)]
    pub min_subclonal_af: f64,
    #[allow(dead_code)]
    pub min_subclonal_alt_reads: u64,
    #[allow(dead_code)]
    pub max_strand_bias_pvalue: f64,
}

pub fn validate_product(cli: &Cli) -> Result<ProductJob, CliError> {
    if cli.starter.is_empty() {
        return Err(CliError::StarterMandatory);
    }
    if cli.batch_manifest.is_none() && cli.edited.is_empty() {
        return Err(CliError::StarterMandatory);
    }
    if cli.refs.is_empty() {
        return Err(CliError::MissingRef);
    }
    let Some(outdir) = cli.outdir.clone() else {
        return Err(CliError::MissingOutdir);
    };
    let Some(raw) = cli.editor.as_deref() else {
        return Err(CliError::BadEditor {
            value: "(missing)".into(),
        });
    };
    let editor = parse_editor(raw)?;
    check_fastq_pairing("--starter", &cli.starter)?;
    if !cli.edited.is_empty() {
        check_fastq_pairing("--edited", &cli.edited)?;
    }

    let is_batch = cli.batch_manifest.is_some();
    let spacer = match (editor, cli.spacer.as_deref()) {
        (Editor::Cas9 | Editor::Cas12a, None) if !is_batch => return Err(CliError::MissingSpacer),
        (Editor::Cas9 | Editor::Cas12a, Some(s)) if s.trim().is_empty() && !is_batch => {
            return Err(CliError::MissingSpacer)
        }
        (_, Some(s)) if !s.trim().is_empty() => Some(validate_spacer(s)?),
        _ => None,
    };

    let pam = match cli.pam.as_deref() {
        Some(p) if !p.trim().is_empty() => Some(validate_pam(p)?),
        _ => None,
    };

    for r in &cli.refs {
        if !r.exists() {
            return Err(CliError::FileNotFound {
                flag: "--ref",
                path: r.clone(),
            });
        }
    }
    for f in &cli.starter {
        if !f.exists() {
            return Err(CliError::FileNotFound {
                flag: "--starter",
                path: f.clone(),
            });
        }
    }
    for f in &cli.edited {
        if !f.exists() {
            return Err(CliError::FileNotFound {
                flag: "--edited",
                path: f.clone(),
            });
        }
    }
    if let Some(ref bm) = cli.batch_manifest {
        if !bm.exists() {
            return Err(CliError::FileNotFound {
                flag: "--batch-manifest",
                path: bm.clone(),
            });
        }
    }
    if let Some(ref p) = cli.intended {
        if !p.exists() {
            return Err(CliError::FileNotFound {
                flag: "--intended",
                path: p.clone(),
            });
        }
    }

    if !cli.min_subclonal_af.is_finite() || cli.min_subclonal_af < 0.01 || cli.min_subclonal_af > 0.50 {
        return Err(CliError::InvalidParameter {
            flag: "--min-subclonal-af",
            reason: format!(
                "value {} out of valid range [0.01, 0.50]",
                cli.min_subclonal_af
            ),
        });
    }
    if cli.min_subclonal_alt_reads < 1 {
        return Err(CliError::InvalidParameter {
            flag: "--min-subclonal-alt-reads",
            reason: "minimum ALT reads must be >= 1".into(),
        });
    }
    if !cli.max_strand_bias_pvalue.is_finite() || cli.max_strand_bias_pvalue < 0.0 || cli.max_strand_bias_pvalue > 1.0 {
        return Err(CliError::InvalidParameter {
            flag: "--max-strand-bias-pvalue",
            reason: format!(
                "value {} out of valid range [0.0, 1.0]",
                cli.max_strand_bias_pvalue
            ),
        });
    }

    Ok(ProductJob {
        starter: cli.starter.clone(),
        edited: cli.edited.clone(),
        batch_manifest: cli.batch_manifest.clone(),
        refs: cli.refs.clone(),
        editor,
        spacer,
        pam,
        threads: cli.threads.max(1),
        outdir,
        keep_bam: cli.keep_bam,
        intended: cli.intended.clone(),
        hypothesis: cli.experimental_hypothesis_annotation && !cli.no_hypothesis,
        schema_version: cli.schema_version,
        offtarget_association_window: cli.offtarget_association_window,
        max_dna_bulge: cli.max_dna_bulge,
        max_rna_bulge: cli.max_rna_bulge,
        jc_min_support_reads: cli.jc_min_support_reads,
        jc_min_frequency: cli.jc_min_frequency,
        subclonal: cli.subclonal,
        min_subclonal_af: cli.min_subclonal_af,
        min_subclonal_alt_reads: cli.min_subclonal_alt_reads,
        max_strand_bias_pvalue: cli.max_strand_bias_pvalue,
    })
}

fn check_fastq_pairing(flag: &'static str, files: &[PathBuf]) -> Result<(), CliError> {
    if files.is_empty() {
        return Ok(());
    }
    if files.len() != 1 && files.len() % 2 != 0 {
        return Err(CliError::OddFastqCount {
            flag,
            count: files.len(),
        });
    }
    Ok(())
}

pub fn validate_evidence(args: &EvidenceArgs) -> Result<(), CliError> {
    if args.refs.is_empty() {
        return Err(CliError::MissingRef);
    }
    if args.fastq.is_empty() {
        return Err(CliError::EmptyFastq);
    }
    check_fastq_pairing("--fastq", &args.fastq)?;
    for r in &args.refs {
        if !r.exists() {
            return Err(CliError::FileNotFound {
                flag: "--ref",
                path: r.clone(),
            });
        }
    }
    for f in &args.fastq {
        if !f.exists() {
            return Err(CliError::FileNotFound {
                flag: "--fastq",
                path: f.clone(),
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static TEST_COUNTER: AtomicUsize = AtomicUsize::new(0);

    fn test_dir(prefix: &str) -> PathBuf {
        let id = TEST_COUNTER.fetch_add(1, Ordering::SeqCst);
        let dir = std::env::temp_dir().join(format!(
            "prokadiff_cli_test_{prefix}_{id}_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn touch(dir: &std::path::Path, name: &str) -> PathBuf {
        let p = dir.join(name);
        std::fs::write(&p, b"").unwrap();
        p
    }

    #[test]
    fn missing_starter_is_mandatory_error() {
        let cli = Cli::try_parse_from([
            "prokdiff", "--edited", "e.fq", "--ref", "r.fa", "--editor", "cas9", "--outdir", "out",
        ])
        .unwrap();
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(err, CliError::StarterMandatory));
        assert!(err.to_string().contains("starter-strain WGS is mandatory"));
    }

    #[test]
    fn schema_version_defaults_to_v2_and_accepts_v1() {
        let default_cli = Cli::try_parse_from(["prokadiff"]).expect("schema version has a default");
        assert_eq!(default_cli.schema_version, SchemaVersion::V2);
        let legacy_cli = Cli::try_parse_from(["prokadiff", "--schema-version", "v1"])
            .expect("v1 schema version parses");
        assert_eq!(legacy_cli.schema_version, SchemaVersion::V1);
    }

    #[test]
    fn missing_edited_is_also_mandatory_error() {
        let cli = Cli::try_parse_from([
            "prokdiff",
            "--starter",
            "s.fq",
            "--ref",
            "r.fa",
            "--editor",
            "cas9",
            "--outdir",
            "out",
        ])
        .unwrap();
        assert!(matches!(
            validate_product(&cli).unwrap_err(),
            CliError::StarterMandatory
        ));
    }

    #[test]
    fn cast_editor_points_to_roadmap() {
        let err = parse_editor("cast").unwrap_err();
        assert!(err.to_string().contains("roadmap"));
        let err = parse_editor("is110").unwrap_err();
        assert!(err.to_string().contains("roadmap"));
    }

    #[test]
    fn evidence_subcommand_does_not_require_starter() {
        let d = test_dir("evidence_ok");
        let r = touch(&d, "r.fa");
        let a = touch(&d, "a.fq");
        let b = touch(&d, "b.fq");
        let cli = Cli::try_parse_from([
            "prokdiff",
            "evidence",
            "--ref",
            r.to_str().unwrap(),
            "--fastq",
            a.to_str().unwrap(),
            "--fastq",
            b.to_str().unwrap(),
            "--outdir",
            "out",
        ])
        .unwrap();
        match cli.command {
            Some(Commands::Evidence(args)) => {
                validate_evidence(&args).unwrap();
                assert_eq!(args.fastq.len(), 2);
            }
            other => panic!("expected evidence, got {other:?}"),
        }
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn generate_script_pins_wgsim_no_extra_mutations() {
        let p =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../testdata/generate.sh");
        let s = std::fs::read_to_string(&p).unwrap_or_else(|_| panic!("missing {}", p.display()));
        assert!(s.contains("wgsim"), "generate.sh must call wgsim");
        assert!(s.contains("-r 0"));
        assert!(s.contains("-R 0"));
        assert!(s.contains("-X 0"));
        assert!(
            s.contains("login"),
            "generate.sh must refuse the login node"
        );
    }

    #[test]
    fn cas9_requires_spacer() {
        let cli = Cli::try_parse_from([
            "prokdiff",
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
        .unwrap();
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(err, CliError::MissingSpacer));
    }

    #[test]
    fn dsb_does_not_require_spacer() {
        let d = test_dir("dsb_ok");
        let s = touch(&d, "s.fq");
        let e = touch(&d, "e.fq");
        let r = touch(&d, "r.fa");
        let cli = Cli::try_parse_from([
            "prokdiff",
            "--starter",
            s.to_str().unwrap(),
            "--edited",
            e.to_str().unwrap(),
            "--ref",
            r.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            "out",
            "--no-hypothesis",
        ])
        .unwrap();
        let job = validate_product(&cli).unwrap();
        assert_eq!(job.editor, Editor::Dsb);
        assert!(job.spacer.is_none());
        assert!(!job.hypothesis);
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn odd_fastq_count_is_rejected() {
        let cli = Cli::try_parse_from([
            "prokdiff",
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
        .unwrap();
        assert!(matches!(
            validate_product(&cli).unwrap_err(),
            CliError::OddFastqCount {
                flag: "--starter",
                count: 3
            }
        ));
    }

    #[test]
    fn intended_path_is_optional_and_stored() {
        let d = test_dir("intended_ok");
        let s = touch(&d, "s.fq");
        let e = touch(&d, "e.fq");
        let r = touch(&d, "r.fa");
        let i = touch(&d, "intended.tsv");
        let cli = Cli::try_parse_from([
            "prokdiff",
            "--starter",
            s.to_str().unwrap(),
            "--edited",
            e.to_str().unwrap(),
            "--ref",
            r.to_str().unwrap(),
            "--editor",
            "cas9",
            "--spacer",
            "acgtacgtacgtacgtacgt",
            "--intended",
            i.to_str().unwrap(),
            "--outdir",
            "out",
        ])
        .unwrap();
        let job = validate_product(&cli).unwrap();
        assert_eq!(job.spacer.as_deref(), Some("ACGTACGTACGTACGTACGT"));
        assert_eq!(job.intended.as_ref(), Some(&i));
        assert!(!job.hypothesis);
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn hypothesis_opt_in_flag_enables_hypothesis() {
        let d = test_dir("hyp_opt_in");
        let s = touch(&d, "s.fq");
        let e = touch(&d, "e.fq");
        let r = touch(&d, "r.fa");
        let cli = Cli::try_parse_from([
            "prokdiff",
            "--starter",
            s.to_str().unwrap(),
            "--edited",
            e.to_str().unwrap(),
            "--ref",
            r.to_str().unwrap(),
            "--editor",
            "cas9",
            "--spacer",
            "acgtacgtacgtacgtacgt",
            "--outdir",
            "out",
            "--experimental-hypothesis-annotation",
        ])
        .unwrap();
        let job = validate_product(&cli).unwrap();
        assert!(job.hypothesis);
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn file_not_found_on_starter() {
        let d = test_dir("fnf_starter");
        let e = touch(&d, "e.fq");
        let r = touch(&d, "r.fa");
        let missing = d.join("missing_starter.fq");
        let cli = Cli::try_parse_from([
            "prokdiff",
            "--starter",
            missing.to_str().unwrap(),
            "--edited",
            e.to_str().unwrap(),
            "--ref",
            r.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            "out",
        ])
        .unwrap();
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(
            err,
            CliError::FileNotFound {
                flag: "--starter",
                ..
            }
        ));
        assert!(err.to_string().contains("not found"));
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn file_not_found_on_ref() {
        let d = test_dir("fnf_ref");
        let s = touch(&d, "s.fq");
        let e = touch(&d, "e.fq");
        let missing = d.join("missing_ref.fa");
        let cli = Cli::try_parse_from([
            "prokdiff",
            "--starter",
            s.to_str().unwrap(),
            "--edited",
            e.to_str().unwrap(),
            "--ref",
            missing.to_str().unwrap(),
            "--editor",
            "dsb",
            "--outdir",
            "out",
        ])
        .unwrap();
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(err, CliError::FileNotFound { flag: "--ref", .. }));
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn invalid_spacer_rejects_uracil_and_bad_chars() {
        assert!(validate_spacer("").is_err());
        let u_err = validate_spacer("AUGCAUGCAUGCAUGCAUGC").unwrap_err();
        assert!(u_err.to_string().contains("use 'T' instead of 'U'"));
        let bad_err = validate_spacer("ACGT123").unwrap_err();
        assert!(bad_err.to_string().contains("invalid character"));
        assert_eq!(
            validate_spacer("acgtacgtacgtacgtacgt").unwrap(),
            "ACGTACGTACGTACGTACGT"
        );
    }

    #[test]
    fn invalid_pam_rejects_uracil_and_bad_chars() {
        assert!(validate_pam("").is_err());
        let u_err = validate_pam("NGGU").unwrap_err();
        assert!(u_err.to_string().contains("use 'T' instead of 'U'"));
        let bad_err = validate_pam("NGG@").unwrap_err();
        assert!(bad_err.to_string().contains("invalid character"));
        assert_eq!(validate_pam("ngg").unwrap(), "NGG");
        assert_eq!(validate_pam("tttv").unwrap(), "TTTV");
    }

    #[test]
    fn quiet_and_verbose_flags_parse_correctly() {
        let cli = Cli::try_parse_from([
            "prokadiff",
            "--starter",
            "s.fq",
            "--edited",
            "e.fq",
            "--ref",
            "r.fa",
            "-q",
            "-vv",
        ])
        .unwrap();
        assert!(cli.quiet);
        assert_eq!(cli.verbose, 2);
    }

    #[test]
    fn cli_parses_batch_manifest_flag() {
        let cli = Cli::try_parse_from([
            "prokadiff",
            "--starter",
            "s.fq",
            "--batch-manifest",
            "batch_manifest.tsv",
            "--ref",
            "r.fa",
            "--editor",
            "cas9",
            "--spacer",
            "GAGTTCATCTACGCCGTGAA",
            "--outdir",
            "out/",
        ])
        .unwrap();
        assert_eq!(
            cli.batch_manifest,
            Some(PathBuf::from("batch_manifest.tsv"))
        );
        assert!(cli.edited.is_empty());
    }

    #[test]
    fn batch_manifest_mode_still_mandates_starter() {
        let cli = Cli::try_parse_from([
            "prokadiff",
            "--batch-manifest",
            "manifest.tsv",
            "--ref",
            "r.fa",
            "--editor",
            "cas9",
            "--spacer",
            "GAGTTCATCTACGCCGTGAA",
            "--outdir",
            "out/",
        ])
        .unwrap();
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(err, CliError::StarterMandatory));
    }

    #[test]
    fn subclonal_parameter_bounds_validation() {
        let d = test_dir("subclonal_bounds");
        let s = touch(&d, "s.fq");
        let e = touch(&d, "e.fq");
        let r = touch(&d, "r.fa");

        let make_cli = |af: &str, alt_reads: &str, sb_p: &str| {
            Cli::try_parse_from([
                "prokadiff",
                "--starter",
                s.to_str().unwrap(),
                "--edited",
                e.to_str().unwrap(),
                "--ref",
                r.to_str().unwrap(),
                "--editor",
                "dsb",
                "--outdir",
                "out/",
                "--subclonal",
                "--min-subclonal-af",
                af,
                "--min-subclonal-alt-reads",
                alt_reads,
                "--max-strand-bias-pvalue",
                sb_p,
            ])
            .unwrap()
        };

        // AF too low (< 0.01)
        let cli = make_cli("0.005", "3", "0.001");
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(
            err,
            CliError::InvalidParameter {
                flag: "--min-subclonal-af",
                ..
            }
        ));

        // AF too high (> 0.50)
        let cli = make_cli("0.55", "3", "0.001");
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(
            err,
            CliError::InvalidParameter {
                flag: "--min-subclonal-af",
                ..
            }
        ));

        // ALT reads < 1
        let cli = make_cli("0.05", "0", "0.001");
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(
            err,
            CliError::InvalidParameter {
                flag: "--min-subclonal-alt-reads",
                ..
            }
        ));

        // Strand bias p-value > 1.0
        let cli = make_cli("0.05", "3", "1.5");
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(
            err,
            CliError::InvalidParameter {
                flag: "--max-strand-bias-pvalue",
                ..
            }
        ));

        // AF is NaN
        let cli = make_cli("NaN", "3", "0.001");
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(
            err,
            CliError::InvalidParameter {
                flag: "--min-subclonal-af",
                ..
            }
        ));

        // Strand bias p-value is NaN
        let cli = make_cli("0.05", "3", "NaN");
        let err = validate_product(&cli).unwrap_err();
        assert!(matches!(
            err,
            CliError::InvalidParameter {
                flag: "--max-strand-bias-pvalue",
                ..
            }
        ));

        let _ = std::fs::remove_dir_all(d);
    }
}
