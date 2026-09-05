//! CLI parsing and validation (no Bowtie2).

use std::path::PathBuf;

use clap::{ArgAction, Parser, Subcommand};

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
    about = "Prokaryotic genome diff for gene-editing QC — starter vs edited WGS",
    subcommand_negates_reqs = true
)]
pub struct Cli {
    /// Starter-strain FASTQ. Repeatable: two consecutive files = one PE pair (R1 then R2).
    #[arg(long, action = ArgAction::Append)]
    pub starter: Vec<PathBuf>,
    /// Edited-strain FASTQ. Same pairing rule as --starter.
    #[arg(long, action = ArgAction::Append)]
    pub edited: Vec<PathBuf>,
    /// Coordinate-skeleton reference (FASTA or GenBank). Repeatable.
    #[arg(long = "ref", action = ArgAction::Append)]
    pub refs: Vec<PathBuf>,
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
    /// Omit the hypothesis column from unintended.tsv.
    #[arg(long, default_value_t = false)]
    pub no_hypothesis: bool,
    /// Directory where diff outputs (unintended.tsv, summary.txt, .gd) will be written.
    #[arg(long)]
    pub outdir: Option<PathBuf>,
    /// Keep intermediate BAM / Bowtie2 index under outdir.
    #[arg(long, default_value_t = false)]
    pub keep_bam: bool,
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
    #[arg(long = "ref", action = ArgAction::Append, required = true)]
    pub refs: Vec<PathBuf>,
    /// FASTQ: one file = SE; two consecutive files = one PE pair.
    #[arg(long = "fastq", visible_alias = "reads", action = ArgAction::Append, required = true)]
    pub fastq: Vec<PathBuf>,
    #[arg(long, default_value_t = 8)]
    pub threads: usize,
    #[arg(long, required = true)]
    pub outdir: PathBuf,
    #[arg(long, default_value_t = false)]
    pub keep_bam: bool,
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
    pub refs: Vec<PathBuf>,
    pub editor: Editor,
    pub spacer: Option<String>,
    pub pam: Option<String>,
    pub threads: usize,
    pub outdir: PathBuf,
    pub keep_bam: bool,
    pub intended: Option<PathBuf>,
    pub hypothesis: bool,
}

pub fn validate_product(cli: &Cli) -> Result<ProductJob, CliError> {
    if cli.starter.is_empty() || cli.edited.is_empty() {
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
    check_fastq_pairing("--edited", &cli.edited)?;

    let spacer = match (editor, cli.spacer.as_deref()) {
        (Editor::Cas9 | Editor::Cas12a, None) => return Err(CliError::MissingSpacer),
        (Editor::Cas9 | Editor::Cas12a, Some(s)) if s.trim().is_empty() => {
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
    if let Some(ref p) = cli.intended {
        if !p.exists() {
            return Err(CliError::FileNotFound {
                flag: "--intended",
                path: p.clone(),
            });
        }
    }

    Ok(ProductJob {
        starter: cli.starter.clone(),
        edited: cli.edited.clone(),
        refs: cli.refs.clone(),
        editor,
        spacer,
        pam,
        threads: cli.threads.max(1),
        outdir,
        keep_bam: cli.keep_bam,
        intended: cli.intended.clone(),
        hypothesis: !cli.no_hypothesis,
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
}
