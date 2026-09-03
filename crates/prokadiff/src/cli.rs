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
    #[arg(long)]
    pub spacer: Option<String>,
    #[arg(long)]
    pub pam: Option<String>,
    #[arg(long, default_value_t = 8)]
    pub threads: usize,
    /// Omit the hypothesis column from unintended.tsv.
    #[arg(long, default_value_t = false)]
    pub no_hypothesis: bool,
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

#[derive(Debug)]
pub enum CliError {
    StarterMandatory,
    EditorRoadmap { value: String },
    BadEditor { value: String },
    MissingRef,
    MissingOutdir,
    EmptyFastq,
    MissingSpacer,
    OddFastqCount { flag: &'static str },
}

impl std::fmt::Display for CliError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::StarterMandatory => f.write_str(STARTER_MANDATORY),
            Self::EditorRoadmap { value } => {
                write!(f, "unsupported --editor {value}: {EDITOR_ROADMAP}")
            }
            Self::BadEditor { value } => {
                write!(
                    f,
                    "unknown --editor {value}; first-period values are cas9 | cas12a | dsb"
                )
            }
            Self::MissingRef => f.write_str("at least one --ref is required"),
            Self::MissingOutdir => f.write_str("--outdir is required"),
            Self::EmptyFastq => f.write_str("at least one FASTQ is required"),
            Self::MissingSpacer => f.write_str(
                "--spacer is required for --editor cas9 and cas12a (omit for --editor dsb)",
            ),
            Self::OddFastqCount { flag } => write!(
                f,
                "{flag}: expected 1 file (SE) or an even number of files (PE pairs, R1 then R2)"
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
    let spacer = cli
        .spacer
        .as_ref()
        .map(|s| s.trim().to_ascii_uppercase())
        .filter(|s| !s.is_empty());
    match editor {
        Editor::Cas9 | Editor::Cas12a if spacer.is_none() => {
            return Err(CliError::MissingSpacer);
        }
        _ => {}
    }
    Ok(ProductJob {
        starter: cli.starter.clone(),
        edited: cli.edited.clone(),
        refs: cli.refs.clone(),
        editor,
        spacer,
        pam: cli.pam.clone(),
        threads: cli.threads.max(1),
        outdir,
        keep_bam: cli.keep_bam,
        intended: cli.intended.clone(),
        hypothesis: !cli.no_hypothesis,
    })
}

fn check_fastq_pairing(flag: &'static str, files: &[PathBuf]) -> Result<(), CliError> {
    if files.is_empty() {
        return Err(CliError::StarterMandatory);
    }
    if files.len() != 1 && files.len() % 2 != 0 {
        return Err(CliError::OddFastqCount { flag });
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
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

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
        let cli = Cli::try_parse_from([
            "prokdiff", "evidence", "--ref", "r.fa", "--fastq", "a.fq", "--fastq", "b.fq",
            "--outdir", "out",
        ])
        .unwrap();
        match cli.command {
            Some(Commands::Evidence(args)) => {
                validate_evidence(&args).unwrap();
                assert_eq!(args.fastq.len(), 2);
            }
            other => panic!("expected evidence, got {other:?}"),
        }
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
        let cli = Cli::try_parse_from([
            "prokdiff",
            "--starter",
            "s.fq",
            "--edited",
            "e.fq",
            "--ref",
            "r.fa",
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
            CliError::OddFastqCount { flag: "--starter" }
        ));
    }

    #[test]
    fn intended_path_is_optional_and_stored() {
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
            "--spacer",
            "acgtacgtacgtacgtacgt",
            "--intended",
            "intended.tsv",
            "--outdir",
            "out",
        ])
        .unwrap();
        let job = validate_product(&cli).unwrap();
        assert_eq!(job.spacer.as_deref(), Some("ACGTACGTACGTACGTACGT"));
        assert_eq!(
            job.intended.as_deref(),
            Some(std::path::Path::new("intended.tsv"))
        );
        assert!(job.hypothesis);
    }
}
