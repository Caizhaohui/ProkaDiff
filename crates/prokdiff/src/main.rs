//! ProkDiff CLI. Engine path B: Bowtie2 + Rust RA/MC/JC. No runtime breseq.
//! Classification (near_homolog / scattered_snv / structural) is stage 3.

mod cli;

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use prokdiff_evidence::align::FastqInput;
use prokdiff_evidence::engine::{run_sample, EngineOptions};
use prokdiff_evidence::fasta::write_combined_fasta;

use cli::{validate_evidence, validate_product, Cli, CliError, Commands, EvidenceArgs, ProductJob};

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

enum RunError {
    Cli(CliError),
    Evidence(prokdiff_evidence::EvidenceError),
    Io(std::io::Error),
}

impl From<CliError> for RunError {
    fn from(e: CliError) -> Self {
        Self::Cli(e)
    }
}

impl From<prokdiff_evidence::EvidenceError> for RunError {
    fn from(e: prokdiff_evidence::EvidenceError) -> Self {
        Self::Evidence(e)
    }
}

impl From<std::io::Error> for RunError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

impl std::fmt::Display for RunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Cli(e) => write!(f, "{e}"),
            Self::Evidence(e) => write!(f, "{e}"),
            Self::Io(e) => write!(f, "{e}"),
        }
    }
}

fn run(cli: Cli) -> Result<(), RunError> {
    match cli.command {
        Some(Commands::Evidence(args)) => run_evidence(args),
        None => {
            let job = validate_product(&cli)?;
            run_product(job)
        }
    }
}

fn run_evidence(args: EvidenceArgs) -> Result<(), RunError> {
    validate_evidence(&args)?;
    let work = args.outdir.join("work");
    std::fs::create_dir_all(&work)?;
    let ref_fa = materialize_ref(&args.refs, &work)?;
    let reads = FastqInput {
        files: args.fastq.clone(),
    };
    let opts = EngineOptions {
        threads: args.threads.max(1),
        keep_bam: args.keep_bam,
        ..EngineOptions::default()
    };
    let gd = run_sample(&ref_fa, &reads, &args.outdir, &opts)?;
    eprintln!("wrote {}", gd.display());
    Ok(())
}

fn run_product(job: ProductJob) -> Result<(), RunError> {
    let _ = job.editor; // stage 3 classify; validated only
    if job.intended.is_some() {
        eprintln!("note: --intended is accepted but not applied until stage 3 classify");
    }
    std::fs::create_dir_all(&job.outdir)?;
    let work = job.outdir.join("work");
    std::fs::create_dir_all(&work)?;
    let ref_fa = materialize_ref(&job.refs, &work)?;
    let opts = EngineOptions {
        threads: job.threads,
        keep_bam: job.keep_bam,
        ..EngineOptions::default()
    };

    let starter_dir = job.outdir.join("starter");
    let edited_dir = job.outdir.join("edited");
    let starter_gd = run_sample(
        &ref_fa,
        &FastqInput {
            files: job.starter.clone(),
        },
        &starter_dir,
        &opts,
    )?;
    let edited_gd = run_sample(
        &ref_fa,
        &FastqInput {
            files: job.edited.clone(),
        },
        &edited_dir,
        &opts,
    )?;

    let starter_out = job.outdir.join("starter.gd");
    let edited_out = job.outdir.join("edited.gd");
    std::fs::copy(&starter_gd, &starter_out)?;
    std::fs::copy(&edited_gd, &edited_out)?;
    eprintln!("wrote {}", starter_out.display());
    eprintln!("wrote {}", edited_out.display());
    eprintln!(
        "note: mutation classification (structural / near_homolog / scattered_snv) is not in this engine MVP"
    );
    Ok(())
}

fn materialize_ref(refs: &[PathBuf], work: &Path) -> Result<PathBuf, RunError> {
    if refs.len() == 1 {
        let p = &refs[0];
        let name = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
        if name.ends_with(".fa")
            || name.ends_with(".fna")
            || name.ends_with(".fasta")
            || name.ends_with(".fa.gz")
        {
            // Uncompressed FASTA can be used directly by bowtie2-build.
            if !name.ends_with(".gz") {
                return Ok(p.clone());
            }
        }
    }
    let dest = work.join("reference.fa");
    write_combined_fasta(refs, &dest)?;
    Ok(dest)
}
