//! ProkDiff CLI. Engine path B: Bowtie2 + Rust RA/MC/JC. No runtime breseq.
//! Product path: both strains → subtract → intended mask → classes (3)→(1)→(2).

mod cli;

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use prokadiff_classify::{
    classify, parse_intended_path, ClassifyOptions, EditorKind, RefContig, DEFAULT_MAX_MISMATCHES,
    DEFAULT_NEAR_DISTANCE,
};
use prokadiff_evidence::align::FastqInput;
use prokadiff_evidence::engine::{run_sample, EngineOptions};
use prokadiff_evidence::fasta::{read_reference, read_references, write_combined_fasta};
use prokadiff_gd::GenomeDiff;
use prokadiff_offtarget::{
    link_mutations_to_sites, scan_genome_bulge, write_mutation_offtarget_links_tsv,
    write_offtarget_sites_tsv, BulgeSearchOptions, NucleaseProfile, PamSide,
};
use prokadiff_report::{write_summary, write_unintended_tsv};

use cli::{
    validate_evidence, validate_product, Cli, CliError, Commands, Editor, EvidenceArgs, ProductJob,
};
use tracing::info;

fn init_logging(quiet: bool, verbose: u8) {
    let default_level = match (quiet, verbose) {
        (true, _) => "warn",
        (false, 0) => "info",
        (false, 1) => "debug",
        (false, _) => "trace",
    };
    let filter = tracing_subscriber::EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new(default_level));

    let _ = tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .without_time()
        .with_writer(std::io::stderr)
        .try_init();
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    init_logging(cli.quiet, cli.verbose);
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
    Evidence(prokadiff_evidence::EvidenceError),
    Gd(prokadiff_gd::GdError),
    Intended(prokadiff_classify::IntendedError),
    Io(std::io::Error),
}

impl From<CliError> for RunError {
    fn from(e: CliError) -> Self {
        Self::Cli(e)
    }
}

impl From<prokadiff_evidence::EvidenceError> for RunError {
    fn from(e: prokadiff_evidence::EvidenceError) -> Self {
        Self::Evidence(e)
    }
}

impl From<prokadiff_gd::GdError> for RunError {
    fn from(e: prokadiff_gd::GdError) -> Self {
        Self::Gd(e)
    }
}

impl From<prokadiff_classify::IntendedError> for RunError {
    fn from(e: prokadiff_classify::IntendedError) -> Self {
        Self::Intended(e)
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
            Self::Gd(e) => write!(f, "{e}"),
            Self::Intended(e) => write!(f, "{e}"),
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

fn load_repeats(refs: &[PathBuf]) -> Vec<prokadiff_evidence::fasta::RepeatRegion> {
    let mut repeats = Vec::new();
    for r in refs {
        if let Ok(reps) = prokadiff_evidence::fasta::parse_genbank_repeats(r) {
            repeats.extend(reps);
        }
    }
    repeats
}

fn run_evidence(args: EvidenceArgs) -> Result<(), RunError> {
    validate_evidence(&args)?;
    let work = args.outdir.join("work");
    std::fs::create_dir_all(&work)?;
    let ref_fa = materialize_ref(&args.refs, &work)?;
    let reads = FastqInput {
        files: args.fastq.clone(),
    };
    let repeats = load_repeats(&args.refs);
    let opts = EngineOptions {
        threads: args.threads.max(1),
        keep_bam: args.keep_bam,
        repeats,
        ..EngineOptions::default()
    };
    let gd = run_sample(&ref_fa, &reads, &args.outdir, &opts)?;
    info!("wrote {}", gd.display());
    Ok(())
}

fn run_product(job: ProductJob) -> Result<(), RunError> {
    std::fs::create_dir_all(&job.outdir)?;
    let work = job.outdir.join("work");
    std::fs::create_dir_all(&work)?;
    let ref_fa = materialize_ref(&job.refs, &work)?;
    let repeats = load_repeats(&job.refs);
    let opts = EngineOptions {
        threads: job.threads,
        keep_bam: job.keep_bam,
        repeats,
        ..EngineOptions::default()
    };

    let starter_dir = job.outdir.join("starter");
    let edited_dir = job.outdir.join("edited");
    let starter_gd_path = run_sample(
        &ref_fa,
        &FastqInput {
            files: job.starter.clone(),
        },
        &starter_dir,
        &opts,
    )?;
    let edited_gd_path = run_sample(
        &ref_fa,
        &FastqInput {
            files: job.edited.clone(),
        },
        &edited_dir,
        &opts,
    )?;

    let starter_out = job.outdir.join("starter.gd");
    let edited_out = job.outdir.join("edited.gd");
    std::fs::copy(&starter_gd_path, &starter_out)?;
    std::fs::copy(&edited_gd_path, &edited_out)?;

    let starter_gd = GenomeDiff::from_path(&starter_out)?;
    let edited_gd = GenomeDiff::from_path(&edited_out)?;
    let intended = match &job.intended {
        Some(p) => parse_intended_path(p)?,
        None => Vec::new(),
    };
    let fasta = read_reference(&ref_fa)?;
    let refs: Vec<RefContig> = fasta
        .into_iter()
        .map(|r| RefContig {
            name: r.name,
            seq: r.seq,
        })
        .collect();
    let editor_kind = match job.editor {
        Editor::Cas9 => EditorKind::Cas9,
        Editor::Cas12a => EditorKind::Cas12a,
        Editor::Dsb => EditorKind::Dsb,
    };
    let classified = classify(
        &edited_gd,
        &starter_gd,
        &intended,
        &refs,
        &ClassifyOptions {
            editor: editor_kind,
            spacer: job.spacer.clone(),
            pam: job.pam.clone(),
            near_distance: DEFAULT_NEAR_DISTANCE,
            max_mismatches: DEFAULT_MAX_MISMATCHES,
            hypothesis: job.hypothesis,
        },
    );

    let tsv = job.outdir.join("unintended.tsv");
    let summary = job.outdir.join("summary.txt");
    write_unintended_tsv(
        &tsv,
        &classified.unintended,
        editor_kind.as_str(),
        job.hypothesis,
        &refs,
    )?;
    write_summary(
        &summary,
        &classified,
        job.intended.is_some(),
        editor_kind.as_str(),
    )?;

    let offtarget_tsv = job.outdir.join("offtarget_sites.tsv");
    let links_tsv = job.outdir.join("mutation_offtarget_links.tsv");

    let (offtarget_sites, offtarget_links) = if let Some(spacer) = &job.spacer {
        let profile = match job.editor {
            Editor::Cas9 => {
                let pam = job.pam.as_deref().unwrap_or("NGG");
                NucleaseProfile::custom("SpCas9", spacer.len(), pam, PamSide::ThreePrime)
            }
            Editor::Cas12a => {
                let pam = job.pam.as_deref().unwrap_or("TTTV");
                NucleaseProfile::custom("Cas12a", spacer.len(), pam, PamSide::FivePrime)
            }
            Editor::Dsb => NucleaseProfile::custom("DSB", spacer.len(), "", PamSide::ThreePrime),
        };
        if job.editor != Editor::Dsb {
            let ref_tuples: Vec<(String, Vec<u8>)> = refs
                .iter()
                .map(|r| (r.name.clone(), r.seq.clone()))
                .collect();
            let bulge_opts = BulgeSearchOptions {
                max_mismatches: DEFAULT_MAX_MISMATCHES,
                max_dna_bulge: job.max_dna_bulge,
                max_rna_bulge: job.max_rna_bulge,
            };
            let sites = scan_genome_bulge(&ref_tuples, spacer, &profile, bulge_opts);
            let unintended_entries: Vec<prokadiff_gd::GdEntry> = classified
                .unintended
                .iter()
                .map(|u| u.entry.clone())
                .collect();
            let links = link_mutations_to_sites(
                &unintended_entries,
                &sites,
                job.offtarget_association_window,
            );
            (sites, links)
        } else {
            (Vec::new(), Vec::new())
        }
    } else {
        (Vec::new(), Vec::new())
    };

    write_offtarget_sites_tsv(&offtarget_sites, &offtarget_tsv)?;
    write_mutation_offtarget_links_tsv(&offtarget_links, &links_tsv)?;

    info!("wrote {}", starter_out.display());
    info!("wrote {}", edited_out.display());
    info!("wrote {}", tsv.display());
    info!("wrote {}", summary.display());
    info!("wrote {}", offtarget_tsv.display());
    info!("wrote {}", links_tsv.display());
    Ok(())
}

fn materialize_ref(refs: &[PathBuf], work: &Path) -> Result<PathBuf, RunError> {
    read_references(refs)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn materialize_ref_rejects_duplicate_id_in_direct_single_fasta() {
        // Given: a direct Bowtie2-compatible FASTA with a duplicate contig ID.
        let dir =
            std::env::temp_dir().join(format!("prokadiff-direct-reference-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let source = dir.join("duplicate.fa");
        let mut file = std::fs::File::create(&source).unwrap();
        writeln!(file, ">chr\nACGT\n>chr\nTGCA").unwrap();

        // When: direct reference materialization is requested.
        let result = materialize_ref(std::slice::from_ref(&source), &dir.join("work"));

        // Then: validation fails before Bowtie2 could receive the direct file.
        assert!(matches!(result, Err(RunError::Evidence(_))));
    }
}
