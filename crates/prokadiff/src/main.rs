//! ProkDiff CLI. Engine path B: Bowtie2 + Rust RA/MC/JC. No runtime breseq.
//! Product path: both strains → subtract → intended mask → classes (3)→(1)→(2).

mod cli;

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use prokadiff_classify::{
    classify, parse_intended_path, ClassifyOptions, EditorKind, RefContig, DEFAULT_MAX_MISMATCHES,
};
use prokadiff_evidence::align::FastqInput;
use prokadiff_evidence::engine::{run_sample, EngineOptions};
use prokadiff_evidence::fasta::{read_reference, read_references, write_combined_fasta};
use prokadiff_gd::GenomeDiff;
use prokadiff_report::{write_product_outputs, ProductOutputError};

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
    BatchManifest(prokadiff_report::BatchManifestError),
    Evidence(prokadiff_evidence::EvidenceError),
    Gd(prokadiff_gd::GdError),
    Intended(prokadiff_classify::IntendedError),
    Differential(prokadiff_classify::DifferentialError),
    AssociationProjection(prokadiff_classify::AssociationProjectionError),
    ProductOutput(ProductOutputError),
    Io(std::io::Error),
}

impl From<CliError> for RunError {
    fn from(e: CliError) -> Self {
        Self::Cli(e)
    }
}

impl From<prokadiff_report::BatchManifestError> for RunError {
    fn from(e: prokadiff_report::BatchManifestError) -> Self {
        Self::BatchManifest(e)
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

impl From<prokadiff_classify::DifferentialError> for RunError {
    fn from(e: prokadiff_classify::DifferentialError) -> Self {
        Self::Differential(e)
    }
}

impl From<prokadiff_classify::AssociationProjectionError> for RunError {
    fn from(e: prokadiff_classify::AssociationProjectionError) -> Self {
        Self::AssociationProjection(e)
    }
}

impl From<ProductOutputError> for RunError {
    fn from(e: ProductOutputError) -> Self {
        Self::ProductOutput(e)
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
            Self::BatchManifest(e) => write!(f, "{e}"),
            Self::Evidence(e) => write!(f, "{e}"),
            Self::Gd(e) => write!(f, "{e}"),
            Self::Intended(e) => write!(f, "{e}"),
            Self::Differential(e) => write!(f, "{e}"),
            Self::AssociationProjection(e) => write!(f, "{e}"),
            Self::ProductOutput(e) => write!(f, "{e}"),
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
        jc_min_support_reads: args.jc_min_support_reads,
        jc_min_frequency: args.jc_min_frequency,
        ..EngineOptions::default()
    };
    let gd = run_sample(&ref_fa, &reads, &args.outdir, &opts)?;
    info!("wrote {}", gd.display());
    Ok(())
}

fn run_product(job: ProductJob) -> Result<(), RunError> {
    if let Some(manifest_path) = job.batch_manifest.clone() {
        return run_batch_product(job, &manifest_path);
    }

    std::fs::create_dir_all(&job.outdir)?;
    let work = job.outdir.join("work");
    std::fs::create_dir_all(&work)?;
    let ref_fa = materialize_ref(&job.refs, &work)?;
    let repeats = load_repeats(&job.refs);
    let opts = EngineOptions {
        threads: job.threads,
        keep_bam: job.keep_bam,
        repeats,
        jc_min_support_reads: job.jc_min_support_reads,
        jc_min_frequency: job.jc_min_frequency,
        ..EngineOptions::default()
    };

    let starter_dir = job.outdir.join("starter");
    let edited_dir = job.outdir.join("edited");
    let (starter_gd_path, starter_pileups) = if job.subclonal {
        let (path, piles) = prokadiff_evidence::run_sample_with_pileup(
            &ref_fa,
            &FastqInput {
                files: job.starter.clone(),
            },
            &starter_dir,
            &opts,
        )?;
        (path, Some(piles))
    } else {
        let path = run_sample(
            &ref_fa,
            &FastqInput {
                files: job.starter.clone(),
            },
            &starter_dir,
            &opts,
        )?;
        (path, None)
    };
    let (edited_gd_path, edited_pileups) = if job.subclonal {
        let (path, piles) = prokadiff_evidence::run_sample_with_pileup(
            &ref_fa,
            &FastqInput {
                files: job.edited.clone(),
            },
            &edited_dir,
            &opts,
        )?;
        (path, Some(piles))
    } else {
        let path = run_sample(
            &ref_fa,
            &FastqInput {
                files: job.edited.clone(),
            },
            &edited_dir,
            &opts,
        )?;
        (path, None)
    };

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
            near_distance: job.offtarget_association_window,
            max_mismatches: DEFAULT_MAX_MISMATCHES,
            max_dna_bulge: job.max_dna_bulge,
            max_rna_bulge: job.max_rna_bulge,
            hypothesis: job.hypothesis,
        },
    )?;

    let sample_meta = prokadiff_classify::SampleMetadata {
        starter_names: job
            .starter
            .iter()
            .map(|p| {
                p.file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_default()
            })
            .collect(),
        edited_names: job
            .edited
            .iter()
            .map(|p| {
                p.file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_default()
            })
            .collect(),
        reference_names: job
            .refs
            .iter()
            .map(|p| {
                p.file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_default()
            })
            .collect(),
        editor: editor_kind.as_str().to_string(),
        spacer: job.spacer.clone(),
        pam: job.pam.clone(),
        threads: job.threads,
    };

    let ref_hashes: Vec<String> = job
        .refs
        .iter()
        .map(|p| compute_file_sha256(p).unwrap_or_else(|_| "ERROR_READING_REF".to_string()))
        .collect();
    let reference_sha256 = if ref_hashes.is_empty() {
        None
    } else {
        Some(ref_hashes.join(";"))
    };

    let bowtie2_version = Some(detect_bowtie2_version());

    let offtarget_search_status = if job.spacer.is_some() && job.editor != Editor::Dsb {
        prokadiff_classify::ValidationStatus::FixtureValidated
            .as_str()
            .to_string()
    } else {
        prokadiff_classify::ValidationStatus::NotRequested
            .as_str()
            .to_string()
    };
    let bulge_search_status = if job.max_dna_bulge > 0 || job.max_rna_bulge > 0 {
        prokadiff_classify::ValidationStatus::ExperimentalUnvalidated
            .as_str()
            .to_string()
    } else {
        "EXACT_UNGAPPED".to_string()
    };

    let provenance = prokadiff_classify::AnalysisProvenance {
        prokadiff_version: env!("CARGO_PKG_VERSION").to_string(),
        git_commit: option_env!("GIT_HASH").unwrap_or("UNKNOWN").to_string(),
        reference_sha256,
        bowtie2_version,
        offtarget_search_status,
        cfd_scoring_status: prokadiff_classify::ValidationStatus::Disabled
            .as_str()
            .to_string(),
        hsu_scoring_status: prokadiff_classify::ValidationStatus::Disabled
            .as_str()
            .to_string(),
        bulge_search_status,
        run_timestamp: utc_now_iso8601(),
    };

    let mut gbk_features = Vec::new();
    for ref_path in &job.refs {
        if let Ok(feats) = prokadiff_evidence::parse_genbank_features(ref_path) {
            gbk_features.extend(feats);
        }
    }
    let features: Vec<prokadiff_classify::AnnotatedFeature> = gbk_features
        .into_iter()
        .map(|f| prokadiff_classify::AnnotatedFeature {
            seq_id: f.seq_id,
            start: f.start,
            end: f.end,
            strand: f.strand,
            feature_type: f.feature_type,
            locus_tag: f.locus_tag,
            gene_name: f.gene_name,
            product: f.product,
        })
        .collect();

    let audit_result = prokadiff_classify::build_complete_audit_result(
        sample_meta,
        job.intended.is_some(),
        classified
            .intended_edit_assessments
            .clone()
            .unwrap_or_default(),
        classified.intended_event_ids.clone(),
        classified.starter_vs_ref,
        job.hypothesis,
        &classified.unintended,
        &classified.differential_events,
        &classified.candidate_search.sites,
        &classified.associations,
        classified.candidate_search.status,
        provenance,
        &features,
        &refs,
    )?;

    let audit_result = if job.subclonal {
        let contig_names: Vec<String> = refs.iter().map(|r| r.name.clone()).collect();
        let s_piles = starter_pileups.unwrap_or_default();
        let e_piles = edited_pileups.unwrap_or_default();
        let candidates = prokadiff_classify::subclonal::extract_subclonal_candidates_from_pileups(
            &contig_names,
            &s_piles,
            &e_piles,
            job.min_subclonal_af,
            job.min_subclonal_alt_reads,
        );
        let subclonal_params = prokadiff_classify::subclonal::SubclonalParams {
            min_depth: 20,
            min_alt_reads: job.min_subclonal_alt_reads,
            min_af: job.min_subclonal_af,
            max_af: 0.7499,
            min_delta_af: 0.05,
            max_strand_bias_pvalue: job.max_strand_bias_pvalue,
            fdr_alpha: 0.01,
            error_floor_snp: 0.005,
            error_floor_indel: 0.010,
        };
        let evaluated = prokadiff_classify::subclonal::evaluate_subclonal_candidates(
            candidates,
            &subclonal_params,
        );
        audit_result
            .with_subclonal_enabled(true)
            .with_subclonal_variants(evaluated)
    } else {
        audit_result
    };

    write_product_outputs(&audit_result, job.schema_version, &job.outdir)?;

    info!("wrote {}", starter_out.display());
    info!("wrote {}", edited_out.display());
    info!("wrote product reports under {}", job.outdir.display());
    Ok(())
}

fn run_batch_product(job: ProductJob, manifest_path: &Path) -> Result<(), RunError> {
    let manifest = prokadiff_report::BatchManifest::parse_file(manifest_path)?;
    std::fs::create_dir_all(&job.outdir)?;
    let work = job.outdir.join("work");
    std::fs::create_dir_all(&work)?;
    let ref_fa = materialize_ref(&job.refs, &work)?;
    let repeats = load_repeats(&job.refs);
    let opts = EngineOptions {
        threads: job.threads,
        keep_bam: job.keep_bam,
        repeats,
        jc_min_support_reads: job.jc_min_support_reads,
        jc_min_frequency: job.jc_min_frequency,
        ..EngineOptions::default()
    };

    // 1. Process Starter Strain
    let starter_dir = job.outdir.join("starter");
    let (starter_gd_path, starter_pileups) = if job.subclonal {
        let (path, piles) = prokadiff_evidence::run_sample_with_pileup(
            &ref_fa,
            &FastqInput {
                files: job.starter.clone(),
            },
            &starter_dir,
            &opts,
        )?;
        (path, Some(piles))
    } else {
        let path = run_sample(
            &ref_fa,
            &FastqInput {
                files: job.starter.clone(),
            },
            &starter_dir,
            &opts,
        )?;
        (path, None)
    };
    let starter_out = job.outdir.join("starter.gd");
    std::fs::copy(&starter_gd_path, &starter_out)?;
    let starter_gd = GenomeDiff::from_path(&starter_out)?;

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

    let ref_hashes: Vec<String> = job
        .refs
        .iter()
        .map(|p| compute_file_sha256(p).unwrap_or_else(|_| "ERROR_READING_REF".to_string()))
        .collect();
    let reference_sha256 = if ref_hashes.is_empty() {
        None
    } else {
        Some(ref_hashes.join(";"))
    };
    let bowtie2_version = Some(detect_bowtie2_version());

    let mut gbk_features = Vec::new();
    for ref_path in &job.refs {
        if let Ok(feats) = prokadiff_evidence::parse_genbank_features(ref_path) {
            gbk_features.extend(feats);
        }
    }
    let features: Vec<prokadiff_classify::AnnotatedFeature> = gbk_features
        .into_iter()
        .map(|f| prokadiff_classify::AnnotatedFeature {
            seq_id: f.seq_id,
            start: f.start,
            end: f.end,
            strand: f.strand,
            feature_type: f.feature_type,
            locus_tag: f.locus_tag,
            gene_name: f.gene_name,
            product: f.product,
        })
        .collect();

    let mut clone_ids = Vec::with_capacity(manifest.entries.len());
    let mut clone_summaries = Vec::with_capacity(manifest.entries.len());
    let mut all_observations = Vec::new();

    let clones_dir = job.outdir.join("clones");
    std::fs::create_dir_all(&clones_dir)?;

    for entry in &manifest.entries {
        clone_ids.push(entry.clone_id.clone());
        let clone_outdir = clones_dir.join(&entry.clone_id);
        std::fs::create_dir_all(&clone_outdir)?;

        let mut clone_reads = vec![entry.r1.clone()];
        if let Some(r2) = &entry.r2 {
            clone_reads.push(r2.clone());
        }

        let clone_evidence = clone_outdir.join("evidence");
        let (clone_gd_path, clone_pileups) = if job.subclonal {
            let (path, piles) = prokadiff_evidence::run_sample_with_pileup(
                &ref_fa,
                &FastqInput {
                    files: clone_reads.clone(),
                },
                &clone_evidence,
                &opts,
            )?;
            (path, Some(piles))
        } else {
            let path = run_sample(
                &ref_fa,
                &FastqInput {
                    files: clone_reads.clone(),
                },
                &clone_evidence,
                &opts,
            )?;
            (path, None)
        };
        let clone_gd_out = clone_outdir.join("edited.gd");
        std::fs::copy(&clone_gd_path, &clone_gd_out)?;
        let clone_starter_gd_out = clone_outdir.join("starter.gd");
        std::fs::copy(&starter_out, &clone_starter_gd_out)?;

        let edited_gd = GenomeDiff::from_path(&clone_gd_out)?;

        // Per-clone intended override or fallback to global intended
        let clone_intended_path = entry.intended_tsv.as_ref().or(job.intended.as_ref());
        let intended = match clone_intended_path {
            Some(p) => parse_intended_path(p)?,
            None => Vec::new(),
        };

        let clone_spacer = entry.spacer.as_ref().or(job.spacer.as_ref()).cloned();
        if matches!(job.editor, Editor::Cas9 | Editor::Cas12a) && clone_spacer.is_none() {
            return Err(RunError::Cli(CliError::InvalidSpacer {
                reason: format!(
                    "clone '{}' is missing a required spacer sequence for editor {:?}",
                    entry.clone_id, job.editor
                ),
            }));
        }
        let clone_pam = entry.pam.as_ref().or(job.pam.as_ref()).cloned();

        let classified = classify(
            &edited_gd,
            &starter_gd,
            &intended,
            &refs,
            &ClassifyOptions {
                editor: editor_kind,
                spacer: clone_spacer.clone(),
                pam: clone_pam.clone(),
                near_distance: job.offtarget_association_window,
                max_mismatches: DEFAULT_MAX_MISMATCHES,
                max_dna_bulge: job.max_dna_bulge,
                max_rna_bulge: job.max_rna_bulge,
                hypothesis: job.hypothesis,
            },
        )?;

        let sample_meta = prokadiff_classify::SampleMetadata {
            starter_names: job
                .starter
                .iter()
                .map(|p| {
                    p.file_name()
                        .map(|n| n.to_string_lossy().to_string())
                        .unwrap_or_default()
                })
                .collect(),
            edited_names: clone_reads
                .iter()
                .map(|p| {
                    p.file_name()
                        .map(|n| n.to_string_lossy().to_string())
                        .unwrap_or_default()
                })
                .collect(),
            reference_names: job
                .refs
                .iter()
                .map(|p| {
                    p.file_name()
                        .map(|n| n.to_string_lossy().to_string())
                        .unwrap_or_default()
                })
                .collect(),
            editor: editor_kind.as_str().to_string(),
            spacer: clone_spacer.clone(),
            pam: clone_pam.clone(),
            threads: job.threads,
        };

        let offtarget_search_status = if clone_spacer.is_some() && job.editor != Editor::Dsb {
            prokadiff_classify::ValidationStatus::FixtureValidated
                .as_str()
                .to_string()
        } else {
            prokadiff_classify::ValidationStatus::NotRequested
                .as_str()
                .to_string()
        };
        let bulge_search_status = if job.max_dna_bulge > 0 || job.max_rna_bulge > 0 {
            prokadiff_classify::ValidationStatus::ExperimentalUnvalidated
                .as_str()
                .to_string()
        } else {
            "EXACT_UNGAPPED".to_string()
        };

        let provenance = prokadiff_classify::AnalysisProvenance {
            prokadiff_version: env!("CARGO_PKG_VERSION").to_string(),
            git_commit: option_env!("GIT_HASH").unwrap_or("UNKNOWN").to_string(),
            reference_sha256: reference_sha256.clone(),
            bowtie2_version: bowtie2_version.clone(),
            offtarget_search_status,
            cfd_scoring_status: prokadiff_classify::ValidationStatus::Disabled
                .as_str()
                .to_string(),
            hsu_scoring_status: prokadiff_classify::ValidationStatus::Disabled
                .as_str()
                .to_string(),
            bulge_search_status,
            run_timestamp: utc_now_iso8601(),
        };

        let audit_result = prokadiff_classify::build_complete_audit_result(
            sample_meta,
            clone_intended_path.is_some(),
            classified
                .intended_edit_assessments
                .clone()
                .unwrap_or_default(),
            classified.intended_event_ids.clone(),
            classified.starter_vs_ref,
            job.hypothesis,
            &classified.unintended,
            &classified.differential_events,
            &classified.candidate_search.sites,
            &classified.associations,
            classified.candidate_search.status,
            provenance,
            &features,
            &refs,
        )?;

        let audit_result = if job.subclonal {
            let contig_names: Vec<String> = refs.iter().map(|r| r.name.clone()).collect();
            let s_piles = starter_pileups.as_deref().unwrap_or_default();
            let c_piles = clone_pileups.unwrap_or_default();
            let candidates =
                prokadiff_classify::subclonal::extract_subclonal_candidates_from_pileups(
                    &contig_names,
                    s_piles,
                    &c_piles,
                    job.min_subclonal_af,
                    job.min_subclonal_alt_reads,
                );
            let subclonal_params = prokadiff_classify::subclonal::SubclonalParams {
                min_depth: 20,
                min_alt_reads: job.min_subclonal_alt_reads,
                min_af: job.min_subclonal_af,
                max_af: 0.7499,
                min_delta_af: 0.05,
                max_strand_bias_pvalue: job.max_strand_bias_pvalue,
                fdr_alpha: 0.01,
                error_floor_snp: 0.005,
                error_floor_indel: 0.010,
            };
            let evaluated = prokadiff_classify::subclonal::evaluate_subclonal_candidates(
                candidates,
                &subclonal_params,
            );
            audit_result
                .with_subclonal_enabled(true)
                .with_subclonal_variants(evaluated)
        } else {
            audit_result
        };

        write_product_outputs(&audit_result, job.schema_version, &clone_outdir)?;

        let obs = prokadiff_report::extract_clone_observations(&entry.clone_id, &audit_result, &refs);
        all_observations.extend(obs);

        clone_summaries.push(prokadiff_report::CloneBatchSummary::from_audit(
            entry.clone_id.clone(),
            entry.r1.display().to_string(),
            entry.r2.as_ref().map(|p| p.display().to_string()),
            &audit_result,
        ));
    }

    let variants_matrix = prokadiff_report::align_batch_variants(&clone_ids, &all_observations);
    let batch_variants_tsv = job.outdir.join("batch_variants.tsv");
    prokadiff_report::write_batch_variants_tsv(&batch_variants_tsv, &variants_matrix, &clone_ids)?;

    let starter_sample_id = job
        .starter
        .first()
        .and_then(|p| p.file_name())
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| "starter".to_string());
    let reference_name = refs
        .first()
        .map(|r| r.name.clone())
        .unwrap_or_else(|| "reference".to_string());

    let batch_summary = prokadiff_report::BatchAuditSummary {
        manifest_path: manifest_path.display().to_string(),
        starter_sample_id,
        starter_r1: job
            .starter
            .first()
            .map(|p| p.display().to_string())
            .unwrap_or_default(),
        starter_r2: job.starter.get(1).map(|p| p.display().to_string()),
        reference_name,
        clones: clone_summaries,
        variants_matrix,
    };

    let batch_report_md = job.outdir.join("batch_report.md");
    prokadiff_report::write_batch_report_md(&batch_report_md, &batch_summary)?;

    info!("wrote batch report {}", batch_report_md.display());
    info!("wrote batch variants {}", batch_variants_tsv.display());
    info!("completed batch audit across {} clones", manifest.entries.len());

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

fn utc_now_iso8601() -> String {
    let now = std::time::SystemTime::now();
    let dur = now
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = dur.as_secs();
    let days = (secs / 86400) as i64 + 719468;
    let era = (if days >= 0 { days } else { days - 146096 }) / 146097;
    let doe = (days - era * 146097) as u32;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = (yoe as i64) + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    let time_of_day = secs % 86400;
    let hours = time_of_day / 3600;
    let minutes = (time_of_day % 3600) / 60;
    let seconds = time_of_day % 60;
    format!("{y:04}-{m:02}-{d:02}T{hours:02}:{minutes:02}:{seconds:02}Z")
}

fn compute_file_sha256(path: &std::path::Path) -> std::io::Result<String> {
    use sha2::{Digest, Sha256};
    let mut file = std::fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = std::io::Read::read(&mut file, &mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn detect_bowtie2_version() -> String {
    match std::process::Command::new("bowtie2")
        .arg("--version")
        .output()
    {
        Ok(output) if output.status.success() => {
            let stdout = String::from_utf8_lossy(&output.stdout);
            stdout
                .lines()
                .next()
                .map(|line| line.trim().to_string())
                .unwrap_or_else(|| "UNAVAILABLE".to_string())
        }
        _ => "UNAVAILABLE".to_string(),
    }
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
