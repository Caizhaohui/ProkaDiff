//! Single-sample consensus engine: Bowtie2 → BAM (noodles) → RA / MC / JC → Genome Diff.

pub mod bam_io;
pub(crate) mod emit;
pub(crate) mod jc_cluster;
#[cfg(test)]
mod tests;

use std::path::{Path, PathBuf};

use prokadiff_gd::GenomeDiff;
use rayon::prelude::*;
use tracing::info;

use crate::align::{align_to_bam, AlignKind, FastqInput};
use crate::error::Result;
use crate::fasta::{read_reference, FastaRecord, RepeatRegion};
use crate::mc::MC_DEL_MIN_LEN;
use crate::pileup::{
    apply_read, place_softclips_batch_with_repeat_evidence, AlignedRead, PlaceIndex, SplitCandidate,
};
use crate::ra::{PileupColumn, RaOptions};
use crate::repeat_ambiguous::{RepeatAmbiguousDiagnostics, RepeatAmbiguousSeed};

use bam_io::{read_aligned_bam, read_primary_bam};
use emit::emit_from_pileup;
use jc_cluster::second_pass_splits;

pub use jc_cluster::JC_CLUSTER_TOL_BP;

#[derive(Clone, Debug)]
pub struct EngineOptions {
    pub threads: usize,
    pub keep_bam: bool,
    pub ra: RaOptions,
    pub mc_min_len: usize,
    /// Unique-mapping + total-depth-0 gap must be at least this long to become
    /// a consensus DEL without JC support (see `docs/parity.md`).
    pub mc_del_min_len: usize,
    /// Known repeat regions (e.g. mobile elements) from reference annotations.
    pub repeats: Vec<RepeatRegion>,
    /// Minimum total supporting reads (split + spanning) required to emit a JC record (default: 3).
    pub jc_min_support_reads: usize,
    /// Minimum allele frequency relative to local coverage (default: 0.05, i.e. 5%).
    pub jc_min_frequency: f64,
}

impl Default for EngineOptions {
    fn default() -> Self {
        Self {
            threads: 8,
            keep_bam: false,
            ra: RaOptions::default(),
            mc_min_len: 3,
            mc_del_min_len: MC_DEL_MIN_LEN,
            repeats: Vec::new(),
            jc_min_support_reads: 3,
            jc_min_frequency: 0.05,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ContigPileup {
    pub columns: Vec<PileupColumn>,
    pub unique_depth: Vec<u32>,
    pub total_depth: Vec<u32>,
    pub splits: Vec<SplitCandidate>,
}

pub(crate) struct PileupResult {
    pub(crate) contigs: Vec<ContigPileup>,
    pub(crate) repeat_seeds: Vec<RepeatAmbiguousSeed>,
}

/// Align `reads` to `ref_fa` and write a Genome Diff to `outdir/output.gd`.
pub fn run_sample(
    ref_fa: &Path,
    reads: &FastqInput,
    outdir: &Path,
    opts: &EngineOptions,
) -> Result<PathBuf> {
    let mut opts_buf;
    let opts = if opts.repeats.is_empty() {
        let mut reps = Vec::new();
        if let Ok(r) = crate::fasta::parse_genbank_repeats(ref_fa) {
            reps.extend(r);
        }
        if reps.is_empty() {
            if let Some(parent) = ref_fa.parent() {
                if let Ok(entries) = std::fs::read_dir(parent) {
                    for entry in entries.flatten() {
                        let p = entry.path();
                        if let Some(ext) = p.extension().and_then(|s| s.to_str()) {
                            if ext.eq_ignore_ascii_case("gbk") || ext.eq_ignore_ascii_case("gb") {
                                if let Ok(r) = crate::fasta::parse_genbank_repeats(&p) {
                                    reps.extend(r);
                                }
                            }
                        }
                    }
                }
            }
        }
        if !reps.is_empty() {
            opts_buf = opts.clone();
            opts_buf.repeats = reps;
            &opts_buf
        } else {
            opts
        }
    } else {
        opts
    };
    std::fs::create_dir_all(outdir)?;
    let work = outdir.join("work");
    std::fs::create_dir_all(&work)?;
    let bam_path = outdir.join("aligned.bam");
    info!("prokadiff: primary alignment");
    align_to_bam(
        ref_fa,
        reads,
        &bam_path,
        opts.threads,
        &work,
        AlignKind::Primary,
    )?;
    let fasta = read_reference(ref_fa)?;
    info!(
        "prokadiff: pileup + junction seeds (place S>={})",
        crate::jc_seq::MIN_CLIP_FOR_PLACE
    );
    let primary_data = read_primary_bam(&bam_path, &fasta)?;
    let pileup = pileup_contigs_with_repeat(&fasta, &primary_data.aligned, opts);
    info!("prokadiff: candidate-junction second pass");
    let extra = second_pass_splits(&fasta, reads, &primary_data, &pileup.contigs, opts, &work)?;
    let gd = emit_from_pileup(&fasta, pileup.contigs, opts, &extra.splits);
    let gd_path = outdir.join("output.gd");
    gd.write_path(&gd_path)?;
    let mut repeat_seeds = pileup.repeat_seeds;
    repeat_seeds.extend(extra.repeat_seeds);
    let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(&fasta, repeat_seeds);
    let sample_identity = reads
        .files
        .iter()
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>()
        .join(",");
    diagnostics.write_tsv(
        &outdir.join("repeat_ambiguous_junctions.tsv"),
        &sample_identity,
    )?;
    if !opts.keep_bam {
        let _ = std::fs::remove_file(&bam_path);
        let _ = std::fs::remove_dir_all(&work);
    }
    Ok(gd_path)
}

pub fn call_from_bam(
    bam_path: &Path,
    fasta: &[FastaRecord],
    opts: &EngineOptions,
) -> Result<GenomeDiff> {
    let aligned = read_aligned_bam(bam_path, fasta)?;
    Ok(call_from_aligned(fasta, &aligned, opts))
}

/// In-memory consensus (no BAM). Used by `call_from_bam` and layer-0 tests.
pub fn call_from_aligned(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
) -> GenomeDiff {
    call_from_aligned_extra(fasta, aligned, opts, &[])
}

/// Like `call_from_aligned`, but also returns repeat-ambiguous diagnostics.
pub fn call_from_aligned_with_diagnostics(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
) -> (GenomeDiff, RepeatAmbiguousDiagnostics) {
    call_from_aligned_with_extra_and_diagnostics(fasta, aligned, opts, &[], Vec::new())
}

/// Like `call_from_aligned_extra`, but accepts both extra split candidates and
/// extra repeat seeds, returning GenomeDiff alongside complete repeat-ambiguous diagnostics.
pub fn call_from_aligned_with_extra_and_diagnostics(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
    extra_splits: &[SplitCandidate],
    extra_repeat_seeds: Vec<RepeatAmbiguousSeed>,
) -> (GenomeDiff, RepeatAmbiguousDiagnostics) {
    let pileup = pileup_contigs_with_repeat(fasta, aligned, opts);
    let mut repeat_seeds = pileup.repeat_seeds;
    repeat_seeds.extend(extra_repeat_seeds);
    let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(fasta, repeat_seeds);
    let gd = emit_from_pileup(fasta, pileup.contigs, opts, extra_splits);
    (gd, diagnostics)
}

/// Like `call_from_aligned`, then merge extra split candidates (second-pass
/// junction hits) before clustering and `accept_junction`.
pub fn call_from_aligned_extra(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
    extra_splits: &[SplitCandidate],
) -> GenomeDiff {
    let pileup = pileup_contigs(fasta, aligned, opts);
    emit_from_pileup(fasta, pileup, opts, extra_splits)
}

pub(crate) fn pileup_contigs(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
) -> Vec<ContigPileup> {
    pileup_contigs_with_repeat(fasta, aligned, opts).contigs
}

fn pileup_contigs_with_repeat(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
) -> PileupResult {
    let pool = if opts.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(opts.threads)
            .build()
            .ok()
    } else {
        None
    };

    let call_contigs = || {
        let (mut results, all_clips): (Vec<_>, Vec<_>) = fasta
            .par_iter()
            .enumerate()
            .map(|(idx, rec)| {
                let n = rec.seq.len();
                let mut columns: Vec<PileupColumn> = rec
                    .seq
                    .iter()
                    .map(|&b| PileupColumn {
                        ref_base: b.to_ascii_uppercase(),
                        observations: Vec::new(),
                        insertions: Vec::new(),
                    })
                    .collect();
                let mut unique_depth = vec![0u32; n];
                let mut total_depth = vec![0u32; n];
                let mut splits = Vec::new();
                let mut clips = Vec::new();
                for read in aligned.iter().filter(|r| r.contig_idx == idx) {
                    apply_read(
                        read,
                        &mut columns,
                        &mut unique_depth,
                        &mut total_depth,
                        &mut splits,
                        &mut clips,
                        Some(&rec.seq),
                    );
                }
                (
                    ContigPileup {
                        columns,
                        unique_depth,
                        total_depth,
                        splits,
                    },
                    clips,
                )
            })
            .unzip();

        let has_any_clips = all_clips.iter().any(|c| !c.is_empty());
        let mut repeat_seeds = Vec::new();
        if has_any_clips {
            let place_index = PlaceIndex::from_fasta(fasta);
            let placed = place_softclips_batch_with_repeat_evidence(&all_clips, &place_index);
            for (cp, outcome) in results.iter_mut().zip(placed) {
                cp.splits.extend(outcome.splits);
                repeat_seeds.extend(outcome.repeat_seeds);
            }
        }

        PileupResult {
            contigs: results,
            repeat_seeds,
        }
    };
    match &pool {
        Some(p) => p.install(call_contigs),
        None => call_contigs(),
    }
}
