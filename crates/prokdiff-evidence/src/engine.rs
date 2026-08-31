//! Single-sample consensus engine: Bowtie2 → BAM (noodles) → RA / MC / JC → Genome Diff.

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind as NoodlesKind;
use prokdiff_gd::{GdEntry, GenomeDiff};
use rayon::prelude::*;

use crate::align::{align_to_bam, FastqInput};
use crate::error::{EvidenceError, Result};
use crate::fasta::{read_reference, FastaRecord};
use crate::jc::{accept_junction, JunctionSupport};
use crate::mc::call_missing_coverage;
use crate::normalize::{right_align_del, right_align_ins};
use crate::pileup::{apply_read, AlignedRead, CigarKind, CigarOp, SplitCandidate};
use crate::ra::{call_consensus, ConsensusCall, PileupColumn, RaOptions};

#[derive(Clone, Debug)]
pub struct EngineOptions {
    pub threads: usize,
    pub keep_bam: bool,
    pub ra: RaOptions,
    pub mc_min_len: usize,
}

impl Default for EngineOptions {
    fn default() -> Self {
        Self {
            threads: 8,
            keep_bam: false,
            ra: RaOptions::default(),
            mc_min_len: 3,
        }
    }
}

/// Align `reads` to `ref_fa` and write a Genome Diff to `outdir/output.gd`.
pub fn run_sample(
    ref_fa: &Path,
    reads: &FastqInput,
    outdir: &Path,
    opts: &EngineOptions,
) -> Result<PathBuf> {
    std::fs::create_dir_all(outdir)?;
    let work = outdir.join("work");
    std::fs::create_dir_all(&work)?;
    let bam_path = outdir.join("aligned.bam");
    align_to_bam(ref_fa, reads, &bam_path, opts.threads, &work)?;
    let fasta = read_reference(ref_fa)?;
    let gd = call_from_bam(&bam_path, &fasta, opts)?;
    let gd_path = outdir.join("output.gd");
    gd.write_path(&gd_path)?;
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
    let name_to_idx: HashMap<&str, usize> = fasta
        .iter()
        .enumerate()
        .map(|(i, r)| (r.name.as_str(), i))
        .collect();

    let mut reader = File::open(bam_path).map(bam::io::Reader::new)?;
    let header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(n, _)| n.to_string())
        .collect();

    let mut aligned = Vec::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
        if rec.flags().is_unmapped() || rec.flags().is_secondary() {
            continue;
        }
        let Some(Ok(rid)) = rec.reference_sequence_id() else {
            continue;
        };
        let contig_name = ref_names.get(rid).map(String::as_str).unwrap_or("");
        let Some(&contig_idx) = name_to_idx.get(contig_name) else {
            continue;
        };
        let Some(Ok(start)) = rec.alignment_start() else {
            continue;
        };
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(0); // Option<MappingQuality>
        let minus = rec.flags().is_reverse_complemented();
        let seq_len = rec.sequence().len();
        let mut seq = Vec::with_capacity(seq_len);
        for i in 0..seq_len {
            if let Some(b) = rec.sequence().get(i) {
                seq.push(b);
            }
        }
        let mut cigar = Vec::new();
        for op in rec.cigar().iter() {
            let op = op.map_err(|e| EvidenceError::Alignment(e.to_string()))?;
            let kind = match op.kind() {
                NoodlesKind::Match | NoodlesKind::SequenceMatch | NoodlesKind::SequenceMismatch => {
                    CigarKind::Match
                }
                NoodlesKind::Insertion => CigarKind::Ins,
                NoodlesKind::Deletion => CigarKind::Del,
                NoodlesKind::SoftClip => CigarKind::SoftClip,
                NoodlesKind::HardClip => CigarKind::HardClip,
                NoodlesKind::Skip => CigarKind::Skip,
                _ => continue,
            };
            cigar.push(CigarOp {
                kind,
                len: op.len(),
            });
        }
        aligned.push(AlignedRead {
            contig_idx,
            ref_start_0: start.get() as i64 - 1,
            minus,
            seq,
            cigar,
            mapq,
        });
    }

    let pool = if opts.threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(opts.threads)
                .build()
                .map_err(|e| EvidenceError::Alignment(e.to_string()))?,
        )
    } else {
        None
    };

    let call_contigs = || {
        fasta
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
                let mut splits = Vec::new();
                for read in aligned.iter().filter(|r| r.contig_idx == idx) {
                    apply_read(read, &mut columns, &mut unique_depth, &mut splits);
                }
                (columns, unique_depth, splits)
            })
            .collect()
    };
    let contig_results: Vec<(Vec<PileupColumn>, Vec<u32>, Vec<SplitCandidate>)> = match &pool {
        Some(p) => p.install(call_contigs),
        None => call_contigs(),
    };

    let mut gd = GenomeDiff::new();
    gd.metadata
        .push(("PROGRAM".into(), "prokdiff-evidence".into()));
    let mut next_id = 1u32;

    for (idx, rec) in fasta.iter().enumerate() {
        let (columns, unique_depth, splits) = &contig_results[idx];
        let mut pending_del: Option<(u64, u8)> = None;
        let flush_del = |gd: &mut GenomeDiff, id: &mut u32, start: u64, size: u8| {
            // Match breseq mutation coordinates (RA evidence may stay 5′).
            let start = right_align_del(&rec.seq, start, u64::from(size));
            gd.entries
                .push(GdEntry::del(*id, rec.name.clone(), start, u64::from(size)));
            *id += 1;
        };

        for (i, col) in columns.iter().enumerate() {
            let pos = i as u64 + 1;
            match call_consensus(col, &opts.ra) {
                ConsensusCall::Snp { alt } => {
                    if let Some((s, sz)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz);
                    }
                    gd.entries.push(GdEntry::snp(
                        next_id,
                        rec.name.clone(),
                        pos,
                        (alt as char).to_string(),
                    ));
                    next_id += 1;
                }
                ConsensusCall::Ins { seq } => {
                    if let Some((s, sz)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz);
                    }
                    let (pos, seq) = right_align_ins(&rec.seq, pos, &seq);
                    let s = String::from_utf8_lossy(&seq).into_owned();
                    gd.entries
                        .push(GdEntry::ins(next_id, rec.name.clone(), pos, s));
                    next_id += 1;
                }
                ConsensusCall::Del { size } => match pending_del {
                    Some((s, sz)) if s + u64::from(sz) == pos && sz + size <= 2 => {
                        pending_del = Some((s, sz + size));
                    }
                    Some((s, sz)) => {
                        flush_del(&mut gd, &mut next_id, s, sz);
                        pending_del = Some((pos, size));
                    }
                    None => pending_del = Some((pos, size)),
                },
                _ => {
                    if let Some((s, sz)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz);
                    }
                }
            }
        }
        if let Some((s, sz)) = pending_del.take() {
            flush_del(&mut gd, &mut next_id, s, sz);
        }

        for span in call_missing_coverage(unique_depth, opts.mc_min_len) {
            if span.len() < 3 {
                continue;
            }
            // Contig ends look like unique-depth 0; breseq reports UN, not DEL.
            if span.start_0 == 0 || span.end_0 == unique_depth.len() {
                continue;
            }
            let start = span.start_1();
            let end = span.end_1();
            let size = end.saturating_sub(start) + 1;
            let mc_id = next_id;
            next_id += 1;
            gd.entries
                .push(GdEntry::mc(mc_id, rec.name.clone(), start, end, 0, 0));
            gd.entries
                .push(GdEntry::del(next_id, rec.name.clone(), start, size));
            if let Some(last) = gd.entries.last_mut() {
                last.parent_ids = vec![mc_id];
            }
            next_id += 1;
        }

        let mut by_key: HashMap<(u64, u64, i64, bool), JunctionSupport> = HashMap::new();
        for sp in splits {
            // Long I/D splits are JC candidates (IS / cassette must not be dropped).
            let key = (sp.side1_pos_1, sp.side2_pos_1, sp.overlap, sp.minus);
            let e = by_key.entry(key).or_default();
            let ov1 = sp.left.len().max(3);
            let ov2 = sp.right.len().max(3);
            if sp.minus {
                e.minus_reads += 1;
                e.minus_side1 = e.minus_side1.max(ov1);
                e.minus_side2 = e.minus_side2.max(ov2);
            } else {
                e.plus_reads += 1;
                e.plus_side1 = e.plus_side1.max(ov1);
                e.plus_side2 = e.plus_side2.max(ov2);
            }
            e.min_overlap_side1 = if e.min_overlap_side1 == 0 {
                ov1
            } else {
                e.min_overlap_side1.min(ov1)
            };
            e.min_overlap_side2 = if e.min_overlap_side2 == 0 {
                ov2
            } else {
                e.min_overlap_side2.min(ov2)
            };
        }
        for ((p1, p2, ov, minus), support) in by_key {
            if !accept_junction(&support) {
                continue;
            }
            let s1 = if minus { "-" } else { "+" };
            gd.entries.push(GdEntry::jc(
                next_id,
                rec.name.clone(),
                p1,
                s1,
                rec.name.clone(),
                p2,
                s1,
                ov,
            ));
            next_id += 1;
        }
    }

    Ok(gd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::{apply_read, UNIQUE_MAPQ};

    #[test]
    fn ra_calls_from_in_memory_pileup_without_bam() {
        let fasta = [FastaRecord {
            name: "chr".into(),
            seq: b"ACGTACGTAC".to_vec(),
        }];
        let n = fasta[0].seq.len();
        let mut columns: Vec<PileupColumn> = fasta[0]
            .seq
            .iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: Vec::new(),
                insertions: Vec::new(),
            })
            .collect();
        let mut depth = vec![0u32; n];
        let mut splits = Vec::new();
        let mut alt = fasta[0].seq.clone();
        alt[2] = b'T';
        for minus in [false, false, false, true, true, true] {
            let read = AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus,
                seq: alt.clone(),
                cigar: vec![CigarOp {
                    kind: CigarKind::Match,
                    len: alt.len(),
                }],
                mapq: UNIQUE_MAPQ,
            };
            apply_read(&read, &mut columns, &mut depth, &mut splits);
        }
        assert_eq!(
            call_consensus(&columns[2], &RaOptions::default()),
            ConsensusCall::Snp { alt: b'T' }
        );
    }
}
