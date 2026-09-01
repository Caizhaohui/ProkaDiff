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
use crate::mc::{
    call_missing_coverage, promote_mc_to_del_with_jc, true_missing_core, MC_DEL_MIN_LEN,
};
use crate::normalize::{right_align_del, right_align_ins};
use crate::pileup::{apply_read, place_softclip, AlignedRead, CigarKind, CigarOp, SplitCandidate};
use crate::ra::{call_consensus, ConsensusCall, PileupColumn, RaOptions};

#[derive(Clone, Debug)]
pub struct EngineOptions {
    pub threads: usize,
    pub keep_bam: bool,
    pub ra: RaOptions,
    pub mc_min_len: usize,
    /// Unique-mapping + total-depth-0 gap must be at least this long to become
    /// a consensus DEL without JC support (see `docs/parity.md`).
    pub mc_del_min_len: usize,
}

type ContigPileup = (Vec<PileupColumn>, Vec<u32>, Vec<u32>, Vec<SplitCandidate>);
type JcKey = (usize, u64, usize, u64, i64);

impl Default for EngineOptions {
    fn default() -> Self {
        Self {
            threads: 8,
            keep_bam: false,
            ra: RaOptions::default(),
            mc_min_len: 3,
            mc_del_min_len: MC_DEL_MIN_LEN,
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
    let aligned = read_aligned_bam(bam_path, fasta)?;
    Ok(call_from_aligned(fasta, &aligned, opts))
}

fn read_aligned_bam(bam_path: &Path, fasta: &[FastaRecord]) -> Result<Vec<AlignedRead>> {
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
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(0);
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
    Ok(aligned)
}

/// In-memory consensus (no BAM). Used by `call_from_bam` and layer-0 tests.
pub fn call_from_aligned(
    fasta: &[FastaRecord],
    aligned: &[AlignedRead],
    opts: &EngineOptions,
) -> GenomeDiff {
    let pool = if opts.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(opts.threads)
            .build()
            .ok()
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
                    );
                }
                for hint in &clips {
                    splits.extend(place_softclip(hint, fasta));
                }
                (columns, unique_depth, total_depth, splits)
            })
            .collect()
    };
    let contig_results: Vec<ContigPileup> = match &pool {
        Some(p) => p.install(call_contigs),
        None => call_contigs(),
    };

    let mut gd = GenomeDiff::new();
    gd.metadata
        .push(("PROGRAM".into(), "prokdiff-evidence".into()));
    let mut next_id = 1u32;

    for (idx, rec) in fasta.iter().enumerate() {
        let (columns, unique_depth, total_depth, splits) = &contig_results[idx];
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

        let mut by_key: HashMap<JcKey, JunctionSupport> = HashMap::new();
        let mut jc_meta: HashMap<JcKey, (bool, bool)> = HashMap::new();
        for sp in splits {
            // Do not key on read strand: plus and minus must aggregate (accept_junction).
            let key = (
                sp.contig_idx,
                sp.side1_pos_1,
                sp.side2_contig_idx,
                sp.side2_pos_1,
                sp.overlap,
            );
            let e = by_key.entry(key).or_default();
            jc_meta
                .entry(key)
                .or_insert((sp.side1_minus, sp.side2_minus));
            let ov1 = sp.left.len().max(3);
            let ov2 = sp.right.len().max(3);
            // Best-read semantics: one unbalanced read must not kill the
            // junction (published accept rules are max-of-min per read).
            let m = ov1.min(ov2);
            e.best_min_overlap = e.best_min_overlap.max(m);
            if sp.minus {
                e.minus_reads += 1;
                e.minus_best_min = e.minus_best_min.max(m);
            } else {
                e.plus_reads += 1;
                e.plus_best_min = e.plus_best_min.max(m);
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
        let mut accepted_jc = Vec::new();
        let mut jc_endpoints_1 = Vec::new();
        for (key, support) in &by_key {
            if !accept_junction(support) {
                continue;
            }
            let (c1, p1, c2, p2, ov) = *key;
            let (m1, m2) = jc_meta[key];
            accepted_jc.push((c1, p1, m1, c2, p2, m2, ov));
            if c1 == idx {
                jc_endpoints_1.push(p1);
            }
            if c2 == idx {
                jc_endpoints_1.push(p2);
            }
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
            let mc_id = next_id;
            next_id += 1;
            gd.entries
                .push(GdEntry::mc(mc_id, rec.name.clone(), start, end, 0, 0));
            if !promote_mc_to_del_with_jc(
                span,
                unique_depth,
                total_depth,
                opts.mc_del_min_len,
                opts.mc_min_len,
                &jc_endpoints_1,
            ) {
                continue;
            }
            let core = true_missing_core(
                span,
                unique_depth,
                total_depth,
                crate::mc::coverage_cutoff(unique_depth),
            )
            .unwrap_or(span);
            let del_start = core.start_1();
            let del_end = core.end_1();
            let size = del_end.saturating_sub(del_start) + 1;
            gd.entries
                .push(GdEntry::del(next_id, rec.name.clone(), del_start, size));
            if let Some(last) = gd.entries.last_mut() {
                last.parent_ids = vec![mc_id];
            }
            next_id += 1;
        }

        for (c1, p1, m1, c2, p2, m2, ov) in accepted_jc {
            let s1 = if m1 { "-" } else { "+" };
            let s2 = if m2 { "-" } else { "+" };
            let n1 = fasta[c1].name.clone();
            let n2 = fasta[c2].name.clone();
            gd.entries
                .push(GdEntry::jc(next_id, n1, p1, s1, n2, p2, s2, ov));
            next_id += 1;
        }
    }

    gd
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
            apply_read(
                &read,
                &mut columns,
                &mut depth,
                &mut vec![0u32; n],
                &mut splits,
                &mut Vec::new(),
            );
        }
        assert_eq!(
            call_consensus(&columns[2], &RaOptions::default()),
            ConsensusCall::Snp { alt: b'T' }
        );
    }

    fn fasta_chr(seq: &[u8]) -> [FastaRecord; 1] {
        [FastaRecord {
            name: "chr".into(),
            seq: seq.to_vec(),
        }]
    }

    fn covering_reads(
        seq: &[u8],
        start: i64,
        len: usize,
        n_plus: usize,
        n_minus: usize,
        mapq: u8,
    ) -> Vec<AlignedRead> {
        let slice = seq[start as usize..start as usize + len].to_vec();
        let mut out = Vec::new();
        for minus in std::iter::repeat_n(false, n_plus).chain(std::iter::repeat_n(true, n_minus)) {
            out.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: start,
                minus,
                seq: slice.clone(),
                cigar: vec![CigarOp {
                    kind: CigarKind::Match,
                    len,
                }],
                mapq,
            });
        }
        out
    }

    fn gd_dels(gd: &GenomeDiff) -> Vec<&prokdiff_gd::GdEntry> {
        gd.entries
            .iter()
            .filter(|e| e.kind == prokdiff_gd::GdKind::Del)
            .collect()
    }

    fn opts_single_thread() -> EngineOptions {
        EngineOptions {
            threads: 1,
            ..EngineOptions::default()
        }
    }

    #[test]
    fn short_internal_coverage_gap_does_not_emit_del() {
        let seq = vec![b'A'; 40];
        let mut reads = covering_reads(&seq, 0, 15, 6, 6, UNIQUE_MAPQ);
        reads.extend(covering_reads(&seq, 20, 20, 6, 6, UNIQUE_MAPQ));
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(
            gd_dels(&gd).is_empty(),
            "short unique+total gap must stay MC, not consensus DEL; got {:?}",
            gd_dels(&gd).iter().map(|e| &e.fields).collect::<Vec<_>>()
        );
    }

    #[test]
    fn contig_terminal_zero_does_not_emit_del() {
        let seq = vec![b'A'; 40];
        let reads = covering_reads(&seq, 10, 30, 6, 6, UNIQUE_MAPQ);
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(gd_dels(&gd).is_empty());
    }

    #[test]
    fn repeat_like_low_mapq_coverage_is_not_del() {
        let seq = vec![b'A'; 40];
        let mut reads = covering_reads(&seq, 0, 10, 6, 6, UNIQUE_MAPQ);
        reads.extend(covering_reads(&seq, 10, 16, 6, 6, 3)); // MAPQ < unique
        reads.extend(covering_reads(&seq, 26, 14, 6, 6, UNIQUE_MAPQ));
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(
            gd_dels(&gd).is_empty(),
            "unique-depth 0 with remaining total coverage must not be DEL; got {:?}",
            gd_dels(&gd).iter().map(|e| &e.fields).collect::<Vec<_>>()
        );
    }

    #[test]
    fn long_internal_true_gap_emits_del() {
        let seq = vec![b'A'; 80];
        let mut reads = covering_reads(&seq, 0, 10, 6, 6, UNIQUE_MAPQ);
        reads.extend(covering_reads(&seq, 70, 10, 6, 6, UNIQUE_MAPQ));
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let dels = gd_dels(&gd);
        assert_eq!(
            dels.len(),
            1,
            "expected one unique-mapping DEL, got {dels:?}"
        );
        let start: u64 = dels[0].fields[1].parse().unwrap();
        let size: u64 = dels[0].fields[2].parse().unwrap();
        assert!(start > 1, "DEL must not be contig-terminal, start={start}");
        assert!(
            size >= 50,
            "unique+total-zero core must be long, size={size}"
        );
    }

    fn gd_jcs(gd: &GenomeDiff) -> Vec<&prokdiff_gd::GdEntry> {
        gd.entries
            .iter()
            .filter(|e| e.kind == prokdiff_gd::GdKind::Jc)
            .collect()
    }

    #[test]
    fn both_strand_long_deletion_emits_jc() {
        let seq = vec![b'A'; 50];
        let mut reads = Vec::new();
        for minus in [false, false, false, true, true, true] {
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus,
                seq: vec![b'A'; 30],
                cigar: vec![
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 15,
                    },
                    CigarOp {
                        kind: CigarKind::Del,
                        len: 8,
                    },
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 15,
                    },
                ],
                mapq: UNIQUE_MAPQ,
            });
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(
            !gd_jcs(&gd).is_empty(),
            "plus+minus long-D splits must aggregate into an accepted JC; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
    }

    fn rc(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                x => x,
            })
            .collect()
    }

    #[test]
    fn softclip_onto_repeat_copy_emits_jc() {
        // Unique flanks + two identical 20 bp IS-like copies (first-period mosaic).
        let mut seq = vec![b'C'; 160];
        let motif = *b"ACGTACGTACGTACGTACGT";
        seq[40..60].copy_from_slice(&motif);
        seq[100..120].copy_from_slice(&motif);
        // Plus read crossing the insert left boundary: match ref 71..90, 3′ clip.
        let mut plus_seq = vec![b'C'; 20];
        plus_seq.extend_from_slice(&motif);
        // Minus read crossing the SAME boundary: SEQ = rc([match][clip]) =
        // [rc(clip)][rc(match)], so the clip leads and the CIGAR is S+M over
        // the same reference match span (biologically consistent minus read).
        let mut minus_seq = rc(&motif);
        minus_seq.extend(std::iter::repeat_n(b'C', 20));
        let mut reads = Vec::new();
        for minus in [false, false, false, true, true, true] {
            let (s, cigar) = if minus {
                (
                    minus_seq.clone(),
                    vec![
                        CigarOp {
                            kind: CigarKind::SoftClip,
                            len: 20,
                        },
                        CigarOp {
                            kind: CigarKind::Match,
                            len: 20,
                        },
                    ],
                )
            } else {
                (
                    plus_seq.clone(),
                    vec![
                        CigarOp {
                            kind: CigarKind::Match,
                            len: 20,
                        },
                        CigarOp {
                            kind: CigarKind::SoftClip,
                            len: 20,
                        },
                    ],
                )
            };
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 70,
                minus,
                seq: s,
                cigar,
                mapq: UNIQUE_MAPQ,
            });
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert!(
            !jcs.is_empty(),
            "softclipped IS motif must place a JC onto a repeat copy; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let side2: u64 = jcs[0].fields[4].parse().unwrap();
        assert!(
            (41..=60).contains(&side2) || (101..=120).contains(&side2),
            "side2 should land in an IS copy, got {side2}"
        );
    }

    #[test]
    fn synth_is_mob_geometry_both_strand_softclips_aggregate() {
        // Regression for the synth_is_mob empty-GD failure (job 2371258):
        // plus reads with 3′ softclips and minus reads with 5′ softclips cross
        // the SAME insert left boundary. Before strand normalization the minus
        // hints keyed at first_match_ref+1 (e.g. 11) instead of the junction
        // base (30), so plus/minus never aggregated and no JC was accepted.
        let motif = *b"TGCAAGTCGATCGTTAGCCA";
        let mut seq = vec![b'C'; 160];
        seq[40..60].copy_from_slice(&motif);
        seq[100..120].copy_from_slice(&motif);
        // Junction: unique flank up to ref 30 (1-based), then motif content.
        let mut plus_seq = vec![b'C'; 20];
        plus_seq.extend_from_slice(&motif);
        let mut minus_seq = rc(&motif);
        minus_seq.extend(std::iter::repeat_n(b'C', 20));
        let mut reads = Vec::new();
        for _ in 0..3 {
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 10, // 20M covers ref 11..30 (1-based)
                minus: false,
                seq: plus_seq.clone(),
                cigar: vec![
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 20,
                    },
                    CigarOp {
                        kind: CigarKind::SoftClip,
                        len: 20,
                    },
                ],
                mapq: UNIQUE_MAPQ,
            });
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 10, // same reference match span, minus strand
                minus: true,
                seq: minus_seq.clone(),
                cigar: vec![
                    CigarOp {
                        kind: CigarKind::SoftClip,
                        len: 20,
                    },
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 20,
                    },
                ],
                mapq: UNIQUE_MAPQ,
            });
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            2,
            "one aggregated JC per motif copy; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let mut side2s = Vec::new();
        for jc in &jcs {
            let side1: u64 = jc.fields[1].parse().unwrap();
            let side2: u64 = jc.fields[4].parse().unwrap();
            assert_eq!(
                side1, 30,
                "both strands must key side1 at the junction base"
            );
            side2s.push(side2);
        }
        side2s.sort_unstable();
        assert_eq!(
            side2s,
            vec![41, 101],
            "side2 must be the two motif copy starts"
        );
    }

    #[test]
    fn jc_supported_cigar_gap_promotes_short_del() {
        let seq = vec![b'A'; 50];
        let mut reads = Vec::new();
        for minus in [false, false, false, true, true, true] {
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus,
                seq: vec![b'A'; 30],
                cigar: vec![
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 15,
                    },
                    CigarOp {
                        kind: CigarKind::Del,
                        len: 8,
                    },
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 15,
                    },
                ],
                mapq: UNIQUE_MAPQ,
            });
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        assert!(
            !gd_dels(&gd).is_empty(),
            "short true gap with both-flank JC must become DEL; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
    }
}
