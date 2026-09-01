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

/// Two split candidates whose side1 AND side2 both land within this many bp
/// are clustered into one junction before `accept_junction`. IS terminal /
/// TSD register ambiguity keys reads of one event a few bp apart (measured
/// in the synth_is_mob 2371412 BAM: plus reads of the same junction end at
/// 601 and 602; breseq's own two JC for the event land at 599 and 602,
/// 3 bp apart), so exact-key aggregation can strand reads of one junction
/// in different buckets. Tied to `jc_tol_bp=5` in
/// `scripts/layer2_gd_compare.py` (oracle-side matching tolerance); see
/// `docs/parity.md`.
pub const JC_CLUSTER_TOL_BP: u64 = 5;

/// One junction cluster: representative coordinates are those of the FIRST
/// member (smallest side1/side2 after the deterministic sort; the
/// breseq-side comparison matches greedily within ±5 bp, so any member
/// coordinate would do — first-member is the simplest deterministic choice).
/// `JunctionSupport` is aggregated over ALL members, merging strands across
/// the few-bp register split.
struct JcCluster {
    rep_side1: u64,
    rep_side2: u64,
    rep_overlap: i64,
    rep_side1_minus: bool,
    rep_side2_minus: bool,
    support: JunctionSupport,
}

impl JcCluster {
    fn new(sp: &SplitCandidate) -> Self {
        let mut c = Self {
            rep_side1: sp.side1_pos_1,
            rep_side2: sp.side2_pos_1,
            rep_overlap: sp.overlap,
            rep_side1_minus: sp.side1_minus,
            rep_side2_minus: sp.side2_minus,
            support: JunctionSupport::default(),
        };
        c.add(sp);
        c
    }

    fn add(&mut self, sp: &SplitCandidate) {
        let ov1 = sp.left.len().max(3);
        let ov2 = sp.right.len().max(3);
        // Best-read semantics: one unbalanced read must not kill the
        // junction (published accept rules are max-of-min per read).
        let m = ov1.min(ov2);
        let e = &mut self.support;
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
}

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

        // Junction clustering (replaces exact-key aggregation). Do not key
        // on read strand: plus and minus must aggregate (accept_junction).
        // Group by (contig, side2 contig, overlap sign) — the exact overlap
        // value is register-dependent, but its sign is biologically distinct
        // (deletion-like < 0, blunt = 0, insertion-like > 0). Within a group,
        // sort by (side1, side2) and greedily cluster: a candidate joins the
        // first cluster whose FIRST member is within JC_CLUSTER_TOL_BP on
        // both sides (single-pass against the first member is deterministic;
        // documented in docs/parity.md).
        let mut groups: HashMap<(usize, usize, i8), Vec<&SplitCandidate>> = HashMap::new();
        for sp in splits {
            groups
                .entry((
                    sp.contig_idx,
                    sp.side2_contig_idx,
                    sp.overlap.signum() as i8,
                ))
                .or_default()
                .push(sp);
        }
        let mut group_keys: Vec<(usize, usize, i8)> = groups.keys().copied().collect();
        group_keys.sort_unstable();
        let mut accepted_jc = Vec::new();
        let mut jc_endpoints_1 = Vec::new();
        for gk in group_keys {
            let mut members: Vec<&SplitCandidate> = groups[&gk].clone();
            // Stable sort: read order breaks ties, keeping the first member
            // (hence the representative) deterministic.
            members.sort_by_key(|sp| (sp.side1_pos_1, sp.side2_pos_1));
            let mut clusters: Vec<JcCluster> = Vec::new();
            for sp in members {
                let mut hit = None;
                for (i, c) in clusters.iter().enumerate() {
                    if sp.side1_pos_1.abs_diff(c.rep_side1) <= JC_CLUSTER_TOL_BP
                        && sp.side2_pos_1.abs_diff(c.rep_side2) <= JC_CLUSTER_TOL_BP
                    {
                        hit = Some(i);
                        break;
                    }
                }
                match hit {
                    Some(i) => clusters[i].add(sp),
                    None => clusters.push(JcCluster::new(sp)),
                }
            }
            for c in clusters {
                if !accept_junction(&c.support) {
                    continue;
                }
                accepted_jc.push((
                    gk.0,
                    c.rep_side1,
                    c.rep_side1_minus,
                    gk.1,
                    c.rep_side2,
                    c.rep_side2_minus,
                    c.rep_overlap,
                ));
                if gk.0 == idx {
                    jc_endpoints_1.push(c.rep_side1);
                }
                if gk.1 == idx {
                    jc_endpoints_1.push(c.rep_side2);
                }
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
            // breseq GD convention: side strands are 1 / -1.
            let s1 = if m1 { "-1" } else { "1" };
            let s2 = if m2 { "-1" } else { "1" };
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

    #[test]
    fn softclip_onto_repeat_copy_emits_jc() {
        // Unique flanks + two identical 20 bp IS-like copies (first-period mosaic).
        let mut seq = vec![b'C'; 160];
        let motif = *b"ACGTACGTACGTACGTACGT";
        seq[40..60].copy_from_slice(&motif);
        seq[100..120].copy_from_slice(&motif);
        // Reads crossing the insert left boundary: match ref 71..90, 3′ clip.
        // SAM/BAM stores SEQ reference-forward (minus reads are reverse-
        // complemented), so plus and minus reads crossing the same boundary
        // share the SAME stored SEQ and CIGAR; only the flag differs.
        let mut read_seq = vec![b'C'; 20];
        read_seq.extend_from_slice(&motif);
        let mut reads = Vec::new();
        for minus in [false, false, false, true, true, true] {
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 70,
                minus,
                seq: read_seq.clone(),
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
        // Regression for the synth_is_mob empty-GD failures (jobs 2371258,
        // 2371412): plus and minus reads with 3′ softclips cross the SAME
        // insert left boundary. SAM/BAM stores SEQ reference-forward, so both
        // strands share the same stored SEQ/CIGAR and must key side1 at the
        // junction base (30). Keying minus reads at the opposite end of the
        // match block (the pre-fix behavior) scattered minus candidates tens
        // of bp away, so plus/minus never aggregated and no JC was accepted.
        let motif = *b"TGCAAGTCGATCGTTAGCCA";
        let mut seq = vec![b'C'; 160];
        seq[40..60].copy_from_slice(&motif);
        seq[100..120].copy_from_slice(&motif);
        // Junction: unique flank up to ref 30 (1-based), then motif content.
        let mut read_seq = vec![b'C'; 20];
        read_seq.extend_from_slice(&motif);
        let mut reads = Vec::new();
        for _ in 0..3 {
            for minus in [false, true] {
                reads.push(AlignedRead {
                    contig_idx: 0,
                    ref_start_0: 10, // 20M covers ref 11..30 (1-based)
                    minus,
                    seq: read_seq.clone(),
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
            }
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
    fn few_bp_strand_split_junctions_cluster_into_one_jc() {
        // Regression for the synth_is_mob job 2371412 zero-JC failure. IS
        // terminal / TSD register ambiguity keys reads of one insertion
        // junction a few bp apart (measured in the 2371412 BAM: plus M+S
        // reads end at 601/602, and breseq's own two JC for the event land
        // at 599/602, 3 bp apart). Plus reads: M+S ending at 599. Minus
        // reads: same reference-forward stored SEQ, M+S ending at 602.
        // Exact-key aggregation put the registers in different buckets;
        // ±5 bp clustering must merge them into exactly ONE accepted JC.
        let motif = *b"TGCAAGTCGATCGTTAGCCA";
        let mut seq = vec![b'C'; 700];
        seq[100..120].copy_from_slice(&motif); // one existing copy -> side2 = 101
                                               // Stored SEQ (reference-forward on both strands): match = unique
                                               // flank, 3′ clip = motif content.
        let mut read_seq = vec![b'C'; 20];
        read_seq.extend_from_slice(&motif);
        let mut reads = Vec::new();
        for _ in 0..3 {
            // Plus reads: match ref 580..599 (1-based) -> side1 = 599.
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 579,
                minus: false,
                seq: read_seq.clone(),
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
            // Minus reads: match ref 583..602 (1-based) -> side1 = 602.
            reads.push(AlignedRead {
                contig_idx: 0,
                ref_start_0: 582,
                minus: true,
                seq: read_seq.clone(),
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
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            1,
            "599/602 strand split must cluster into exactly one JC; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let side1: u64 = jcs[0].fields[1].parse().unwrap();
        let side2: u64 = jcs[0].fields[4].parse().unwrap();
        assert_eq!(
            side1, 599,
            "representative = first member after the (side1, side2) sort"
        );
        assert_eq!(side2, 101, "clip must place onto the motif copy start");
    }

    #[test]
    fn junctions_beyond_cluster_tolerance_do_not_merge() {
        // Two junctions 10 bp apart (> JC_CLUSTER_TOL_BP), each with both
        // strands, must stay two clusters -> two JC.
        let motif = *b"TGCAAGTCGATCGTTAGCCA";
        let mut seq = vec![b'C'; 700];
        seq[100..120].copy_from_slice(&motif);
        let mut read_seq = vec![b'C'; 20];
        read_seq.extend_from_slice(&motif);
        let mut reads = Vec::new();
        // Junction A: plus ends 599 / minus ends 602.
        // Junction B: plus ends 609 / minus ends 612 (10 bp from A).
        for (plus_start0, minus_start0) in [(579i64, 582i64), (589, 592)] {
            for _ in 0..3 {
                for (start0, minus) in [(plus_start0, false), (minus_start0, true)] {
                    reads.push(AlignedRead {
                        contig_idx: 0,
                        ref_start_0: start0,
                        minus,
                        seq: read_seq.clone(),
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
                }
            }
        }
        let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
        let jcs = gd_jcs(&gd);
        assert_eq!(
            jcs.len(),
            2,
            "junctions 10 bp apart must not merge; entries={:?}",
            gd.entries
                .iter()
                .map(|e| (e.kind, e.fields.clone()))
                .collect::<Vec<_>>()
        );
        let mut side1s: Vec<u64> = jcs.iter().map(|jc| jc.fields[1].parse().unwrap()).collect();
        side1s.sort_unstable();
        assert_eq!(side1s, vec![599, 609]);
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
