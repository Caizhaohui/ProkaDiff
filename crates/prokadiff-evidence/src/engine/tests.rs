use std::collections::HashSet;

use prokadiff_gd::GenomeDiff;

use super::bam_io::*;
use super::*;
use crate::fasta::FastaRecord;
use crate::jc::SubAlignment;
use crate::pileup::{
    apply_read, AlignedRead, CigarKind, CigarOp, SplitCandidate, SplitOrigin, UNIQUE_MAPQ,
};
use crate::ra::{call_consensus, ConsensusCall, PileupColumn, RaOptions};

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
            deletions: Vec::new(),
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
            molecule_id: 0,
        };
        apply_read(
            &read,
            read.molecule_id,
            &mut columns,
            &mut depth,
            &mut vec![0u32; n],
            &mut splits,
            &mut Vec::new(),
            Some(&fasta[0].seq),
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
            molecule_id: 0,
        });
    }
    out
}

fn gd_dels(gd: &GenomeDiff) -> Vec<&prokadiff_gd::GdEntry> {
    gd.entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Del)
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
fn adjacent_snps_merge_into_sub() {
    // Reference 40 bp
    let mut ref_seq = vec![b'A'; 40];
    ref_seq[15] = b'C';
    ref_seq[16] = b'G';

    // Reads have TT at index 15..17 (1-based 16..17, size 2)
    let mut mut_seq = ref_seq.clone();
    mut_seq[15] = b'T';
    mut_seq[16] = b'T';

    let reads = covering_reads(&mut_seq, 0, 40, 6, 6, UNIQUE_MAPQ);
    let gd = call_from_aligned(&fasta_chr(&ref_seq), &reads, &opts_single_thread());

    let subs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Sub)
        .collect();
    assert_eq!(subs.len(), 1, "expected 1 SUB entry, got {:?}", gd.entries);
    assert_eq!(subs[0].fields[0], "chr");
    assert_eq!(subs[0].fields[1], "16"); // 1-based pos
    assert_eq!(subs[0].fields[2], "2"); // size
    assert_eq!(subs[0].fields[3], "TT"); // new_seq
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

fn gd_jcs(gd: &GenomeDiff) -> Vec<&prokadiff_gd::GdEntry> {
    gd.entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Jc)
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
            molecule_id: 0,
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
            molecule_id: 0,
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
fn same_reads_placed_on_two_copies_fold_to_one_jc() {
    // Multi-copy placement folding (job 2372015). Regression base for
    // the synth_is_mob empty-GD failures (jobs 2371258, 2371412): plus
    // and minus reads with 3′ softclips cross the SAME insert left
    // boundary. SAM/BAM stores SEQ reference-forward, so both strands
    // share the same stored SEQ/CIGAR and must key side1 at the
    // junction base (30). The SAME reads' clips hit both motif copies,
    // so after folding exactly ONE JC remains (not one per copy); with
    // no copy-anchored reads the two placements tie on support and the
    // deterministic tie-break keeps the smallest side2 (first copy).
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
                molecule_id: 0,
            });
        }
    }
    let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
    let jcs = gd_jcs(&gd);
    assert_eq!(
        jcs.len(),
        1,
        "same reads placed on two copies must fold to one JC; entries={:?}",
        gd.entries
            .iter()
            .map(|e| (e.kind, e.fields.clone()))
            .collect::<Vec<_>>()
    );
    let side1: u64 = jcs[0].fields[1].parse().unwrap();
    let side2: u64 = jcs[0].fields[4].parse().unwrap();
    assert_eq!(
        side1, 30,
        "both strands must key side1 at the junction base"
    );
    assert_eq!(
        side2, 41,
        "support tie keeps the lexicographically smallest side2"
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
            molecule_id: 0,
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
            molecule_id: 0,
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
                    molecule_id: 0,
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

/// 80 bp IS-like motif with 2 bp register coincidences against the
/// flank markers used in the folding tests: the first 2 bp equal the
/// right-flank start ("TT") and the last 2 bp equal the left-flank end
/// ("AA"), mirroring the measured synth_is_mob motif/ref coincidences
/// (job 2371412 diagnostic: motif prefix GT == ref 601-602, motif
/// suffix GC == ref 599-600).
fn is_motif_80() -> Vec<u8> {
    let mut m = Vec::new();
    m.extend_from_slice(b"TTGCAAGTCGATCGTTAGCC"); // [0..20]
    m.extend_from_slice(b"GATCCGATGCTAGCTTACGA"); // [20..40]
    m.extend_from_slice(b"ACGTGGATCCTTAGCGATCA"); // [40..60]
    m.extend_from_slice(b"GCTAGGCCATGATCGATAAA"); // [60..80]
    m
}

fn softclip_read(
    ref_start_0: i64,
    minus: bool,
    seq: &[u8],
    cigar: &[(CigarKind, usize)],
) -> AlignedRead {
    AlignedRead {
        contig_idx: 0,
        ref_start_0,
        minus,
        seq: seq.to_vec(),
        cigar: cigar
            .iter()
            .map(|&(kind, len)| CigarOp { kind, len })
            .collect(),
        mapq: UNIQUE_MAPQ,
        molecule_id: 0,
    }
}

#[test]
fn synth_is_mob_2372015_multicopy_and_reverse_duplicates_fold() {
    // End-to-end geometry of job 2372015 (4000 bp genome, 80 bp motif
    // copies at 1201-1280 / 2481-2560, third copy inserted at 601).
    // Each insert boundary is reported (a) once per existing motif copy
    // (multi-copy placement of the same clipped reads) and (b) once in
    // the reverse direction by reads aligned on a copy with the unique
    // flank clipped. Before folding this emitted 6 JC for 2 physical
    // junctions (4 redundant over-JC vs breseq: 599->1278, 602->1203,
    // 2480->599, 2560->601); after folding exactly 2 canonical JC must
    // remain.
    let motif = is_motif_80();
    // Unique flank markers (NOT homopolymers: copy-anchored reads clip
    // the flank, and place_softclip searches both strands — an A/T
    // homopolymer pair would reverse-complement onto each other).
    // lmark ends in "AA" and rmark starts with "TT", mirroring the
    // measured 2 bp register coincidences with the motif ends.
    let lmark = *b"ACGATCGGTTAACTGCGGAA"; // genome 581..600 (1-based)
    let rmark = *b"TTGCGATCGATCCGGAATCG"; // genome 601..620 (1-based)
    let mut seq = vec![b'C'; 4000];
    seq[580..600].copy_from_slice(&lmark);
    seq[600..620].copy_from_slice(&rmark);
    seq[1200..1280].copy_from_slice(&motif); // copy1, 1-based 1201..1280
    seq[2480..2560].copy_from_slice(&motif); // copy2, 1-based 2481..2560

    // Stored SEQ (reference-forward on both strands).
    let mut left_junction_seq = lmark.to_vec();
    left_junction_seq.extend_from_slice(&motif[0..20]);
    let mut right_junction_seq = Vec::new();
    right_junction_seq.extend_from_slice(&motif[60..80]);
    right_junction_seq.extend_from_slice(&rmark);

    let mut reads = Vec::new();
    // Left junction, flank-anchored: M on the left flank, trailing clip
    // = motif prefix. Register 600 (20M+20S) and register 602 (22M+18S;
    // the M block extends 2 bp into the motif because motif[0..2]=="TT"
    // equals the right-flank start).
    for minus in [false, false, true, true] {
        reads.push(softclip_read(
            580,
            minus,
            &left_junction_seq,
            &[(CigarKind::Match, 20), (CigarKind::SoftClip, 20)],
        ));
    }
    for minus in [false, true] {
        reads.push(softclip_read(
            580,
            minus,
            &left_junction_seq,
            &[(CigarKind::Match, 22), (CigarKind::SoftClip, 18)],
        ));
    }
    // Left junction, copy-anchored (bowtie2 placed the motif block on
    // copy2): leading clip = left-flank suffix.
    for minus in [false, true] {
        reads.push(softclip_read(
            2480,
            minus,
            &left_junction_seq,
            &[(CigarKind::SoftClip, 20), (CigarKind::Match, 20)],
        ));
    }
    // Right junction, flank-anchored: leading clip = motif suffix, M on
    // the right flank. Register 601 (20S+20M) and register 599
    // (18S+22M; motif[78..80]=="AA" equals the left-flank end).
    for minus in [false, false, true, true] {
        reads.push(softclip_read(
            600,
            minus,
            &right_junction_seq,
            &[(CigarKind::SoftClip, 20), (CigarKind::Match, 20)],
        ));
    }
    for minus in [false, true] {
        reads.push(softclip_read(
            598,
            minus,
            &right_junction_seq,
            &[(CigarKind::SoftClip, 18), (CigarKind::Match, 22)],
        ));
    }
    // Right junction, copy-anchored: M on copy2 end, trailing clip =
    // right-flank prefix.
    for minus in [false, true] {
        reads.push(softclip_read(
            2540,
            minus,
            &right_junction_seq,
            &[(CigarKind::Match, 20), (CigarKind::SoftClip, 20)],
        ));
    }

    let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
    let jcs = gd_jcs(&gd);
    assert_eq!(
        jcs.len(),
        2,
        "6 redundant reports must fold to exactly 2 canonical JC; entries={:?}",
        gd.entries
            .iter()
            .map(|e| (e.kind, e.fields.clone()))
            .collect::<Vec<_>>()
    );
    // Canonical coordinates: the copy2 placements carry the extra
    // copy-anchored reads, so they win the multi-copy fold.
    assert_eq!(
        jcs[0].fields,
        vec!["chr", "599", "1", "chr", "2558", "-1", "0"],
        "left junction canonical JC"
    );
    assert_eq!(
        jcs[1].fields,
        vec!["chr", "600", "-1", "chr", "2481", "1", "0"],
        "right junction canonical JC"
    );
}

#[test]
fn distinct_junctions_apart_beyond_tol_are_not_folded() {
    // Two genuinely distinct junctions 40 bp apart (> tol), each
    // supported by its own read set whose clips place on both motif
    // copies. Each junction's two copy placements fold (same reads);
    // the two junctions must NOT fold into each other.
    let motif = is_motif_80();
    let mut seq = vec![b'C'; 4000];
    seq[580..600].fill(b'A');
    seq[620..640].fill(b'G');
    seq[1200..1280].copy_from_slice(&motif);
    seq[2480..2560].copy_from_slice(&motif);
    let mut reads = Vec::new();
    for (flank_base, start0) in [(b'A', 580i64), (b'G', 620)] {
        let mut read_seq = vec![flank_base; 20];
        read_seq.extend_from_slice(&motif[0..20]);
        for minus in [false, false, false, true, true, true] {
            reads.push(softclip_read(
                start0,
                minus,
                &read_seq,
                &[(CigarKind::Match, 20), (CigarKind::SoftClip, 20)],
            ));
        }
    }
    let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
    let jcs = gd_jcs(&gd);
    assert_eq!(
        jcs.len(),
        2,
        "junctions 40 bp apart must not fold; entries={:?}",
        gd.entries
            .iter()
            .map(|e| (e.kind, e.fields.clone()))
            .collect::<Vec<_>>()
    );
    let mut side1s: Vec<u64> = jcs.iter().map(|jc| jc.fields[1].parse().unwrap()).collect();
    side1s.sort_unstable();
    assert_eq!(side1s, vec![600, 640]);
    for jc in &jcs {
        // Equal support on both copies -> smallest side2 (copy1).
        assert_eq!(jc.fields[4], "1201");
    }
}

#[test]
fn distinct_read_sets_same_flank_are_not_folded() {
    // Same side1 (within tol) but DISJOINT supporting read sets placing
    // on two different copies: not the same underlying reads, so the
    // multi-copy rule must not fold them (identity gate).
    let motif_p = is_motif_80();
    let mut motif_q = Vec::new();
    motif_q.extend_from_slice(b"GCGTTAACCGGATCGTACGT");
    motif_q.extend_from_slice(b"CCGGAATTAACCGGTTAACC");
    motif_q.extend_from_slice(b"TTGGCCAATTGGCCAATTGG");
    motif_q.extend_from_slice(b"AACCGGTTCCTTAAGGCCAA");
    let mut seq = vec![b'C'; 4000];
    seq[580..600].fill(b'A');
    seq[1200..1280].copy_from_slice(&motif_p);
    seq[2480..2560].copy_from_slice(&motif_q);
    let mut reads = Vec::new();
    for motif in [&motif_p, &motif_q] {
        let mut read_seq = vec![b'A'; 20];
        read_seq.extend_from_slice(&motif[0..20]);
        for minus in [false, false, false, true, true, true] {
            reads.push(softclip_read(
                580,
                minus,
                &read_seq,
                &[(CigarKind::Match, 20), (CigarKind::SoftClip, 20)],
            ));
        }
    }
    let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
    let jcs = gd_jcs(&gd);
    assert_eq!(
        jcs.len(),
        2,
        "disjoint read sets must not fold even at the same side1; entries={:?}",
        gd.entries
            .iter()
            .map(|e| (e.kind, e.fields.clone()))
            .collect::<Vec<_>>()
    );
    let mut side2s: Vec<u64> = jcs.iter().map(|jc| jc.fields[4].parse().unwrap()).collect();
    side2s.sort_unstable();
    assert_eq!(side2s, vec![1201, 2481]);
}

#[test]
fn cigar_split_junctions_are_never_multicopy_folded() {
    // Two long-D deletions sharing side1 (within tol) but reaching
    // different side2 positions: distinct events. CIGAR splits have
    // unique identity and must survive as 2 JC.
    let seq = vec![b'A'; 400];
    let mut reads = Vec::new();
    for del_len in [8usize, 100] {
        for minus in [false, false, false, true, true, true] {
            reads.push(softclip_read(
                0,
                minus,
                &[b'A'; 30],
                &[
                    (CigarKind::Match, 15),
                    (CigarKind::Del, del_len),
                    (CigarKind::Match, 15),
                ],
            ));
        }
    }
    let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
    let jcs = gd_jcs(&gd);
    assert_eq!(
        jcs.len(),
        2,
        "CIGAR-split junctions must never be multi-copy folded; entries={:?}",
        gd.entries
            .iter()
            .map(|e| (e.kind, e.fields.clone()))
            .collect::<Vec<_>>()
    );
    let mut side2s: Vec<u64> = jcs.iter().map(|jc| jc.fields[4].parse().unwrap()).collect();
    side2s.sort_unstable();
    assert_eq!(side2s, vec![24, 116]);
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
            molecule_id: 0,
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

fn short_clip_is_seed_reads() -> (Vec<u8>, Vec<AlignedRead>) {
    let mut seq = vec![b'C'; 160];
    let motif = *b"ACGTACGTAC";
    seq[40..50].copy_from_slice(&motif);
    seq[100..110].copy_from_slice(&motif);
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
                    len: 10,
                },
            ],
            mapq: UNIQUE_MAPQ,
            molecule_id: 0,
        });
    }
    (seq, reads)
}

#[test]
fn ten_bp_clips_seed_but_do_not_accept_without_second_pass() {
    let (seq, reads) = short_clip_is_seed_reads();
    let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
    assert!(
        gd_jcs(&gd).is_empty(),
        "10 bp clips are below the 14 bp accept floor; got {:?}",
        gd_jcs(&gd).iter().map(|e| &e.fields).collect::<Vec<_>>()
    );
}

#[test]
fn second_pass_hits_promote_short_clip_seed_to_jc() {
    let (seq, reads) = short_clip_is_seed_reads();
    // Seed keys side1 at last match (90). Extra hits model spanning
    // second-pass alignments (18 bp into each side, both strands).
    let extra: Vec<SplitCandidate> = [false, false, false, true, true, true]
        .into_iter()
        .map(|minus| SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: true,
            side2_minus: false,
            left: SubAlignment {
                read_start: 0,
                read_end: 18,
            },
            right: SubAlignment {
                read_start: 18,
                read_end: 36,
            },
            side1_pos_1: 90,
            side2_pos_1: 41,
            overlap: 0,
            origin: SplitOrigin::Softclip(10_000),
            seed_reads: 0,
        })
        .collect();
    let gd = call_from_aligned_extra(&fasta_chr(&seq), &reads, &opts_single_thread(), &extra);
    assert!(
        !gd_jcs(&gd).is_empty(),
        "second-pass spanning hits must lift the seed past accept_junction; entries={:?}",
        gd.entries
            .iter()
            .map(|e| (e.kind, e.fields.clone()))
            .collect::<Vec<_>>()
    );
}

// --- candidate-junction second pass (`parse_junction_bam`) ---------------
//
// Fixture BAMs are authored as SAM text and converted with the production
// `sam_to_sorted_bam`, so these cover the real decode path.

fn write_fixture_bam(dir: &str, sq: &[(&str, usize)], recs: &[String]) -> PathBuf {
    let d = std::env::temp_dir().join(dir);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    let mut text = String::from("@HD\tVN:1.6\n");
    for (name, len) in sq {
        text.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    for r in recs {
        text.push_str(r);
        text.push('\n');
    }
    let sam = d.join("in.sam");
    std::fs::write(&sam, text).unwrap();
    let bam = d.join("in.bam");
    crate::align::sam_to_sorted_bam(&sam, &bam).unwrap();
    bam
}

fn sam_rec(qname: &str, flag: u16, rname: &str, pos: usize, cigar: &str, qlen: usize) -> String {
    format!(
        "{qname}\t{flag}\t{rname}\t{pos}\t44\t{cigar}\t*\t0\t0\t{}\t{}",
        "A".repeat(qlen),
        "I".repeat(qlen)
    )
}

/// bowtie2-style record carrying its own alignment score. The junction
/// second pass runs with `-k 2000`, so one molecule gets one record per
/// competing construct, each with the aligner's `AS:i:` for THAT hit.
fn sam_rec_as(
    qname: &str,
    flag: u16,
    rname: &str,
    pos: usize,
    cigar: &str,
    qlen: usize,
    as_score: i32,
) -> String {
    format!(
        "{}\tAS:i:{as_score}",
        sam_rec(qname, flag, rname, pos, cigar, qlen)
    )
}

/// Candidate paired with the breakpoint of ITS construct.
fn jc_cand(breakpoint: usize) -> (crate::jc_seq::CandidateJunction, usize) {
    (
        crate::jc_seq::CandidateJunction {
            side1_contig: 0,
            side1_pos_1: 90,
            side1_minus: true,
            side2_contig: 0,
            side2_pos_1: 41,
            side2_minus: false,
            overlap: 0,
            seed_reads: 3,
        },
        breakpoint,
    )
}

#[test]
fn second_pass_uses_per_construct_breakpoint_not_flank() {
    // Construct seeded at a contig start: left flank clamped to 10 bp, so
    // the breakpoint is 10 even though `flank` was 36. A 30M hit at
    // construct pos 1 puts 10 bp left and 20 bp right of the true
    // breakpoint — genuine spanning evidence.
    let bam = write_fixture_bam(
        "prokdiff-jcbam-breakpoint",
        &[("jc0", 30)],
        &[sam_rec("r1", 0, "jc0", 1, "30M", 30)],
    );
    let scores = PrimaryScores::default();

    let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(10)], &scores, &[], &[]).unwrap();
    assert_eq!(stats.considered, 1);
    assert_eq!(stats.kept, 1);
    assert_eq!(extra.len(), 1);
    assert_eq!(
        extra[0].left.len(),
        10,
        "side 1 extension is breakpoint-cstart"
    );
    assert_eq!(extra[0].right.len(), 20);

    // Same record against the pre-fix constant breakpoint (`flank` = 36):
    // the whole hit looks left-of-breakpoint and the read is discarded.
    // This is the JC=0 mechanism for every construct near a contig edge.
    let (extra_flank, stats_flank) =
        parse_junction_bam(&bam, &[jc_cand(36)], &scores, &[], &[]).unwrap();
    assert!(extra_flank.is_empty());
    assert_eq!(stats_flank.not_spanning, 1);
}

#[test]
fn second_pass_rejects_two_bp_overhang_instead_of_clamping_it_up() {
    // 36M at construct pos 35 lands 2 bp past a breakpoint of 36. The
    // engine used to admit this and clamp both extensions to >=3, which
    // fabricated evidence and made `accept_junction` rule 4 unreachable.
    let bam = write_fixture_bam(
        "prokdiff-jcbam-overhang",
        &[("jc0", 72)],
        &[sam_rec("r1", 0, "jc0", 35, "36M", 36)],
    );
    let (extra, stats) =
        parse_junction_bam(&bam, &[jc_cand(36)], &PrimaryScores::default(), &[], &[]).unwrap();
    assert!(extra.is_empty(), "2 bp past the breakpoint is not evidence");
    assert_eq!(stats.not_spanning, 1);
    assert_eq!(stats.kept, 0);
}

#[test]
fn second_pass_reports_real_overlaps_with_no_floor() {
    // Exactly 3 bp on the short side: admitted, and reported as 3 — not
    // rounded, not clamped.
    let bam = write_fixture_bam(
        "prokdiff-jcbam-exact3",
        &[("jc0", 72)],
        &[sam_rec("r1", 0, "jc0", 34, "36M", 36)],
    );
    let (extra, stats) =
        parse_junction_bam(&bam, &[jc_cand(36)], &PrimaryScores::default(), &[], &[]).unwrap();
    assert_eq!(stats.kept, 1);
    assert_eq!(extra[0].left.len(), 3);
    assert_eq!(extra[0].right.len(), 33);
}

#[test]
fn second_pass_counts_short_cover_rejects() {
    // 20 aligned query bases is below JC_MIN_COVER_BASES (28).
    let bam = write_fixture_bam(
        "prokdiff-jcbam-shortcover",
        &[("jc0", 72)],
        &[sam_rec("r1", 0, "jc0", 25, "20M16S", 36)],
    );
    let (extra, stats) =
        parse_junction_bam(&bam, &[jc_cand(36)], &PrimaryScores::default(), &[], &[]).unwrap();
    assert!(extra.is_empty());
    assert_eq!(stats.short_cover, 1);
    assert_eq!(stats.not_spanning, 0, "cover is checked before spanning");
}

#[test]
fn second_pass_counts_primary_score_gate_separately_from_spanning() {
    // A spanning 36M hit scores 36. The gate keeps it against a weaker
    // primary alignment and drops it against a stronger one; the counter
    // distinguishes the two so a Clonal job can tell which gate ate the
    // reads without re-running.
    let bam = write_fixture_bam(
        "prokdiff-jcbam-scoregate",
        &[("jc0", 72)],
        &[sam_rec("r1", 0, "jc0", 19, "36M", 36)],
    );
    let cands = [jc_cand(36)];

    let weaker = PrimaryScores::from_single_ended([("r1", 30)]);
    let (extra, stats) = parse_junction_bam(&bam, &cands, &weaker, &[], &[]).unwrap();
    assert_eq!(stats.kept, 1);
    assert_eq!(stats.worse_than_primary, 0);
    assert_eq!(extra.len(), 1);

    let stronger = PrimaryScores::from_single_ended([("r1", 40)]);
    let (extra, stats) = parse_junction_bam(&bam, &cands, &stronger, &[], &[]).unwrap();
    assert!(extra.is_empty());
    assert_eq!(stats.worse_than_primary, 1);
    assert_eq!(stats.not_spanning, 0);
}

#[test]
fn pe_mate_scores_are_isolated_so_perfect_mate_does_not_gate_junction_mate() {
    // In paired-end sequencing, mate 1 may have aligned perfectly to the
    // reference genome (score 150), while mate 2 crosses a structural junction
    // and had a weak primary alignment (score 20).
    //
    // The second-pass junction alignment for mate 2 scores 36 on the junction
    // construct. With per-mate score isolation, mate 2 compares only against its
    // own primary score (20), so 36 >= 20 keeps the junction read.
    // Without mate isolation, mate 1's score of 150 would erroneously gate out mate 2.
    let bam = write_fixture_bam(
        "prokdiff-jcbam-pe-mate-isolation",
        &[("jc0", 72)],
        &[sam_rec("read_pe", 129, "jc0", 19, "36M", 36)], // flag 129: paired + mate 2
    );
    let cands = [jc_cand(36)];

    let scores = PrimaryScores::from_iter([
        (("read_pe", Mate::First), 150),
        (("read_pe", Mate::Last), 20),
    ]);

    let (extra, stats) = parse_junction_bam(&bam, &cands, &scores, &[], &[]).unwrap();
    assert_eq!(
        stats.kept, 1,
        "mate 2 is accepted against its own primary score"
    );
    assert_eq!(stats.worse_than_primary, 0);
    assert_eq!(extra.len(), 1);

    // If mate 2's own primary score had been superior (e.g. 40), it would be rejected.
    let scores_higher = PrimaryScores::from_iter([
        (("read_pe", Mate::First), 150),
        (("read_pe", Mate::Last), 40),
    ]);
    let (extra_high, stats_high) =
        parse_junction_bam(&bam, &cands, &scores_higher, &[], &[]).unwrap();
    assert!(extra_high.is_empty());
    assert_eq!(stats_high.worse_than_primary, 1);
}

#[test]
fn second_pass_origins_are_per_read_and_never_alias_softclip_ids() {
    let bam = write_fixture_bam(
        "prokdiff-jcbam-origins",
        &[("jc0", 72)],
        &[
            sam_rec("r1", 0, "jc0", 19, "36M", 36),
            sam_rec("r2", 16, "jc0", 20, "36M", 36),
        ],
    );
    let (extra, stats) =
        parse_junction_bam(&bam, &[jc_cand(36)], &PrimaryScores::default(), &[], &[]).unwrap();
    assert_eq!(stats.kept, 2);
    // Two reads, two identities: the multi-copy fold must be able to tell
    // them apart.
    let ids: HashSet<SplitOrigin> = extra.iter().map(|c| c.origin).collect();
    assert_eq!(ids.len(), 2);
    // Variant-tagged, so a qname hash can never collide with a first-pass
    // `SoftclipHint` index in the supporting-read identity set.
    assert!(extra
        .iter()
        .all(|c| matches!(c.origin, SplitOrigin::SecondPass(_))));
    // Read strand is carried through; side strands come from the candidate.
    assert!(extra.iter().any(|c| c.minus));
    assert!(extra.iter().any(|c| !c.minus));
}

#[test]
fn second_pass_ignores_records_on_unknown_constructs() {
    // A reference name that is not `jc<idx>`, and an index past the end of
    // the retained candidate list: neither may be counted or promoted.
    let bam = write_fixture_bam(
        "prokdiff-jcbam-unknown",
        &[("chr", 72), ("jc7", 72)],
        &[
            sam_rec("r1", 0, "chr", 19, "36M", 36),
            sam_rec("r2", 0, "jc7", 19, "36M", 36),
        ],
    );
    let (extra, stats) =
        parse_junction_bam(&bam, &[jc_cand(36)], &PrimaryScores::default(), &[], &[]).unwrap();
    assert!(extra.is_empty());
    assert_eq!(stats.considered, 0);
}

#[test]
fn is150_microhomology_sliding_recovers_mob_from_dup2() {
    let mut seq = vec![b'G'; 2500];
    // Target sequence around 1700:
    // 0-based 1700..1706 is 1-based 1701..1706.
    // 1-based:
    // 1701: 'A'
    // 1702: 'C'
    // 1703: 'A'
    // If p_plus = 1702, p_minus = 1703, dup = 2.
    // Base at 1701 ('A') == base at 1703 ('A'), so p_plus slides to 1701 and dup becomes 3.
    seq[1700..1706].copy_from_slice(b"ACAACA");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS150".into(),
        strand: 1,
        start: 500,
        end: 1500,
    }];
    let opts = EngineOptions {
        repeats,
        ..Default::default()
    };
    let mut splits = Vec::new();
    // Junction 1: IS150 end (1500, minus) to target (1703, minus)
    for (i, minus) in [false, false, false, true, true, true]
        .into_iter()
        .enumerate()
    {
        splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: true,
            side2_minus: true,
            left: SubAlignment {
                read_start: 0,
                read_end: 18,
            },
            right: SubAlignment {
                read_start: 18,
                read_end: 36,
            },
            side1_pos_1: 1500,
            side2_pos_1: 1703,
            overlap: 0,
            origin: SplitOrigin::Softclip(100 + i),
            seed_reads: 0,
        });
    }
    // Junction 2: IS150 start (500, plus) to target (1702, plus)
    for (i, minus) in [false, false, false, true, true, true]
        .into_iter()
        .enumerate()
    {
        splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: false,
            side2_minus: false,
            left: SubAlignment {
                read_start: 0,
                read_end: 18,
            },
            right: SubAlignment {
                read_start: 18,
                read_end: 36,
            },
            side1_pos_1: 500,
            side2_pos_1: 1702,
            overlap: 0,
            origin: SplitOrigin::Softclip(200 + i),
            seed_reads: 0,
        });
    }
    let contig_results = vec![ContigPileup {
        columns: fasta[0]
            .seq
            .iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: vec![],
                insertions: vec![],
                deletions: vec![],
            })
            .collect(),
        unique_depth: vec![10; 2500],
        total_depth: vec![10; 2500],
        splits,
    }];
    let gd = emit_from_pileup(&fasta, contig_results, &opts, &[]);
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(mobs.len(), 1, "Expected 1 MOB, got: {:?}", mobs);
    assert_eq!(mobs[0].fields[0], "chr");
    assert_eq!(mobs[0].fields[1], "1701");
    assert_eq!(mobs[0].fields[2], "IS150");
    assert_eq!(mobs[0].fields[4], "3");
    assert_eq!(
        mobs[0].parent_ids.len(),
        2,
        "MOB must link to constituent JCs: {:?}",
        mobs[0].parent_ids
    );
}

#[test]
fn test_depth_aware_jc_frequency_filter_rejects_noise_and_keeps_true_junction() {
    // Models RW-002: At 400x coverage, 4 reads (1% frequency) is chimeric noise and
    // must be filtered by jc_min_frequency = 0.05. A 40-read junction (10% frequency)
    // at the same depth passes. At 30x depth, a 4-read junction (13.3% frequency) also passes.
    let seq = vec![b'A'; 2000];
    let fasta = fasta_chr(&seq);
    let opts = EngineOptions {
        threads: 1,
        jc_min_support_reads: 3,
        jc_min_frequency: 0.05,
        ..EngineOptions::default()
    };

    // Case 1: High depth (400x) with 4 supporting reads (2 plus, 2 minus).
    // Clears candidate accept_junction (both strands, >=14bp overlap), but frequency = 4/400 = 1% < 5%.
    let mut noise_splits = Vec::new();
    for (i, &minus) in [false, false, true, true].iter().enumerate() {
        noise_splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: true,
            side2_minus: false,
            left: SubAlignment {
                read_start: 0,
                read_end: 20,
            },
            right: SubAlignment {
                read_start: 20,
                read_end: 40,
            },
            side1_pos_1: 500,
            side2_pos_1: 1500,
            overlap: 0,
            origin: SplitOrigin::Softclip(100 + i),
            seed_reads: 0,
        });
    }

    let contig_results_noise = vec![ContigPileup {
        columns: vec![
            PileupColumn {
                ref_base: b'A',
                observations: vec![],
                insertions: vec![],
                deletions: vec![],
            };
            2000
        ],
        unique_depth: vec![400; 2000],
        total_depth: vec![400; 2000],
        splits: noise_splits,
    }];
    let gd_noise = emit_from_pileup(&fasta, contig_results_noise, &opts, &[]);
    assert_eq!(
        gd_jcs(&gd_noise).len(),
        0,
        "4 reads at 400x depth (1% freq) must be filtered by 5% threshold"
    );

    // Case 2: High depth (400x) with 40 supporting reads (20 plus, 20 minus) = 10% freq.
    let mut true_splits = Vec::new();
    for i in 0..40 {
        true_splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus: i % 2 == 0,
            side1_minus: true,
            side2_minus: false,
            left: SubAlignment {
                read_start: 0,
                read_end: 20,
            },
            right: SubAlignment {
                read_start: 20,
                read_end: 40,
            },
            side1_pos_1: 500,
            side2_pos_1: 1500,
            overlap: 0,
            origin: SplitOrigin::Softclip(1000 + i),
            seed_reads: 0,
        });
    }

    let contig_results_true = vec![ContigPileup {
        columns: vec![
            PileupColumn {
                ref_base: b'A',
                observations: vec![],
                insertions: vec![],
                deletions: vec![],
            };
            2000
        ],
        unique_depth: vec![400; 2000],
        total_depth: vec![400; 2000],
        splits: true_splits,
    }];
    let gd_true = emit_from_pileup(&fasta, contig_results_true, &opts, &[]);
    assert_eq!(
        gd_jcs(&gd_true).len(),
        1,
        "40 reads at 400x depth (10% freq) must pass"
    );

    // Case 3: Moderate depth (30x) with 4 supporting reads (2 plus, 2 minus) = 13.3% freq.
    let mut low_depth_splits = Vec::new();
    for (i, &minus) in [false, false, true, true].iter().enumerate() {
        low_depth_splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: true,
            side2_minus: false,
            left: SubAlignment {
                read_start: 0,
                read_end: 20,
            },
            right: SubAlignment {
                read_start: 20,
                read_end: 40,
            },
            side1_pos_1: 500,
            side2_pos_1: 1500,
            overlap: 0,
            origin: SplitOrigin::Softclip(2000 + i),
            seed_reads: 0,
        });
    }

    let contig_results_low = vec![ContigPileup {
        columns: vec![
            PileupColumn {
                ref_base: b'A',
                observations: vec![],
                insertions: vec![],
                deletions: vec![],
            };
            2000
        ],
        unique_depth: vec![30; 2000],
        total_depth: vec![30; 2000],
        splits: low_depth_splits,
    }];
    let gd_low = emit_from_pileup(&fasta, contig_results_low, &opts, &[]);
    assert_eq!(
        gd_jcs(&gd_low).len(),
        1,
        "4 reads at 30x depth (13.3% freq >= 5%) must pass"
    );
}

#[test]
fn mob_lineage_unlinked() {
    // 1) Test with emit_from_pileup: MOB is promoted from constituent JCs and links to them.
    // An unrelated third JC (not flanking this MOB) is emitted independently and NOT in MOB parent_ids.
    let mut seq = vec![b'G'; 2500];
    seq[1700..1706].copy_from_slice(b"ACAACA");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS150".into(),
        strand: 1,
        start: 500,
        end: 1500,
    }];
    let opts = EngineOptions {
        repeats,
        jc_min_support_reads: 3,
        jc_min_frequency: 0.05,
        ..Default::default()
    };
    let mut splits = Vec::new();
    // Constituent JC 1: IS150 end (1500, minus) to target (1703, minus)
    for (i, minus) in [false, false, false, true, true, true]
        .into_iter()
        .enumerate()
    {
        splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: true,
            side2_minus: true,
            left: SubAlignment {
                read_start: 0,
                read_end: 18,
            },
            right: SubAlignment {
                read_start: 18,
                read_end: 36,
            },
            side1_pos_1: 1500,
            side2_pos_1: 1703,
            overlap: 0,
            origin: SplitOrigin::Softclip(100 + i),
            seed_reads: 0,
        });
    }
    // Constituent JC 2: IS150 start (500, plus) to target (1702, plus)
    for (i, minus) in [false, false, false, true, true, true]
        .into_iter()
        .enumerate()
    {
        splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: false,
            side2_minus: false,
            left: SubAlignment {
                read_start: 0,
                read_end: 18,
            },
            right: SubAlignment {
                read_start: 18,
                read_end: 36,
            },
            side1_pos_1: 500,
            side2_pos_1: 1702,
            overlap: 0,
            origin: SplitOrigin::Softclip(200 + i),
            seed_reads: 0,
        });
    }
    // Unrelated JC 3: at positions 100 and 300 (not near repeat)
    for (i, minus) in [false, false, false, true, true, true]
        .into_iter()
        .enumerate()
    {
        splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: true,
            side2_minus: false,
            left: SubAlignment {
                read_start: 0,
                read_end: 18,
            },
            right: SubAlignment {
                read_start: 18,
                read_end: 36,
            },
            side1_pos_1: 100,
            side2_pos_1: 300,
            overlap: 0,
            origin: SplitOrigin::Softclip(300 + i),
            seed_reads: 0,
        });
    }
    let contig_results = vec![ContigPileup {
        columns: fasta[0]
            .seq
            .iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: vec![],
                insertions: vec![],
                deletions: vec![],
            })
            .collect(),
        unique_depth: vec![10; 2500],
        total_depth: vec![10; 2500],
        splits,
    }];
    let gd = emit_from_pileup(&fasta, contig_results, &opts, &[]);
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(mobs.len(), 1);
    let mob = mobs[0];
    assert_eq!(
        mob.parent_ids.len(),
        2,
        "MOB must link to exactly 2 constituent JCs"
    );

    let jcs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Jc)
        .collect();
    assert_eq!(jcs.len(), 3, "Expected 2 constituent JCs + 1 unrelated JC");

    let unrelated_jc = jcs
        .iter()
        .find(|j| j.fields[1] == "100" || j.fields[4] == "100")
        .expect("unrelated JC exists");
    assert!(
        !mob.parent_ids.contains(&unrelated_jc.id),
        "MOB parent_ids must NOT contain unrelated JC"
    );
    assert_eq!(
        unrelated_jc.attrs.get("mob_evidence"),
        None,
        "Unrelated JC must not have mob_evidence=1"
    );

    for pid in &mob.parent_ids {
        let constituent = jcs.iter().find(|j| j.id == *pid).expect("parent JC exists");
        assert_eq!(
            constituent.attrs.get("mob_evidence"),
            Some(&"1".to_string())
        );
    }

    // 2) Verify unlinked/historical MOB without parent_ids remains unlinked:
    let unlinked_mob = prokadiff_gd::GdEntry::mob(10, "chr", 1701, "IS150", "1", 3);
    assert!(unlinked_mob.parent_ids.is_empty());
}

#[test]
fn ra_metrics_absent() {
    // 1) Test an entry not called by RA (e.g. an MC-promoted DEL).
    // Set up a genome with a true gap that promotes to DEL via MC.
    let seq = vec![b'A'; 80];
    let mut reads = covering_reads(&seq, 0, 10, 6, 6, UNIQUE_MAPQ);
    reads.extend(covering_reads(&seq, 70, 10, 6, 6, UNIQUE_MAPQ));
    let gd = call_from_aligned(&fasta_chr(&seq), &reads, &opts_single_thread());
    let dels: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Del)
        .collect();
    assert_eq!(dels.len(), 1);
    let del = dels[0];
    // MC promoted DEL has parent_ids containing the MC record ID:
    assert!(!del.parent_ids.is_empty());
    // But MUST NOT have call-time pd_ra_* metrics:
    assert_eq!(del.attrs.get("pd_ra_depth"), None);
    assert_eq!(del.attrs.get("pd_ra_support_reads"), None);
    assert_eq!(del.attrs.get("pd_ra_frequency"), None);

    // 2) An RA-called SNP DOES have pd_ra_* metrics:
    let fasta = [FastaRecord {
        name: "chr".into(),
        seq: b"ACGTACGTAC".to_vec(),
    }];
    let mut alt = fasta[0].seq.clone();
    alt[2] = b'T';
    let reads = covering_reads(&alt, 0, 10, 6, 6, UNIQUE_MAPQ);
    let gd_snp = call_from_aligned(&fasta, &reads, &opts_single_thread());
    let snps: Vec<_> = gd_snp
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Snp)
        .collect();
    assert_eq!(snps.len(), 1);
    let snp = snps[0];
    assert_eq!(snp.attrs.get("pd_ra_depth"), Some(&"12".to_string()));
    assert_eq!(
        snp.attrs.get("pd_ra_support_reads"),
        Some(&"12".to_string())
    );
    assert_eq!(snp.attrs.get("pd_ra_frequency"), Some(&"1".to_string()));
}

#[test]
fn gd_core_fields_unchanged() {
    // Assert that adding attributes and parent_ids does not alter
    // the core positional fields or subtract_key().
    let snp_raw = prokadiff_gd::GdEntry::snp(1, "chr", 100, "T");
    let mut snp_with_attrs = prokadiff_gd::GdEntry::snp(1, "chr", 100, "T");
    snp_with_attrs
        .attrs
        .insert("pd_ra_depth".into(), "30".into());
    snp_with_attrs
        .attrs
        .insert("pd_ra_support_reads".into(), "30".into());
    snp_with_attrs
        .attrs
        .insert("pd_ra_frequency".into(), "1".into());

    assert_eq!(snp_raw.fields, snp_with_attrs.fields);
    assert_eq!(snp_raw.subtract_key(), snp_with_attrs.subtract_key());

    let mob_raw = prokadiff_gd::GdEntry::mob(2, "chr", 1701, "IS150", "1", 3);
    let mut mob_with_parents = prokadiff_gd::GdEntry::mob(2, "chr", 1701, "IS150", "1", 3);
    mob_with_parents.parent_ids = vec![5, 6];

    assert_eq!(mob_raw.fields, mob_with_parents.fields);
    assert_eq!(mob_raw.subtract_key(), mob_with_parents.subtract_key());

    let jc_raw = prokadiff_gd::GdEntry::jc(3, "chr", 500, "1", "chr", 1702, "1", 0);
    let mut jc_with_attrs = prokadiff_gd::GdEntry::jc(3, "chr", 500, "1", "chr", 1702, "1", 0);
    jc_with_attrs
        .attrs
        .insert("mob_evidence".into(), "1".into());
    jc_with_attrs
        .attrs
        .insert("pd_support_reads".into(), "10".into());

    assert_eq!(jc_raw.fields, jc_with_attrs.fields);
    assert_eq!(jc_raw.subtract_key(), jc_with_attrs.subtract_key());
}

// ---------------------------------------------------------------------------
// G1 Clonal JC parity fixes (RC1/RC2/RC3/RC4) — synthetic geometry only,
// no sample-specific coordinates.
// ---------------------------------------------------------------------------

#[test]
fn second_pass_equal_score_is_rejected_for_unique_candidates() {
    // Published second-pass rule: a junction alignment must score STRICTLY
    // higher than the read's best primary reference alignment. A read the
    // reference already fully explains (36 = 36) carries no junction
    // information; admitting it pads repeat-side constructs with
    // reference-explained reads, which both inflates spurious IS-end hub
    // junctions and ties the multi-copy fold (RC1/RC2).
    let bam = write_fixture_bam(
        "prokdiff-jcbam-equal-score",
        &[("jc0", 72)],
        &[sam_rec("r1", 0, "jc0", 19, "36M", 36)],
    );
    let equal = PrimaryScores::from_single_ended([("r1", 36)]);

    let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &equal, &[], &[]).unwrap();
    assert!(
        extra.is_empty(),
        "equal score must not be junction evidence"
    );
    assert_eq!(stats.worse_than_primary, 1);
}

#[test]
fn second_pass_keeps_equal_score_for_repeat_repeat_candidates() {
    // Repeat↔repeat junction reads can at best TIE the reference (both
    // flanks are reference sequence — the IS-mediated rearrangements
    // breseq emits as prediction=unknown / frequency=NA). When BOTH sides
    // of the candidate fall inside annotated repeats, an equal score is
    // kept so those junctions stay recoverable (RC3).
    let bam = write_fixture_bam(
        "prokdiff-jcbam-equal-score-repeat",
        &[("jc0", 72)],
        &[sam_rec("r1", 0, "jc0", 19, "36M", 36)],
    );
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 20,
            end: 100,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 30,
            end: 120,
        },
    ];
    let names = vec!["chr".to_string()];
    let equal = PrimaryScores::from_single_ended([("r1", 36)]);

    let (extra, stats) =
        parse_junction_bam(&bam, &[jc_cand(36)], &equal, &repeats, &names).unwrap();
    assert_eq!(stats.kept, 1);
    assert_eq!(extra.len(), 1);
    // Stage-2 seed strength travels on the promoted split candidate.
    assert_eq!(extra[0].seed_reads, 3);
}

/// A candidate junction that differs from `jc_cand` so two competing
/// constructs in one fixture are distinguishable by geometry.
fn jc_cand_other(breakpoint: usize) -> (crate::jc_seq::CandidateJunction, usize) {
    (
        crate::jc_seq::CandidateJunction {
            side1_contig: 0,
            side1_pos_1: 9_000,
            side1_minus: false,
            side2_contig: 0,
            side2_pos_1: 8_000,
            side2_minus: true,
            overlap: 0,
            seed_reads: 2,
        },
        breakpoint,
    )
}

#[test]
fn second_pass_rejects_suboptimal_competing_junction_alignment() {
    // Wave 11: paralogous constructs differ by 1–2 bp, so the SAME
    // molecule aligns to both — best on the true paralog (AS 36), worse
    // on the cross-mapping target (AS 26). bowtie2 reports one record per
    // construct (`-k 2000`); the sub-optimal record must not be credited
    // to the paralog it lost to, or it fabricates support there (the
    // 53-read `jc194` cross-map that built the FP MOB).
    let bam = write_fixture_bam(
        "prokdiff-jcbam-competing",
        &[("jc0", 72), ("jc1", 72)],
        &[
            sam_rec_as("r1", 0, "jc0", 19, "36M", 36, 36),
            sam_rec_as("r1", 0, "jc1", 19, "36M", 36, 26),
            sam_rec_as("r2", 16, "jc0", 20, "36M", 36, 35),
            sam_rec_as("r2", 16, "jc1", 20, "36M", 36, 27),
        ],
    );
    let cands = [jc_cand(36), jc_cand_other(36)];

    let (extra, stats) =
        parse_junction_bam(&bam, &cands, &PrimaryScores::default(), &[], &[]).unwrap();
    assert_eq!(stats.considered, 4);
    assert_eq!(
        stats.kept, 2,
        "only the best-scoring construct per molecule is credited"
    );
    assert_eq!(stats.worse_than_competing, 2);
    assert_eq!(extra.len(), 2);
    // Every credited alignment belongs to construct A (jc0).
    assert!(extra
        .iter()
        .all(|c| c.side1_pos_1 == 90 && c.side2_pos_1 == 41));
}

#[test]
fn second_pass_preserves_tied_scores_across_identical_repeat_copies() {
    // Invariance (Codex Sol Recommendation 2): identical repeat copies
    // give a molecule the SAME score on every construct. The competition
    // filter must keep ALL tied best placements — destroying the tie
    // would delete genuine multi-copy evidence and turn an unresolved
    // ambiguity into a random copy assignment.
    let bam = write_fixture_bam(
        "prokdiff-jcbam-competing-tie",
        &[("jc0", 72), ("jc1", 72)],
        &[
            sam_rec_as("r1", 0, "jc0", 19, "36M", 36, 36),
            sam_rec_as("r1", 0, "jc1", 19, "36M", 36, 36),
            sam_rec_as("r2", 16, "jc0", 20, "36M", 36, 36),
            sam_rec_as("r2", 16, "jc1", 20, "36M", 36, 36),
        ],
    );
    // Both constructs are repeat↔repeat, so the equal primary score is
    // also kept (RC3) before the competition stage runs.
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 20,
            end: 100,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 30,
            end: 120,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 7_950,
            end: 8_050,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 8_950,
            end: 9_050,
        },
    ];
    let names = vec!["chr".to_string()];
    let equal = PrimaryScores::from_single_ended([("r1", 36), ("r2", 36)]);
    let cands = [jc_cand(36), jc_cand_other(36)];

    let (extra, stats) = parse_junction_bam(&bam, &cands, &equal, &repeats, &names).unwrap();
    assert_eq!(stats.worse_than_competing, 0, "ties are not losers");
    assert_eq!(stats.kept, 4, "both repeat copies keep both molecules");
    assert_eq!(extra.len(), 4);
    let on_a = extra.iter().filter(|c| c.side1_pos_1 == 90).count();
    let on_b = extra.iter().filter(|c| c.side1_pos_1 == 9_000).count();
    assert_eq!((on_a, on_b), (2, 2));
}

#[test]
fn multicopy_fold_support_tie_prefers_nearer_copy_over_higher_seed() {
    // Identical termini share spanning reads, so a higher Stage-2 seed on
    // the farther copy is not copy-discriminating. On an own-support tie
    // the displayed copy is the one closer to the unique flank. A support
    // gap above the tie window still wins, and that case is
    // `support_gap_above_tie_window_keeps_stronger_copy`.
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'A'; 30_000],
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 2_443,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 20_000,
            end: 21_443,
        },
    ];
    let mk = |rep_pos: u64, seed: usize, id: usize| {
        super::jc_cluster::AcceptedJc::from_parts(
            0,
            5_000,
            true,
            0,
            rep_pos,
            false,
            0,
            crate::jc::JunctionSupport {
                plus_reads: 5,
                minus_reads: 5,
                best_min_overlap: 18,
                plus_best_min: 18,
                minus_best_min: 16,
                min_overlap_side1: 16,
                min_overlap_side2: 15,
            },
            HashSet::from([(0, SplitOrigin::Softclip(id))]),
            false,
            seed,
        )
    };

    // Shared reads make the two termini one junction. Seed 20 on the far
    // copy must not hide the nearer end. Family seed stays the max.
    let folded = super::jc_cluster::fold_multicopy_placements(
        vec![mk(1_000, 5, 1), mk(20_000, 20, 1)],
        &repeats,
        &fasta,
    );
    assert_eq!(
        folded.len(),
        1,
        "shared-evidence copies must fold to one JC"
    );
    assert_eq!(
        folded[0].p2, 1_000,
        "equal support must display the copy nearer the unique flank"
    );
    assert_eq!(folded[0].own_seed_reads, 5);
    assert_eq!(folded[0].seed_reads, 20);
    assert_eq!(
        folded[0].support.plus_reads + folded[0].support.minus_reads,
        10,
        "shared observations must not be summed into the family support"
    );

    // Same repeat name, no shared reads, and different unique flanks:
    // still two junctions. The family name is not equivalence.
    let other_flank = super::jc_cluster::AcceptedJc::from_parts(
        0,
        25_000,
        false,
        0,
        20_000,
        false,
        0,
        crate::jc::JunctionSupport {
            plus_reads: 5,
            minus_reads: 5,
            best_min_overlap: 18,
            plus_best_min: 18,
            minus_best_min: 16,
            min_overlap_side1: 16,
            min_overlap_side2: 15,
        },
        HashSet::from([(0, SplitOrigin::Softclip(2))]),
        false,
        5,
    );
    let named_only = super::jc_cluster::fold_multicopy_placements(
        vec![mk(1_000, 20, 1), other_flank],
        &repeats,
        &fasta,
    );
    assert_eq!(named_only.len(), 2);

    // Equal seeds, shared reads: the copy closer to the unique flank wins.
    let tied = super::jc_cluster::fold_multicopy_placements(
        vec![mk(1_000, 7, 3), mk(20_000, 7, 3)],
        &repeats,
        &fasta,
    );
    assert_eq!(tied.len(), 1);
    assert_eq!(tied[0].p2, 1_000);
}

#[test]
fn normalize_snaps_repeat_side_to_terminus_and_absorbs_overlap() {
    // RC4: repeat-side registers that over-count microhomology past a copy
    // terminus slide back onto the terminus with the overlap fully
    // absorbed; near-terminus sides snap; interior and unique sides are
    // untouched.
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS150".into(),
        strand: 1,
        start: 1_000,
        end: 2_443,
    }];
    let side = |pos: u64, minus: bool| super::jc_cluster::JcSideRef {
        pos,
        minus,
        contig: "chr",
    };
    let norm = |p1: u64, m1: bool, p2: u64, m2: bool, ov: i64| {
        super::jc_cluster::normalize_jc_coords(side(p1, m1), side(p2, m2), ov, &repeats)
    };

    // 9 bp past the copy end with 9 bp counted microhomology: slide −9
    // lands exactly on the terminus → absorbed, overlap 0.
    assert_eq!(norm(5_000, false, 2_452, true, 9), (5_000, 2_443, 0));
    // Same geometry on the side-1 repeat flank: a plus-strand side slides
    // +k, so 9 bp before the copy START absorbs onto the start.
    assert_eq!(norm(991, false, 9_000, true, 9), (1_000, 9_000, 0));
    // Terminus snap without overlap: 2 bp inside the end lands on it.
    assert_eq!(norm(5_000, false, 2_441, true, 0), (5_000, 2_443, 0));
    // Interior repeat side: the slide does not land on a terminus, so the
    // overlap and coordinates are preserved.
    assert_eq!(norm(5_000, false, 2_000, true, 9), (5_000, 2_000, 9));
    // Unique sides: no repeat normalization applies.
    assert_eq!(norm(5_000, false, 9_000, true, 9), (5_000, 9_000, 9));
}

#[test]
fn interrepeat_junction_gate_requires_copy_anchoring() {
    // RC3: repeat↔repeat junctions are no longer blanket-skipped.
    //  - both sides at termini of DIFFERENT copies → emitted (element-
    //    boundary event; previously skipped outright);
    //  - one side at the terminus of the SAME copy containing the other
    //    side → emitted with pd_repeat_junction=1;
    //  - two INTERIOR points of the SAME copy → Stage-2 intra-copy
    //    pairing noise → dropped.
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'A'; 30_000],
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 2_443,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 20_000,
            end: 21_443,
        },
    ];
    let opts = EngineOptions {
        repeats,
        ..EngineOptions::default()
    };
    let mut splits = Vec::new();
    let push =
        |splits: &mut Vec<SplitCandidate>, base: usize, p1: u64, m1: bool, p2: u64, m2: bool| {
            for (i, minus) in [false, false, false, true, true, true]
                .into_iter()
                .enumerate()
            {
                splits.push(SplitCandidate {
                    contig_idx: 0,
                    side2_contig_idx: 0,
                    minus,
                    side1_minus: m1,
                    side2_minus: m2,
                    left: SubAlignment {
                        read_start: 0,
                        read_end: 18,
                    },
                    right: SubAlignment {
                        read_start: 18,
                        read_end: 36,
                    },
                    side1_pos_1: p1,
                    side2_pos_1: p2,
                    overlap: 0,
                    origin: SplitOrigin::Softclip(base + i),
                    seed_reads: 0,
                });
            }
        };
    push(&mut splits, 100, 2_443, true, 20_000, false); // terminus↔terminus, different copies
    push(&mut splits, 200, 1_200, false, 2_443, true); // interior↔terminus, same copy
    push(&mut splits, 300, 1_300, true, 1_400, false); // interior↔interior, same copy

    let contig_results = vec![ContigPileup {
        columns: fasta[0]
            .seq
            .iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: vec![],
                insertions: vec![],
                deletions: vec![],
            })
            .collect(),
        unique_depth: vec![10; 30_000],
        total_depth: vec![10; 30_000],
        splits,
    }];
    let gd = emit_from_pileup(&fasta, contig_results, &opts, &[]);
    let jcs = gd_jcs(&gd);
    let mut seen: Vec<(u64, u64)> = jcs
        .iter()
        .map(|j| {
            (
                j.fields[1].parse::<u64>().unwrap(),
                j.fields[4].parse::<u64>().unwrap(),
            )
        })
        .collect();
    seen.sort_unstable();
    assert_eq!(
        seen,
        vec![(1_200, 2_443), (2_443, 20_000)],
        "anchored repeat↔repeat junctions must be emitted, intra-copy noise dropped"
    );
    assert!(jcs
        .iter()
        .all(|j| j.attrs.get("pd_repeat_junction") == Some(&"1".to_string())));
}

#[test]
fn mob_promotion_rejects_repeat_side_targets() {
    // A MOB target is a unique-mapping site (docs/parity.md red line). Two
    // junctions whose shared "target" flank lies INSIDE an annotated
    // repeat copy are two repeat-side placements of an IS-boundary event
    // and must not promote a MOB.
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'G'; 5_000],
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 500,
            end: 1_500,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 4_000,
            end: 5_000,
        },
    ];
    let opts = EngineOptions {
        repeats,
        ..EngineOptions::default()
    };
    let mut splits = Vec::new();
    let push =
        |splits: &mut Vec<SplitCandidate>, base: usize, p1: u64, m1: bool, p2: u64, m2: bool| {
            for (i, minus) in [false, false, false, true, true, true]
                .into_iter()
                .enumerate()
            {
                splits.push(SplitCandidate {
                    contig_idx: 0,
                    side2_contig_idx: 0,
                    minus,
                    side1_minus: m1,
                    side2_minus: m2,
                    left: SubAlignment {
                        read_start: 0,
                        read_end: 18,
                    },
                    right: SubAlignment {
                        read_start: 18,
                        read_end: 36,
                    },
                    side1_pos_1: p1,
                    side2_pos_1: p2,
                    overlap: 0,
                    origin: SplitOrigin::Softclip(base + i),
                    seed_reads: 0,
                });
            }
        };
    // Shared in-repeat "target" flanks at 1200/1198 (opposing strands,
    // would-be dup = 3), other sides on copy termini of the same family.
    push(&mut splits, 100, 1_200, true, 4_000, false);
    push(&mut splits, 200, 1_198, false, 5_000, false);

    let contig_results = vec![ContigPileup {
        columns: fasta[0]
            .seq
            .iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: vec![],
                insertions: vec![],
                deletions: vec![],
            })
            .collect(),
        unique_depth: vec![10; 5_000],
        total_depth: vec![10; 5_000],
        splits,
    }];
    let gd = emit_from_pileup(&fasta, contig_results, &opts, &[]);
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert!(mobs.is_empty(), "in-repeat target must not promote a MOB");
}

fn push_flank(
    splits: &mut Vec<SplitCandidate>,
    n_per_strand: usize,
    origin: usize,
    p1: u64,
    m1: bool,
    p2: u64,
    m2: bool,
) {
    for i in 0..n_per_strand {
        for (k, minus) in [false, true].into_iter().enumerate() {
            splits.push(SplitCandidate {
                contig_idx: 0,
                side2_contig_idx: 0,
                minus,
                side1_minus: m1,
                side2_minus: m2,
                left: SubAlignment {
                    read_start: 0,
                    read_end: 18,
                },
                right: SubAlignment {
                    read_start: 18,
                    read_end: 36,
                },
                side1_pos_1: p1,
                side2_pos_1: p2,
                overlap: 0,
                origin: SplitOrigin::Softclip(origin + i * 2 + k),
                seed_reads: 0,
            });
        }
    }
}

fn emit_prepared(
    seq: Vec<u8>,
    repeats: Vec<RepeatRegion>,
    splits: Vec<SplitCandidate>,
) -> prokadiff_gd::GenomeDiff {
    let n = seq.len();
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let opts = EngineOptions {
        repeats,
        ..EngineOptions::default()
    };
    let contig_results = vec![ContigPileup {
        columns: fasta[0]
            .seq
            .iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: vec![],
                insertions: vec![],
                deletions: vec![],
            })
            .collect(),
        unique_depth: vec![10; n],
        total_depth: vec![10; n],
        splits,
    }];
    emit_from_pileup(&fasta, contig_results, &opts, &[])
}

fn emit_flank_gd(
    seq: Vec<u8>,
    repeats: Vec<RepeatRegion>,
    flanks: &[(usize, usize, u64, bool, u64, bool)],
) -> prokadiff_gd::GenomeDiff {
    let mut splits = Vec::new();
    for &(n_per_strand, origin, p1, m1, p2, m2) in flanks {
        push_flank(&mut splits, n_per_strand, origin, p1, m1, p2, m2);
    }
    emit_prepared(seq, repeats, splits)
}

fn jc_has_side(gd: &prokadiff_gd::GenomeDiff, pos: u64, strand: &str) -> bool {
    gd.entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Jc)
        .any(|j| {
            (j.fields[1] == pos.to_string() && j.fields[2] == strand)
                || (j.fields[4] == pos.to_string() && j.fields[5] == strand)
        })
}

#[test]
fn mob_boundary_shortfall_survives_terminus_snap() {
    // Repeat side sits 2 bp before a forward-copy start. Terminus snap must
    // still land the emitted JC on the start, but only AFTER the equal
    // shortfall has shifted p_plus. Snapping first yields dup = 1 and drops
    // the MOB. Bases around the unshifted register are unequal so the IS150
    // 1–2 bp rescue cannot invent the duplication on its own.
    let mut seq = vec![b'A'; 2_500];
    seq[1_698] = b'C';
    seq[1_699] = b'A';
    seq[1_700] = b'G';
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS150".into(),
        strand: 1,
        start: 500,
        end: 1_500,
    }];
    let gd = emit_flank_gd(
        seq,
        repeats,
        &[
            (3, 100, 498, false, 1_700, false),
            (3, 200, 1_500, true, 1_700, true),
        ],
    );
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(
        mobs.len(),
        1,
        "pre-terminus shortfall must still promote a MOB"
    );
    assert_eq!(mobs[0].fields[1], "1698");
    assert_eq!(mobs[0].fields[2], "IS150");
    assert_eq!(mobs[0].fields[3], "-1");
    assert_eq!(mobs[0].fields[4], "3");
    assert!(
        jc_has_side(&gd, 500, "1"),
        "emitted JC must still snap onto the terminus"
    );
    assert!(
        !jc_has_side(&gd, 498, "1"),
        "pre-snap coordinate must not be emitted"
    );
}

#[test]
fn is150_dup5_register_shifts_plus_side_onto_dup3() {
    // Plus and minus targets are 4 bp apart (duplication 5). That is the
    // 2 bp overlap register of a 3 bp IS150 target site. The MOB and the
    // plus-side JC move forward by 2 bp. A duplication of 4 is left alone.
    let seq = vec![b'A'; 3_000];
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS150".into(),
        strand: 1,
        start: 500,
        end: 1_500,
    }];
    let gd = emit_flank_gd(
        seq,
        repeats.clone(),
        &[
            (3, 100, 500, false, 2_000, false),
            (3, 200, 1_500, true, 2_004, true),
        ],
    );
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(mobs.len(), 1);
    assert_eq!(mobs[0].fields[1], "2002");
    assert_eq!(mobs[0].fields[2], "IS150");
    assert_eq!(mobs[0].fields[4], "3");
    assert!(
        jc_has_side(&gd, 2_002, "1"),
        "plus JC must follow the corrected register, entries={:?}",
        jc_fields(&gd)
    );
    assert!(!jc_has_side(&gd, 2_000, "1"));

    let mut seq4 = vec![b'A'; 3_000];
    seq4[1_999] = b'C';
    seq4[2_002] = b'T';
    let gd4 = emit_flank_gd(
        seq4,
        repeats,
        &[
            (3, 100, 500, false, 2_000, false),
            (3, 200, 1_500, true, 2_003, true),
        ],
    );
    let mob4: Vec<_> = gd4
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(mob4.len(), 1, "dup 4 must still promote");
    assert_eq!(mob4[0].fields[1], "2000");
    assert_eq!(mob4[0].fields[4], "4");
    assert!(jc_has_side(&gd4, 2_000, "1"));
}

#[test]
fn mob_flanks_follow_stronger_copy() {
    // Opposite element ends on different copies. Wave 8 Target C: the MOB
    // pair itself demonstrates that both flanks belong to one insertion
    // event, so the weaker 5' is displayed on the stronger copy's opposite
    // end even when that placement was never a retained alternative. The
    // unique side is untouched and the MOB duplication is unchanged.
    // Synthetic intervals only.
    let seq = vec![b'G'; 10_000];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 2_443,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: -1,
            start: 7_000,
            end: 8_443,
        },
    ];
    let gd = emit_flank_gd(
        seq,
        repeats,
        &[
            (3, 100, 8_443, true, 4_000, false),
            (6, 400, 2_443, true, 4_002, true),
        ],
    );
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(mobs.len(), 1);
    assert_eq!(mobs[0].fields[1], "4000");
    assert_eq!(mobs[0].fields[4], "3");
    assert!(
        jc_has_side(&gd, 1_000, "1"),
        "weaker flank adopts the stronger copy's opposite (5') end"
    );
    assert!(
        !jc_has_side(&gd, 8_443, "-1"),
        "the weaker flank's original copy must not remain displayed"
    );
    assert!(
        jc_has_side(&gd, 2_443, "-1"),
        "stronger 3' must stay on its own copy"
    );
    assert!(
        jc_links(&gd, 1_000, "1", 4_000, "1"),
        "weak flank keeps its unique target side"
    );
}

#[test]
fn mob_weak_flank_tie_support_not_synthesized() {
    // Wave 9: synthesis is reserved for genuinely weaker flanks. When the
    // weak flank's support is within COPY_SUPPORT_TIE_GAP of the strong
    // flank, the pair does not demonstrate that the strong flank's copy is
    // the inserted one, so the weak flank keeps its observed placement —
    // no synthesized move onto the strong copy (and no retained
    // alternative exists here either).
    let seq = vec![b'G'; 10_000];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 2_443,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: -1,
            start: 7_000,
            end: 8_443,
        },
    ];
    // Weak flank support 8 vs strong flank support 10: inside the tie
    // window, unlike mob_flanks_follow_stronger_copy's 6 vs 12.
    let gd = emit_flank_gd(
        seq,
        repeats,
        &[
            (4, 100, 8_443, true, 4_000, false),
            (5, 400, 2_443, true, 4_002, true),
        ],
    );
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(
        mobs.len(),
        1,
        "the MOB pair itself is unaffected by the display guard"
    );
    assert!(
        jc_has_side(&gd, 8_443, "-1"),
        "tie-supported weak flank keeps its observed copy"
    );
    assert!(
        !jc_has_side(&gd, 1_000, "1"),
        "no synthesized move onto the stronger copy's opposite end"
    );
    assert!(
        jc_links(&gd, 8_443, "-1", 4_000, "1"),
        "weak flank stays linked to its unique target"
    );
}

#[test]
fn mob_flanks_keep_copies_when_support_ties() {
    let seq = vec![b'G'; 10_000];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 2_443,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: -1,
            start: 7_000,
            end: 8_443,
        },
    ];
    let gd = emit_flank_gd(
        seq,
        repeats,
        &[
            (3, 100, 8_443, true, 4_000, false),
            (3, 400, 2_443, true, 4_002, true),
        ],
    );
    assert!(
        jc_has_side(&gd, 8_443, "-1"),
        "equal support must not reassign the repeat copy"
    );
    assert!(!jc_has_side(&gd, 1_000, "1"));
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(mobs.len(), 1);
    assert_eq!(mobs[0].fields[4], "3");
}

#[test]
fn nearer_copy_register_shift_keeps_presnap_shortfall_mob() {
    // Fold, shortfall, and snap interact. Two placements of the minus
    // flank share reads, so they are one junction: the farther copy uses
    // the plus flank's unique register, the nearer copy is 1 bp lower.
    // Own scores tie and the nearer copy is displayed while geometry is
    // measured. The plus repeat side is still 2 bp before its start.
    // Rejecting p_plus > p_minus before that shortfall drops the MOB;
    // snapping first would erase the shortfall too. Bases block a
    // leftward IS150 slide and allow one rightward step, so the
    // duplication lands on the shortfall position with dup = 3.
    // Deduped minus support is below the plus flank, and the farther
    // terminus was retained, so the displayed flank may move onto that
    // alternative after the MOB is recorded. Synthetic intervals only.
    let mut seq = vec![b'A'; 4_000];
    seq[1_996] = b'C';
    seq[1_998] = b'G';
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 500,
            end: 1_500,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_600,
            end: 1_800,
        },
    ];
    let gd = emit_flank_gd(
        seq,
        repeats,
        &[
            (6, 100, 498, false, 2_000, false),
            (3, 200, 1_500, true, 2_000, true),
            (3, 200, 1_800, true, 1_999, true),
        ],
    );
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(
        mobs.len(),
        1,
        "a 1 bp nearer-copy register shift must not drop the shortfall MOB"
    );
    assert_eq!(mobs[0].fields[1], "1998");
    assert_eq!(mobs[0].fields[2], "IS150");
    assert_eq!(mobs[0].fields[3], "-1");
    assert_eq!(mobs[0].fields[4], "3");
    assert!(
        jc_has_side(&gd, 1_500, "-1"),
        "displayed minus flank must be the retained alternative on the stronger copy"
    );
    assert!(
        !jc_has_side(&gd, 1_800, "-1"),
        "nearer copy is an alternative used for geometry, not the final display"
    );
    assert!(
        jc_has_side(&gd, 500, "1"),
        "emitted plus flank must still snap onto the terminus"
    );
    assert!(
        !jc_has_side(&gd, 498, "1"),
        "pre-snap shortfall coordinate must not be emitted"
    );
}

#[test]
fn mob_flank_relocates_only_onto_retained_alternative() {
    // The weaker 5' was also observed on the stronger copy. That
    // alternative is retained through the fold, so the displayed flank
    // may move there. The unobserved coordinate is not what is written.
    let seq = vec![b'G'; 10_000];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 2_443,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: -1,
            start: 7_000,
            end: 8_443,
        },
    ];
    let gd = emit_flank_gd(
        seq,
        repeats,
        &[
            (3, 100, 8_443, true, 4_000, false),
            (1, 100, 1_000, false, 4_000, false),
            (6, 400, 2_443, true, 4_002, true),
        ],
    );
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(mobs.len(), 1);
    assert_eq!(mobs[0].fields[1], "4000");
    assert_eq!(mobs[0].fields[4], "3");
    assert!(
        jc_has_side(&gd, 1_000, "1"),
        "retained alternative on the stronger copy may be displayed"
    );
    assert!(
        !jc_has_side(&gd, 8_443, "-1"),
        "the lower-scoring copy is not the displayed flank"
    );
}

#[allow(clippy::too_many_arguments)]
fn scored_copy(
    unique: u64,
    repeat_pos: u64,
    repeat_minus: bool,
    plus: usize,
    minus: usize,
    seed: usize,
    ids: Vec<usize>,
    overlap: i64,
) -> super::jc_cluster::AcceptedJc {
    super::jc_cluster::AcceptedJc::from_parts(
        0,
        unique,
        false,
        0,
        repeat_pos,
        repeat_minus,
        overlap,
        crate::jc::JunctionSupport {
            plus_reads: plus,
            minus_reads: minus,
            best_min_overlap: 18,
            plus_best_min: 18,
            minus_best_min: 16,
            min_overlap_side1: 16,
            min_overlap_side2: 15,
        },
        ids.into_iter()
            .map(|id| (0, SplitOrigin::Softclip(id)))
            .collect(),
        false,
        seed,
    )
}

fn is150_pair() -> (Vec<FastaRecord>, Vec<RepeatRegion>) {
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'A'; 30_000],
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 2_443,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 20_000,
            end: 21_443,
        },
    ];
    (fasta, repeats)
}

#[test]
fn partitioned_unique_flank_keeps_nearer_end_despite_far_seed() {
    // One unique flank, two identical 5' ends, disjoint read ids. The far
    // copy's seed is larger. The flank is still one junction, displayed on
    // the nearer end. This is the equal-support copy choice that dropped
    // a protected unique↔terminus match when seed outranked distance.
    let (fasta, repeats) = is150_pair();
    let near = scored_copy(5_000, 1_000, false, 4, 4, 1, vec![1, 2, 3, 4], 0);
    let far = scored_copy(5_000, 20_000, false, 4, 4, 40, vec![5, 6, 7, 8], 0);
    let folded = super::jc_cluster::fold_multicopy_placements(vec![far, near], &repeats, &fasta);
    assert_eq!(folded.len(), 1, "partitioned copies of one flank must fold");
    assert_eq!(folded[0].p1, 5_000);
    assert_eq!(
        folded[0].p2, 1_000,
        "farther higher-seed copy must not replace the nearer end"
    );
}

#[test]
fn each_unique_flank_keeps_its_own_nearer_copy() {
    // Two flanks of the same element end. Each keeps the copy closer to
    // its own unique side. A higher seed on the other copy must not pull
    // both flanks onto one reference copy.
    let (fasta, repeats) = is150_pair();
    let folded = super::jc_cluster::fold_multicopy_placements(
        vec![
            scored_copy(4_500, 20_000, false, 4, 4, 40, vec![1, 2, 3, 4], 0),
            scored_copy(4_500, 1_000, false, 4, 4, 1, vec![1, 2, 3, 4], 0),
            scored_copy(18_000, 1_000, false, 4, 4, 40, vec![11, 12, 13, 14], 0),
            scored_copy(18_000, 20_000, false, 4, 4, 1, vec![11, 12, 13, 14], 0),
        ],
        &repeats,
        &fasta,
    );
    assert_eq!(folded.len(), 2, "distinct unique flanks must stay distinct");
    let mut shown: Vec<(u64, u64)> = folded.iter().map(|j| (j.p1, j.p2)).collect();
    shown.sort_unstable();
    assert_eq!(shown, vec![(4_500, 1_000), (18_000, 20_000)]);
}

#[test]
fn shortfall_register_beats_nearer_exact_terminus_on_support_tie() {
    // The farther placement is the 2 bp pre-terminus shortfall MOB
    // promotion needs. An equal-support nearer exact terminus with a
    // higher seed must not replace it.
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'A'; 6_000],
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 500,
            end: 1_500,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 3_600,
            end: 4_600,
        },
    ];
    let folded = super::jc_cluster::fold_multicopy_placements(
        vec![
            scored_copy(4_000, 3_600, false, 4, 4, 40, vec![1, 2, 3, 4], 0),
            scored_copy(4_000, 498, false, 4, 4, 1, vec![1, 2, 3, 4], 0),
        ],
        &repeats,
        &fasta,
    );
    assert_eq!(folded.len(), 1);
    assert_eq!(
        folded[0].p2, 498,
        "pre-terminus shortfall must stay the displayed register"
    );
}

#[test]
fn same_element_end_pair_does_not_steal_forward_mob() {
    // High-support junctions put the same element end on both sides of
    // the target. That pair must not become the MOB: it records the
    // opposite strand and then blocks the opposite-end pair. The kept
    // MOB uses the opposite ends, strand +1, duplication 6. Synthetic
    // intervals only.
    let seq = vec![b'A'; 10_000];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS186".into(),
            strand: 1,
            start: 1_000,
            end: 2_342,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS186".into(),
            strand: 1,
            start: 4_000,
            end: 5_342,
        },
    ];
    let gd = emit_flank_gd(
        seq,
        repeats,
        &[
            (10, 10, 1_000, false, 8_000, false),
            (10, 500, 4_000, false, 8_005, true),
            (4, 1_000, 2_342, true, 8_005, true),
            (8, 1_500, 2_342, true, 8_000, false),
            (8, 2_000, 1_000, false, 8_005, true),
        ],
    );
    let mobs: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
        .collect();
    assert_eq!(mobs.len(), 1, "same-end flanks must not add a second MOB");
    assert_eq!(mobs[0].fields[1], "8000");
    assert_eq!(mobs[0].fields[2], "IS186");
    assert_eq!(mobs[0].fields[3], "1");
    assert_eq!(mobs[0].fields[4], "6");
    let parents: Vec<_> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Jc && mobs[0].parent_ids.contains(&e.id))
        .collect();
    assert!(
        parents.iter().any(|j| {
            (j.fields[1] == "8000"
                && j.fields[2] == "1"
                && j.fields[4] == "2342"
                && j.fields[5] == "-1")
                || (j.fields[4] == "8000"
                    && j.fields[5] == "1"
                    && j.fields[1] == "2342"
                    && j.fields[2] == "-1")
        }),
        "plus target must be parented on the 3' end, parents={:?}",
        parents.iter().map(|j| j.fields.clone()).collect::<Vec<_>>()
    );
}

#[test]
fn emitted_flanks_display_each_nearer_copy_when_far_seed_is_higher() {
    // Emitted coordinates, not only the fold struct. Each unique flank
    // shows the nearer 5' even though the other copy has the higher seed.
    let mut splits = Vec::new();
    push_flank(&mut splits, 4, 0, 4_500, false, 2_000, false);
    push_flank(&mut splits, 4, 0, 4_500, false, 8_000, false);
    for sp in &mut splits {
        if sp.side1_pos_1 == 8_000 || sp.side2_pos_1 == 8_000 {
            sp.seed_reads = 40;
        }
    }
    let before = splits.len();
    push_flank(&mut splits, 4, 100, 7_200, false, 2_000, false);
    push_flank(&mut splits, 4, 100, 7_200, false, 8_000, false);
    for sp in splits.iter_mut().skip(before) {
        if sp.side1_pos_1 == 2_000 || sp.side2_pos_1 == 2_000 {
            sp.seed_reads = 40;
        }
    }
    let gd = emit_prepared(
        vec![b'A'; 12_000],
        vec![
            RepeatRegion {
                seq_id: "chr".into(),
                name: "IS150".into(),
                strand: 1,
                start: 2_000,
                end: 3_442,
            },
            RepeatRegion {
                seq_id: "chr".into(),
                name: "IS150".into(),
                strand: 1,
                start: 8_000,
                end: 9_442,
            },
        ],
        splits,
    );
    assert!(
        jc_links(&gd, 4_500, "1", 2_000, "1"),
        "first flank must keep the nearer end, entries={:?}",
        jc_fields(&gd)
    );
    assert!(
        !jc_links(&gd, 4_500, "1", 8_000, "1"),
        "higher seed must not move the first flank onto the far end"
    );
    assert!(
        jc_links(&gd, 7_200, "1", 8_000, "1"),
        "second flank must keep the end nearer to itself, entries={:?}",
        jc_fields(&gd)
    );
    assert!(!jc_links(&gd, 7_200, "1", 2_000, "1"));
}

#[test]
fn modest_support_gap_keeps_nearer_identical_end() {
    // One extra read on the farther identical 5' is not copy-discriminating.
    // The unique flank stays on the nearer end. Disjoint ids, synthetic
    // intervals only.
    let (fasta, repeats) = is150_pair();
    let near = scored_copy(
        5_000,
        1_000,
        false,
        5,
        5,
        1,
        vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        0,
    );
    let far = scored_copy(
        5_000,
        20_000,
        false,
        6,
        5,
        40,
        vec![11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
        0,
    );
    let folded = super::jc_cluster::fold_multicopy_placements(vec![far, near], &repeats, &fasta);
    assert_eq!(folded.len(), 1);
    assert_eq!(folded[0].p2, 1_000);
}

#[test]
fn support_gap_above_tie_window_keeps_stronger_copy() {
    let (fasta, repeats) = is150_pair();
    let near = scored_copy(5_000, 1_000, false, 5, 5, 1, (1..=10).collect(), 0);
    let far = scored_copy(5_000, 20_000, false, 10, 6, 4, (101..=116).collect(), 0);
    let folded = super::jc_cluster::fold_multicopy_placements(vec![near, far], &repeats, &fasta);
    assert_eq!(folded.len(), 1);
    assert_eq!(
        folded[0].p2, 20_000,
        "a support gap above the tie window still selects that copy"
    );
}

#[test]
fn identical_constructs_keep_nearer_copy_and_drop_far_hub() {
    // Same unique flank, same element end, identical terminal sequence.
    // The far copy's seed is larger. Only the nearer construct is kept, so
    // the second-pass aligner cannot report the far copy as a second hub.
    let mut seq = vec![b'A'; 4_000];
    seq[99..107].copy_from_slice(b"ACGTACGT");
    seq[999..1007].copy_from_slice(b"TTTTTTTT");
    seq[2999..3007].copy_from_slice(b"TTTTTTTT");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 1_100,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 3_000,
            end: 3_100,
        },
    ];
    let mk = |repeat_pos: u64, seed: usize| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: 100,
        side1_minus: false,
        side2_contig: 0,
        side2_pos_1: repeat_pos,
        side2_minus: false,
        overlap: 0,
        seed_reads: seed,
    };
    let far = mk(3_000, 40);
    let near = mk(1_000, 2);
    assert_eq!(
        crate::jc_seq::junction_construct(&fasta, &far, 8).sequence(),
        crate::jc_seq::junction_construct(&fasta, &near, 8).sequence(),
        "fixture must be one junction sequence"
    );
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![far, near],
        &fasta,
        &repeats,
        8,
    );
    assert_eq!(
        collapsed.len(),
        1,
        "far identical end must not stay a second candidate"
    );
    assert_eq!(collapsed[0].side2_pos_1, 1_000);
    assert_eq!(collapsed[0].seed_reads, 40);
}

#[test]
fn identical_repeat_repeat_constructs_collapse_to_one_representative() {
    // Both sides are annotated repeats, and a uniform genome makes the two
    // constructs the same sequence. The aligner cannot distinguish them:
    // keeping both only splits the second-pass reads across identical
    // targets (Wave 8 Target A — the BR113 dilution failure). One
    // representative survives and carries the group's strongest seed.
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'C'; 4_000],
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 1_000,
            end: 1_050,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 2_000,
            end: 2_050,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 3_000,
            end: 3_050,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 3_500,
            end: 3_550,
        },
    ];
    let mk = |p1: u64, p2: u64, seed: usize| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: p1,
        side1_minus: false,
        side2_contig: 0,
        side2_pos_1: p2,
        side2_minus: true,
        overlap: 0,
        seed_reads: seed,
    };
    let left = mk(1_000, 3_000, 3);
    let right = mk(2_000, 3_500, 9);
    assert_eq!(
        crate::jc_seq::junction_construct(&fasta, &left, 8).sequence(),
        crate::jc_seq::junction_construct(&fasta, &right, 8).sequence(),
        "fixture must be one junction sequence"
    );
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![left, right],
        &fasta,
        &repeats,
        8,
    );
    assert_eq!(
        collapsed.len(),
        1,
        "identical repeat↔repeat constructs must not dilute reads"
    );
    assert_eq!(
        collapsed[0].seed_reads, 9,
        "representative carries the group's max seed"
    );
}

/// Fixture genome whose four loci hold `L`, `R`, `rc(R)`, `rc(L)` so two
/// candidates can be built whose constructs are exact reverse complements.
fn rc_pair_fixture() -> (Vec<FastaRecord>, Vec<u8>, Vec<u8>) {
    let mut seq = vec![b'A'; 4_000];
    let l = b"ACGTTGCA"; // left flank of construct 1
    let r = b"TTGCAACG"; // right flank of construct 1
    seq[200..208].copy_from_slice(l);
    seq[400..408].copy_from_slice(r);
    // rc(r) = CGTTGCAA, rc(l) = TGCAACGT
    seq[600..608].copy_from_slice(b"CGTTGCAA");
    seq[800..808].copy_from_slice(b"TGCAACGT");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    // cand1: side1 minus @208 → flank seq[200..208) = L;
    //        side2 plus  @401 → flank seq[400..408) = R.
    // cand2: side1 minus @608 → flank seq[600..608) = rc(R);
    //        side2 plus  @801 → flank seq[800..808) = rc(L).
    // construct1 = L + R; construct2 = rc(R) + rc(L) = rc(L + R).
    let construct1 = crate::jc_seq::junction_construct(
        &fasta,
        &crate::jc_seq::CandidateJunction {
            side1_contig: 0,
            side1_pos_1: 208,
            side1_minus: true,
            side2_contig: 0,
            side2_pos_1: 401,
            side2_minus: false,
            overlap: 0,
            seed_reads: 0,
        },
        8,
    )
    .sequence();
    let construct2 = crate::jc_seq::junction_construct(
        &fasta,
        &crate::jc_seq::CandidateJunction {
            side1_contig: 0,
            side1_pos_1: 608,
            side1_minus: true,
            side2_contig: 0,
            side2_pos_1: 801,
            side2_minus: false,
            overlap: 0,
            seed_reads: 0,
        },
        8,
    )
    .sequence();
    (fasta, construct1, construct2)
}

#[test]
fn reverse_complement_constructs_share_one_canonical_group() {
    // Wave 11: bowtie2 indexes both strands, so a construct and its
    // reverse complement are ONE double-stranded reference. Keying groups
    // by the raw sequence emitted both orientations into `junctions.fa`,
    // and the second-pass reads split across them. The strand-canonical
    // key (min(seq, rc(seq))) must put them in one group — and a
    // repeat↔repeat group then collapses to a single representative.
    let (fasta, construct1, construct2) = rc_pair_fixture();
    assert_ne!(construct1, construct2, "fixture must not be palindromic");
    assert_eq!(
        construct2,
        crate::jc_seq::rc_dna(&construct1),
        "fixture: the two constructs are reverse complements"
    );
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 200,
            end: 260,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 400,
            end: 460,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 600,
            end: 660,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 800,
            end: 860,
        },
    ];
    let mk = |s1: u64, s2: u64, seed: usize| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: s1,
        side1_minus: true,
        side2_contig: 0,
        side2_pos_1: s2,
        side2_minus: false,
        overlap: 0,
        seed_reads: seed,
    };
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![mk(208, 401, 4), mk(608, 801, 9)],
        &fasta,
        &repeats,
        8,
    );
    assert_eq!(
        collapsed.len(),
        1,
        "a construct and its reverse complement are one reference"
    );
    assert_eq!(collapsed[0].seed_reads, 9, "representative keeps max seed");
}

#[test]
fn reverse_complement_pair_with_distinct_unique_flanks_collapses() {
    // The production geometry (Wave 11): an ldr↔IS construct and a second
    // IS↔ldr construct at a DIFFERENT locus are reverse-complement
    // identical — the aligner cannot tell them apart, so keeping both
    // splits the molecule hits between them. They must collapse to one
    // representative even though each has its own distinct unique flank
    // (the one-flank guard alone would keep both).
    let (fasta, construct1, construct2) = rc_pair_fixture();
    assert_eq!(construct2, crate::jc_seq::rc_dna(&construct1));
    // Annotated copies start exactly at the repeat-side coordinates, so
    // each candidate classifies as unique↔element-end.
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 401,
            end: 460,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 801,
            end: 860,
        },
    ];
    let mk = |s1: u64, s2: u64, seed: usize| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: s1,
        side1_minus: true,
        side2_contig: 0,
        side2_pos_1: s2,
        side2_minus: false,
        overlap: 0,
        seed_reads: seed,
    };
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![mk(208, 401, 9), mk(608, 801, 4)],
        &fasta,
        &repeats,
        8,
    );
    assert_eq!(
        collapsed.len(),
        1,
        "reverse-complement duplicates must not split second-pass reads"
    );
    assert_eq!(collapsed[0].side1_pos_1, 208, "stronger seed is the rep");
    assert_eq!(collapsed[0].seed_reads, 9);
}

#[test]
fn distinct_unique_flanks_with_the_same_sequence_do_not_collapse() {
    let mut seq = vec![b'A'; 4_000];
    seq[99..107].copy_from_slice(b"ACGTACGT");
    seq[399..407].copy_from_slice(b"ACGTACGT");
    seq[999..1007].copy_from_slice(b"TTTTTTTT");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS150".into(),
        strand: 1,
        start: 1_000,
        end: 1_100,
    }];
    let mk = |unique: u64| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: unique,
        side1_minus: false,
        side2_contig: 0,
        side2_pos_1: 1_000,
        side2_minus: false,
        overlap: 0,
        seed_reads: 4,
    };
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![mk(100), mk(400)],
        &fasta,
        &repeats,
        8,
    );
    assert_eq!(collapsed.len(), 2, "two unique positions are not one flank");
}

#[test]
fn interior_repeat_side_member_does_not_block_identical_group_dedup() {
    // One member of an identical-sequence group lands its repeat side a
    // few bp INSIDE the copy, so the all-or-nothing
    // nearer_unique_repeat_representative classification fails for the
    // whole group (Wave 8: the BR96/1270660 group kept >20 identical
    // constructs this way). The group still collapses to one
    // representative anchored on the nearer annotated terminus.
    let mut seq = vec![b'A'; 4_000];
    seq[99..107].copy_from_slice(b"ACGTACGT");
    seq[999..1007].copy_from_slice(b"TTTTTTTT");
    seq[2999..3007].copy_from_slice(b"TTTTTTTT");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 1_100,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 3_000,
            end: 3_100,
        },
    ];
    let mk = |repeat_pos: u64, seed: usize| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: 100,
        side1_minus: false,
        side2_contig: 0,
        side2_pos_1: repeat_pos,
        side2_minus: false,
        overlap: 0,
        seed_reads: seed,
    };
    let far = mk(3_000, 40);
    let near = mk(1_000, 2);
    // Interior member: repeat side 30 bp inside the far copy. Its construct
    // is identical only if the terminal 8 bp repeat; force it by placing it
    // at the far copy start with the same 8-mer window.
    let interior = mk(3_030, 5);
    let seqs_equal = |a: &crate::jc_seq::CandidateJunction,
                      b: &crate::jc_seq::CandidateJunction| {
        crate::jc_seq::junction_construct(&fasta, a, 8).sequence()
            == crate::jc_seq::junction_construct(&fasta, b, 8).sequence()
    };
    // The interior member must NOT classify (repeat_containing path).
    assert!(seqs_equal(&far, &near));
    assert!(
        !seqs_equal(&interior, &far),
        "fixture: interior member is a different sequence"
    );
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![far.clone(), near.clone(), interior],
        &fasta,
        &repeats,
        8,
    );
    let same_seq: Vec<_> = collapsed
        .iter()
        .filter(|cand| seqs_equal(cand, &near))
        .collect();
    assert_eq!(
        same_seq.len(),
        1,
        "identical group with an interior-side sibling still collapses to one representative"
    );
    assert_eq!(same_seq[0].side2_pos_1, 1_000);
}

#[test]
fn unique_flank_registers_within_tolerance_cluster_into_one_construct() {
    // Wave 8 Target B: the same unique flank enumerated at 1-bp-shifted
    // registers joins the same element end. One construct must remain so
    // the second-pass reads are not split across near-identical targets.
    let mut seq = vec![b'A'; 4_000];
    seq[99..107].copy_from_slice(b"ACGTACGT");
    seq[999..1007].copy_from_slice(b"TTTTTTTT");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS150".into(),
        strand: 1,
        start: 1_000,
        end: 1_100,
    }];
    let mk = |unique: u64, seed: usize| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: unique,
        side1_minus: false,
        side2_contig: 0,
        side2_pos_1: 1_000,
        side2_minus: false,
        overlap: 0,
        seed_reads: seed,
    };
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![mk(100, 7), mk(101, 12), mk(103, 3)],
        &fasta,
        &repeats,
        8,
    );
    let same_family: Vec<_> = collapsed
        .iter()
        .filter(|cand| cand.side2_pos_1 == 1_000)
        .collect();
    assert_eq!(
        same_family.len(),
        1,
        "1–5 bp-shifted registers of one flank join one construct"
    );
    assert_eq!(
        same_family[0].side1_pos_1, 101,
        "strongest seed keeps its register"
    );
    assert_eq!(same_family[0].seed_reads, 12, "cluster carries max seed");
}

#[test]
fn flank_register_clusters_require_same_element_end() {
    // Two registers of one unique flank joining OPPOSITE element ends are
    // the two flanks of an insertion (a MOB pair), not one junction. They
    // must stay separate constructs.
    let mut seq = vec![b'A'; 4_000];
    seq[99..107].copy_from_slice(b"ACGTACGT");
    seq[989..1007].copy_from_slice(b"TTTTTTTTTTTTTTTTTT");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS150".into(),
        strand: 1,
        start: 1_000,
        end: 1_100,
    }];
    let five_prime = crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: 100,
        side1_minus: false,
        side2_contig: 0,
        side2_pos_1: 1_000,
        side2_minus: false,
        overlap: 0,
        seed_reads: 10,
    };
    let three_prime = crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: 101,
        side1_minus: true,
        side2_contig: 0,
        side2_pos_1: 1_100,
        side2_minus: true,
        overlap: 0,
        seed_reads: 10,
    };
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![five_prime, three_prime],
        &fasta,
        &repeats,
        8,
    );
    assert_eq!(
        collapsed.len(),
        2,
        "opposite element ends are a MOB pair, not one junction"
    );
}

#[test]
fn flank_register_clusters_require_same_strand_and_distance() {
    // Opposite-strand registers and registers farther than the cluster
    // tolerance are distinct junction evidence and stay separate.
    let mut seq = vec![b'A'; 8_000];
    seq[99..107].copy_from_slice(b"ACGTACGT");
    seq[599..607].copy_from_slice(b"ACGTACGT");
    seq[999..1007].copy_from_slice(b"TTTTTTTT");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS150".into(),
        strand: 1,
        start: 1_000,
        end: 1_100,
    }];
    let mk = |unique: u64, minus: bool, seed: usize| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: unique,
        side1_minus: minus,
        side2_contig: 0,
        side2_pos_1: 1_000,
        side2_minus: false,
        overlap: 0,
        seed_reads: seed,
    };
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![mk(100, false, 5), mk(101, true, 5), mk(600, false, 5)],
        &fasta,
        &repeats,
        8,
    );
    assert_eq!(
        collapsed.len(),
        3,
        "strand and 5 bp distance guards keep distinct flanks separate"
    );
}

#[test]
fn flank_register_clusters_require_same_repeat_copy() {
    // Wave 9 (Wave 8 Target B fallback): two registers of one unique
    // flank — same strand, within the cluster tolerance, same family,
    // same element end — that join DIFFERENT copies of the family are
    // distinct junctions, one construct per copy. Cross-copy merging
    // deletes a copy's true junction while a byte-identical sibling
    // captures its reads.
    let mut seq = vec![b'A'; 8_000];
    seq[99..107].copy_from_slice(b"ACGTACGT");
    seq[999..1007].copy_from_slice(b"TTTTTTTT");
    seq[4_999..5_007].copy_from_slice(b"TTTTTTTT");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 1_100,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 5_000,
            end: 5_100,
        },
    ];
    let mk = |unique: u64, terminus: u64, seed: usize| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: unique,
        side1_minus: false,
        side2_contig: 0,
        side2_pos_1: terminus,
        side2_minus: false,
        overlap: 0,
        seed_reads: seed,
    };
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![mk(100, 1_000, 5), mk(102, 5_000, 9)],
        &fasta,
        &repeats,
        8,
    );
    assert_eq!(
        collapsed.len(),
        2,
        "registers joining different repeat copies stay separate constructs"
    );
    let has = |terminus: u64| collapsed.iter().any(|cand| cand.side2_pos_1 == terminus);
    assert!(has(1_000), "copy 1000..1100 keeps its own junction");
    assert!(
        has(5_000),
        "copy 5000..5100 keeps its own junction; the stronger seed must not capture it"
    );
}

#[test]
fn interleaved_repeat_copy_registers_merge_same_copy_cleanly() {
    // Wave 10 / Codex Sol [P2] static counterexample. One unique flank is
    // enumerated at three registers: copy A at 100, copy B at 101, copy A
    // again at 102. When the cluster sort ran by flank POSITION first, the
    // copy-B register sat between the two copy-A registers; the linear
    // scan compares only the cluster's last member, so it flushed at the
    // interloper and produced THREE constructs (one per offset), not two.
    // Copy identity must sort before position so all registers joining
    // one copy are contiguous and merge within JC_CLUSTER_TOL_BP.
    let mut seq = vec![b'A'; 8_000];
    seq[99..107].copy_from_slice(b"ACGTACGT");
    seq[999..1007].copy_from_slice(b"TTTTTTTT");
    seq[2999..3007].copy_from_slice(b"TTTTTTTT");
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 1_100,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 3_000,
            end: 3_100,
        },
    ];
    let mk = |unique: u64, terminus: u64, seed: usize| crate::jc_seq::CandidateJunction {
        side1_contig: 0,
        side1_pos_1: unique,
        side1_minus: false,
        side2_contig: 0,
        side2_pos_1: terminus,
        side2_minus: false,
        overlap: 0,
        seed_reads: seed,
    };
    // Deliberately interleaved by unique position: A(100), B(101), A(102).
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        vec![mk(100, 1_000, 5), mk(101, 3_000, 6), mk(102, 1_000, 7)],
        &fasta,
        &repeats,
        8,
    );
    assert_eq!(
        collapsed.len(),
        2,
        "same-copy registers must merge even when another copy interleaves"
    );
    let copy_a: Vec<_> = collapsed
        .iter()
        .filter(|cand| cand.side2_pos_1 == 1_000)
        .collect();
    let copy_b: Vec<_> = collapsed
        .iter()
        .filter(|cand| cand.side2_pos_1 == 3_000)
        .collect();
    assert_eq!(copy_a.len(), 1, "copy A keeps exactly one representative");
    assert_eq!(copy_b.len(), 1, "copy B keeps exactly one representative");
    assert_eq!(
        copy_a[0].side1_pos_1, 102,
        "copy A representative is the strongest-seed register"
    );
    assert_eq!(
        copy_a[0].seed_reads, 7,
        "copy A representative carries the cluster's max seed"
    );
    assert_eq!(copy_b[0].side1_pos_1, 101);
}

#[test]
fn mob_weak_flank_stays_when_stronger_copy_is_farther() {
    // The synthesized move only shortens the unique↔terminus distance.
    // When the weaker flank already sits on the copy nearest its target,
    // moving it onto the stronger flank's copy would lengthen that
    // distance and the display stays put.
    let seq = vec![b'G'; 30_000];
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq,
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 1_000,
            end: 2_443,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS150".into(),
            strand: 1,
            start: 20_000,
            end: 21_443,
        },
    ];
    let opts = EngineOptions {
        repeats,
        jc_min_support_reads: 3,
        jc_min_frequency: 0.05,
        ..Default::default()
    };
    let mut splits = Vec::new();
    // Strong flank: copy A 3' end (2_443, minus) to far target
    // (25_002, minus).
    for (i, minus) in [false; 8].into_iter().chain([true; 8]).enumerate() {
        splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: true,
            side2_minus: true,
            left: SubAlignment {
                read_start: 0,
                read_end: 18,
            },
            right: SubAlignment {
                read_start: 18,
                read_end: 36,
            },
            side1_pos_1: 2_443,
            side2_pos_1: 25_002,
            overlap: 0,
            origin: SplitOrigin::Softclip(100 + i),
            seed_reads: 0,
        });
    }
    // Weak flank: copy B 5' end (20_000, plus) to target (25_000, plus) —
    // already on the copy NEAREST the target.
    for (i, minus) in [false; 3].into_iter().chain([true; 3]).enumerate() {
        splits.push(SplitCandidate {
            contig_idx: 0,
            side2_contig_idx: 0,
            minus,
            side1_minus: false,
            side2_minus: false,
            left: SubAlignment {
                read_start: 0,
                read_end: 18,
            },
            right: SubAlignment {
                read_start: 18,
                read_end: 36,
            },
            side1_pos_1: 20_000,
            side2_pos_1: 25_000,
            overlap: 0,
            origin: SplitOrigin::Softclip(200 + i),
            seed_reads: 0,
        });
    }
    let contig_results = vec![ContigPileup {
        columns: fasta[0]
            .seq
            .iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: vec![],
                insertions: vec![],
                deletions: vec![],
            })
            .collect(),
        unique_depth: vec![10; 30_000],
        total_depth: vec![10; 30_000],
        splits,
    }];
    let gd = emit_from_pileup(&fasta, contig_results, &opts, &[]);
    assert_eq!(
        gd.entries
            .iter()
            .filter(|e| e.kind == prokadiff_gd::GdKind::Mob)
            .count(),
        1
    );
    assert!(
        jc_links(&gd, 20_000, "1", 25_000, "1"),
        "weak flank on the nearer copy stays displayed there"
    );
    assert!(
        !jc_links(&gd, 1_000, "1", 25_000, "1"),
        "no synthesized move onto the farther strong copy"
    );
}

fn jc_fields(gd: &prokadiff_gd::GenomeDiff) -> Vec<Vec<String>> {
    gd.entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Jc)
        .map(|e| e.fields.clone())
        .collect()
}

fn jc_links(gd: &prokadiff_gd::GenomeDiff, p1: u64, s1: &str, p2: u64, s2: &str) -> bool {
    let a = p1.to_string();
    let b = p2.to_string();
    gd.entries.iter().any(|j| {
        j.kind == prokadiff_gd::GdKind::Jc
            && ((j.fields[1] == a && j.fields[2] == s1 && j.fields[4] == b && j.fields[5] == s2)
                || (j.fields[4] == a && j.fields[5] == s1 && j.fields[1] == b && j.fields[2] == s2))
    })
}

#[test]
fn own_score_outranks_aggregated_family_support() {
    // Three copies of one demonstrated junction. The first two would sum
    // to 20 under the old absorber and hide the third placement's own
    // support of 15. Ranking uses the frozen scores, and shared
    // observations are counted once.
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 1_000,
            end: 1_100,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 2_000,
            end: 2_100,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 3_000,
            end: 3_100,
        },
    ];
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'A'; 5_000],
    }];
    let ids = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    let mk = |pos: u64, plus: usize, minus: usize| {
        scored_copy(4_500, pos, false, plus, minus, 1, ids.clone(), 0)
    };
    let forward = super::jc_cluster::fold_multicopy_placements(
        vec![mk(1_000, 5, 5), mk(2_000, 5, 5), mk(3_000, 8, 7)],
        &repeats,
        &fasta,
    );
    let reverse = super::jc_cluster::fold_multicopy_placements(
        vec![mk(3_000, 8, 7), mk(2_000, 5, 5), mk(1_000, 5, 5)],
        &repeats,
        &fasta,
    );
    assert_eq!(forward.len(), 1);
    assert_eq!(reverse.len(), 1);
    assert_eq!(forward[0].p2, 3_000);
    assert_eq!(reverse[0].p2, 3_000);
    assert_eq!(forward[0].own_support_reads, 15);
    let family = forward[0].support.plus_reads + forward[0].support.minus_reads;
    assert!(
        family < 35,
        "duplicate observations must not be summed, family={family}"
    );
    assert_eq!(forward[0].alternatives.len(), 3);
}

#[test]
fn incompatible_end_or_overlap_does_not_fold() {
    let repeats = vec![RepeatRegion {
        seq_id: "chr".into(),
        name: "IS".into(),
        strand: 1,
        start: 1_000,
        end: 2_443,
    }];
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'A'; 10_000],
    }];
    let ids = vec![1, 2, 3, 4];
    let five = scored_copy(5_000, 1_000, false, 4, 4, 1, ids.clone(), 0);
    let three = scored_copy(5_002, 2_443, true, 4, 4, 1, ids.clone(), 0);
    let ends = super::jc_cluster::fold_multicopy_placements(vec![five, three], &repeats, &fasta);
    assert_eq!(
        ends.len(),
        2,
        "opposite element ends are different junctions"
    );

    let blunt = scored_copy(5_000, 1_000, false, 4, 4, 1, ids.clone(), 0);
    let gapped = scored_copy(5_000, 1_000, false, 4, 4, 1, ids, 12);
    let overlaps =
        super::jc_cluster::fold_multicopy_placements(vec![blunt, gapped], &repeats, &fasta);
    assert_eq!(overlaps.len(), 2, "unequal overlap is not one junction");
}

#[test]
fn reversed_copy_same_element_end_folds_with_shared_evidence() {
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 1_000,
            end: 2_000,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: -1,
            start: 7_000,
            end: 8_000,
        },
    ];
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'A'; 10_000],
    }];
    let ids = vec![1, 2, 3, 4];
    let forward = scored_copy(5_000, 1_000, false, 4, 4, 2, ids.clone(), 0);
    let reverse = super::jc_cluster::AcceptedJc::from_parts(
        0,
        5_002,
        false,
        0,
        8_000,
        true,
        0,
        crate::jc::JunctionSupport {
            plus_reads: 4,
            minus_reads: 4,
            best_min_overlap: 18,
            plus_best_min: 18,
            minus_best_min: 16,
            min_overlap_side1: 16,
            min_overlap_side2: 15,
        },
        ids.into_iter()
            .map(|id| (0, SplitOrigin::Softclip(id)))
            .collect(),
        false,
        2,
    );
    let folded =
        super::jc_cluster::fold_multicopy_placements(vec![forward, reverse], &repeats, &fasta);
    assert_eq!(
        folded.len(),
        1,
        "inverted copies of one element end may fold"
    );
    assert_eq!(folded[0].alternatives.len(), 2);
}

#[test]
fn equal_seed_unique_flank_prefers_nearer_repeat_copy() {
    // Support and Stage-2 seeds tie (identical termini). The reported copy
    // is the one closer to the unique flank, not the lower coordinate.
    // Circular distance: a copy just before the origin beats a farther
    // linear neighbor when the chromosome wraps.
    let fasta = vec![FastaRecord {
        name: "chr".into(),
        seq: vec![b'A'; 10_000],
    }];
    let repeats = vec![
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 500,
            end: 1_900,
        },
        RepeatRegion {
            seq_id: "chr".into(),
            name: "IS".into(),
            strand: 1,
            start: 9_950,
            end: 9_980,
        },
    ];
    let mk = |rep_pos: u64, seed: usize, id: usize| {
        super::jc_cluster::AcceptedJc::from_parts(
            0,
            200,
            false,
            0,
            rep_pos,
            false,
            0,
            crate::jc::JunctionSupport {
                plus_reads: 4,
                minus_reads: 4,
                best_min_overlap: 18,
                plus_best_min: 18,
                minus_best_min: 16,
                min_overlap_side1: 16,
                min_overlap_side2: 15,
            },
            HashSet::from([(0, SplitOrigin::Softclip(id))]),
            false,
            seed,
        )
    };
    let folded = super::jc_cluster::fold_multicopy_placements(
        vec![mk(500, 5, 1), mk(9_950, 5, 1)],
        &repeats,
        &fasta,
    );
    assert_eq!(folded.len(), 1);
    assert_eq!(
        folded[0].p2, 9_950,
        "circularly nearer terminus must win a support/seed tie"
    );

    // Unequal seeds on a support tie must not move the flank onto the
    // farther terminus. Seed remains only the last tie-break.
    let seeded = super::jc_cluster::fold_multicopy_placements(
        vec![mk(500, 9, 3), mk(9_950, 5, 3)],
        &repeats,
        &fasta,
    );
    assert_eq!(seeded.len(), 1);
    assert_eq!(
        seeded[0].p2, 9_950,
        "on a support tie, circular distance beats a stronger Stage-2 seed"
    );
    assert_eq!(seeded[0].own_seed_reads, 5);
    assert_eq!(seeded[0].seed_reads, 9);
}

#[test]
fn distant_interior_repeat_echo_is_dropped_nested_boundary_is_kept() {
    // Interior of one element × terminus of a distant copy is an echo of
    // one chimera on every identical copy. Drop it. The same interior
    // sitting inside the interval that owns the terminus is the nested
    // boundary and stays, tagged pd_repeat_junction.
    let far = emit_flank_gd(
        vec![b'A'; 25_000],
        vec![
            RepeatRegion {
                seq_id: "chr".into(),
                name: "ISa".into(),
                strand: 1,
                start: 1_000,
                end: 2_443,
            },
            RepeatRegion {
                seq_id: "chr".into(),
                name: "ISb".into(),
                strand: 1,
                start: 20_000,
                end: 21_443,
            },
        ],
        &[
            (3, 100, 1_500, false, 20_000, false),
            (3, 300, 5_000, false, 1_000, false),
        ],
    );
    assert!(
        !jc_has_side(&far, 1_500, "1"),
        "distant interior↔terminus echo must not be emitted"
    );
    assert!(
        jc_has_side(&far, 5_000, "1"),
        "unique↔terminus junction must stay"
    );

    let nested = emit_flank_gd(
        vec![b'A'; 6_000],
        vec![
            RepeatRegion {
                seq_id: "chr".into(),
                name: "ISnested".into(),
                strand: -1,
                start: 1_500,
                end: 2_500,
            },
            RepeatRegion {
                seq_id: "chr".into(),
                name: "ISouter".into(),
                strand: 1,
                start: 1_000,
                end: 2_000,
            },
        ],
        &[(3, 100, 1_100, false, 2_500, true)],
    );
    let jcs = gd_jcs(&nested);
    assert_eq!(jcs.len(), 1, "nested boundary must be emitted");
    assert_eq!(
        jcs[0].attrs.get("pd_repeat_junction").map(String::as_str),
        Some("1")
    );
}

/// Diagnostic replay of the G1 Clonal parity run (frozen job artifacts +
/// CURRENT engine code). Read-only: never writes into the run directory.
///
/// Gated behind `PROKADIFF_G1_REPLAY=<run_dir>` (a directory containing
/// `testdata`-style layout: `prokadiff/aligned.bam`, `prokadiff/work/…`,
/// `breseq/output/output.gd`) plus `PROKADIFF_G1_REF=<REL606.gbk>`.
/// Skips silently when unset so ordinary `cargo test` is unaffected.
#[test]
#[ignore = "set PROKADIFF_G1_REPLAY + PROKADIFF_G1_REF to replay frozen Clonal artifacts"]
fn g1_replay_clonal_jc_gate() {
    use std::collections::HashMap;
    use std::path::PathBuf;

    let Some(run_dir) = std::env::var_os("PROKADIFF_G1_REPLAY").map(PathBuf::from) else {
        return;
    };
    let Some(ref_gbk) = std::env::var_os("PROKADIFF_G1_REF").map(PathBuf::from) else {
        eprintln!("PROKADIFF_G1_REPLAY set but PROKADIFF_G1_REF missing; skipping");
        return;
    };

    let repeats = crate::fasta::parse_genbank_repeats(&ref_gbk).unwrap_or_default();
    assert!(!repeats.is_empty(), "repeats table must parse");
    // Replay fidelity: use the reference FASTA the frozen bowtie2 run aligned
    // against (LOCUS-derived id `REL606`). The current `read_reference`
    // prefers the VERSION accession (RW-005), which would not match the
    // frozen BAM header.
    let fasta = crate::fasta::read_reference(run_dir.join("prokadiff/work/reference.fa"))
        .expect("reference fasta");
    let opts = EngineOptions {
        threads: 1,
        repeats,
        ..EngineOptions::default()
    };

    let bam = run_dir.join("prokadiff/aligned.bam");
    let primary = super::bam_io::read_primary_bam(&bam, &fasta).expect("primary bam");
    let pileup = super::pileup_contigs_with_repeat(&fasta, &primary.aligned, &opts);
    let work = run_dir.join("prokadiff/work");

    // Replicate `second_pass_splits` candidate pipeline against frozen
    // Stage-2 SAMs, then reuse the frozen junction alignment BAM.
    let flank = primary
        .aligned
        .iter()
        .map(|r| r.seq.len())
        .max()
        .unwrap_or(36)
        .max(28);
    let mut all_cands = super::jc_cluster::seed_candidates_from_pileup(&pileup.contigs);
    let contig_names: Vec<String> = fasta.iter().map(|r| r.name.clone()).collect();
    for sam_name in [
        "stage2_r1.sam",
        "stage2_r2.sam",
        "stage2_pe.sam",
        "stage2_se.sam",
        "stage2.sam",
    ] {
        let sam_path = work.join(sam_name);
        if sam_path.is_file() {
            let res = crate::split_seed::extract_candidate_junctions_with_repeat_evidence_from_sam(
                &sam_path,
                &contig_names,
            )
            .expect("stage2 sam");
            all_cands.extend(res.candidates);
        }
    }
    type JcMergeKey = (usize, u64, bool, usize, u64, bool, i64);
    let mut merged: HashMap<JcMergeKey, crate::jc_seq::CandidateJunction> = HashMap::new();
    for cand in all_cands {
        let (c1, p1, m1, c2, p2, m2) =
            if (cand.side1_contig, cand.side1_pos_1) <= (cand.side2_contig, cand.side2_pos_1) {
                (
                    cand.side1_contig,
                    cand.side1_pos_1,
                    cand.side1_minus,
                    cand.side2_contig,
                    cand.side2_pos_1,
                    cand.side2_minus,
                )
            } else {
                (
                    cand.side2_contig,
                    cand.side2_pos_1,
                    !cand.side2_minus,
                    cand.side1_contig,
                    cand.side1_pos_1,
                    !cand.side1_minus,
                )
            };
        merged
            .entry((c1, p1, m1, c2, p2, m2, cand.overlap))
            .and_modify(|e| e.seed_reads += cand.seed_reads)
            .or_insert(crate::jc_seq::CandidateJunction {
                side1_contig: c1,
                side1_pos_1: p1,
                side1_minus: m1,
                side2_contig: c2,
                side2_pos_1: p2,
                side2_minus: m2,
                overlap: cand.overlap,
                seed_reads: cand.seed_reads,
            });
    }
    let ref_len: usize = fasta.iter().map(|r| r.seq.len()).sum();
    eprintln!(
        "replay: primary.aligned={} merged_cands={}",
        primary.aligned.len(),
        merged.len()
    );
    let cands = crate::jc_seq::rank_and_cap(merged.into_values().collect(), ref_len, flank);
    eprintln!("replay: rank_and_cap kept {}", cands.len());
    let mut recs = Vec::new();
    let mut kept: Vec<(crate::jc_seq::CandidateJunction, usize)> = Vec::new();
    for cand in &cands {
        let construct = crate::jc_seq::junction_construct(&fasta, cand, flank);
        if construct.len() < crate::jc_seq::JC_MIN_COVER_BASES {
            continue;
        }
        recs.push(FastaRecord {
            name: format!("jc{}", kept.len()),
            seq: construct.sequence(),
        });
        kept.push((cand.clone(), construct.breakpoint()));
    }
    eprintln!("replay: {} second-pass candidates kept", kept.len());

    // Run the junction alignment ourselves (current code must be
    // self-consistent; the frozen jc/aligned.bam indices belonged to the
    // frozen candidate order). Outputs go to scratch only.
    let scratch = std::env::var_os("PROKADIFF_G1_SCRATCH")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("opencode/g1_replay"));
    std::fs::create_dir_all(&scratch).expect("scratch dir");
    let jc_dir = scratch.join("jc");
    std::fs::create_dir_all(&jc_dir).expect("jc scratch dir");
    let fa = jc_dir.join("junctions.fa");
    crate::fasta::write_records(&recs, &fa).expect("write constructs");
    let pass_reads = crate::align::FastqInput {
        files: vec![
            work.join("jc/reads/pass_0.fastq"),
            work.join("jc/reads/pass_1.fastq"),
        ],
    };
    let jc_bam = jc_dir.join("aligned.bam");
    crate::align::align_to_bam(
        &fa,
        &pass_reads,
        &jc_bam,
        1,
        &jc_dir,
        crate::align::AlignKind::Junction,
    )
    .expect("junction alignment");

    let (extra, stats) = super::bam_io::parse_junction_bam(
        &jc_bam,
        &kept,
        &primary.primary_scores,
        &opts.repeats,
        &contig_names,
    )
    .expect("junction bam");
    eprintln!("replay: second-pass stats {:?}", stats);

    let gd = super::emit::emit_from_pileup(&fasta, pileup.contigs, &opts, &extra);

    // Mini accepted-JC gate vs the breseq oracle in the same run dir.
    struct JcRow {
        id: String,
        pos1: u64,
        pos2: u64,
        s1: i8,
        s2: i8,
        ov: i64,
        rejected: bool,
    }
    let parse_gd_jc = |path: &PathBuf| -> Vec<JcRow> {
        std::fs::read_to_string(path)
            .expect("gd file")
            .lines()
            .filter_map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                if f.len() < 10 || f[0] != "JC" {
                    return None;
                }
                Some(JcRow {
                    id: f[1].to_string(),
                    pos1: f[4].parse().ok()?,
                    pos2: f[7].parse().ok()?,
                    s1: f[5].parse().ok()?,
                    s2: f[8].parse().ok()?,
                    ov: f[9].parse().ok()?,
                    rejected: f[10..].iter().any(|a| a.starts_with("reject=")),
                })
            })
            .collect()
    };
    let oracle = parse_gd_jc(&run_dir.join("breseq/output/output.gd"));
    let oracle_acc: Vec<&JcRow> = oracle.iter().filter(|r| !r.rejected).collect();

    let mut pd_rows: Vec<&prokadiff_gd::GdEntry> = gd
        .entries
        .iter()
        .filter(|e| e.kind == prokadiff_gd::GdKind::Jc)
        .collect();
    pd_rows.sort_by_key(|e| e.position().unwrap_or(0));

    let side_close = |p: u64, s: i8, q: u64, t: i8| s == t && p.abs_diff(q) <= 5;
    let mut matched = 0usize;
    let mut over = Vec::new();
    let mut used = vec![false; oracle_acc.len()];
    for e in &pd_rows {
        let f = &e.fields;
        let (p1, s1, p2, s2, ov): (u64, i64, u64, i64, i64) = (
            f[4].parse().unwrap(),
            f[5].parse().unwrap(),
            f[7].parse().unwrap(),
            f[8].parse().unwrap(),
            f[9].parse().unwrap(),
        );
        let mut hit = None;
        for (i, o) in oracle_acc.iter().enumerate() {
            if used[i] || o.ov != ov {
                continue;
            }
            let (a1, a2, b1, b2) = if p1 <= p2 {
                (p1, p2, o.pos1, o.pos2)
            } else {
                (p2, p1, o.pos2, o.pos1)
            };
            let (sa, sb, ta, tb) = if p1 <= p2 {
                (s1, s2, o.s1, o.s2)
            } else {
                (s2, s1, o.s2, o.s1)
            };
            if side_close(a1, sa as i8, b1, ta) && side_close(a2, sb as i8, b2, tb) {
                hit = Some(i);
                break;
            }
        }
        match hit {
            Some(i) => {
                used[i] = true;
                matched += 1;
                println!("MATCHED  PD{}\tBR{}", f[1], oracle_acc[i].id);
            }
            None => {
                let sup = e.attrs.get("pd_support_reads").cloned().unwrap_or_default();
                let best = e
                    .attrs
                    .get("pd_best_min_overlap")
                    .cloned()
                    .unwrap_or_default();
                over.push(format!(
                    "OVER  PD{}\t{}({})->{}({})\tov={}\tsup={}\tbest={}",
                    f[1], f[3], f[4], f[6], f[7], ov, sup, best
                ));
            }
        }
    }
    let under: Vec<String> = oracle_acc
        .iter()
        .enumerate()
        .filter(|(i, _)| !used[*i])
        .map(|(_, o)| format!("UNDER BR{}", o.id))
        .collect();
    for l in &over {
        println!("{l}");
    }
    for l in &under {
        println!("{l}");
    }
    println!(
        "G1_REPLAY_SUMMARY pd_jc={} matched={} over={} under={}",
        pd_rows.len(),
        matched,
        over.len(),
        under.len()
    );
}

/// TEMPORARY Wave-11 diagnostic (to be removed): faithful replay of the W11
/// Clonal run from frozen artifacts, pinpointing the stage where the
/// BR95/BR113/BR114-class junctions die.
#[test]
#[ignore = "set PROKADIFF_W11_DIAG=<run_dir> PROKADIFF_G1_REF=<gbk>"]
fn w11_diag_replay() {
    use std::collections::HashMap;
    use std::path::PathBuf;

    let Some(run_dir) = std::env::var_os("PROKADIFF_W11_DIAG").map(PathBuf::from) else {
        return;
    };
    let Some(ref_gbk) = std::env::var_os("PROKADIFF_G1_REF").map(PathBuf::from) else {
        eprintln!("PROKADIFF_G1_REF missing");
        return;
    };
    let repeats = crate::fasta::parse_genbank_repeats(&ref_gbk).unwrap_or_default();
    let fasta =
        crate::fasta::read_reference(run_dir.join("prokadiff/work/reference.fa")).unwrap();
    let opts = EngineOptions {
        threads: 1,
        repeats,
        ..EngineOptions::default()
    };
    let primary =
        super::bam_io::read_primary_bam(&run_dir.join("prokadiff/aligned.bam"), &fasta).unwrap();
    let mut pileup = super::pileup_contigs_with_repeat(&fasta, &primary.aligned, &opts);
    let work = run_dir.join("prokadiff/work");
    let flank = primary
        .aligned
        .iter()
        .map(|r| r.seq.len())
        .max()
        .unwrap_or(36)
        .max(28);
    let mut all_cands = super::jc_cluster::seed_candidates_from_pileup(&pileup.contigs);
    let contig_names: Vec<String> = fasta.iter().map(|r| r.name.clone()).collect();
    let contig_lens: Vec<u64> = fasta.iter().map(|r| r.seq.len() as u64).collect();
    for sam_name in [
        "stage2_r1.sam",
        "stage2_r2.sam",
        "stage2_pe.sam",
        "stage2_se.sam",
        "stage2.sam",
    ] {
        let sam_path = work.join(sam_name);
        if sam_path.is_file() {
            let res =
                crate::split_seed::extract_candidate_junctions_with_repeat_evidence_from_sam_anchored(
                    &sam_path,
                    &contig_names,
                    &opts.repeats,
                    &contig_lens,
                )
                .unwrap();
            all_cands.extend(res.candidates);
        }
    }
    type JcMergeKey = (usize, u64, bool, usize, u64, bool, i64);
    let mut merged: HashMap<JcMergeKey, crate::jc_seq::CandidateJunction> = HashMap::new();
    for cand in all_cands {
        let (c1, p1, m1, c2, p2, m2) =
            if (cand.side1_contig, cand.side1_pos_1) <= (cand.side2_contig, cand.side2_pos_1) {
                (
                    cand.side1_contig,
                    cand.side1_pos_1,
                    cand.side1_minus,
                    cand.side2_contig,
                    cand.side2_pos_1,
                    cand.side2_minus,
                )
            } else {
                (
                    cand.side2_contig,
                    cand.side2_pos_1,
                    !cand.side2_minus,
                    cand.side1_contig,
                    cand.side1_pos_1,
                    !cand.side1_minus,
                )
            };
        merged
            .entry((c1, p1, m1, c2, p2, m2, cand.overlap))
            .and_modify(|e| e.seed_reads += cand.seed_reads)
            .or_insert(crate::jc_seq::CandidateJunction {
                side1_contig: c1,
                side1_pos_1: p1,
                side1_minus: m1,
                side2_contig: c2,
                side2_pos_1: p2,
                side2_minus: m2,
                overlap: cand.overlap,
                seed_reads: cand.seed_reads,
            });
    }
    let collapsed = super::jc_cluster::collapse_identical_unique_repeat_copies(
        merged.into_values().collect(),
        &fasta,
        &opts.repeats,
        flank,
    );
    let ref_len: usize = fasta.iter().map(|r| r.seq.len()).sum();
    let cands = crate::jc_seq::rank_and_cap(collapsed, ref_len, flank);
    let mut kept: Vec<(crate::jc_seq::CandidateJunction, usize)> = Vec::new();
    for cand in &cands {
        let construct = crate::jc_seq::junction_construct(&fasta, cand, flank);
        if construct.len() < crate::jc_seq::JC_MIN_COVER_BASES {
            continue;
        }
        kept.push((cand.clone(), construct.breakpoint()));
    }
    // Faithfulness: rebuilt constructs must equal the frozen junctions.fa.
    let frozen = crate::fasta::read_fasta(work.join("jc/junctions.fa")).unwrap();
    assert_eq!(
        kept.len(),
        frozen.len(),
        "rebuilt candidate count differs from frozen junctions.fa"
    );
    for (i, (cand, _)) in kept.iter().enumerate() {
        let seq = crate::jc_seq::junction_construct(&fasta, cand, flank).sequence();
        assert_eq!(
            seq,
        frozen[i].seq,
            "construct {i} differs from frozen junctions.fa"
        );
    }
    let (extra, stats) = super::bam_io::parse_junction_bam(
        &work.join("jc/aligned.bam"),
        &kept,
        &primary.primary_scores,
        &opts.repeats,
        &contig_names,
    )
    .unwrap();
    eprintln!(
        "diag: parse stats considered={} kept={} short_cover={} not_spanning={} worse_primary={} worse_competing={}",
        stats.considered, stats.kept, stats.short_cover, stats.not_spanning,
        stats.worse_than_primary, stats.worse_than_competing
    );

    // --- accept stage ---
    for sp in &extra {
        if sp.contig_idx < pileup.contigs.len() {
            pileup.contigs[sp.contig_idx].splits.push(sp.clone());
        }
    }
    let targets: [(u64, u64); 3] = [
        (3652533, 3901931), // BR114
        (3651091, 3741531), // BR113
        (589886, 590471),   // BR95
    ];
    let near = |p1: u64, p2: u64| {
        targets.iter().any(|(t1, t2)| {
            p1.abs_diff(*t1) <= 40 && p2.abs_diff(*t2) <= 40
                || p1.abs_diff(*t2) <= 40 && p2.abs_diff(*t1) <= 40
        })
    };
    let mut accepted_all: Vec<super::jc_cluster::AcceptedJc> = Vec::new();
    for (idx, _rec) in fasta.iter().enumerate() {
        let splits = std::mem::take(&mut pileup.contigs[idx].splits);
        for ((c1, c2, _sgn), c) in super::jc_cluster::cluster_key_splits(&splits) {
            let acc = crate::jc::accept_junction(&c.support);
            if near(c.rep_side1, c.rep_side2) {
                eprintln!(
                    "diag cluster c1={c1} {}({}) -> c2={c2} {}({}) ov={} plus={} minus={} best_min={} plus_best={} minus_best={} min_s1={} min_s2={} accept={acc}",
                    fasta[c1].name, c.rep_side1, fasta[c2].name, c.rep_side2, c.rep_overlap,
                    c.support.plus_reads, c.support.minus_reads,
                    c.support.best_min_overlap, c.support.plus_best_min,
                    c.support.minus_best_min, c.support.min_overlap_side1,
                    c.support.min_overlap_side2,
                );
            }
            if !acc {
                continue;
            }
            accepted_all.push(super::jc_cluster::AcceptedJc::from_parts(
                c1, c.rep_side1, c.rep_side1_minus, c2, c.rep_side2, c.rep_side2_minus,
                c.rep_overlap, c.support, c.ids, c.has_cigar, c.seed_reads,
            ));
        }
        pileup.contigs[idx].splits = splits;
    }
    eprintln!("diag: accepted_all={}", accepted_all.len());
    let folded = super::jc_cluster::fold_multicopy_placements(
        super::jc_cluster::fold_reverse_duplicates(accepted_all),
        &opts.repeats,
        &fasta,
    );
    eprintln!("diag: folded={}", folded.len());
    for j in &folded {
        if near(j.p1, j.p2) {
            eprintln!(
                "diag folded {}({}) -> {}({}) ov={} plus={} minus={}",
                fasta[j.c1].name, j.p1, fasta[j.c2].name, j.p2, j.overlap,
                j.support.plus_reads, j.support.minus_reads
            );
        }
    }
}
