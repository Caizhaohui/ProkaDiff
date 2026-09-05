use std::collections::{HashMap, HashSet};

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
    let scores = HashMap::new();

    let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(10)], &scores).unwrap();
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
    let (extra_flank, stats_flank) = parse_junction_bam(&bam, &[jc_cand(36)], &scores).unwrap();
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
    let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
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
    let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
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
    let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
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

    let weaker = HashMap::from([("r1".to_string(), 30)]);
    let (extra, stats) = parse_junction_bam(&bam, &cands, &weaker).unwrap();
    assert_eq!(stats.kept, 1);
    assert_eq!(stats.worse_than_primary, 0);
    assert_eq!(extra.len(), 1);

    let stronger = HashMap::from([("r1".to_string(), 40)]);
    let (extra, stats) = parse_junction_bam(&bam, &cands, &stronger).unwrap();
    assert!(extra.is_empty());
    assert_eq!(stats.worse_than_primary, 1);
    assert_eq!(stats.not_spanning, 0);
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
    let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
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
    let (extra, stats) = parse_junction_bam(&bam, &[jc_cand(36)], &HashMap::new()).unwrap();
    assert!(extra.is_empty());
    assert_eq!(stats.considered, 0);
}
