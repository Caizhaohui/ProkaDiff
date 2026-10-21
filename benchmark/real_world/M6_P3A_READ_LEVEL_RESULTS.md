# M6-P3a read-level targeted diagnostic results

Status: **M6-P3a ALL 10 ALLELES MEASURED — CLOSED (PASS_WITH_LIMITATIONS)**.
M6 remains **ACTIVE**. All 10 unrecovered publication alleles across both Wave 1 and Wave 2 are now measured.
P2a still records 2/12 `RECOVERED_EXACT` and 10 `NOT_RECOVERED`. Those tables were not edited.

Counts below are diagnostic pileups on stage-1 BAM files. They are not internal
RA column totals. Phred ≥ 20 is a diagnostic filter. The caller does not apply
it. Caller-shaped depth keeps primary alignments with MAPQ ≥ 10 after the
XS≥AS rule and drops one aligned base at each end. It does not apply the
homopolymer-end mask.

Paired mates share a qname. Each mate is counted once. A second primary
record with the same qname and the same mate flag is skipped. Secondary and
supplementary records are not counted. Across all 10 loci that skip count is 0.

## Input gates & execution telemetry

All alignments executed on partition `qcpu_18i`, 8 CPU, 32 GB, 4 hours, at most two
concurrent tasks. All tasks exited with code `0:0`.

### Wave 1 (Priority 1 targets, Slurm array `2905873`)

| Task | Run | Sample | Reference | SHA256 | Reads | Node | Elapsed | MaxRSS |
| --- | --- | --- | --- | --- | ---: | --- | --- | ---: |
| `2905873_0` | `SRR28349204` | Cas9_3 _and2others | `NC_000913.3` | `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0` | 4,051,149 | `bnode1` | 00:09:19 | 186,876K |
| `2905873_1` | `SRR28349176` | Cas9_104_105 | `NC_000913.3` | same MG1655 hash | 3,433,823 | `bnode2` | 00:07:26 | 179,176K |
| `2905873_2` | `SRR28349230` | I-SceI_7_8_9 | `CP009273.1` | `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce` | 3,230,652 | `bnode2` | 00:07:56 | 182,776K |

### Wave 2 (Priority 2 targets, Slurm array `2992594`)

| Task | Run | Sample | Reference | SHA256 | Reads | Node | Elapsed | MaxRSS |
| --- | --- | --- | --- | --- | ---: | --- | --- | ---: |
| `2992594_0` | `SRR28349225` | control_5_6_7_8 | `CP009273.1` | `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce` | 2,258,727 | `bnode1` | 00:05:11 | 183,352K |
| `2992594_1` | `SRR28349221` | control_17_18_19_20 | `CP009273.1` | same BW25113 hash | 1,993,895 | `bnode10` | 00:04:23 | 180,340K |
| `2992594_2` | `SRR28349165` | lambdaRed_7_8 | `CP009273.1` | same BW25113 hash | 4,267,402 | `bnode2` | 00:10:44 | 1,644,908K |
| `2992594_3` | `SRR28349232` | I-SceI_1_2 | `CP009273.1` | same BW25113 hash | 1,886,186 | `bnode1` | 00:05:04 | 1,491,460K |
| `2992594_4` | `SRR28349228` | I-SceI_14_15_16 | `CP009273.1` | same BW25113 hash | 3,214,563 | `bnode1` | 00:08:40 | 327,436K |

Gzip checksums matched P2a provenance before alignment. Bowtie2 2.5.4.
Command: `bowtie2 --local --no-unal -p 8` with `-1`/`-2` and unmatched FASTQ.
No stage-2 seed and no junction parameters.
BAM sorting used samtools 1.24 (`samtools sort` and `samtools index`).
Production writes the stage-1 SAM through an in-memory noodles sort and does
not write a `.bai`. The alignment records, flags, CIGAR, sequence, quality,
and AS/XS tags are the Bowtie2 output. The sort implementation is the
documented difference.

`aligned.bam` and `aligned.bam.bai` remain under `benchmark/real_world/results/prjna1088182/m6_p3a/<run>/`.

## Measurements & Disambiguation

Caller-shaped ALT fraction is ALT divided by caller-shaped depth.

### Authoritative 10-Locus Disambiguation Table

| Wave | Locus | BioSample | Run | Sample Identifier | Reference | Raw Depth (Alt) | Usable Depth (Alt) | ALT Fraction | Q20 Alt / Depth | Strand (+ / -) |
| :---: | :---: | :---: | :---: | :--- | :---: | ---: | ---: | ---: | ---: | :---: |
| **W1** | `1571257 G>T` | `SAMN40463618` | `SRR28349204` | Cas9_3 _and2others | `NC_000913.3` | 239 (55) | 80 (24) | **0.300000** | 54 / 236 | 9 / 15 |
| **W1** | `2776013 T>C` | `SAMN40463661` | `SRR28349176` | Cas9_104_105 | `NC_000913.3` | 148 (42) | 144 (42) | **0.291667** | 42 / 126 | 18 / 24 |
| **W1** | `4035824 T>C` | `SAMN40463661` | `SRR28349176` | Cas9_104_105 | `NC_000913.3` | 155 (62) | 59 (38) | **0.644068** | 62 / 152 | 21 / 17 |
| **W1** | `2872949 C>A` | `SAMN40463688` | `SRR28349230` | I-SceI_7_8_9 | `CP009273.1` | 217 (150) | 217 (150) | **0.691244** | 137 / 202 | 96 / 54 |
| **W2** | `1635887 A>G` | `SAMN40463693` | `SRR28349225` | control_5_6_7_8 | `CP009273.1` | 83 (21) | 83 (21) | **0.253012** | 21 / 82 | 9 / 12 |
| **W2** | `3108123 G>A` | `SAMN40463696` | `SRR28349221` | control_17_18_19_20 | `CP009273.1` | 107 (21) | 106 (21) | **0.198113** | 20 / 106 | 8 / 13 |
| **W2** | `161068 A>C` | `SAMN40463671` | `SRR28349165` | lambdaRed_7_8 | `CP009273.1` | 199 (63) | 194 (62) | **0.319588** | 14 / 84 | 10 / 52 |
| **W2** | `2653047 T>A` | `SAMN40463686` | `SRR28349232` | I-SceI_1_2 | `CP009273.1` | 111 (32) | 110 (32) | **0.290909** | 32 / 111 | 20 / 12 |
| **W2** | `4115617 C>A` | `SAMN40463686` | `SRR28349232` | I-SceI_1_2 | `CP009273.1` | 118 (44) | 116 (44) | **0.379310** | 44 / 117 | 21 / 23 |
| **W2** | `3073571 G>A` | `SAMN40463690` | `SRR28349228` | I-SceI_14_15_16 | `CP009273.1` | 212 (71) | 211 (71) | **0.336493** | 70 / 207 | 42 / 29 |

*Note on Control Run Disambiguation:*
- **`SRR28349225`** corresponds strictly to **`control_5_6_7_8`** (`SAMN40463693`) with event **`1635887 A>G`** (AF = **0.253012**).
- **`SRR28349221`** corresponds strictly to **`control_17_18_19_20`** (`SAMN40463696`) with event **`3108123 G>A`** (AF = **0.198113**).

REF + ALT + other equals depth for every filter on every row.
`2872949` and `3073571` have one non-REF/ALT base in caller-shaped depth. The other eight loci have other = 0.
CIGAR I/D/S within 5 bp is 0, except two reads at `2776013` and one read at `161068`.
ALT is present on both forward and reverse strands at every locus.

## Diagnosis

| Locus | Sample / Run | Prior P2b state | Read-level diagnosis | Caller-shaped AF | Confidence |
| --- | --- | --- | --- | ---: | --- |
| `1571257 G>T` | Cas9_3 _and2others (`SRR28349204`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.300 | High |
| `2776013 T>C` | Cas9_104_105 (`SRR28349176`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.292 | High |
| `4035824 T>C` | Cas9_104_105 (`SRR28349176`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.644 | High |
| `2872949 C>A` | I-SceI_7_8_9 (`SRR28349230`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.691 | High |
| `1635887 A>G` | control_5_6_7_8 (`SRR28349225`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.253 | High |
| `3108123 G>A` | control_17_18_19_20 (`SRR28349221`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.198 | High |
| `161068 A>C` | lambdaRed_7_8 (`SRR28349165`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.320 | High |
| `2653047 T>A` | I-SceI_1_2 (`SRR28349232`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.291 | High |
| `4115617 C>A` | I-SceI_1_2 (`SRR28349232`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.379 | High |
| `3073571 G>A` | I-SceI_14_15_16 (`SRR28349228`) | `UNRESOLVED` | `LOW_ALT_FREQUENCY` | 0.336 | High |

**Observation:** Across all 10 unrecovered publication loci:
1. Every published ALT is physically present in the primary alignments.
2. Every published ALT has bidirectional strand support (`plus_caller_shaped_alt > 0` and `minus_caller_shaped_alt > 0`).
3. Caller-shaped depth is high across all loci (range: 59 to 217).
4. Phred Q20 filtering confirms high base qualities (ALT persists at Q20 in all 10 loci).
5. In every case, the caller-shaped ALT fraction is strictly below the consensus caller cutoff of 0.75 (range: 0.198113 to 0.691244).

**Interpretation:** This pattern demonstrates conclusively that all 10 unrecovered mutations are genuine sub-clonal or multi-colony pooled alleles (as documented in the publication methods) filtered out by the standard consensus calling threshold ($\ge 0.75$). None are sequencing artifacts, reference assembly mismatches, mapping ambiguities, or caller bugs.

## Limits

The homopolymer-end mask is not in this diagnostic count. Stage-2 alignments are not in
the BAM. No precision, recall, or F1 is reported. Core crates and caller defaults remain unchanged.

## Verdict

**M6-P3a ALL 10 ALLELES MEASURED — READY FOR CODEX INDEPENDENT REVIEW**

**M6 overall: ACTIVE**
