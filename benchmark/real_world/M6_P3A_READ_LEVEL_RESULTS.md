# M6-P3a first-wave read-level results

Status: **M6-P3a FIRST WAVE COMPLETE — READY FOR SCIENTIFIC REVIEW**.
M6 remains **ACTIVE**. The other six unrecovered alleles were not measured.
P2a still records 2/12 `RECOVERED_EXACT` and 10 `NOT_RECOVERED`. Those tables
were not edited.

Counts below are diagnostic pileups on a stage-1 BAM. They are not internal
RA column totals. Phred ≥ 20 is a diagnostic filter. The caller does not apply
it. Caller-shaped depth keeps primary alignments with MAPQ ≥ 10 after the
XS≥AS rule and drops one aligned base at each end. It does not apply the
homopolymer-end mask.

Paired mates share a qname. Each mate is counted once. A second primary
record with the same qname and the same mate flag is skipped. Secondary and
supplementary records are not counted. In these four loci that skip count is 0.

## Input gate

Slurm array `2905873` on `qcpu_18i`, 8 CPU, 32 GB, 4 hours, at most two
concurrent tasks. All three tasks exited 0.

| Task | Run | Sample | Reference | SHA256 | Reads |
| --- | --- | --- | --- | --- | ---: |
| `2905873_0` | `SRR28349204` | Cas9_3 _and2others | `NC_000913.3` | `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0` | 4,051,149 |
| `2905873_1` | `SRR28349176` | Cas9_104_105 | `NC_000913.3` | same MG1655 hash | 3,433,823 |
| `2905873_2` | `SRR28349230` | I-SceI_7_8_9 | `CP009273.1` | `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce` | 3,230,652 |

Gzip checksums matched the P2a provenance before alignment. Bowtie2 2.5.4.
Command: `bowtie2 --local --no-unal -p 8` with `-1`/`-2` and unmatched FASTQ.
No stage-2 seed and no junction parameters.

BAM sorting used samtools 1.24 (`samtools sort` and `samtools index`).
Production writes the stage-1 SAM through an in-memory noodles sort and does
not write a `.bai`. The alignment records, flags, CIGAR, sequence, quality,
and AS/XS tags are the Bowtie2 output. The sort implementation is the
documented difference.

Elapsed time and batch MaxRSS: 00:09:19 and 186,876K; 00:07:26 and 179,176K;
00:07:56 and 182,776K. Batch CPU time was 01:10:35, 00:54:44, and 00:58:52.
No task approached the 32 GB request. `aligned.bam` and `aligned.bam.bai`
remain under `benchmark/real_world/results/prjna1088182/m6_p3a/<run>/`.

## Measurements

Caller-shaped ALT fraction is ALT divided by caller-shaped depth.

| Locus | Raw depth | Raw ALT | MAPQ≥10 depth | Caller-shaped | ALT fraction | Q20 ALT / Q20 depth | ALT plus / minus |
| --- | ---: | ---: | ---: | --- | ---: | --- | --- |
| `1571257 G>T` | 239 | 55 | 81 | 24/80 | 0.300 | 54/236 | 9/15 |
| `2776013 T>C` | 148 | 42 | 148 | 42/144 | 0.292 | 42/126 | 18/24 |
| `4035824 T>C` | 155 | 62 | 62 | 38/59 | 0.644 | 62/152 | 21/17 |
| `2872949 C>A` | 217 | 150 | 217 | 150/217 | 0.691 | 137/202 | 96/54 |

REF + ALT + other equals the depth of that same filter for every row.
`2872949` has one non-REF/ALT base in the caller-shaped set. The other three
loci have other = 0. CIGAR I/D/S within 5 bp is 0, except two reads at
`2776013`. ALT is present on both strands at every locus.

## Diagnosis

| Locus | Prior state | This measurement | Confidence |
| --- | --- | --- | --- |
| `1571257 G>T` | `NOT_RECOVERED` / `UNKNOWN` / P2b `UNRESOLVED` | `LOW_ALT_FREQUENCY` | High for this pileup |
| `2776013 T>C` | same | `LOW_ALT_FREQUENCY` | High for this pileup |
| `4035824 T>C` | same | `LOW_ALT_FREQUENCY` | High for this pileup; 0.644 is the closest of the four |
| `2872949 C>A` | same | `LOW_ALT_FREQUENCY` | High for this pileup; 0.691 is below 0.75 on the raw, MAPQ, and Q20 fractions as well |

Observation: each published ALT is present, on both strands, with caller-shaped
depth from 59 to 217, and the ALT fraction is below 0.75. Q20 does not remove
the ALT. Absence and insufficient coverage do not fit these counts.

Interpretation: that pattern is consistent with a consensus rule that requires
frequency 0.75, which is the current RA default. It does not show that the
RA engine applied that rule to these reads, and it is not evidence of a caller
defect. `4035824` would need about seven more caller-shaped ALT bases at the
same depth to reach 0.75. `2872949` would need about thirteen. Mixture
percentages were not used.

## Limits

The homopolymer-end mask is not in this count. Stage-2 alignments are not in
the BAM. The six priority-2 alleles are still unmeasured. No precision, recall,
or F1 is reported.

## Verdict

**M6-P3a FIRST WAVE COMPLETE — READY FOR SCIENTIFIC REVIEW**

**M6 overall: ACTIVE**
