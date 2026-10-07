# M6-P1a FASTQ staging

Status: **M6-P1a NOT READY FOR M6-P1b** (2026-10-07).

M6-P0 Tier A remains frozen at `06f6259bd3892d2e9a518d5a5450e6af2322b229`.
M6 stays ACTIVE. Tier B stays NOT_EVALUABLE. This task did not run ProkaDiff,
Bowtie2, or breseq.

## Scope

Four independent smoke runs, execution reference `NC_000913.3`:

| SampleName | Run | Publication truth |
| --- | --- | --- |
| Cas9_1 | SRR28349244 | `301690 C>T`, eligible for a later recovery test |
| Cas9_1 | SRR28357593 | same truth; this run did not become paired FASTQ |
| Cas9_2 | SRR28349243 | `NONE` (not a mutation-free genome and not an FPR control) |
| Cas9_2 | SRR28357592 | same `NONE` limit; this run did not become paired FASTQ |

## Method

SRA Toolkit 3.4.1 was installed into gitignored
`data/tools/sra-tools` because the login image had no `prefetch`.
Staging used `prefetch`, `vdb-validate`, `fasterq-dump --split-files`, then
`pigz -n` (pigz 2.8) when both mates existed.

Slurm array `2894179` on `qcpu_18i`: 4 CPU, 16G, 2 hours. That size covers the
largest object, about 155 MB and 1.6 million spots, not a BL21 analysis.
Scratch stayed under `data/prjna1088182/tmp/<Run>` on the project filesystem.
Logs: `benchmark/real_world/results/m6_p1a/slurm_2894179_*.log`.

Local `data/references/MG1655.gbk` ORIGIN was checked before each download:
length 4,641,652, SHA256
`e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0`.

## Per-run result

| Run | Array task | In-job id | Wall | Exit | Peak RSS (batch) | Status |
| --- | --- | --- | --- | --- | --- | --- |
| SRR28349244 | 2894179_0 | 2894180 | 00:20:09 | 0 | 13052K recorded by sacct | READY_FOR_EXECUTION |
| SRR28357593 | 2894179_1 | 2894181 | 00:02:21 | 1 | 111160K | BLOCKED_PAIRING |
| SRR28349243 | 2894179_2 | 2894182 | 00:12:39 | 0 | 359844K | READY_FOR_EXECUTION |
| SRR28357592 | 2894179_3 | 2894179 | 00:01:37 | 1 | 111352K | BLOCKED_PAIRING |

The two `SRA1824580` runs produced R1 and R2 with read counts equal to SRA
spots (`MATCH`). The two `SRA1825115` runs have `spots_with_mates=0`.
`fasterq-dump --split-files` wrote only `<Run>_1.fastq`, and the written read
count equals the spot count. Those singleton files were kept. They are not
paired FASTQ, so those runs are not ready for P1b.

Machine-readable rows:
`manifests/prjna1088182_m6_p1a_fastq_provenance.tsv`.
FASTQ bytes stay under gitignored `data/prjna1088182/fastq/`.
