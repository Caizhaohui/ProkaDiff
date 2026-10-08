# M6-P1a FASTQ staging

Status: **M6-P1a RESOLVED — two paired runs ready for initial P1b** (2026-10-08).

The two `SRA1825115` runs are valid public objects with one biological read
per spot. They stay in the inventory and are not eligible for the initial
paired-end smoke. M6-P1b has not been started.

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
spots (`MATCH`). Reported `LibraryLayout` remains `PAIRED` for all four runs.
The two `SRA1825115` runs also report `spots_with_mates=0`.

Job 2897413 (`qcpu_18i`, SRA Toolkit 3.4.1) inspected the local `.sra` objects.
`sra-stat -x` gives `spot_count_mates=0` and `base_count_bio` equal to
`base_count` for `SRR28357593` (1,899,109 spots) and `SRR28357592`
(1,730,336 spots). Each SEQUENCE row has one `READ_LEN` (150) and one
`READ_TYPE` (`SRA_READ_TYPE_BIOLOGICAL` plus `FORWARD` or `REVERSE`). That
orientation flag is not a second read. No `SRA_READ_TYPE_TECHNICAL` value
appears. `fastq-dump -X 8 --split-files` writes only `_1.fastq`; `--split-3`
writes one unsuffixed FASTQ. Those singleton files stay on disk.

The paired controls differ. `SRR28349244` has `spot_count_mates=1600722`.
Its first SEQUENCE rows are two biological reads, for example `151, 149` and
spot length 300. `fastq-dump -X 8 --split-files` writes `_1` and `_2`.
`SRR28349243` has the same two-biological-read pattern.

Classification: both `SRA1825115` runs are
`SINGLE_BIOLOGICAL_READ_PER_SPOT` and
`NOT_ELIGIBLE_INITIAL_PAIRED_SMOKE`. They are not download failures.
Initial paired P1b, not yet started, can use `SRR28349244` and `SRR28349243`
only. Tier-A truth rules are unchanged.

Machine-readable rows:
`manifests/prjna1088182_m6_p1a_fastq_provenance.tsv`.
FASTQ bytes stay under gitignored `data/prjna1088182/fastq/`.
