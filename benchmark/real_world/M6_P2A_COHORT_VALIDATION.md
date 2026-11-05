# M6-P2a eligible paired cohort

Status: **M6-P2a EXECUTION COMPLETE — SCIENTIFIC DISCREPANCIES REQUIRE REVIEW**.
M6 remains **ACTIVE**. M6-P3 has not started. Tier B remains **NOT_EVALUABLE**.

This is reference-relative publication-event recovery for the frozen 25-SampleName
development subset of PRJNA1088182. It is not a Cas9 causality, specificity, or
F1 benchmark, and it is not a result for the other 64 SampleNames in the BioProject.

## Cohort

Source tables: `prjna1088182_cohort.tsv`, `prjna1088182_relationships.tsv`,
and `prjna1088182_sra_runinfo.csv`. The execution accession is
`relationships.tsv::execution_reference`. `m6_reference_policy.py --check` passed.
`avgLength` was not a selection rule.

| Class | Count | Rule |
| --- | ---: | --- |
| Development SampleNames | 25 | Not the 89-SampleName BioProject |
| SRA runs considered | 50 | Both public runs of each SampleName |
| `ELIGIBLE_PAIRED` | 25 | `SRA1824580`, `spots_with_mates` equal to spots, reported layout `PAIRED` |
| `BLOCKED_READ_STRUCTURE` | 25 | `SRA1825115`, this run's `spots_with_mates=0` |

`SRR28357593` and `SRR28357592` remain `SINGLE_BIOLOGICAL_READ_PER_SPOT` from the
P1a `sra-stat` inspection. The other 23 `SRA1825115` runs are
`RUNINFO_SPOTS_WITH_MATES_ZERO`. Their biological-read layout was not inferred
from those two inspections. They were not downloaded again and are not technical
replicates.

Machine-readable cohort:
`manifests/prjna1088182_m6_p2a_execution_cohort.tsv`.

## References

| Strain | Execution accession | Length | Sequence SHA256 | Local file |
| --- | --- | ---: | --- | --- |
| MG1655 | `NC_000913.3` (`U00096.3` equivalent) | 4,641,652 | `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0` | `data/references/MG1655.gbk` |
| BW25113 | `CP009273.1` | 4,631,469 | `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce` | `data/references/m6_identity/CP009273.1.fasta` |

Both hashes were checked from the local files before staging. No run mixed
strain and reference. BW25113 did not use `NC_000913.3`.

The BW25113 file is the sequence-verified FASTA. It has no GenBank
`repeat_region` features. MG1655 used the GenBank file, so repeat annotations
were loaded. Every completed BW25113 run has `MOB=0`. Every MG1655 run has
2 or 3 MOB records. That difference is an input-format limitation. It is not a
retune of the caller.

## Staging

SRA Toolkit 3.4.1, `prefetch`, `vdb-validate`, `sra-stat -x`, `fasterq-dump
--split-files`, `pigz -n` 2.8. Slurm array `2899261` on `qcpu_18i`, 4 CPU, 16G,
6 h, at most 6 concurrent tasks. The two P1a paired runs were not downloaded
again.

All 23 new runs finished exit 0. Each new log contains `sra-stat` with
`spot_count_mates` equal to spots. R1 and R2 counts match those spots.
Checksums are in `manifests/prjna1088182_m6_p2a_fastq_provenance.tsv`.

Staging wall clock was 16:02:10 to 17:18:41 on 2026-10-08. Per-task elapsed
ranged from 00:02:09 (`SRR28349177`) to 00:47:03 (`SRR28349227`). The slow
tasks spent almost all of that time in HTTPS `prefetch`:

| Run | Spots | Download window (UTC) | Task elapsed |
| --- | ---: | --- | --- |
| `SRR28349227` | 3,038,119 | 08:31:42–09:16:59 | 00:47:03 |
| `SRR28349204` | 4,051,149 | 08:02:14–08:29:20 | 00:29:47 |
| `SRR28349231` | 2,617,115 | 08:27:09–08:47:40 | 00:21:56 |

`fasterq-dump`, validation, and compression after a finished download were
about one to three minutes. FASTQ conversion was not the slow stage. Evidence
waited on `afterany:2899261`, so the last slow download held the cohort.

## Execution

Command: `prokadiff evidence` with the strain-specific reference, the paired
FASTQ, and `--threads 8`. No `--starter`, no intended-edit table, no threshold
change, no breseq.

Slurm array `2899282`, `qcpu_18i`, 8 CPU, 32G, 16 h, at most 4 concurrent
tasks. Calendar span 2026-10-08 17:18:41 to 2026-10-09 05:10:00. All 23 tasks
completed exit 0. The binary SHA256 for every executed run, including the two
reused P1b runs, is
`02077c761b1c9bd83cf1b626a3e2e955d8d2a1c63abb3e35852421e0dbc20bfd`.
A later `cargo build --release` relinked a different artifact and was reverted
before this array. Array tasks did not run cargo.

`SRR28349244` and `SRR28349243` are `REUSED_VERIFIED_P1B`. Their telemetry
commands are `prokadiff evidence` against `MG1655.gbk` with the same FASTQ and
without `--starter`. P1b commit `7025051` and this array's commit `13fd4bc`
have no `crates/` difference.

P1b reference point, unchanged: about 1.6 million spots, evidence wall 1:35:20
and 1:33:46, peak RSS about 2.9 GB. New-run `/usr/bin/time` walls track sacct
elapsed within about a minute, so FASTQ re-hashing is not the evidence cost.
The slowest evidence tasks are the largest lambda Red libraries:

| Run | Sample | Spots | sacct elapsed | telemetry wall | telemetry RSS |
| --- | --- | ---: | --- | --- | ---: |
| `SRR28349166` | lambdaRed_5_6 | 4,362,744 | 03:27:14 | 3:26:43 | 9,642,832K |
| `SRR28349164` | lambdaRed_9_10 | 4,093,992 | 03:11:21 | 3:10:57 | 8,913,572K |
| `SRR28349165` | lambdaRed_7_8 | 4,267,402 | 02:57:13 | 2:56:38 | 9,107,604K |
| `SRR28349168` | lambdaRed_1_2 | 4,960,886 | 02:40:33 | 2:40:04 | 10,506,824K |

Peak RSS stayed under the 32G request. The highest was `SRR28349168` at
10,650,888K in sacct. Relative to P1b, about 5.0 million spots took 2:40 and
about 10 GB, against 1.6 million spots at 1:35 and 2.9 GB. Wall time did not
grow only with spot count: `SRR28349231` has 2,617,115 spots and took 2:51:56.
The overnight calendar is the four-wide array, not a hung task. The longest
single evidence task was 3:27.

## Publication-event recovery

Truth is the frozen coordinate audit, not the GD files. Duplicate
`2872949 C>A` is one allele. `lambdaRed_3_4` `3112908 +AAT` stays
`NOT_EVALUABLE` and is outside the exact-recovery denominator. `NONE` rows
are additional-call inventories, not false positives.

| Denominator | Value |
| --- | ---: |
| Development cohort SampleNames | 25 |
| Paired-run eligible | 25 |
| Staged and execution-ready | 25 |
| Successful executions | 25 (23 new + 2 reused) |
| Unique eligible publication alleles | 12 |
| Unique alleles recovered in at least one run | 2 |
| Per-run eligible observations | 12 |
| Exact recoveries | 2 |
| Equivalent recoveries | 0 |
| Representation differences | 0 |
| Non-recoveries | 10 |
| Excluded truth events | 1 (`+AAT`) |

Exact recoveries:

- `Cas9_1` / `SRR28349244`: `301690 C>T`, `RECOVERED_EXACT`, GD record id 3, DV1 EventId `NOT_EVALUABLE`
- `Cas9_103` / `SRR28349177`: `3536100 G>T`, `RECOVERED_EXACT`

The other 10 eligible alleles are `NOT_RECOVERED` with diagnostic `UNKNOWN`.
This reconciliation did not reopen BAM pileups, so coverage, quality, mapping,
and caller-threshold classes were not assigned. Thresholds were not changed.

`Cas9_2` remains `NONE` and is not an FPR denominator. `301690` is absent there,
as in P1b.

## Structural pattern for review

JC burden is uneven. `Cas9_110_111` / `SRR28349173` has 505 JC records
(support reads: minimum 6, median 16, maximum 2742; 449 records have at least
10). Its repeat-diagnostic table has 34,450 data rows. `Cas9_104_105` has 140
JC records. P1b `Cas9_1` had 49. These counts are reported for review. They
are not an automatic rejection and are not publication JC truth.

MC records and repeat-diagnostic rows are not mutation identities.

## DV1 and M4

Evidence GD numeric ids are not DV1 EventIds. No legitimate parent comparator
was available, so DV1 EventId checks are `NOT_EVALUABLE`. The seven M4
differential products are `NOT_EVALUABLE_BY_TIER_A`.

## Limits

No biological parent. No breseq 0.40.2 oracle. Publication rows stay
`PUBLICATION_REPORTED`. Results are not genome-wide sensitivity and do not
extend to the 89-SampleName BioProject. `VALIDATION_REGISTRY.tsv` was not
updated.
