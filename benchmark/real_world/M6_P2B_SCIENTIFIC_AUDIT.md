# M6-P2b scientific discrepancy audit

Status: **M6-P2b REVIEW COMPLETE — READY FOR M6-P3 PLANNING**.
M6 remains **ACTIVE**. M6-P3 has not started. This audit did not rerun the
cohort, change caller thresholds, or edit `crates/`.

P2a remains **EXECUTION COMPLETE — SCIENTIFIC DISCREPANCIES REQUIRE REVIEW**.
Its tables were not rewritten. The prior `UNKNOWN` diagnostic is preserved in
`p2a_diagnostic` of `manifests/prjna1088182_m6_p2b_event_diagnostics.tsv`.

## A. P2a reconciliation integrity

The frozen coordinate audit has 13 eligible rows and 12 unique alleles after
dropping the repeated `I-SceI_7_8_9` `2872949 C>A` row. The excluded row is
still `lambdaRed_3_4` `3112908 +AAT`.

Independent checks against the local chromosomes and the evidence GD files:

| Check | Result |
| --- | --- |
| Unique eligible alleles | 12 |
| Execution reference base equals the published reference allele | 12/12 |
| `Cas9_1` `301690 C>T` | `SNP` GD id 3, `NC_000913.3`, new base `T`, `RECOVERED_EXACT` |
| `Cas9_103` `3536100 G>T` | `SNP` GD id 16, depth 154, frequency 1, `RECOVERED_EXACT` |
| Other 10 alleles | no SNP, SUB, INS, or DEL at that coordinate |
| Duplicate `2872949 C>A` | one diagnostic row |
| `+AAT` | `NOT_EVALUABLE`, not converted to a zero |
| GD numeric ids | recorded as GD ids; DV1 EventId stays `NOT_EVALUABLE` |

`Cas9_1` used the P1b GD. `Cas9_103` used the P2a GD. Both references match
the relationship table: MG1655 `NC_000913.3`, BW25113 `CP009273.1`.

## B. Ten unrecovered alleles

`aligned.bam` is absent for every run. `keep_bam` defaults to false, and
`run_sample` deletes `aligned.bam` and `work/` after writing the GD. Depth,
ALT count, allele fraction, base quality, mapping quality, and strand balance
are **NOT_MEASURED**. They are not estimated.

No MC interval covers any of the ten coordinates, so the GD does not record
those sites as missing-coverage spans. No exact-position GD mutation of
another class is present, so this is not a representation difference in the
retained output.

| Sample | Run | Event | Reference base | Mixture | P2a state | P2b diagnosis |
| --- | --- | --- | --- | --- | --- | --- |
| Cas9_3 _and2others | SRR28349204 | `1571257 G>T` | G | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |
| Cas9_104_105 | SRR28349176 | `2776013 T>C` | T | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |
| Cas9_104_105 | SRR28349176 | `4035824 T>C` | T | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |
| control_5_6_7_8 | SRR28349225 | `1635887 A>G` | A | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |
| control_17_18_19_20 | SRR28349221 | `3108123 G>A` | G | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |
| lambdaRed_7_8 | SRR28349165 | `161068 A>C` | A | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |
| I-SceI_1_2 | SRR28349232 | `2653047 T>A` | T | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |
| I-SceI_1_2 | SRR28349232 | `4115617 C>A` | C | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |
| I-SceI_7_8_9 | SRR28349230 | `2872949 C>A` | C | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |
| I-SceI_14_15_16 | SRR28349228 | `3073571 G>A` | G | documented pooled | `NOT_RECOVERED` / `UNKNOWN` | `UNRESOLVED` |

Observation: the publication allele is not an exact GD record, and the MC
table does not mark the base as uncovered.

Interpretation, not a measurement: RA consensus emits a SNP only when coverage
is at least 5 and the majority base frequency is at least 0.75, unless strand
bias rejects it. Every documented pool component for these ten alleles has a
recorded percentage below 75% (25%, 33/67%, 35/65%, or 20/35/45%). A
colony-private allele could sit under that rule. That hypothesis is not
proven, because the reads were not re-counted.

Confidence in the cause is **LOW**. Confidence that the retained GD lacks the
exact allele is **HIGH**.

## C. Pooled versus mixture-not-recorded

`mixed=` in `prjna1088182_relationships.tsv` is the experimental field.
`n=` is a colony index, not a pool size (`n=103` is Cas9_103). SampleName and
`manifest_sample_type` were not used as clonality.

| Mixture field | Eligible alleles | Exact recovery |
| --- | ---: | ---: |
| `mixed=` lists components (`DOCUMENTED_POOLED`) | 10 | 0 |
| `mixed=N/A` (`UNRESOLVED_MIXTURE`) | 2 | 2 |

The two recoveries are `Cas9_1` and `Cas9_103`. Their mixture field is `N/A`.
That does not document a single clone. The builder label `single_clone` is
not treated as an assay.

All ten non-recoveries are in documented pools. With only two alleles on the
other side, this split is a description of the table. It is not a selected
biological rate.

## D. MG1655 versus BW25113

| Strain | Eligible alleles | Exact recovery | Non-recovery |
| --- | ---: | ---: | ---: |
| MG1655 `NC_000913.3` | 5 | 2 | 3 |
| BW25113 `CP009273.1` | 7 | 0 | 7 |

The three unrecovered MG1655 alleles are all on documented pools
(`Cas9_3 _and2others`, `Cas9_104_105`). The two MG1655 recoveries have
`mixed=N/A`. All seven BW25113 eligible alleles are on documented pools and
were not recovered. Strain and pooling are not separable in this subset.
Reference bases match on both chromosomes, so this is not a coordinate
mismatch between `U00096.3` / `NC_000913.3` and `CP009273.1`.

## E. BW25113 MOB capability

Read-only path:

- `load_repeats` calls `parse_genbank_repeats`.
- That parser reads GenBank `repeat_region` and `mobile_element` features.
- A FASTA file yields no features.
- If the loaded repeat list is empty, `run_sample` looks for a sibling `.gbk`
  next to the reference path. `data/references/m6_identity/` has no `.gbk`.
- MOB promotion requires `find_repeat_at` on both flanking junction sides.

MG1655 used `MG1655.gbk`. That file contains `mobile_element` features, and
every MG1655 run emitted 2 or 3 MOB records. Every BW25113 run emitted
`MOB=0`. The FASTA execution path could not observe an annotated MOB.
Classification: **INPUT_FORMAT_INCOMPATIBILITY**. This is not
`TRUE_ZERO_SUPPORTED` and is not a biological absence of insertion sequences.

## F. Cas9_110_111 JC review

`SRR28349173` has 505 JC records, 2 MOB, 4 SNP, 1 INS, 1 DEL, and 85 MC.
The repeat-diagnostic table has 34,450 data rows: 34,447 `CANDIDATE` and 3
`SUPPORTED_AMBIGUOUS`.

JC support reads: minimum 6, median 16, 449 at least 10, 25 at least 100,
maximum 2742. Tiers: 56 with 6–9 reads, 266 with 10–19, 158 with 20–99, 25
with 100 or more.

The 505 side-pairs are unique. None is an exact duplicate of another JC.
Median distance between the two sides is 731 bp, and 282 junctions are within
1 kb. Four records carry `mob_evidence=1`. Of 50 `mobile_element` spans parsed
from `MG1655.gbk`, 27 junctions have one side inside a span and 5 have both
sides inside. That is not a single repeat-anchored hotspot.

Observation: the count is many distinct, mostly local junctions with modest
support, plus a large candidate diagnostic table.

Interpretation: a high JC count is not, by itself, a false-positive call.
Thresholds were not changed.

For scale, other MG1655 JC totals in this cohort are 49 (`Cas9_1`), 61
(`Cas9_103`), and 140 (`Cas9_104_105`).

## G. Runtime bottlenecks

Existing Slurm and `/usr/bin/time` records only. No new job was submitted.

**Data transfer.** Staging array `2899261` spent its long tasks in HTTPS
`prefetch`. `SRR28349227` took 00:47:03, with download from 08:31:42 to
09:16:59 UTC for about 310 MB. `SRR28349204` took 00:29:47, mostly download.
After the download, `fasterq-dump` and `pigz` were about one to three minutes.
`SRR28349177` finished staging in 00:02:09. Evidence used
`afterany:2899261`, so the slowest download held the whole evidence array
until 17:18.

**Scheduler concurrency.** Evidence array `2899282` allowed 4 tasks at a time.
The 23 task elapsed times sum to 43.0 hours. Divided by 4, that is 10.8 hours.
The calendar span was 17:18 on 2026-10-08 to 05:10 on 2026-10-09, about 11.9
hours. The gap between 10.8 and 11.9 hours is queue packing, not a hung task.

**Scientific compute.** Telemetry CPU percent on 8 cores was 767–789%, so the
cores were busy. The longest evidence wall was `SRR28349166` at 3:26:43
(4,362,744 spots, peak RSS 9,642,832K). `SRR28349168` peaked at 10,506,824K
for 4,960,886 spots. sacct elapsed and telemetry wall differ by about a
minute, so FASTQ re-hashing is not the evidence cost. P1b at about 1.6
million spots took about 1:35 and 2.9 GB. Wall time grew with library size
and was not only proportional to spot count.

Possible later changes, not applied here: keep the evidence array from waiting
on every download, raise the concurrency cap if the partition is free, and
retain locus BAMs when a publication site needs pileup. None of those were
submitted.

## H. Confirmed scientific defects

None. Exact recovery still holds for the two alleles that P2a recovered.
Reference bases match. The ten missing calls are unexplained at read level.
BW25113 `MOB=0` follows from the FASTA input, not from a patched caller.

## I. Remaining unknowns

- Allele fraction and coverage at the ten unrecovered sites.
- Whether any of those alleles were present below the RA 0.75 frequency rule.
- Whether `Cas9_1` and `Cas9_103` are single clones. `mixed=N/A` does not
  establish that.
- Which of the 505 JC records would remain after an independent oracle.
- Whether BW25113 would emit MOB records from a GenBank file of `CP009273.1`
  with the same chromosome sequence. That file was not staged.

## J. Files and tests

Added:

- `benchmark/real_world/M6_P2B_SCIENTIFIC_AUDIT.md`
- `benchmark/real_world/manifests/prjna1088182_m6_p2b_event_diagnostics.tsv`
- `benchmark/real_world/scripts/audit_m6_p2b.py`
- `benchmark/real_world/scripts/test_m6_p2b_audit.py`

P2a result files were left in place. `crates/` and
`VALIDATION_REGISTRY.tsv` were not modified. `test_m6_p2b_audit.py`: 5 passed.

## K. Verdict

**M6-P2b REVIEW COMPLETE — READY FOR M6-P3 PLANNING**

**M6 overall: ACTIVE**
