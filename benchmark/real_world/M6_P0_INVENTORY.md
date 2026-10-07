# M6-P0 inventory — PRJNA1088182

Status: **M6-P0 FROZEN FOR TIER A** (2026-10-07). M6 overall is **ACTIVE**. M6-P1 is **NOT STARTED**. Tier B is **NOT_EVALUABLE**.

P0.1 tables, which do not replace the cohort manifest:

- `manifests/prjna1088182_relationships.tsv`
- `manifests/prjna1088182_truth_coordinate_audit.tsv`

This file audits versioned inputs and local metadata. It does not execute
M6-P1. Canonical machine-readable inputs stay:

- `manifests/prjna1088182_cohort.tsv` (the 25-sample declaration; `reference` is the execution accession and must match `prjna1088182_relationships.tsv`)
- `manifests/prjna1088182_sra_runinfo.csv` (SRA run table, 178 runs)
- `truth/prjna1088182_publication_truth.tsv` (publication-reported rows)

`manifests/prjna1088182.tsv` is a conflicting 6-row legacy file. It is not
the M6 cohort. Its run accessions are absent from the versioned runinfo.

Local-only, gitignored inputs used for this audit:

- `benchmark/real_world/metadata/media-1.xlsx` (supplementary table)
- `data/references/MG1655.gbk` and `MG1655.fna`
- `benchmark/real_world/results/prjna1088182/` (not ProkaDiff calls)

## A. Scope

| Milestone | Dataset | This audit |
| --- | --- | --- |
| M6 | PRJNA1088182 (Widney & Copley cohort) | In scope |
| M7 | PRJNA884016 assembly truth | Out of scope except the handoff note below |

M7 handoff only: `VALIDATION_REGISTRY.tsv` row `prjna884016 / SY22A1` is
`PENDING` with all validation fields `false`. Assembly GD/MOB/JC files are
versioned. No ProkaDiff structural comparison is versioned. M7 is not started.

Publication source retrieved 2026-10-07:

- Title: *CRISPR-Cas9-assisted genome editing in E. coli elevates the frequency of unintended mutations*
- bioRxiv DOI: `10.1101/2024.03.19.584922` (`https://www.biorxiv.org/content/10.1101/2024.03.19.584922v1`)
- A PMC HTML record `PMC10983943` was also retrieved. This audit did not
  re-extract coordinates from the HTML. The versioned truth TSV was checked
  against local `media-1.xlsx`, not against a re-keyed journal table.

## B. Cohort_25 selection

`build_prjna1088182_cohort.py` hard-codes 25 `SampleName` values (10 Cas9,
5 control, 5 lambda Red, 5 I-SceI). The runinfo has 89 sample names and 178
runs (every name has two runs and two experiments). The builder keeps the
run with `avgLength > 200` and drops the shorter companion. No sample name
has two runs above that cutoff, so the kept run is unique.

Why these 25: a curated mix of names the script comments call “clean” or
“mutated”, plus group coverage. It is not a random or complete BioProject
sample. Selection bias is therefore present by construction: Cas9 is 10/25,
22/25 are labeled `pooled` by a three-name exception (`Cas9_1`, `Cas9_2`,
`Cas9_103` only), and 64 other sample names are absent.

Script comments disagree with the truth TSV and must not be used as truth:

- Comment says `Cas9_103` is clean. Truth row is `3536100 G>T`.
- Comment says `I-SceI_14_15_16` mutation `4068404` and `I-SceI_17_18_19`
  mutation `3073571`. Truth has `3073571 G>A` on `I-SceI_14_15_16` and
  `NONE` on `I-SceI_17_18_19`.

`single_clone` versus `pooled` is a name rule in the builder, not a
clonality assay. Spreadsheet `gDNA` mixture rows show multi-colony libraries
for most names. That supports “not a single colony” for those libraries. It
does not prove the three `single_clone` labels.

## C. Accession join (25 kept runs)

Every kept `sample_id` occurs once in the runinfo. `SampleName` matches the
cohort name. BioProject is `PRJNA1088182`. Layout is `PAIRED`. Platform is
`ILLUMINA` / `Illumina NovaSeq X Plus`. Taxon strings are not uniform:

- 10 Cas9 runs: `Escherichia coli str. K-12 substr. MG1655`
- 15 control / lambda Red / I-SceI runs: `Escherichia coli BW25113`

| sample_id | sample_name | group | BioSample | Experiment | companion run (not in cohort) | SRA organism |
| --- | --- | --- | --- | --- | --- | --- |
| SRR28349244 | Cas9_1 | Cas9 | SAMN40463616 | SRX23955191 | SRR28357593 | MG1655 |
| SRR28349243 | Cas9_2 | Cas9 | SAMN40463617 | SRX23955192 | SRR28357592 | MG1655 |
| SRR28349177 | Cas9_103 | Cas9 | SAMN40463660 | SRX23955258 | SRR28357554 | MG1655 |
| SRR28349204 | Cas9_3 _and2others | Cas9 | SAMN40463618 | SRX23955231 | SRR28357581 | MG1655 |
| SRR28349178 | Cas9_101_102 | Cas9 | SAMN40463659 | SRX23955257 | SRR28357555 | MG1655 |
| SRR28349176 | Cas9_104_105 | Cas9 | SAMN40463661 | SRX23955259 | SRR28357553 | MG1655 |
| SRR28349175 | Cas9_106_107 | Cas9 | SAMN40463662 | SRX23955260 | SRR28357552 | MG1655 |
| SRR28349174 | Cas9_108_109 | Cas9 | SAMN40463663 | SRX23955261 | SRR28357551 | MG1655 |
| SRR28349173 | Cas9_110_111 | Cas9 | SAMN40463664 | SRX23955262 | SRR28357550 | MG1655 |
| SRR28349172 | Cas9_112_113 | Cas9 | SAMN40463665 | SRX23955263 | SRR28357549 | MG1655 |
| SRR28349226 | control_1_2_3_4 | control | SAMN40463692 | SRX23955209 | SRR28357517 | BW25113 |
| SRR28349225 | control_5_6_7_8 | control | SAMN40463693 | SRX23955210 | SRR28357516 | BW25113 |
| SRR28349224 | control_9_10_11_12 | control | SAMN40463694 | SRX23955211 | SRR28357515 | BW25113 |
| SRR28349223 | control_13_14_15_16 | control | SAMN40463695 | SRX23955212 | SRR28357514 | BW25113 |
| SRR28349221 | control_17_18_19_20 | control | SAMN40463696 | SRX23955214 | SRR28357512 | BW25113 |
| SRR28349168 | lambdaRed_1_2 | lambdaRed | SAMN40463668 | SRX23955267 | SRR28357545 | BW25113 |
| SRR28349167 | lambdaRed_3_4 | lambdaRed | SAMN40463669 | SRX23955268 | SRR28357544 | BW25113 |
| SRR28349166 | lambdaRed_5_6 | lambdaRed | SAMN40463670 | SRX23955269 | SRR28357543 | BW25113 |
| SRR28349165 | lambdaRed_7_8 | lambdaRed | SAMN40463671 | SRX23955270 | SRR28357542 | BW25113 |
| SRR28349164 | lambdaRed_9_10 | lambdaRed | SAMN40463672 | SRX23955271 | SRR28357541 | BW25113 |
| SRR28349232 | I-SceI_1_2 | I-SceI | SAMN40463686 | SRX23955203 | SRR28357523 | BW25113 |
| SRR28349231 | I-SceI_3_4_5_6 | I-SceI | SAMN40463687 | SRX23955204 | SRR28357522 | BW25113 |
| SRR28349230 | I-SceI_7_8_9 | I-SceI | SAMN40463688 | SRX23955205 | SRR28357521 | BW25113 |
| SRR28349228 | I-SceI_14_15_16 | I-SceI | SAMN40463690 | SRX23955207 | SRR28357519 | BW25113 |
| SRR28349227 | I-SceI_17_18_19 | I-SceI | SAMN40463691 | SRX23955208 | SRR28357518 | BW25113 |

Companion runs share `SampleName` and `BioSample` with the kept run and have
a different `Experiment`. They are not represented in the cohort. Sharing a
BioSample does not by itself prove they are technical replicates of one
library.

FASTQ paths in the cohort (`data/prjna1088182/<RUN>_1.fastq.gz` and `_2`) are
**MISSING**. `data/prjna1088182/` exists and is empty. Size and checksum:
not applicable until files exist. Missing is not zero.

SRA `size_MB` / `spots` in the runinfo are archive metadata for the kept
run, not local FASTQ measurements.

## D. Reference and biological parent

Local reference files (gitignored `data/references/`):

| File | Identity in file header | Bases | SHA256 |
| --- | --- | --- | --- |
| `MG1655.gbk` | `ACCESSION NC_000913` / `VERSION NC_000913.3` | locus length 4641652 | `313752bfe10d4aaf04abdb541a1309bce9a6b66f059e75fe215898e736a5a8a5` |
| `MG1655.fna` | `>NC_000913.3` | 4641652 | `6b195feda4c66140f6762742eb8b30c2652f02b45878b174f5b00ef85ecc95d7` |

The cohort `reference` column is now the execution accession (`NC_000913.3`
or `CP009273.1`), matching `prjna1088182_relationships.tsv`. The local file
`data/references/MG1655.gbk` is the sequence for `NC_000913.3` only. It is not
the BW25113 execution reference.

`starter_id` is `SRR28349226` for all 25 rows, including that control itself.
The builder comment calls this control the “Reference starter”. Spreadsheet
parent fields do not support that collapse.

Documented parents from `media-1.xlsx` (local, gitignored), and whether that
parent BioSample is in the versioned runinfo:

| Cohort library | Spreadsheet parent | Parent BioSample or accession | In runinfo |
| --- | --- | --- | --- |
| Cas9_1 | MG1655* fbaB syn intm | SAMN40533628 | no |
| Cas9_2 | MG1655* pgi352 syn intm | SAMN40533629 | no |
| Cas9_103, Cas9_101_102 | MG1655* acs328::mCherry intm C | SAMN40533622 | no |
| Cas9_3 _and2others | MG1655* cyaA499::mCherry intm | SAMN40533627 | no |
| Cas9_104_105 through Cas9_112_113 | MG1655* acs610::mCherry intm | SAMN40533623 | no |
| five control pools | BW25113* | SAMN40533614 | no |
| five lambda Red pools | BW25113 | GenBank `CP009273.1` (not a BioSample in this table) | no |
| I-SceI_1_2 | BW25113* yifE Q100* | SAMN40533615 | no |
| I-SceI_3_4_5_6 | BW25113* gapA* ΔybhA Δpgl | SAMN40533617 | no |
| I-SceI_7_8_9 | BW25113* gapA* | SAMN40533616 | no |
| I-SceI_14_15_16 | BW25113* | SAMN40533614 | no |
| I-SceI_17_18_19 | BW25113* gapA* ΔybhA Δpgl rpoC* | SAMN40533618 | no |

No audited parent accession is a row of `prjna1088182_sra_runinfo.csv`.
`SRR28349226` is a four-colony control pool (`C1-1`…`C1-4`, BioSample
`SAMN40463692`), not those parental records. Biological parenthood is
therefore **documented in the supplement and absent from the run table**. It
is not established by the manifest `starter_id`.

`intended_tsv` is `.` for every cohort row. Spacers exist only where the
builder split a Cas9 protospacer string. That is guide metadata, not an
intended-edit allele table. `editor=dsb` on lambda Red / I-SceI / control
rows is a builder default, not a verified statement that lambda Red is a DSB
editor.

## E. Truth versus oracle

Truth TSV: 28 rows, 25 sample IDs, all inside the cohort.

| Class | Count | Meaning |
| --- | ---: | --- |
| `NONE` | 14 | Spreadsheet recorded no unintended mutation for that library |
| non-`NONE` | 14 | Publication-reported alleles |

Non-`NONE` types in the file: synonymous 4, nonsynonymous 3, nonsense 2,
ribosomal RNA 1, snp_intergenic 1, snp_nonsynonymous 1, small indel 1,
intergenic 1. There is no SNP/DEL/JC/MOB GD type column.

Duplicate key: `SRR28349230` / `2872949` / `C` / `A` appears twice
(`I-SceI_7_8_9`). The builder comment says two identical colonies. Count it
once. Do not treat the duplicate as two events.

Provenance class for every truth row: **PUBLICATION_REPORTED**.

| Source class | Status |
| --- | --- |
| PUBLICATION_REPORTED | PRESENT_VERIFIED against local `media-1.xlsx` for the audited rows; the xlsx itself is gitignored |
| INDEPENDENT_ORACLE | MISSING. No breseq/gdtools GD is versioned for this cohort |
| ASSEMBLY_DERIVED | NOT_APPLICABLE for M6 (that is M7 / PRJNA884016) |
| MANUALLY_CURATED | NOT the class of this TSV |
| ProkaDiff calls used as truth | Forbidden. None are versioned |

`mutation_recovery.tsv` under ignored `results/prjna1088182/` labels rows
`detected_by_prokadiff=CONFIRMED` and `RECOVERED` without a ProkaDiff GD or
BAM. Those files annotate the publication table. They are not an independent
oracle and not a validation measurement. The registry note already says this.
This audit does not delete them and does not promote them.

Coordinate convention recorded in the TSV: integer `pos`, `ref_allele`,
`alt_allele`, no contig, no strand, no zero/one-based flag. `evaluate_prjna1088182.py`
uses genome length `4641652`, equal to local `NC_000913.3`. Exact recovery
against that reference for BW25113-labeled runs is **NOT_EVALUABLE** until
the coordinate backbone is explicitly chosen. A missing denominator is not 0%.

## F. Calls and comparisons

| Artifact | Class |
| --- | --- |
| Cohort manifest and runinfo | PRESENT_VERIFIED as tables; FASTQ bytes MISSING |
| Publication truth TSV | PRESENT_VERIFIED as a table of reported alleles |
| Local MG1655 `NC_000913.3` | PRESENT_VERIFIED as files; fitness as the analysis reference for BW25113 runs is UNFROZEN |
| Parental SRA/FASTQ | MISSING from the versioned runinfo |
| ProkaDiff `.gd` / M4 products | MISSING |
| breseq oracle | MISSING |
| Same-job wall/RSS comparison | MISSING |
| Ignored `group_summary.tsv` and `REPORT.md` | PRESENT_UNVERIFIED as publication annotation only |

## G. Registry claim audit

`benchmark/real_world/VALIDATION_REGISTRY.tsv`

| dataset | sample | truth_level | variant | intended | mob | jc | report | status | Justified? |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| prjna1088182 | cohort_25 | PUBLICATION_REPORTED | false | false | false | false | false | PENDING | Yes. Manifest and publication truth exist. Calls and comparisons do not. Do not set any field `true`. |

No registry edit is proposed. M5 BL21 rows and the M7 `prjna884016` row were
not changed.

## H. Gap report

Blockers before any FASTQ job:

1. Parental libraries named in `media-1.xlsx` are not in the versioned runinfo. `starter_id=SRR28349226` is not those parents.
2. Fifteen cohort runs are SRA-labeled BW25113. The only local reference is MG1655 `NC_000913.3`.
3. FASTQ files are absent. Checksums do not exist.
4. No independent call oracle is pinned to these samples. `testdata/VERSIONS.txt` pins breseq 0.40.2 for the project environment; that pin is not an executed oracle here.
5. Publication positions have no contig column. Coordinate eligibility is unresolved for BW25113-labeled runs.
6. `intended_tsv` is empty, so M2 status must stay **NOT_EVALUABLE**, not `MISSING`.
7. Legacy `prjna1088182.tsv` (`SRR28385261`–`SRR28385265`, `SRR28385285`) is unsupported by the runinfo.
8. Companion shorter runs are excluded by an `avgLength > 200` rule that is not yet an accepted scientific choice.
9. `editor=dsb` for lambda Red is not a verified editor class.
10. Ignored result TSVs must not be cited as ProkaDiff recovery.

## I. M6-P1 readiness

| Gate | State |
| --- | --- |
| Cohort membership described | Yes, as a curated 25, not as the whole BioProject |
| Run identity for the 25 kept accessions | Yes, against versioned runinfo |
| Parent/control semantics | No. Documented parents are missing from the run table |
| Reference identity | Local file verified; analysis role for BW25113 runs unfrozen |
| Truth source attribution | Yes: PUBLICATION_REPORTED |
| Independent oracle | Not selected for execution |
| Event matching | Project ±5 bp rules stay as-is; dataset-specific tolerances forbidden. Exact publication-integer match is UNFROZEN until coordinates are accepted |
| Metrics | See contract. Thresholds that lack a denominator stay NOT_EVALUABLE |
| Registry | Correctly PENDING |
| Output location | Future jobs must use a new gitignored results prefix, not overwrite `results/prjna1088182/` annotation files without labeling them |
| Resources | Uncertain proposal only: `qcpu_18i`, 16 CPU, 32G, 4 h per smoke library. BL21 peer jobs used 16 CPU and about 16 GiB RSS for larger private WGS. These SRA `size_MB` values are smaller. This is not a measurement |
| Stop conditions | Stop if parent is required but absent, if the reference accession is not the one named in the job record, if FASTQ checksums are missing, or if a product would emit M2 `COMPLETE`/`MISSING` without an intended table |

## L. P0.1.1 accession equivalence and run policy

NCBI assembly report `GCF_000005845.2_ASM584v2` line:

`U00096.3`, Relationship `=`, `NC_000913.3`, length 4,641,652.

Chromosome sequence SHA256
`e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0`
is shared by NCBI `U00096.3`, NCBI `NC_000913.3`, and local
`data/references/MG1655.gbk`. The 10 Cas9 samples are
`VERIFIED_EQUIVALENT_ACCESSION`: publication reference `U00096.3`, execution
reference `NC_000913.3`.

`CP009273.1` is a different chromosome: length 4,631,469, sequence SHA256
`3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce`.
Every BW25113 substitution reference allele matches that sequence at the
1-based coordinate. Those rows stay `COORDINATE_VERIFIED` and eligible.
Execution reference is `CP009273.1`, not `NC_000913.3`. The `+AAT` row stays
`COORDINATE_UNRESOLVED` and is excluded from the denominator. It does not
remove the verified substitutions.

Both WGS runs of a SampleName are `INDEPENDENT_WGS_RUN`. Execution unit is
the SRA Run. Publication truth is attached to SampleName/BioSample. Runs are
not merged and not dropped by read length.

Downloaded FASTA files used for the hashes are under gitignored
`data/references/m6_identity/`. They are not cohort FASTQ.

## K. P0.1 identity resolution

Methods sentence used as the coordinate source (bioRxiv `10.1101/2024.03.19.584922`, Whole-genome sequencing and analysis): reads were analyzed with breseq 0.35.4 and aligned to MG1655 `U00096.3` or BW25113 `CP009273.1`. The project oracle pin remains breseq 0.40.2. The publication’s 0.35.4 calls are not reclassified as this project’s independent oracle.

### Run pairs

Each of the 25 `SampleName` values has two runs. Both are `WGS` / `RANDOM` / `GENOMIC` / `PAIRED` / `ILLUMINA` / `Illumina NovaSeq X Plus`, same `LibraryName` and BioSample, and different Experiment accessions.

| Submission | Release | avgLength | spots_with_mates | Class |
| --- | --- | --- | --- | --- |
| `SRA1824580` | 2024-03-14 | about 292–298 | equal to spots | `INDEPENDENT_WGS_RUN` |
| `SRA1825115` | 2024-03-18 | about 143–148 | 0 in runinfo | `INDEPENDENT_WGS_RUN` |

`avgLength > 200` is not a selection rule. Neither run is `PRIMARY_WGS`.
They are not technical replicates. Future execution runs each SRA Run
separately and keeps Experiment and Submission on the record.

### Strain and reference

| Group | SRA organism | Publication reference | Execution reference | Assignment |
| --- | --- | --- | --- | --- |
| 10 Cas9 | MG1655, taxid 511145 | `U00096.3` | `NC_000913.3` | `VERIFIED_EQUIVALENT_ACCESSION` |
| 15 control / lambda Red / I-SceI | BW25113, taxid 679895 | `CP009273.1` | `CP009273.1` | `VERIFIED_CP009273.1` |

The two MG1655 accession strings are one chromosome. BW25113 jobs must not use the manifest’s MG1655 path.

### Coordinates

Fourteen non-`NONE` truth rows are in the coordinate audit. Thirteen substitution rows are `COORDINATE_VERIFIED` and eligible (12 unique alleles; `2872949 C>A` is duplicated and counts once). The lambda Red `+AAT` row is `COORDINATE_UNRESOLVED` and is excluded from the denominator. Script comments that disagree with `media-1.xlsx` are rejected and listed in the relationship table.

### Parents and cohort role

Every documented parent pointer is `NO_PUBLIC_RUN_FOUND`. `starter_id=SRR28349226` stays a computational comparator. `cohort_role=DEVELOPMENT_SUBSET`. `single_clone` / `pooled` in the manifest is the builder’s name rule; mixture columns in `media-1.xlsx` are the experimental record and are copied into `mixture_evidence`.

## J. Proposed Tier A smoke subset (do not run)

Contract-eligible. Not submitted. FASTQ staging is M6-P1.

| SampleName | Runs, each independent | Truth | Execution reference |
| --- | --- | --- | --- |
| Cas9_1 | `SRR28349244` and `SRR28357593` | `301690 C>T`, eligible | `NC_000913.3` |
| Cas9_2 | `SRR28349243` and `SRR28357592` | `NONE`, not an FPR control | `NC_000913.3` |

`Cas9_104_105` can be a later smoke sample. It is mixed gDNA. BW25113 smoke
uses `CP009273.1` only, and the `+AAT` row stays out of the denominator.
`SRR28349226` is not the biological parent.
