# M6 validation contract — PRJNA1088182

Status: **M6-P0 FROZEN FOR TIER A** (2026-10-07).

M6 overall remains **ACTIVE**. M6-P1 is **NOT STARTED**. Tier B is
**NOT_EVALUABLE**. FASTQ acquisition and checksums belong to M6-P1 staging.

Tier B is **NOT_EVALUABLE** for the current cohort and does not block Tier A.

Identity evidence:

- `M6_P0_INVENTORY.md`
- `manifests/prjna1088182_relationships.tsv`
- `manifests/prjna1088182_truth_coordinate_audit.tsv`

The cohort manifest, SRA runinfo, and publication-truth TSV stay the sample
lists. `manifests/prjna1088182.tsv` is not the cohort. PRJNA884016 is out of
scope. `starter_id` is not a biological parent.

## References

Assembly `GCF_000005845.2` (`ASM584v2`) reports one chromosome:

`U00096.3` Relationship `=` `NC_000913.3`, length 4,641,652.

Sequence SHA256 of the chromosome bases
`e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0`
is identical for NCBI `U00096.3`, NCBI `NC_000913.3`, and the sequence in
local `data/references/MG1655.gbk` / `MG1655.fna`.

The 10 MG1655 / Cas9 samples are `VERIFIED_EQUIVALENT_ACCESSION`:

- publication reference = `U00096.3`
- execution reference = `NC_000913.3`

They are not different biological references.

BW25113 execution uses only `CP009273.1`.

| Accession | Length | Sequence SHA256 |
| --- | ---: | --- |
| `CP009273.1` | 4,631,469 | `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce` |

`NC_000913.3` must not be the execution reference for BW25113 samples.

`manifests/prjna1088182_cohort.tsv` now stores the same accessions in
`publication_reference`, `execution_reference`, and `reference`. The
`reference` column is the execution accession, not a file path.
`prjna1088182_relationships.tsv` column `execution_reference` remains the
authoritative selector. `scripts/m6_reference_policy.py --check` fails closed
when a requested accession contradicts that table. A local sequence file may
be used only after the job record names this accession. For MG1655 that file
is `data/references/MG1655.gbk`. For BW25113 it is the `CP009273.1` chromosome,
not the MG1655 file.

## Run policy

`execution_unit = SRA Run`.

`publication_truth_unit = SampleName/BioSample`.

For each SampleName with two WGS runs:

- future M6-P1/P2 executes each run independently
- reads are not merged by default
- neither run is discarded because of `avgLength`
- Experiment and Submission stay on the job record
- recovery is compared between runs descriptively

Neither run is `PRIMARY_WGS`. They are not technical replicates. Both classes
in the relationship table are `INDEPENDENT_WGS_RUN`. A later metadata finding
may change the interpretation. It does not rewrite which run was which.

## Tiers

### Tier A — reference-relative publication recovery

Eligible when the run, the execution reference above, and a
`COORDINATE_VERIFIED` truth row are the ones named in the job. A sequenced
parent is not required.

Event eligibility:

| Coordinate status | Tier A |
| --- | --- |
| `COORDINATE_VERIFIED` | Eligible for exact publication-event recovery. Included in the denominator. |
| `COORDINATE_UNRESOLVED` | `NOT_EVALUABLE`. Excluded from the denominator. |
| `COORDINATE_CONFLICT` | `NOT_EVALUABLE`. Excluded until resolved. |

Not every publication row must be eligible. Current denominator is the 13
`COORDINATE_VERIFIED` substitution rows, with duplicate `2872949 C>A` counted
once (12 unique alleles). The `+AAT` row stays out.

`publication NONE` means the current truth source reports no mutation for
that sample. It does not mean a mutation-free genome. `NONE` samples may be
executed for pipeline output, additional-call inventory, and run-to-run
comparison. They are not complete-negative truth for specificity or FPR.

No numeric pass threshold. No precision, recall, or F1. No new M1 tolerance.
Publication integers are not matched with the structural ±5 bp rule.

A missing eligible call, once a run exists, is classified as one of
`TRUE NON-RECOVERY`, `REFERENCE MISMATCH`, `COORDINATE MISMATCH`,
`REPRESENTATION DIFFERENCE`, `INSUFFICIENT EVIDENCE`, or `NOT_EVALUABLE`.
Representation difference is not an automatic false negative. M1 is not
retuned to this dataset.

Tier A must not claim editing causality, a parent-derived mutation, or an
editing-attributable unexpected-mutation rate.

### Tier B — parent-relative genome-edit audit

**NOT_EVALUABLE** for all 25 samples. Documented parents have
`NO_PUBLIC_RUN_FOUND`. This does not fail Tier A. M2 intended status stays
`NOT_EVALUABLE`, not `MISSING`.

## Metrics

| Metric | Truth | Denominator | If data are absent |
| --- | --- | --- | --- |
| Exact publication-event recovery | `COORDINATE_VERIFIED` substitutions | those rows, duplicate counted once | `NOT_EVALUABLE`, not 0 |
| Representation disagreement | same rows | same | reported separately |
| breseq 0.40.2 parity | independent oracle, not the paper’s 0.35.4 calls | none yet | `NOT_EVALUABLE` |
| Specificity / FPR from `NONE` | not defined | `NONE` is not that denominator | `NOT_EVALUABLE` |
| M2 intended status | intended TSV is `.` | none | `NOT_EVALUABLE` |
| JC/MOB versus the paper | no such truth rows | none | `NOT_EVALUABLE` |
| Guide search | parsed Cas9 spacer only | rows with a spacer | a hit is not causality |
| DV1 / seven M4 products | M4 V2 internal checks | one output directory | failure is failure |
| Wall time, peak RSS | same-job record | none until a job | `NOT_EVALUABLE`; no speed claim |

## Smoke proposal (not executed)

| SampleName | Runs, each independent | Truth | Reference |
| --- | --- | --- | --- |
| Cas9_1 | `SRR28349244`, `SRR28357593` | `301690 C>T` | `NC_000913.3` (`U00096.3` equivalent) |
| Cas9_2 | `SRR28349243`, `SRR28357592` | `NONE` (not a negative control for FPR) | same |

`Cas9_104_105` is optional later smoke. It is mixed gDNA and has two verified
substitutions. BW25113 smoke waits on a staged `CP009273.1` job record; the
reference policy itself is already frozen.

## Initial paired-smoke amendment (2026-10-08)

Tier-A truth eligibility is unchanged. Each SRA Run remains an execution unit.
`sra-stat -x` and SEQUENCE rows show that `SRR28357593` and `SRR28357592`
(`SRA1825115`) store one biological read per spot (`spot_count_mates=0`).
Reported `LibraryLayout` remains `PAIRED`. Those runs stay in the inventory
and are not eligible for the first paired-end P1b smoke. That smoke, not yet
started, is `SRR28349244` and `SRR28349243` only.

`VALIDATION_REGISTRY.tsv` stays `PENDING` until execution evidence exists.
