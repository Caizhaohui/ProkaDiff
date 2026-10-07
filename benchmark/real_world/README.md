# ProkaDiff Real-World Validation Framework

This directory contains the automated framework for benchmarking and validating ProkaDiff against real-world prokaryotic editing datasets, published literature cohorts, and assembly-confirmed ground truth.

## Directory Structure

```text
benchmark/real_world/
├── README.md
├── manifests/
│   ├── bl21_user.tsv          # Real Cas9-edited BL21 clones (B21_2_1, B21_3_1, B21_4_1)
│   ├── prjna1088182.tsv       # Widney/Copley 2024 deep-resequencing cohort
│   ├── prjna1088182_cohort.tsv # Curated 25-sample cohort declaration
│   ├── prjna884016.tsv        # IS/transposon structural ground truth cohort
│   └── rel606.tsv             # REL606 long-term evolution & mobile element truth
├── truth/
│   └── bl21_curated_truth.tsv # Curated breakpoints and SV annotations for BL21 clones
│   ├── prjna1088182_publication_truth.tsv # Reported, not WGS-measured, variants
│   └── prjna884016/            # Assembly-derived GD/MOB/JC structural truth
├── scripts/
│   ├── download_dataset.sh    # SRA fetch utility for public accession datasets
│   ├── run_prokadiff.sh       # Runner for ProkaDiff differential audit
│   ├── run_breseq.sh          # Runner for breseq oracle comparison
│   ├── compare_variant_calls.py   # Evaluates variant precision/recall vs oracle/truth
│   ├── compare_intended_edits.py  # Checks M4 intended-edit rows; peer runs do not preset a status
│   ├── compare_mob_truth.py       # Evaluates mobile element insertion identification
│   ├── validate_reports.py        # Validates markdown reports and neutrality checks
│   └── aggregate_results.py       # Summarizes benchmark metrics across all datasets
└── results/
    ├── bl21_user/
    ├── prjna1088182/
    ├── prjna884016/
    └── rel606/
```

## Validation evidence contract

`VALIDATION_REGISTRY.tsv` is the release-facing source of truth. A `true` validation
field requires all of the following to be versioned or directly reproducible from
versioned inputs:

1. A manifest and declared truth source.
2. ProkaDiff calls from the declared sample pair.
3. A comparison result that identifies the exact calls and truth used.
4. For runtime or breseq-parity claims, the same-job evidence required by
   [`docs/parity.md`](../../docs/parity.md).

The PRJNA1088182 publication table and the PRJNA884016 assembly files are validation
inputs. They are not measured ProkaDiff results. Their registry rows remain `PENDING`
until the queue-only validation milestones write and preserve a comparison record.
Generated `results/` directories and raw FASTQ/BAM remain local and ignored.

## BL21 M5 validation inputs and products

M5 (BL21 real-data validation) is **COMPLETE**. The authoritative closure
record is [`M5_CLOSEOUT.md`](M5_CLOSEOUT.md). Single-sample
source-separated reconciliation remains in
[`M5_BL21_P1_RECONCILIATION.md`](M5_BL21_P1_RECONCILIATION.md). Generated
`results/` stay local and gitignored; the closeout summarizes observed
job identities and product validation, not versioned FASTQ/BAM payloads.

`manifests/bl21_user.tsv` is versioned as `BL21_INPUTS_V2`. Each row explicitly
declares `comparison_role=PEER_COMPARATOR` and
`biological_parent_verified=false`; `starter_id` identifies the computational
comparator only and does not establish biological parenthood. The role vocabulary
also supports `REFERENCE_ONLY` and `MATCHED_PARENT` for future validation cases.
It lists all three
sample/starter pairs, the six private FASTQ paths, the GCF_013167015.1 FASTA,
the shared intended-edit table, and both guides as separate named spacer
columns (`BL21_G1`, `BL21_G2`). `scripts/validate_bl21_registry.py` checks
sample and pair names, file presence, reference accession, FASTA/GenBank
contig sequence identity, intended coordinates, curated structural calls,
and both spacer declarations. It stats FASTQ files but does not read them.

The BL21 product jobs use schema v2 and validate the M4 outputs directly:
`summary.txt`, `report.md`, `edit_outcomes.tsv`, `post_edit_variants.tsv`,
`unintended.tsv`, `offtarget_sites.tsv`, and `mutation_offtarget_links.tsv`.
Each sample/comparator comparison runs once per named guide into separate product
directories. Reference-relative single-sample review and peer-comparator
assessment are separate questions. A clone can show an abnormal target locus
relative to `GCF_013167015.1` while the clone-to-clone differential assessment
is `MISSING`, `COMPLETE`, `PARTIAL`, or `UNEXPECTED_STRUCTURE`, according to
the differential events that remain. Peer-comparator `MISSING` is valid when no
intended-related differential event remains, both DV1 edit lists are empty, and
no product row contradicts that status. It does not mean the target locus is
normal relative to the reference. `COMPLETE` is valid when the differential
evidence supports the expected edit. Neither outcome is preset from curated
single-sample structure. Both guide runs continue independently. Single-sample
target review is a separate validation TSV with source-specific ProkaDiff,
breseq, and curated records; it does not generate M4 events or assessments.
`repeat_ambiguous_junctions.tsv` stays diagnostic-only and is not an EventId
source. Raw bidirectional `gdtools SUBTRACT A B` and `gdtools SUBTRACT B A`
residuals are retained before any tolerance/whitelist classification.
Submit `scripts/run_bl21_validation.sbatch` to `qcpu_18i`. The focused
`scripts/run_p1_3_diff.sbatch` accepts the three frozen directions
`B21_3_1`/`B21_4_1`, `B21_4_1`/`B21_3_1`, and `B21_2_1`/`B21_3_1`, and runs G1
and G2 independently for the selected direction. Neither script is run on the
login node.

Run `bash scripts/check_validation_assets.sh` after changing a cohort manifest, truth asset,
or registry row. It is a pure metadata check: it does not read FASTQ, run Bowtie2, or
claim a validation result.

## Key Release Metrics

- **FALSE_COMPLETE == 0**: Crucial release gate. If a strain suffered an unintended aberrant rearrangement or large deletion at the target locus, ProkaDiff must NEVER report `COMPLETE`.
- **Variant Precision / Recall**: Evaluated against oracle `breseq 0.40.2` and assembly truth.
- **IS / MOB Precision**: Accurate identification of transposon insertion sites and families.
- **Report Audit Neutrality**: Distal small variants must not receive speculative causal labels.
