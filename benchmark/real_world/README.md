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
│   ├── compare_intended_edits.py  # Evaluates intended edit status (PASS/PARTIAL/UNEXPECTED)
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

Run `bash scripts/check_validation_assets.sh` after changing a cohort manifest, truth asset,
or registry row. It is a pure metadata check: it does not read FASTQ, run Bowtie2, or
claim a validation result.

## Key Release Metrics

- **FALSE_COMPLETE == 0**: Crucial release gate. If a strain suffered an unintended aberrant rearrangement or large deletion at the target locus, ProkaDiff must NEVER report `COMPLETE`.
- **Variant Precision / Recall**: Evaluated against oracle `breseq 0.40.2` and assembly truth.
- **IS / MOB Precision**: Accurate identification of transposon insertion sites and families.
- **Report Audit Neutrality**: Distal small variants must not receive speculative causal labels.
