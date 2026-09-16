# ProkaDiff Real-World Validation Framework

This directory contains the automated framework for benchmarking and validating ProkaDiff against real-world prokaryotic editing datasets, published literature cohorts, and assembly-confirmed ground truth.

## Directory Structure

```text
benchmark/real_world/
├── README.md
├── manifests/
│   ├── bl21_user.tsv          # Real Cas9-edited BL21 clones (B21_2_1, B21_3_1, B21_4_1)
│   ├── prjna1088182.tsv       # Widney/Copley 2024 deep-resequencing cohort
│   ├── prjna884016.tsv        # IS/transposon structural ground truth cohort
│   └── rel606.tsv             # REL606 long-term evolution & mobile element truth
├── truth/
│   └── bl21_curated_truth.tsv # Curated breakpoints and SV annotations for BL21 clones
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

## Key Release Metrics

- **FALSE_COMPLETE == 0**: Crucial release gate. If a strain suffered an unintended aberrant rearrangement or large deletion at the target locus, ProkaDiff must NEVER report `COMPLETE`.
- **Variant Precision / Recall**: Evaluated against oracle `breseq 0.40.2` and assembly truth.
- **IS / MOB Precision**: Accurate identification of transposon insertion sites and families.
- **Report Audit Neutrality**: Distal small variants must not receive speculative causal labels.
