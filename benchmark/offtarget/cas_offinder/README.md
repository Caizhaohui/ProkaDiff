# Cas-OFFinder Validation Harness

This harness evaluates ProkaDiff exact mismatch-only search against Cas-OFFinder 2.x across both strands, multiple contigs, and mismatches 0..5.

## Test Dataset
- Reference: `testdata/offtarget/cas_offinder/reference.fa` (generated via `generate_ref.py`)
  - Contig `chr` (1000 bp)
  - Contig `plasmid` (500 bp)
- Matrix: 12 controlled cases defined in `cases.tsv`.

## Execution
```bash
# Run Cas-OFFinder oracle (if binary installed)
bash scripts/run_cas_offinder.sh

# Run ProkaDiff standalone search
bash scripts/run_prokadiff_search.sh

# Compare outputs
python3 scripts/compare_sites.py --oracle cas_offinder.raw.tsv --prokadiff prokadiff.tsv
```

## Parity Acceptance Gate
- Precision: 1.000000 (0 false positives)
- Recall: 1.000000 (0 false negatives)
