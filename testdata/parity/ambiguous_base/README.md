# Mixed-N Denominator Parity (FIX-014)

## Issue

`crates/prokadiff-evidence/src/ra.rs` uses `observations.len()` as the coverage
denominator for allele frequency (AF) computation. When reads carry ambiguous
bases (`N`), those observations are counted in the denominator but **do not
enter any specific allele bucket** (A/C/G/T).

This may cause AF underestimation relative to breseq:

```
N observations inflate denominator → AF = alt_count / (alt_count + ref_count + N_count)
```

breseq's actual behavior must be confirmed by running a controlled oracle experiment.

## Status

```
BLOCKED_ORACLE_RUN
```

An oracle run on Slurm partition `qcpu_18i` with the synthetic ambiguous-base
FASTQ in `testdata/parity/ambiguous_base/cases/` is required before modifying `ra.rs`.

**No change to `ra.rs` may be made without first observing breseq's denominator
behavior on the controlled test cases.**

## Test Cases

Located in `testdata/parity/ambiguous_base/cases/`:

| Case | Reads                          | Expected denominator (TBD) |
|------|--------------------------------|---------------------------|
| A    | 10 G + 0 N                    | TBD (oracle run required)  |
| B    | 10 G + 1 N                    | TBD                       |
| C    | 10 G + 5 N                    | TBD                       |
| D    | 10 G + 10 N                   | TBD                       |
| E    | 5 G + 5 ref + 10 N            | TBD                       |

## Decision Procedure

1. Submit Slurm job with `testdata/parity/ambiguous_base/run_oracle.sh` to `qcpu_18i`
2. Capture breseq RA evidence (`output/evidence/evidence.gd`)
3. Record denominator in `oracle_results.tsv`
4. If breseq excludes N from denominator:
   ```rust
   let informative_cov: u32 = counts.iter().sum(); // A+C+G+T only
   ```
5. If breseq includes N, keep `observations.len()`
6. Add regression tests

## Files

- `run_oracle.sh` — Slurm submission script
- `cases/` — synthetic read FASTQ (to be generated; see `generate_cases.sh`)
- `oracle_results.tsv` — to be filled after oracle run
