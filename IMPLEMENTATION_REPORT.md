# Implementation Report: ProkaDiff Validation & Scientific-Correctness Sprint

- **Validated Commit**: `40eda5a`
- **Baseline Commit**: `ec6be8c0a57377117ffed17c53d4ced34ef0a572`
- **Report Generated At**: 2026-09-15
- **Toolchain**: rustc 1.89.0 / cargo 1.89.0 / Bowtie2 2.5.4 / breseq 0.40.2 / gdtools 0.40.2

---

## 1. Executive Summary

This report documents the comprehensive execution of the **Validation & Scientific-Correctness Sprint** (Tasks FIX-001 through FIX-020).

The primary objective of this sprint was **scientific rigor and external oracle alignment over feature novelty**:
1. **Disabled unvalidated heuristic scores**: Removed custom approximate CFD penalty models (`base_penalty`, `position_factor`); marked CFD as `disabled` and Hsu2013 as `experimental` pending real external oracle golden datasets.
2. **Oracle Harness Infrastructure**: Established standardized external oracle testing harnesses for **Cas-OFFinder 2.x/3.x**, **FlashFry 1.15**, and **CRISPRitz 2.8.1**, enforcing clean segregation between raw external tool outputs and ProkaDiff normalized schemas.
3. **GenomeDiff Parity Hardening**: Upgraded `compare_gd.py` to use `gdtools SUBTRACT` as the primary strict mutation parity gate (`gdtools_strict_over`/`gdtools_strict_under`), added normalized parity accounting, and introduced an automated JC discrepancy classifier.
4. **Intended Edit Semantics**: Shifted intended mutation masking evaluation from raw event counts to per-edit assessments (`IntendedEditAssessment`: Complete, Partial, Missing), eliminating false partial-call artifacts on multi-junction cassette insertions.
5. **Documentation & Provenance**: Removed causal biological language ("off-target cleavage", "repair stress response") in favor of descriptive, scientifically neutral terminology; added validation status metadata to `summary.txt`.

---

## 2. Component Reports

### 2.1 Variant Calling Engine (Oracle: breseq 0.40.2 + Bowtie2 2.5.4)

#### Mutation Parity
| Dataset | breseq Mutations | ProkaDiff Mutations | Strict Over | Strict Under | Parity Status |
| :--- | :---: | :---: | :---: | :---: | :--- |
| `ambiguous_n` | 0 | 0 | 0 | 0 | **100% Strict Match** |
| `synthetic_del` | 1 (DEL) | 1 (DEL) | 0 | 0 | **100% Strict Match** |
| `synth_snp_indel` | 7 | 7 | 0 | 0 | **100% Strict Match** |
| `sub_adjacent` | 1 (SUB) | 1 (SUB) | 0 | 0 | **100% Strict Match** |
| `synth_is_mob` | 2 (MOB) | 2 (MOB) | 0 | 0 | **100% Tolerant Match** (TSD microhomology sliding) |
| `clonal` (REL606, 7.6M reads) | 8 | 8 | 0 | 0 | **100% Red-line Mutation Recall** (`over_red = 0`) |

#### Junction (JC) Calling Parity
- **Exact Matches**: 5 canonical junctions match breseq within $\pm 5$ bp on `clonal`.
- **Candidate Discrepancies**: ProkaDiff emits 29 unpromoted candidate split-read junctions on 36 bp reads that breseq filters via read-position hash entropy and coverage skew tests.
- **Diagnostics**: `scripts/parity/analyze_jc_discrepancies.py` operationalized to classify coordinate jitter, strand mismatches, duplicates, and low-support candidates.

#### Performance & Resource Utilization (Measured on Slurm partition `qcpu_18i`)
| Dataset | breseq Wall Time | ProkaDiff Wall Time | Speedup | breseq Peak RSS | ProkaDiff Peak RSS |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `clonal` (REL606) | 3,098.01 s | **544.55 s** | **5.69×** | 1,856,948 kB (1.77 GB) | 3,432,960 kB (3.27 GB) |
| `synthetic_del` | 8.81 s | **2.32 s** | **3.80×** | 97,364 kB (95.1 MB) | 104,188 kB (101.7 MB) |
| `synth_is_mob` | 16.88 s | **1.51 s** | **11.18×** | 97,308 kB (95.0 MB) | 107,196 kB (104.7 MB) |
| `synth_snp_indel` | 5.78 s | **0.62 s** | **9.32×** | 97,180 kB (94.9 MB) | 98,760 kB (96.4 MB) |
| `ambiguous_n` | 8.77 s | **0.49 s** | **17.90×** | 97,240 kB (95.0 MB) | 99,504 kB (97.2 MB) |
| `sub_adjacent` | NA | **0.58 s** | NA | NA | 99,548 kB (97.2 MB) |

*(All values measured; unmeasured values explicitly denoted NA per Rule 68).*

---

### 2.2 CRISPR Off-Target Engine

#### Cas-OFFinder 2.x Mismatch Search Parity
- **Status**: **PASS** (12/12 controlled test cases verified against `testdata/offtarget/cas_offinder/reference.fa`).
- **Precision**: 1.000000 (0 false positives).
- **Recall**: 1.000000 (0 false negatives).
- **Scope**: Plus/minus strand, 0–5 mismatches, PAM variants, multi-contig, chromosome + plasmid, contig boundary coordinates.
- **Parser Architecture**: Segregated into `oracle/cas_offinder_v2.rs` (6-column standard mismatch) and `oracle/cas_offinder_v3.rs` (experimental bulge extension).

#### FlashFry 1.15 Scoring Parity (CFD / Hsu2013)
- **Status**: **BLOCKED_EXTERNAL_DEPENDENCY**
- **External Dependency**: FlashFry 1.15 not installed in compute environment.
- **Harness Infrastructure**: Completed in `benchmark/offtarget/flashfry/` (`cases.tsv`, `run_flashfry.sh`, `normalize_flashfry.py`, `compare_scores.py`).
- **Production Guard**: `calculate_cfd_score()` unconditionally returns `None`; `HSU2013_VALIDATED = false`. TSV exports output `NA` for CFD/Hsu scores until golden data verification passes $\ge 100$ cases within $10^{-6}$ tolerance.
- **Parser**: Updated to native FlashFry 1.15 column schema (`contig`, `start`, `stop`, `orientation`, `Doench2016CFDScore`, `Hsu2013`) and 0-based half-open to 1-based inclusive conversion.

#### CRISPRitz 2.8.1 Bulge Search Parity
- **Status**: **BLOCKED_EXTERNAL_DEPENDENCY**
- **External Dependency**: CRISPRitz 2.8.1 not installed in compute environment.
- **Harness Infrastructure**: Completed in `benchmark/offtarget/crispritz/` (`cases.tsv`, `run_crispritz.sh`, `normalize_crispritz.py`, `compare_bulges.py`).
- **Production Guard**: Bulge search flags (`--max-dna-bulge`, `--max-rna-bulge`) remain labeled `EXPERIMENTAL` and default to `0`.

---

## 3. Test Suite Verification

- `cargo fmt --all --check`: **PASS** (0 differences)
- `cargo clippy --workspace --all-targets -- -D warnings`: **PASS** (0 warnings)
- `cargo test --workspace`: **PASS** (201 passed; 0 failed; 1 ignored for Slurm cluster)
  - `prokadiff`: 15 passed
  - `cli_e2e`: 15 passed (1 ignored)
  - `prokadiff-classify`: 28 passed
  - `prokadiff-evidence`: 96 passed
  - `prokadiff-gd`: 16 passed
  - `prokadiff-offtarget`: 20 passed
  - `prokadiff-report`: 11 passed
