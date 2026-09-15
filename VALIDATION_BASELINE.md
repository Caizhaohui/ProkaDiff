# ProkaDiff Validation Baseline

- **baseline_commit**: `ec6be8c0a57377117ffed17c53d4ced34ef0a572`
- **timestamp**: `2026-09-15T09:21:00+08:00`
- **recorded_by**: Gemini (ProkaDiff Senior Rust/Bioinformatics Engineer)

---

## 1. Environment and Tool Versions

| Component | Status | Version | Location |
| :--- | :--- | :--- | :--- |
| **Git Commit** | Verified | `ec6be8c0a57377117ffed17c53d4ced34ef0a572` | `.` |
| **Rust compiler** | Verified | `rustc 1.89.0 (29483883e 2025-08-04)` | `/hpcfs/fhome/caizhh/.cargo/bin/rustc` |
| **Cargo** | Verified | `cargo 1.89.0 (c24e10642 2025-06-23)` | `/hpcfs/fhome/caizhh/.cargo/bin/cargo` |
| **Bowtie2** | Verified | `2.5.4` | `/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/bowtie2` |
| **breseq** | Verified | `0.40.2` | `/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/breseq` |
| **gdtools** | Verified | `0.40.2` | `/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/gdtools` |
| **Cas-OFFinder** | Missing | `NOT AVAILABLE` | Not found in active environment |
| **FlashFry** | Missing | `NOT AVAILABLE` | Not found in active environment |
| **CRISPRitz** | Missing | `NOT AVAILABLE` | Not found in active environment |
| **OS / Architecture** | Verified | `Linux 5.14.0-570.12.1.el9_6.x86_64` | `x86_64 GNU/Linux` |

---

## 2. Test Suite Baseline (`cargo test --workspace`)

- `cargo fmt --all --check`: **PASS**
- `cargo clippy --workspace --all-targets -- -D warnings`: **PASS**
- `cargo test --workspace`: **PASS**

### Breakdown by Crate:
- `prokadiff-evidence`: 96 passed, 0 failed, 0 ignored
- `prokadiff-gd`: 16 passed, 0 failed, 0 ignored
- `prokadiff-offtarget`: 17 passed, 0 failed, 0 ignored
- `prokadiff-report`: 11 passed, 0 failed, 0 ignored
- **Total**: 140 passed, 0 failed, 0 ignored

---

## 3. Parity Status Baseline

- `ambiguous_n`: 0 over, 0 under (exact match)
- `sub_adjacent`: 0 over, 0 under (exact match)
- `synthetic_del`: 0 del_over, 0 del_under (exact match on DEL 501 50); supporting JC exact canonical match
- `synth_snp_indel`: 0 over, 0 under
- `synth_is_mob`: 2 over, 2 under (coordinate jitter $\pm 5$ bp tolerant match)
- `clonal`: red-line mutations exact match; 34 JC emitted with 5 matching within $\pm 5$ bp

---

## 4. Current Identified Scientific Discrepancies & Actions Needed

1. **CFD Scoring**: Current implementation in `scoring/cfd.rs` uses an ad-hoc heuristic penalty (`base_penalty * pos_factor`) rather than the true published Doench 2016 20x4 mismatch matrix + PAM weights. **Action (FIX-002)**: Disable unvalidated CFD in production path immediately until official matrix/FlashFry parity is established.
2. **Hsu Scoring**: Hsu 2013 calculation is uncalibrated against real FlashFry outputs. **Action (FIX-003)**: Mark `HSU2013_VALIDATED = false` and disable in production by default.
3. **Cas-OFFinder Parser**: Current parser attempts to parse bulges and mismatches simultaneously without explicit version segregation. **Action (FIX-005)**: Segregate v2 (mismatch-only) and v3.
4. **Parity Comparison Logic**: `compare_gd.py` strictly checks Python tuples instead of official `gdtools SUBTRACT` output. **Action (FIX-012/013)**: Align with `gdtools SUBTRACT` for strict and normalized parity.
5. **Ambiguous Base Denominator**: `ra.rs` counts $N$ in coverage denominator, potentially deflating allele frequency. **Action (FIX-014)**: Evaluate against breseq oracle and adjust.
6. **Intended Edits Semantics**: Edit records vs observed events are conflated. **Action (FIX-015)**: Separate edit-level and event-level tracking.
