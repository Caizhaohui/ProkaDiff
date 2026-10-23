# Milestone M6 Formal Closure Record

**Document ID:** `M6_CLOSEOUT.md`
**Repository:** `Caizhaohui/ProkaDiff` (Public GitHub)
**Git Baseline:** [`c4dc0015e6de43d29479f0334e7a6cab3e89e270`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff)
**Formal Milestone Status:** **`COMPLETE (2026-10-23)`** (Human-Authorized)
**Scientific Closeout Designation:** **`CLOSED_WITH_LIMITATIONS`**
**Milestone State in Governance:** `M6: COMPLETE (2026-10-23)` (per local [DEVELOPMENT_PLAN.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/DEVELOPMENT_PLAN.md))
**Validation Registry State:** [`benchmark/real_world/VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv) line 6 remains strictly **`PENDING`**

---

## 1. Executive Summary & Purpose

This document provides the formal closure record for **Milestone M6 (Real-World Resequencing Benchmark — BioProject `PRJNA1088182`)**.

All computational, algorithmic, and engineering deliverables defined in [DEVELOPMENT_PLAN.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/DEVELOPMENT_PLAN.md) and [M6_VALIDATION_CONTRACT.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/M6_VALIDATION_CONTRACT.md) have been completed, verified on Slurm compute nodes, and audited by independent reviewers.

Crucially, **milestone completion is strictly distinguished from biological variant validation**:
- The computational engineering objectives of Milestone M6 are **COMPLETE**.
- Independent biological truth for unconfirmed candidate variants remains **`NOT_EVALUABLE`**.
- Milestone M6 is formally closed under the scientific designation **`CLOSED_WITH_LIMITATIONS`**, with human project owner authorization granted on 2026-10-23.

---

## 2. Computational Objectives Completed

1. **Dataset Inventory & Provenance (M6-P0):**
   - 25 SRA runs cataloged and mapped to reference strains (*E. coli* MG1655 `NC_000913.3` and BW25113 `CP009273.1`).
   - Two runs (`SRR28357593`, `SRR28357592`) with non-standard single-read spot structures were rigorously excluded (`BLOCKED_READ_STRUCTURE`).
2. **FASTQ Staging & Pipeline Execution (M6-P1):**
   - 23 eligible paired-end runs staged with verified byte counts and gzip SHA-256 digests.
   - Successful execution on Slurm compute partition `qcpu_18i`.
3. **Reference Sequence Bit-Identity:**
   - Both ProkaDiff and breseq aligned reads against bit-identical references:
     - MG1655 (`NC_000913.3`): `4,641,652` bp, sequence SHA-256 `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0`.
     - BW25113 (`CP009273.1`): `4,631,469` bp, sequence SHA-256 `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce`.
4. **Reference Repeat Masking Validation (M6-P3b):**
   - A/B experiment demonstrated that ingesting 401 repeat annotations from `CP009273.1.gbk` filtered out 7 repeat-associated junction candidates (`ANNOTATION_FILTERED_JC`) while preserving point mutations perfectly ($\Delta\text{SNP} = 0, \Delta\text{DEL} = 0$).

---

## 3. Tier A: Reference-Relative Publication-Event Recovery

- **Frozen Definition:** Tier A assesses reference-relative publication-event recovery against published variant lists. Structural breseq comparison belongs to M6-P3c, not to the Tier A truth definition.
- **Publication Ground Truth Denominator:** Assessed 12 unique eligible publication substitution alleles documented in Widney & Copley (2024).
- **Exact Consensus Recovery ($\ge 0.75$ AF cutoff):** Exactly **2 / 12** unique publication alleles were recovered:
  - `NC_000913.3:301690 C>T` in `Cas9_1` (`SRR28349244`): ProkaDiff AF = 0.984; breseq AF = 0.984.
  - `NC_000913.3:3536100 G>T` in `Cas9_103` (`SRR28349177`): ProkaDiff AF = 0.992; breseq AF = 0.992.
- **Uncalled Published Alleles (10 / 12):**
  - Read-level pileup diagnostics across all 10 uncalled loci in [`benchmark/real_world/manifests/prjna1088182_m6_p3a_read_level.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/manifests/prjna1088182_m6_p3a_read_level.tsv) confirmed that all 10 published alternate alleles have read-level support on both strands below the 0.75 consensus threshold.
  - The caller-shaped ALT fraction range across all 10 uncalled loci is **0.198113–0.691244**.
  - Under ProkaDiff's clonal consensus decision model ($\ge 0.75$ AF), these alternate alleles remain unpromoted by design. Biological subclonality remains unverified. They do not indicate alignment failure or caller bugs.
- **Orthogonal Biological Validation:** Orthogonal biological validation remains unavailable. Tier A measures recovery against publication tables, not wet-lab ground truth.

---

## 4. Tier B: Parent-Relative Genome-Edit Audit

- **Frozen Definition:** Tier B assesses parent-relative genome-edit audit requiring matched starter-strain WGS.
- **Status:** **`NOT_EVALUABLE`**
- **Justification:** Tier B is `NOT_EVALUABLE` because verified biological-parent WGS is unavailable. The source publication (Widney & Copley, 2024) did not deposit parental strain WGS runs in public repositories (`NO_PUBLIC_RUN_FOUND`). In accordance with AGENTS.md Rule 6 (mandatory starter WGS), parent-relative differential subtraction could not be performed.
- **Contractual Relationship:** Tier B being `NOT_EVALUABLE` does not block or fail Tier A publication-event recovery.

---

## 5. Structural Variant Corroboration (M6-P3c): Findings & Limitations

Across the 5-sample Benchmark Quintet (`SRR28349173`, `SRR28349176`, `SRR28349244`, `SRR28349231`, `SRR28349226`), evaluated against `breseq 0.40.2` on Slurm `qcpu_18i`. Structural breseq comparison is an independent computational comparator under M6-P3c, not a Tier A publication truth definition and not wet-lab biological truth:

1. **Junction Accounting (1,319 Bidirectional Reconciliation Accounting Rows):**
   - ProkaDiff emitted **806 JC records**; breseq emitted **639 JC evidence records**.
   - These records yield exactly **1,319 bidirectional reconciliation accounting rows**, comprising:
     - **126 paired comparison rows**
     - **1,193 one-sided accounting rows**
   - Partitioned across the 6 frozen, mutually exclusive accounting categories:
     - **`EXACT_MATCH` ($D = 0$ bp):** 91 paired comparison rows (91 ProkaDiff JCs, 91 breseq JCs).
     - **`NEAR_MATCH_WITHIN_TOLERANCE` ($0 < \max(d_1, d_2) \le 5$ bp):** 10 paired comparison rows (10 ProkaDiff JCs, 10 breseq JCs).
     - **Dimension 3 Corroborated Subtotal:** Exactly **101 paired comparison rows** ($91 + 10 = 101$).
     - **`REPEAT_FAMILY_MATCH` (Case 6 paralogous repeat breakends):** 15 paired comparison rows (15 ProkaDiff JCs, 15 breseq JCs).
     - **`NEAR_MISS_BOUNDARY` ($6 \le \max(d_1, d_2) \le 10$ bp, microhomology slips):** 10 paired comparison rows (10 ProkaDiff JCs, 10 breseq JCs; Dimension 3 `NOT_CORROBORATED`).
     - **`BRESEQ_ONLY_REJECTED`:** 513 one-sided accounting rows (0 ProkaDiff JCs, 513 breseq JCs; suppressed by internal breseq coverage/frequency cutoffs).
     - **`PROKADIFF_ONLY`:** 680 one-sided accounting rows (680 ProkaDiff JCs, 0 breseq JCs; sub-threshold candidates and rRNA operon chimeras).
   - *Sum Conservation Invariant:*
     - Paired comparison rows: $91 + 10 + 15 + 10 = 126$.
     - One-sided accounting rows: $513 + 680 = 1,193$.
     - Total accounting rows: $126 + 1,193 = 1,319$.
     - ProkaDiff record sum: $126 + 680 = 806$.
     - breseq record sum: $126 + 513 = 639$.
   - *Epistemic Invariant:* These 1,319 rows represent bidirectional computational reconciliation accounting rows between callers, **not** 1,319 matched graph edges, and **not** 1,319 biological structural events.
2. **Decoupled MOB Concordance:**
   - **`BRESEQ_REPORTED_MOB_COVERAGE`:** **7 / 7 (100.0%)**.
   - **`PROKADIFF_MOB_COMPARATOR_AGREEMENT`:** **7 / 8 (87.5%)**.
   - 6 mutual MOBs match with exact TSD duplication lengths (`EXACT_MOB_MATCH`).
   - 1 mutual MOB exhibits a 1 bp microhomology slip (`CONCORDANT_MOB_1BP_TSD_RESOLVED`).
   - 1 ProkaDiff-exclusive candidate (`MOB 247` at 2,386,696) is strictly classified as an unresolved computational prediction.
   - Both BW25113 samples independently exhibited zero MOB calls (`BOTH_CALLERS_NO_MOB_REPORTED`).

---

## 6. Epistemic Invariants & Scientific Rigor

1. **No Orthogonal Biological Ground Truth:** Dataset PRJNA1088182 lacks orthogonal PCR, Sanger sequencing, or long-read confirmation. Dimension 4 (Biological Validation) is strictly maintained as **`NOT_EVALUABLE`**.
2. **Prohibition of Diagnostic Claims:** All ungrounded epidemiological claims (`BIOLOGICAL_SENSITIVITY`, `BIOLOGICAL_SPECIFICITY`, `BIOLOGICAL_PRECISION`, `TRUE_NEGATIVE`) have been excised.
3. **In-Vivo Equivalence Boundary:** Dual-caller computational agreement reflects shared alignment evidence from short reads and is never equated with validated in-vivo transposition or physical ground truth.
4. **P3b Repeat Filter Boundary:** Seven repeat-associated JC candidates were removed by annotation-dependent filtering. Their biological false-positive status remains unverified; repeat-associated events are not automatically assumed to be artifacts.
5. **Telemetry Observations:** Wall-clock and CPU differences are descriptive cluster measurements. Candidate graph construction overhead is an unprofiled engineering hypothesis, not a proven root cause.

---

## 7. Production Codebase Freeze & Test Verification

- **Production Engine (`crates/**`):** 100% frozen (zero modifications throughout M6).
- **Automated Regression Suites:**
  - `pytest benchmark/real_world/scripts/test_m6_p3c_reconcile.py`: **20 / 20 PASSED** (0.35s).
  - Full Python benchmark suite: **108 / 108 PASSED** (47.2s).
  - Cargo workspace gate: **231 PASSED, 0 FAILED** (exit code 0).

---

## 8. Validation Registry Disposition

- **Disposition:** **`MUST REMAIN STRICTLY PENDING`**
- Line 6 of [`benchmark/real_world/VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv) (`prjna1088182 cohort_25`) remains **`PENDING`** with all validation boolean flags set to `false`.
- A benchmark dataset without orthogonal wet-lab truth cannot legitimately be toggled to validated.
- **Registry Enum Integrity:** The `VALIDATION_REGISTRY.tsv` schema defines only `ACTIVE` and `PENDING`. `CLOSED_WITH_LIMITATIONS` is not an authorized registry enum value and must not be injected into the TSV.

---

## 9. Formal Milestone Closure & Vocabulary Verification

- **Governance Verification:** In project governance ([AGENTS.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/AGENTS.md) and [DEVELOPMENT_PLAN.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/DEVELOPMENT_PLAN.md)), lifecycle milestone states are `ACTIVE`, `COMPLETE (YYYY-MM-DD)`, or `NOT STARTED`.
- **Status Nature:** `CLOSED_WITH_LIMITATIONS` is a **scientific closeout designation** and review verdict; it is **not** an authorized registry status enum in `VALIDATION_REGISTRY.tsv` (where valid values are `ACTIVE` and `PENDING`).
- **Milestone State in Governance:** `M6: COMPLETE (2026-10-23)` (per local [DEVELOPMENT_PLAN.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/DEVELOPMENT_PLAN.md)).
- **Formal Scientific Closeout Designation:** **`CLOSED_WITH_LIMITATIONS`**.
- **Human Authorization:** Formally authorized by Human Project Owner on 2026-10-23. Milestone M6 is closed.
