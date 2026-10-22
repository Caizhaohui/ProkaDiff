# Milestone M6 Proposed Closure Record

**Document ID:** `M6_CLOSEOUT.md`
**Repository:** `Caizhaohui/ProkaDiff` (Public GitHub)
**Git Baseline:** [`c4dc0015e6de43d29479f0334e7a6cab3e89e270`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff)
**Proposed Status:** **`CLOSED_WITH_LIMITATIONS — PENDING HUMAN AUTHORIZATION`**
**Current Milestone State:** `M6 overall: ACTIVE`
**Validation Registry State:** [`benchmark/real_world/VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv) line 6 remains strictly **`PENDING`**

---

## 1. Executive Summary & Purpose

This document provides the formal closure record for **Milestone M6 (Real-World Resequencing Benchmark — BioProject `PRJNA1088182`)**.

All computational, algorithmic, and engineering deliverables defined in [DEVELOPMENT_PLAN.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/DEVELOPMENT_PLAN.md) and [M6_VALIDATION_CONTRACT.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/M6_VALIDATION_CONTRACT.md) have been completed, verified on Slurm compute nodes, and audited by independent reviewers.

Crucially, **milestone completion is strictly distinguished from biological variant validation**:
- The computational engineering objectives of Milestone M6 are **COMPLETE**.
- Independent biological truth for unconfirmed candidate variants remains **`NOT_EVALUABLE`**.
- The proposed status is **`CLOSED_WITH_LIMITATIONS`**, awaiting explicit human project owner authorization.

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

## 3. Tier A Reference-Relative Benchmark: Findings & Limitations

- **Publication Ground Truth Denominator:** 12 unique evaluable substitution alleles documented in Widney & Copley (2024).
- **Consensus Clonal Recoveries ($\ge 0.75$ AF cutoff):** Exactly **2 / 12** alleles were recovered:
  - `NC_000913.3:301690 C>T` in `Cas9_1` (`SRR28349244`): ProkaDiff AF = 0.984; breseq AF = 0.984.
  - `NC_000913.3:3536100 G>T` in `Cas9_103` (`SRR28349177`): ProkaDiff AF = 0.992; breseq AF = 0.992.
- **Uncalled Consensus Alleles (10 / 12):**
  - Deep BAM pileup diagnostics in M6-P3a confirmed robust physical read support on both strands for uncalled loci (e.g. Priority-1 loci at depths 59–217).
  - The allele fractions for these variants range between **$0.292$ and $0.691$**, reflecting the pooled multi-colony sequencing design of the source publication.
  - Under ProkaDiff's clonal consensus decision model ($\ge 0.75$ AF), these sub-clonal alleles are unpromoted by design. They do not indicate alignment failure or caller bugs.

---

## 4. Tier B Biological-Parent-Relative Verification: Status

- **Status:** **`NOT_EVALUABLE`**
- **Justification:** The source publication (Widney & Copley, 2024) did not deposit parental strain WGS runs in public repositories (`NO_PUBLIC_RUN_FOUND`). In accordance with AGENTS.md Rule 6 (mandatory starter WGS), parent-relative differential subtraction could not be performed. Tier B is permanently recorded as `NOT_EVALUABLE` for PRJNA1088182.

---

## 5. Structural Variant Corroboration: Findings & Limitations

Across the 5-sample Benchmark Quintet (`SRR28349173`, `SRR28349176`, `SRR28349244`, `SRR28349231`, `SRR28349226`), evaluated against `breseq 0.40.2` on Slurm `qcpu_18i`:

1. **Bijective Junction Accounting (1,319 Accounting Rows):**
   - ProkaDiff emitted **806 JC records**; breseq emitted **639 JC evidence records**.
   - All records are partitioned into 6 frozen, mutually exclusive accounting categories:
     - **`EXACT_MATCH` ($D = 0$ bp):** 91 pairs.
     - **`NEAR_MATCH_WITHIN_TOLERANCE` ($0 < \max(d_1, d_2) \le 5$ bp):** 10 pairs.
     - **Dimension 3 Corroborated Total:** Exactly **101 pairs** ($91 + 10 = 101$).
     - **`REPEAT_FAMILY_MATCH` (Case 6 paralogous repeat breakends):** 15 pairs.
     - **`NEAR_MISS_BOUNDARY` ($6 \le \max(d_1, d_2) \le 10$ bp, microhomology slips):** 10 pairs (Dimension 3 `NOT_CORROBORATED`).
     - **`BRESEQ_ONLY_REJECTED`:** 513 records (suppressed by internal breseq coverage/frequency cutoffs).
     - **`PROKADIFF_ONLY`:** 680 records (sub-threshold candidates and rRNA operon chimeras).
   - *Sum Conservation:* ProkaDiff $= 91 + 10 + 15 + 10 + 680 = 806$; breseq $= 91 + 10 + 15 + 10 + 513 = 639$.
   - *Epistemic Note:* These 1,319 rows represent computational edge accounting pairings between callers, **not** 1,319 validated biological mutations.
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

---

## 9. Proposed Milestone State

- **Current State:** `M6 overall: ACTIVE`
- **Proposed Closure Designation:** **`CLOSED_WITH_LIMITATIONS`**
- **Action:** This record is prepared and ready for closure. Activation of the closed status awaits explicit authorization from the Human Project Owner.
