# Milestone M6 Final Scientific Close-Gate Review Report

**Document ID:** `M6_FINAL_CLOSE_GATE_REVIEW.md`
**Repository:** `Caizhaohui/ProkaDiff` (Public GitHub)
**Synchronized Git Commit:** [`c4dc0015e6de43d29479f0334e7a6cab3e89e270`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff)
**Target Milestone:** Milestone M6 — Real-World Empirical Resequencing Benchmark (BioProject `PRJNA1088182` Cohort)
**Milestone State:** `M6 overall: ACTIVE`
**Validation Registry Status:** [`benchmark/real_world/VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv) line 6 remains strictly **`PENDING`**
**Final Close-Gate Status:** **`M6_CLOSE_GATE_READY_FOR_HUMAN_DECISION`**

---

## A. Verified Repository State

- **Local Git HEAD:** `c4dc0015e6de43d29479f0334e7a6cab3e89e270` (`docs: sanitize pre-push structural claims and performance caveats`).
- **Remote `origin/main`:** `c4dc0015e6de43d29479f0334e7a6cab3e89e270` (100% bit-identical, verified via `git fetch origin && git rev-parse HEAD origin/main`).
- **Ancestry Status:** Strict fast-forward descendant (`0 ahead, 0 behind`).
- **Production Engine Code Freeze:** Directory `crates/**` contains **0 modifications** and has remained 100% frozen throughout all M6 benchmark phases.
- **Rule & Roadmap Integrity:** Files `AGENTS.md` and `DEVELOPMENT_PLAN.md` remain untouched.
- **Validation Registry Status:** Line 6 of `benchmark/real_world/VALIDATION_REGISTRY.tsv` (`prjna1088182 / cohort_25`) is preserved strictly in its preregistered **`PENDING`** state with all validation boolean flags set to `false`.

---

## B. Agent Execution Ledger

In strict accordance with the project multi-agent governance rules, four specialized agents participated in the close-gate review with zero unauthorized substitutions:

| Role | Agent | Platform / Model | Reasoning Level | Execution Mode | Scope & Responsibilities |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **Chief Orchestrator** | **Antigravity** | Gemini 3.8 Flash | `MEDIUM` / `FAST` | Coordination / Governance | Protocol enforcement, forensic verification, cross-agent synthesis, reporting |
| **Scientific Reviewer** | **Claude** | Opus 5.5 | `HIGH` / `STANDARD` | `READ_ONLY` | Phase 1 scientific architecture review across all M6 sub-phases (P0–P3c) |
| **Independent Auditor** | **Codex** | GPT-6.1 | `EXTRA HIGH` | `READ_ONLY` | Phase 2 independent audit, error correction, gate-by-gate verification, registry defense |
| **Implementation Engineer** | **Grok** | Cursor Grok 4.7 | `HIGH` / `FAST` | `READ_ONLY` (Standby) | Standby for complex engineering remediation; no production code touched |

---

## C. M6 Phase-by-Phase Acceptance Matrix

| Sub-Phase | Core Objectives | Verified Empirical Evidence on Disk | Epistemic Boundary & Limitations | Codex Gate Verdict |
| :--- | :--- | :--- | :--- | :---: |
| **M6-P0: Inventory & Contract** | - Map 25 SRA samples.<br>- Audit library layouts.<br>- Assign reference contigs.<br>- Freeze Tier A denominator. | - `prjna1088182_cohort.tsv` maps 25 samples.<br>- `sra-stat -x` audit confirmed `SRR28357593` and `SRR28357592` have 0 mate spots (`BLOCKED_READ_STRUCTURE`).<br>- MG1655 (`NC_000913.3`) & BW25113 (`CP009273.1`) assigned.<br>- Tier A frozen at 12 unique evaluable substitution alleles. | Tier B marked `NOT_EVALUABLE` due to `NO_PUBLIC_RUN_FOUND` for documented parental strains. | **COMPLETE** |
| **M6-P1: FASTQ Staging & Smoke** | - Stage 23 eligible paired runs.<br>- Compute SHA-256 checksums.<br>- Slurm smoke runs on `qcpu_18i`. | - `prjna1088182_m6_p1a_fastq_provenance.tsv` records exact byte counts and gzip SHA-256 for all 23 runs.<br>- Smoke jobs `2894179` and `2894180` completed with exit code 0. | None. | **COMPLETE** |
| **M6-P2b: Published Allele Recovery** | - Reference-relative recovery against 12 unique published alleles.<br>- Clonal consensus calling at $\ge 0.75$ AF. | - Exact consensus recovery: 2/12 (`301690 C>T` in `SRR28349244` AF=0.984; `3536100 G>T` in `SRR28349177` AF=0.992).<br>- 10/12 uncalled by consensus filter.<br>- Verified read-level support below 0.75 consensus threshold (caller-shaped ALT fraction range 0.198113–0.691244); biological subclonality remains unverified. | Clonal consensus caller does not recover alleles below consensus threshold diluted in pooled gDNA libraries. | **COMPLETE_WITH_LIMITATIONS** |
| **M6-P3a: Read-Level Diagnostics** | - Measure raw and caller-shaped read support at Priority-1 loci.<br>- Invariant conservation checks. | - Diagnostic pileups across uncalled loci confirm read support on both strands; caller-shaped ALT fraction range across all 10 uncalled loci is 0.198113–0.691244.<br>- All loci show physical read support with depths 59–239.<br>- Mathematical invariants hold: raw depth and caller-shaped depth sums balance perfectly. | Validates presence of alleles in reads without claiming caller defect or altering consensus threshold. | **COMPLETE** |
| **M6-P3b: Reference Repeat Masking** | - A/B experiment: plain FASTA vs GenBank repeat features.<br>- BW25113 triad evaluation (`SRR28349226`, `SRR28349231`, `SRR28349230`). | - 401 repeat features parsed from `CP009273.1.gbk`.<br>- $\Delta\text{SNP} = 0, \Delta\text{DEL} = 0$ (point mutations invariant).<br>- Seven repeat-associated JC candidates were removed by annotation-dependent filtering. Their biological false-positive status remains unverified.<br>- Both callers report 0 MOBs in both conditions. | Eliminates repeat-associated junction candidates without altering point mutation calls or implying that repeat-associated events are automatically artifacts. | **COMPLETE** |
| **M6-P3c / P3c.3a: Structural Reconciliation** | - Cohort structural comparison across 5 runs.<br>- Kuhn maximum-cardinality bipartite matching.<br>- 6-category bidirectional accounting of 1,319 rows (126 paired comparison rows, 1,193 one-sided accounting rows).<br>- Decoupled MOB concordance. | - All 5 Slurm comparator runs (`breseq 0.40.2`) exited 0 on `qcpu_18i`.<br>- Kuhn bipartite matching proved greedy matching globally optimal ($d_{max} = 1$).<br>- 1,319 rows partitioned: 91 EXACT, 10 NEAR, 15 REPEAT_FAM, 10 NEAR_MISS, 513 BRESEQ_REJ, 680 PK_ONLY.<br>- Dimension 3 corroborated = 101 paired rows.<br>- MOB concordance: breseq coverage = 7/7 (100%), PK agreement = 7/8 (87.5%). | PRJNA1088182 lacks orthogonal PCR/Sanger ground truth; Dimension 4 is strictly `NOT_EVALUABLE`. | **COMPLETE_WITH_LIMITATIONS** |

---

## D. Scientific Evidence Summary

### 1. Data Provenance & Reference Identity
- **Biological Source:** BioProject PRJNA1088182 (Widney & Copley, 2024).
- **Run Allocation:** 25 total runs. Two companion runs (`SRR28357593`, `SRR28357592`) were excluded during M6-P0 due to `BLOCKED_READ_STRUCTURE` (0 mate spots).
- **Reference Sequences:**
  - *E. coli* K-12 MG1655 (`NC_000913.3`): Length `4,641,652` bp; normalized uppercase sequence SHA-256: `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0`.
  - *E. coli* BW25113 (`CP009273.1`): Length `4,631,469` bp; normalized uppercase sequence SHA-256: `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce`.
  - Bit-identity between FASTA and GenBank ORIGIN sequences is 100% verified across both callers.

### 2. Point Mutation Recovery — Tier A: Reference-Relative Publication-Event Recovery
- **Frozen Contract Definition:** Tier A assesses reference-relative publication-event recovery against published variant lists. Structural breseq comparison belongs to M6-P3c, not to the Tier A truth definition. Tier B (parent-relative genome-edit audit) is permanently `NOT_EVALUABLE` due to missing parental WGS (`NO_PUBLIC_RUN_FOUND`). Orthogonal biological validation remains unavailable.
- **Publication Ground Truth Denominator:** Assessed 12 unique eligible publication substitution alleles documented in Widney & Copley (2024).
- **Exact Consensus Recovery ($\ge 0.75$ AF):** Exactly **2 / 12** published alleles were recovered:
  - `NC_000913.3:301690 C>T` in `Cas9_1` (`SRR28349244`): ProkaDiff AF = 0.984; breseq AF = 0.984.
  - `NC_000913.3:3536100 G>T` in `Cas9_103` (`SRR28349177`): ProkaDiff AF = 0.992; breseq AF = 0.992.
- **Uncalled Published Alleles (10 / 12):** Investigated via read-level pileup diagnostics in M6-P3a (`prjna1088182_m6_p3a_read_level.tsv`). Verified that all 10 published alternate alleles have read-level support on both strands below the 0.75 consensus threshold, with a caller-shaped ALT fraction range of **0.198113–0.691244**. Biological subclonality remains unverified. Under ProkaDiff's clonal consensus decision model ($\ge 0.75$ AF), these alleles remain unpromoted by design.

### 3. Annotation Dependency & Repeat Masking
- Seven repeat-associated JC candidates were removed by annotation-dependent filtering. Their biological false-positive status remains unverified; repeat-associated structural events are not automatically artifacts. Point mutations were preserved perfectly ($\Delta\text{SNP} = 0, \Delta\text{DEL} = 0$).
- Both callers independently called zero MOB events in BW25113 controls (`BOTH_CALLERS_NO_MOB_REPORTED`), with biological truth strictly preserved as `NOT_EVALUABLE`.

### 4. Cohort Structural Event Accounting (M6-P3c Computational Corroborator)
Across the 5 benchmark runs (`SRR28349173`, `SRR28349176`, `SRR28349244`, `SRR28349231`, `SRR28349226`), all 806 ProkaDiff junction records and 639 breseq junction records were reconciled into exactly **1,319 bidirectional reconciliation accounting rows**, comprising:
- **126 paired comparison rows**
- **1,193 one-sided accounting rows**

Preserving the original 6 frozen, mutually exclusive accounting categories:
- **`EXACT_MATCH` ($D = 0$):** 91 paired comparison rows (91 ProkaDiff JCs, 91 breseq JCs).
- **`NEAR_MATCH_WITHIN_TOLERANCE` ($0 < \max(d_1, d_2) \le 5$ bp):** 10 paired comparison rows (10 ProkaDiff JCs, 10 breseq JCs).
- **Dimension 3 Corroborated Subtotal:** Exactly **101 paired comparison rows** ($91 + 10 = 101$).
- **`REPEAT_FAMILY_MATCH` (Case 6 paralogous repeat breakends):** 15 paired comparison rows (15 ProkaDiff JCs, 15 breseq JCs).
- **`NEAR_MISS_BOUNDARY` ($6 \le \max(d_1, d_2) \le 10$ bp, microhomology slips):** 10 paired comparison rows (10 ProkaDiff JCs, 10 breseq JCs; Dimension 3 `NOT_CORROBORATED`).
- **`BRESEQ_ONLY_REJECTED`:** 513 one-sided accounting rows (0 ProkaDiff JCs, 513 breseq JCs; suppressed by internal breseq coverage/frequency cutoffs).
- **`PROKADIFF_ONLY`:** 680 one-sided accounting rows (680 ProkaDiff JCs, 0 breseq JCs; sub-threshold candidates and rRNA operon chimeras).
- **Sum Checks:**
  - Paired comparison rows $= 91 + 10 + 15 + 10 = \mathbf{126}$.
  - One-sided accounting rows $= 513 + 680 = \mathbf{1,193}$.
  - Total accounting rows $= 126 + 1,193 = \mathbf{1,319}$.
  - ProkaDiff $= 126 + 680 = \mathbf{806}$.
  - breseq $= 126 + 513 = \mathbf{639}$.

*Epistemic Invariant:* These 1,319 rows represent bidirectional computational reconciliation accounting rows between callers, **not** 1,319 matched graph edges, and **not** 1,319 biological structural events.

### 5. Decoupled MOB Concordance
- **Breseq Coverage:** $7 / 7$ ($100.0\%$).
- **ProkaDiff Agreement:** $7 / 8$ ($87.5\%$).
- 6 mutual MOBs match with exact TSD duplication lengths (`EXACT_MOB_MATCH`).
- 1 mutual MOB exhibits a 1 bp microhomology slip (`CONCORDANT_MOB_1BP_TSD_RESOLVED`).
- 1 ProkaDiff-exclusive candidate (`MOB 247` at 2,386,696) is strictly classified as an unresolved computational candidate.

---

## E. Unresolved Limitations

1. **Clonal Consensus Horizon:** ProkaDiff operates under a clonal consensus model requiring $\ge 0.75$ allele frequency. Published alternate alleles with caller-shaped ALT fraction below threshold (observed range: 0.198113–0.691244) remain uncalled by design; biological subclonality remains unverified.
2. **Missing Parental WGS (Tier B `NOT_EVALUABLE`):** The lack of publicly deposited parental WGS runs for PRJNA1088182 permanently prevents parent-relative intended edit verification on this cohort.
3. **Algorithmic Structural Divergence:** The 680 ProkaDiff-only and 513 breseq-only junction records reflect divergent split-read clustering and filtering heuristics. Determining their physical validity requires orthogonal long-read sequencing (PacBio/ONT) or targeted PCR.
4. **Dimension 4 Ground Truth:** PRJNA1088182 lacks orthogonal experimental ground truth; Dimension 4 is strictly **`NOT_EVALUABLE`**.
5. **Cluster Telemetry Limitations:** Wall-clock and CPU overhead differences are descriptive empirical measurements; candidate graph building overhead remains an unprofiled hypothesis.

---

## F. Benchmark Reproducibility Assessment

- **Execution Environment:** Slurm partition `qcpu_18i` on compute nodes `bnode1` and `bnode2` using pinned toolchains (`bowtie2 2.5.4`, `breseq 0.40.2`).
- **Deterministic Test Suites:**
  - `python3 -m pytest benchmark/real_world/scripts/test_m6_p3c_reconcile.py`: **20 / 20 PASSED** (0.35s).
  - `python3 -m pytest benchmark/real_world/scripts/`: **108 / 108 PASSED** (47.2s).
  - `cargo test --workspace`: **231 PASSED, 0 FAILED** (exit code 0).
- **Reproducibility Manifests:** Complete telemetry, bidirectional accounting matrices, and reconciliation tables are versioned and committed in [`benchmark/real_world/manifests/`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/manifests).

---

## G. Independent Codex Review Verdict

- **Reviewer:** Codex (GPT-6.1 EXTRA HIGH, Mode: `READ_ONLY`)
- **Formal Verdict:** **`PASS_WITH_LIMITATIONS`**
- **Findings:**
  - All computational, algorithmic, and engineering deliverables specified for Milestone M6 are complete and reproducible.
  - Zero modifications to production crates (`crates/**` 100% frozen).
  - All known scientific limitations are fully documented, mathematically conserved, and rigorously bounded.
  - All unsupported claims of biological ground truth, sensitivity, specificity, and precision were successfully eliminated.

---

## H. Explicit Milestone Closure Recommendation

### Recommendation: **`READY_FOR_CLOSURE_WITH_LIMITATIONS`** (Awaiting Human Authorization)

The engineering deliverables and computational benchmark objectives of Milestone M6 are complete. The scientific findings are verified against physical artifacts, and biological ground truth is honestly bounded as `NOT_EVALUABLE`.

Milestone M6 is recommended to be closed under the narrative closeout designation **`CLOSED_WITH_LIMITATIONS`**, subject to explicit authorization by the Human Project Owner. Note that in project governance ([AGENTS.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/AGENTS.md) and [DEVELOPMENT_PLAN.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/DEVELOPMENT_PLAN.md)), `CLOSED_WITH_LIMITATIONS` functions as a narrative closeout designation and scientific review verdict, rather than an authorized registry enum value in `VALIDATION_REGISTRY.tsv`. Milestone state remains `M6 overall: ACTIVE` until human authorization.

---

## I. Validation Registry Disposition

### Disposition: **MUST REMAIN STRICTLY `PENDING`**

Codex and the Orchestrator explicitly reject marking row 6 of `benchmark/real_world/VALIDATION_REGISTRY.tsv` as validated:
1. **Schema Integrity:** The registry schema only recognizes `PENDING` and `ACTIVE`. `CLOSED_WITH_LIMITATIONS` is not an authorized registry enum value and must not be injected into the TSV.
2. **Absence of Biological Ground Truth:** Marking `variant_validated`, `intended_validated`, or `mob_validated` as `true` would constitute an ungrounded scientific assertion because:
   - Only 2/12 published consensus SNPs were recovered (10 published alternate alleles had read-level support below the 0.75 consensus threshold, with caller-shaped ALT fraction range 0.198113–0.691244; biological subclonality remains unverified).
   - Tier B is permanently `NOT_EVALUABLE` due to missing parental WGS.
   - PRJNA1088182 lacks orthogonal experimental validation for structural variants.
3. **Contractual Invariant:** Row 6 of [`benchmark/real_world/VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv) must remain **`PENDING`** with all validation booleans set to `false`.

---

## J. Recommendations for Milestone M7

1. **Sub-Clonal / Heterogeneous Mixture Caller:** Design an opt-in sub-clonal variant caller mode (e.g. `--subclonal` with configurable allele fraction down to $\sim 0.10$ and beta-binomial error modeling) to address pooled multi-colony resequencing libraries.
2. **Junction Graph Profiling & Performance Tuning:** Generate CPU flame graphs and line-level memory profiles of the split-read indexing and candidate junction clustering engine to systematically address runtime overhead.
3. **Long-Read Ground Truth Integration:** Integrate synthetic or public long-read datasets (Oxford Nanopore / PacBio HiFi) to benchmark structural variant resolution across multi-copy ribosomal operons (`rrnA`–`rrnH`) and paralogous transposable elements.
4. **Calibrated Candidate Filtering:** Formulate principled coverage-evenness, strand-bias, and homopolymer-overhang filters to separate true sub-clonal rearrangements from PCR/library chimeras in low-support candidates ($N < 10$).

---

**Final Status:** **`M6_CLOSE_GATE_READY_FOR_HUMAN_DECISION`**
*(No automatic push, tag, release, registry update, or M7 execution has been performed).*
