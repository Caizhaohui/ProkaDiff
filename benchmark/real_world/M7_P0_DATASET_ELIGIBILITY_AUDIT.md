# Milestone M7-P0: PRJNA884016 Assembly Truth Eligibility Audit

**Document ID:** `M7_P0_DATASET_ELIGIBILITY_AUDIT.md`
**Target Milestone:** Milestone M7 — PRJNA884016 Assembly Truth Structural Verification
**BioProject:** `PRJNA884016` (NCBI SRA / GenBank)
**Primary Publication:** Sheng et al., *"Insertion sequence transposition inactivates CRISPR-Cas immunity"*, *Nature Communications* 14, 4366 (July 2023), [PMC10359306](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10359306/)
**Audit Date:** 2026-10-23
**Status:** **`M7_P0_READY_FOR_HUMAN_REVIEW`**

---

## 1. Executive Summary & Core Forensic Discovery

Milestone M7 was planned in [DEVELOPMENT_PLAN.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/DEVELOPMENT_PLAN.md) to evaluate ProkaDiff's structural variant (MOB/JC) calling against closed assembly truth from BioProject `PRJNA884016` (*E. coli* strain `SY22A1` vs `SY221`). Line 7 of [`benchmark/real_world/VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv) tracked `SY22A1` under `truth_level = ASSEMBLY_CONFIRMED`, based on a versioned truth asset ([`SY22A1_truth.gd`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/truth/prjna884016/SY22A1_truth.gd)) reporting 4 mobile element insertions (MOBs) and 7 breakpoint junctions (JCs) between `CP107005.1` and `CP106998.1`.

A comprehensive multi-agent forensic audit (Antigravity Orchestrator + Claude Opus 5.5 Scientific Architect) incorporating direct sequence retrieval from NCBI Assembly has revealed a **fundamental scientific discrepancy**:

1. **100.000% Chromosome Bit-Identity:**
   Direct download and pairwise sequence alignment of the 12 complete NCBI chromosome assemblies deposited under PRJNA884016 (`GCA_030345415.1` – `GCA_030345635.1`, accessions `CP106994.1` – `CP107005.1`, all length 4,687,476 bp) demonstrates that **all 11 CRISPR-tolerant mutant chromosomes (including `SY22A1` / `CP107005.1`) are 100.000% bit-identical to the uninduced control chromosome `SY221` (`CP106998.1`) across all 4,687,476 base pairs (`diff_count = 0`, zero SNPs, zero indels, zero structural variants)**.
2. **Concordance with Published Findings:**
   This bit-identity directly confirms the primary finding of Sheng et al. (2023):
   > *"Sequence alignment of SY221 with the CRISPR-tolerant mutants revealed 100% sequence identity at 100% chromosome coverage, suggesting that the transpositions of ISs into cas genes did not result in simultaneous transpositions into other genomic loci or trigger genomic rearrangements."*
3. **Plasmids Were Not Assembled or Deposited:**
   The paper's experimentally validated IS1 and IS10 transpositions occurred exclusively into *cas* genes on **extrachromosomal cas-expression plasmids** (`pCasY-3-ΔSIM-kan`), NOT on the chromosome. However, **zero plasmid sequences were assembled or deposited in NCBI under BioProject PRJNA884016**.
4. **Root Cause of Local Truth Fabrication:**
   The local file `data/prjna884016/assemblies/SY22A1.fna` was truncated at base 4,492,785 (omitting 194,691 bp). `build_prjna884016_truth.py` parsed MUMmer `rdiff` against this truncated tail, misidentified unaligned coordinates starting at `4,492,786` as biological rearrangements, and heuristically fabricated `SY22A1_truth.gd`. In physical reality and in NCBI accession `CP107005.1`, no such truncation or structural variants exist.

---

## 2. BioProject PRJNA884016 Provenance & Replicon Audit

BioProject PRJNA884016 comprises 12 BioSamples and 24 SRA sequencing runs:

| BioSample | Strain Name | Strain Phenotype / Biological Role | Complete Chromosome Accession | Chromosome Length | Replicon Molecules in Assembly | SRA WGS Runs (Illumina / Long-Read) |
| :--- | :--- | :--- | :---: | :---: | :---: | :--- |
| `SAMN31007788` | `SY221` | Uninduced control (Cas-plasmid without sgRNA) | `CP106998.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706776` (ILMN) / `SRR21706786` (PacBio) |
| `SAMN31007783` | `SY22A1` | CRISPR-tolerant escaper mutant | `CP107005.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706795` (ILMN) / `SRR21706791` (PacBio) |
| `SAMN31007784` | `SY22A3` | CRISPR-tolerant escaper mutant | `CP107004.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706794` (ILMN) / `SRR21706790` (PacBio) |
| `SAMN31007787` | `SY22A4` | CRISPR-tolerant escaper mutant | `CP106999.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706777` (ILMN) / `SRR21706787` (PacBio) |
| `SAMN31007792` | `SY22A6` | CRISPR-tolerant escaper mutant | `CP106994.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706772` (ILMN) / `SRR21706781` (ONT) |
| `SAMN31007794` | `SY22B1` | CRISPR-tolerant escaper mutant | `CP107003.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706792` (ILMN) / `SRR21706779` (ONT) |
| `SAMN31007786` | `SY22B2` | CRISPR-tolerant escaper mutant | `CP107000.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706778` (ILMN) / `SRR21706788` (PacBio) |
| `SAMN31007791` | `SY22B3` | CRISPR-tolerant escaper mutant | `CP106995.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706773` (ILMN) / `SRR21706782` (ONT) |
| `SAMN31007785` | `SY22B5` | CRISPR-tolerant escaper mutant | `CP107001.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706783` (ILMN) / `SRR21706789` (PacBio) |
| `SAMN31007789` | `SY22B6` | CRISPR-tolerant escaper mutant | `CP106997.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706775` (ILMN) / `SRR21706785` (PacBio) |
| `SAMN31007793` | `SY22C3` | CRISPR-tolerant escaper mutant | `CP107002.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706793` (ILMN) / `SRR21706780` (ONT) |
| `SAMN31007790` | `SY22D8` | CRISPR-tolerant escaper mutant | `CP106996.1` | 4,687,476 bp | 1 (Chromosome only) | `SRR21706774` (ILMN) / `SRR21706784` (PacBio) |

---

## 3. Candidate Truth & Eligibility Matrix

To rigorously distinguish evidence categories, each potential target locus is evaluated under four epistemological tiers:
- **`ASSEMBLY_SUPPORTED`**: Confirmed by closed, high-quality, reciprocal whole-genome assembly alignment.
- **`EXPERIMENTALLY_VALIDATED`**: Supported by targeted physical laboratory assays (PCR amplicon, Sanger sequencing).
- **`COMPUTATIONAL_CANDIDATE`**: Inferred by algorithmic heuristics without physical confirmation.
- **`NOT_EVALUABLE`**: Lacking essential sequence records or reference assets.

| Target Sample | BioSample | SRA Run (PE) | Assembly Accn | Target Replicon | Reference Accn | Evidence Source | Variant Type | Reference Coordinates | Alternate Coordinates | IS Family | TSD Length | Biological Relationship | Truth Eligibility Classification |
| :--- | :--- | :--- | :---: | :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`SY22A1`** | `SAMN31007783` | `SRR21706795` | `CP107005.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22A1`** | `SAMN31007783` | `SRR21706795` | *Not Deposited* | Cas-Plasmid | *Not Deposited* | PCR / Sanger (*Nat Commun*) | MOB (IS1 / IS10) | *Unmapped* | *Unmapped* | IS1 / IS4 | 9 bp / 9 bp | Isolate vs Uninduced Control | **`NOT_EVALUABLE`** (No plasmid assembly) |
| **`SY22A3`** | `SAMN31007784` | `SRR21706794` | `CP107004.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22A4`** | `SAMN31007787` | `SRR21706777` | `CP106999.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22A6`** | `SAMN31007792` | `SRR21706772` | `CP106994.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22B1`** | `SAMN31007794` | `SRR21706792` | `CP107003.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22B2`** | `SAMN31007786` | `SRR21706778` | `CP107000.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22B3`** | `SAMN31007791` | `SRR21706773` | `CP106995.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22B5`** | `SAMN31007785` | `SRR21706783` | `CP107001.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22B6`** | `SAMN31007789` | `SRR21706775` | `CP106997.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22C3`** | `SAMN31007793` | `SRR21706793` | `CP107002.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| **`SY22D8`** | `SAMN31007790` | `SRR21706774` | `CP106996.1` | Chromosome | `CP106998.1` | NCBI Assembly Alignment | **NONE** (0 variants) | None | None | None | None | Isolate vs Uninduced Control | **`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** |
| *Fabricated* | `SAMN31007783` | `SRR21706795` | `CP107005.1` (Truncated) | Chromosome | `CP106998.1` | Truncated fasta + MUMmer | 4 MOB, 7 JC (`SY22A1_truth.gd`) | 4,492,786–4,687,476 | Truncation tail | IS1, IS10 | 9 bp | File Truncation Artifact | **`FABRICATED_COMPUTATIONAL_ARTIFACT`** |

---

## 4. Formal Acceptance Decision

### Decision: **C. `ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`** (and **D. `INSUFFICIENT_GROUND_TRUTH`** for Positive SV Calibration)

#### Justifications:
1. **Rejection of Option A (`ELIGIBLE_FOR_CHROMOSOMAL_SV_TRUTH`):**
   Claiming that PRJNA884016 contains positive chromosomal structural variants would constitute scientific falsification. The true assemblies in NCBI demonstrate 100.000% identity across all 4,687,476 bases.
2. **Rejection of Option B (`ELIGIBLE_FOR_PLASMID_SV_TRUTH`):**
   No plasmid assemblies or reference accessions were deposited in NCBI under PRJNA884016. Assembly-confirmed plasmid truth is impossible without deposited sequences.
3. **Acceptance of Option C (`ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY`):**
   Because the complete circular chromosomes of `SY22A1` and `SY221` are bit-identical, PRJNA884016 provides a clean benchmark for **Negative Specificity / False Positive Resistance** (verifying that differential callers emit zero unintended calls between identical genomic backgrounds).
4. **Acceptance of Option D (`INSUFFICIENT_GROUND_TRUTH` for Positive Recall):**
   Because there are zero true positive chromosomal structural variants, PRJNA884016 cannot measure sensitivity, recall, breakpoint tolerance, or TSD concordance.

---

## 5. Recommended Alternative Public Datasets for Positive Structural Truth

To rigorously calibrate positive structural variant calling (MOB, JC, deletions, TSD resolution), ProkaDiff must employ datasets with verified, published chromosomal structural mutations:

1. ***E. coli* B REL606 Long-Term Evolution Experiment (LTEE) Isolates (Milestone M8):**
   - Evolved isolates (Tenaillon et al. 2016 *Nature*, Blount et al. 2012 *Nature*), including `SRR030257` (20,000 generations) already pre-registered in Line 5 of [`VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv).
   - Contains authentic chromosomal IS transpositions (IS150, IS1, IS186) and large chromosomal deletions with complete hybrid assemblies.
   - Primary reference benchmark used by the Barrick laboratory for breseq structural validation.
2. **Consuegra et al. 2021 (*Nat Commun* 12, 980):**
   - High-throughput short- and long-read resequencing of *E. coli* undergoing active chromosomal IS mobilization, with Sanger-verified transposition coordinates.
3. **Synthetic Dual-Strain Ground Truth Fixtures (Layer 1 / Layer 2):**
   - Fully controllable synthetic reads generated from edited reference genomes with known single-base breakpoint truth.

---

## 6. Recommended M7 Execution Scope & Protocol

1. **Retire Fabricated Truth Asset:**
   - Immediately deprecate and remove [`benchmark/real_world/truth/prjna884016/SY22A1_truth.gd`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/truth/prjna884016/SY22A1_truth.gd).
   - Replace it with a verified **EMPTY** GenomeDiff ground truth file representing 0 chromosomal mutations.
2. **M7-P1 Execution Scope (Negative Specificity Gate):**
   - Execute ProkaDiff and breseq 0.40.1 on Slurm partition `qcpu_18i` comparing `SRR21706795` (SY22A1) against `SRR21706776` (SY221) aligned to `CP106998.1`.
   - Pass Criteria: The differential caller must emit **0 true chromosomal calls**, demonstrating resistance to false positive structural calls on identical backgrounds.
3. **Validation Registry Invariant:**
   - Keep line 7 of [`benchmark/real_world/VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv) as **`PENDING`** with all validation booleans `false`.
   - Update `truth_level` from `ASSEMBLY_CONFIRMED` to `NEGATIVE_CONTROL_ASSEMBLY_CONFIRMED`.
## 7. Independent Codex Review Verdict

- **Reviewer:** Codex (GPT-6.1, Mode: `READ_ONLY`), Independent Verification and Acceptance Reviewer
- **Formal Verification Verdict:** **`M7_P0_VERIFICATION_PASSED`**
- **Governance Review Verdict:** **`PASS_WITH_LIMITATIONS`**
- **Auditor Conclusions:**
  1. *100.000% Chromosome Bit-Identity:* Reconciled with Sheng et al. (2023), confirming zero chromosomal mutations across all 12 NCBI chromosome assemblies (`CP106994.1`–`CP107005.1` all identical to `CP106998.1`).
  2. *Plasmid Localization:* Confirmed that genuine IS1 and IS10 transpositions occurred into plasmid-borne *cas* genes, with zero plasmid sequences deposited in PRJNA884016.
  3. *Truth Fabrication Forensic:* Proved that `SY22A1_truth.gd` was an artifact of an artificial 194.7 kb file truncation (`4,492,785 bp` vs `4,687,476 bp`).
  4. *Acceptance Classification:* Fully concurred with `ELIGIBLE_FOR_NEGATIVE_CONCORDANCE_ONLY` and `INSUFFICIENT_GROUND_TRUTH` for positive structural recall.
  5. *Alternative Public Datasets:* Endorsed pivoting positive structural calibration to *E. coli* B REL606 LTEE evolved isolates (Milestone M8) and Consuegra et al. (2021).
  6. *Boundaries Respected:* Zero Slurm jobs submitted; production engine `crates/**` 100% frozen; `VALIDATION_REGISTRY.tsv` line 7 preserved as `PENDING`.

---

## 8. Summary of Findings & Next Phase Recommendation

Milestone M7-P0 is complete and conclusively documented. BioProject PRJNA884016 is scientifically eligible **only** as a gold-standard negative control for false-positive resistance (target: 0 calls). Formal deprecation of `SY22A1_truth.gd` and scheduling of M7-P1 Slurm negative concordance smoke testing await human review and authorization.
