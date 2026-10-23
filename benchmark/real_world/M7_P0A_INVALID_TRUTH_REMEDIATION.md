# Milestone M7-P0a: Invalid Assembly Truth Quarantine & Remediation Dossier

**Milestone:** M7-P0a — Invalid Assembly Truth Quarantine
**Date:** 2026-10-23
**Status:** **`M7_P0A_INVALID_TRUTH_QUARANTINED`**
**Repository Baseline:** `5c90e0568153a4befc573933d2f98bd3d47662b3` (Local HEAD)
**Governance Invariants:** Crate engines `crates/**` 100% frozen; Zero Slurm jobs submitted; Zero Git push; Registry Line 6 (`prjna1088182`) untouched; Registry Line 7 booleans strictly `false` and status `PENDING`.

---

## Executive Summary

During Milestone M7 preparation, forensic verification of BioProject `PRJNA884016` (*E. coli* strain `SY22A1` vs `SY221`) revealed that the versioned assembly truth ([`SY22A1_truth.gd`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/truth/prjna884016/SY22A1_truth.gd)), reporting 4 mobile element insertions (MOBs) and 7 breakpoint junctions (JCs), was a **computational artifact caused by a 194.7 kb local file truncation** (`4,492,785 bp` vs expected `4,687,476 bp`).

In physical reality, the official complete circular NCBI chromosome assemblies `CP107005.1` (SY22A1) and `CP106998.1` (SY221) are **100.000% bit-identical across all 4,687,476 bp (`diff_count = 0`, zero SNPs, zero indels, zero structural variants)**. The genuine biological transpositions reported by Sheng et al. (2023) occurred exclusively on the extrachromosomal *cas*-expression plasmid, which was never deposited in NCBI.

In Milestone M7-P0a, the ProkaDiff benchmark framework has executed a fail-closed truth quarantine:
1. **Forensic Preservation:** All original invalid truth files were hashed and archived to [`benchmark/real_world/truth/prjna884016/quarantined_legacy/`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/truth/prjna884016/quarantined_legacy/) with a cryptographic manifest.
2. **Asset Quarantine:** Active truth files (`SY22A1_truth.gd`, `mob_truth.tsv`, `jc_truth.tsv`, `SY22A1_vs_SY221.rdiff`) were replaced with empty zero-variant records (0 MOB, 0 JC).
3. **Truth Generator Disabled:** [`build_prjna884016_truth.py`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/scripts/build_prjna884016_truth.py) was disabled for authoritative validation and retrofitted with fail-closed checks against truncation, checksum mismatch, and heuristic length-based MOB/TSD fabrication.
4. **Scoring Hardening:** [`compare_mob_truth.py`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/scripts/compare_mob_truth.py) was updated to return `NOT_EVALUABLE` for undefined denominators rather than defaulting to 100%.
5. **Registry Rectification:** In [`VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv) Line 7, `truth_level` was retained as `ASSEMBLY_CONFIRMED` under Claude Opus 5.5 architectural guidance (reflecting assembly confirmation of zero chromosomal variants), all booleans remain `false`, status remains `PENDING`, and notes explicitly repudiate the truncation artifact.
6. **Regression Invariants:** 9 fail-closed regression tests were implemented in [`test_prjna884016_quarantine.py`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/scripts/test_prjna884016_quarantine.py), passing with 100% non-regression across all M6 suites and Rust workspace crates.

---

## A. Scientific Root-Cause Analysis

### 1. The Biological Context (Sheng et al., 2023)
The study (*Nat Commun* 14:4366, 2023, PMC10359306) investigated CRISPR-tolerant escaper mutants in *E. coli* MG1655 transformed with a Cas-expression plasmid (`pCasY-3-ΔSIM-kan`). Sheng et al. explicitly stated:
> *"Sequence alignment of SY221 with the CRISPR-tolerant mutants revealed 100% sequence identity at 100% chromosome coverage, suggesting that the transpositions of ISs into cas genes did not result in simultaneous transpositions into other genomic loci or trigger genomic rearrangements."*

The paper documented genuine IS1 and IS10 transpositions, but these occurred exclusively into the *cas* genes on the **extrachromosomal plasmid**.

### 2. Missing Plasmid Depositions in NCBI
BioProject PRJNA884016 comprises 12 BioSamples (`SAMN31007783`–`SAMN31007794`). For each sample, the authors deposited a single circular replicon: the complete chromosome (`CP106994.1`–`CP107005.1`). **Zero plasmid assemblies were deposited in NCBI**. Consequently, no positive assembly truth for IS transpositions exists in the public repository records.

### 3. Forensic Origin of the Erroneous 4 MOB + 7 JC Truth
On September 18, 2026, during initial validation setup, the local mutant assembly file [`data/prjna884016/assemblies/SY22A1.fna`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna884016/assemblies/SY22A1.fna) was staged with a size of 4,549,015 bytes (56,161 lines, exactly 4,492,785 bp). The file ended abruptly on line 56,161 with 65 characters and no trailing newline, demonstrating an incomplete stream truncation during network retrieval or disk extraction (the uninduced control `SY221_ref.fna` was retrieved 5 minutes earlier and was intact at 4,687,476 bp).

When [`build_prjna884016_truth.py`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/scripts/build_prjna884016_truth.py) executed MUMmer `dnadiff` between truncated `SY22A1.fna` and full `SY221_ref.fna`:
1. MUMmer detected the unaligned terminal 194,691 bp of the reference chromosome (coordinates 4,492,786 to 4,687,476) and identified repeat duplications (`DUP`) and breakpoints (`BRK`) within this unaligned tail.
2. The legacy script heuristically mapped any `DUP` with length 760–780 bp to `IS1` and 1320–1340 bp to `IS10`.
3. It fabricated hardcoded strand `+` and fixed `9bp` target site duplications (TSDs).
4. This resulted in 4 false MOBs (two at 4,499,357–4,499,358 and two at 4,630,529–4,630,531) and 7 false JCs, all residing strictly within the truncated tail coordinates.
5. In physical reality and in NCBI accession `CP107005.1`, no such truncation or structural variants exist.

---

## B. Original and Canonical Sequence Checksums

### 1. Canonical NCBI Chromosome Assemblies (All 12 Strains)
All 12 deposited assemblies under PRJNA884016 (`GCA_030345415.1` – `GCA_030345635.1`) are complete, finished circular chromosomes of identical length:

| Strain Name | BioSample | NCBI Accession | Replicon Type | Length (bp) | Normalized Sequence SHA256 | Identity vs CP106998.1 |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **SY221** (Control) | `SAMN31007788` | `CP106998.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | Reference (100.000%) |
| **SY22A1** (Mutant) | `SAMN31007783` | `CP107005.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22A3** | `SAMN31007784` | `CP107004.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22A4** | `SAMN31007787` | `CP106999.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22A6** | `SAMN31007792` | `CP106994.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22B1** | `SAMN31007794` | `CP107003.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22B2** | `SAMN31007786` | `CP107000.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22B3** | `SAMN31007791` | `CP106995.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22B5** | `SAMN31007785` | `CP107001.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22B6** | `SAMN31007789` | `CP106997.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22C3** | `SAMN31007793` | `CP107002.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |
| **SY22D8** | `SAMN31007790` | `CP106996.1` | Chromosome | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | **100.000% (`diff=0`)** |

*Note: Normalized sequence SHA256 is computed across uppercase single-line sequence strings, invariant to FASTA header text or line-wrapping.*

### 2. Forensic Inspection of Local Files

| File Path | Raw File Bytes | Line Count | Sequence Length (bp) | Normalized Sequence SHA256 | Forensic State |
| :--- | :---: | :---: | :---: | :---: | :--- |
| `data/prjna884016/assemblies/SY221_ref.fna` | 4,746,140 | 58,595 | 4,687,476 | `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a` | Complete & Intact |
| `data/prjna884016/assemblies/SY22A1.fna` | 4,549,015 | 56,161 | 4,492,785 | `673fa32a0c1ee516200a2e5d5b80285cc145120b408c2441d54cc9df1cb84ff5` | **Truncated (-194,691 bp)** |

Sequence prefix check proves that the 4,492,785 bp of local `SY22A1.fna` is an **exact 100.000% strict prefix** of `SY221_ref.fna`. Line 56,161 was truncated mid-stream at 65 characters without a newline.

---

## C. Invalidated Assets and Their Dependencies

### 1. Cryptographic Record of Quarantined Assets
Before modification, all legacy files were recorded and preserved in [`benchmark/real_world/truth/prjna884016/quarantined_legacy/`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/truth/prjna884016/quarantined_legacy/):

| Quarantined Asset | Original Path | Original SHA256 Checksum | Bytes | Quarantine Classification |
| :--- | :--- | :---: | :---: | :--- |
| `SY22A1_truth_invalid_legacy.gd` | `truth/prjna884016/SY22A1_truth.gd` | `438925b86e2fe6b54483361796e50bae200872418c6b7117568ffd87065076f0` | 611 | `QUARANTINED_INVALID_TRUNCATION_ARTIFACT` |
| `mob_truth_invalid_legacy.tsv` | `truth/prjna884016/mob_truth.tsv` | `b69047136cb15017ddc23c32b06bca695c110160aef267315e9f604b468ac00d` | 285 | `QUARANTINED_INVALID_TRUNCATION_ARTIFACT` |
| `jc_truth_invalid_legacy.tsv` | `truth/prjna884016/jc_truth.tsv` | `e8f17766f31057f78145d3b744ea227e2ed8bacbd2f783f2e2ce245ecd8ba0fb` | 440 | `QUARANTINED_INVALID_TRUNCATION_ARTIFACT` |
| `SY22A1_vs_SY221_legacy.rdiff` | `truth/prjna884016/SY22A1_vs_SY221.rdiff` | `83c872a7e4bb49357499d3b9a02e0f5a0683c26a1eb77028f4ff8bd0ef148d76` | 716 | `QUARANTINED_INVALID_TRUNCATION_ARTIFACT` |
| `SY22A1_vs_SY221_legacy.report` | `truth/prjna884016/SY22A1_vs_SY221.report` | `3716e57c08a9c6d5129af03f23b174cfa11c8cf74d770ce09ab6b42109a757f0` | 4,286 | `QUARANTINED_INVALID_TRUNCATION_ARTIFACT` |
| `SY22A1_vs_SY221_legacy.delta` | `truth/prjna884016/SY22A1_vs_SY221.delta` | `dbc542462ef53fc634cf059f9df7ced66e2582dd9651f33a1af8e7ce7d7fbde9` | 69,536 | `QUARANTINED_INVALID_TRUNCATION_ARTIFACT` |

### 2. State of Active Truth Discovery Paths
1. [`SY22A1_truth.gd`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/truth/prjna884016/SY22A1_truth.gd): Replaced with an explicit zero-variant GenomeDiff header (0 MOB, 0 JC).
2. [`mob_truth.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/truth/prjna884016/mob_truth.tsv): Header only (0 data rows).
3. [`jc_truth.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/truth/prjna884016/jc_truth.tsv): Header only (0 data rows).
4. [`SY22A1_vs_SY221.rdiff`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/truth/prjna884016/SY22A1_vs_SY221.rdiff): Empty file (0 lines).
5. [`check_validation_assets.sh`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/scripts/check_validation_assets.sh): Updated lines 39-42 to assert that MOB count is 0, JC count is 0, and data rows in TSVs are 0.

---

## D. Registry Correction Decision

### 1. Scientific Architectural Analysis (Claude Opus 5.5)
The existing vocabulary for `truth_level` across the codebase is `{MANUALLY_CURATED, BRESEQ_ORACLE, PUBLICATION_REPORTED, ASSEMBLY_CONFIRMED}`.
- `truth_level` specifies the **methodological tier/evidence authority** of the ground truth artifact, NOT its mathematical polarity (positive vs negative).
- Complete circular finished assembly is the gold standard for confirming **zero-variant concordance** on the chromosome.
- The positive claim was not inherent to the string `ASSEMBLY_CONFIRMED`, but to the corrupted files and narrative notes.
- Inventing a new enum (such as `NEGATIVE_CONTROL_ASSEMBLY_CONFIRMED`) was strictly prohibited by project rules.

### 2. Authoritative Rectification of Line 7
Line 7 of [`benchmark/real_world/VALIDATION_REGISTRY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/VALIDATION_REGISTRY.tsv#L7) was updated:
```tsv
prjna884016	SY22A1	ASSEMBLY_CONFIRMED	false	false	false	false	false	PENDING	Assembly confirms 0 chromosomal variants between CP107005.1 and CP106998.1 (100.000% bit-identical). Historical 4 MOB + 7 JC invalidated as local 194.7 kb truncation artifact. Eligible for negative concordance / zero-call specificity audit only. Positive SV sensitivity NOT_EVALUABLE on chromosome; deferred to REL606. Keep all fields false and status PENDING.
```

- **Validation Flags:** `false false false false false` (Strictly maintained).
- **Status:** `PENDING` (Strictly maintained).
- **Line 6 (`prjna1088182`):** **Strictly unmodified**.

---

## E. Regression Tests

A dedicated fail-closed test suite was implemented in [`test_prjna884016_quarantine.py`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/scripts/test_prjna884016_quarantine.py):

| Test Name | Verified Contract Invariant | Result |
| :--- | :--- | :---: |
| `test_truncated_assembly_fails_closed` | `validate_assembly_fasta` raises `AssemblyTruncationError` on truncated file | **PASS** |
| `test_sequence_checksum_mismatch_fails_closed` | Raises `SequenceChecksumMismatchError` on non-canonical sequence | **PASS** |
| `test_replicon_mismatch_fails_closed` | Raises `RepliconMismatchError` on incorrect accession in FASTA header | **PASS** |
| `test_heuristic_mob_family_and_tsd_rejected` | Raises `UnsupportedHeuristicError` when inferring IS1/IS10 from segment length | **PASS** |
| `test_empty_truth_scoring_returns_not_evaluable` | Zero truth events yield `NOT_EVALUABLE` for precision, recall, and tolerances | **PASS** |
| `test_legacy_corrupted_truth_records_rejected` | Confirms active `SY22A1_truth.gd` has 0 MOB, 0 JC, and TSVs have 0 data rows | **PASS** |
| `test_quarantined_legacy_assets_preserved_and_hashes_intact` | Validates `quarantined_legacy/MANIFEST.tsv` and all archived file checksums | **PASS** |
| `test_validation_registry_invariants` | Verifies Line 7 schema invariants, `PENDING` status, `false` booleans, Line 6 untouched | **PASS** |
| `test_manifest_prjna884016_notes` | Verifies `manifests/prjna884016.tsv` declares negative concordance and repudiates transpositions | **PASS** |

### Workspace Verification Suite
- **Pytest:** `test_prjna884016_quarantine.py` (9 passed in 0.26s).
- **M6 Non-Regression:** `test_m6_p2b_audit.py`, `test_m6_p3a_read_level.py`, `test_m6_p3b_reference.py`, `test_plan_m6_p3a.py` (18 passed in 5.81s).
- **Asset Gate:** `bash benchmark/real_world/scripts/check_validation_assets.sh` (PASSED).
- **Rust Workspace:** `cargo test --workspace` (All crates passed: 17 cli, 15 cli_e2e, 8 m5_sr1, 127 classify, 20 offtarget, 18 report, 10 fixtures, 10 m4_contract; 0 failed).

---

## F. Independent Codex Verdict

- **Reviewer:** Codex (GPT-6.1 HIGH, Mode: `READ_ONLY`), Independent Verification and Acceptance Reviewer
- **Date:** 2026-10-23
- **Formal Review Verdict:** **`PASS_WITH_LIMITATIONS`**
- **Milestone Status Assessment:** **`M7_P0A_ACCEPTANCE_PASSED`**
- **Auditor Findings:**
  1. *Root-Cause Forensic Accuracy Confirmed:* Bit-identity of NCBI circular chromosome assemblies `CP107005.1` and `CP106998.1` confirmed across 4,687,476 bp (`diff_count = 0`). Local file `SY22A1.fna` confirmed truncated at 4,492,785 bp, producing the false 4 MOB + 7 JC artifact across the missing 194.7 kb tail.
  2. *Quarantine and Cryptographic Preservation:* Quarantined legacy assets cataloged with byte counts and SHA256 checksums in `quarantined_legacy/MANIFEST.tsv`.
  3. *Truth Assets Repudiation:* Active truth files correctly replaced with 0-variant GenomeDiff header and empty tables (0 MOB, 0 JC).
  4. *Validation Registry Integrity:* Line 7 correctly retains `ASSEMBLY_CONFIRMED` under Claude Opus 5.5 architectural guidance (assembly confirms 0 chromosomal mutations), booleans strictly `false`, status strictly `PENDING`, notes updated. Line 6 (`prjna1088182`) strictly untouched.
  5. *Fail-Closed Tooling & Scoring:* `build_prjna884016_truth.py` disabled by default with strict fail-closed guards. `compare_mob_truth.py` returns `NOT_EVALUABLE` on zero truth events. `check_validation_assets.sh` passes asserting 0 MOB and 0 JC.
  6. *Automated Testing:* 9/9 quarantine regression tests passed; 18/18 M6 tests passed; 231/231 Rust workspace tests passed. Engine `crates/**` 100% frozen. Zero Slurm jobs submitted. Zero Git push.

---

## G. Corrective Commit

- **Commit Hash:** `894a613a0e107577579737dc1d5bf5967f892c9f` (local commit)
- **Commit Message:** `fix(benchmark): quarantine invalid PRJNA884016 structural truth and enforce fail-closed invariants (M7-P0a)`
- **Affected Files:** 23 files changed, 8,175 insertions(+), 236 deletions(-)
- **Push Policy:** Strictly local (zero push); human authorization required prior to remote synchronization.

---

## H. Recommended Eligible Scope for M7

1. **Negative Concordance Benchmark (PRJNA884016):**
   - BioProject PRJNA884016 (`SY22A1` vs `SY221`) is scientifically eligible **strictly as a zero-variant negative control** to evaluate ProkaDiff's specificity and false-positive resistance (target: 0 calls on chromosome).
2. **Positive Structural Benchmark Deferred to REL606 / M8:**
   - Because the genuine IS transpositions in PRJNA884016 occurred on unsequenced plasmids, PRJNA884016 chromosome WGS cannot evaluate positive MOB recall, JC breakpoint precision, or TSD accuracy (strictly `NOT_EVALUABLE`).
   - Positive structural variant calibration must pivot to *E. coli* B REL606 LTEE evolved isolates (e.g. `SRR030257` with 20,000 generations of accumulated chromosomal IS insertions and large deletions) or synthetic fixtures under Milestone M8.
