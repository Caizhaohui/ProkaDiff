# M6-P3c Structural Computational Evidence Assessment Plan

**Milestone:** `M6-P3c` (Structural Computational Corroboration & Repeat Ambiguity Audit)
**Status:** **PLAN READY FOR CODEX INDEPENDENT DESIGN REVIEW**
**Core Codebase:** Frozen (`crates/**` read-only).
**HPC Compute Target:** Slurm partition `qcpu_18i`.

---

## 1. Executive Summary & Objective

Phase M6-P2a established that structural event reporting in PRJNA1088182 varies dramatically across strains and biological conditions:
- **MG1655 cohort (10 runs against `NC_000913.3`):** 45 to 505 JC records per run (total: 1,230 JCs, mean: 123.00), with 2 to 3 MOB records per run (total: 23 MOBs, mean: 2.30).
- **BW25113 cohort (15 runs against `CP009273.1`):** 35 to 72 JC records per run (total: 719 JCs, mean: 47.93), with 0 MOB records per run (`MOB_NOT_CALLED`).

The primary objective of Phase M6-P3c is to establish a reproducible, independent computational evidence assessment for reported structural events in the PRJNA1088182 benchmark, addressing:
1. How many reported JC events have reproducible read-level support?
2. How should repeat-associated and repeat-ambiguous junctions be classified?
3. Can independently generated computational evidence corroborate reported structural events?
4. Which MOB calls have sufficient mechanistic evidence for a candidate classification?
5. Which results remain biologically `NOT_EVALUABLE` without orthogonal experimental validation?

---

## 2. Primary Audit Target: `SRR28349173` (`Cas9_110_111`)

The primary audit target for structural evidence evaluation is `SRR28349173` (`Cas9_110_111`, MG1655 against `NC_000913.3`), which exhibits the highest structural complexity in the PRJNA1088182 cohort:
- **505 JC records** in `output.gd` (IDs 91 through 596)
- **2 MOB records** in `output.gd` (IDs 597 & 598)
- **4 SNP records** in `output.gd` (IDs 1 through 4) + 1 INS record (ID 5)
- **85 MC records** in `output.gd` (IDs 6 through 90)
- **34,450 repeat-ambiguity diagnostic entries** in `repeat_ambiguous_junctions.tsv` (34,447 `CANDIDATE` and 3 `SUPPORTED_AMBIGUOUS`)

This target is designated as a priority evidence-audit target to evaluate structural sensitivity and repeat ambiguity under extreme structural load, not as a predefined false-positive sample.

---

## 3. Scientific Design Requirements (A through L)

### A. Structural Evidence Taxonomy

Every structural record emitted or analyzed is classified into a single, mutually exclusive taxonomic category:

1. `MOB_PROMOTED`: A mobile element insertion event supported by a reciprocal junction pair, matching an annotated IS element boundary ($\le 20$ bp), opposite strands, target distance $\le 10$ bp, and target site duplication ($dup \ge 3$ bp).
2. `DEL_CORROBORATED`: A genomic deletion supported by both a high-confidence deletion junction (`overlap < 0` or collinear $p_1 < p_2$) and an overlapping missing-coverage (MC) span.
3. `JC_UNIQUE_ANCHORED`: A junction where both breakpoints land in non-repeat genomic sequence with high mapping confidence ($AS > XS$, MAPQ $\ge 10$, flank overlap $\ge 14$ bp on both sides).
4. `JC_REPEAT_ANCHORED`: A junction where one breakpoint is unique genomic sequence and the other lands at an annotated repeat boundary, but lacks a matching reciprocal partner for MOB promotion.
5. `JC_INTER_REPEAT`: A junction where both breakpoints land inside annotated reference repeats. (Treated as putative mapping artifacts between multi-copy elements and suppressed when repeat annotations are provided).
6. `REPEAT_AMBIGUOUS_DIAGNOSTIC`: Reads whose split or softclip alignments map to $\ge 2$ reference repeat copies without bridging unique sequence (retained in `repeat_ambiguous_junctions.tsv`, never emitted as mutation calls).

### B. Distinction Between JCs, Biological Events, Ambiguity, and MOBs

- **Junction Call (JC):** An alignment-level hypothesis reporting that sequence from $(c_1, p_1, s_1)$ connects to $(c_2, p_2, s_2)$. A single physical biological event often generates multiple JCs (e.g., a transposon insertion produces two flanking JCs; an inversion produces two reciprocal JCs).
- **Biological Structural Event:** The physical genomic alteration in the chromosome.
- **Repeat Ambiguity:** An epistemic limitation of short reads ($L \le 151$ bp) when spanning repeats longer than the read length.
- **MOB Promotion:** An evidence-layer synthesis that resolves two complementary JCs into a unified biological transposition call with an explicit target site duplication (TSD).

### C. Read-Level Support and Mapping Quality Criteria

A junction is classified as having **Robust Read-Level Support** if and only if:
1. Total supporting reads $N_{\text{reads}} \ge 10$ (`pd_support_reads`).
2. Bilateral strand support: $N_{\text{plus}} \ge 2$ AND $N_{\text{minus}} \ge 2$ (`pd_plus_reads`, `pd_minus_reads`).
3. Spanning flank depth: $min\_overlap\_side1 \ge 14$ bp AND $min\_overlap\_side2 \ge 14$ bp.
4. Anchor uniqueness: Alignments covering the junction breakpoint do not have $XS \ge AS$ across the entire breakpoint window.

### D. Handling of Repeated Elements and Microhomology

1. **Microhomology Sliding:** Junction alignments often exhibit microhomology (1–5 bp) at the breakpoint. Breseq and ProkaDiff resolve this by reporting the leftmost coordinate or right-aligning across the microhomology window (`JC_CLUSTER_TOL_BP = 5`).
2. **Boundary Displacement Adjustment:** When a junction breakpoint on an IS element is displaced by $\le 5$ bp from the annotated repeat boundary due to terminal microhomology, ProkaDiff applies an equal-and-opposite displacement to the target breakpoint.
3. **Multi-Copy Placement Folding:** Split reads mapping to identical repeat copies are collapsed into a single candidate with a multi-copy placement flag, preventing artificial multiplication of junction candidates.

### E. Representative Sample Selection Strategy (The Benchmark Quintet)

To bound compute while maintaining comprehensive structural coverage, Phase M6-P3c selects a representative 5-run cohort from PRJNA1088182:

| Sample Name | Run Accession | Strain / Ref | P2a JCs | P2a MOBs | Role & Evaluation Focus |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Cas9_110_111** | `SRR28349173` | MG1655 (`NC_000913.3`) | **505** | **2** | **Primary Audit Target:** Extreme structural complexity, 34k repeat diagnostics, 2 high-support MOBs. |
| **Cas9_104_105** | `SRR28349176` | MG1655 (`NC_000913.3`) | **140** | **3** | **Intermediate MG1655:** Multi-colony pool, verified P3a SNV target, moderate JC count. |
| **Cas9_1** | `SRR28349244` | MG1655 (`NC_000913.3`) | **49** | **3** | **Baseline MG1655:** Reused P1b baseline, standard clonal profile. |
| **I-SceI_3_4_5_6** | `SRR28349231` | BW25113 (`CP009273.1`) | **72** | **0** | **Max BW25113 Complexity:** Positive endonuclease DSBs, highest BW25113 JC count. |
| **control_1_2_3_4** | `SRR28349226` | BW25113 (`CP009273.1`) | **40** | **0** | **Negative Control:** Negative parental baseline, false-positive control. |

### F. Independent Computational Corroboration Protocol (breseq Oracle)

Per AGENTS.md Rule 4, 5, and 12:
1. `breseq 0.40.x` is invoked strictly as an offline computational oracle. It is never called as a runtime engine.
2. Binary: `breseq` version pinned in `testdata/VERSIONS.txt`.
3. Execution: Submitted to Slurm partition `qcpu_18i` with identical inputs:
   ```bash
   breseq -j 8 -r <reference.gbk> -o <outdir>/breseq <run>_1.fastq.gz <run>_2.fastq.gz
   ```
4. Output Extraction: Convert `output.gd` from breseq and ProkaDiff into standardized JC comparison tables:
   - Breakpoint match window: $|p_{1,A} - p_{1,B}| \le 5$ bp AND $|p_{2,A} - p_{2,B}| \le 5$ bp.
   - Strand concordance: $s_{1,A} == s_{1,B}$ AND $s_{2,A} == s_{2,B}$.

### G. Policy Distinguishing Computational Corroboration from Ground Truth

1. **Computational Corroboration $\ne$ Biological Ground Truth:**
   - Concordance between ProkaDiff and breseq proves that both algorithms independently identify the same chimeric sequence configuration in the raw reads.
   - If an Illumina library contains an in-vitro PCR chimera or optical artifact, both tools will predictably detect it. Computational agreement is therefore **evidence of algorithmic reproducibility, not wet-lab biological reality**.
2. **Ground Truth Requirements:** Biological ground truth strictly requires orthogonal, non-sequencing assays (e.g. PCR amplification spanning the breakpoint, Sanger sequencing across the insertion, or PacBio/ONT long-read single-molecule validation).

### H. Evidence-Tier Classification System

Each structural event is assigned to an evidence tier:

| Tier | Classification Label | Operational Criteria |
| :---: | :--- | :--- |
| **Tier 1** | `CORROBORATED_CANDIDATE` | Supporting reads $\ge 10$, bilateral strand support ($+ \ge 2, - \ge 2$), verified by independent oracle (breseq $\pm 5$ bp), unique genomic flanks. |
| **Tier 2** | `ALGORITHMICALLY_CONSISTENT` | ProkaDiff support $\ge 10$, bilateral strands, unique flanks, but uncalled by oracle due to different frequency/coverage cutoffs. |
| **Tier 3** | `REPEAT_AMBIGUOUS` | Junction lands in multi-copy repeats; candidate placements cannot be resolved to a single physical locus. |
| **Tier 4** | `ARTEFACT_SUSPECT` | Extreme strand bias ($+ == 0$ or $- == 0$), support $< 10$ reads, low flank overlap ($< 14$ bp), or inter-repeat chimera without genomic anchors. |
| **Tier 5** | `BIOLOGICALLY_NOT_EVALUABLE` | Structural calls lacking starter-strain comparison, or repeat-internal variants lacking long-read resolution. |

### I. Scientifically Valid Acceptance Metrics

1. **Oracle Concordance on Unique High-Support JCs:** $\ge 80\%$ of Tier-1 ProkaDiff JCs must match breseq within $\pm 5$ bp.
2. **MOB Promotion Consistency:** 100% of promoted MOBs must possess valid parent JCs with canonical target site duplications ($dup \ge 3$ bp).
3. **Absorbed JC Accounting:** Every promoted MOB must absorb and tag exactly two constituent JCs (`mob_evidence=1`).
4. **Small-Variant Invariance:** SNP and indel calls must remain 100% invariant between standalone and oracle runs.

### J. Failure, Uncertainty, and NOT_EVALUABLE Criteria

A structural call must be declared `NOT_EVALUABLE` under any of the following conditions:
1. **Unsubtracted Reference Comparison:** Calls made without an isogenic starter strain cannot be evaluated for Cas9 causality.
2. **Repeat-Internal Chimeras:** JCs where both sides land within identical repeat families without external genomic anchors.
3. **Marginal Coverage ($< 5$ Reads):** Loci with insufficient depth to distinguish low-frequency variants from background Illumina error.

### K. Estimated HPC Resource & Storage Requirements

Per AGENTS.md Rule 11, execution is restricted to Slurm partition `qcpu_18i`:

- **Compute Allocation:**
  - Tasks: 5 representative runs.
  - CPUs per Task: 8 cores (`--cpus-per-task=8`).
  - Memory: 32 GB RAM (`--mem=32G`).
  - Wall Time: 04:00:00 per task (ProkaDiff takes $\sim 1.5$–3h; breseq oracle takes $\sim 2.5$–4h).
  - Array Concurrency: `--array=0-4%2` (at most 2 tasks running concurrently).
- **Storage Footprint:**
  - Fastq inputs: already staged in `data/prjna1088182/fastq/` (read-only, no re-download).
  - Output scratch per run: $\sim 6$ GB (BAMs, SAMs, and breseq intermediate files).
  - Total scratch required: $\sim 30$ GB (well within the $\ge 40$ GB storage gate).

### L. Rules for Preserving the Original Frozen Benchmark

1. `data/prjna1088182/` FASTQ and reference files are **STRICTLY READ-ONLY**.
2. P2a baseline results in `manifests/prjna1088182_m6_p2a_run_results.tsv` and `VALIDATION_REGISTRY.tsv` remain **COMPLETELY FROZEN**.
3. All P3c outputs must be written to isolated directories under `benchmark/real_world/results/prjna1088182/m6_p3c/<run>/`.
4. No modification of core engine code (`crates/**`) is permitted.
