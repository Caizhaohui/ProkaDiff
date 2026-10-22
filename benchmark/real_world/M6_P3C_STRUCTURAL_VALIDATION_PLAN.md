# M6-P3c Structural Computational Evidence Assessment Plan

**Milestone:** `M6-P3c` (Structural Computational Corroboration & Repeat Ambiguity Audit)
**Status:** **PLAN REVISED & HARDENED — READY FOR CODEX INDEPENDENT DESIGN REVIEW**
**Core Codebase:** Frozen (`crates/**` read-only).
**HPC Compute Target:** Slurm partition `qcpu_18i`.

---

## 1. Executive Summary & Objective

Phase M6-P2a established that structural event reporting in PRJNA1088182 varies dramatically across strains and biological conditions:
- **MG1655 cohort (10 runs against `NC_000913.3`):** 45 to 505 JC records per run (total: 1,230 JCs, mean: 123.00), with 2 to 3 MOB records per run (total: 23 MOBs, mean: 2.30).
- **BW25113 cohort (15 runs against `CP009273.1`):** 35 to 72 JC records per run (total: 719 JCs, mean: 47.93), with 0 MOB records per run (`MOB_NOT_CALLED` due to plain FASTA reference input, resolved in M6-P3b).

The primary objective of Phase M6-P3c is to establish a reproducible, independent computational evidence assessment for reported structural events in the PRJNA1088182 benchmark, addressing:
1. How many reported JC events have reproducible read-level support?
2. How should repeat-associated and repeat-ambiguous junctions be classified without confounding placement ambiguity with confirmed artifact status?
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

This target is designated as an evidence-audit benchmark to evaluate structural sensitivity, caller concordance, and repeat ambiguity under extreme structural load.

---

## 3. Scientific Design Requirements (A through L)

### A. Structural Evidence Taxonomy

Every structural record emitted or analyzed is classified into a single, mutually exclusive taxonomic category based on its topological and genomic configuration:

1. `MOB_PROMOTED`: A mobile element insertion event supported by a reciprocal junction pair, matching an annotated IS element boundary ($\le 20$ bp), opposite strands, target distance $\le 10$ bp, and canonical target site duplication ($dup \ge 3$ bp).
2. `DEL_CORROBORATED`: A genomic deletion supported by both a high-confidence deletion junction (`overlap < 0` or collinear $p_1 < p_2$) and an overlapping missing-coverage (MC) span.
3. `JC_UNIQUE_ANCHORED`: A junction where both breakends land in non-repeat genomic sequence with high mapping confidence ($AS > XS$, MAPQ $\ge 10$, flank overlap $\ge 14$ bp on both sides).
4. `JC_REPEAT_ANCHORED`: A junction where exactly one breakend lands in an annotated repeat boundary ($\le 20$ bp), while the other breakend is uniquely anchored in non-repeat sequence.
5. `JC_INTER_REPEAT`: A junction where both breakends land inside annotated reference repeats.
6. `REPEAT_AMBIGUOUS_DIAGNOSTIC`: Reads whose split or softclip alignments map to $\ge 2$ reference repeat copies without bridging unique sequence (retained in `repeat_ambiguous_junctions.tsv`, never emitted as mutation calls).

---

### B. Distinction Between JCs, Biological Events, Ambiguity, and MOBs

- **Junction Call (JC):** An alignment-level hypothesis reporting that sequence from $(c_1, p_1, s_1)$ connects to $(c_2, p_2, s_2)$. A single physical biological event often generates multiple JCs (e.g., an insertion produces two flanking JCs; an inversion produces two reciprocal JCs).
- **Biological Structural Event:** The physical genomic alteration in the chromosome.
- **Repeat Ambiguity:** An epistemic limitation of short reads ($L \le 151$ bp) when spanning repeats longer than the read length.
- **MOB Promotion:** An evidence-layer synthesis that resolves two complementary JCs into a unified biological transposition call with an explicit target site duplication (TSD).

---

### C. Read-Level Support and Mapping Quality Criteria

Operational classification of read-level support relies on pre-registered quantitative criteria:
1. **Total Supporting Reads:** $N_{\text{reads}} \ge 10$ (`pd_support_reads`).
2. **Bilateral Strand Support:** $N_{\text{plus}} \ge 2$ AND $N_{\text{minus}} \ge 2$ (`pd_plus_reads`, `pd_minus_reads`).
3. **Flank Extension Depth:** $min\_overlap\_side1 \ge 14$ bp AND $min\_overlap\_side2 \ge 14$ bp.
4. **Anchor Mapping Quality:** Uniquely mapped reads without $XS \ge AS$ across the breakpoint.
5. **Data Availability Standard:** If any specific metric (e.g., per-read Phred quality at breakpoint) is not natively recorded in the output GD or diagnostic TSV, it must be reported as `NOT_MEASURED` rather than estimated or invented.

---

### D. Handling of Repeated Elements and Microhomology

1. **Microhomology Sliding:** Junction alignments often exhibit microhomology (1–5 bp) at the breakpoint. Breseq and ProkaDiff resolve this by reporting the leftmost coordinate or right-aligning across the microhomology window (`JC_CLUSTER_TOL_BP = 5`).
2. **Boundary Displacement Adjustment:** When a junction breakpoint on an IS element is displaced by $\le 5$ bp from the annotated repeat boundary due to terminal microhomology, ProkaDiff applies an equal-and-opposite displacement to the target breakpoint.
3. **Multi-Copy Placement Folding:** Split reads mapping to identical repeat copies are collapsed into a single candidate with a multi-copy placement flag, preventing artificial multiplication of junction candidates.
4. **Non-Artifactual Nature of Repeats:** Repeat-associated and repeat-ambiguous junctions must **never** be labeled as "confirmed artefacts" without orthogonal experimental evidence. Genuine biological transpositions, homologous recombinations, and inversions naturally occur at repetitive loci.

---

### E. Representative Sample Selection Strategy (The Benchmark Quintet)

To bound compute while maintaining comprehensive structural coverage, Phase M6-P3c selects a representative 5-run cohort from PRJNA1088182:

| Sample Name | Run Accession | Strain / Ref | P2a JCs | P2a MOBs | Role & Evaluation Focus |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Cas9_110_111** | `SRR28349173` | MG1655 (`NC_000913.3`) | **505** | **2** | **Primary Audit Target:** Extreme structural complexity, 34k repeat diagnostics, 2 high-support MOBs. |
| **Cas9_104_105** | `SRR28349176` | MG1655 (`NC_000913.3`) | **140** | **3** | **Intermediate MG1655:** Multi-colony pool, verified P3a SNV target, moderate JC count. |
| **Cas9_1** | `SRR28349244` | MG1655 (`NC_000913.3`) | **49** | **3** | **Baseline MG1655:** Reused P1b baseline, standard clonal profile. |
| **I-SceI_3_4_5_6** | `SRR28349231` | BW25113 (`CP009273.1`) | **72** | **0** | **Max BW25113 Complexity:** Positive endonuclease DSBs, highest BW25113 JC count. |
| **control_1_2_3_4** | `SRR28349226` | BW25113 (`CP009273.1`) | **40** | **0** | **Negative Control:** Negative parental baseline, false-positive control. |

---

### F. Independent Computational Corroboration Protocol (breseq Oracle)

Per AGENTS.md Rule 4, 5, and 12:
1. `breseq 0.40.x` is invoked strictly as an offline computational comparator. It is never called as a runtime engine.
2. Execution: Submitted to Slurm partition `qcpu_18i` with identical inputs:
   ```bash
   breseq -j 8 -r <reference.gbk> -o <outdir>/breseq <run>_1.fastq.gz <run>_2.fastq.gz
   ```
3. Standalone invocation ensures that breseq generates an independent Genome Diff file for matching against ProkaDiff's `output.gd`.

---

### G. Policy Distinguishing Computational Corroboration from Ground Truth

1. **Computational Corroboration $\ne$ Biological Ground Truth:**
   - Concordance between ProkaDiff and breseq demonstrates that both algorithms independently identify the same chimeric sequence configuration in the raw reads.
   - If an Illumina library contains an in-vitro PCR chimera, optical artifact, or template switch, both tools will predictably detect it. Computational agreement is **evidence of algorithmic reproducibility, not wet-lab biological reality**.
2. **Ground Truth Requirements:** Biological ground truth strictly requires orthogonal, non-sequencing assays (e.g., PCR amplification spanning the breakpoint, Sanger sequencing across the junction, or single-molecule long-read sequencing).

---

### H. Orthogonal Evidence Classification Model

The evaluation replaces collapsed 1D tiers with a **4-dimensional orthogonal evidence classification model**. Every structural call is evaluated across all four independent axes:

```
                          [Structural Call]
                                  │
         ┌────────────────────────┼────────────────────────┐
         ▼                        ▼                        ▼
   [Dimension 1]            [Dimension 2]            [Dimension 3]
  Genomic Placement         Read Evidence             Corroboration
  - UNIQUE_ANCHOR           - STRONG                 - CORROBORATED
  - REPEAT_ANCHORED         - LIMITED                - NOT_CORROBORATED
  - REPEAT_AMBIGUOUS        - UNRESOLVED             - NOT_EVALUABLE
  - UNRESOLVED                                             │
                                                           ▼
                                                     [Dimension 4]
                                                     Biological Truth
                                                     - VALIDATED
                                                     - NOT_EVALUABLE
```

#### Dimension 1: Genomic Placement
- `UNIQUE_ANCHOR`: Both breakends map uniquely to non-repeat genomic reference regions ($AS > XS$, MAPQ $\ge 10$, no repeat overlap within 20 bp).
- `REPEAT_ANCHORED`: Exactly one breakend lands in an annotated repeat region ($\le 20$ bp of boundary), and the other lands in unique sequence.
- `REPEAT_AMBIGUOUS`: Read evidence connects multiple potential loci, or both breakends land in multi-copy repeats where the specific physical copy cannot be uniquely identified from short reads ($L_{\text{read}} \le 151$ bp).
- `UNRESOLVED`: Insufficient placement information or extraction failure.

#### Dimension 2: Read Evidence
- `STRONG`: Satisfies pre-registered operational criteria: supporting reads $N \ge 10$, bilateral strand support ($N_{\text{plus}} \ge 2 \land N_{\text{minus}} \ge 2$), and flank extension $\ge 14$ bp on both sides.
- `LIMITED`: Read support present ($N \ge 3$) but fails one or more strong criteria ($3 \le N < 10$, unilateral strand support, or flank overlap $< 14$ bp).
- `UNRESOLVED`: Supporting read metrics missing or below detection threshold ($N < 3$).

#### Dimension 3: Independent Computational Corroboration
- `CORROBORATED`: An independent clean-room caller (`breseq 0.40.x`) calls an equivalent junction matching reference contig, both breakends ($\le 5$ bp microhomology tolerance), and topology under the formal matching contract.
- `NOT_CORROBORATED`: Comparator executed on the same inputs, but no matching junction was called.
- `NOT_EVALUABLE`: Comparator was not run, failed execution, or reference/read inputs were incompatible.

#### Dimension 4: Biological Validation
- `VALIDATED`: Orthogonal wet-lab experimental evidence exists (PCR breakpoint verification, Sanger sequencing, or PacBio/ONT long-read single-molecule validation).
- `NOT_EVALUABLE`: In PRJNA1088182 reference-relative runs without orthogonal assays, all structural variants are biologically `NOT_EVALUABLE` regarding true in-vivo physical presence or Cas9 causality.

*(Optional Composite Confidence Label: For reporting convenience, a call with `UNIQUE_ANCHOR` + `STRONG` + `CORROBORATED` may be summarized as `HIGH_CONFIDENCE_COMPUTATIONAL_CANDIDATE`, while strictly retaining all 4 underlying dimensional labels).*

---

### I. Formal Independent JC Matching Contract (ProkaDiff vs breseq)

A rigorous mathematical matching contract governs whether a ProkaDiff junction $J_A$ and a breseq junction $J_B$ represent the same structural event:

#### 1. Contract Specification
Two junctions $J_A = (c_{1,A}, p_{1,A}, s_{1,A}, c_{2,A}, p_{2,A}, s_{2,A})$ and $J_B = (c_{1,B}, p_{1,B}, s_{1,B}, c_{2,B}, p_{2,B}, s_{2,B})$ match if and only if:
1. **Contig Identity:** Breakends reside on identical reference contigs.
2. **Dual-Breakend Equivalence Under Symmetric Edge Comparison:**
   Junctions are undirected genomic connections. A match exists under either collinear or swapped breakend ordering:
   $$\text{Direct Match:} \quad \Big( c_{1,A} == c_{1,B} \land |p_{1,A} - p_{1,B}| \le 5 \land s_{1,A} == s_{1,B} \Big) \land \Big( c_{2,A} == c_{2,B} \land |p_{2,A} - p_{2,B}| \le 5 \land s_{2,A} == s_{2,B} \Big)$$
   $$\text{Swapped Match:} \quad \Big( c_{1,A} == c_{2,B} \land |p_{1,A} - p_{2,B}| \le 5 \land s_{1,A} == s_{2,B} \Big) \land \Big( c_{2,A} == c_{1,B} \land |p_{2,A} - p_{1,B}| \le 5 \land s_{2,A} == s_{1,B} \Big)$$
3. **Topology Concordance:** The orientation combination ($s_1, s_2$) dictates event topology. A $\pm 5$ bp coordinate match with differing strand signs represents an entirely different structural topology (e.g. deletion vs inversion) and is strictly **NOT A MATCH**.
4. **Deterministic Greedy Bipartite Matching:**
   - Step 1: Exact matches ($D = |p_1 - p_1| + |p_2 - p_2| == 0$, identical strands) are linked first.
   - Step 2: Near matches within microhomology tolerance ($D \le 10$, individual $|p - p| \le 5$) are ranked by minimal combined delta distance $D$.
   - Tie-breaker: Highest supporting read count, then lowest record ID.
   - Once paired, both records are removed from the candidate pool (strict 1:1 matching).

---

### J. Concrete Matching Examples (Synthetic Cases 1 through 7)

```
Case 1: Exact Match              Case 2: Swapped Breakends        Case 3: Orientation Mismatch
ProkaDiff: 1000(+) ~ 2500(-)    ProkaDiff: 1500(+) ~ 3000(-)    ProkaDiff: 1000(+) ~ 2500(-) [DEL]
breseq:    1000(+) ~ 2500(-)    breseq:    3000(-) ~ 1500(+)    breseq:    1000(+) ~ 2500(+) [INV]
Decision:  EXACT_MATCH          Decision:  SWAPPED_MATCH        Decision:  TOPOLOGY_MISMATCH

Case 4: Single Breakpoint Match  Case 5: Competing Candidates     Case 6: Repeat Ambiguity
ProkaDiff: 1000(+) ~ 2500(-)    ProkaDiff: 1000(+) ~ 2000(-)    ProkaDiff: 5000(+) ~ IS1_copyA
breseq:    1000(+) ~ 8500(-)    breseq 1:  1002(+) ~ 2001(-)    breseq:    5000(+) ~ IS1_copyB
Decision:  DISCORDANT_BREAKEND  breseq 2:  1004(+) ~ 2000(-)    Decision:  AMBIGUOUS_FAMILY_MATCH
                                Decision:  Matches B1 (D=3 < 4)
```

1. **Case 1: Exact Structural Match:**
   - $J_A$: `NC_000913.3  1000  1  NC_000913.3  2500  -1`
   - $J_B$: `NC_000913.3  1000  1  NC_000913.3  2500  -1`
   - *Outcome:* `EXACT_MATCH` ($D = 0$).
2. **Case 2: Swapped Breakend Order:**
   - $J_A$: `NC_000913.3  1500  1  NC_000913.3  3000  -1`
   - $J_B$: `NC_000913.3  3000  -1  NC_000913.3  1500  1`
   - *Outcome:* `EQUIVALENT_SWAPPED_MATCH` ($D = 0$, topology identical).
3. **Case 3: Orientation / Strand Mismatch:**
   - $J_A$: `NC_000913.3  1000  1  NC_000913.3  2500  -1` (Deletion topology)
   - $J_B$: `NC_000913.3  1000  1  NC_000913.3  2500  1` (Inversion topology)
   - *Outcome:* `TOPOLOGY_MISMATCH` (No match; strand orientation is non-equivalent).
4. **Case 4: Only One Matching Breakpoint:**
   - $J_A$: `NC_000913.3  1000  1  NC_000913.3  2500  -1`
   - $J_B$: `NC_000913.3  1000  1  NC_000913.3  8500  -1`
   - *Outcome:* `DISCORDANT_BREAKEND` (No match; different structural target).
5. **Case 5: Competing Candidates Within Microhomology Window:**
   - $J_A$: `NC_000913.3  1000  1  NC_000913.3  2000  -1`
   - $J_{B,1}$: `NC_000913.3  1002  1  NC_000913.3  2001  -1` ($D = |1000-1002| + |2000-2001| = 3$)
   - $J_{B,2}$: `NC_000913.3  1004  1  NC_000913.3  2000  -1` ($D = |1000-1004| + |2000-2000| = 4$)
   - *Outcome:* $J_A$ matches $J_{B,1}$ ($D=3 < 4$). $J_{B,2}$ remains unmatched (`NOT_CORROBORATED`).
6. **Case 6: Repeat-Copy Placement Ambiguity:**
   - $J_A$: `NC_000913.3  5000  1  NC_000913.3  250000  -1` (lands in IS1 copy A)
   - $J_B$: `NC_000913.3  5000  1  NC_000913.3  1400000  -1` (lands in IS1 copy B)
   - *Outcome:* `REPEAT_FAMILY_MATCH_AMBIGUOUS_LOCUS` (Dimension 1: `REPEAT_AMBIGUOUS`; Dimension 3: `CORROBORATED_AT_FAMILY_LEVEL_UNRESOLVED_AT_LOCUS`). No 1:1 coordinate match.
7. **Case 7: No Independent Comparator Match:**
   - $J_A$: `NC_000913.3  345000  1  NC_000913.3  345800  -1`
   - Comparator: No junction called within 100 kb.
   - *Outcome:* `NOT_CORROBORATED` (Dimension 3).

---

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

---

### L. Rules for Preserving the Original Frozen Benchmark

1. `data/prjna1088182/` FASTQ and reference files are **STRICTLY READ-ONLY**.
2. P2a baseline results in `manifests/prjna1088182_m6_p2a_run_results.tsv` and `VALIDATION_REGISTRY.tsv` remain **COMPLETELY FROZEN**.
3. All P3c outputs must be written to isolated directories under `benchmark/real_world/results/prjna1088182/m6_p3c/<run>/`.
4. No modification of core engine code (`crates/**`) is permitted.
5. In accordance with M1 contracts, differential `EventId`s must **never** be fabricated for single-sample evidence records (`JC`, `MOB`).
