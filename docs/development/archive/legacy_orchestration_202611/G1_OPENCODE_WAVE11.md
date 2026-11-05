# Wave 11 Task Specification: Paralog Disambiguation, Reverse-Complement Canonicalization, and Repeat-Repeat Junction Recovery

**TASK_ID:** `G1_OPENCODE_WAVE11`  
**Developer Role:** Developer P1 Advanced (`zhipuai-coding-plan/glm-5.3`)  
**Controller:** Antigravity CLI (`gemini-3.8-flash-medium`)  
**Target Worktree:** `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff` (main repo)  

---

## Mandate & Policy

Read and follow `~/.config/ai-collaboration/POLICY.md` before acting. It is not auto-loaded. If it cannot be read, stop and report `POLICY_UNREADABLE`.

---

## Scientific Root Cause Analysis & Empirical Evidence

Detailed audit of the Clonal sample alignment artifacts (Slurm Job 3066886, `aligned.sam`, and `junctions.fa`) established the following biological and algorithmic facts:

1. **Paralogous Cross-Mapping & FP MOB @ 1269590:**
   - In Clonal sample, IS150 inserted at 1270660 in `ldrC` (breseq JC 96: `664688(1) -> 1270660(1)`, ProkaDiff PD 331: `1270660(1) -> 2775877(-1)`).
   - In REL606, `ldrA/B` (around 1269590) and `ldrC` (around 1270660) are paralogous tandem repeats differing by only 2 bp over the junction window.
   - 53 reads in `aligned.sam` mapped to `jc194` (1269590). **Every single one of these 53 reads also mapped to `jc112/jc113` (1270660)**.
   - For 49 of these reads, the alignment score was 1 to 10 points **higher on 1270660 than on 1269590** (`AS:i:36` vs `AS:i:26..32`). The remaining 4 reads did not even span the breakpoint (starting after pos 36).
   - `parse_junction_bam` in `bam_io.rs` currently iterates over alignments independently without comparing scores across candidate junctions. It accepted all 53 sub-optimal alignments to `jc194`.
   - Crediting sub-optimal cross-mapping reads created JC 330 (`1269590(1) -> 2775877(-1)`).
   - JC 330 paired with JC 326 (`666130(-1) -> 1269593(-1)`) to form FP MOB 350 (`1269590 / IS150 / -1 / dup 4`), which bypassed secondary frequency filtering.

2. **Reverse-Complement Duplicate Construct Escaping Collapse:**
   - `3630716(-1) <-> 3652533(-1)` (BR 112, an `ldrD <-> IS150` junction) and `666130(-1) <-> 1269593(-1)` (JC 326, an `IS150 <-> ldrA/B` junction) share identical 72 bp construct sequences in reverse complement: `cand1 == rc_dna(&cand2)`.
   - `collapse_identical_unique_repeat_copies` in `jc_cluster.rs:1254` keyed `groups` by raw sequence `seq`. Because `seq != rc_dna(&seq)`, reverse-complement constructs were not collapsed, emitting duplicate opposite-strand constructs into `junctions.fa` and splitting reads.

3. **Repeat-Repeat Junctions Silently Filtered (BR 95 and BR 113):**
   - BR 113 (`3651091(1) -> 3741531(-1)`) was present in `junctions.fa` (`jc182`, `jc437`) with 39-42 balanced reads.
   - It was silently dropped in `emit.rs:1035` because `interior(r2, 3741531)` evaluated true for IS1, misclassifying a genuine 42-read structural rearrangement as an IS-end chimera.
   - BR 95 (`589886(1) -> 590471(-1)`) was present in `junctions.fa` (`jc174..jc178`) with 15-16 reads. It was dropped in `emit.rs:1064` because `local_depth * 0.05` was computed against the inflated multi-copy sum of IS150 + IS1 repeats (~400x). In breseq, repeat-repeat junctions have `frequency=NA` and are exempt from single-copy frequency filtering.

---

## Required Implementations

### 1. `crates/prokadiff-evidence/src/engine/bam_io.rs`
In `parse_junction_bam`:
- Add competing-junction alignment score resolution:
  - When decoding junction BAM records, evaluate qualified spanning alignments (`q_cover >= JC_MIN_COVER_BASES` and `spans_breakpoint(cstart, cend, breakpoint, JC_MIN_SPAN_EACH_SIDE)` and `jc_score >= rs`).
  - Compute candidate alignment score. Prefer `AS:i:` tag if present (falling back to `cigar_match_indel_score(m, i, d)` or subtracting `NM:i:` if AS is absent, preserving unit test compatibility).
  - Track for each molecule `(norm_qname, mate)` the maximum score achieved across candidate junctions: `best_scores: HashMap<(String, Mate), i32>`.
  - Filter out alignments whose score is strictly less than `best_scores[&(norm_qname, mate)]` with counter `stats.worse_than_competing += 1` (or record in `SecondPassReject`).
  - **Invariance:** Tied scores across identical repeat copies (`jc_score == best_score`) MUST BE PRESERVED to maintain multi-copy evidence and unresolved ambiguity (Codex Sol Recommendation 2).

### 2. `crates/prokadiff-evidence/src/engine/jc_cluster.rs`
In `collapse_identical_unique_repeat_copies`:
- Canonicalize the grouping key:
  ```rust
  let rc_seq = crate::jc_seq::rc_dna(&seq);
  let key_seq = if seq <= rc_seq { seq } else { rc_seq };
  groups.entry(key_seq).or_default().push(cand);
  ```
  Ensure helper functions in `jc_cluster.rs` handle reverse-complement members when selecting representatives.

### 3. `crates/prokadiff-evidence/src/engine/emit.rs`
- In secondary filtering (`emit.rs:1064`):
  - When `pd_repeat_junction` is true (both sides are in annotated repeats):
    Exempt from `opts.jc_min_frequency` (because repeat depth is multi-copy; matching breseq's `frequency=NA` convention). Require only `support_reads >= opts.jc_min_support_reads`.
- In repeat-repeat filtering (`emit.rs:1035`):
  - Do not drop repeat-repeat junctions if `j.support.plus_reads >= 5 && j.support.minus_reads >= 5` and `accept_junction(&j.support)` is satisfied.

### 4. Deterministic Unit Tests (`crates/prokadiff-evidence/src/engine/tests.rs`)
- Add unit test: `second_pass_rejects_suboptimal_competing_junction_alignment` verifying that a read aligning to construct A with score 36 and construct B with score 26 is only credited to construct A.
- Add unit test verifying that tied scores on identical repeat copies are both retained.
- Add unit test verifying reverse-complement canonical grouping in `collapse_identical_unique_repeat_copies`.
- Verify all existing 187+ unit tests continue to pass.

---

## Acceptance Criteria
1. `cargo test -p prokadiff-evidence --lib` exits 0 (all 187+ tests pass + new regression tests).
2. Clonal rerun on Slurm compute node (`qcpu_18i`):
   - Strict gate: matched >= 10 at cost=0, extra FP MOB @ 1269590 eliminated (0 FP MOBs), 5/5 exact oracle MOBs preserved.
   - Standard gate: matched 18/18 (under=0), OVER significantly reduced.
3. No hardcoding of coordinates or strain names. Clean-room Path B only.
