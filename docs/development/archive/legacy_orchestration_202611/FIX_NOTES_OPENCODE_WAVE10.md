# FIX NOTES — OpenCode Wave 10 (G1_CLONAL_JC_OPENCODE_WAVE10)

- **Date:** 2026-10-29
- **Developer:** OpenCode `zhipuai-coding-plan/glm-5.3` (Developer P1 Advanced)
- **Scope:** Codex Sol [P2] interleaved-copy clustering defect; duplicate-register elimination; deterministic unit test; FP MOB @ 1269590 audit.
- **Repo state:** working tree of `04_ProkaDiff` @ `1a0c7b5` (dirty — carries uncommitted Waves 8/9 G1 work; Wave 10 delta described below). Single writer respected; no git push; no commits made.

## Target A — interleaved-copy clustering defect (jc_cluster.rs)

**Defect (Codex Sol [P2], confirmed):** `cluster_unique_flank_registers` sorted
`keyed` by `(contig, unique_pos, strand, family, copy_span, end)`. Because the
flank POSITION sorted before the repeat-copy identity, a register joining copy B
could sit between two ≤5 bp registers of copy A. The subsequent linear scan
compares only the cluster's LAST member and flushes on any key mismatch, so the
copy-B interloper split copy A's registers back into one construct per offset —
exactly the duplicate-register defect Wave 8 Target B was meant to remove.

**Fix:** sort key reordered to
`(contig, strand, repeat_family, repeat_span, end_kind, unique_pos)` — copy
identity first, position last. All registers joining one copy are now
contiguous; inside a run they are position-sorted, so every register within
`JC_CLUSTER_TOL_BP` (5 bp) of the run merges. Cross-copy merging remains
impossible: `repeat_span` is still a required-equal extend key (Wave 9
invariant preserved — see test `flank_register_clusters_require_same_repeat_copy`).
Function doc comment updated to record the contiguity requirement.

**Wave 10 delta:** `crates/prokadiff-evidence/src/engine/jc_cluster.rs`
(sort comparator + comments only; no behavioral change anywhere else).
Downstream ordering is unchanged: the function re-sorts its output by
`(side1_contig, side1_pos_1, side2_contig, side2_pos_1)` as before.

## Target B — deterministic unit test (tests.rs)

New test `interleaved_repeat_copy_registers_merge_same_copy_cleanly`:
3 candidates — unique 100→copy A (1000..1100, seed 5), unique 101→copy B
(3000..3100, seed 6), unique 102→copy A (seed 7). Asserts: exactly 2 constructs;
copy A keeps ONE representative at register 102 carrying max seed 7; copy B
keeps ONE representative at 101.

**Red-proof (job 3066797, qfree):** with the pre-fix sort order temporarily
restored, the new test FAILS (`test result: FAILED. 0 passed; 1 failed`,
rc=101). With the fix (job 3066790): passes. The test discriminates the defect
deterministically.

## Target C — FP MOB @ 1269590 audit (Slurm 3065896 evidence)

**Emitted (3065896):** `MOB 352 327,331 REL606.6 1269590 IS150 -1 4` from
JC327 `666130(-1)→1269593(-1)` (support 25) + JC331 `1269590(1)→2775877(-1)`
(support 28). Breseq oracle: NO MOB and NO junction at 1269590.

**Root cause (measured, not hypothesized):** REL606 carries a paralogous block
around 1269590 that is 115/120 identical over a ±60 bp window and 236/250 over
a ±(200+50) bp window to the block around 1270660. The Clonal sample's true
event is an IS150-boundary junction at 1270660 — breseq JC96
`664688(1)→1270660(1)` (note: breseq calls NO MOB there either; it is a
side_2=repeat junction, not an insertion). ProkaDiff's stage-2 alignment maps a
share of the same split reads onto the paralogous displaced locus, producing
three over-called junctions: JC322 `590471→1269591` (25), JC327
`666130→1269593` (25), JC331 `1269590→2775877` (28) — all strict-gate
ACCEPTED_OVER in the 3065896 reconcile TSV. ProkaDiff ALSO calls the correct
paralog placement (JC332 `1270660→2775877`, support 35 ≡ breseq JC96 under
repeat-copy equivalence).

**Why the MOB promoter fires:** JC327's repeat side anchors the 3' end of the
+ copy 664688..666130; JC331's repeat side anchors the 5' end of the − copy
2774435..2775877. Targets 1269593(−)/1269590(+) are opposite strands 3 bp
apart at a unique-mapping site with adequate unique coverage. Every gate in
`emit.rs` MOB promotion legitimately passes (opposite element ends, same family
IS150, dup=4 after no microhomology rewrite). The pair is geometrically
indistinguishable from a real insertion *given the two input junctions*.

**Why constituent JC filtering cannot remove it:** MOB constituents are exempt
from secondary support/frequency filtering (`mob_constituent_jcs`, emit.rs), and
even unexempt JC327/JC331 far exceed the thresholds (support 25/28 ≥ 3; local
frequency ≫ 0.05×depth). The reads are real — they are displaced onto the
paralog, not noise.

**Does Target A eliminate it?** No — structurally impossible. JC327 (minus
flank @1269593 joining copy 664688..666130) and JC331 (plus flank @1269590
joining copy 2774435..2775877) are different (strand, copy) clusters; Target A
only merges same-copy ≤5 bp registers and can only strengthen survivors. (Rerun
confirmation below.)

**Clean-room remedy direction (proposal for Controller / Wave 11 — NOT
implemented here):** general paralog-flank suppression among accepted
unique↔repeat junctions sharing repeat family+element end — when two accepted
junctions' unique-flank windows are near-identical (e.g. ≥90% identity over the
construct flank), keep the stronger-supported placement and suppress the
paralog-displaced sibling(s). This is a general rule (no Clonal coordinates, no
strain whitelist); it would suppress JC327/JC331/JC322 (~25–28 support) in
favor of JC332 @1270660 (35). Risks a Wave-11 counterexample audit against
oracle far/near copy keepers (BR96/107/109/113/114 class) before adoption.

## Verification evidence

| Item | Result | Evidence |
|---|---|---|
| Unit tests | 187 passed, 0 failed, 1 ignored | job 3066790 (qfree, vnode2), `benchmark/logs/w10_tests_3066790_cargo.log` |
| Red-proof (pre-fix sort) | new test FAILS rc=101 → discriminating | job 3066797, `benchmark/logs/w10_redproof_3066797_cargo.log` |
| Post-restore unit tests | 187 passed | job 3066810, `benchmark/logs/w10_verify_3066810_cargo.log` |
| `cargo fmt --check -p prokadiff-evidence` | clean (empty log) | job 3066810, `benchmark/logs/w10_verify_3066810_fmt.log` |
| Clippy (project convention, `--lib`) | rc=0, 3 lib warnings — identical set to Wave 9 (jc_cluster.rs:161 too_many_arguments 11/7; if_same_then_else in `classify_unique_repeat_candidate`; split_seed.rs:659 needless lifetimes; + 8/7 in tests). None introduced by Wave 10 (delta touches no signatures/branches) | job 3066820, `benchmark/logs/w10_clippy_3066820_clippy.log` vs `benchmark/logs/w9_clippy_3065897_clippy.log` |

## G1 rerun (Slurm 3066886, qcpu_18i)

Slurm 3066886 · qcpu_18i · bnode5 · 8 threads · 9m24s · MaxRSS 3.57 GB ·
git `1a0c7b5` dirty=44 · breseq 0.40.2 · bowtie2 2.5.4 (run_header.txt).
Command: `sbatch benchmark/real_world/scripts/run_g1_clonal_jc_rerun.sbatch`.
Job exit 2:0 AFTER the strict gate completed — identical post-gate
script/layer2 behavior as Wave 9's 3065896 (also FAILED 2:0, no layer2
files); not a Wave 10 regression. All application outputs and gate runs are
complete.

| Metric (strict, no whitelist) | Wave 9 (3065896) | Wave 10 (3066886) | Δ |
|---|---|---|---|
| matched (cost=0) | 10 | **10** (89,90,91,94,100,102,108,110,112,114 all cost=0) | floor held |
| accepted_over | 26 | **24** | −2 |
| accepted_under | 8 | **8** | = |
| breseq_rejected | 9 | **9** | = |
| MOB exact (5 true oracle) | 5/5 | **5/5** exact (position, element, strand, dup) | held |
| extra FP MOB @1269590 | 1 | 1 (MOB 350 = JC326 666130→1269593 + JC330 1269590→2775877) | persists (Target C) |
| JC records | 36 | **34** | −2 duplicates |

| Metric (standard gate, repeat-copy equivalence) | Wave 9 | Wave 10 | Δ |
|---|---|---|---|
| matched exact | 10 | **10** | = |
| equiv_repeat_copy | 6 | **6** | = |
| total equivalent | 16/18 | **16/18** | = (BR95, BR113 remain under) |
| accepted_over | 20 | **18** | −2 |
| accepted_under | 2 | **2** (BR95, BR113) | = |

Standard-gate evidence recomputed post-hoc with the same command shape as
Wave 9's controller run (--allow-repeat-copy-equivalence
--normalize-contig-version) and persisted at
benchmark/results/parity_closure/g1_clonal_jc/slurm/3066886/gate_standard/.

Duplicate-register eliminations (set-diff of emitted JCs, verified):
exactly two constructs removed, zero added, zero lost:

1. `REL606:588495(1)->1733649(1)` (Wave 9 JC321, support 46): a
   duplicate-register survivor of the plus-flank family at 1733648/49. The
   interleaved-copy split had left both a copy-C-anchored (588495) and a
   copy-A-anchored (664688) construct of the same flank; post-fix the
   family consolidates and the surviving plus flank is
   `JC324 664688(1)->1733648(1)` (support 43), paired with
   `JC327 666130(-1)->1733648(-1)` (31) in MOB 348 @1733647 — same MOB
   position/element/strand/dup as the oracle.
2. `REL606:2774196(-1)->4524525(-1)` (Wave 9 JC339, support 43):
   IS186-region duplicate register eliminated.

Red line: SNP=28, DEL=3, INS=3 — field-level identical diff vs 3065896
(REDLINE_IDENTICAL).

## Protected floor check (strict gate, no whitelist)

- Matched JCs >= 10: **PASS** — 89, 90, 91, 94, 100, 102, 108, 110, 112, 114
  all MATCHED at cost=0 (reconcile rows PD 312/317/318/320/335/337/339/342/343/344).
- 5 true oracle MOBs exact: **PASS** — `4524522 IS186 1 6`, `3015771 IS150 -1 4`,
  `1733647 IS150 -1 3`, `1821525 IS150 -1 3`, `16972 IS150 -1 3`.
- Strict gate overall verdict remains `pass=False` because
  accepted_over=24 > 19 (same failing dimension as Wave 9's 26; improved).
- **MOB_PARITY still NOT passed**: `MOB_COVERAGE_5_OF_5 != MOB_PARITY_PASS`;
  the FP MOB @ 1269590 (dup 4) persists — Wave 10 Target A was not expected
  to and did not remove it (Target C analysis above; rerun confirms
  JC326/JC330 supports 25/28 and pair geometry unchanged).

## Bottom line for Controller

1. **Wave 10 strictly dominates Wave 9**: identical protected floor (10
   matched cost=0, 5/5 MOB exact, 16/18 equivalent, red line untouched),
   with 2 fewer strict OVER and 2 fewer standard OVER, via removal of two
   genuine duplicate-register artifacts. No new calls, no losses.
2. The FP MOB @1269590 is a **paralog-displacement over-call** (measured
   94–96% identity block vs the true 1270660 locus; breseq JC96 is a
   repeat-side junction there, NOT a MOB). It cannot be removed by register
   clustering (different strand/copy clusters) nor by constituent-JC
   filters (exempt by design; supports 25/28 far above thresholds). The
   clean-room remedy candidate is the general paralog-flank suppression
   rule proposed in Target C — needs a Wave-11 counterexample audit
   (BR96/107/109/113/114 far/near copy keepers) or an L3 ruling.
3. Remaining known gaps unchanged: BR95/BR113 under (architecture-level
   copy-display decision, per Wave 9 notes), strict over=24.

## Constraints honored

- No hardcoded Clonal coordinates; no strain-specific whitelist expansions.
- All 5 true oracle MOBs and protected matched JCs re-checked in rerun (below).
- Single writer: no other Developer process touched the tree; no `git push`.
- Compute: all cargo/evidence work ran inside Slurm allocations (qfree for
  unit/lint jobs; qcpu_18i for the WGS rerun), ≤2 running batch jobs at a time
  under this Developer's management.
