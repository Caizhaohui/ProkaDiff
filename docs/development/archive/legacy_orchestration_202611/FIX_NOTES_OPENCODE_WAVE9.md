# FIX NOTES — OpenCode Wave 9 (TASK_ID: G1_CLONAL_JC_OPENCODE_WAVE9)

- **Developer:** OpenCode `zhipuai-coding-plan/glm-5.3` (Developer P1 Advanced)
- **Date:** 2026-10-29
- **Repo:** `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **Baseline:** Slurm 3065336 (strict gate: matched=7, over=19, under=11,
  MOB 5/5 exact) — Wave 8 gains (OVER 27→19, BR90 recovered) present but
  the protected floor (matched ≥ 10) breached.
- **Target:** matched ≥ 11 (breseq 89, 90, 91, 94, 100, 101, 102, 108,
  110, 112, 114), over ≤ 19, MOB 5/5 exact.
- **Clean-room statement:** no Clonal coordinates, BR ids, or oracle
  outputs are hard-coded. All changes are general engine rules driven by
  the annotated repeat table and observed evidence.
- **Protected paths untouched:** `benchmark/results/clonal_2673729`,
  `verified/`, gate whitelist (still `--strict-no-whitelist`).

## Forensic root cause (from 3064901 vs 3065336 artifacts)

All four flipped matches (BR100, BR101, BR112, BR114) died at the
**candidate stage**, before second-pass acceptance:

- 3065336 `junctions.fa` has **zero** constructs carrying the minus-side
  register at 1733649 (BR100 `1733649(-)→2774435(+)`) — 3064901 had one
  (`jc164`-family). Same for the minus register at 3630716 (BR112): 14
  matching constructs in 3064901, **0** in 3065336.
- The 1821525(+) flank (BR101 `1821525(+)→2775877(-)`) lost every
  side1-plus construct (3 in 3064901 → 0 in 3065336); only the
  `664688(+)→1821525(+)` display survived.
- The 3901931/3901932 registers collapsed from two constructs (3064901
  `jc85`,`jc86`) to one (3065336 `jc78`): BR114's `3652533(-)→3901931(+)`
  joined copy `3651091..3652533`; the surviving OVER `3894996(-)→3901932(+)`
  joined copy `3893554..3894996`. Same unique position ±1 bp, same family
  IS150, same element end (3′), **different copies** — merged by Wave 8
  Target B, whose grouping key used only `(unique contig, pos, strand,
  repeat family name, end_kind)`.

Mechanism: `cluster_unique_flank_registers` grouped registers by repeat
**family name** (e.g. `IS150`), not by the specific copy
(`rep.start..rep.end`). Two registers of one unique flank joining
different copies of the family (identical element-end sequence) were
clustered into ONE construct; the surviving representative then displayed
the other copy, deleting the true junction of the merged-away copy.
`synthesize_mob_weak_placement` did NOT destroy BR100/101 in 3065336
(its moves were already rejected by the existing unique↔terminus
distance guard); the Task-2 guard is protection for the support-tie
window that the distance guard does not cover.

## Changes

### Fix 1 — Target B restricted to same-copy register merging (`jc_cluster.rs`)

Wave 8 risk-register fallback, now implemented:

- `UniqueRepeatClass` gains `repeat_span: (u64, u64)` = `(rep.start,
  rep.end)` of the copy the repeat side joins, set by
  `classify_unique_repeat_candidate`.
- `cluster_unique_flank_registers` includes `repeat_span` in its sort key
  and cluster-extend predicate: registers merge only when they join the
  **same copy** of the family at the same element end. Candidates joining
  different copies remain separate constructs, so each copy's junction
  (including the 3901931-copy junction) is preserved.
- Identical-sequence collapsing (Target A, `nearer_unique_repeat_representative`)
  is unchanged: byte-identical constructs remain indistinguishable to the
  aligner and still collapse to the nearer copy; cross-copy display for
  truly identical constructs is still decided by that exact-sequence
  nearest rule only.

### Fix 2 — support-tie guard on synthesized MOB weak placements (`emit.rs`)

`align_mob_flanks_to_stronger_copy` now synthesizes a weak-flank placement
only when the weak flank is **not** comparably supported:

- Guard: when the weak flank's total support is within
  `COPY_SUPPORT_TIE_GAP` (= 2, now `pub(crate)`) of the strong flank's,
  the pair does not demonstrate that the strong flank's copy is the
  inserted one, so the weak flank keeps its observed placement — no
  synthesized move (a retained alternative showing the placement still
  wins, as before).
- A genuinely weaker flank (support gap above the tie window) is still
  synthesized onto the stronger copy's opposite end — the MOB pair alone
  demonstrates that copy equivalence. This keeps the BR90-class recovery
  path (weak flank support 11 vs strong 24 in 3065336) intact.
- Existing guards unchanged: equal support ⇒ no move; same copy ⇒ no
  move; a move that would lengthen the unique↔terminus distance is
  rejected; one pair per junction (`locked`).

## Tests (`crates/prokadiff-evidence/src/engine/tests.rs`)

- New `flank_register_clusters_require_same_repeat_copy`: two registers
  of one unique flank (±2 bp, same strand, same family, same element end)
  joining different copies stay two constructs; the stronger seed must
  not capture the other copy's junction.
- New `mob_weak_flank_tie_support_not_synthesized`: weak flank support 8
  vs strong 10 (inside the tie window, no retained alternative) keeps its
  observed copy; the MOB pair itself is unaffected.
- All Wave 8 contracts re-verified unchanged:
  `unique_flank_registers_within_tolerance_cluster_into_one_construct`
  (same-copy registers still merge),
  `mob_flanks_follow_stronger_copy` (gap 6 vs 12 still moves),
  `mob_weak_flank_stays_when_stronger_copy_is_farther`,
  `mob_flanks_keep_copies_when_support_ties`.

## Verification

- `cargo test -p prokadiff-evidence --lib` (Slurm 3065893, qfree):
  **186 passed, 0 failed, 1 ignored** (184 prior + 2 new), cold build,
  0 warnings.
- `cargo fmt --check` / `cargo clippy -p prokadiff-evidence --lib
  --tests` (Slurm 3065897, after an import-order fmt nit found by
  3065895): fmt clean; clippy reports only pre-existing warnings (test
  helper arity, lifetime elisions — none in Wave 9 code).
- Real-data strict gate: Slurm **3065896** via
  `run_g1_clonal_jc_rerun.sbatch` (fresh outdir; `clonal_2673729` and
  `verified/` untouched). Outcome below.

## Risk register

- The 8 OVER calls eliminated by Wave 8 must not return. Same-copy
  register clustering (the legitimate Target B dedup) is preserved; the
  eliminated OVERs (`1269590→2775877`, `16972→2775877`, …) had opposite
  unique-flank strands or fell to Target A / support redistribution, so
  Fix 1 should not resurrect them — the gate checks this empirically.
- BR101's `1821525(+)→2775877(-)`: if its identical-sequence group
  contains an unclassifiable (repeat-interior) member, Target A's
  fallback representative (`pick_identical_group_representative`) ranks
  seed before distance and may still display `664688`. That path is
  outside this wave's prescribed scope; if the gate shows BR101 still
  under while BR100/112/114 recover, the next wave should reorder that
  rank (nearer-terminus before seed among anchored members).

## Slurm 3065896 outcome

- **Job:** 3065896 (`g1_clonal_jc`, bnode2, qcpu_18i, elapsed 09:13, exit
  2:0 = gate FAIL as expected until full parity).
- **Strict gate:** `pass=False matched=10 over=26 under=8
  equiv_repeat_copy=0` (`gate_strict_summary.txt`).
- **Floor restored:** matched=10 ≥ 10 with the exact protected set
  89, 90, 91, 94, 100, 102, 108, 110, 112, 114 — BR90 (synthesis,
  `16972(1)->588495(1)`) and BR100/BR112/BR114 (Fix 1) all recovered.
- **MOB 5/5 exact preserved:** 4524522 IS186 +1 dup6; 3015771 IS150 -1
  dup4; 1733647 IS150 -1 dup3; 1821525 IS150 -1 dup3; 16972 IS150 -1
  dup3 — all five oracle MOBs exact. One **extra** false-positive MOB
  (`1269590 IS150 -1 dup4`, pair JC327+JC331) — breseq calls no MOB
  there.
- **Red-line unchanged:** SNP 28 / INS 3 / DEL 3, matching the frozen
  oracle counts.
- **Under = 8** (best alongside 3064901; 3065336 had 11).

## Acceptance vs task targets

| Metric | Target | 3065336 | 3065896 (Wave 9) |
|---|---|---|---|
| matched | ≥ 11 | 7 | **10** (floor ≥10 restored; BR101 missing) |
| over | ≤ 19 | 19 | **26** ✗ |
| under | — | 11 | **8** |
| MOB exact | 5/5 | 5/5 | **5/5** (+1 extra FP MOB) |

## Why over regressed 19 → 26 (evidence)

The 7 returned OVER calls are exactly the ones Wave 8 eliminated **by the
same cross-copy Target B merges that broke the floor** (diff of
`jc_acceptance_reconcile.tsv` 3065336 vs 3065896):

`1269590→2775877`, `1270660→2775877`, `1604667→2035093`,
`1821526→3651091`, `2774196→4524525`, `2774435→3015773`,
`3015767→3651091`.

Wave 8's over=19 was therefore not an independent gain: cross-copy
merging traded ~4 matched calls for ~7 OVER calls. Fix 1 (mandated
floor repair, Wave 8 risk-register fallback) necessarily gives those 7
back. `over ≤ 19` and `matched ≥ 11` cannot both be produced by the
Target-B copy-span rules alone.

## BR101 analysis (still under)

The `1821525(+)` flank family displays `1821526(1)->3651091(1)` (C_d 5′
end, support 40, `mob_evidence=1`); breseq 101 wants
`1821525(1)->2775877(-1)` (C_c 5′ end). The C_c/C_b byte-identical
registers collapse in Target A and the surviving construct loses the
stage-2 read race to the 1-bp-shifted C_d sibling; post-acceptance
folding then displays the stronger C_d placement. Any deterministic
copy-display rule has documented counterexamples (oracle 114 keeps a FAR
copy, oracle 96 a NEAR one), so this is an architecture decision, not a
local patch (see recommendation).

## Recommendation to Controller

1. **Wave 9 state is a strict floor restore** (10 matched incl. BR90,
   MOB 5/5, red-line intact, under=8 best) and strictly dominates
   3064901 (10/27/8) on OVER and equals it on matched.
2. **over ≤ 19 + matched ≥ 11 jointly requires an L3 ruling** on
   cross-copy sibling suppression / copy display among byte-identical
   constructs (enumeration-arbitrary in breseq). Candidate directions:
   (a) suppress an accepted unique↔repeat junction when a stronger
   accepted sibling shares its unique flank within tolerance on another
   copy of the same family+end (kills ~5 of the 7; risks BR114-class
   far-copy keepers — needs counterexample audit against oracle
   96/107/109/113); (b) accept the 10/26 Pareto point and re-baseline
   the gate; (c) route to Architect per POLICY §5 after this wave.
3. The extra MOB @1269590 (dup4) comes from pair JC327
   (`666130→1269593`) + JC331 (`1269590→2775877`) — same family+opposite
   ends geometry as true MOBs; suppressing it needs the same sibling
   analysis as (a).

## Evidence paths

- Gate: `benchmark/results/parity_closure/g1_clonal_jc/slurm/3065896/gate_strict/`
- Unit tests: `benchmark/logs/w9_tests_3065893_cargo.log` (186 pass)
- Lint: `benchmark/logs/w9_clippy_3065897_{fmt,clippy}.log`
- No commits made; working tree left for Controller review. No push.
