# FIX NOTES — OpenCode Wave 8 (TASK_ID: G1_CLONAL_JC_OPENCODE_WAVE8)

- **Developer:** OpenCode `zhipuai-coding-plan/glm-5.3` (Developer P1 Advanced)
- **Date:** 2026-10-29
- **Repo:** `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **Baseline:** Slurm 3064901 (`run_g1_clonal_jc_rerun.sbatch`, strict gate:
  matched=10, over=27, under=8, MOB 5/5 exact)
- **Protected floor:** matched JCs ≥ 10 (breseq 89, 91, 94, 100, 101, 102,
  108, 110, 112, 114), MOB 5/5 exact (4524522 IS186 +1 dup6; 3015771 IS150
  -1 dup4; 1733647 IS150 -1 dup3; 1821525 IS150 -1 dup3; 16972 IS150 -1
  dup3), red-line SNP/INS/DEL unchanged.
- **Clean-room statement:** no Clonal coordinates, BR ids, or oracle outputs
  are hard-coded. All changes are general engine rules driven by the
  annotated repeat table and observed evidence.

## Changes

### Target A — deduplicate identical construct sequences (`jc_cluster.rs`)

`collapse_identical_unique_repeat_copies` now collapses EVERY
identical-sequence candidate group to one representative, not only groups
that pass the strict unique↔repeat classification:

- Groups where some member's repeat side landed in a copy **interior**
  (`repeat_containing` ⇒ classification `None`) previously fell through the
  all-or-nothing `nearer_unique_repeat_representative` and kept every
  member. On 3064901 this left >20 byte-identical constructs per family in
  `junctions.fa` (1102 constructs total, 57 duplicated sequences, up to 26
  copies of one 72-bp both-repeat sequence), diluting stage-2 bowtie2 reads
  below `accept_junction` (the BR113 failure mode).
- New guard `identical_group_has_one_flank`: a group collapses when its
  classifiable members share one unique flank within `JC_CLUSTER_TOL_BP`
  (unclassifiable members must not introduce a different flank), or when
  every member is repeat-anchored (repeat↔repeat groups). Groups pairing
  two distinct far-apart unique flanks remain intact (distinct junctions).
- Representative (`pick_identical_group_representative`): prefer a member
  anchored on an annotated terminus, then strongest `seed_reads`, then
  nearer unique↔terminus distance, then deterministic coordinates. The
  representative carries the group's **max** seed.

### Target B — cluster 1–5 bp-shifted unique-flank registers (`jc_cluster.rs`)

New pass `cluster_unique_flank_registers` runs after the identical-sequence
collapse. Unique↔repeat representatives whose unique flanks share contig +
strand and sit within `JC_CLUSTER_TOL_BP`, joining the same repeat family
at the same element end (5' or 3'), are one junction whose register was
enumerated at several offsets (e.g. `1821527` vs `1821528`). Each cluster
keeps one representative (strongest seed, then nearer copy, then
coordinates) with the cluster's max seed. Opposite element ends (a MOB's
two flanks), opposite strands, and flanks farther than 5 bp never cluster,
so an insertion's plus/minus flank pair cannot be merged.

### Target C — MOB weak-flank copy disambiguation (`emit.rs`)

`align_mob_flanks_to_stronger_copy` previously moved the weaker flank onto
the stronger copy only when a retained alternative showed that placement;
otherwise the weak flank stayed on its (arbitrary) copy. The MOB pair
itself demonstrates the copy equivalence (same family, opposite ends, one
insertion event), so when no alternative exists the placement is now
**synthesized** (`synthesize_mob_weak_placement`): only the weak flank's
repeat side moves onto `place_element_end(strong.rep, weak_end_kind)`, the
unique side and frozen scores are untouched. Existing guards retained:
equal support ⇒ no move; same copy ⇒ no move; a move that would lengthen
the unique↔terminus distance is rejected; one pair per junction (`locked`).

**Deliberately NOT implemented:** a blanket "nearest annotated copy"
display rule for unique↔repeat junctions. Analysis of oracle 114
(`3652533(-1)→3901931(1)`, a protected match) shows breseq keeps a FAR copy
there (nearest 3'-end copy would be 3894996, ~7 kb away), while oracle 96
uses the NEAR copy — breseq's copy choice among byte-identical copies is
enumeration-arbitrary and not reproducible by any distance rule. A blanket
nearest rule would break the floor (114 and 109 both regress). BR96's
`664688` selection is instead expected from Targets A+B: the 1270660-flank
identical-construct group contains observed members on multiple copies, and
the nearest-observed representative rule selects `664688`.

## Tests (`crates/prokadiff-evidence/src/engine/tests.rs`)

- `identical_repeat_repeat_constructs_collapse_to_one_representative`
  (replaces wave-7 `..._are_not_collapsed`, whose contract Target A
  reverses): repeat↔repeat identical group ⇒ 1 representative, max seed.
- `interior_repeat_side_member_does_not_block_identical_group_dedup`:
  interior-side sibling no longer keeps the whole group alive.
- `unique_flank_registers_within_tolerance_cluster_into_one_construct`:
  1–3 bp-shifted registers ⇒ one construct, strongest-seed register, max
  seed (Target B).
- `flank_register_clusters_require_same_element_end`: 5'/3' pair (a MOB's
  two flanks) never merges.
- `flank_register_clusters_require_same_strand_and_distance`: opposite
  strands and >5 bp flanks stay separate.
- `mob_flanks_follow_stronger_copy` (updated to wave-8 contract): the weak
  flank now adopts the stronger copy's opposite end without a retained
  alternative; MOB position/dup unchanged; unique side intact.
- `mob_weak_flank_stays_when_stronger_copy_is_farther`: synthesized moves
  that would lengthen the unique↔terminus distance are rejected.

Housekeeping: `crates/prokadiff/tests/m5_sr1_evidence_only.rs` was
pre-broken in the working tree by wave 7's uncommitted
`RepeatAmbiguousSeed::failure_stage` field (missing initializer field);
added `failure_stage: None`.

## Verification

- `cargo test -p prokadiff-evidence --lib`: **184 passed, 0 failed**
  (inside Slurm allocation 3048499, bnode1).
- `cargo test --workspace`: **all crates pass**.
- `cargo clippy -p prokadiff-evidence --lib --tests`: no new warnings from
  this wave (remaining warnings are pre-existing wave-7 code).
- `cargo fmt -p prokadiff-evidence`: applied.
- Real-data strict gate: Slurm **3065336** submitted via
  `run_g1_clonal_jc_rerun.sbatch` (fresh outdir `slurm/3065336/`;
  `clonal_2673729` and `verified/` untouched). Result: see
  "Slurm 3065336 outcome" below once complete.

## Risk register (floor-critical, watch in 3065336)

- Target B merges register pairs like {JC339-family, JC323-family}
  (unique flanks 1821527/1821528) and {JC347-family, JC348-family}
  (3901931/3901932). The surviving representative is chosen by strongest
  seed; if the seed ordering is adverse, a protected match (102 or 114)
  could flip. The strict-gate run is the acceptance check; if the floor
  breaks, fallback is to restrict Target B to same-copy register merging
  (cross-copy display then decided by the exact-sequence nearest rule only).
- Target A changes which constructs enter `junctions.fa` (~1102 → fewer),
  so support distributions shift; MOB 5/5 and red-line SNP/INS/DEL must be
  re-checked on the gate output.

## Slurm 3065336 outcome

_Pending — fill in after job completion._
