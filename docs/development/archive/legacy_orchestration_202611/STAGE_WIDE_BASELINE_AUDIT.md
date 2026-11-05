# Stage-wide baseline audit

**Date:** 2026-11-03  
**Controller:** Cursor CLI / Grok 4.7  
**HEAD:** `1a0c7b553753bfa1fd6f0c0741c0a8bd43d4cf3c` on `main` (dirty)  
**Allocation:** Slurm `3095377` on `bnode2`

Authoritative stage boundary is [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md). M0–M6 have in-repo closeouts. M7 is closed with limitations on the dedicated worktree, not by a registry promotion. M8 is the active milestone. M9 has not started.

## Verified dispositions

| Stage | Disposition | Evidence |
| --- | --- | --- |
| M0–M4 | COMPLETE | Plan entries and recorded `qcpu_18i` jobs |
| M5 | COMPLETE | [benchmark/real_world/M5_CLOSEOUT.md](benchmark/real_world/M5_CLOSEOUT.md) |
| M6 | CLOSED_WITH_LIMITATIONS | [benchmark/real_world/M6_CLOSEOUT.md](benchmark/real_world/M6_CLOSEOUT.md). Registry `cohort_25` stays `PENDING` |
| M7 | COMPLETE_WITH_LIMITATIONS | Closeout only at `../04_ProkaDiff-m7-worktree/benchmark/real_world/M7_CLOSEOUT.md`. `sacct`: 3049255 `0:0` 6 min; 3049330 `1:0` 1:40:13; 3051116 `0:0` 1:06:44. Registry `prjna884016` stays `PENDING`. Positive chromosomal sensitivity is `NOT_EVALUABLE` |
| M8 / G1 | BLOCKED_AUTHORIZATION | Strict gate job 3066886: matched 10, OVER 24, UNDER 8, `pass=false`. `docs/parity.md` forbids a distance rule. No production attribution edit. No new G1 job |
| M9 | NOT_STARTED | Depends on M8 |

## This-turn check

`python3 -m pytest benchmark/real_world/scripts/test_m7_p0b_scoring.py -q` on job 3095377: **19 passed**, exit 0. The fail-closed M7 scorer on this worktree still emits `NOT_EVALUABLE` and does not promote the registry.

A second pass on the same allocation, `test_prjna884016_quarantine.py`, `test_compare_variant_calls.py`, and `test_m6_reference_policy.py`: **19 passed**, exit 0. These lock the quarantine, the `NOT_EVALUABLE` sentinel, and the M6 reference policy. They do not open M8 production work.

The remaining real-world script tests on the same allocation: **94 passed, 4 failed**, exit 1. The failures are `test_contracts_and_crates_unchanged`, `test_contracts_unchanged`, `test_versioned_manifest_matches_builder_and_contracts_are_unchanged`, and `test_gd_ids_are_not_written_as_event_ids_and_registry_is_untouched`. Each asserts that `git diff --stat -- crates benchmark/real_world/VALIDATION_REGISTRY.tsv` is empty. `VALIDATION_REGISTRY.tsv` itself is clean. The diff is nine files, 5389 insertions and 143 deletions: `bam_io.rs`, `emit.rs`, `jc_cluster.rs`, `tests.rs`, `jc_seq.rs`, `pileup.rs`, `repeat_ambiguous.rs`, `split_seed.rs`, and `crates/prokadiff/tests/m5_sr1_evidence_only.rs`. The freeze checks are doing what they were written to do. They were not weakened, and the engine diff was not reverted.

The uncommitted `docs/parity.md` adds both instructions in one diff. Items 7 and 10 tell the fold to keep the nearer copy when support ties, and item 10 puts distance ahead of seed. The Wave 8 paragraph in the same diff says not to add a general nearest-copy display rule, because oracle 114 keeps the far copy and oracle 96 keeps the near copy. The file was not edited to pick one of those sentences. That choice is an authorization change, not an implementation fix.

## G1 resource record (not a speed claim)

`sacct` for the existing jobs, exit `2:0` is the strict gate's non-zero exit:

| Job | Elapsed | MaxRSS | CPUs | Partition |
| --- | --- | --- | --- | --- |
| 3066886 | 00:09:24 | 3,570,756 K | 8 | qcpu_18i |
| 3072550 | 00:10:17 | 3,579,356 K | 8 | qcpu_18i |

Breseq wall time is not a separate Slurm step in these jobs. No parity or speed claim is made.

## Not eligible now

- G1 production edit or a new G1 Slurm acceptance job.
- Promoting `VALIDATION_REGISTRY.tsv` for PRJNA1088182 or PRJNA884016.
- Restoring the quarantined 4 MOB + 7 JC truth.
- A second supervisor.
