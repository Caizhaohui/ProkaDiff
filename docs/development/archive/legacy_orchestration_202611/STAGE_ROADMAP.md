# Stage roadmap

**Active milestone:** M8 REL606 engine parity. M7 is closed with limitations. M9 is not open.

| Workstream | Deliverables | Tests | Review | Status |
| --- | --- | --- | --- | --- |
| M0–M6 | In-repo closeouts | Recorded cluster jobs | Historical | COMPLETE or CLOSED_WITH_LIMITATIONS |
| M7 negative control | Closeout on the m7 worktree; registry stays PENDING | `test_m7_p0b_scoring.py` 19 passed, 2026-11-03, job 3095377 | Closeout already reviewed on that worktree | COMPLETE_WITH_LIMITATIONS |
| M8 G1 strict JC gate | Job 3066886 gate | matched 10, OVER 24, UNDER 8, pass false | Claude: no production change. Codex D2: wording PASS | BLOCKED_AUTHORIZATION |
| M8 same-job breseq wall/RSS | Not separately accounted | `sacct` has only the combined job | None | BLOCKED_EXTERNAL |
| M9 report review | Not started | — | Human review required by the plan | NOT_STARTED |

## Task matrix

| Task | Workstream | Status | Why |
| --- | --- | --- | --- |
| `m7_p0b_regression` | M7 | COMPLETE | 19 passed on 3095377 |
| `m8_g1_production` | M8 | BLOCKED_AUTHORIZATION | Distance rule forbidden; Claude `PRODUCTION_CHANGE: NO` |
| `m8_g1_new_slurm` | M8 | BLOCKED_AUTHORIZATION | No new G1 acceptance job |
| `m8_breseq_split_timing` | M8 | BLOCKED_EXTERNAL | Existing jobs do not split breseq wall time |
| `m6_crate_freeze_tests` | M6 | FAILED_RETRYABLE | Four tests fail because `crates/prokadiff-evidence` is dirty. Do not weaken the tests. Do not revert the engine diff in this session. |
| `m9_report_review` | M9 | BLOCKED_AUTHORIZATION | Plan places it after M8 |

G1 investigation files under `benchmark/results/parity_closure/g1_clonal_jc/` stay read-only evidence. They do not change the strict gate.
