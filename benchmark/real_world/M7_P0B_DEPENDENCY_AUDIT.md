# M7-P0b Dependency Audit

Milestone: Ground-Truth Dependency Audit and Fail-Closed Hardening.

Final status: **M7_P0B_PUSH_READY**

This ledger records completed CLI executions. Planned or failed attempts are marked separately. No GitHub push, tag, or release was created.

## Git state

| Ref | Hash |
| --- | --- |
| `origin/main` | `5c90e0568153a4befc573933d2f98bd3d47662b3` |
| Local HEAD before this commit | `70e299d845b8204d63c0ae3e6906fa173e8ab424` |

`70e299d` is the unpushed M7-P0a quarantine commit. This audit is part of the local M7-P0b commit that follows it. That commit is not pushed.

M6 remains `CLOSED_WITH_LIMITATIONS`. `VALIDATION_REGISTRY.tsv` row `prjna884016` remains `ASSEMBLY_CONFIRMED`, five booleans `false`, status `PENDING`. `crates/**` has no diff. Pre-existing M6-P2A worktree files, including the modified `benchmark/real_world/scripts/stage_m6_p1a_fastq.sh`, are excluded from the M7-P0b commit.

## Orchestrator

- CLI: Cursor
- Model: Grok 4.7 (fallback; Composer 2.5 Standard was not the session model)
- Mode: orchestration only
- Task: preflight, dispatch, stage gates, ledger, narrow local commit
- Files written by the orchestrator: `M7_P0B_DEPENDENCY_INVENTORY.md`, `M7_P0B_SCIENTIFIC_REVIEW.md`, this audit. Implementation edits were made by ZCode.
- Tests after the final code and local-result revisions: `python3 -m pytest -q` on `test_m7_p0b_scoring.py`, `test_prjna884016_quarantine.py`, and `test_compare_variant_calls.py` — 31 passed, exit 0. `bash benchmark/real_world/scripts/check_validation_assets.sh` — exit 0.
- Verdict: accepted for a local commit. Human approval is still required before GitHub synchronization.

## CLI preflight

All three CLIs were callable. Versions observed: ZCode 0.17.0, Codex CLI 0.163.0, Claude Code 2.1.216. ZCode has no `--model` flag. `GLM-5.3` was selected with a temporary `ZCODE_PERSONAL_PROVIDER_CONFIG_FILE` whose `modelId` is `GLM-5.3` on provider `account:bigmodel-individual-coding-plan`. `GLM-5.3-Flash` is the default personal model and was invoked with that variable unset. DeepSeek was not selected. Codex catalog slug used: `gpt-6.1-sol` with `model_reasoning_effort="high"`. Claude model id used: `claude-opus-5-5` with `--effort high`.

| CLI | Provider | Exact model | Mode | Task | Invocation | Result |
| --- | --- | --- | --- | --- | --- | --- |
| zcode | bigmodel-individual-coding-plan | GLM-5.3-Flash | plan | smoke | `zcode --mode plan --disallowed-tools Write,Edit --json -p "/model"` | completed, session `sess_e974a203-617e-43fa-8674-5cd9e286b35e`. `/model` is not a headless model switch. |
| zcode | bigmodel-individual-coding-plan | GLM-5.3 | plan | smoke | `ZCODE_PERSONAL_PROVIDER_CONFIG_FILE=…/provider_config_glm53.json zcode --mode plan` | completed, session `sess_fb0105d3-6b0c-4e01-80bb-bad45ab8aadb` |
| codex | openai | gpt-6.1-sol | read-only, effort high | smoke | `codex exec -m gpt-6.1-sol -c 'model_reasoning_effort="high"' -s read-only --ephemeral` | `CODEX_SMOKE_OK`, exit 0, session `01a16750-85b8-7d53-b1d9-97d208938465`. An earlier `codex exec -a` failed because `-a` is not valid on `exec`. |
| claude | anthropic | claude-opus-5-5 | plan, effort high | smoke | `claude -p --model claude-opus-5-5 --effort high --permission-mode plan --tools "" --no-session-persistence` | `CLAUDE_SMOKE_OK`, exit 0 |

## Workflow A — dependency inventory

- CLI: zcode 0.17.0
- Provider / model: `account:bigmodel-individual-coding-plan` / `GLM-5.3-Flash`
- Mode: READ_ONLY, `--mode plan`, Write and Edit disallowed
- Session: `sess_ce6e2624-4ffc-4282-9273-efc45b171255`
- Invocation: `zcode --mode plan --cwd <repo> --disallowed-tools Write,Edit --no-color --json -p <workflow-a-prompt>`
- Exit: 0
- Files modified by the model: none. Orchestrator saved `benchmark/real_world/M7_P0B_DEPENDENCY_INVENTORY.md`.
- Tests: none (read-only)
- Verdict: inventory completed. It is not the scientific gate.

## Workflow B — scientific review

| Attempt | Status |
| --- | --- |
| `claude -p --model claude-opus-5-5 --effort high --permission-mode plan --tools "Read,Grep,Glob,Bash" --output-format json` | Not the review. Session `8546df6d-3c7b-4267-b690-641b2a42d5a4` produced no usable review and was stopped. `workflow-b.json` records `error_during_execution`. |

Completed review:

- CLI: claude 2.1.216
- Provider / model: Anthropic / `claude-opus-5-5`
- Mode: READ_ONLY. `--effort high`, `--permission-mode dontAsk`, allowed tools `Read,Grep,Glob,Bash`, disallowed `Edit,Write,NotebookEdit`, `--output-format stream-json`
- Session: `80f44b04-2177-4101-850e-58e5168d0a40`
- Exit: 0. Result subtype `success`, 12 turns, about 119 s. Tools actually used were Read and Bash. The model did not edit the worktree.
- Files: orchestrator saved `benchmark/real_world/M7_P0B_SCIENTIFIC_REVIEW.md`
- Verdict: **SCIENTIFIC_DESIGN_PASS**

Contract used for implementation: PRJNA884016 has no valid positive chromosomal structural truth. Chromosome assemblies support descriptive zero-difference concordance and do not by themselves establish tool specificity. With zero positive truth, recall, sensitivity, breakpoint accuracy, and positive-event TSD accuracy are `NOT_EVALUABLE`. Precision follows `tp+fp`. An empty call set makes precision `NOT_EVALUABLE`; unmatched calls against an assembly-confirmed zero-difference locus are false positives. Plasmid IS evidence stays separated. Registry row 7 stays `ASSEMBLY_CONFIRMED` with all validation flags false and status `PENDING`. REL606 remains a breseq comparator, not independent assembly truth.

## Workflow C — implementation

| Attempt | Model | Mode | Session | Outcome |
| --- | --- | --- | --- | --- |
| First | GLM-5.3 | `--mode build` | cancelled before a JSON session record (`workflow-c.json` is empty; stderr: turn cancelled) | Edit/Write denied. No worktree change. Not the implementation. |
| Second | GLM-5.3 | `--mode yolo` | `sess_fe0749cf-9a52-48ae-b180-5fa275c05697` | Exit 0. Scoring, quarantine guard, README, MUMmer sidecar quarantine, and `test_m7_p0b_scoring.py`. |
| Third | GLM-5.3 | `--mode yolo` | `sess_8d54d280-6efc-4752-b82b-0ae3185985f1` | Exit 0. Aggregation alias and fail-closed merge after Codex `NEEDS_REVISION`. |
| Fourth | GLM-5.3 | `--mode yolo` | `sess_28f9422d-daec-4ef9-bb65-6bb97b69fc2f` | Exit 0. Gitignored live result files only. No tracked-file change. |

GLM-5.3 invocations set `ZCODE_PERSONAL_PROVIDER_CONFIG_FILE` to the temporary GLM-5.3 provider config and used `zcode --mode yolo --cwd <repo> --disallowed-tools WebSearch,WebFetch --no-color --json -p …`. Runtime logs show only `account:bigmodel-individual-coding-plan/GLM-5.3`.

Tracked files changed:

- `benchmark/real_world/scripts/compare_variant_calls.py`
- `benchmark/real_world/scripts/compare_mob_truth.py`
- `benchmark/real_world/scripts/aggregate_results.py`
- `benchmark/real_world/scripts/check_validation_assets.sh`
- `benchmark/real_world/scripts/test_m7_p0b_scoring.py`
- `benchmark/real_world/README.md`
- `benchmark/real_world/truth/prjna884016/quarantined_legacy/MANIFEST.tsv`
- quarantined `SY22A1_vs_SY221_legacy.qdiff` and `.snps`
- deletion of active `SY22A1_vs_SY221.{delta,report,1delta,mdelta,1coords,mcoords,qdiff,snps}`

Gitignored local republication, not in the commit:

- `benchmark/real_world/results/prjna884016/summary.tsv`
- `benchmark/real_world/results/prjna884016/mob_metrics.tsv`
- `benchmark/real_world/results/prjna884016/STRUCTURAL_DISCREPANCIES.md`

Forensic copies `*.invalidated-m7p0a` retain the withdrawn text. Live score cells are `NOT_EVALUABLE`. A search of the three live files found no `1.0000`, `100%`, or `tp=4`.

## Workflow D — regression checks

| Attempt | Model | Mode | Session | Tests actually run by ZCode | Verdict |
| --- | --- | --- | --- | --- | --- |
| First | GLM-5.3-Flash | `--mode plan` | `sess_15e25df8-5c1d-47da-9fae-406aebc470f1` | no | plan only; not `CHECKS_PASS` |
| Second | GLM-5.3-Flash | `--mode build`, Edit/Write disallowed | `sess_92c5e5c3-3224-4e2d-b19f-883f881355a8` | no; Bash permission denied for pytest and the asset gate | `INCONCLUSIVE` |

Flash did not redesign scoring rules and did not execute the regression commands. The orchestrator ran them. After the aggregation revision the same three test modules were 31 passed, exit 0, and the asset gate was exit 0. The same pair was re-run after the gitignored result rewrite: 31 passed, gate exit 0.

## Workflow E — independent audit

Codex 0.163.0, provider openai, model `gpt-6.1-sol`, sandbox `read-only`, reasoning effort `high`. Invocation form:

`codex exec -m gpt-6.1-sol -c 'model_reasoning_effort="high"' -s read-only --ephemeral -C <repo> -o <file> <prompt>`

| Review | Session | Exit | Verdict | File-level finding |
| --- | --- | --- | --- | --- |
| 1 | `01a167ae-5d37-7e60-bae5-03e047066bcf` | 0 | `NEEDS_REVISION` | `aggregate_results.py` did not alias comparator keys; a later `1.0000` could replace `NOT_EVALUABLE`. Live `results/prjna884016/summary.tsv` still published perfect scores. |
| 2 | `01a167bb-a32d-7d50-9883-5347baeb8c67` | 0 | `NEEDS_REVISION` | Alias and merge defects accepted. Live `STRUCTURAL_DISCREPANCIES.md`, `summary.tsv`, and `mob_metrics.tsv` still published the withdrawn PASS / `tp=4` / perfect breakpoint claims. No additional blocking code defect. |
| 3 | `01a167c5-8888-7c20-af7d-42ead63a7fc4` | 0 | `PASS` | Live files now use `NOT_EVALUABLE` and withdraw the PASS claims. Tracked code still fail-closes empty positive truth. No diff to `crates/**` or `VALIDATION_REGISTRY.tsv`. Gitignored results will not travel with the commit; Codex judged that a verification limit, not a block. |

Codex did not re-run pytest inside the read-only sandbox. The orchestrator pytest exit codes above are the test evidence.

## Gate summary

| Gate | Result |
| --- | --- |
| Scientific design | `SCIENTIFIC_DESIGN_PASS` |
| Empty-positive metrics | `NOT_EVALUABLE`; precision keeps its own denominator |
| Quarantined truth | rejected before scoring |
| Valid positive-truth numeric scoring | preserved by regression tests |
| Registry flags | unchanged, all false, `PENDING` |
| Production crates | unchanged |
| Frozen M6 | not included in this commit |
| M7-P1 / Slurm / FASTQ download | not started |
| Engineering audit | `PASS` |

**M7_P0B_PUSH_READY**

Local commit only. Do not push until a human approves GitHub synchronization.
