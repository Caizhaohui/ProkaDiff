# 自主编排测试报告

命令在 `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff` 执行，主机 `bnode2`，Slurm `3095377`。

```text
python3 -m unittest benchmark/results/parity_closure/supervisor/test_orchestration.py -v
```

最终一次：退出码 0，`Ran 14 tests in 5.372s`，`OK`。时间约 2026-11-02T13:52+08:00。

前一次退出码 1。原因是 `--sleep` 被声明为 `int`，`0.05` / `0.1` 被 argparse 拒绝并返回 2；测试 12 在 CLI 探测结束前读取了尚不存在的日志。修复后 `--sleep` 为 `float`，队列路径不再先运行 `agent -v`，测试 12 等到 `QUEUE_IDLE`。修复后的这一次是验收结果。

## 哈希

| 文件 | sha256 |
| --- | --- |
| `parity_closure_supervisor.py` | `9c7008a890097a0a3b26c9bd1925f114c569b65f57cd2f39c6ad7211f028f44b` |
| `orchestration_runtime.py` | `3077ff07f0a95a6698935e6d96ca30a65c00aee5fc1b753acf42208bdb8f2065` |
| `test_orchestration.py` | `9dd642cbd5aa001a6a3e159942c8bee17223b1ed8b73a2056b267c2196a60297` |
| `dispatch_g1_opencode_wave11b0.sh` | `aae54d51756fdd7c0276ae677822228108af3b6c9c2742587d2c0e7c9f1948b8` |
| `dispatch_g1_opencode_wave11b2.sh` | `9e88d15b3142075b2aa0c7835a3ea1b31b71fcbc1952b5792817474f0e1c0579` |
| `dispatch_codex_wave11b_review.sh` | `fc408c2efaf815fa040521f09e81097a512ce23df4044fa540bc34b980abf2b1` |

## 十二项对照

测试使用临时仓库和假命令，不调用 OpenCode、Codex 或真实 Slurm。

| 要求 | 测试 | 结果 |
| --- | --- | --- |
| 1 worker 完成触发下一轮 | `test_1_worker_completion_starts_next_cycle` | 通过。A、B 均为 `COMPLETE`，`CONTROLLER_CYCLE` 至少 2 次 |
| 2 Claude 不可用仍能开始 B0 | `test_2_and_3_claude_block_does_not_stop_ready_work` | 通过。决策为 `WAVE11B0_SHADOW_AS` |
| 3 科学阻断不挡住无关 `READY` | 同上 | 通过。B0 与 B2 可运行，B1 与 E 不可运行 |
| 4 Slurm 完成触发收集和下一决策 | `test_4_slurm_completion_collects_and_dispatches` | 通过。stub 作业 42 `COMPLETED 0:0` |
| 5 主控崩溃不重复派发 | `test_5_crash_does_not_duplicate_worker` | 通过。重启后 `ADOPT`，`SPAWN` 一次 |
| 6 过期锁文本不是所有权 | `test_6_stale_lock_text_is_not_ownership` | 通过 |
| 7 多次调度后循环仍在 | `test_7_loop_survives_multiple_dispatches` | 通过。三条链完成后进程仍在，直到 `SIGTERM` |
| 8 重启从检查点继续 | `test_8_restart_resumes_checkpoint` | 通过 |
| 9 不启动 Gemini | `test_9_cursor_argv_never_selects_gemini` 与 `test_9_process_refuses_gemini_without_launching_it` | 通过。`agy` 返回 3，无 `AGENT_START` |
| 10 两路不能同持 writer 锁 | `test_10_second_writer_lease_fails` | 通过 |
| 11 失败 worker 不会卡死队列 | `test_11_failed_worker_is_terminal_and_does_not_stick` | 通过。`/bin/false` 为 `FAILED_TERMINAL`，兄弟任务 `COMPLETE` |
| 12 全部外部阻断时进程仍在 | `test_12_all_external_blocks_keep_the_loop_alive` | 通过。`QUEUE_IDLE`，无 `SPAWN` |

余额测试 `test_stale_claude_balance_is_not_inferred_from_spend` 同时通过：磁盘上的 2026-10-28 记录被判为 stale。

## 生产循环的现场证据

安装命令退出码 0：

```text
python3 benchmark/results/parity_closure/supervisor/parity_closure_supervisor.py --repo "$PWD" --install-queue
```

随后唯一循环 pid `3655965`。2026-11-02T13:53:46+08:00 日志：

```text
QUEUE_LOOP_START pid=3655965 once=False
CONTROLLER_CYCLE n=1 dispatch=WAVE11B0_SHADOW_AS ready=WAVE11B0_SHADOW_AS,WAVE11B2_3072550_REGRESSION
SPAWN WAVE11B0_SHADOW_AS pid=3655967 worker=opencode/glm-5.3
```

13:54 复查：`3655965`、wrapper `3655967`、`opencode` `3655968` 仍在；两把锁都是 `LOCK_HELD`；`squeue` 仍只有 `3095377` 与 `3095279`。进程树里没有 `agy`。OpenCode 会话 `ses_e652dd2e3ffeT7xDTx0NDUpZ5Y` 已读策略文件和 3066886 证据。

## 尚未完成的部分

合成多周期测试已经通过，生产循环在第一次调度后没有退出。

## 2026-11-02T15:29 复查（真实链路实况）

链路在 B2 处中断，不是代码缺陷，是 GLM 配额。实测：

| 时间 | 事件 |
| --- | --- |
| 13:53:46 | `CONTROLLER_CYCLE n=1` 派出 B0 |
| 14:17:46 | `WORKER_EXIT WAVE11B0_SHADOW_AS exit=0 status=COMPLETE` |
| 14:17:46 | `CONTROLLER_CYCLE n=2` 自动派出 B2（依赖链正常） |
| 14:36:16 | B2 exit 1，被判 `FAILED_RETRYABLE`，同一秒重派 |
| 14:38:16 | B2 exit 1，被判 `FAILED_TERMINAL` |
| 14:38–15:29 | `QUEUE_IDLE` 空转，Codex 依赖永不满足 |

OpenCode 原始错误：`provider.rate-limit`，`status 429`，`code 1308`，"已达到 5 小时的使用上限。您的限额将在 2026-11-02 16:49:43 重置"。

B0 已交付四个文件（`WAVE11B0_AS_SHADOW_RESULTS.tsv` 34 行、`WAVE11B0_READ_PROVENANCE.tsv`、`WAVE11B0_TIE_AUDIT.tsv`、`WAVE11B0_REGRESSION_REPORT.md`）和确定性分析器。关键行：JC330 28→0（`remain NO`）、JC331 35→29、JC344 28/31（`remain YES`）、JC326/321/343 均 25→25。与 Wave 11A Codex 复算一致。

## 限流分类修复（第二轮）

缺陷：普通失败和配额限流走同一条 `classify_exit()`，429 在 2 分钟内烧掉全部重试预算变成 `FAILED_TERMINAL`，一个上游限流就永久切断下游依赖。

修复：新增 `detect_rate_limit()`、`classify_rate_limited()`、`refresh_rate_limits()`。任务可声明 `rate_limit_logs`。非零退出时先检测限流标记（`provider.rate-limit`、`code 1308`、`使用上限`、`429`），命中则：

- 状态置 `BLOCKED_EXTERNAL`，`retry_count` 不增加
- `resume_condition` 写入 `provider_quota_reset_at=<从日志解析的重置时间>`
- 循环每轮调用 `refresh_rate_limits()`，到点自动回到 `READY` 并写 `REACTIVATE_RATE_LIMITED`

普通失败仍按原规则累加重试。

新增 `test_13_provider_rate_limit_is_external_block_not_terminal`：断言限流不消耗重试预算、阻断期间不可调度、重置时间前不复活、到点复活，同时普通 `/bin/false` 仍为 `FAILED_TERMINAL`。全量 15 项 `OK`，4.997 秒，退出码 0。

## 线上纠正与重启

- B2 从 `FAILED_TERMINAL` 纠正为 `BLOCKED_EXTERNAL`，`retry_count` 归 0，`resume_after=2026-11-02T16:49:43+08:00`
- 旧 supervisor 用 `SIGTERM` 优雅退出：`LOOP_SIGNAL signum=15` → `QUEUE_LOOP_END pid=3655965 cycles=194`，没有强杀
- 新 supervisor pid `3668146`，父进程 1，单一实例
- `orchestration_runtime.py` 新哈希 `0e3eef4ff41b74ba09a93230ecab57a9241433688f183314225ecb6a943bf7a0`
- `supervisor.lock` `LOCK_HELD`，`writer.lock` `LOCK_FREE`（当前无 writer）

## 真实链路的最终状态

- B0 完成，**未做** Codex 审查
- B2 待 16:49:43 配额重置后自动重派
- Codex 审查因依赖 B2 尚未开始
- 因此“Claude 不可用时自动完成 B0、Codex 审查和 3072550 调查并回到队列”仍未闭环。G1 仍是 `PROJECT_ACCEPTANCE_PENDING`。

## 多模型降级（第三轮，对齐 fastBact）

owner 指出 GLM 额度耗尽后应调用 Grok 4.7 继续开发。fastBact 依据：`AGENTS.md` 要求选模型前先查 `quota_manager.py status`；`F27_WORKER_DISPATCH_LEDGER.md` 的 `DISP-F27-002`/`003` 明确标注 *failover from GLM rate limit*，用 `cursor-agent` + `grok-4.7-high-fast` 完成生产任务；配额状态跨项目共享于 `~/.local/state/ai-routing/quota-state.json`。

新增能力：

- `fallback_worker` / `fallback_command` / `fallback_rate_limit_logs` / `escalation_history` 任务字段
- `escalate_worker()`：限流后把同一任务转给 Grok，不消耗 `retry_count`，记 `ESCALATE_WORKER`
- `quota_pool_status()` / `glm_pool_blocked()`：读共享配额状态，冷却未过即判 GLM 不可用
- `pre_escalate_if_pool_blocked()`：派发前预检，已知限流直接派 Grok
- 防回弹：`refresh_rate_limits()` 跳过 `fallback_worker_active`；同一 fallback 只升级一次

Grok 调度脚本用 `cursor-agent -p --output-format json --model grok-4.7-high-fast --trust --auto-review --workspace <repo>`，**不使用 `--force`**（AGENTS.md 与 V5 均禁止）。smoke 验证在 `/tmp/grok_smoke` 生成 `b.txt`（`SMOKE_OK`，退出 0），确认无 `--force` 也能写。

新增测试 14/15/16。全量 `Ran 18 tests in 4.668s`，`OK`，退出码 0。

现场：`quota-state.json` 为 `bigmodel-coding-plan: RATE_LIMITED`，`next_probe 2026-11-02T08:49:43Z`（=16:49:43+08）。B2 已立即改派 Grok，supervisor `3669202` 于 15:36:54 记录 `SPAWN WAVE11B2_3072550_REGRESSION worker=cursor/grok-4.7-high-fast`，cursor-agent pid `3669205` 运行中，进程树内无 `agy`。

哈希：

| 文件 | sha256 |
| --- | --- |
| `orchestration_runtime.py` | `c653b016d79f07dcf2c9d87d014a4b10827d80acb8bc691b1ec5b3d2a367e6dd` |
| `test_orchestration.py` | `70b53aff91163cd1967e2052a7649df1ab942fbf2f89219847062b03796936b9` |
| `dispatch_g1_grok_wave11b0.sh` | `c617d0952b6bf6165792fef7131c2b0d90d09dc7f10dd6e1c42149e72f11de2f` |
| `dispatch_g1_grok_wave11b2.sh` | `4197ac3d085accbdc05fa34ef2cd009c93e3f71ea6a3509313c27e8268f184da` |

此轮之后 B2 不再等 16:49:43，而是由 Grok 立即接手。Codex 仍等待 B0 与 B2 均完成。G1 仍为 `PROJECT_ACCEPTANCE_PENDING`。
