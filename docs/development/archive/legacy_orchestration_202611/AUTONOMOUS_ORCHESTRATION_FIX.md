# 自主编排修复

修复范围只在 supervisor、队列运行时、调度脚本和状态文件。没有修改 `crates/` 里的生产归因，没有 `sbatch`，没有改验收阈值。

## 运行时

新增 `benchmark/results/parity_closure/supervisor/orchestration_runtime.py`。`parity_closure_supervisor.py` 的 `--loop`、`--once`、`--install-queue` 都进入 `run_queue_loop()`。无这些参数时返回 2。`--sleep` 改为浮点。锁冲突返回 75，不再把 `SystemExit(0)` 当成成功。

循环在持有 `supervisor.lock` 的 `flock` 期间运行。一个 worker 结束后回到：

`reconcile` → 重建可运行队列 → 记录 `CONTROLLER_CYCLE` → 只派出一个 writer → `sleep` → 再评估。

`QUEUE_LOOP_END` 只出现在 `stop`、信号、`PARITY_MAX_CYCLES` 或 `--once`。生产循环不设置 `PARITY_MAX_CYCLES`。暂停时写 `PAUSED` 并继续睡眠，不退出。

## 任务字段与状态

每个任务含 `task_id`、`task_type`、`worker`、`status`、`dependencies`、`blocked_reason`、`resume_condition`、`allowed_while_other_tasks_blocked`、`retry_count`、`maximum_retries`、`execution_pid`、`output_paths`、`source_hash`、`completion_timestamp`，另有调度用的 `priority`、`command`、`exit_code`、`slurm_job_id`。

状态：`READY`、`RUNNING`、`WAITING_WORKER`、`WAITING_SLURM`、`BLOCKED_EXTERNAL`、`REVIEW_REQUIRED`、`COMPLETE`、`FAILED_RETRYABLE`、`FAILED_TERMINAL`。

可运行条件：状态为 `READY` 或 `FAILED_RETRYABLE`，依赖均为 `COMPLETE`，命令非空，重试次数未到上限。同工作树同时只持有一把 `writer.lock`。没有命令的任务不会被自动执行，因此 Claude 任务和 Slurm 任务目前不会被误派。

退出分类：0 为 `COMPLETE`；非 0 增加 `retry_count`，达到 `maximum_retries` 为 `FAILED_TERMINAL`，否则 `FAILED_RETRYABLE`。默认 `maximum_retries` 为 2。

## 已安装队列

`--install-queue` 于 2026-11-02T13:53:40+08:00 退出 0。安装前状态哈希 `8f5c886e1d3004dcf0164e41d8c33d85a9a04f5dc9fa5fa4719b9481711a7474`。`wave11a_codex_verdict=MECHANISM_SUPPORTED` 与 `claude_dispatch=NOT_INVOKED_BALANCE_STALE` 保留。治理标记改为 `PROKADIFF_MULTI_AI_V5_QUEUE`。

| 任务 | 优先级 | 初始状态 | 说明 |
| --- | --- | --- | --- |
| `WAVE11B0_SHADOW_AS` | 10 | `READY` | OpenCode GLM-5.3，只写 `wave11b0/` |
| `WAVE11B2_3072550_REGRESSION` | 20 | `READY` | 等 B0 释放 writer 锁后运行 |
| `WAVE11B1_SCIENTIFIC_REVIEW` | 30 | `BLOCKED_EXTERNAL` | 余额原因 `claude_balance_stale_age_hours=122.9`，无命令 |
| `WAVE11D_CODEX_REVIEW` | 40 | `READY` | 依赖 B0 与 B2，只读 Codex |
| `WAVE11E_PRODUCTION_AS` | 50 | `BLOCKED_EXTERNAL` | 依赖 B0、B1、D，无生产编辑命令 |
| `WAVE11F_SLURM_VALIDATION` | 60 | `READY` | 依赖 E，命令为空，不会自动 `sbatch` |

Claude 判断只读 `~/.local/state/ai-routing/claude-budget.json`。该文件 `last_verified_balance_at_utc=2026-10-28T03:00:00Z`，余额 `45.886322`，保留金 `10`。超过 24 小时，不使用本地支出账本推断余额。只阻断 B1。

## 锁与单写者

`try_lease()` 使用 `LOCK_EX|LOCK_NB`。过期 pid 文本不能取得所有权。2026-11-02T13:54 实测 `supervisor.lock` 与 `writer.lock` 均为 `LOCK_HELD`，文本 pid `3655965`。恢复时收养活着的 `execution_pid`，不杀 worker，不第二次 `SPAWN`。

调度脚本不再改 `parity_closure_state.json`，也不删除 `writer.lock`。OpenCode 与 Codex 使用绝对路径。

## Cursor 与 Gemini

队列循环不调用 `agy`。`PARITY_AGENT_CLI=agy` 或模型名含 `gemini` 时进程返回 3，且日志没有 `AGENT_START`。确定性队列不在每个周期再启动一个 `agent` 进程；需要显式询问时，唯一构造函数是 `cursor_launch_argv()`，模式为 `plan`，输出为 `json`。

## 生产循环

2026-11-02T13:53:46+08:00 将 `pause` 设为 false 后启动唯一循环：

`setsid nohup python3 benchmark/results/parity_closure/supervisor/parity_closure_supervisor.py --repo <repo> --loop --sleep 30`

- supervisor pid `3655965`，父进程 1，主机 `bnode2`
- `QUEUE_LOOP_START` 后 `CONTROLLER_CYCLE n=1 dispatch=WAVE11B0_SHADOW_AS`
- `SPAWN` pid `3655967`，OpenCode pid `3655968`，模型 `zhipuai-coding-plan/glm-5.3`
- 同时刻 B2 仍是 `READY`，B1 与 E 仍是 `BLOCKED_EXTERNAL`
- 日志没有 `QUEUE_LOOP_END`，没有新的 Slurm 作业

交互中的 Cursor 会话不再派发 worker。不要再启动第二个 `--loop`。

## 状态刷新（2026-11-02T15:31+08:00）

真实链路在 B2 处被 GLM 配额阻断，不是代码问题。OpenCode 原始错误 `provider.rate-limit` / `status 429` / `code 1308`："已达到 5 小时的使用上限。您的限额将在 2026-11-02 16:49:43 重置"。

旧 `classify_exit()` 把限流当成普通失败，14:36:16 与 14:38:16 两次失败烧完全部重试预算，B2 变成 `FAILED_TERMINAL`，Codex 依赖永远不满足，队列空转约 50 分钟。这与 V5 要求的“一次上游限流不得永久切断下游依赖”冲突。

修复内容：

- `detect_rate_limit()` 读任务声明的 `rate_limit_logs` 尾部，匹配 `provider.rate-limit`、`code 1308`、`使用上限`、`429 Too Many Requests`，再解析重置时间
- `classify_rate_limited()` 置 `BLOCKED_EXTERNAL`，`retry_count` 不增，写 `resume_condition=provider_quota_reset_at=...` 与 `resume_after`
- `refresh_rate_limits()` 每轮检查，到点写 `REACTIVATE_RATE_LIMITED` 并回到 `READY`
- 普通失败仍走 `classify_exit()` 累加重试
- `task_template()` 新增 `rate_limit_logs`、`resume_after` 字段；B0/B2 任务登记各自 stdout 路径

测试：新增 `test_13_provider_rate_limit_is_external_block_not_terminal`。全量 `Ran 15 tests in 4.997s`，`OK`，退出码 0。运行时哈希 `0e3eef4ff41b74ba09a93230ecab57a9241433688f183314225ecb6a943bf7a0`，测试哈希 `45f5d2493b2ad258d4ec2c7ba28f346817c60c9986c5f7cfbb220033e861decc`。

线上纠正：

- B2 `FAILED_TERMINAL` → `BLOCKED_EXTERNAL`，`retry_count` 归 0，`resume_after=2026-11-02T16:49:43+08:00`
- 旧 supervisor `3655965` 以 `SIGTERM` 优雅退出：`LOOP_SIGNAL signum=15`，`QUEUE_LOOP_END cycles=194`，未强杀
- 新 supervisor pid `3668146`，父进程 1，唯一实例
- 15:31 实测 `supervisor.lock` `LOCK_HELD`，`writer.lock` `LOCK_FREE`

B0 测量（影子分析，不改生产）：JC330 28→0、JC331 35→29、JC326/321/343 25→25、JC344 28/31 保留。与 Wave 11A Codex 复算一致。

Codex 审查尚未开始，真实链路未闭环。G1 仍为 `PROJECT_ACCEPTANCE_PENDING`。

## 多模型降级（对齐 fastBact）

owner 指出：GLM 额度耗尽后应调用 Grok 4.7 继续开发，而不是空等配额窗口。fastBact 的做法可直接对照：

- `AGENTS.md`：选择模型前先查 `python3 ~/.local/share/ai-routing/quota_manager.py status`，`UNKNOWN` 视为未确认
- `F27_WORKER_DISPATCH_LEDGER.md`：`DISP-F27-002` / `DISP-F27-003` 明确标注 *failover from GLM rate limit*，用 `cursor-agent` + `grok-4.7-high-fast` 完成生产任务
- 配额状态由跨项目共享的 `~/.local/state/ai-routing/quota-state.json` 管理，不按项目各存一份

ProkaDiff 原先只有"限流转外部阻断 + 等重置"，缺了降级这一层。现补齐：

1. `task_template()` 增加 `fallback_worker`、`fallback_command`、`fallback_rate_limit_logs`、`escalation_history`
2. `escalate_worker()`：命中限流后**把同一任务转给 Grok**，不消耗 `retry_count`，记 `ESCALATE_WORKER`，并写入升级历史（含原因与配额重置时间）
3. `quota_pool_status()` / `glm_pool_blocked()`：读共享 `quota-state.json`，冷却未过即判 GLM 不可用
4. `pre_escalate_if_pool_blocked()`：**派发前**预检，已知 GLM 限流就直接派 Grok，不浪费一次调度
5. 防止回弹：`refresh_rate_limits()` 跳过 `fallback_worker_active` 的任务，`escalate_worker()` 对同一 fallback 只升级一次

新增脚本：

- `dispatch_g1_grok_wave11b0.sh`
- `dispatch_g1_grok_wave11b2.sh`

调用形式为 `cursor-agent -p --output-format json --model grok-4.7-high-fast --trust --auto-review --workspace <repo>`。**不使用 `--force`**：AGENTS.md 与 V5 指令均禁止；实测 `--trust --auto-review` 已能在无 `--force` 下写文件（smoke test 在 `/tmp/grok_smoke` 生成 `b.txt`，内容 `SMOKE_OK`，exit 0）。

新增测试：

- `test_14_glm_quota_escalates_to_grok_instead_of_waiting`
- `test_15_production_queue_declares_grok_failover`
- `test_16_known_glm_rate_limit_skips_glm_before_dispatch`

全量 `Ran 18 tests in 4.668s`，`OK`，退出码 0。

现场验证：`quota-state.json` 当前 `bigmodel-coding-plan: RATE_LIMITED`，`next_probe 2026-11-02T08:49:43Z`（= 16:49:43+08）。B2 已立即改派 Grok，supervisor `3669202` 于 15:36:54 记录 `SPAWN WAVE11B2_3072550_REGRESSION worker=cursor/grok-4.7-high-fast`，cursor-agent pid `3669205` 正在运行，进程树内无 `agy`。

哈希：

| 文件 | sha256 |
| --- | --- |
| `orchestration_runtime.py` | `c653b016d79f07dcf2c9d87d014a4b10827d80acb8bc691b1ec5b3d2a367e6dd` |
| `test_orchestration.py` | `70b53aff91163cd1967e2052a7649df1ab942fbf2f89219847062b03796936b9` |
| `dispatch_g1_grok_wave11b0.sh` | `c617d0952b6bf6165792fef7131c2b0d90d09dc7f10dd6e1c42149e72f11de2f` |
| `dispatch_g1_grok_wave11b2.sh` | `4197ac3d085accbdc05fa34ef2cd009c93e3f71ea6a3509313c27e8268f184da` |
