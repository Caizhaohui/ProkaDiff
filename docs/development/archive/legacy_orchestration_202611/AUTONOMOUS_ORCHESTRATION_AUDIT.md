# 自主编排运行时审计

审计时间：2026-11-02T13:53+08:00。主机 `bnode2`，Slurm `3095377`（`qcpu_18i`，1 CPU）。HEAD `1a0c7b553753bfa1fd6f0c0741c0a8bd43d4cf3c`。证据来自 `parity_closure_supervisor.py`、`supervisor.log` 和当时的控制文件，不以旧的 PASS 记录代替。

通用策略文件仍是 `~/.config/ai-collaboration/POLICY.md` V3.2，默认主控写的是 Antigravity / Gemini。本仓库的所有者指令把智能主控改为 Cursor CLI / Grok 4.7 High Fast。冲突记在 `HANDOFF.md`，科学契约未改。

## 1. Supervisor 为什么停着

`supervisor_control.json` 在 2026-11-02T12:20+08:00 被写成 `pause: true`，原因 `cursor_v4_takeover_lifecycle_only`。当时没有存活的 `parity_closure_supervisor.py`。暂停是为了阻止已死亡的旧循环重新拉起 `agy`。`supervisor.lock` 文本停在 2026-10-29，不是内核锁。日志里 2026-10-29T14:48:37 起反复出现 `DECIDE action=PAUSE` 然后 `PAUSED operator control`。

## 2. 循环在调度之后是否还在

旧生产路径是 `supervisor_once()`。它做完一次 `DECIDE` 就 `return 0`。`--loop` 若还在，只是睡醒后再进同一函数；进程本身不保证还在跑。2026-11-02 接管时没有 supervisor 进程，所以 worker 完成后没有人继续。

## 3. `RUN_CONTROLLER` 是否提前返回

是。`supervisor_once()` 在派出 Wave 7/8 或跑完一轮 controller 后直接返回。日志：

- 2026-10-29T14:08:19 `GROK_DISPATCH pid=1375722`
- 2026-10-29T14:13:19 状态仍是 `WAITING_FOR_WORKER`，动作却是 `RUN_CONTROLLER`，并启动了 `controller_round_1_1793254399`

一次调度的返回值结束了这一轮，没有回到“worker 完成后的下一任务”。

## 4. `--once` 和 `--loop` 是否分路

旧实现里，不带 `--loop` 就是一次 `supervisor_once()`。`--loop` 只是在同一次决策函数外包了睡眠；`PAUSE`、`WAIT_WORKER`、`TERMINAL`、`STALLED` 都返回 0，失败的 controller 返回非 0 时旧循环退出。二者不是两条持久语义。当前代码里，无 `--loop`、`--once`、`--install-queue` 的调用写 `REFUSING_IMPLICIT_ONESHOT` 并返回 2。`--loop` 与 `--once` 进入同一个 `run_queue_loop()`，只有 `once` 标志不同。

## 5. Cursor CLI 是否以非交互方式调用

历史日志实际启动过 `model=gemini-3.8-flash-medium`（2026-10-29T14:46:00 和 14:49:38）。旧默认 CLI 是 `agy`。当前 `run_agent_controller()` 和 `cursor_launch_argv()` 构造的是 `agent -p --output-format json --mode plan|ask --model grok-4.7-high-fast`，没有 `--force`。`agy`、`gemini`、`antigravity` 在执行前被拒绝并返回 3。队列循环本身按优先级做确定性选择，不在每个周期再起一个 Cursor 进程，避免与交互会话形成两个主控。

## 6. JSON 是否解析正确

旧解析只认逐行 `type=result` 或同时含 `conversation_id` 与 `response` 的事件。当前 `parse_agent_stdout()` 先接受整份 JSON 文档，再接受逐行事件，并归一化 `result` / `response` 与 `session_id` / `conversation_id`。`test_json_result_objects_parse` 通过。

## 7. Worker 完成是否唤醒下一轮

旧全局状态做不到。`worker_pids` 为空时，即使状态是 `WAITING_FOR_WORKER`，`decide_next_action()` 也会落到 `RUN_CONTROLLER`（日志第 9–10 行）。新队列在同一进程里 `reconcile()` 后重新 `eligible_tasks()`。测试 1 中 A 完成后自动调度 B，日志里至少两次 `CONTROLLER_CYCLE`。

## 8. Slurm 完成是否唤醒下一轮

旧逻辑只有一个全局 `slurm_id`。`sacct` 结果既非 running 也非 terminal 时，`WAITING_FOR_SLURM` 仍会变成 `RUN_CONTROLLER`。日志 2026-10-29T14:49:38 就是这个状态。新队列把任务标成 `WAITING_SLURM`，只有 `sacct` 或 `PARITY_SLURM_STUB` 给出终态才 `COLLECT_SLURM` 并进入下一决策。测试 4 通过。本阶段没有提交新的 G1 作业。

## 9. 一个 BLOCKED 任务是否停掉全部工作

旧 `decide_next_action()` 把 `HARD_BLOCKED`、`ACCEPTED`、`PROJECT_SUCCESS` 变成 `TERMINAL`，并把连续无进展变成 `ENGINEERING_STALLED`。唯一的 `next_action` 被科学阻断占用后，其他工作没有独立队列。Wave 11A 之后的 `CLAUDE_BALANCE_REFRESH_BEFORE_WAVE11B` 就是这个瓶颈。新队列里 `WAVE11B1_SCIENTIFIC_REVIEW` 单独为 `BLOCKED_EXTERNAL` 时，`WAVE11B0_SHADOW_AS` 与 `WAVE11B2_3072550_REGRESSION` 仍可入选。测试 2、3、12 通过。

## 10. 失败 Worker 是否永久停在 `WAITING_FOR_WORKER`

会。旧代码只在 `worker_pids_alive()` 时等待；pid 文本消失就改去启动 controller，pid 仍显示存活则永远等待。Wave 7/8 的 writer 锁只包住 `Popen`，worker 活着时锁已经释放。新代码用 `fcntl.flock`，worker 退出写入 `state/task_exits/<id>.exit`，`classify_exit()` 分成 `COMPLETE`、`FAILED_RETRYABLE` 或 `FAILED_TERMINAL`。崩溃恢复会收养仍活着的 `execution_pid`，不第二次派出。测试 5 和 11 通过。

## 结论

停止的直接原因是：进程已经不在，控制文件仍暂停，而旧决策函数又把“一个任务结束或被阻断”当成整条流水线的终点。队列、锁和 Claude 余额规则已按本次审计改到 `orchestration_runtime.py`。G1 仍是 `PROJECT_ACCEPTANCE_PENDING`。
