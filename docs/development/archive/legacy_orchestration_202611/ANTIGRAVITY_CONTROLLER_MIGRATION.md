# ProkaDiff — Antigravity Multi-AI Controller Migration Report

- **Date:** 2026-10-29
- **Controller Platform:** Antigravity CLI (`agy` v1.4.10)
- **Controller Model:** Gemini 3.8 Flash Medium (`gemini-3.8-flash-medium`)
- **Authority / Policy:** Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **Source Reference:** fastBact (`/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact`)
- **Target Project:** ProkaDiff (`/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`)
- **Active Task:** G1 Clonal JC Parity Closure

---

## 1. fastBact 多 AI 协同体系的实际结构审查

经对源项目 `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact` 及其引用的全局规范 `~/.config/ai-collaboration/POLICY.md`、`~/.local/share/ai-routing/`、`HANDOFF.md` 与历史提交审计，核心机制如下：

1. **Gemini 3.8 Flash 担任 Controller**：
   - 依据 `POLICY.md` 规定，Antigravity CLI（`agy`）运行 `gemini-3.8-flash-medium` 作为默认唯一的 `AUTONOMOUS_EXECUTION_CONTROLLER`。
   - 负责任务规划、分解、Developer/Reviewer 严格顺序调度、门禁组织与持续接续推进，原则上不直接编写底层生产代码。
2. **多 AI CLI 统一调度规范**：
   - **Developer P1**：OpenCode CLI，`opencode run --standalone -m zhipuai-coding-plan/glm-5.3-flash`（常规）或 `glm-5.3`（复杂算法），`--format json`。
   - **Developer P2**：Cursor CLI，`agent -p --output-format json --model grok-4.7-high-fast` 或 `grok-4.7-xhigh-fast`，仅当 GLM 额度用尽或遇特定重大瓶颈时由 Controller 调度。
   - **Independent Reviewer**：Codex CLI，`codex exec -m gpt-6.1-sol`（核心正确性）或 `gpt-6.1-terra`（常规），必须加 `-s read-only --ephemeral`。
   - **Formal Validator**：Codex CLI，`codex exec -m gpt-6-luna -s read-only --ephemeral`，仅基于退出码、真实测试和输出产物做不可篡改的法定验收。
   - **Scientific Architect**：Claude Code CLI，`claude -p --model claude-opus-5-5 --effort high`，受严格财务保护（任务软预算 $5.00，网关保留金 $10.00）。
3. **Developer 与 Reviewer 的绝对独立性**：
   - 写入 Worker（GLM/Grok）严禁担任审查者；审查者（Codex Sol/Terra）与验收者（Luna）严格以只读沙箱运行，绝不直接修改代码。
   - 审查结论按缺陷等级（P0/P1/P2）形式化输出并关联可验证反例与哈希。
4. **长周期持续推进机制**：
   - 状态机驱动：`PLAN → IMPLEMENT → REVIEW → VALIDATE → NEXT_MILESTONE`。
   - 普通门禁通过后自动进入下一任务，不输出挂起会话的提问，仅在预算耗尽、硬权限/科学契约变更或目标全部达成时停止。
5. **Slurm 计算任务管理**：
   - 登录节点仅用于轻量元数据管理与作业提交调度。
   - 编译、测试、WGS 对拍基准全部通过 Slurm 提交计算节点（`qfree` 或 `qcpu_18i`）。
   - 记录 `SLURM_JOB_ID`，采用 `sacct` 异步轮询非阻塞状态机。
6. **中断恢复与断点接续**：
   - 每次状态转移通过原子写更新 `HANDOFF.md`、`ACTIVE_TASK.md` 与状态 JSON。
   - 恢复启动时，先检查锁与活动 PID，校验 Git SHA 与 dirty 状态，读取最新 checkpoint，无缝继续 `next_action`。
7. **单写入者契约（Single-Writer Invariant）**：
   - 同一 worktree 同一时间只能有一个活动写入者。
   - 通过 OS 文件锁（`fcntl.flock`）配合 PID 活性核验，避免并发双写冲突。
8. **客观验收与反欺骗门禁**：
   - 严禁采信 Worker 自我汇报（Self-Reported PASS）。
   - 必须基于真实测试套件退出码 `0:0`、独立 Oracle 对拍、产物哈希及 Codex Luna 只读验收。

---

## 2. 迁移文件映射与最小改动原则

| 源 / 目标文件 | 迁移处理 | 说明 |
| :--- | :--- | :--- |
| `~/.config/ai-collaboration/POLICY.md` | 引用生效 | 目标仓库通过 `AGENTS.md` 显式加载，作为通用自主协作标准 V3.2。 |
| `AGENTS.md` | 原位最小合并 | 移除历史 Cursor 限制，确立 Antigravity Gemini 3.8 Flash 为 Controller，更新多模型分工表，完全保留 ProkaDiff 15 条科学纪律与 clean-room 约束。 |
| `HANDOFF.md` | 原位更新 | 治理标识更新为 `PROKADIFF_MULTI_AI_V2_0`，确立 `agy` Controller 身份，记录 G1 Wave 7 完成状态与后续动作。 |
| `ACTIVE_TASK.md` | 原位更新 | 更新当前 Controller 为 Antigravity Gemini 3.8 Flash，记录任务从 Wave 7 完成转入 Slurm 对拍重跑。 |
| `parity_closure_supervisor.py` | 最小适配增强 | 默认 `AGENT_CLI="agy"`，`CONTROLLER_MODEL="gemini-3.8-flash-medium"`，适配 `agy` 命令行参数与 JSON 输出解析，保留对 Cursor 的向后兼容。 |
| `parity_closure_state.json` | 原子更新 | 记录 Controller 平台迁移，清空已退出 Worker PID，推进下一步 `CONTROLLER_SBATCH_G1_RERUN`。 |

**未迁移内容及原因：**
- 未复制 fastBact 的基因组组装测试数据、QUAST 评测流程及 Rust 组装器代码（科学领域不同）。
- 未复制 fastBact 历史任务 ID、Slurm 作业编号及特定 C0-02R 状态。
- 未复制任何 API Key、Token 或认证文件（全部复用既有系统级环境配置）。

---

## 3. 旧 Cursor Supervisor 安全交接核查

在执行迁移前进行了严格的只读核查：
1. **旧 Supervisor 进程**：发现 PID 1375608（`parity_closure_supervisor.py --loop`）曾处于运行态。
2. **安全交接操作**：立即调用受支持的控制命令 `bash benchmark/results/parity_closure/supervisor/parity_closure_ctl.sh pause`，将 `supervisor_control.json` 设置为 `pause: true`。
3. **活动写入者保护**：核查发现 Grok 4.7 XHigh Fast 正在执行 Wave 7（PID 1375794），严格遵守单写入者原则，未做强杀，耐心等待其完成。
4. **安全交接点到达**：
   - Grok Wave 7 于 14:35:37 顺利完成退出，提交并在 Slurm 计算节点完成单元测试（Job 3064893，exit 0:0，179 pass / 0 fail），生成 `FIX_NOTES_GROK_WAVE7.md`。
   - 退出脚本将 `parity_closure_state.json` 更新为 `WORKER_COMPLETED`，`next_action: CONTROLLER_SBATCH_G1_RERUN`。
   - 由于 Supervisor 处于 PAUSED 状态，无任何新进程竞争启动，实现了零冲突的安全平稳交接。

---

## 4. 迁移验证结果矩阵

### Test A — Model Dispatch CLI Verification

所有模型与工具 CLI 均在运行环境中完成真实探测：
- **Gemini Flash Controller** (`agy`): v1.4.10，`gemini-3.8-flash-medium` 探测成功，JSON 输出测试通过（exit 0）。
- **OpenCode Developer** (`opencode`): v2.2.1，支持 `zhipuai-coding-plan/glm-5.3-flash` 与 `glm-5.3`。
- **Cursor Grok Developer** (`agent`): v2026.10.22，支持 `grok-4.7-high-fast` 与 `grok-4.7-xhigh-fast`。
- **Codex Reviewer / Validator** (`codex`): v0.163.0，配置支持 `gpt-6-luna`、`gpt-6.1-sol`、`gpt-6.1-terra`。
- **Scientific Architect** (`claude`): v2.1.216 (Claude Code)，支持 `claude-opus-5-5`。

### Test B — Checkpoint Recovery

成功完整恢复当前状态：
- Task: `G1_CLONAL_JC`
- Git SHA: `1a0c7b553753bfa1fd6f0c0741c0a8bd43d4cf3c`
- Status: `WAITING_FOR_SLURM`
- Protected Floor: Slurm `3063624` (matched=6, over=22, under=12, MOB 5/5)
- Previous Run: Slurm `3064609` (matched=11, over=28, under=7)
- Worker: Wave 7 completed with exit code 0
- Next Action: `CONTROLLER_SBATCH_G1_RERUN`

### Test C — Autonomous Continuation

`parity_closure_supervisor.py` 已适配 Antigravity CLI：
- `--verify-cli` 成功返回 `agent: agy`，`version: 1.4.10`，`supports_print: true`，`supports_json: true`。
- Supervisor 接续测试与循环控制逻辑验证完毕。

### Test D — PAUSE / RESUME Gate

- `parity_closure_ctl.sh pause` 成功将控制标志置为暂停，决策层正确返回 `PAUSE`。
- `parity_closure_supervisor.py --pause-test` 返回 `PAUSE_RESUME_TEST_PASS`（exit 0）。

### Test E — Single Writer Gate

- 验证了内核级 `flock` 锁文件（`supervisor.lock`, `controller.lock`, `writer.lock`）互斥机制。
- 确认当有活动 Worker 或 Controller 时，互斥锁立即生效，防止双写。

---

## 5. 迁移结论与门禁判定

```
ANTIGRAVITY_CONTROLLER_MIGRATION_PASS
AUTO_CONTINUATION_PASS
SINGLE_WRITER_GATE_PASS
CHECKPOINT_RECOVERY_PASS
```

**正式接管状态：** Antigravity CLI（Gemini 3.8 Flash Medium）已正式接管 ProkaDiff 项目的 Controller 职责。

---

## 6. 回滚程序说明

若未来因特定环境原因需回滚至 Cursor Composer 2.5 Controller：
1. 设置环境变量：`export PARITY_AGENT_CLI="agent"`，`export PARITY_CONTROLLER_MODEL="composer-2.5"`。
2. 将 `AGENTS.md` 中的 Controller 默认项切换回 `agent / composer-2.5`。
3. `parity_closure_supervisor.py` 内建多 CLI 兼容层，将自动识别并切回 Cursor 命令行语法。
