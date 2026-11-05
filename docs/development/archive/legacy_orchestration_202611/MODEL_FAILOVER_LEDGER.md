# 模型降级台账（V5.1）

记录 Worker 选择、配额判断与升级事实。只写可核实的内容；不记录任何 API key、OAuth 凭据或 token。

## 路由规则（对齐 fastBact）

| 顺序 | 角色 | CLI / 模型 | 用途 |
| --- | --- | --- | --- |
| 1 | Primary Developer | `opencode` / `zhipuai-coding-plan/glm-5.3` | Rust、取证、回归、基准准备、文档 |
| 2 | Fallback Developer | `cursor-agent` / `grok-4.7-high-fast` | GLM 限流时接管；或任务需要更强推理 |
| 3 | Independent Auditor | `codex exec` / `gpt-6.1-sol` | 只读独立审查 |
| 4 | Scientific Reviewer | `claude` / `claude-opus-5-5` | 条件性，需核实余额 |

共享配额状态：`~/.local/state/ai-routing/quota-state.json`，由 `~/.local/share/ai-routing/quota_manager.py` 管理。不另建配额系统。

## 触发条件

派发前（`pre_escalate_if_pool_blocked`）：

1. 查 `bigmodel-coding-plan` 状态。
2. `AVAILABLE` → 用 GLM。
3. `RATE_LIMITED` 且冷却未过 → **直接派 Grok**，不浪费一次调度。
4. `UNKNOWN` → 不当作可用。
5. 两个都不可用 → 仅该任务 `BLOCKED_EXTERNAL`。

运行中（`detect_rate_limit` + `escalate_worker`）：

- 识别 `provider.rate-limit`、`code 1308`、`使用上限`、`429 Too Many Requests`。
- 命中后：保留原 task_id 与规格 → 标 `worker_rate_limit` → **不消耗 `retry_count`** → 同一任务升级到 Grok → 记录 `ESCALATE_WORKER` → 继续执行，不等 GLM 恢复。
- `refresh_rate_limits()` 跳过 `resume_condition=fallback_worker_active`，防止过早切回 GLM。
- Grok 不会再升级到自身（`escalate_worker` 检查 `escalation_history`）。

## 实际升级记录

### ESC-001 — `WAVE11B2_3072550_REGRESSION`

| 项 | 值 |
| --- | --- |
| 原 Worker | `opencode/glm-5.3` |
| 升级 Worker | `cursor/grok-4.7-high-fast` |
| 原因 | `provider_rate_limit` |
| 触发证据 | `provider.rate-limit` / `status 429` / `code 1308`，"已达到 5 小时的使用上限。您的限额将在 2026-11-02 16:49:43 重置" |
| 发生时间 | 2026-11-02T15:36:36+08:00 |
| 原 Worker 尝试次数 | 2 次（14:36:16、14:38:16），均 exit 1 |
| 重试预算消耗 | 已修正为 0；原误判为 `FAILED_TERMINAL` 后人工纠正 |
| Grok 结果 | exit 0，`2026-11-02T15:51:24+08:00`，`duration_ms=835307` |
| 产出 | `wave11b2/WAVE11B2_REGRESSION_REPORT.md`（8740 字节） |
| 结论 | breseq 114 消失在多拷贝代表位点选择，不是打分/聚类/频率过滤 |

### ESC-002 — 无

`WAVE11B0_SHADOW_AS` 全程在 GLM 上完成，无需升级。

### ESC-003 — 无

`WAVE11B0R_ANALYZER_REPAIR` 于 2026-11-03T11:41 派发时 GLM 已 `AVAILABLE`，直接使用 GLM。

## 配额快照

| 时间 | bigmodel-coding-plan | antigravity-gemini |
| --- | --- | --- |
| 2026-11-02T14:54 | `RATE_LIMITED`（`next_probe 08:49:43Z`） | `UNKNOWN` |
| 2026-11-03T11:39 | `AVAILABLE`（`evidence=RECOVERY`） | `UNKNOWN` |

`antigravity-gemini` 在本架构中不使用。V5 明确禁止 Antigravity/Gemini Controller。

## Claude

| 项 | 值 |
| --- | --- |
| CLI | `/hpcfs/fhome/caizhh/.local/share/omniakey/bin/claude` 2.1.216 |
| `auth status` | `loggedIn: false`，`authMethod: none` |
| 最后核实余额 | `45.886322 USD` @ `2026-10-28T03:00:00Z`（已过期） |
| 保留金 | `10.0 USD` |
| 本地累计支出 | `4.113678 USD`（不是余额） |
| 判断 | 余额无法核实 → `WAVE11B1_SCIENTIFIC_REVIEW` 保持 `BLOCKED_EXTERNAL` |
| 影响范围 | 仅阻断 B1 与依赖它的 E；B0R 等无关任务照常执行 |

`loggedIn: false` 单独不作为唯一判据（若能返回可核实的 `total_cost_usd` 仍可记录），但本次余额时间戳过期且无可用登录态，因此不发起付费调用。

## 未使用的能力

- 未使用 `--force`。Grok 走 `--trust --auto-review`，smoke 验证可写文件（`/tmp/grok_smoke/b.txt` = `SMOKE_OK`，exit 0）。
- 未启用 Antigravity/Gemini。
- 未启用 Pi / Pi Durable。
