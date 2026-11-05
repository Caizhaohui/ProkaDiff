# Milestone M8: Worker Dispatch Ledger

> Active Governing Policy: **Multi-AI Collaboration V3.2** (`~/.config/ai-collaboration/POLICY.md`)  
> Operating Strategy: [`docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md`](docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md)  
> Active Controller: **Antigravity CLI / Gemini 3.8 Flash Medium** (`AUTONOMOUS_EXECUTION_CONTROLLER`)  
> Updated: 2026-11-03T17:00:00+08:00  

---

## 外部 Worker 调度台账 (Worker Dispatch Ledger)

| 调度编号 (DISP_ID) | 关联任务 (TASK_ID) | 调度角色与模型 | CLI 入口与工作目录 | 授权范围 / 变更集 | 执行结果 (Exit Code) | 产出物 / 耗费 | 审核与验收状态 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **`DISP-M8-001`** | `M8-T01` | Controller (Gemini 3.8 Flash) | Antigravity CLI | 全项目环境与锁清理 | 0 (SUCCESS) | 清除旧锁，停用 supervisor，重构规范 | Controller 自核验证通过 |
| **`M8-DISP-002`** | `M8-T02` | Controller / Slurm bnode2 | `cargo test --workspace`, `pytest` | 单元与回归测试基准审计 | 0 / 1 (4 test diff checks) | 307/307 Rust tests pass; 140/144 Python tests pass | Controller 审计通过并识别脆弱断言 |
| **`M8-DISP-003`** | `M8-T03` | Developer P1 (GLM-5.3) / Reviewer (Codex Sol) | `opencode run`, `codex exec` | `scripts/parity/audit_rel606_non_jc.py`, `test_audit_rel606_non_jc.py` | 0 (SUCCESS) | 39/39 严格一致 (0 over, 0 under, delta=0); 19/19 tests pass | Codex Sol 独立审查 PASS |
| **`M8-DISP-004`** | `M8-T04` | Developer P1 (GLM-5.3) / Reviewer (Codex Sol) | `opencode run`, `codex exec` | `benchmark/real_world/scripts/test_m6_*.py` | 0 (SUCCESS) | `git status --porcelain` 守护；163/163 Python 测试全绿 | Codex Sol 独立审查 PASS |
| **`M8-DISP-005`** | `M8-T05` | Developer P1 (GLM-5.3) / Reviewer (Codex Sol) | `opencode run`, `codex exec` | `scripts/parity/generate_m8_perf_summary.py`, `test_m8_perf_summary.py` | 0 (SUCCESS) | `m8_perf_summary.tsv` 提取真实测得 5.29x 加速与内存比；34/34 tests pass | Codex Sol 独立审查 PASS |
| **`M8-DISP-006`** | `M8-T06` | Controller (Antigravity) / Reviewer (Codex Sol) | Slurm `bnode2` / `codex exec` | `benchmark/results/parity_closure/g1_clonal_jc/review/WAVE11E_HOLD.md` | 0 (SUCCESS) | G1 严格门重算 MATCHED=10, OVER=24, UNDER=8, PASS=false 完全一致，维持受控固化 | 独立审查确认重算一致，指出生产源码残留距离函数需如实记入限制 |
| **`M8-DISP-007`** | `M8-T07` | Reviewer (Codex Sol `gpt-6.1-sol`) | `codex exec -m gpt-6.1-sol -s read-only` | M8 全部交付件、测试、台账与审计代码 (Round 1) | 0 (SUCCESS) | 独立审查报告：Gate C/D PASS；Gate E/F FAIL (指出 Gate E 必须诚实报告为受控未通过，禁止粉饰 PASS，对齐测试计数，补全台账闭环) | 审查完成，提出整改意见并由 Controller 落实修复 |
| **`M8-DISP-008`** | `M8-T07` | Reviewer (Codex Sol `gpt-6.1-sol`) | `codex exec -m gpt-6.1-sol -s read-only` | M8 全部交付件与审计代码 (Round 2) | 0 (SUCCESS) | 独立审查报告：确认 Gate E 披露正确；指出 Cargo 准确计数为 441 passed / 3 ignored，发现报告中存在重复段落与过早记录 PASS | 审查完成，Controller 已完成清理与精确对齐 |
| **`M8-DISP-009`** | `M8-T07` | Reviewer (Codex Sol `gpt-6.1-sol`) | `codex exec -m gpt-6.1-sol -s read-only` | M8 全部交付件与审计代码 (Round 4 Final) | 0 (SUCCESS) | 独立审查最终裁决：**PASS WITH RECORDED LIMITATIONS**。确认测试计数精准、台账闭环一致、G1 限制保留明晰 | 最终独立签收通过，Milestone M8 正式闭环收口 |
