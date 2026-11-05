# Milestone M8: Task Ledger

> Active Governing Policy: **Multi-AI Collaboration V3.2** (`~/.config/ai-collaboration/POLICY.md`)  
> Operating Strategy: [`docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md`](docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md)  
> Active Controller: **Antigravity CLI / Gemini 3.8 Flash Medium** (`AUTONOMOUS_EXECUTION_CONTROLLER`)  
> Updated: 2026-11-03T17:00:00+08:00  

---

## 任务执行台账 (Task Ledger)

| 任务编号 | 任务名称 | 责任角色 / 分配模型 | 前置依赖 | 验收条件 | 产出物 / 证据路径 | 当前状态 |
| --- | --- | --- | --- | --- | --- | --- |
| **`M8-T01`** | **治理体系核验与环境清理** | Controller (Antigravity Gemini 3.8) | 无 | 确认唯一主控；无未受控后台守护进程；模型与预算就绪；Slurm 正常。 | `benchmark/results/parity_closure/state/supervisor_control.json` (DEACTIVATED), lock archive | **COMPLETE** |
| **`M8-T02`** | **基线代码与测试套件状态审计** | Controller (Antigravity Gemini 3.8) | `M8-T01` | Cargo 测试全过；排查 Python 4 项失败根因；确认 G1 严格门状态。 | `cargo test --workspace` 307/307 pass; 4 failures identified as working tree diff checks | **COMPLETE** |
| **`M8-T03`** | **REL606 非 JC 突变全量对拍审计与报告** | Developer P1 (OpenCode GLM-5.3) | `M8-T02` | 39 个非 JC 突变 (SNP/INS/DEL/MOB) 100% 严格一致；0 over, 0 under；版本规范化；Kuhn 迭代式增广路最大二分匹配；19/19 单元测试全过；Codex Sol 独立审查 PASS。 | `benchmark/results/m8_rel606_non_jc_parity.tsv`, `scripts/parity/audit_rel606_non_jc.py`, `scripts/parity/test_audit_rel606_non_jc.py` | **COMPLETE** |
| **`M8-T04`** | **M6 回归测试套件工作树隔离对齐** | Developer P1 (OpenCode GLM-5.3) | `M8-T02` | 优化 `test_m6_*.py` 中对 `crates/` 的脆弱检查，使用 `git status --porcelain` 严格守护 `VALIDATION_REGISTRY.tsv`；`pytest` 163/163 全绿；Codex Sol 审查 PASS。 | `benchmark/real_world/scripts/test_m6_*.py` | **COMPLETE** |
| **`M8-T05`** | **ProkaDiff vs breseq 独立性能基准脚本与测量** | Developer P1 (OpenCode GLM-5.3) / Controller | `M8-T01` | 解析 Slurm 2769249 `/usr/bin/time -v` 日志生成 `m8_perf_summary.tsv`，精准记录 Wall-clock、CPU% 与 Peak RSS，证明 5.29x 加速；34/34 测试全过；Codex Sol 审查 PASS。 | `benchmark/results/m8_perf_summary.tsv`, `scripts/parity/generate_m8_perf_summary.py` | **COMPLETE** |
| **`M8-T06`** | **G1 结构变异现状冻结与限制陈述** | Controller / Codex Sol | `M8-T02` | 保持严格门；确认无未授权生产篡改；如实保留 OVER=24, UNDER=8 限制说明（状态保持 PASS=false / FROZEN_HOLD）。 | `benchmark/results/parity_closure/g1_clonal_jc/review/WAVE11E_HOLD.md` | **COMPLETE** (严格受控固化，PASS=false) |
| **`M8-T07`** | **M8 阶段集成与独立审查验收** | Codex Sol (Reviewer) / Controller | `M8-T03`, `M8-T04`, `M8-T05`, `M8-T06` | 全门禁 (Gate A–Gate F) 审计收口；Codex Sol 审查反馈落实；Gate E 诚实记录为未通过/受控阻断；Codex Sol Round 4 签收。 | `M8_FINAL_INTEGRATION_REPORT.md`, `docs/development/HANDOFF.md` | **COMPLETE (PASS WITH RECORDED LIMITATIONS)** |
