# Milestone M11: Task Ledger
## 亚克隆低频变异与混群审计模式任务台账

> Active Governing Policy: **Multi-AI Collaboration V3.2** (`~/.config/ai-collaboration/POLICY.md`)  
> Operating Strategy: [`docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md`](docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md)  
> Active Controller: **OpenCode CLI / GLM-5.3-Flash** (acting Controller per user-authorized handoff; former Controller Antigravity / Gemini 3.8 Flash Medium exhausted weekly quota)  
> Updated: 2026-11-05T11:00:00+08:00  

---

## 任务执行台账 (Task Ledger)

| 任务编号 | 任务名称 | 责任角色 / 分配模型 | 前置依赖 | 验收条件 | 产出物 / 证据路径 | 当前状态 |
| --- | --- | --- | --- | --- | --- | --- |
| **`M11-T01`** | **M11 科学契约制定与低频统计模型设计** | Scientific Architect (Claude Opus 5.5) / Controller | M10 闭环 | 确立低频变异二项分布显著性检验、Wilson 置信区间、链偏倚 Fisher Exact 检验、背景错误率估计与 CLI 参数契约；初始化台账。 | `M11_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`, `M11_TASK_LEDGER.md`, `M11_WORKER_DISPATCH_LEDGER.md` | **COMPLETE** |
| **`M11-T02`** | **低频 Pileup 统计与质量过滤引擎实现** | Developer P1 (OpenCode GLM-5.3) | `M11-T01` | 实现链特异读段解析、本地错误率背景估计、Fisher 链偏倚计算与 Wilson CI；单元测试全绿。 | `crates/prokadiff-evidence/src/subclonal.rs` | **COMPLETE** |
| **`M11-T03`** | **SubclonalVariant 领域模型与分类器聚合集成** | Developer P1 (OpenCode GLM-5.3) | `M11-T02` | 将 `SubclonalMetrics` 挂载至变异模型，支持默认共识模式与 opt-in 低频模式双模解析，实现亲本背景扣除；测试全绿。 | `crates/prokadiff-classify/src/subclonal.rs` | **COMPLETE** |
| **`M11-T04`** | **非处方性报告与 TSV 输出扩展** | Developer P1 (OpenCode GLM-5.3) | `M11-T03` | `post_edit_variants.tsv` 新增频率与 CI 列；`report.md` 渲染亚克隆小节，严格附带非处方性与不确定度警示；发布门禁通过。 | `crates/prokadiff-report/src/subclonal.rs`, `tables.rs` | **COMPLETE** |
| **`M11-T05`** | **合成梯度稀释测试套件 (Titration Test Suite)** | Developer P1 (OpenCode GLM-5.3) | `M11-T04` | 构建 5%, 10%, 25%, 50% 梯度稀释合成夹具；验证检出敏感度、假阳性控制及默认模式无回归。 | `crates/prokadiff/tests/subclonal_titration.rs` | **COMPLETE** |
| **`M11-T06`** | **M11 阶段集成与独立审查验收** | Reviewer (Codex Sol `gpt-6.1-sol`) / Controller | `M11-T05` | Gate A–Gate E 全门禁通过；`cargo test --workspace` 全绿；发布门禁 CLI 0 findings；Codex Sol 给出正式签收裁决；输出最终交接。 | `M11_FINAL_INTEGRATION_REPORT.md`, `HANDOFF.md` | **IN_PROGRESS** — R1 REJECTED (9 findings, 09:43) → 修复 → R2 REJECTED (4×P1+2×P2, 10:19, `.prokadiff/codex_m11_r2_review.txt`) → 6 项 R2 修复已落地 (10:20–10:27, 未复审) → Controller 接管恢复审计：全工作区 481 passed / 0 failed / 3 ignored (bnode2, alloc 3095377, `.prokadiff/m11_recovery_full_workspace_test.log`)；待 Codex Sol Round 3 复审当前快照 |
