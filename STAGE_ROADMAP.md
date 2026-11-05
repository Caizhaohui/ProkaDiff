# ProkaDiff Stage Roadmap

> Active Governing Policy: **Multi-AI Collaboration V3.2** (`~/.config/ai-collaboration/POLICY.md`)  
> Operating Strategy: [`docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md`](docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md)  
> Active Controller: **Antigravity CLI / Gemini 3.8 Flash Medium** (`AUTONOMOUS_EXECUTION_CONTROLLER`)  
> Updated: 2026-11-03T16:59:00+08:00  

---

## 1. 里程碑全局路线图 (Milestone Roadmap Overview)

| 里程碑 | 核心目标 | 状态 | 验证与科学结论 |
| --- | --- | --- | --- |
| **M0** | 验证声明与可复核输入收口 | **COMPLETE** (2026-09-18) | Manifest、truth 与 SRA 元数据版本化纳管；缺失结果退回 PENDING。 |
| **M1** | Evidence → Differential Variant | **COMPLETE** (2026-09-22) | `DV1 EventId` 确定性生成；MCMF 容差亲本扣除；全工作区测试与 Slurm 验证全过。 |
| **M2** | Differential Variant → Intended Verification | **COMPLETE** (2026-09-23) | 逐编辑评估与异常结构判定；生物拟合度决胜；Layer 1 端到端双株 FASTQ 差分验证通过。 |
| **M3** | Intended Verification → Mechanistic Association | **COMPLETE** (2026-09-24) | 单次候选位点扫描；类型化错误隔离；空间/序列机制关联与 Cas 切割因果解耦。 |
| **M4** | Mechanistic Association → AuditResult 与全部产品输出 | **COMPLETE** (2026-09-24) | V2 默认追加列、V1 字节兼容；单一 AuditResult 聚合与原子化发布/回滚测试通过。 |
| **M5** | BL21 真实重测序数据复核 | **COMPLETE** (2026-10-07) | SR1 重复接头恢复、BL21 单样本参考复核及 6 格 PEER_COMPARATOR 对照完成。 |
| **M6** | PRJNA1088182 队列复核 | **CLOSED_WITH_LIMITATIONS** (2026-10-23) | Tier A 发表突变共识检出 2/12（其余亚克隆 ALT fraction 0.20–0.69）；无正交生物学金标准；登记行保持 PENDING。 |
| **M7** | PRJNA884016 组装真值复核 | **COMPLETE_WITH_LIMITATIONS** (2026-10-28) | 染色体组装 100% 比对相同；历史 4 MOB + 7 JC 判定为局部截断假象并隔离；阳性 SV 灵敏度设为 NOT_EVALUABLE；登记行保持 PENDING。 |
| **M8** | **REL606 引擎扩展对拍与性能基准** | **CLOSED_WITH_LIMITATIONS** (2026-11-04) | 39 个非 JC 突变 100% 一致 (0 over, 0 under)；独立测定 5.29x 加速；G1 严格门保持 FROZEN_HOLD (PASS=false)；Codex Sol 签收通过。 |
| **M9** | **报告审阅与发布门禁** | **COMPLETE** (2026-11-04) | 确定性发布门禁工具与 195 项测试全绿；全库 0 违规因果词汇；M5–M8 真实报告与 TSV 100% 具备有效 DV1 EventId 回指；Codex Sol 独立签收通过。 |
| **M10** | **多克隆工程批次审计模式** | **COMPLETE** (2026-11-04) | 依据 `docs/batch_mode_spec.md` 实现 batch manifest 输入、一祖先多克隆审计聚合、线协议变异复现性分类、Taskbook §66 非处方性批次报告与 TSV 矩阵；Codex Sol 独立签收通过。 |

---

## 2. 历史阶段闭环状态
- **M0–M7**: 均已收口（M6 CLOSED_WITH_LIMITATIONS, M7 COMPLETE_WITH_LIMITATIONS）。
- **M8**: REL606 引擎对拍与性能基准已收口（CLOSED_WITH_LIMITATIONS, 39 个核心突变 100% 一致，5.29x 加速，G1 严格门受控固化保持 PASS=false）。
- **M9**: 报告全面审阅与发布门禁已收口（COMPLETE, 195 门禁测试/392 Python 测试/441 Rust 测试全绿，Codex Sol 独立审查签收）。

---

## 3. 已完成阶段 (Milestone M9: 报告审阅与发布门禁) 任务矩阵

| 任务 ID | 工作包内容 | 依赖关系 | 执行角色 | 验证工具 / 门禁 | 最终状态 |
| --- | --- | --- | --- | --- | --- |
| **`M9-T01`** | **M9 契约制定与治理矩阵初始化** | M8 闭环 | Controller (Antigravity) | `M9_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md` | **COMPLETE** |
| **`M9-T02`** | **真实数据集报告全量证据链审计** | `M9-T01` | Codex Sol / Controller | `scripts/audit_release_gates.py` Gate 2 (PASS, 0 findings) | **COMPLETE** |
| **`M9-T03`** | **科学措辞与因果越界防护审计** | `M9-T01` | OpenCode (P1) / Reviewer | `scripts/audit_release_gates.py` Gate 1 (PASS, 0 findings) | **COMPLETE** |
| **`M9-T04`** | **跨里程碑能力与限制总台账构建** | `M9-T02`, `M9-T03` | Controller / OpenCode (P1) | `docs/CAPABILITIES_AND_LIMITATIONS_REGISTRY.md` | **COMPLETE** |
| **`M9-T05`** | **自动化防御性发布门禁工具开发** | `M9-T02`, `M9-T03` | OpenCode GLM-5.3 (P1) | `scripts/audit_release_gates.py`, `test_audit_release_gates.py` (195 passed in 41.40s) | **COMPLETE** |
| **`M9-T06`** | **M9 阶段集成与独立审查验收** | `M9-T04`, `M9-T05` | Codex Sol (Reviewer) / Controller | `M9_FINAL_INTEGRATION_REPORT.md` (Codex Sol 签收通过) | **COMPLETE (PASS)** |

---

## 4. 当前阶段 (Milestone M10: 多克隆工程批次审计模式) 工作包矩阵

依据 `docs/batch_mode_spec.md`，将批次审计功能分解为如下独立工作包：

| 任务 ID | 工作包内容 | 依赖关系 | 推荐执行角色 | 验证工具 / 门禁 | 当前状态 |
| --- | --- | --- | --- | --- | --- |
| **`M10-T01`** | **M10 契约制定与多克隆数据模型确立**：依据 `docs/batch_mode_spec.md` 明确批次清单格式、多克隆审计聚合模型 `BatchAuditResult`、复现性分类 `RecurrenceStatus` 及非推荐性原则（Taskbook §66）。 | M9 闭环 | Controller (Antigravity) | `M10_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`, `M10_TASK_LEDGER.md`, `M10_WORKER_DISPATCH_LEDGER.md` | **READY** |
| **`M10-T02`** | **Batch Manifest 解析器与 CLI 参数拓展**：实现 `--batch-manifest <path>` 解析与校验，检查克隆 ID 唯一性、FASTQ 路径存在性、可选每克隆 spacer/intended 覆盖，缺少 starter 或 ref 报错 fail-closed。 | `M10-T01` | Developer P1 (OpenCode GLM-5.3) | 单元测试 + CLI 语法与错误拦截测试 | **QUEUED** |
| **`M10-T03`** | **多克隆变异对齐与复现性分类器 (Recurrence Engine)**：在保持各克隆内部 `DV1 EventId` 确定性的前提下，跨克隆对齐差分事件，将变异精确分类为 `CloneSpecific`、`SharedSubset`、`SharedAll`。 | `M10-T01` | Developer P1 (OpenCode GLM-5.3) | 变异对齐算法单元测试 + 守恒性校验 | **QUEUED** |
| **`M10-T04`** | **批次报告生成器 (`batch_report.md` & `batch_variants.tsv`)**：生成执行摘要、克隆摘要矩阵、共享/特异变异表、结构变异概览与局限说明，严禁自动评选“最佳克隆”。 | `M10-T02`, `M10-T03` | Developer P1 (OpenCode GLM-5.3) | 黄金报告夹具测试 + Release Gate 语法审计 | **QUEUED** |
| **`M10-T05`** | **多克隆批次模式端到端集成测试**：构建包含 3+ 克隆的合成批次夹具，验证从 `--batch-manifest` 到最终批次报告与各克隆独立报告的完整一致性。 | `M10-T04` | Developer P1 (OpenCode GLM-5.3) / Reviewer | `cargo test --workspace`, 集成测试通过 | **QUEUED** |
| **`M10-T06`** | **M10 阶段集成与独立审查验收**：运行全量 Cargo / Python 测试套件与发布门禁，由 Codex Sol 独立签收并完成阶段闭环交接。 | `M10-T05` | Codex Sol (Reviewer) / Controller | `M10_FINAL_INTEGRATION_REPORT.md`, `HANDOFF.md` | **QUEUED** |

---

## 3. 执行规则与状态推进原则

1. **单任务顺序实施**: 严格按照 `READY` 任务逐个推进，优先调用 **OpenCode GLM-5.3** 进行实施。
2. **三轮修复上限**: 针对具体工程缺陷，同一 Worker 连续修复不超过 3 轮；若未通过测试门禁，自动上报 Controller 并考虑模型升级（Grok 或 Astra/Opus）。
3. **独立审计与验证**: 任何重要代码改动完成后，实施者不得自行判定通过，必须由 **Codex Sol/Terra** 独立审查，并在 Slurm 计算节点上由 **Codex Luna + 确定性脚本** 完成验证。
4. **G1 严格红线**: 保持严格门不变，禁止在未获用户授权前修改底层生产归因或提交新的未授权 G1 Slurm 作业。
