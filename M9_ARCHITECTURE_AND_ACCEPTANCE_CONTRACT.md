# Milestone M9: Architecture and Acceptance Contract
## 报告全面审阅与发布门禁阶段契约 (Report Audit & Release Gate Contract)

- **Milestone ID**: `M9_REPORT_AUDIT_AND_RELEASE_GATE`
- **Controller**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **Reviewer**: OpenAI Codex CLI / `gpt-6.1-sol`
- **Developer P1**: OpenCode CLI / GLM-5.3
- **Operating Policy**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **Strategy Document**: [`docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md`](docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md)
- **Status**: `M9_CONTRACT_FORMALIZED`
- **Date**: 2026-11-04

---

## 1. 阶段目标与适用范围 (Objectives & Scope)

### 1.1 核心目标
1. **全量报告证据链双向可追溯性审计 (Bidirectional Traceability Audit)**:
   - 系统核对 ProkaDiff 产品报告（Markdown 与 TSV）中的每一项 Finding 是否 100% 具备双向回指证据：
     `Evidence (GD / BAM) ↔ Differential Variant (DV1 EventId) ↔ Intended Evaluation ↔ Mechanistic Association ↔ AuditResult / Final Report`.
   - 杜绝报告中出现任何凭空捏造、缺少底表证据或无法回溯原始测序支持的孤立断言。
2. **科学措辞与因果越界防护审计 (Scientific Language & Causality Guarding)**:
   - 依据 [`docs/scientific_model.md`](docs/scientific_model.md)，严格执行“观察优先、因果最后”准则。
   - 杜绝将“空间/几何关联”夸大为“Cas 切割因果”，统一使用 `candidate guide-dependent off-target`，严禁使用 `confirmed off-target`。
   - 显式声明成对 WGS 分析的局限性：*Spatial association is consistent with a possible guide-dependent event but does not establish cleavage causality.*
3. **跨里程碑能力与已知局限性权威注册 (Capabilities & Limitations Registry)**:
   - 整合 M5 (BL21 真实克隆)、M6 (PRJNA1088182 队列)、M7 (PRJNA884016 组装真值)、M8 (REL606 单样本对拍与性能) 的实际测定结论与边界。
   - 建立权威机器可读与人类可读的已验证能力台账，将未经正交生物学验证、高拷贝重复解析受限等已知科学局限如实固定。
4. **防御性发布门禁程序与自动化检查 (Automated Release Gate)**:
   - 编写确定性发布门禁检查工具（`scripts/audit_release_gates.py`），对生成的报告、TSV 及文本进行静态与动态完整性校验，违背科学契约时严格 fail-closed。

---

## 2. 核心技术不变量与科学红线 (Invariants & Red Lines)

1. **证据链不可断裂原则 (M4 Invariant)**:
   - 所有产品报告必须且仅能从统一的 `AuditResult` 实例派生。
   - `DV1 EventId` 为公开差分变异的唯一权威身份标识；任何衍生报告行必须包含或可直接追溯至有效 `EventId`。
2. **科学因果红线 (Rule 7, Rule 8 & Scientific Model)**:
   - 严禁在缺乏体外切割或独立正交拯救实验的情况下，宣称某突变为 Cas9/Cas12a 切割所致的直接突变。
   - 严禁在未独立解析重复拷贝时发出拷贝特异性结构调用。
   - 严禁擅自修改已冻结的 G1 结构变异现状（保持 `FROZEN_HOLD`，`PASS=false`）。
3. **环境与执行红线 (Rule 11)**:
   - 登录节点仅用于轻量元数据管理与调度；所有自动化发布门禁检查与测试必须在 Slurm 分区（`qfree` 或 `qcpu_18i`）执行。
4. **Git 与工作树安全红线**:
   - 严格单写入者；严禁自动 push、merge、tag 或 release；最终发布与合并权限属于人类用户。

---

## 3. 验收门禁矩阵 (Mandatory Acceptance Gates)

| 门禁 ID | 门禁名称 | 验收标准 | 验证工具 / 证据 |
| --- | --- | --- | --- |
| **Gate A** | **代码与回归测试套件全绿** | 全工作区 Cargo 测试 441 项全过 (3 ignored)；Python 测试 264+ 项全绿。 | `cargo test --workspace`, `pytest` on `bnode2` |
| **Gate B** | **双向证据链回指审计** | M5–M8 产品报告中所有 finding 100% 回指有效 `EventId` 与 GD 记录，无悬空断言。 | `test_audit_report_fixtures.rs`, 审计验证脚本 |
| **Gate C** | **科学措辞与因果越界防护** | 全文扫描 0 处违规因果词汇（如 `confirmed off-target`）；空间关联均附带因果限定声明。 | 确定性文本正则与 AST 审计扫描器 |
| **Gate D** | **权威能力与局限总台账** | 形成 `CAPABILITIES_AND_LIMITATIONS_REGISTRY.md`，如实收口 M5–M8 真实边界。 | 审查归档文档 |
| **Gate E** | **自动化防御性发布门禁** | 编写 `scripts/audit_release_gates.py`，全量测试用例验证通过，具备 fail-closed 防护力。 | `pytest scripts/test_audit_release_gates.py` |
| **Gate F** | **独立同行审查与发布签收** | Codex Sol 独立审查通过，形成完整的 M9 交付与发布准备报告。 | `M9_FINAL_INTEGRATION_REPORT.md` |

---

## 4. 交付件清单 (Deliverables)

1. `docs/CAPABILITIES_AND_LIMITATIONS_REGISTRY.md`: M5–M8 全面能力、边界与已知限制权威台账。
2. `scripts/audit_release_gates.py`: 报告证据链与措辞防御性自动化发布门禁工具。
3. `scripts/test_audit_release_gates.py`: 发布门禁工具的正反向测试用例（含因果越界拦截负向用例）。
4. `M9_FINAL_INTEGRATION_REPORT.md`: Milestone M9 最终集成审计与发布准备报告。
5. 完整的 M9 台账系统（`STAGE_ROADMAP.md`, `M9_TASK_LEDGER.md`, `M9_WORKER_DISPATCH_LEDGER.md`, `HANDOFF.md`, `ACTIVE_TASK.md`）。
