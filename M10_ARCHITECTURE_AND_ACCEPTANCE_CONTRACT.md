# Milestone M10: Architecture and Acceptance Contract
## 多克隆工程批次审计模式架构与验收契约 (Multi-Clone Engineering Batch Audit Mode)

- **Milestone ID**: `M10`
- **Controller**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **Governing Policy**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **Specification Source**: [`docs/batch_mode_spec.md`](docs/batch_mode_spec.md) & Taskbook §64–67
- **Activation Date**: 2026-11-04
- **Status**: `ACTIVE`

---

## 1. 阶段目标与业务场景 (Objective & Scope)

在微生物基因组工程实践中，实验人员通常会从同一批转化或突变实验中挑选多个克隆进行并行测序分析。
当前 ProkaDiff 的标准操作为成对审计（一对一）：
$$\text{Starter Strain WGS} \longleftrightarrow \text{Edited Strain WGS}$$

Milestone M10 扩展并实现**一祖先对多克隆（One-starter vs Many-clones）工程批次模式**：
```text
                     Starter Strain WGS
                             │
     ┌───────────────┬───────┴───────┬───────────────┐
     ▼               ▼               ▼               ▼
 Edited Clone 01  Edited Clone 02  Edited Clone 03  Edited Clone N
```

核心交付目标：
1. **统一批次清单驱动**：通过 `--batch-manifest <tsv>` 批量组织与校验各克隆测序数据（FASTQ、每克隆可选 intended/spacer 覆盖）。
2. **多克隆并行/顺序审计聚合**：对清单中各克隆独立执行标准化差分审计，产出各克隆独立审计子目录。
3. **跨克隆变异对齐与复现性分类 (Recurrence Engine)**：对齐各克隆差分事件，将突变分类为 `CloneSpecific`（克隆特异）、`SharedSubset`（部分克隆共享）、`SharedAll`（全克隆共享）。
4. **批次审计综合报告与数据表**：生成 `batch_report.md` 与 `batch_variants.tsv`，汇总全批次编辑状态与共性变异。

---

## 2. 核心架构设计与数据模型 (Domain Models)

### 2.1 批次清单规范 (`batch_manifest.tsv`)
- 列定义：
  - `clone_id` (必填): 唯一克隆字符串标识符（如 `C01`, `Clone_A`）。
  - `r1` (必填): R1 FASTQ 路径（支持 `.fq`, `.fastq`, `.gz`）。
  - `r2` (选填): R2 FASTQ 路径（双端测序填入，单端留空）。
  - `intended_tsv` (选填): 每克隆可选预期编辑表覆盖（未填则沿用全局 `--intended`）。
  - `spacer` (选填): 每克隆可选 spacer 序列覆盖。
  - `pam` (选填): 每克隆可选 PAM 序列覆盖。

### 2.2 复现性分类模型 (`RecurrenceStatus`)
```rust
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum RecurrenceStatus {
    CloneSpecific,
    SharedSubset { count: usize, total: usize },
    SharedAll { total: usize },
}
```

### 2.3 批次变异记录 (`BatchVariantRecord`)
```rust
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BatchVariantRecord {
    pub batch_variant_id: String,
    pub seq_id: String,
    pub position: u64,
    pub end: u64,
    pub gd_type: String,
    pub recurrence: RecurrenceStatus,
    pub clone_presence: HashMap<String, Option<String>>, // clone_id -> Option<EventId>
}
```

### 2.4 批次聚合结果 (`BatchAuditResult`)
```rust
pub struct BatchAuditResult {
    pub manifest_path: PathBuf,
    pub starter_info: SampleInfo,
    pub clone_entries: Vec<CloneAuditEntry>,
    pub variants_matrix: Vec<BatchVariantRecord>,
}
```

---

## 3. 严格科学与安全契约 (Binding Invariants)

1. **严格非推荐性原则 (Strict Non-Prescriptive Policy, Taskbook §66)**:
   - ProkaDiff 是客观科学审计工具，**严禁**在批次报告或日志中宣称任何克隆为“最佳 (best)”、“最安全 (safest)”或“推荐 (recommended)”。
   - 克隆优选决策取决于生物学表型、表达滴度及工程下游目标，属于用户科研裁决范畴。
2. **EventId 与单克隆审计真实性不变量**:
   - 每个克隆内部保持完整的独立 `AuditResult` 与 `DV1 EventId` 生成机制，绝不伪造全局虚拟事件覆盖各克隆真实变异。
   - 各克隆标准产品集（`report.md`, `post_edit_variants.tsv`, `unintended.tsv`, `edit_outcomes.tsv`, `summary.txt`）完整生成在其独立子目录下。
3. **出发株 WGS 强制性 (Rule 6)**:
   - 批次模式依然必须提供公共 `--starter` 与 `--ref`，缺少任一参数直接退出并报错。
4. **发布门禁与科学用语 100% 遵从**:
   - `batch_report.md` 必须通过 `scripts/audit_release_gates.py` 全部 3 个发布门禁，绝不允许引入未经限定的因果词汇或虚假断言。

---

## 4. 阶段验收门禁 (Quality Gates)

| 验收门禁 | 门禁要求 | 验证方式 |
|---|---|---|
| **Gate A: 数据模型与 Manifest 解析** | 准确解析合法 TSV manifest，严格拦截克隆 ID 重复、FASTQ 文件不存在、列数不齐等非法输入。 | 纯内存单元测试 |
| **Gate B: 跨克隆变异对齐与分类** | 确定性对齐多克隆变异；验证 `CloneSpecific`, `SharedSubset`, `SharedAll` 状态分类准确，支持不同克隆同一突变回指其独立 `EventId`。 | 核心算法单元测试 + 守恒性测试 |
| **Gate C: 批次报告确定性输出** | 输出结构完整的 `batch_report.md` 与 `batch_variants.tsv`；包含执行摘要、克隆指标表、共性突变分析；0 违规推荐用语。 | 黄金报告渲染测试 |
| **Gate D: 端到端批次集成测试** | 针对包含多个克隆的合成批次数据集，CLI 成功运行并输出完整批次与分克隆产品目录。 | 集成测试 (`tests/batch_e2e.rs`) |
| **Gate E: 全工作区回归与发布门禁** | `cargo test --workspace` 全绿；`audit_release_gates.py` 3 个门禁 0 findings 全绿。 | 计算节点 `bnode2` 实测 |
| **Gate F: 独立评审与签收** | Codex Sol 审查通过并出具正式签收结论。 | `codex exec -m gpt-6.1-sol` 审查 |
