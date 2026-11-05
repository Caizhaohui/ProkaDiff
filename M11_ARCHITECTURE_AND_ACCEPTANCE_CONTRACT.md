# Milestone M11: Architecture and Acceptance Contract
## 亚克隆低频变异与混群审计模式架构与验收契约 (Subclonal / Pool-seq Variant Audit Mode)

- **Milestone ID**: `M11`
- **Authorized by**: User (2026-11-05)
- **Governing Policy**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **Controller**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **Scientific Architect**: Claude Opus 5.5 (`claude-opus-5-5`)
- **Developer P1**: OpenCode CLI / GLM-5.3
- **Independent Reviewer**: OpenAI Codex CLI / `gpt-6.1-sol`
- **Status**: `ACTIVATED`

---

## 1. 核心目标与科学边界 (Scope & Invariants)

### 1.1 核心问题
在原核生物基因编辑（如 CRISPR/Cas9/Cas12a 或混群富集实验）中，样本常呈现菌群异质性或不完全敲除，存在等位基因频率低于 0.75 共识阈值的真实生物学亚克隆变异。M11 旨在为 ProkaDiff 引入**亚克隆低频变异与混群审计模式**，回答：
> **“在共识固定变异之外，该工程菌群中是否存在统计显著、高于测序背景噪声的低频编辑或脱靶事件？”**

### 1.2 严格科学约束与铁律
1. **默认共识调用完全不变（Opt-in Only）**：
   - 默认模式下，ProkaDiff 维持现行 breseq 对拍与高纯度共识变异调用行为（$\ge 0.75$）；
   - 亚克隆审计必须通过显式命令行参数开启（如 `--subclonal` 或 `--min-subclonal-af <FLOAT>`，默认关闭）；
   - 默认模式下的所有现有 460+ 项测试与输出产物**逐字节不发生任何回归改变**。
2. **变异类型严格受限**：
   - 本期**仅授权**：单核苷酸变异（SNP）、替换（SUB）、小型插入与删除（$\le 2\text{ bp}$ 的 INS/DEL）；
   - **严格排除**：移动元件（MOB）、接合点（JC）、基因重复扩增（AMP）及复杂结构变异的低频调用（此类变异的假阳性与重复区比对歧义极高，不在本期授权范围内）。
3. **目标研究范围 vs 灵敏度声明**：
   - 目标分析区间为 ALT fraction 0.05–0.70；
   - **严禁**直接将 0.05–0.70 宣称为“已验证的生物学灵敏度区间”；报告中必须明确说明低频检测灵敏度受测序深度（Depth）、背景测序错误率（Error Rate）、链偏倚（Strand Bias）及比对歧义（Mapping Quality）的物理制约。
4. **统计模型与置信度要求**：
   - 必须基于链特异支持读段（Forward/Reverse ALT count）计算链偏倚检验（如 Fisher Exact Test $p$-value 或 Strand Odds Ratio）；
   - 必须提供双侧二项分布（Binomial / Wilson Score）置信区间（如 95% CI）；
   - 必须设定最小绝对支持读段数（如 ALT reads $\ge 3$）与位点总深度门槛（如 Depth $\ge 20\times$）。
5. **历史科学结论与严格门固化**：
   - G1 Clonal JC 严格门受控固化保持 `MATCHED=10, OVER=24, UNDER=8, PASS=false`；
   - PRJNA1088182 与 PRJNA884016 登记行保持 `PENDING`；
   - 严禁在未经授权前发起大规模计算节点真实 WGS 重跑作业。

---

## 2. 统计模型数学推导与契约规范 (Statistical Formulations)

### 2.1 连续校正 Wilson Score 置信区间 (Wilson Score with Continuity Correction)
样本等位基因比例 $\hat{p} = A / N$，其中 $A$ 为候选突变等位基因（ALT）支持读段数，$N = R + A + O$ 为有效位点深度。
显著性水平 $\alpha = 0.05$（双侧 $z = 1.959964$）。
为消除极小样本/极低频率下的覆盖偏差，规范采用连续校正（Newcombe 1998, Method 4）：
$$L = \max\left(0.0, \frac{2A + z^2 - 1 - z \sqrt{z^2 - 2 - \frac{1}{N} + 4 \frac{A}{N}(N - A + 1)}}{2(N + z^2)}\right)$$
$$U = \min\left(1.0, \frac{2A + z^2 + 1 + z \sqrt{z^2 + 2 - \frac{1}{N} + 4 \frac{A}{N}(N - A - 1)}}{2(N + z^2)}\right)$$
- 边界约束：$A = 0 \implies L = 0$；$A = N \implies U = 1$；若根号内因浮点误差小于 0 则置为 0；$N = 0$ 返回 `Err(StatsError::ZeroDepth)`。
- 测试验证基准向量：$A = 0, N = 20 \implies U \approx 0.2005$（标准 Wilson $U = 0.16113$）；CC 区间对任意 $(A, N)$ 恒包含标准 Wilson 区间。

### 2.2 链偏倚双侧 Fisher 精确检验 (Two-Tailed Fisher's Exact Test)
列联表定义：
$$\begin{pmatrix} R_f & R_r \\ A_f & A_r \end{pmatrix}, \quad r_1 = R_f + R_r, \quad r_2 = A_f + A_r, \quad c_1 = R_f + A_f, \quad c_2 = R_r + A_r, \quad n = r_1 + r_2$$
在超几何分布零假设下，观测值 $x = R_f$ 的概率为：
$$P(x) = \frac{\binom{r_1}{x} \binom{r_2}{c_1 - x}}{\binom{n}{c_1}}$$
双侧 $p$ 值按 R 语言 `fisher.test` 准则求和：
$$p = \sum_{x: P(x) \le P(x_{\text{obs}}) \times (1 + 10^{-7})} P(x)$$
- **防溢出实现**：在对数空间利用对数阶乘 $\ln(k!)$ 计算 $\ln P(x) = \ln \binom{r_1}{x} + \ln \binom{r_2}{c_1 - x} - \ln \binom{n}{c_1}$。$k < 256$ 使用预计算阶乘表，$k \ge 256$ 使用 Stirling 高阶渐近展开，保证误差 $< 10^{-14}$，严禁直接浮点阶乘运算。
- 测试基准向量：$\begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix} \implies p = 0.4857143$。

### 2.3 测序背景噪声单侧二项检验 (Upper-Tail Binomial Error Test)
$H_0: p \le e_0$ 对择 $H_1: p > e_0$。
$$p_{\text{err}} = P(X \ge A \mid N, e_0) = \sum_{k=A}^{N} \binom{N}{k} e_0^k (1 - e_0)^{N - k}$$
- 本地背景错误率：$e_0 = \max(e_{\text{floor}}, \hat{e}_{\text{window}})$，其中 $e_{\text{floor}} = 0.005$（SNP/SUB）或 $0.010$（INDEL）；$\hat{e}_{\text{window}}$ 为 $\pm 50\text{ bp}$ 窗口内非 REF、非候选 ALT 噪声错配率。
- 测试基准向量：$A = 3, N = 20, e_0 = 0.005 \implies p_{\text{err}} \approx 1.34 \times 10^{-4}$。
- 多重假设检验校正：全样本候选集采用 Benjamini–Hochberg 获得 $q_{\text{err}}$，门槛为 $q_{\text{err}} \le 0.01$。

### 2.4 亲本背景差分与亚克隆分类器 (Differential Calling & Classification)
差分列联表：$\begin{pmatrix} A_e & R_e \\ A_s & R_s \end{pmatrix}$，应用双侧 Fisher 检验获得 $p_{\text{diff}}$，经 BH 校正获得 $q_{\text{diff}}$。
差分频率 $\Delta = \text{AF}_e - \text{AF}_s$，其置信区间采用 Newcombe 混合评分法：
$$\Delta_L = \Delta - \sqrt{(\text{AF}_e - L_e)^2 + (U_s - \text{AF}_s)^2}$$
$$\Delta_U = \Delta + \sqrt{(U_e - \text{AF}_e)^2 + (\text{AF}_s - L_s)^2}$$
当 $A_s = 0$ 时，记录亲本排除上限 $U_s$（在当前深度下未检出亲本突变的置信上限）。

#### 状态转移与分类决策表 (严格按顺序判定)
1. `INDETERMINATE_STARTER_DEPTH`: $N_s < 20$。亲本深度不足，无法判断是否已存在。
2. `EDITED_SPECIFIC`: 编辑样本检出、亲本样本未检出、$q_{\text{diff}} \le 0.01$、$\Delta \ge 0.05$ 且 $\Delta_L > 0$。
3. `EDITED_ENRICHED`: 编辑与亲本均检出、$q_{\text{diff}} \le 0.01$、$\Delta \ge 0.05$ 且 $\Delta_L > 0$（富集突变，非 de novo）。
4. `PRE_EXISTING_IN_STARTER`: 亲本已检出且不属于富集突变（作为亲本背景从编辑突变列表中扣除）。
5. `DEPLETED_IN_EDITED`: 亲本检出、$q_{\text{diff}} \le 0.01$ 且 $\Delta \le -0.05$（编辑样本中显著稀释/缺失）。
6. `AMBIGUOUS_STARTER_TRACE`: 编辑检出、亲本有痕量读段（$A_s \ge 1$）但差分不显著。
7. `NOT_SUBCLONAL`: 其他情况，不予报告。

---

## 3. 领域数据模型扩展 (Data Models)

### 3.1 `SubclonalMetrics`
```rust
/// Detailed quantitative metrics for subclonal (low-frequency) variant observations.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SubclonalMetrics {
    /// Observed alternate allele frequency (ALT reads / total depth).
    pub alt_fraction: f64,
    /// 95% Wilson Score confidence interval [lower, upper].
    pub ci_95_lower: f64,
    pub ci_95_upper: f64,
    /// Forward and reverse supporting reads for ALT allele.
    pub alt_forward_reads: u32,
    pub alt_reverse_reads: u32,
    /// Total forward and reverse coverage at locus.
    pub total_forward_reads: u32,
    pub total_reverse_reads: u32,
    /// Strand bias p-value (Fisher exact test).
    pub strand_bias_pvalue: f64,
    /// Estimated local background sequencing error rate.
    pub background_error_rate: f64,
}
```

### 3.2 CLI 接口规范
- `--subclonal`: 启用低频亚克隆变异审计模式（默认关闭）。
- `--min-subclonal-af <FLOAT>`: 最低等位基因频率过滤阈值（默认 `0.05`，有效范围 `0.01..=0.50`）。
- `--min-subclonal-alt-reads <INT>`: 最少支持读段数（默认 `3`）。
- `--max-strand-bias-pvalue <FLOAT>`: 链偏倚过滤 $p$ 值阈值（默认 `0.001`，低于此值标记 `STRAND_BIAS` 软警示）。


---

## 4. 非处方性报告与不确定度警示规范 (Reporting & Language Constraints)

### 4.1 报告中强制包含的科学不确定度警示（Verbatim Caveat Block）
亚克隆模式报告必须在小节首部原样包含以下说明，严禁因果化断言：
> **Subclonal audit mode reports read-level ALT fractions at sites outside the consensus threshold. These are audit observations, not validated variant calls. The AF reporting threshold (default 0.05) is a configurable setting, not an experimental sensitivity floor. Effective detection limits depend on sequencing depth, strand balance, mapping quality, and the local sequencing error rate, and may be well above the configured threshold. Absence of an allele in the starter sample is reported with its exclusion bound at the observed depth and does not establish absence. Statistical tests assume independent reads; PCR duplicates, library artifacts, and systematic errors can violate this assumption.**

### 4.2 严格禁止的因果宣称词汇 (Prohibited Wording)
- 严禁宣称：“Detects variants down to 5%”、“sensitivity of 5%”、“validated/confirmed subclonal variant”、“true heterogeneity”、“population frequency”、“mutation rate”；
- 严禁直接断言编辑样本中特异的低频突变“必然是由本次编辑工具诱发”。

---

## 5. 工作包分解矩阵 (Work Packages)

| 任务编号 | 工作包名称 | 责任模型 | 交付物 | 验收条件 |
|---|---|---|---|---|
| **`M11-T01`** | **M11 科学架构设计与低频统计模型契约** | Claude Opus 5.5 / Controller | `M11_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`, 统计模型数学推导 | 确立链特异二项检验、Wilson CI、链偏倚 Fisher 检验与假阳性控制算法；完成台账初始化。 |
| **`M11-T02`** | **低频 Pileup 统计与质量过滤引擎实现** | Developer P1 (OpenCode GLM-5.3) | `crates/prokadiff-evidence/src/engine/subclonal.rs` | 实现链特异读段计数、本地错误率估计、Fisher 链偏倚计算与 Wilson CI；单元测试全绿。 |
| **`M11-T03`** | **SubclonalVariant 领域模型与分类器聚合集成** | Developer P1 (OpenCode GLM-5.3) | `crates/prokadiff-classify/src/subclonal.rs` | 将 `SubclonalMetrics` 挂载至变异模型，区分共识变异与低频变异，支持亲本低频背景比对扣除。 |
| **`M11-T04`** | **非处方性报告与 TSV 输出扩展** | Developer P1 (OpenCode GLM-5.3) | `crates/prokadiff-report/src/subclonal.rs`, `tables.rs` | `post_edit_variants.tsv` 新增频率与 CI 列；`report.md` 渲染亚克隆小节，严格附带非处方性与不确定度警示。 |
| **`M11-T05`** | **合成梯度稀释测试套件 (Titration Test Suite)** | Developer P1 (OpenCode GLM-5.3) | `crates/prokadiff/tests/subclonal_titration.rs` | 构建 5%, 10%, 25%, 50% 梯度稀释合成夹具；验证检出敏感度、假阳性控制及默认模式无回归。 |
| **`M11-T06`** | **M11 阶段集成与 Codex 独立审查验收** | Codex Sol (`gpt-6.1-sol`) / Controller | `M11_FINAL_INTEGRATION_REPORT.md` | Gate A–Gate G 全绿；工作区测试全过；Release Gate 3 gates PASS；Codex 独立审查签收。 |

---

## 6. 验收门禁 (Quality Gates)

- **Gate A (Non-Regression)**: 默认模式（未开启 `--subclonal`）下，全工作区现有 463+ 测试 100% 通过，输出逐字节无变化。
- **Gate B (Opt-in Precision)**: 开启 `--subclonal` 后，能准确识别合成梯度的低频突变（0.05, 0.10, 0.25, 0.50），且 0% 阴性对照位点不发生误报。
- **Gate C (Scientific Language)**: 发布门禁 `scripts/audit_release_gates.py` 严格通过，0 违规因果词汇，显式标注低频检测的物理与测序限制。
- **Gate D (HPC Execution)**: 所有 Rust 构建、测试与合成数据比对仅在 Slurm `qfree` 运行，严禁在未授权前提交 `qcpu_18i` 大规模分析。
- **Gate E (Independent Review)**: 经过 Codex Sol (`gpt-6.1-sol`) 独立代码与契约审查签收。

---

## 7. Gate A 验收基准修订 (User-Authorized Amendment — 2026-11-05)

- **状态**: `AUTHORIZED_BY_USER` (2026-11-05, FINAL_AUTHORIZATION_AUTHORITY)。
- **背景**: 原 Gate A 要求默认模式输出与 pre-M11 状态"逐字节无变化"。M8–M11 全部按治理规则以未提交工作树方式开发（BASE_MAIN_SHA `1a0c7b5` 之后无任何提交），故不存在可比较的 pre-M11 二进制或源树：GPFS 无用户文件系统快照；OpenCode snapshot 存储中不存在携带 M10-final 标记（`M10_FINAL_INTEGRATION_REPORT.md` 存在且 `M11_TASK_LEDGER.md` 不存在）的根树；Codex 会话仅包含 M11 之后的文件读取。Codex Sol Round 5 复审（Gate B/C/D PASS）明确裁定：要么提供可信 pre-M11 对比，要么"obtain explicit authorization to revise Gate A's acceptance basis"。
- **用户裁决**: 2026-11-05 用户（`FINAL_AUTHORIZATION_AUTHORITY`）在 Controller 呈报完整证据与三个选项后，明确选择 **"Revise Gate A basis"**（修订 Gate A 验收基准）。
- **修订后的 Gate A 验收基准**（取代上方 §6 Gate A 的"输出逐字节无变化"字面要求）:
  1. **Opt-in 隔离的代码级验证**：默认模式代码路径不含任何 subclonal 处理（`product_output_names` 在 `subclonal_enabled=false` 时排除 `subclonal_audit.tsv`；`write_markdown_report` 仅在启用时渲染亚克隆小节）；
  2. **全部 pre-existing 测试无改动且全绿**（mtime 证据 + 全量套件通过）；
  3. **committed pre-M11 oracle 套件通过**（58 parity + 3 layer-2 + 195 release-gate，oracle 粒度字节比较）；
  4. **当前二进制默认输出 provenance 哈希固定**（`.prokadiff/m11_gate_a_default_outputs.txt`：命令行、二进制 sha256、Slurm job、10 产品哈希、`subclonal_audit.tsv` 缺席验证）；
  5. **显式记录残余局限**：与 M10-final 二进制的字面 byte-diff 不可重建，如实记录于验证日志。
- **基线提交义务**: 本次修订同时授权将当前快照提交为**首个内容寻址基线**（本地提交，不推送）；自 M12 起，每个里程碑必须在新开发开始前提交内容寻址基线，使字面 byte-diff 比较成为标准实践。
- **记录文件**: `.prokadiff/m11_gate_a_authorization.md`、`.prokadiff/m11_r4_repair_validation.log`（ERRATUM 附注）。
