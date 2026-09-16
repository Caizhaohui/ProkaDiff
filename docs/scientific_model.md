# ProkaDiff 科学解释模型与证据规范
## Scientific Interpretation Model & Evidence Hierarchy

> **版本：** v0.2.0-draft  
> **适用范围：** 原核工程菌株全基因组编辑审计（Post-edit Genome Audit & Edit Verification）  
> **核心准则：** Observation first, attribution second, causality last.

---

## 1. 核心定位与回答的核心问题

ProkaDiff 并不是一个先入为主的“CRISPR 脱靶预测器”（CRISPR off-target predictor），而是：
> **面向原核工程菌株的高性能、WGS-first 的编辑验证与全基因组完整性审计框架。**  
> *(A high-performance, WGS-first framework for post-edit genome auditing and edit verification in engineered prokaryotes)*

在工程菌株（如大肠杆菌、放线菌、枯草芽孢杆菌等）基因编辑质控中，核心回答两个基本问题：
1. **这个工程菌株是否获得了预期的编辑？结构是否完整合规？**
2. **除此之外，它的整个基因组还发生了什么变化？与编辑系统是否存在序列或机制上的关联？**

---

## 2. 观察优先的分析顺序（Pipeline Order）

严禁“先根据 gRNA 预测全基因组脱靶位点，再去测序数据中强行寻找支持”的伪闭环。ProkaDiff 必须执行严格的**观察优先（Observation-First）**分析流向：

```text
[1. Observed Genome Changes]
    ├── Starter WGS vs. Edited WGS 严格集合差
    └── 检出所有真实存在的新增变异 (SNVs, indels, large deletions, novel junctions, MOB)
            ↓
[2. Edit Verification]
    ├── 检验声明的预期编辑 (Intended edits)
    └── 判定状态 (Complete / Partial / Missing / UnexpectedStructure)
            ↓
[3. Genome-wide Annotation]
    ├── 物理尺度分类 (Small variants vs. Structural variants)
    └── 基因组元件注释 (CDS, tRNA, rRNA, 移动元件 IS, 重复序列)
            ↓
[4. Mechanistic Association]
    ├── 候选向导依赖位点关联 (Candidate guide-dependent association)
    └── 移动元件活化、重组断点与供体/载体序列关联
            ↓
[5. Human-readable Interpretation]
    ├── 按复核优先级 (ReviewPriority) 组织报告
    └── 附带明确的科学局限性 (Limitations) 与技术追溯 (Provenance)
```

---

## 3. 科学证据四层级（Evidence Hierarchy）

在配对 WGS 差分分析中，不同结论具有不同的科学确定性，必须在数据结构与报告中严格区分：

| 层级 | 概念名称 | 科学内涵与支撑依据 | ProkaDiff 报告规范 |
| :--- | :--- | :--- | :--- |
| **Level 0** | **Direct Observation**<br>(直接观察事实) | 编辑株中高置信度检出，且亲本出发株中确认缺失的基因组变异。由 Illumina 读段覆盖与断点接头直接支撑。 | 统一称为 **Post-edit Differential Variant**。<br>不能直接标注为“编辑引起的突变”。 |
| **Level 1** | **Contextual Association**<br>(序列与空间关联) | 变异位点与特定基因组上下文存在几何或同源关联（如：邻近向导同源区、位于已知 IS 元件插入热点、处于反向重复区）。 | 标注为关联特征（如 `distance_to_site`, `pam`, `tsd`）。<br>表明相关性，不表明因果性。 |
| **Level 2** | **Mechanistic Hypothesis**<br>(生物学机制假说) | 基于文献机理提出的可能归因假说（如：向导依赖误切候选、转座子应激活化、DSB 修复伴生变异）。 | 统一使用 **Candidate** 或 **Hypothesis** 前缀。<br>严禁使用“Confirmed Off-target”。 |
| **Level 3** | **Validated Causality**<br>(实验证实的因果) | 通过独立阴性对照、多克隆生物学重复、反向互补拯救或体外切割实验证实变异直接由编辑系统造成。 | **ProkaDiff 默认不能单凭 Pair WGS 给出此结论**，保留给用户后续实验复核。 |

---

## 4. 核心实体解释规范

### 4.1 预期编辑（Intended Edits）
对用户通过 `--intended` 声明的编辑，系统不再简单进行布尔匹配（Mask / Not mask），而是生成多维验证评估：
- **`Complete`**：预期的等位基因或结构完整检出，断点精准，无非预期异常结构。
- **`Partial`**：仅检测到部分预期事件（如整合盒仅检出左侧接头、大片段缺失断点偏差过大）。
- **`Missing`**：目标位点完全维持亲本型，未检测到任何预期编辑事件。
- **`UnexpectedStructure`**：目标位点检出复杂重排、异常外源片段插入或二次接头，提示靶向编辑异常。

### 4.2 向导依赖候选误切（Candidate Guide-dependent Off-target）
- **术语规范：** 一律使用 `candidate guide-dependent off-target`，严禁使用 `confirmed off-target`。
- **判定标准：** 仅当编辑株新增小变异（SNP/indel）满足：
  1. 空间距离在向导候选靶点窗口内（默认 $\le 50\text{ bp}$）；
  2. 靶点与 spacer 错配数符合阈值（默认 $\le 4$）；
  3. 具有匹配的有效 PAM；
  此时赋予候选关联，并在报告中显式注明：“*Spatial association is consistent with a possible guide-dependent event but does not establish cleavage causality.*”

### 4.3 附带基因组变异（Collateral Genome Changes）
对于远离向导靶点的全基因组非预期差异，统一称为 **Collateral Genome Changes**，可能来源包括：
- 菌株培养传代漂移（Culture drift）；
- 电转/热激等转化物理应激；
- 抗生素或致死筛选压力；
- 质粒复制与外源蛋白表达负担（Metabolic burden）；
- 重组酶系统（如 $\lambda$-Red）诱发的复制叉停滞；
- 自发突变（Spontaneous mutation）。

**铁律：严禁在单克隆测序中自动标注“该散在突变由 Cas9 DSB / SOS 诱发”。** 假说列仅可作为可选辅助参考，不可作为既定事实。

### 4.4 移动元件事件（Mobile-element / IS Events）
在原核工程菌（如大肠杆菌 MG1655、BL21）中，插入序列（IS 元件，如 IS1、IS5 等）的跳跃往往比微小点突变具有更显著的表型破坏性。
- ProkaDiff 重点提取接头（JC）、靶位点重复（TSD）与元件家族。
- 客观描述事件本身及其是否为编辑后独有，不妄断转座诱发诱因。

---

## 5. 人工复核优先级（ReviewPriority）与无打分准则

### 5.1 复核优先级（ReviewPriority）
为了辅助实验人员快速排查风险，系统提供三级排查优先级，**仅表达人工核对的建议顺序，不表达菌株的安全性质**：
1. **`HighAttention`**：
   - 预期靶位点出现异常结构（`UnexpectedStructure`）；
   - 新增移动元件插入（`MOB`）；
   - 大片段结构变异（大型缺失、染色体重排、异常接头）；
   - 高度同源的向导依赖候选误切事件。
2. **`Review`**：
   - 编码区非同义突变；
   - 散在点突变及小 indel；
   - 部分达成的预期编辑（`Partial`）。
3. **`Info`**：
   - 完整达成的预期编辑（`Complete`）；
   - 良性同义突变或非关键区域变异。

### 5.2 严禁伪科学打分
ProkaDiff 明确**禁止**开发或输出类似：
- `0–100 Genome Quality Score`（基因组质量综合分）
- `Genome Safety Score`（基因组安全指数）
- `Safe / Unsafe`（安全/不安全二元标签）

菌株是否可用、是否能进入下一轮发酵或功能测试，取决于下游具体工程目标，生物信息工具不得替代生物学专家妄下定论。

---

## 6. 算法门禁机制（CFD / Hsu / Bulge Policy）

依据 Clean-room 原则及科学诚信规范：
1. **真实 Oracle 门禁：** 任何评分算法（如 CFD 特异性评分、Hsu 矩阵打分、Bulge 允许匹配）在未接入公认外部基准（如标准 Cas-OFFinder/CFD 官方实现）并完成跨平台双盲对拍前，**必须处于未激活（Disabled / NA）状态**。
2. **严禁自创加权公式：** 绝对禁止开发任何未经同行评审实证的启发式伪打分算法（Heuristic seed scores）。
3. **透明声明：** 若评分功能未激活，输出数据与报告必须显式注明 `NA` 或 `DISABLED_UNVALIDATED_ORACLE`，并在报告《局限性》章节主动说明原因。
