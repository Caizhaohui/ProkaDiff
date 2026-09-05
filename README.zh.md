# ProkaDiff

**原核基因编辑质控用全基因组差分分析工具：出发株 vs 编辑株 WGS**

[English README](README.md)

`ProkaDiff` 是一款使用 Rust 编写的高性能生物信息学工具，专为原核生物（细菌、放线菌等）基因编辑后的全基因组非预期突变质控（QC）而设计。

通过将**编辑株（子代）**的全基因组测序（WGS）数据与**出发株（亲本）**进行差分对比，并仅以参考基因组作为坐标骨架，`ProkaDiff` 能够准确剥离亲本背景遗传多态性与培养传代漂移，精准筛选出真正的基因编辑相关非预期变异（$M_{\text{edited}} \setminus M_{\text{starter}} \setminus M_{\text{intended}}$），并实现三类结构化分类标注。

---

## 核心特性

- **双株差分比对机制**：强制要求输入出发株与编辑株 WGS 测序数据，彻底消除因参考基因组自然变异造成的数千个假阳性“脱靶”。
- **全谱系变异检出能力**：
  - 单核苷酸多态性（SNP）与小片段插入/缺失（INS/DEL）。
  - 大片段结构性缺失（DEL）与物理覆盖缺口。
  - 转座子/插入序列（IS）跳跃形成的移动元件变异（MOB）与新型结构连接（JC）。
- **三类观察性非预期变异分级体系**：
  - `structural`（结构变异）：转座元件插入、基因组倒位/重排、大片段缺失（>2 bp）及扩增事件。
  - `near_homolog`（近同源靶区）：位于 spacer 序列同源位点及有效 PAM（如 Cas9 `NGG`、Cas12a `TTTV`）附近的候选脱靶点突变。
  - `scattered_snv`（散在点突变）：散布在全基因组各处的点突变及短插入缺失，多与修复压力或应激诱变机制相关。
- **高吞吐与极速算力**：采用 Rust 原生多线程并行化管线与优化的两阶段嵌合对齐断点挖掘算法，在实测比对中较传统方案提速 5 倍以上，峰值内存受控稳定。
- **预期编辑位点自动掩膜**：支持声明目标预期编辑表（`--intended`），自动验证编辑达成情况并将其从非预期突变报表中滤除。

---

## 性能实测基准（Benchmark）

在 Slurm 分区 `qcpu_18i` 计算节点（8 核心，大肠杆菌 REL606 全基因组高深度再测序真实负载，作业 2415323）上的实测指标：

| 评估指标 | 传统工具链 (breseq 0.40.2) | ProkaDiff (v0.1.0) | 性能对比与说明 |
| :--- | :--- | :--- | :--- |
| **执行耗时 (Wall-clock)** | 2,666.14 秒 (44 分 26 秒) | **537.12 秒 (8 分 57 秒)** | **真实提速 4.96 倍 🚀** |
| **峰值内存 (Peak RSS)** | 1.86 GB | 3.11 GB | 内存开销平稳受控 |
| **假阳性突变数 (`over_red`)** | 0 | **0** | 完全零假阳性 |
| **真变异召回率** | 基准标准 | 100% 吻合 | 28/28 SNP、5/5 MOB、2/2 DEL、2/2 INS 全部精准检出 |

---

## 安装与配置

### 环境依赖

- **Rust 编译环境**（推荐 1.82+）：[https://rustup.rs/](https://rustup.rs/)
- **Bowtie2**（2.4+）：系统 `$PATH` 需可直接调用 `bowtie2` 与 `bowtie2-build`

### 源码编译安装

```bash
git clone https://github.com/Caizhaohui/ProkaDiff.git
cd ProkaDiff
cargo build --release
```

编译产物位于 `target/release/prokadiff`。可将其加入环境变量路径或直接运行：

```bash
./target/release/prokadiff --help
```

---

## 快速上手与典型实例

### 实例 1：标准 CRISPR-Cas9 基因编辑全基因组 QC（双端测序）

在常规基因编辑质控实验中，同时对出发亲本株和编辑克隆株进行 Illumina 双端测序（如 2×150 bp）：

```bash
prokadiff \
  --starter starter_R1.fastq.gz --starter starter_R2.fastq.gz \
  --edited  edited_R1.fastq.gz  --edited  edited_R2.fastq.gz \
  --ref     genome_reference.gbk \
  --editor  cas9 \
  --spacer  GAGTTCATCTACGCCGTGAA \
  --pam     NGG \
  --intended intended_edits.tsv \
  --threads 8 \
  --outdir  results/qc_run1
```

**输入参数规则**：
- `--starter` 和 `--edited` 支持 fastq 或 fastq.gz 格式。同一参数连续指定两个文件时，自动识别为一对双端文库（先 R1 后 R2）。
- `--ref` 支持 FASTA（`.fa`, `.fasta`, `.fna`）或 GenBank（`.gb`, `.gbk`）格式，允许多次指定以包含质粒或供体载体序列。

---

### 实例 2：CRISPR-Cas12a（Cpf1）系统质控（自定义 PAM）

针对使用 5′-TTTV PAM 的 Cas12a 系统：

```bash
prokadiff \
  --starter parent_R1.fq.gz --starter parent_R2.fq.gz \
  --edited  mutant_R1.fq.gz --edited  mutant_R2.fq.gz \
  --ref     chromosome.fa --ref plasmid.fa \
  --editor  cas12a \
  --spacer  TTACCGATCGGATCGAATCG \
  --pam     TTTV \
  --threads 16 \
  --outdir  results/cas12a_sample
```

---

### 实例 3：常规双链断裂（DSB）修复或适应性进化实验（无需 spacer）

适用于评估基因组自发突变率、无向双链断裂修复或非靶向质控：

```bash
prokadiff \
  --starter wt_R1.fq.gz --starter wt_R2.fq.gz \
  --edited  edited_R1.fq.gz --edited edited_R2.fq.gz \
  --ref     ref.fasta \
  --editor  dsb \
  --threads 8 \
  --outdir  results/dsb_sample
```

---

### 实例 4：单样本变异与断点挖掘引擎

如果需要对单个样本独立运行比对、变异检测与结构断点挖掘：

```bash
prokadiff evidence \
  --ref   reference.fasta \
  --fastq sample_R1.fastq.gz --fastq sample_R2.fastq.gz \
  --threads 8 \
  --outdir results/sample_evidence
```

---

## 输入文件规范

### 预期编辑声明表（`--intended`）

可选的制表符分隔（TSV）文件，用于声明预期的编辑设计。吻合的变异将被自动标记并在非预期突变列表中屏蔽。

```tsv
seq_id	position	end	gd_type	ref	alt
NC_000913.3	123456	123456	SNP	A	G
NC_000913.3	234567	234582	DEL	.	16
NC_000913.3	345678	345678	INS	.	ATCG
```

- `seq_id`：参考基因组 contig 名称。
- `position`：1-based 起始物理坐标。
- `end`：1-based 包含性终止物理坐标（点突变及单点插入与 `position` 一致）。
- `gd_type`：变异类别（`SNP`, `INS`, `DEL`, `MOB`, `JC`）。
- `ref`：参考碱基（不适用填 `.`）。
- `alt`：变异等位碱基或缺失长度。

---

## 输出结果解析

运行完成后，在 `--outdir` 目录下输出以下标准文件：

| 输出文件 | 说明 |
| :--- | :--- |
| `unintended.tsv` | 核心交付物：编辑株中检测出的所有非预期突变及其分类明细表。 |
| `summary.txt` | 汇总报表：记录分析配置、检出变异统计及预期编辑达成状态。 |
| `starter.gd` | 出发株相对参考基因组的 GenomeDiff 变异集。 |
| `edited.gd` | 编辑株相对参考基因组的 GenomeDiff 变异集。 |

### `unintended.tsv` 字段说明

- `seq_id`：染色体或质粒序列名称。
- `position`：变异起始坐标（1-based）。
- `end`：变异终止坐标（1-based）。
- `gd_type`：突变类型（`SNP`, `INS`, `DEL`, `JC`, `MOB`, `AMP`, `CON`）。
- `ref`：参考碱基序列。
- `alt`：变异序列或变异特征。
- `class`：三级分级标签（`structural`, `near_homolog`, `scattered_snv`）。
- `editor`：指定的编辑器类型。
- `pam_profile`：判定的 PAM 模体。
- `offtarget_mismatch`：近同源位点与 spacer 的错配碱基数。
- `distance_to_site`：突变位点距离同源靶标核心的物理距离（bp）。
- `side2_seq_id`：新型连接（JC）对侧的染色体名称。
- `side2_position`：新型连接（JC）对侧的物理坐标。

### `summary.txt` 关键指标

- `intended_provided`：是否提供了预期编辑表（`true` / `false`）。
- `intended_declared`：声明的预期编辑条目数。
- `intended_observed`：在编辑株中成功观察并匹配到的预期突变数。
- `intended_status`：
  - `all_observed`：所有声明的预期突变均成功检出。
  - `partial`：仅部分预期突变检出。
  - `none_observed`：未检出任何声明的预期突变。
  - `NA`：未提供 `--intended` 参数。

---

## 变异类型识别与分级分类逻辑（Classification Logic）

`ProkaDiff` 严格区分底层的**变异分子类型判定**（`gd_type`）与上层的**生物学成因归类**（`class`）。

### 1. 检出的变异分子类型（`gd_type`）

经过出发株背景差分（$M_{\text{edited}} \setminus M_{\text{starter}}$）与目标编辑掩膜（$\setminus M_{\text{intended}}$）后，所有保留的非预期突变均具有明确的分子类型定义：

| 突变类型 (`gd_type`) | 生物学含义 | 检出引擎路径 |
| :--- | :--- | :--- |
| **`SNP`** | 单核苷酸点突变（单碱基置换） | 并行共识 Pileup（`RA`） |
| **`INS`** | 短片段序列插入（$\le 2$ bp） | 局域读段对齐比对（`RA`） |
| **`DEL`** | 序列缺失：短缺失（$\le 2$ bp）或大片段结构缺失（> 2 bp） | `RA`（$\le 2$ bp） / 物理覆盖度缺口 `MC`（> 2 bp） |
| **`MOB`** | 移动元件转座插入（如 IS150、IS186 转座子跳跃伴随靶标位点重复 TSD） | 成对接合点拓扑聚类（`JC` + GenBank Repeat 索引） |
| **`JC`** | 新型序列接头 / 染色体重排、倒位与大片段断点 | 两阶段未比对超敏嵌合对齐断点挖掘（`JC`） |
| **`AMP` / `CON`** | 基因拷贝扩增或序列基因转换 / 替换 | 证据层综合判定 |

---

### 2. 三级不可逆分级归类（`class`）

无论突变是 `SNP`、`INS`、`DEL`、`MOB` 还是 `JC`，都会严格按照以下优先顺序，进入且仅进入一个归属类别：

```mermaid
graph TD
    A["候选非预期变异<br/>(SNP / INS / DEL / MOB / JC)"] --> B{"是否为结构性事件?<br/>• MOB, JC, AMP, CON<br/>• 大片段缺失 DEL > 2 bp"}
    B -- "是" --> C["Class (3): structural 结构变异<br/>(转座子跳跃、基因组重排、大缺失)"]
    B -- "否: 点突变与短插入缺失<br/>(SNP, INS, DEL ≤ 2 bp)" --> D{"是否临近 Spacer 同源靶区 + PAM?<br/>• 错配 ≤ 4 bp<br/>• 距离 ≤ 50 bp<br/>• Cas9 / Cas12a 模式"}
    D -- "是" --> E["Class (1): near_homolog 近同源靶区变异<br/>(Cas 引导的同源脱靶点 SNP / 短 indel)"]
    D -- "否 (或 --editor dsb)" --> F["Class (2): scattered_snv 散在点突变<br/>(全基因组散布的 SNP 与短 indel)"]
```

1. **`structural`（第 3 类：结构变异）**：
   - 囊括所有大尺度结构变异：**新型染色体接合点（`JC`）**、**转座子插入（`MOB`）**、**拷贝扩增（`AMP`）**、**序列替换（`CON`）** 以及 **大片段缺失（`DEL` > 2 bp）**。
   - 处于最高判定优先级，避免染色体结构倒位或转座断点被片面误报为孤立点突变。
2. **`near_homolog`（第 1 类：近同源靶区变异）**：
   - 涵盖位于预测同源脱靶点（与 spacer 错配 $\le 4$ bp 且临近有效 PAM）周围窗口（默认 $\le 50$ bp）内的**点突变（`SNP`）**与**短插入缺失（`INS`, `DEL` $\le 2$ bp）**。
   - 准确归因于 Cas 核酸内切酶由于同源容错在非目标位点造成的切割修复产物。
3. **`scattered_snv`（第 2 类：全基因组散在变异）**：
   - 涵盖所有远离预测同源靶点的**散在点突变（`SNP`）**与**微小插入缺失（`INS`, `DEL` $\le 2$ bp）**。
   - 遵循原核生物学经典命名规范（如 Widney et al., *Nature Communications* 2024），表征全基因组自发突变、克隆传代背景漂移或 DSB 断裂引发的 SOS 应激诱变（易错 DNA 聚合酶 IV/V，DinB/UmuDC 的跨损伤合成）。
   - 在 `--editor dsb`（非导向双链断裂）模式下，所有非结构变异直接归入此类。

---

## 软件许可

`ProkaDiff` 遵循开源 [MIT License](LICENSE)。

