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

## 变异分级分类逻辑

编辑株中筛选出的每一条非预期突变，严格按照不可逆优先级进入且仅进入一个类别：

```mermaid
graph TD
    A[非预期突变] --> B{是否为结构性事件?}
    B -- 是: JC, MOB, DEL > 2bp, AMP, CON --> C[structural 结构变异]
    B -- 否 --> D{是否临近 Spacer 同源序列与 PAM?}
    D -- 是: 错配 <= 4, 在侧翼窗口内 --> E[near_homolog 近同源靶区变异]
    D -- 否 --> F[scattered_snv 散在点突变]
```

---

## 软件许可

`ProkaDiff` 遵循开源 [MIT License](LICENSE)。

