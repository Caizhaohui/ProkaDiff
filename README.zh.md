# ProkDiff

[English README](README.md)

**原核基因编辑质控用全基因组差分：出发株 vs 编辑株 WGS**

面向原核基因编辑质控的全基因组非预期变异检测器（Rust）。出发株与编辑株都做 Illumina WGS；NCBI / 近缘基因组只作**坐标骨架**；非预期 = 编辑株相对出发株多出来的突变，再按三类标注。

这不是人细胞 GUIDE-seq，也不是「相对 MG1655 的全部 SNP 清单」。

## 定位


| 项 | 说明 |
| --- | --- |
| 展示名 | **ProkDiff**（CLI：`prokdiff`） |
| 第一期读长 | Illumina 短读，克隆 consensus |
| 引擎 | 外部 **Bowtie2** + Rust 实现 RA / MC / JC（清洁室对标 [breseq methods](https://breseq.barricklab.org/0.40.1/methods/)）。系统 `breseq` **仅作 oracle**，不是运行时引擎。 |
| 差分 | `M_edited \ M_parent \ {intended}` |
| 编辑器（v0.1） | `--editor cas9 \| cas12a \| dsb`（`cast` / `is110` 拒绝并指向路线图） |

## 双样本强制

`--starter` 与 `--edited` FASTQ **都必填**。缺任一株 → 非零退出，并写明必须提供出发株 WGS。不得把相对 NCBI 的全部差异当成非预期编辑。

- 同一选项可重复：**连续两个文件 = 一对 PE**（先 R1 后 R2）；单个文件 = SE。
- `--intended` 可选；缺省则差分后全部进入非预期列表。
- `--ref` 可重复（染色体 + 质粒 + donor 骨架）。

## 三类标签（观察定义）

判定顺序锁定：**结构 (3) → 近同源 (1) → 散在 (2)**。一条突变只进一类。**不得**因「该编辑器不做 DSB」关闭散在突变。

1. **`near_homolog`** — 落在 spacer 近同源 + 该编辑器 PAM 附近（Cas9 默认 NGG；Cas12a 默认 TTTV；最多 4 个错配；默认 50 bp）。注释，非唯一证据。`--editor dsb` 不打此类。
2. **`scattered_snv`** — 其余点突变 / 短 indel。所有编辑器都要报。可选假说 `sos_widney2014`（`--no-hypothesis` 可关）。
3. **`structural`** — MOB / JC、大片段缺失（>2 bp）、AMP / CON。

## CLI

产品入口（两株必填；差分后再分类）：

```bash
prokdiff \
  --starter starter_R1.fastq.gz --starter starter_R2.fastq.gz \
  --edited  edited_R1.fastq.gz  --edited  edited_R2.fastq.gz \
  --ref     NC_000913.3.fa \
  --intended intended.tsv \
  --editor cas9 \
  --spacer  NNNNNNNNNNNNNNNNNNNN \
  --threads 8 \
  --outdir  results/prokdiff
```

写出 `starter.gd`、`edited.gd`、`unintended.tsv`、`summary.txt`。

`summary.txt` 记录是否提供了 `--intended`（`intended_provided`）、声明行数（`intended_declared`）、掩码命中条数（`intended_observed`），以及汇总状态 `intended_status`（`all_observed` / `partial` / `none_observed` / `NA`）。TSV 中 SNP 的 `ref` 为参考碱基；JC 另有 `side2_seq_id` / `side2_position`。完整契约见 [docs/schema.md](docs/schema.md)。

单样本引擎（oracle / 对拍任务）：

```bash
prokdiff evidence \
  --ref genome.fa \
  --fastq sample_R1.fastq --fastq sample_R2.fastq \
  --threads 8 \
  --outdir results/sample
```

输入输出契约见 [docs/schema.md](docs/schema.md)。与 breseq 的结果 / 墙钟 / 峰值 RSS 对拍规则见 [docs/parity.md](docs/parity.md)。

## v0.1 明确不做

- CAST / bRNA-IS110 专用近同源与 donor 解释
- 混群 polymorphism、长读长、GUIDE-seq
- 运行时包装系统 `breseq` 当引擎
- 内置组装；重写 Bowtie2 / HMMER / ESM
- 把完整大肠深度 FASTQ 或 Widney/Copley SRA（PRJNA1088182）推进 git
- 把相对 NCBI 的全部差异报成编辑脱靶

## 当前进度

阶段 4 报告层：`summary.txt` 含目的编辑声明/观察/状态；`unintended.tsv` 写出 SNP 参考碱基与 JC 第二端坐标。引擎 MVP（`prokdiff evidence`）不变。层 0 测试为纯内存。层 1 / 层 2 引擎作业必须提交 Slurm **`qcpu_18i`**。

## 许可

MIT。方法学：[breseq 0.40.1](https://breseq.barricklab.org/0.40.1/methods/)。科学语义：Widney et al. 2024 bioRxiv [10.1101/2024.03.19.584922](https://doi.org/10.1101/2024.03.19.584922)。
