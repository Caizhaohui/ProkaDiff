# ProkDiff

**Prokaryotic genome diff for gene-editing QC — starter vs edited WGS**

面向原核基因编辑质控的全基因组非预期变异检测器（Rust）。出发株与编辑株都做 Illumina WGS；NCBI / 近缘基因组只作**坐标骨架**；非预期 = 编辑株相对出发株多出来的突变，再按三类标注。

当前仓库：**阶段 0–1**（Cargo workspace + Genome Diff I/O + 层 0 RA/JC 单元测试）。引擎 / 层 1 对拍见阶段 2。

这不是人细胞 GUIDE-seq，也不是「相对 MG1655 的全部 SNP 清单」。

## 定位

| 项 | 说明 |
| --- | --- |
| 展示名 | **ProkDiff** |
| CLI | `prokdiff` |
| 第一期读长 | Illumina 短读，克隆 consensus |
| 引擎 | 外部 **Bowtie2** + Rust 实现 RA / MC / JC（对标 breseq consensus；breseq `.gd` 仅作 oracle） |
| 差分 | `M_edited \ M_parent \ {intended}` |
| 编辑器 | `--editor cas9 \| cas12a \| dsb`（`cast` / `is110` 第一期拒绝） |

## 三类标签（观察定义）

判定顺序锁定：**结构 (3) → 近同源 (1) → 散在 (2)**。一条突变只进一类。**不得**因「编辑器不做 DSB」关闭散在突变。

1. **`near_homolog`** — 落在 spacer 近同源 + PAM 附近（Cas9 默认 NGG；Cas12a 默认 TTTV）。注释，非唯一证据。
2. **`scattered_snv`** — 其余点突变 / 短 indel。所有编辑器都要报。可选假说列可关（`--no-hypothesis`）。
3. **`structural`** — MOB / JC、大片段、质粒相关等。

## 快速开始（计划中的 CLI）

阶段 0 脚手架尚未落地；下列为第一期约定接口：

```bash
prokdiff \
  --starter starter_R1.fastq.gz --starter starter_R2.fastq.gz \
  --edited  edited_R1.fastq.gz  --edited  edited_R2.fastq.gz \
  --ref     NC_000913.3.gbk \
  --ref     plasmid.gbk \
  --intended intended.tsv \
  --editor cas9 \
  --spacer  NNNNNNNNNNNNNNNNNNNN \
  --pam     NGG \
  --threads 8 \
  --outdir  results/prokdiff
```

- `--starter` 与 `--edited` **都必填**。缺出发株 → 退出，不会把相对 NCBI 的差异当成脱靶。
- 同一选项连续两个文件 = 一对 PE（R1 后 R2）；单文件 = SE。
- `--intended` 可选；缺省则差分后全部进入非预期列表。

输出：`unintended.tsv`、`starter.gd` / `edited.gd`、摘要。契约见 [docs/schema.md](docs/schema.md)。

## 明确不做（第一期）

- CAST / bRNA-IS110 专用近同源与 donor 解释
- 混群 polymorphism、长读长、GUIDE-seq
- 运行时包装系统 `breseq` 当引擎（breseq 只作对拍 oracle）
- 内置组装；重写 Bowtie2 / HMMER / ESM
- 把完整大肠深度 FASTQ 或 Widney SRA（PRJNA1088182）推进 git

## 架构（简图）

```text
starter FASTQ ──┐
                ├── Bowtie2 → RA / MC / JC → starter.gd ─┐
edited  FASTQ ──┘         (same refs)      edited.gd  ─┼─ subtract → classify → unintended.tsv
ref GBK/FASTA ─────────────────────────────────────────┘
intended.tsv (optional) ────────────────────────────────┘
```

## 开发与测试（HPC）

登录节点只做编辑与 `fmt` / `clippy` / `check`。合成夹具、真实 WGS、以及相对 breseq 的 **结果 + wall-clock + peak RSS** 对照一律提交 **`qcpu_18i`**，并给足 CPU / 内存。细节见 [DEVELOPMENT.md](DEVELOPMENT.md) 与 [docs/parity.md](docs/parity.md)。

## 文档索引

| 文件 | 用途 |
| --- | --- |
| [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) | 完整规格与分阶段任务 |
| [AGENTS.md](AGENTS.md) | 给 agent 的硬规则 |
| [DEVELOPMENT.md](DEVELOPMENT.md) | 工具链、测试分层、HPC |
| [docs/schema.md](docs/schema.md) | CLI / TSV / GD 契约 |
| [docs/parity.md](docs/parity.md) | 与 breseq 结果 / 时间 / 资源对拍规则 |
| [原核基因编辑脱靶检测-开发计划.md](原核基因编辑脱靶检测-开发计划.md) | 旧单文件规格 → 指针 |

## 许可与致谢

建议许可证：**MIT**（清洁室实现；不复制 breseq 源码）。方法学参考 [breseq methods](https://breseq.barricklab.org/0.40.1/methods/)；科学语义参考 Widney et al. 2024 bioRxiv [10.1101/2024.03.19.584922](https://doi.org/10.1101/2024.03.19.584922)。
