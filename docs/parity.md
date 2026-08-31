# ProkDiff ↔ breseq 对拍规则

**Oracle + 性能对照：** 钉死版本的 **breseq 0.40.x** + 配套 bowtie2（写入 `testdata/VERSIONS.txt`）。同一参考、同一 FASTQ、**同一计算节点**。

对拍三维（缺一不可）：

1. **结果** — Genome Diff 子集 + `gdtools SUBTRACT`
2. **时间** — wall-clock（总时间；可选分步 RA / MC / JC）
3. **资源** — peak RSS（可选 CPU 时间、线程数）

**运行时：** ProkDiff **不**调用 `breseq`。breseq 只作 oracle 与对照跑次。

**执行环境：** 禁止在登录节点跑对拍 / 计时。提交 Slurm **`qcpu_18i`**，给足 `cpus-per-task` 与内存，使 ProkDiff `--threads` 与 breseq `-j`（或等价）对齐。实测写入 `benchmark/results/`；禁止臆造数字。

实现方法依据：[breseq 0.40.1 methods](https://breseq.barricklab.org/0.40.1/methods/)。清洁室：不复制 breseq 源码。

## 单样本引擎对拍（结果）

```text
gdtools SUBTRACT rust.gd breseq.gd   # 应为空（无多报）
gdtools SUBTRACT breseq.gd rust.gd   # 应为空（无漏报）
```

在下方允许差异规范化之后，两个方向都应为空（或仅剩 parity 白名单内条目）。

## 时间与资源记录

建议每行一条 workload（CSV 或表）：

```text
workload,tool,threads,wall_s,max_rss_kb,node,breseq_version,bowtie2_version,notes
synth_snp_indel,prokdiff,8,<sec>,<kb>,<hostname>,…,…,
synth_snp_indel,breseq,8,<sec>,<kb>,<hostname>,…,…,
```

规则：

- 两边同一节点、同一输入、同一线程设定；先后顺序固定并在 notes 中写明。
- 申请资源 ≥ 所用线程；勿欠申请导致排队抖动或 OOM。
- 先有实测再写 README 加速比；规格不写死「快 N×」。
- 缺少路径或二进制时：停止该层，不伪造表。

## 允许的差异（结果白名单）

下列差异 **不**构成结果失败；应在测试中规范化或忽略：

| 类别 | 说明 |
| --- | --- |
| UN / 未知区 | `UN` 证据不对拍失败 |
| HTML / R 图 | breseq 报告产物；ProkDiff 第一期不要求 |
| BAM 文件名 / 路径 | 中间比对文件命名可不同 |
| 重复区 JC 拷贝歧义 | 「哪一个拷贝」与 breseq 一样可有 unassigned JC；允许 junction 坐标 ± 数 bp 的规范化（须在夹具说明中写明阈值） |
| consensus 频率舍入 | `frequency` / 相关浮点字段在 consensus 下接近 1 即可，允许舍入差 |

## 不允许的差异（结果红线）

下列任一出现即结果失败：

| 类别 | 说明 |
| --- | --- |
| SNP 位点 / 等位基因不一致 | 同一位点 ref/alt 与 oracle 不符 |
| unique-mapping MOB 漏报 | 明确可 unique 映射的移动元件插入漏报 |
| 大 DEL 漏报 | 明确的大片段缺失漏报（MC/JC 路径） |
| 白名单外的系统性多报 | `SUBTRACT rust.gd breseq.gd` 留下非白名单条目 |

时间 / RSS 无「必须更快」硬阈值；但 **必须有可复核记录**。若 ProkDiff 明显更慢或吃内存异常，记入 notes，不得静默省略。

## 产品层（双样本）

breseq 本身不做「出发株 vs 编辑株」差分。产品层金标准：

```text
gdtools SUBTRACT edited.gd starter.gd
```

再减去 `intended`（若有），与 `unintended.tsv` 变异集合一致（允许本文件白名单内的规范化差异）。

合成夹具 `synth_parent_child`：出发株相对骨架的历史 SNP **不得**出现在非预期列表。

引擎单样本仍须做时间 / 资源对照；产品层差分本身可不与 breseq 比 wall（breseq 无等价流水线）。

## 版本钉死

- 不因文献复现改钉 breseq **0.35.4**。
- 比对器版本与 oracle 一致；否则 RA/JC 易对不齐。
- 发版前记录 `breseq --version` 与 `bowtie2 --version` 哈希到 `testdata/VERSIONS.txt`。

## 不纳入 oracle / testdata

- Widney/Copley 2024 SRA（PRJNA1088182）及补充表
- 第一期不做 `Population_Sample.tgz`（混群）必过
