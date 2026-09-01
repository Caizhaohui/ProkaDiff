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
| RA 短 indel | 引擎在写出 GD 前将 INS/DEL **3′ 对齐**（tandem repeat 上滑动并旋转插入序列），以便与 breseq / `gdtools SUBTRACT` 键一致。若仍有残留 ±1，视为实现缺陷而非白名单。 |
| consensus 频率舍入 | `frequency` / 相关浮点字段在 consensus 下接近 1 即可，允许舍入差 |

## 不允许的差异（结果红线）

下列任一出现即结果失败：

| 类别 | 说明 |
| --- | --- |
| SNP 位点 / 等位基因不一致 | 同一位点 ref/alt 与 oracle 不符 |
| unique-mapping MOB 漏报 | 明确可 unique 映射的移动元件插入漏报 |
| 大 DEL 漏报 | 明确的大片段缺失漏报（MC/JC 路径） |
| 白名单外的系统性多报 | `SUBTRACT rust.gd breseq.gd` 留下非白名单条目 |

### 第一期 MC→DEL（作业 2369858 收口）

Clonal leftover（breseq **0.39.0**，见 `benchmark/results/clonal_leftover.md`）：134 条 over DEL 均伴 MC，与 breseq 0 精确交集。原因：unique-mapping 深度 0 在重复区 / IS 拷贝上也会出现（低 MAPQ 读段仍覆盖），不等于缺失。

第一期规则（层 0 可测，写入 `prokdiff-evidence::mc`）：

1. 覆盖缺口**始终写 MC**（contig 末端仍不写 DEL / 不把末端 0-depth 当缺失）。
2. **不**再把每个内部 MC 自动升格为 consensus DEL。
3. 升格 DEL 仅当 unique-0 且 **total（任意 MAPQ）深度低于 cutoff** 的核心满足：
   - 两端有 unique 覆盖（允许 1 bp pileup end-trim 边距）；
   - 无 JC：核心长度 ≥ 50 bp（`MC_DEL_MIN_LEN`）；
   - 或 JC 端点落在两侧翼 ±20 bp 内：长度 ≥ 3 bp。
4. 大 DEL（如 REL606 `3894997` / 6934 bp IS 介导）应走 **两侧 JC 或 MC+JC**，不靠「unique 缺口一律 DEL」。

测试：`repeat_like_unique_gap_is_not_promoted_to_del`、`short_internal_true_gap_is_not_promoted_to_del`、`contig_terminal_zero_depth_is_not_promoted_to_del`、`long_internal_zero_total_gap_promotes_to_del`、`jc_supported_short_gap_promotes_to_del`；引擎层 `short_internal_coverage_gap_does_not_emit_del`、`repeat_like_low_mapq_coverage_is_not_del`、`long_internal_true_gap_emits_del`。

### 第一期 JC / MOB

作业 2369858：ProkDiff JC=0、MOB=0；breseq 27 JC、5 unique-mapping MOB。根因（清洁室，非 breseq 源码）：(1) 正负链 JC 被拆进不同聚合键，`accept_junction` 永远缺一链；(2) CIGAR softclip 未对照参考放置第二端。

第一期实现：

- 长 CIGAR I/D（>2 bp）仍作 mosaic 候选；**聚合键不含 read strand**，正负链合并后再 `accept_junction`。
- ≥14 bp softclip：在参考上精确匹配（含反向互补），跳过本比对的平凡延续；重复拷贝允许 unassigned JC。
- **MOB：** 第一期不写 `IS150` 等注释名（无 GBK repeat 表）。Clonal 的 5 条 unique-mapping MOB **仍可能只以 JC 出现**；等 breseq **0.40.x** 层 2 再验。这落在上表「重复区 JC 拷贝歧义」白名单，**不**把 breseq `.gd` 抄进运行时。
- 层 0 测例：`both_strand_long_deletion_emits_jc`、`softclip_onto_repeat_copy_emits_jc`、`jc_supported_cigar_gap_promotes_short_del`、`synth_is_mob_geometry_both_strand_softclips_aggregate`。
- **未在登录节点重跑 Clonal FASTQ。** 不得声称 27 条 JC 已与 0.39.0 leftover 对齐。
- **Softclip 链向规范化 + best-read 接受规则（作业 2371258 收口）：** synth_is_mob 层 1 红线类型双向 leftover=0，但 ProkDiff 吐出 empty GD，breseq 报 4 JC（200→401、200→881、201→480、201→960；夹具在 201 处插入第三个 80 bp motif 拷贝，已有拷贝在 401-480、881-960）。根因二条：(1) softclip hint 的 `aligned_pos_1` / `clip_is_left` 此前按 **read 方向**记录，minus 链接合点键错位（5′ clip 键在 151 而非 200），正负链无法聚合；现统一为**参考方向**（`ref_clip_is_left = clip_on_read_left XOR minus`；left→`first_match_ref+1`，right→`last_match_ref+1`），`clip_seq` 仍保持 read 方向（`place_softclip` 双链搜索），side2 公式相应改为 `clip_is_left == (is_rc XOR minus)`。(2) `accept_junction` 的 ≥14 bp 规则此前取所有支持读段的 **min**，一条不平衡读段（如单侧 5 bp）即可杀死 junction；现按发表方法改为**最佳单读段**每侧 ≥14 bp（max over reads of min(ov1,ov2)），且**每链**有读段每侧 ≥9 bp（per-strand max-of-min），每侧最小延伸 ≥3 bp 保留为文档性检查（`.max(3)` 钳位使其恒真）。

层 2 等 0.40.x 再重跑；不要重提作业 2369858。

时间 / RSS 无「必须更快」硬阈值；但 **必须有可复核记录**。若 ProkDiff 明显更慢或吃内存异常，记入 notes，不得静默省略。

## 产品层（双样本）

breseq 本身不做「出发株 vs 编辑株」差分。产品层金标准：

```text
gdtools SUBTRACT edited.gd starter.gd
```

再减去 `intended`（若有），与 `unintended.tsv` 变异集合一致（允许本文件白名单内的规范化差异）。

合成夹具 `synth_parent_child`：出发株相对骨架的历史 SNP **不得**出现在非预期列表。

**产品层 subtract 与引擎 JC 白名单的张力：** 单样本引擎允许 junction 坐标 ± 数 bp 的规范化（见上表「重复区 JC 拷贝歧义」）。产品层 `M_edited \ M_parent` 目前与 `gdtools SUBTRACT` 一样用**精确字段键**（`prokdiff-gd::GdEntry::subtract_key`）。因此出发株与编辑株两次独立调用若对同一 JC 报出差 1 bp 的坐标，该 JC 会留在非预期里并标 `structural`。这不是引擎漏报，是差分语义；第一期不把模糊匹配写进运行时。层 3 配对数据若实测到这类假阳性，在本表记下容差后再改键，禁止未测先改。MC 衍生的大 DEL 边界同此风险。

`summary.txt` 报 `intended_provided` / `intended_declared` / `intended_observed` / `intended_status`（`all_observed` \| `partial` \| `none_observed` \| `NA`）以及 QC 行 `intended_missing`。字段契约见 [docs/schema.md](schema.md)。不做逐条长篇报告。

引擎单样本仍须做时间 / 资源对照；产品层差分本身可不与 breseq 比 wall（breseq 无等价流水线）。

## 版本钉死

- 不因文献复现改钉 breseq **0.35.4**。
- 比对器版本与 oracle 一致；否则 RA/JC 易对不齐。
- 发版前记录 `breseq --version` 与 `bowtie2 --version` 哈希到 `testdata/VERSIONS.txt`。

## 不纳入 oracle / testdata

- Widney/Copley 2024 SRA（PRJNA1088182）及补充表
- 第一期不做 `Population_Sample.tgz`（混群）必过
