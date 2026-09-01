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
- 层 0 测例：`both_strand_long_deletion_emits_jc`、`softclip_onto_repeat_copy_emits_jc`、`jc_supported_cigar_gap_promotes_short_del`、`same_reads_placed_on_two_copies_fold_to_one_jc`（原 `synth_is_mob_geometry_both_strand_softclips_aggregate`，见下方折叠规则）。
- **未在登录节点重跑 Clonal FASTQ。** 不得声称 27 条 JC 已与 0.39.0 leftover 对齐。
- **Softclip 链向规范化 + best-read 接受规则（作业 2371258 收口）：** synth_is_mob 层 1 红线类型双向 leftover=0，但 ProkDiff 吐出 empty GD，breseq 报 4 JC（200→401、200→881、201→480、201→960；夹具在 201 处插入第三个 80 bp motif 拷贝，已有拷贝在 401-480、881-960）。根因二条：(1) softclip hint 的 `aligned_pos_1` / `clip_is_left` 此前按 **read 方向**记录，minus 链接合点键错位（5′ clip 键在 151 而非 200），正负链无法聚合；现统一为**参考方向**（`ref_clip_is_left = clip_on_read_left XOR minus`；left→`first_match_ref+1`，right→`last_match_ref+1`），`clip_seq` 仍保持 read 方向（`place_softclip` 双链搜索），side2 公式相应改为 `clip_is_left == (is_rc XOR minus)`。(2) `accept_junction` 的 ≥14 bp 规则此前取所有支持读段的 **min**，一条不平衡读段（如单侧 5 bp）即可杀死 junction；现按发表方法改为**最佳单读段**每侧 ≥14 bp（max over reads of min(ov1,ov2)），且**每链**有读段每侧 ≥9 bp（per-strand max-of-min），每侧最小延伸 ≥3 bp 保留为文档性检查（`.max(3)` 钳位使其恒真）。

- **synth_is_mob 夹具几何修正（作业 2371383 BAM 诊断收口）：** (a) 旧夹具 1440 bp（left400+motif80+mid400+motif80+right400），插入点在 201；wgsim 插入片段 ~500 bp（`-d 500` 默认），跨接合点片段起点只能 ≤280，接合点恒落在片段左半，只有正链 left mate 能跨过（实测 263 条接合点 softclip 读段全部为正链），双链 `accept_junction` 规则因此吐出 0 JC。引擎本身记录 softclip 无误（位置精确在 200/201，MAPQ 全部 ≥10），属夹具几何缺陷而非引擎缺陷。(b) breseq 0.40.2 对同一数据报 4 条 JC，全部单链且全部带 `reject=COVERAGE_EVENNESS_SKEW`；ProkDiff 保留已发表的双链接受规则，不为对齐 oracle 放宽。`scripts/layer2_gd_compare.py` 现把未匹配 JC 按 `reject=` 拆分单列（`*_under_jc_oracle_rejected`），oracle 拒收条数不再混入普通 unmatched。(c) 夹具改为 ~4000 bp（left1200+motif80+mid1200+motif80+right1440），插入点移到 601（0-based 600），距基因组两端均 ≥600 bp，双链跨接合点读段在几何上存在；wgsim 参数不变。

- **minus 链 softclip 键位修正 + ±5 bp 接合点聚类（作业 2371412 BAM 诊断收口，诊断作业 2371899）：** 修正后夹具下 ProkDiff 仍吐 0 JC，breseq 0.40.2 报 2 条 JC（`599 (+1) -> 2558 (-1)`、`602 (-1) -> 2483 (+1)`，frequency ~0.6，**未**被 oracle 拒收）。对 `synth_is_mob_2371412/prokdiff/aligned.bam` 的只读诊断（`scripts/inspect_bam_jc.py`，窗口 595-610，clip ≥14）实测：
  - plus 读段键位紧密聚在 599（74 条 leading S）与 602（62+1 条 trailing S）——寄存器歧义属实：motif 前 2 bp `GT` == ref 601-602，motif 末 2 bp `GC` == ref 599-600，breseq 自己的两条 JC（同一插入事件的左、右两个接合点）也分落 599 与 602（差 3 bp）。
  - **但 minus 读段根本没键在接合点附近**：68 条 minus trailing-S（比对块终于 602）被旧规则键在 first match（散布 517-580），59 条 minus leading-S（始于 599）被键在 last match（散布 631-682）。根因：SAM/BAM 对 minus 读段的 SEQ 存**参考正向**（反向互补后），clip 在参考上的左右侧与 read 链向无关；旧的 `clip_on_read_left XOR minus` 规则把 minus 候选键到匹配块的错误一端。实测佐证：这些 minus 读段的 clip 序列精确落在 motif 拷贝上（1203/1278/2483/2558 = 拷贝两端 ±2），若 SEQ 为原始方向则应落在 503-516 唯一区。
  - 修复（层 0 可测）：(a) `pileup.rs` softclip 键位改为链向无关——leading clip → first match、`clip_is_left=true`；trailing clip → last match、`clip_is_left=false`；`place_softclip` 的 `content_rc_hit = is_rc`（不再 XOR minus）。side2 落位公式在两处错误相互抵消下数值不变，故放置结果无回归。(b) `engine.rs` 用**确定性聚类**取代精确键聚合：按 (contig, side2 contig, overlap 符号) 分组，组内按 (side1, side2) 稳定排序后贪心聚类——候选与某簇**首成员**两侧坐标均差 ≤ `JC_CLUSTER_TOL_BP`（=5，与 `layer2_gd_compare.py` 的 `jc_tol_bp=5` 对齐）即入簇；`JunctionSupport` 跨簇内全体成员聚合（正负链合并），簇代表坐标 = 排序后**首成员**（最简确定选择；oracle 侧 ±5 bp 贪心匹配下簇内任选均可）。每簇接受后吐 1 条 JC。(c) JC 链字段按 breseq 约定写 `1`/`-1`，链向由几何推导（该侧序列离接合点的延伸方向：side1 = 与 clip 侧相反；side2 = 接合点在 hit 末端则为 -1），不再取代表读段的 read 链向（原来随 BAM 读序抛硬币）。已实测复现 synth 两条 JC 与 Clonal 缺失侧翼 JC（`3289961 -1 -> 3289978 1`）的 oracle 链向；反向（is_rc）放置与 INS 接合点的链向约定**未经实测**，属 best-effort。
  - 层 0 测例：`few_bp_strand_split_junctions_cluster_into_one_jc`（599/602 跨链聚成 1 条 JC）、`junctions_beyond_cluster_tolerance_do_not_merge`（10 bp 不并簇）、`minus_strand_leading_clip_keys_at_first_match`、`minus_strand_trailing_clip_keys_at_last_match`；既有 softclip/JC 测例的 minus 读段模型改为参考正向存储。

- **Clonal 2371389 头条（oracle breseq 0.40.2，新引擎）：** 134 条 MC→DEL 过报**消除**（over_red 139 → 5，仅剩 5 SNP；DEL 134→0）。JC 仍为 0（本轮跑在聚类与 minus 键位修复**之前**），27 条 oracle JC 全未匹配（其中 9 条 oracle 自带 `reject=`）。wall / RSS（同节点单次，**非正式**记录）：ProkDiff 137.23 s / 2.97 GB，breseq 2667.67 s / 1.86 GB。详见 `benchmark/results/clonal_leftover.md` 2371389 节。

- **JC 冗余折叠（作业 2372015 收口）：** 2372015 实测 ProkDiff 吐 6 条 JC、breseq 0.40.2 吐 2 条（`599 (+1) -> 2558 (-1)`、`602 (-1) -> 2483 (+1)`，jc_matched_tol=2、under=0）；4 条 over-JC 全是同一物理事件的冗余报告：
  1. 多拷贝放置——同一批 softclip 读段的 clip 同时精确命中两个 motif 拷贝（1201-1280 / 2481-2560），多出 `599 (+1) -> 1278 (-1)` 与 `602 (-1) -> 1203 (+1)`（copy1 放置；breseq 只保留 copy2 2483/2558）。
  2. 反向重复——`2560 (-1) -> 601 (+1)` 与 `2480 (+1) -> 599 (-1)` 是两条已匹配 JC 从另一侧报告的方向重复（寄存器漂移：599+2=601、2558+2=2560；602−3=599、2483−3=2480）。

  第一期折叠规则（`prokdiff-evidence::engine`，层 0 可测；先反向、后多拷贝——先折多拷贝会把 copy 侧报告残留为额外 JC）：
  1. **反向规范化**：交换两侧后在 ±`JC_CLUSTER_TOL_BP`（5 bp）内匹配（contig 相同、各侧链向相同——链向随参考位置固有，**不**随报告方向翻转；2372015 实测反向报告的两侧链向与原报告逐位置一致，仅坐标差 ±2/3 bp）。每个物理事件只吐一条：代表 = 支持读段最多者，并列取 canonical 键（forward vs 交换侧后的字典序）最小者。
  2. **多拷贝放置折叠**：`place_softclip` 给每个候选打上来源 softclip 的身份（hint 序号；一次命中多拷贝的读段每个拷贝贡献一个候选、共享同一身份）。共享一侧（同 contig/链向、±5 bp，任一方向）且支持读段身份交集 > 较小集合 50% 的已接受 JC 并为一条；代表 = 支持读段最多的拷贝（2372015 实测 breseq 保留的正是带额外 copy-anchored 读段的拷贝），并列取 canonical 键最小者。CIGAR 长 I/D 候选身份唯一，**永不**被本规则折叠。

  写出坐标为代表簇原坐标、原方向（±5 bp 寄存器 / 方向歧义在「重复区 JC 拷贝歧义」白名单内，对拍脚本自行做侧规范化）。层 0 测例：`synth_is_mob_2372015_multicopy_and_reverse_duplicates_fold`（2372015 全几何 6→2）、`distinct_junctions_apart_beyond_tol_are_not_folded`、`distinct_read_sets_same_flank_are_not_folded`、`cigar_split_junctions_are_never_multicopy_folded`；`synth_is_mob_geometry_both_strand_softclips_aggregate` 相应改为 `same_reads_placed_on_two_copies_fold_to_one_jc`（同一批读段命中两个拷贝 → 1 条 JC）。

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
