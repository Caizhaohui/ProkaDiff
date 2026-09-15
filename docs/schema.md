# ProkaDiff 输入 / 输出契约

本文是 CLI / Genome Diff / TSV 的落地契约。实现与测试以本文件为准。

**引擎子命令（单样本，oracle / 对拍）：**

```text
prokadiff evidence --ref genome.fa --fastq R1.fastq --fastq R2.fastq --threads 8 --outdir out/
```

写出 `out/output.gd`。产品层分类走无子命令的入口（两株 FASTQ 必填）。

## CLI（第一期）

```text
prokadiff \
  --starter starter_R1.fastq.gz --starter starter_R2.fastq.gz \
  --edited  edited_R1.fastq.gz  --edited  edited_R2.fastq.gz \
  --ref     NC_000913.3.gbk \
  --ref     plasmid.gbk \
  --intended intended.tsv \
  --editor cas9 \
  --spacer  NNNNNNNNNNNNNNNNNNNN \
  --pam     NGG \
  --threads 8 \
  --outdir  results/prokadiff
```

| 选项 | 必填 | 说明 |
| --- | --- | --- |
| `--starter` | 是 | 可重复。同一选项连续两个文件 = 一对 PE（先 R1 后 R2）；单个文件 = SE；可多 lane（多对）。 |
| `--edited` | 是 | 同上。缺任一株 → 非零退出，错误信息写明必须提供出发株 WGS。 |
| `--ref` | 是 | 可重复。染色体 + 质粒 + donor 骨架进入同一套参考。所有输入中的 contig ID 必须全局唯一；同一文件内或不同文件间重名会在写合并 FASTA 或调用 Bowtie2 前以非零错误拒绝，并显示重复 ID 与两个冲突来源路径。 |
| `--intended` | 否 | 目的编辑表。缺省则差分后全部进非预期。坐标相对骨架参考。 |
| `--editor` | 是 | `cas9` \| `cas12a` \| `dsb`。`cast` / `is110` → 明确错误并指向路线图。 |
| `--spacer` | cas9/cas12a 是 | 引导 RNA 的 DNA 字母（T 而非 U）。`dsb` 不要求。 |
| `--pam` | 否 | 默认：cas9 → `NGG`；cas12a → `TTTV`。可覆盖。近同源扫描允许最多 **4** 个 spacer 错配；突变距最近位点 ≤ **50 bp** 标 `near_homolog`。 |
| `--threads` | 否 | 传给 Bowtie2 / 并行 pileup。 |
| `--no-hypothesis` | 否 | 输出不含假说列。 |
| `--outdir` | 是 | 结果目录。写出 `starter.gd` / `edited.gd` / `unintended.tsv` / `summary.txt`。 |

缺出发株时 **不得** 静默把相对 NCBI 的全部 SNP 当非预期。

## `intended.tsv`（用户提供，可选）

| 列 | 含义 |
| --- | --- |
| `seq_id` | 参考序列名（与 GBK/FASTA 一致） |
| `start` | 1-based 闭区间起点 |
| `end` | 1-based 闭区间终点 |
| `ref` | 出发等位（SNP/短 indel）；结构编辑可空 |
| `alt` | 目的等位；大片段插入写 `INS` + 另附 FASTA id |
| `kind` | `snp` / `indel` / `del` / `ins` / `cassette` |

掩码规则：`M_edited` 中与某一行坐标重叠且等位基因匹配（或 JC 落在 cassette 预期连接上）的条目从非预期列表去掉，改记入「目的编辑已观察到」。

实现（`prokadiff-classify`）：

- `snp` / `ins`：坐标重叠且 `alt` 匹配（`alt` 为 `.` 或空则只按坐标）。SNP 的 GD `ref` 不另加校验。
- `del` / `indel`：仅匹配 GD `DEL`，且观测大小必须等于声明的 `end - start + 1`；GD 的 `DEL` 不带 alt，因此该行的 `alt` 必须为空或 `.`（否则不掩码）。不校验 `ref` / 具体序列。
- `cassette`：可匹配落在声明区间上的 `JC` / `MOB` / `INS` / `CON`。
- 其他 `kind` 不掩码任何条目（保守起见）。

## 引擎中间：GD 兼容子集

每个样本写出 `starter.gd` / `edited.gd`。字段与 breseq Genome Diff 对齐，至少支持：

| GD 类型 | 含义 | 产品层用途 |
| --- | --- | --- |
| `SNP` | 单碱基 | RA；差分后进 (1) 或 (2) |
| `INS` / `DEL` | 短 indel（第一期对标 breseq：RA 侧重 ≤2 bp；更长走 MC/JC）。写出 GD 前将 INS/DEL **3′ 对齐**（tandem repeat 上滑动并旋转插入序列），与 breseq / `gdtools SUBTRACT` 键一致 | 同 SNP |
| `MOB` | 移动元件插入 | class (3) |
| `JC` / `UN` | 新连接 / 未知区 | (3) 或覆盖缺口；UN 不对拍失败 |
| `AMP` / `CON` | 扩增 / 替换（若引擎能报则保留） | class (3) |

对拍命令（单样本引擎）：

```text
gdtools SUBTRACT rust.gd breseq.gd   # 应为空（无多报）
gdtools SUBTRACT breseq.gd rust.gd   # 应为空（无漏报）
```

产品层金标准（双样本）：

```text
gdtools SUBTRACT edited.gd starter.gd
```

再减去 intended（若有），与 `unintended.tsv` 的变异集合一致（允许 [parity.md](parity.md) 中列出的规范化差异）。

集合差语义（`prokadiff-gd`）：匹配键为 **坐标 + GD 类型 + 等位基因**。为处理出发株与编辑株独立比对聚类时的代表位点微漂移，`subtract` 对 JC、MOB 以及结构缺失 DEL（`size > 2`，起点与终点边界双端差 ≤5 bp）实施 ±5 bp 规范化容差消去（双端 contig、链向一致且坐标差 ≤5 bp 即视为同一物理事件减除），彻底避免先祖结构变异假阳性冒出；提供 `subtract_exact` 供严格文本精确比对。非 JC/MOB/DEL 变异与短 indel（≤2 bp）保持精确字段匹配。


## `unintended.tsv`（产品输出）

列顺序锁定（`--no-hypothesis` 时省略最后一列 `hypothesis`）：

```text
seq_id	position	end	gd_type	ref	alt	class	editor	pam_profile	offtarget_mismatch	distance_to_site	side2_seq_id	side2_position	hypothesis
```

| 列 | 必填 | 说明 |
| --- | --- | --- |
| `seq_id` | 是 | 骨架参考序列（JC 为 side 1） |
| `position` | 是 | 1-based（JC 为 side 1） |
| `end` | 是 | 闭区间（SNP/INS/JC 与 `position` 相同；DEL 为起点 + size − 1） |
| `gd_type` | 是 | `SNP` / `INS` / `DEL` / `MOB` / `JC` / … |
| `ref` | SNP 是 | SNP：参考该 1-based 位点的碱基。INS / DEL / MOB / JC 第一期为 `.` |
| `alt` | SNP/INS 是 | 新等位；DEL / JC 为 `.` |
| `class` | 是 | `structural` \| `near_homolog` \| `scattered_snv`（对应 (3)(1)(2)） |
| `editor` | 是 | CLI 传入值 |
| `pam_profile` | class (1) 是 | 如 `NGG` / `TTTV`；`dsb` 为空 |
| `offtarget_mismatch` | class (1) 建议 | 相对 spacer 的错配数 |
| `distance_to_site` | class (1) 是 | 距预测近同源位点 bp；默认阈值 50 |
| `side2_seq_id` | JC 是 | JC 第二端序列名；非 JC 行为空 |
| `side2_position` | JC 是 | JC 第二端 1-based 坐标；非 JC 行为空 |
| `hypothesis` | 否 | `sos_widney2014` / `unknown_global` / 空；`--no-hypothesis` 时省略整列 |
| `coverage` / `frequency` | 建议（后续） | consensus 下频率应接近 1；第一期 TSV **不写这两列** |

另写 `summary.txt`（TSV 风格 `key\tvalue`，一行一项）：

| 键 | 未提供 `--intended` | 提供了 `--intended` | 说明 |
| --- | --- | --- | --- |
| `editor` | CLI 值 | 同左 | 编辑器类型 |
| `intended_provided` | `no` | `yes` | 是否提供预期编辑表 |
| `intended_edits_declared` | `NA` | 声明行数（`intended.len()`） | FIX-015: 编辑级声明数 |
| `intended_edits_complete` | `NA` | 匹配到的预期编辑行数 | FIX-015: 完整匹配的编辑数 |
| `intended_edits_partial` | `NA` | 部分匹配的预期编辑行数 | FIX-015: 部分匹配的编辑数 |
| `intended_edits_missing` | `NA` | 未匹配到的预期编辑行数 | FIX-015: 缺失的预期编辑数 |
| `intended_events_observed` | `NA` | 匹配到的实际突变条目数 | FIX-015: 观测到的突变事件总目数 |
| `intended_declared` | `NA` | 同 `intended_edits_declared` | （已弃用，向后兼容） |
| `intended_observed` | `NA` | 掩码命中的突变条数 | （已弃用，向后兼容） |
| `intended_status` | `NA` | `all_observed` / `partial` / `none_observed` | FIX-015: 基于编辑级状态判定 |
| `intended_missing` | `NA` | 声明行数 − 观测突变数 | （已弃用，向后兼容） |
| `structural` / `near_homolog` / `scattered_snv` | 非预期三类计数 | 同左 | 各类别非预期突变数 |
| `starter_vs_ref_mutations` | 出发株相对骨架的变异总数（质控，不列入非预期） | 同左 | 亲本背景变异计数 |
| `offtarget_search_validation_status` | `validated_via_self_test` | 同左 | FIX-018 验证状态元数据 |
| `cfd_validation_status` | `disabled` | 同左 | FIX-018 验证状态元数据 |
| `hsu_validation_status` | `experimental` | 同左 | FIX-018 验证状态元数据 |
| `bulge_validation_status` | `experimental` | 同左 | FIX-018 验证状态元数据 |

`intended_events_observed` 计的是**被掩码的突变条目**，而 `intended_edits_*` 按**声明的编辑行**为单位进行评估。一条 cassette 声明匹配两条 JC 时，在 FIX-015 语义下该编辑行被评为 `complete`，`intended_events_observed` 计为 2，不会再错误地因为 `observed > declared` 导致判定偏差。

