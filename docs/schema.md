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
| `--experimental-hypothesis-annotation` | 否 | 默认关闭。启用时输出机制假说列（如需纯净审计输出请勿开启）。 |
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

M2 掩码规则：先通过 M1 的 `DifferentialEvent` 集合完成亲本扣除，再逐行评估预期编辑。仅评估为预期组成事件且未在其他声明中被判为不一致或异常的差分事件从非预期列表移出；单纯坐标重叠不能掩码异常事件。

实现（`prokadiff-classify`）：

- `snp` / `ins`：坐标重叠且 `alt` 匹配（`alt` 为 `.` 或空则只按坐标）。SNP 的 GD `ref` 不另加校验。
- `del` / `indel`：仅匹配 GD `DEL`，且观测大小必须等于声明的 `end - start + 1`；GD 的 `DEL` 不带 alt，因此该行的 `alt` 必须为空或 `.`（否则不掩码）。不校验 `ref` / 具体序列。
- `cassette`：在声明断点附近评估 `JC` / `MOB` 接头；`INS` / `CON` 等附近事件保留为异常观察，不因重叠而掩码。
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

M1 产品差分由 `prokadiff-classify::build_differential_events` 统一生成规范事件集：普通变异按 **规范坐标 + GD 类型 + 等位基因** 精确消去；JC、MOB 与双方长度均 >2 bp 的 DEL 才可进入下述确定性 ±5 bp 匹配。`prokadiff-gd::GenomeDiff::subtract_exact` 保留用于严格文本对拍。非 JC/MOB/结构 DEL 变异与短 indel（≤2 bp）始终保持精确匹配。


## `unintended.tsv`（产品输出）

列顺序锁定（默认未启用 `--experimental-hypothesis-annotation` 时省略最后一列 `hypothesis`）：

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
| `hypothesis` | 否 | 仅在启用显式参数 `--experimental-hypothesis-annotation` 时输出；默认省略整列以保持严格的观察性审计 |
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

`intended_events_observed` 计的是**被掩码的不同差分事件数**，而 `intended_edits_*` 按**声明的编辑行**为单位进行评估。一条 cassette 声明匹配两条 JC 时，该编辑行被评为 `complete`，`intended_events_observed` 计为 2；同一差分事件支持两条声明时事件只计一次。

## 差分事件与稳定标识符契约（M1）

### 1. 规范化差分事件与 DV1 标识符

为实现数据链可追溯性，M1 建立基于内容寻址的稳定事件身份 `DV1_<32位小写十六进制SHA-256截断哈希>`。

```rust
pub struct EventId(String); // 校验正则: ^DV1_[0-9a-f]{32}$
```

#### CanonicalEvent 变异表示

| 变异类型 | 规范化表示与规则 | 容错与非法输入 |
| --- | --- | --- |
| `SNP` | 确切参考 contig 标识及序列 SHA-256；1-based 坐标；单字节大写 ASCII DNA/IUPAC 碱基。 | 缺失参考、坐标 0、等位长度 != 1 或字符不在 `ACGTRYSWKMBDHVN` 返回错误。 |
| `INS` | 确切 contig 标识；大写插入序列；经 3′ 右对齐滑动与序列旋转规范化；1-based 插入前碱基。 | 空序列、位置 0、缺失参考、坐标越界或字符不在 `ACGTRYSWKMBDHVN` 返回错误。 |
| `DEL` | 确切 contig 标识；正整数大小；经 3′ 右对齐滑动规范化；1-based 第一个删除碱基。 | 大小 0、缺失参考、区间越界返回错误。 |
| `SUB` | 确切 contig 标识；1-based 起点、正整数大小、大写替换序列。不拆分为 SNP/INS/DEL，不去除公共前缀后缀。 | 大小 0、空替换、缺失参考或字符不在 `ACGTRYSWKMBDHVN` 返回错误。 |
| `MOB` | 确切 contig 标识；官方首个重复靶位 1-based 坐标；repeat_name 保持精确字节；链向规范化为 `1` 或 `-1`；保留带符号重复区长度。 | 未知链向、缺失参考、坐标 0 返回错误。 |
| `AMP` | 确切 contig 标识；1-based 起点、正整数大小、正整数 `new_copy_number`。 | 非法数字、大小或拷贝数为 0、区间越界返回错误。 |
| `CON` | 确切 contig 标识；1-based 目标起点、正整数大小、精确 UTF-8 `region`。 | 空 region、缺失参考、区间越界返回错误。 |
| `INV` | 确切 contig 标识；1-based 起点、正整数大小。反向互补由参考序列与区间唯一定义，不另存序列。 | 非法数字、大小 0、缺失参考、区间越界返回错误。 |
| `JC` | 链向规范化为 `1` 或 `-1`；两侧 `(contig digest, seq_id, position, strand)` 按字典序排序使得双向报告一致；保留带符号 overlap。 | 未知链向、坐标 0、缺失参考返回错误。 |

- **参考序列身份**：Contig 序列规范化大写后计算 32 字节 SHA-256，避免文件换行、注释差异影响。
- **等位字母表**：SNP、SUB、INS 仅接受大小写不敏感的 DNA/IUPAC `ACGTRYSWKMBDHVN`，进入规范事件前统一转为大写。
- **真实参考要求**：`classify()` 必须接收调用所涉及的真实 `RefContig`；缺失 contig 返回 `DifferentialError`，不得按事件坐标构造全 `N` 伪参考。所有坐标和区间终点使用 checked arithmetic 验证后才进入规范化。
- **3′ 对齐共享**：`INS` 和 `DEL` 在 `prokadiff-evidence` 与 `prokadiff-classify` 间共享同一种参考序列 3′ 右对齐实现。

#### DV1 二进制编码与哈希 Preimage

- **Preimage 格式**：
  1. 域前缀固定字节：`ProkaDiff\0DV1\0`。
  2. 单字节类型码：SNP `01`, SUB `02`, INS `03`, DEL `04`, MOB `05`, AMP `06`, CON `07`, INV `08`, JC `09`。
  3. Contig 编码：大端 `u32` 长度前缀 UTF-8 `seq_id` + 32 字节序列 SHA-256。
  4. 数值字段：8 字节大端整数（`u64` / `i64`）。
  5. 序列/名称字段：大端 `u32` 长度前缀字节。
  6. SNP 碱基：单字节 ASCII（无长度前缀）。
  7. 链向：单字节 `01`（`+1`）或 `02`（`-1`）。
  8. JC：按字典序先后写入规范化 side 1、side 2，最后写入 8 字节大端 `overlap`。
- **哈希算法与截断**：
  对 Preimage 计算 SHA-256，取前 16 字节（128 bits），小写十六进制编码为 32 个字符，拼接 `DV1_` 前缀。
- **排除字段**：
  原始 GD 数字 ID、`parent_ids`、证据记录、`pd_*` 属性、合并源 ID、样本名、文件路径、编辑器类型、输出顺序及无关事件严格排除在 Preimage 之外。
- **碰撞防护**：
  内存构建期使用 `BTreeMap<EventId, Vec<u8>>` 记录 Preimage，如遇不同 Preimage 产生相同 128 位截断哈希，立即返回 `DifferentialError::EventIdCollision`，严禁静默加盐或重编号。

### 2. 亲本扣除与结构容差

- **执行顺序**：
  1. 独立解析并规范化编辑株与出发株记录；
  2. 样本内合并完全相同的规范化重复事件（exact duplicates）；
  3. 先行一对一扣除规范化完全匹配（exact match）的亲本事件；
  4. 随后在允许的结构变异类型间按确定性排序规则执行一对一 ±5 bp 容差匹配消去；
  5. 生成仅属于编辑株的 `DifferentialEvent` 集合。
- **容差适用类型限定**：
  - 容差匹配规则见 [parity.md](parity.md#双样本差分亲本容差消去)。
  - **仅允许**以下三类事件适用 ±5 bp 坐标容差消去：
    - `JC`：双端 contig、链向及带符号 `overlap` 一致，且两侧坐标差均 ≤ 5 bp；
    - `MOB`：靶位 contig、repeat_name、链向及 `duplication_size` 一致，且靶位坐标差 ≤ 5 bp；
    - 结构型 `DEL`：出发株与编辑株两方删除长度均 > 2 bp，同一 contig，且起点及终点坐标差均 ≤ 5 bp。
  - **严禁**将 ±5 bp 容差扩展至 SNP、SUB、INS、短 DEL（≤2 bp）、AMP、CON、INV。
- **确定性消去原则**：
  在 exact canonical subtraction 后，对剩余候选执行确定性的最小费用最大二分匹配：先最大化有效一对一匹配数，再最小化累计 `(max_coordinate_delta, sum_coordinate_delta)`；同费用路径按编辑株和出发株 canonical bytes 顺序稳定裁决。输入行顺序不参与结果。容差仅用于消除继承背景变异，绝不修改编辑株事件坐标或哈希。

### 3. 重复物理事件合并

- 仅在规范化后的 DV1 Preimage 字节完全相同时才进行去重合并，避免因 ±5 bp 近邻非传递性导致不同物理事件被错误坍塌。
- 选代表项（`representative: GdEntry`）采用确定性全序元组：
  `(kind_tag, normalized_core_fields, raw_core_fields, raw_id, sorted_parent_ids, sorted_attribute_key_value_pairs)`。
- 所有来源的原始 GD ID 归入 `merged_source_ids: Vec<u32>`（排序且去重）。
- MOB 的组装 JC 仅在引擎显式通过 `parent_ids` 链接时吸收为 MOB 的证据，不进行坐标推测性合并；该吸收在出发株和编辑株 GD 上对称执行。

### 4. 证据引用与真实性法则（No-Inference Rule）

```rust
pub struct EvidenceReferences {
    pub gd_parent_ids: Option<Vec<u32>>,
    pub ra: Option<Vec<GdEvidenceRef>>,
    pub mc: Option<Vec<GdEvidenceRef>>,
    pub jc: Option<Vec<GdEvidenceRef>>,
    pub pd_attributes: std::collections::BTreeMap<String, Vec<String>>,
    pub unresolved_parent_ids: Vec<u32>,
}
```

- **禁止变异类型推断证据**：彻底废除“SNP 必然意味着 RA”或“大 DEL 必然意味着 MC”的硬编码推断。未观测到实际证据时，引用字段必须为 `None`，在输出边界格式化为字面值 `NA`。
- **证据记录解析**：仅将 `parent_ids` 中实际指向的 RA/MC/JC 记录解析为对应类型的 `GdEvidenceRef`；悬空 ID 保留在 `unresolved_parent_ids`。
- **兼容层 `EvidenceSummary`**：`ra` / `mc` / `jc` 字段由 `bool` 改为 `Option<bool>`。只有观测到对应证据记录时才为 `Some(true)`，无证据时为 `None`。`format_brief()` 输出格式严格为 `RA=1;MC=NA;JC=NA`（未观测写 `NA`）。
- **审计传递**：`AuditResult` 和报告兼容层只能从对应 `DifferentialEvent.evidence` 生成 `EvidenceSummary`；不得从 `GdEntry.kind` 重建或推断证据。

### 5. 公开兼容性边界

- M1 维持现有 7 个公开 TSV（`unintended.tsv`、`edit_outcomes.tsv`、`post_edit_variants.tsv`、`offtarget_sites.tsv`、`mutation_offtarget_links.tsv`、`provenance.tsv`）及 `summary.txt` 的列名、列序与现有行结构完全不变。
- M1 不在公开 TSV 中增加 `event_id` 列，现有数字 ID、`mut_<id>` 和 `VAR_####` 保留由 `representative` 映射兼容。
- `DifferentialEvent` 在 M1 作为 Rust 内存领域模型完成闭环，跨文件统一公开暴露出入将在 M2–M4 逐阶段演进。

## M2 预期编辑验证契约

- `assess_intended_edits` 直接接收亲本扣除后的 `DifferentialEvent`。`IntendedEditAssessment` 内部以 `EventId` 引用差分事件，并用 `IntendedEventRelationship` 区分预期组成、部分/不一致观察和靶位异常。一个编辑可关联多个事件，一个事件可参加多个编辑声明的检查；关联不复制差分事件或证据。
- 验证先于掩码。只有 `ExpectedConstituent` 且不在任何声明中标为部分或异常的事件被移出 `unintended.tsv`。异常事件仍可通过 `EventId` 回指原 `DifferentialEvent.evidence`。等价候选按坐标/大小偏差及稳定事件身份裁决，GD 行顺序不参与身份或决胜。
- MC 是独立证据诊断，不产生 `EventId`，也不能仅凭靶位附近的 MC 判定编辑株新发结构突变。已观测且非同记录亲本 MC 的诊断记为 `Observed`；未有可证明的诊断时为 `Unknown`，不推断证据缺失或生物学缺失。RA/MC/JC 证据仍只来自 `DifferentialEvent.evidence`，缺失输出 `NA`。
- M2 保持所有既有公开 TSV 列名、列序及 `summary.txt` 键不变。`edit_outcomes.tsv` 的 `matched_event_ids` / `unexpected_events` 在写出边界把内部 `EventId` 映射为代表 GD 数字 ID，维持现有单元格格式；M2 不新增公开 DV1 列。公开 ID 全面协调留待 M4。
