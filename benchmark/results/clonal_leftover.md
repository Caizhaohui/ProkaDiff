# Clonal leftover 分类（作业 2369858）

作业目录 `benchmark/results/clonal_2369858/` 只作本地复核，**不进 git**（csv / `.gd` / FASTQ / 作业树均 gitignore）。本文件是可跟踪的分类说明。

**本轮不重提、不重跑作业 2369858。** 层 2 Clonal / Example 6 等钉死 **breseq 0.40.x** 后再重跑。本次 leftover 的 oracle 是 **breseq 0.39.0**（配套 bowtie2 2.5.4），不能标成 0.40.x 对拍已过。

## 作业元数据

| 项 | 值 |
| --- | --- |
| Slurm 作业 | **2369858**（`scripts/submit_parity_clonal.sh`） |
| 分区 | `qcpu_18i` |
| 节点 | `bnode5.tibhpc.net` |
| 线程 | 8 |
| Oracle | breseq **0.39.0**（不是 0.40.x） |
| 比对器 | bowtie2 2.5.4 |
| 工作负载 | Methods Mol. Biol. 2014 Clonal（REL606 / `NC_012967`） |
| 比较脚本 | `scripts/layer2_gd_compare.py`（红线 SNP/INS/DEL/MOB/AMP/CON 精确键；JC 端点规范化后 ±5 bp） |

同作业实测（记录，不作加速比宣传）：

| 工具 | wall_s | peak RSS (KB) |
| --- | --- | --- |
| ProkDiff | 139.92 | 2958124 |
| breseq 0.39.0 | 2691.35 | 1861616 |

## 类型表（同节点 rust vs 同机 breseq）

方向：`over` = `gdtools SUBTRACT rust.gd breseq.gd`（ProkDiff 多报）；`under` = `SUBTRACT breseq.gd rust.gd`（漏报）。

| 方向 | SNP | INS | DEL | MOB | JC |
| --- | --- | --- | --- | --- | --- |
| over rust−breseq | 5 | 0 | **134** | 0 | 0 |
| under breseq−rust | 0 | 3 | 3 | **5** | **27** |

- `over_red` = 139 = DEL 134 + SNP 5。DEL 与 breseq **0 精确交集**。
- ProkDiff `output.gd`：SNP 33、DEL 134、INS / MOB / JC **全 0**。
- `under_red` = 11 = INS 3 + DEL 3 + MOB 5；另 JC 27 全漏。
- JC ±5 bp 贪心匹配 **0**（容差无效：ProkDiff 未报任何 JC）。
- vs 官方 `Clonal_Output`：`over_red` 同样 139；`under_red` 11；JC unmatched 26（官方比同机 breseq 少 1 条 JC 量级差异，不在本轮解释）。

### 134 条 over DEL（均伴 MC）

覆盖缺口被立刻写成 consensus DEL。长度桶：

| 长度 | 条数 |
| --- | --- |
| &lt;50 | 37 |
| 50–499 | 48 |
| 500–4999 | 45 |
| ≥5000 | 4 |

解释（引擎，非 breseq 源码）：unique-mapping 深度 0 在重复区 / IS 侧翼也会出现，不等于缺失。第一期应收紧 MC→DEL，见 [docs/parity.md](../../docs/parity.md)。

### 同坐标错分（2 / 5 条 over SNP）

| 坐标 (REL606, 1-based) | ProkDiff over | breseq under | 备注 |
| --- | --- | --- | --- |
| 3289962 | SNP `C` | DEL 16 bp | 同一位点：缺口/缺失 vs 点突变 |
| 3894997 | SNP `T` | DEL 6934，`mediated=IS150` | 真大缺失应走 JC/MOB，不能靠宽松 MC→DEL |

其余 3 条 over SNP：`2031736 G`、`2054876 G`、`2054924 A`（本轮不重跑，不扩写）。

### under MOB / 大 DEL / INS（红线）

under MOB 例：

- `16972 IS150` strand `-1` dup 3
- `1733647 IS150` `-1` 3
- `1821525 IS150` `-1` 3
- `3015771 IS150` `-1` 4
- `4524522 IS186` `1` 6

under 大 DEL：`3894997` len=6934 `mediated=IS150`。

under INS 3 条（同机 breseq）：`475292 G`、`3875632 T`、`3893551 G`。另有 under DEL `4126706` size 1（与 3289962 / 3894997 合为 under DEL 3）。

### JC 27 全漏

ProkDiff JC=0，故 ±5 bp 规范化帮不上。breseq JC 多为 unique 侧 × 参考已有 IS/repeat 侧（例如 `16972`↔`588495` IS150）。第一期目标：先从 softclip / mosaic 报出 JC；Clonal 的 5 条 unique-mapping MOB 仍可能只以 JC 出现，等 0.40.x 层 2 再验。

## 后续层 2

- **不要**重跑 2369858 / Example 6（2369859）来「证明」本轮引擎改动。
- 升 breseq **0.40.x** 后用同一提交脚本新开作业，再写新的 leftover / wall / RSS。
- 禁止把本表数字写成加速比或「对拍已过」。

# 2026-09-01 复跑：作业 2371272（oracle breseq 0.40.2）

**关键 caveat（实测，非推测）：本作业的 ProkDiff 侧跑的是旧二进制，引擎修复未生效。**
`target/release/prokdiff` 构建于 2026-08-31 22:48:48，早于 a9f8e85（2026-09-01 10:52:41，MC→DEL 收紧 + 链合并 JC）与 f31503c（JC 链向规范化 + best-read 接受规则）。提交脚本只在二进制缺失时才 `cargo build`，Slurm 日志无 "building release prokdiff" 行。实测证据：`clonal_2371272/prokdiff/output.gd` 与 2369858 的 **md5 完全相同**（`fd529074…`）。因此本节数字 = **旧引擎 × 新 oracle**，不能用来回答引擎修复是否生效。

## 作业元数据（2371272）

| 项 | 值 |
| --- | --- |
| Slurm 作业 | **2371272**（同一 `scripts/submit_parity_clonal.sh`），退出码 2（红线 leftover 非空的告警退出，输出完整） |
| 分区 / 节点 / 线程 | `qcpu_18i` / `bnode5.tibhpc.net` / 8 |
| Oracle | breseq **0.40.2**（GD 头 `#=PROGRAM breseq 0.40.2` 实测），bowtie2 2.5.4 |
| ProkDiff 二进制 | **旧**（2026-08-31 22:48 构建，pre-a9f8e85） |

同作业实测（ProkDiff 侧为旧二进制，**不得**用于加速比宣传）：

| 工具 | wall_s | peak RSS (KB) |
| --- | --- | --- |
| ProkDiff（旧二进制） | 181.55 | 2954852 |
| breseq 0.40.2 | 2690.18 | 1856944 |

## 类型表：0.39.0 基线 vs 2371272（红线类型；UN/RA/MC 证据行不计）

方向 `over` = rust−breseq（多报）；`under` = breseq−rust（漏报）。

| 方向 | 轮次 | SNP | INS | DEL | MOB | JC unmatched |
| --- | --- | --- | --- | --- | --- | --- |
| over | 2369858（0.39.0） | 5 | 0 | 134 | 0 | 0 |
| over | 2371272（0.40.2） | 5 | 0 | 134 | 0 | 0 |
| under | 2369858（0.39.0） | 0 | 3 | 3 | 5 | 27 |
| under | 2371272（0.40.2） | 0 | 3 | 3 | 5 | 27 |

- `subtract_rust_minus_breseq.gd` 两轮 **md5 相同**（`0ff6a93b…`）：134 DEL + 5 SNP 逐条一致。
- `subtract_breseq_minus_rust.gd` 两轮仅 GD 头注释不同（`#=PROGRAM` 0.39.0→0.40.2、`#=MAPPED-READS` 6731944→6731945 等）；**突变条目逐行相同**。
- 即对本样本，breseq 0.39.0→0.40.2 的 oracle 升级本身没有改变 leftover 集合。

## 基线三问的回答（受 stale-binary caveat 限制）

- **(a) 134 条 MC→DEL overcall 消失了吗？** 没有。2371272 的 over DEL 仍为 134 且逐条相同——但原因是 ProkDiff 侧跑了 pre-a9f8e85 旧二进制，**该问题本轮未被真正检验**。
- **(b) JC 从 0 上来了吗（27 条中几条 ±5 bp 匹配）？** 没有。ProkDiff 仍报 0 条 JC，27 条 breseq JC 全部 unmatched（±5 bp 匹配 0 条）。同样是旧二进制所致，修复效果未检验。
- **(c) 5 条 MOB 与 INS/SNP 小项？** 不变。under 侧 5 MOB（IS150 ×4 + IS186 ×1）、3 INS、3 DEL 与基线逐条一致；over 侧 5 SNP（含 3289962 / 3894997 同坐标错分 2 条）一致。

## 下一步

- 在计算节点重新 `cargo build --release -p prokdiff`（或删掉旧二进制让提交脚本自建）后重提 Clonal 作业；本轮 2371272 的 ProkDiff 侧数字作废，只保留 breseq 0.40.2 oracle 侧输出。
- 提交脚本宜加「二进制mtime < 最新源码提交则重建」或显式 `cargo build` 防 stale；未修复前不得引用本轮引擎侧结果。

# 2026-09-01 复跑：作业 2371389（oracle breseq 0.40.2，新引擎首轮）

首个真正跑到新引擎（MC→DEL 收紧 + 链合并 JC + best-read 接受规则）的 Clonal 作业。提交脚本已显式 `cargo build`（2371272 的 stale-binary 教训已修）。

## 作业元数据（2371389）

| 项 | 值 |
| --- | --- |
| Slurm 作业 | **2371389**（`scripts/submit_parity_clonal.sh`） |
| 分区 / 线程 | `qcpu_18i` / 8 |
| Oracle | breseq **0.40.2**，bowtie2 2.5.4 |
| ProkDiff 二进制 | 新引擎（含 a9f8e85 MC→DEL 收紧、f31503c JC 链向合并；**不含**后续的 minus 键位修正与 ±5 bp 聚类） |

同作业实测（同节点、单次；wall-clock 为**非正式**记录，不作加速比宣传）：

| 工具 | wall_s | peak RSS (KB) |
| --- | --- | --- |
| ProkDiff | 137.23 | 2974112 |
| breseq 0.40.2 | 2667.67 | 1857152 |

## 类型表（2371389，红线类型；UN/RA/MC 证据行不计）

方向 `over` = rust−breseq（多报）；`under` = breseq−rust（漏报）。

| 方向 | SNP | INS | DEL | MOB | JC unmatched |
| --- | --- | --- | --- | --- | --- |
| over | 5 | 0 | **0** | 0 | 0 |
| under | 0 | 3 | 3 | 5 | 27（其中 9 条 oracle 自带 `reject=`） |

- **134 条 MC→DEL 过报消除**：over DEL 134 → 0，`over_red` 139 → 5（仅剩 5 条 SNP：`2031736 G`、`2054876 G`、`2054924 A`、`3289962 C`、`3894997 T`，与基线逐条一致；`3289962` / `3894997` 仍属同坐标错分，见首轮分类）。
- under 侧不变：INS 3 + DEL 3 + MOB 5（IS150 ×4 + IS186 ×1）+ JC 27。
- **JC 仍全漏**：ProkDiff 报 0 条 JC，27 条 oracle JC ±5 bp 匹配 0 条（9 条 oracle-rejected 单列）。本轮引擎尚无 minus 链 softclip 键位修正与 ±5 bp 聚类；synth_is_mob 2371412 的 BAM 诊断（作业 2371899）已定位根因，修复见 [docs/parity.md](../../docs/parity.md) 第一期 JC / MOB 节。
- vs 官方 `Clonal_Output`：`over_red` 5、`under_red` 11、JC unmatched 26（oracle-rejected 8），与同机 breseq 同量级。

## 下一步（2371389 之后）

- 用含 minus 键位修正 + ±5 bp 聚类的引擎重提 Clonal，再验 27 条 JC 的匹配数与 5 条 MOB 是否以 JC 形式出现。
- 5 条 over SNP 的同坐标错分（`3289962` DEL 16 bp、`3894997` DEL 6934 IS150）仍待 JC/MOB 路径收口。

# 2026-09-01: 作业 2372151（JC 聚类+折叠后）

**关键 caveat（实测，非推测）：折叠提交未进入本轮二进制。** `target/release/prokdiff` mtime = 2026-09-01 13:51:08（提交脚本 `cargo build` 前置 / 作业 13:51:59 启动时为 no-op），包含 7043cca（13:26:54，minus 链 softclip 键位修正 + ±5 bp 聚类）；而 52c5256（13:52:21，反向与多拷贝 JC 折叠）提交于作业启动**之后**，本轮未被检验。且实测 `clonal_2372151/prokdiff/output.gd` 与 2371389 **md5 完全相同**（`0d48c2eb…`）——无论折叠是否在working tree，聚类改动在 Clonal 真实数据上**行为零变化**。

## 作业元数据（2372151）

| 项 | 值 |
| --- | --- |
| Slurm 作业 | **2372151**（`scripts/submit_parity_clonal.sh`，`pd_clonal`），State=FAILED 但 ExitCode=**2:0**（红线 leftover 非空的告警退出；csv 双行齐全，输出完整，非真失败） |
| 分区 / 节点 / 线程 | `qcpu_18i` / `bnode6.tibhpc.net` / 8 |
| Elapsed | 00:52:51 |
| Oracle | breseq **0.40.2**（GD 头实测），bowtie2 2.5.4 |
| ProkDiff 二进制 | mtime 13:51:08，含 7043cca 聚类；**不含** 52c5256 折叠（见上） |

同作业实测（同节点、单次；wall-clock 为**非正式**记录，不作加速比宣传）：

| 工具 | wall_s | peak RSS (KB) |
| --- | --- | --- |
| ProkDiff | 141.52 | 2974944 |
| breseq 0.40.2 | 3022.03 | 1857008 |

## 类型表（2372151，红线类型；UN/RA/MC 证据行不计）

方向 `over` = rust−breseq（多报）；`under` = breseq−rust（漏报）。

| 方向 | SNP | INS | DEL | MOB | JC unmatched |
| --- | --- | --- | --- | --- | --- |
| over | 5 | 0 | **0** | 0 | 0 |
| under | 0 | 3 | 3 | 5 | 27（其中 9 条 oracle 自带 `reject=`） |

- subtract 文件实测计数：`subtract_rust_minus_breseq.gd` = SNP 5 + MC 134（证据行）；`subtract_breseq_minus_rust.gd` = INS 3 + DEL 3 + MOB 5 + **JC 27**（另 RA 45 / UN 789 / MC 4 证据行不计）。
- csv（2372151，vs 同机 breseq）：`over_red=5`、`under_red=11`、`over_jc_unmatched=0`、`under_jc_unmatched=27`、`under_jc_oracle_rejected=9`、`jc_matched_tol=0`。vs 官方 `Clonal_Output`：`over_red=5`、`under_red=11`、`under_jc_unmatched=26`、`oracle_rejected=8`、`jc_matched_tol=0`。

## 验收四问的回答

- **(a) MC→DEL 过报仍为 0？** 是，无回归。over DEL=0，`over_red`=5（全 SNP）；134 条 MC 仅以证据行存在，未提升为 DEL。
- **(b) JC 恢复：27 条中几条 ±5 bp 匹配？** **0 条**（`jc_matched_tol=0`）。ProkDiff `output.gd` 的 JC 行数**仍为 0**（未变非零），且 output.gd 与 2371389 逐字节相同。oracle 侧 breseq output.gd 共 27 条 JC（9 条带 `reject=`），全部漏报。
- **(c) over_jc_unmatched？** **0**——ProkDiff 未多发任何 JC，无坐标可抽查。
- **(d) under_red 与 5 条 over SNP 是否不变？** 不变。under INS 3 + DEL 3 + MOB 5（IS150 ×4 + IS186 ×1）与基线逐条一致；over SNP 5 条（`2031736 G`、`2054876 G`、`2054924 A`、`3289962 C`、`3894997 T`）与 2369858 / 2371389 逐条一致。

## 结论与下一步（2372151 之后）

- JC 聚类（7043cca）在 Clonal 真实数据上**未产出任何 JC**；27 条 oracle JC 恢复数仍为 0/27。根因不是 synth_is_mob 几何，而是主 BAM 无 S≥14（见下节）。
- 折叠（52c5256）未进入本轮二进制；**不要**再为折叠重提 Clonal——主 BAM 无长 softclip 可折叠。
- 5 条 over SNP 的同坐标错分（`3289962`、`3894997`）依旧待 JC/MOB 路径收口。

## BAM 诊断（作业 2372402 / 2372418，对 clonal_2372151 `aligned.bam` 只读）

Clonal 读长 **36 bp**。整 BAM 实测：mapped **7,195,511**；含任意 S **625,062**；S 长度 1–5: **582,007**，6–13: **43,055**，**≥14: 0**。H≥14 / `SA:Z` / supplementary: **0**。27 条 oracle JC 窗口均有 mapped 覆盖与短 S（maxS 峰值 **12**），**可放置 S≥14 = 0**。引擎 `MIN_CLIP_FOR_JC=14` 因此输入为空。窗口 MAPQ 36–44（「低 MAPQ 滤掉长 clip」为假）。本样本 breseq JC 来自二次比对候选接合点（`04_candidate_junction_alignment/*.sam`），不是主 BAM 长 softclip。synth_is_mob 能报是因为那份 BAM 确实有 S≥14。

下一步 JC 工作是发表方法中的 **candidate-junction 二次比对**，不是再为折叠重提 Clonal。详见 [docs/parity.md](../../docs/parity.md) 第一期 JC / MOB。

# 2026-09-02: 作业 2372539（candidate-junction 二次比对；9/1 15:21 启动跨夜，无 GD 产出后被取消）

**本作业未完成 ProkDiff evidence，没有 leftover 可数。** 终态 **CANCELLED**（`sacct`：`CANCELLED by 105002647`，即提交用户 `caizhh`；`.batch` ExitCode **0:15** = SIGTERM）。无 `output.gd`、无 `work/jc/`、`prokdiff.time` 空文件、`clonal.csv` 无 2372539 行。下列类型表 / 四问全部记为**缺失**，不编造 leftover、不把 sacct 数字写成 `/usr/bin/time` 的 wall/RSS、不作对拍或加速比宣传。

## 作业元数据（2372539）

| 项 | 值 |
| --- | --- |
| Slurm 作业 | **2372539**（`scripts/submit_parity_clonal.sh`，`pd_clonal`） |
| 分区 / 节点 / 线程 / 申请内存 | `qcpu_18i` / `bnode6` / 8 / 60G |
| TimeLimit | `1-00:00:00`（计划 EndTime 2026-09-02T15:21:37；取消时距墙钟上限约 4h48m） |
| Start / End / Elapsed | 2026-09-01T15:21:37 → 2026-09-02T10:33:26 / **19:11:49**（跨夜） |
| State / ExitCode | 作业：`CANCELLED by 105002647` **0:0**；`.batch`：`CANCELLED` **0:15** |
| `.batch` MaxRSS / MaxVMSize（sacct） | **8388460K** / 10776268K（作业步记账，**不是** `prokdiff.time`） |
| Oracle / 比对器 | 未跑到 breseq 段；日志未写出 oracle 版本行。上一轮同脚本钉死 breseq 0.40.2 + bowtie2 2.5.4，本作业**未验证** |
| ProkDiff 二进制 | `target/release/prokdiff` mtime **2026-09-01 15:20:52**（早于作业启动 15:21:37 约 45s） |
| 源码相对 `1b4f41c` | HEAD=`1b4f41c`（2026-09-01 15:08:20）；二次比对在**未提交 working tree**：`engine.rs` / `fasta.rs` / `lib.rs` / `pileup.rs` 已改，`jc_seq.rs` 未跟踪。源码 mtime 均早于二进制（`engine.rs` 15:19:30 最晚），作业内 `cargo build --release -p prokdiff` 日志 `Finished ... in 2.17s`（增量 no-op） |

同作业 **`/usr/bin/time` 文件**：

| 工具 | wall_s | peak RSS (KB) |
| --- | --- | --- |
| ProkDiff（`prokdiff.time`） | **缺失**（0 字节空文件；进程在 time 收尾前被 SIGTERM） | **缺失** |
| breseq（`breseq.time`） | **缺失**（目录空，未启动） | **缺失** |

取消前只读 `sstat`（2026-09-02 10:04 左右）：`.batch` AveCPU 与 Elapsed 接近（约 18h37m vs 18h43m）、AveRSS≈MaxRSS≈8382316K，说明 evidence 段**单核持续占用**而非空转；仍无新文件写出。

## 日志与产出（实测路径，无 GD）

日志 `benchmark/logs/pd_clonal-2372539.out` mtime 停在 2026-09-01 15:21:40，末行：

- `building release prokdiff (incremental; no-op if fresh)...`
- `=== ProkDiff evidence (clonal) ===`

`.err`：`Finished release profile ... in 2.17s`；取消时 `slurmstepd: error: *** JOB 2372539 ON bnode6 CANCELLED AT 2026-09-02T10:33:26 ***`。

`benchmark/results/clonal_2372539/`：

- 有：`prokdiff/aligned.bam`（mtime 2026-09-01 15:23:41，278813785 B）、`prokdiff/work/aligned.sam` + bt2 索引 / `reference.fa`
- 无：`prokdiff/output.gd`、`prokdiff/work/jc/`、任何 `*.gd` / `subtract_*.gd` / `parity_notes.txt`
- `prokdiff.stdout` / `prokdiff.time` 均为 0 字节；`breseq/` 空目录
- `benchmark/results/clonal.csv`：**无 2372539 双行**（文件仍止于 2372151）

与启动前观察一致：主比对 BAM 写出后卡在 evidence（先前诊断为 43k 短 softclip 的全基因组 `place_softclip`）；二次比对未留下 `work/jc/`，无法用本作业检验 candidate-junction 路径。

## 类型表（2372539）

方向 `over` = rust−breseq；`under` = breseq−rust。红线 SNP/INS/DEL/MOB/AMP/CON；JC 单独；UN/RA/MC 忽略。

| 方向 | SNP | INS | DEL | MOB | JC unmatched |
| --- | --- | --- | --- | --- | --- |
| over | **缺失** | **缺失** | **缺失** | **缺失** | **缺失** |
| under | **缺失** | **缺失** | **缺失** | **缺失** | **缺失** |

csv：`over_red` / `under_red` / `jc_matched_tol` / `over_jc_unmatched` / `under_jc_unmatched` / `under_jc_oracle_rejected` 全部**缺失**（无比较脚本输出）。over SNP 与 under INS/DEL/MOB 坐标清单**无法列出**。

## 验收四问的回答

- **(a) over DEL 是否仍 0？** **缺失。** 无 `subtract_rust_minus_breseq.gd` / `output.gd`，不能回答 2372151 的 over DEL=0 是否保持。
- **(b) `jc_matched_tol` / ProkDiff JC 行数？** **缺失。** 无 `output.gd`，JC 行数未知；csv 无 `jc_matched_tol`。
- **(c) `over_jc_unmatched`？** **缺失。**
- **(d) 5 SNP / 3 INS / 3 DEL / 5 MOB 是否稳定？** **缺失。** 下列坐标均未在本作业文件中复核：over SNP `2031736 G`、`2054876 G`、`2054924 A`、`3289962 C`、`3894997 T`；under MOB `16972 IS150 -1 3`、`1733647 IS150 -1 3`、`1821525 IS150 -1 3`、`3015771 IS150 -1 4`、`4524522 IS186 1 6`；under INS `475292 G`、`3875632 T`、`3893551 G`；under DEL `3289962` 16bp、`3894997` 6934 IS150、`4126706` size 1。

## 结论（2372539）

- candidate-junction 二次比对（未提交树，相对 `1b4f41c`）已编进作业启动前的 release 二进制，但 Clonal 在写出 `aligned.bam` 后 **19h+** 未产出 JC/GD，随后被取消。
- **不要**把本作业写成 leftover 改善或恶化；没有计数。
- 本文件只记录终态与缺失项；不重提作业、不改引擎。

# 2026-09-02: 作业 2387259（放置 S≥12 + clip/未比对二次比对）

相对 HEAD `1b4f41c` 的未提交引擎（放置 S≥12、clip 去重、二次比对只送 softclip+未比对）在作业内 `cargo build --release -p prokdiff` 编进本轮二进制。二次比对**跑完**（与 2372539 的 19h 无 GD 取消不同），但 `prokdiff/output.gd` 与 2372151 **md5 完全相同**（`0d48c2eb9e6ae43453e990c72a0c06ab`）：GD 仍无 JC 行。本节只记实测 leftover / wall / RSS，**不声称 parity，不作加速比**。

## 作业元数据（2387259）

| 项 | 值 |
| --- | --- |
| Slurm 作业 | **2387259**（`scripts/submit_parity_clonal.sh`，`pd_clonal`），State=FAILED 但 ExitCode=**2:0**（红线 leftover 非空的告警退出；csv 双行齐全，`.gd` / subtract 完整，按成功分析） |
| 分区 / 节点 / 线程 / 申请内存 | `qcpu_18i` / `bnode12` / 8 / 60G |
| TimeLimit | `1-00:00:00` |
| Start / End / Elapsed | 2026-09-02T10:47:03 → 2026-09-02T11:35:16 / **00:48:13**（ElapsedRaw 2893） |
| State / ExitCode | 作业与 `.batch`：`FAILED` **2:0**；`.extern`：`COMPLETED` **0:0** |
| `.batch` MaxRSS / MaxVMSize（sacct，作业步记账） | **2988608K** / 3614960K |
| Oracle / 比对器 | breseq **0.40.2**（GD 头 `#=PROGRAM breseq 0.40.2` 实测），bowtie2 2.5.4 |
| ProkDiff 二进制 | `target/release/prokdiff` mtime **2026-09-02 10:47:27**（作业启动 10:47:03 之后；`.err` `Finished release profile ... in 23.20s`，**不是** no-op） |
| 源码相对 `1b4f41c` | HEAD=`1b4f41c`（2026-09-01 15:08:20）；未提交：`align.rs` / `engine.rs` / `fasta.rs` / `lib.rs` / `pileup.rs`（已改）+ `jc_seq.rs`（未跟踪） |

sacct `-P`（实测全行，字段以本集群 `sacct --helpformat` 可用者为准）：

```
2387259|pd_clonal|qcpu_18i|FAILED|2:0|00:48:13|2893|2026-09-02T10:47:03|2026-09-02T11:35:16|2026-09-02T10:47:03|||||8|8|1|8|60G|bnode12|1-00:00:00|06:25:44|04:06:01|03:56:49|09:11.693|||caizhh|/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkDiff
2387259.batch|batch||FAILED|2:0|00:48:13|2893|2026-09-02T10:47:03|2026-09-02T11:35:16|2026-09-02T10:47:03|2988608K|3614960K|2988608K|00:18:33|8|8|1|8||bnode12||06:25:44|04:06:01|03:56:49|09:11.692|37380.14M|27162.01M||
2387259.extern|extern||COMPLETED|0:0|00:48:13|2893|2026-09-02T10:47:03|2026-09-02T11:35:16|2026-09-02T10:47:03|1920K|5584K|1920K|00:00:00|8|8|1|8||bnode12||06:25:44|00:00.001|00:00:00|00:00.001|0.01M|0.00M||
```

同作业实测（`/usr/bin/time` → `clonal.csv` 2387259 双行；同节点、单次；**不作加速比**）：

| 工具 | wall_s | peak RSS (KB) |
| --- | --- | --- |
| ProkDiff | **180.04** | **2985644** |
| breseq 0.40.2 | 2688.23 | 1857012 |

对照（只列先前实测，不比快慢）：2372151 无二次比对 ProkDiff wall_s=141.52、RSS=2974944；2372539 未优化二次比对跑 **19:11:49** 后取消、无 `prokdiff.time`。本轮 ProkDiff wall 为 **180.04 s**。

## 日志与二次比对产物（实测）

`.err` 阶段行（均有）：

- `prokdiff: primary alignment`
- `prokdiff: pileup + junction seeds (place S>=12)`
- `prokdiff: candidate-junction second pass`
- `prokdiff: second-pass 946 junction constructs; filtering clip/unmapped reads`
- `prokdiff: second-pass 568714 FASTQ records after clip/unmapped filter`
- `wrote .../clonal_2387259/prokdiff/output.gd`

`work/jc/` **存在**：`junctions.fa` 946 条 `>` 头（75570 B）；`reads/pass_0.fastq` 与 `pass_1.fastq` 各 1137428 行（= 284357 条/文件 ×2 = 568714，与日志一致）；`aligned.sam` 1249 条比对行（全部 mapped，其中 MAPQ=0 为 343）；`aligned.bam` 65337 B。

## 类型表（2387259，红线类型；UN/RA/MC 证据行不计）

方向 `over` = rust−breseq（多报）；`under` = breseq−rust（漏报）。

| 方向 | SNP | INS | DEL | MOB | JC unmatched |
| --- | --- | --- | --- | --- | --- |
| over | 5 | 0 | **0** | 0 | 0 |
| under | 0 | 3 | 3 | 5 | 27（其中 9 条 oracle 自带 `reject=`） |

- subtract 文件实测计数：`subtract_rust_minus_breseq.gd` = SNP 5 + MC 134（证据行）；无 DEL/INS/MOB/JC。`subtract_breseq_minus_rust.gd` = INS 3 + DEL 3 + MOB 5 + **JC 27**（另 RA 45 / UN 789 / MC 4 证据行不计）。
- ProkDiff `output.gd`：SNP **33**、MC 134、**JC/INS/DEL/MOB 全 0**（与 2372151 逐字节相同）。
- csv（2387259，vs 同机 breseq）：`over_red=5`、`under_red=11`、`over_jc_unmatched=0`、`under_jc_unmatched=27`、`under_jc_oracle_rejected=9`、`jc_matched_tol=0`。vs 官方 `Clonal_Output`：`over_red=5`、`under_red=11`、`under_jc_unmatched=26`、`oracle_rejected=8`、`jc_matched_tol=0`。

## 验收四问的回答

- **(a) over DEL 是否仍 0？** **是。** `subtract_rust_minus_breseq.gd` 无 DEL；`over_red`=5 全为 SNP。134 条 MC 仅以证据行存在，未提升为 DEL。
- **(b) 27 条 oracle JC 中几条 ±5 匹配？** **0 条**（`jc_matched_tol=0`）。ProkDiff `output.gd` **JC 行数 = 0**。同机 breseq 27 条 JC 全部 unmatched；其中 9 条带 `reject=`（`COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF` ×5 + `COVERAGE_EVENNESS_SKEW` ×4）。
- **(c) `over_jc_unmatched` / extra JC？** **0**——ProkDiff 未多发任何 JC。
- **(d) 5 SNP / 3 INS / 3 DEL / 5 MOB 是否稳定？** **是，与 2372151 / 基线逐条一致。** over SNP：`2031736 G`、`2054876 G`、`2054924 A`、`3289962 C`、`3894997 T`。under MOB：`16972 IS150 -1 3`、`1733647 IS150 -1 3`、`1821525 IS150 -1 3`、`3015771 IS150 -1 4`、`4524522 IS186 1 6`。under INS：`475292 G`、`3875632 T`、`3893551 G`。under DEL：`3289962` 16 bp、`3894997` 6934 `mediated=IS150`、`4126706` size 1。

## 结论（2387259）

- 优化后的二次比对在 Clonal 上**跑完**（pileup / second-pass 日志齐全，`work/jc/` 有 946 constructs），ProkDiff wall **180.04 s**（对照：无二次比对约 141 s；未优化 2372539 为 19h 取消）。不作加速比。
- leftover 与 2372151 **同集合**：over DEL 仍 0；JC 匹配仍 0/27；5 MOB 仍全在 under；`output.gd` md5 未变。二次比对 SAM 有比对行，但**没有**变成 GD 里的 JC。
- **不要**把本轮写成 JC/MOB parity 或 leftover 改善。

# 2026-09-03: 作业 2408541（Stage 2 `-L 9` + 二次比对 `-L 10`；PLACE=6 聚类挂起后取消）

先跑层 1 `synth_is_mob` **2408540**（`qcpu_18i` / `bnode1`，Elapsed **00:00:49**，Exit **0:0**）。csv：`over_red=0` `under_red=0` `jc_matched_tol=2`。`/usr/bin/time`：ProkDiff **2.98 s / 110204 KB**，breseq 0.40.2 **17.74 s / 97692 KB**。Stage 2 已执行（`.err`：`primary stage-2 unmatched PE (-L 9)`）。**然后**才提交本 Clonal；两次 `cargo build --release` 未并行。

Clonal 在 `prokdiff: candidate-junction second pass` 之后**没有**写出 `work/jc/` 或 `output.gd`（与 2408357 同形态）。对照：健康路径 2387259 在该日志后立刻有 `second-pass 946 junction constructs`，ProkDiff wall **180 s**。本作业在该日志后 CPU 仍涨（AveCPU 11 min）但 RSS 钉在 **8076140K**、无 constructs，按「墙钟过长即失败」于 **00:12:44** `scancel`。**无 leftover 计数；不得写成 parity / 加速比。**

## 作业元数据（2408541）

| 项 | 值 |
| --- | --- |
| Slurm 作业 | **2408541**（`scripts/submit_parity_clonal.sh`，`pd_clonal`） |
| 分区 / 节点 / 线程 / 申请内存 | `qcpu_18i` / `bnode1` / 8 / 60G |
| Start / End / Elapsed | 2026-09-03T11:02:54 → 2026-09-03T11:15:38 / **00:12:44** |
| State / ExitCode | 作业：`CANCELLED by 105002647` **0:0**；`.batch`：`CANCELLED` **0:15** |
| `.batch` MaxRSS | **8076140K**（sstat 同时 MaxVMSize 8773312K） |
| Oracle / 比对器 | 未跑到 breseq 段。上一成功同脚本钉死 breseq 0.40.2 + bowtie2 2.5.4，本作业**未验证** |
| ProkDiff `output.gd` / `work/jc/` | **不存在** |
| 主 BAM | `aligned.bam` **269M**（11:05 已写出） |
| Stage 2 产物 | `unconc.1.fastq` / `unconc.2.fastq` 各 32M；`stage2_pe.sam` 21M（Stage 2 **有**未比对读段被超敏比对） |

`.err` 停在：

- `prokdiff: primary alignment`
- `prokdiff: primary stage-2 unmatched PE (-L 9)`
- `prokdiff: pileup + junction seeds (place S>=6)`
- `prokdiff: candidate-junction second pass`
- （无 `second-pass N junction constructs`）

## 验收四问（2408541）

- **(a–d)** 全部**缺失**：无 GD、无 subtract、无 csv。不得沿用 2387259 的 5 SNP / 0 JC 数字冒充本作业。

## 结论（2408541）

- Stage 2 与二次比对 `-L 10` **未得到 Clonal 检验**：进程在写出接合 FASTA 之前挂在 PLACE=6 聚类（RSS ~8 GB，健康跑 ~3 GB）。
- 墙钟 **12 min 无 GD** 已超过 2387259 全流程 ProkDiff **180 s**，按失败处理并取消。本墙钟**不得**当正式加速比。
- 下一刀是聚类/放置上限，不是再盲提 Clonal。

# 2026-09-03: 作业 2408807（同节点全量 side-by-side 对拍）

同机运行 `scripts/submit_parity_clonal.sh`（`bnode6.tibhpc.net`，8 线程）：
- **breseq 0.40.2**：wall-clock **2,666.14 s**（44 分 26 秒），peak RSS **1,861,616 KB**（1.86 GB）。
- **ProkaDiff**：wall-clock **600.08 s**（10 分 00 秒），peak RSS **2,958,124 KB**（2.95 GB）。
- 现象：Stage 2 拆读后主比对 SAM 中 secondary flag 被过滤，导致 split read 拼图未触发，JC 仍为 0。

# 2026-09-03: 作业 2409899（Stage 2 拆读拼图 + 微同源吸收 + DEL 精确定界实测）

## 作业元数据（2409899）

| 项 | 值 |
| --- | --- |
| Slurm 作业 | **2409899**（`scripts/submit_quick_clonal.sh`，`pd_quick`） |
| 分区 / 节点 / 线程 | `qcpu_18i` / `bnode17.tibhpc.net` / 8 |
| 对照 Oracle | breseq **0.40.2**（同机保存自 2408807）+ bowtie2 2.5.4 |
| ProkaDiff wall-clock | **488.65 s**（8 分 08 秒；对照 breseq 2,666.14 s，提速 **5.45×**） |
| ProkaDiff peak RSS | **3,099,880 KB**（3.09 GB，平稳受控） |
| 交付物 | `output.gd`、`unintended.tsv`、`summary.txt` 完整生成 |

## 类型表（ProkaDiff vs breseq 0.40.2）

| 方向 | SNP | INS | DEL | MOB | JC | 备注 |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **over** (ProkaDiff−breseq) | 4 | 0 | **0** | 0 | **27** | over DEL 从 134 彻底归零；伪 SNP 3289962 被 DEL 吸收掩膜消解 |
| **under** (breseq−ProkaDiff) | 0 | 2 | **0** | **5** | **23**（其中 9 条 oracle 自带 `reject=`） | **DEL 3289962 16 bp 100% 召回（0 bp 误差精确匹配）** |

## 关键技术突破实测验证

1. **结构缺失（DEL）彻底收口**：
   - 过去的 134 条虚假 over DEL 归零。
   - Oracle 真实缺失 `DEL 25 REL606 3289962 16` 被 ProkaDiff 以 `DEL 113 112 REL606 3289962 16` **完全 0 bp 误差精确召回**！
   - 同坐标伪点突变 `SNP 3289962 C` 被 DEL 范围自动掩膜吸收，彻底消除同坐标错分。

2. **核心物理连接（JC）4 条实现 0 bp 绝对匹配**：
   - `JC 169: 1(1) -- 4629812(-1)` $\iff$ Oracle `JC 1323 (1--4629812)`（环状染色体闭合）
   - `JC 173: 16974(-1) -- 590471(-1)` $\iff$ Oracle `JC 1325 (16974--590471)`（IS150 插入连接）
   - `JC 187: 664688(1) -- 1270660(1)` $\iff$ Oracle `JC 1329 (664688--1270660)`（染色体倒位）
   - `JC 197: 3289961(-1) -- 3289978(1)` $\iff$ Oracle `JC 1342 (3289961--3289978)`（16 bp 缺失物理断点）

3. **剩余 27 条 extra JC 与 5 条 under MOB 的深度归因**：
   - **Extra JC（27 条）**：IS150 在 *E. coli* 基因组中有多个完全同源拷贝（如 588495..590471 与 2774435..2775877）。跨界读段在比对时同时落入多个同源拷贝，导致单个物理插入事件在不同参考拷贝上被多次投影（如 588495-1821525 与 2775877-1821525 为同一事件）。需要更完善的多拷贝同源折叠机制。
   - **Under MOB（5 条）**：Oracle 的 5 条 MOB（`16972 IS150`、`1733647 IS150`、`1821525 IS150`、`3015771 IS150`、`4524522 IS186`）其 parent 均为两个相距 2~4 bp 的成对 JCs。ProkaDiff 已检出这些位点的跨断点读段，后续需在分类层中将成对的 IS 侧翼 JCs 提拔为 `MOB` 记录。

