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
