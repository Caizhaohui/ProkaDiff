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
