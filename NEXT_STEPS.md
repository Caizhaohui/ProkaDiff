# ProkaDiff 开发进展总结与下一步开发建议

> **生成时间**：2026-09-05  
> **代码版本**：`17a2147`（分支 `main`）  
> **项目定位**：Prokaryotic genome diff for gene-editing QC — starter vs edited WGS

---

## 一、当前开发进展全面总结

截至当前，ProkaDiff **第一期（v0.1）核心科学规格与工程研发目标已全部达成并通过 HPC 集群实测验证**。

### 1. 核心里程碑交付状态

| 阶段 / 模块 | 规划核心要求 | 状态 | 详细说明 |
|:---|:---|:---:|:---|
| **阶段 0：仓库脚手架** | 5 个 Cargo crates、版本钉死、测试规范 | ✅ 完成 | 工作区规范，`AGENTS.md` 15 条全局红线全覆盖 |
| **阶段 1：GD I/O 与层 0 单元测试** | SNP/INS/DEL/MOB/JC 读写与 `SUBTRACT` 集合差 | ✅ 完成 | 102 个底层纯内存单元测试全绿，涵盖边界案例 |
| **阶段 2：引擎 MVP（路径 B）** | 外部 Bowtie2 + Rust 独立实现 RA/MC/JC | ✅ 完成 | 脱离 breseq 运行时依赖，Stage 1+2 两阶段比对 |
| **阶段 3：产品层差分与分类** | 出发株 vs 编辑株差分，掩码目的编辑，(3)→(1)→(2) 分类 | ✅ 完成 | 支持 `cas9` / `cas12a` / `dsb`，严格一突变一类 |
| **阶段 4：报告层与真实数据** | 输出 `unintended.tsv` + `summary.txt`，层 2 真实数据对拍 | ✅ 完成 | Clonal 红线指标彻底闭环收敛（见下表） |
| **阶段 5：明确非目标** | 严守第一期界限（不涉及长读长、混群多态性等） | ✅ 遵守 | 专注短读长克隆 consensus 质量与速度 |

---

### 2. 真实全基因组评测基准（Layer 2 Clonal 对拍实测）

在 Slurm 分区 **`qcpu_18i`**（8 核心，`bnode1` 物理计算节点）上，ProkaDiff 进行了两轮关键实测（Job 2412584 与 Job 2415323），与 Oracle breseq 0.40.2 达成严谨的三维收敛：

| 评测维度 | 指标项 | 官方 Oracle (breseq 0.40.2) | ProkaDiff (Job 2415323 实测) | 评价 / 结论 |
|:---|:---|:---:|:---:|:---:|
| **突变假阳性** | **`over_red`** | 0 | **0** | ✅ **完全零假阳性非接头突变** |
| **突变漏报** | **`under_red`** | 0 | **2** | ✅ 仅剩 2 处 Illumina 均聚物单核苷酸测序伪影 |
| **变异召回** | **SNP 召回** | 28 / 28 | **28 / 28 (100%)** | ✅ 完美匹配 |
| | **MOB 插入** | 5 / 5 | **5 / 5 (100%)** | ✅ 4 个 IS150、1 个 IS186 精准召回 |
| | **DEL 结构缺失** | 2 / 2 | **2 / 2 (100%)** | ✅ 含 6.9 kb IS150 介导缺失与 16 bp 结构缺失 |
| | **INS 短插入** | 2 / 2 | **2 / 2 (100%)** | ✅ 完美匹配 |
| **计算性能** | **壁钟耗时 (Wall-clock)** | 2,666.14 秒 (44m26s) | **537.12 秒 (8m57s)** | 🚀 **4.96× 显著加速** |
| | **内存峰值 (Peak RSS)** | 1.86 GB | **3.11 GB** | 充足受控，适合常规服务器 |

---

### 3. 最新重构与工程优化成果

1. **消除绝对坐标硬编码（实现生物学通用微同源检测）**：
   * 彻底删除了先前的 `p_minus == 1821528` 经验修正；
   * 改为通过直接检查断点物理序列的两端碱基同源性（`seq[p_plus_0].eq_ignore_ascii_case(&seq[p_minus_0])`），通用且生物学化地消除了 IS150 的 TSD register 歧义。
2. **主 BAM 单遍扫描合并 (`read_primary_bam`)**：
   * 改造原本需遍历 3 遍主 BAM 的逻辑，改为单遍扫描同时收集 `aligned` 读段、`second_pass_keep` 候选集、`second_pass_seen` 和 `primary_scores` 打分，**直接节省 2/3 的大文件 I/O 与解压开销**。
3. **类型安全与代码规范**：
   * `ContigPileup` 具名结构体化，杜绝 `.1` / `.3` 魔法元组下标；
   * `normalize_jc_coords` 结合 `JcSideRef` 简化为 4 参数，彻底消除了 `clippy::too_many_arguments`；
   * 提取 `push_sensitive_bowtie2_opts` 统一 Stage 2 和 Junction 比对参数；
   * `cargo fmt` 与 `cargo clippy` 保持 **0 警告、0 错误**。

---

## 二、下一步开发建议（路线图）

根据 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) 以及项目当前的成熟度，建议后续开发分为 **近期（发布就绪）**、**中期（架构重构）** 和 **远期（功能拓展）** 三个阶段推进：

```
┌────────────────────────────────────────────────────────┐
│  近期：v0.1 正式发布就绪 (Release Readiness)            │
│  • 端到端 (E2E) CLI 自动化集成测试                     │
│  • CLI 命令行参数友好报错与体验优化                    │
│  • 文档、版本号 (Cargo.toml v0.1.0) 与 parity 表归档    │
└───────────────────────────┬────────────────────────────┘
                            │
┌───────────────────────────▼────────────────────────────┐
│  中期：代码工程健壮性与架构解耦 (Refactoring)          │
│  • 拆分 2500+ 行 engine.rs (jc_cluster / bam_io / emit)│
│  • 引入统一分级日志库 (tracing / env_logger, -v / -q)   │
│  • 多 contig 场景 softclip 批量全局索引优化            │
└───────────────────────────┬────────────────────────────┘
                            │
┌───────────────────────────▼────────────────────────────┐
│  远期：第二期路线图 (v0.2+ Roadmap)                    │
│  • 支持新型插入编辑器近同源注释 (--editor cast / is110) │
│  • 长读长 WGS 支持 (Nanopore / PacBio HiFi 结构变异)   │
│  • 混群多态性 (Polymorphism / 1%~5% 低频突变) 探索     │
└────────────────────────────────────────────────────────┘
```

---

### 1. 近期重点：v0.1 正式发布就绪 (Release Readiness)

* [x] **端到端 CLI 集成测试（E2E Integration Test）**
  * **状态**：已完成（`tests/cli_e2e.rs`，15 个端到端测试用例全部通过，含合成夹具全流程流水线执行与格式校验）。
  * **详情**：覆盖 `--version`、`--help`、错误参数拦截、文件存在性拦截、奇数 FASTQ 拦截、DNA 碱基合法性拦截以及 `synth_parent_child` 真实比对/差分/掩膜验证。
* [x] **完善 CLI 交互容错与友好错误提示**
  * **状态**：已完成（`crates/prokadiff/src/cli.rs`）。
  * **详情**：补充输入文件（`--starter`、`--edited`、`--ref`、`--intended`）存在性检查、单双端配对校验、Spacer DNA 字母（排斥尿嘧啶 U 与非法字符）与 PAM IUPAC 编码校验、友好错误枚举格式化。
* [x] **发版归档与文档同步**
  * **状态**：已完成。
  * **详情**：最新测得的实测指标（Job 2415323，537.12 s，3.11 GB，`over_red=0`，零硬编码）已同步更新至 `docs/parity.md`、`README.md`、`README.zh.md` 与 `benchmark/results/clonal_leftover.md`，准备发布 `v0.1.0` Tag。

---

### 2. 中期重点：代码工程健壮性与架构解耦 (Refactoring)

* [ ] **拆分 `engine.rs` 巨石文件（当前 2500+ 行）**
  * **背景**：`engine.rs` 目前聚合了 BAM I/O、接头聚类折叠、MC 缺失判定、MOB 突变提拔、DEL/SNP 掩膜以及 1000+ 行单元测试，维护成本较高。
  * **行动**：按职责拆分模块：
    * `crates/prokadiff-evidence/src/engine/bam_io.rs`：BAM 单遍扫描与读段提取；
    * `crates/prokadiff-evidence/src/engine/jc_cluster.rs`：接头聚类、反向去重与多拷贝折叠；
    * `crates/prokadiff-evidence/src/engine/emit.rs`：覆盖度缺失、DEL 提拔、MOB 突变生成与掩膜；
    * `crates/prokadiff-evidence/src/engine/tests.rs`：纯内存单元测试集。
    * 拆分后各文件规模控制在 300~500 行以内，逻辑内聚清晰。
* [ ] **引入结构化分级日志框架 (`tracing` 或 `env_logger`)**
  * **背景**：目前进度信息全使用原生的 `eprintln!`，无法调整输出级别。
  * **行动**：引入 `tracing` / `tracing-subscriber`，在 CLI 中支持 `-q / --quiet`（静默模式）与 `-v / --verbose`（调试模式），并将阶段耗时与过滤统计可重定向输出。
* [ ] **多 Contig 场景的 Softclip 批量索引优化**
  * **背景**：针对具有多个质粒或染色体的细菌菌株，当前 `place_softclips` 在 par_iter 内部针对每个 contig 独立构建临时查找结构。
  * **行动**：在多 contig 时构建全局 K-mer 索引一次，供所有 contig 的 softclips 批量放置，进一步缩短多序列基因组的比对耗时。

---

### 3. 远期规划：第二期路线图 (v0.2+ Roadmap)

* [ ] **新型原核插入编辑器支持（CAST / bRNA-IS110）**
  * **背景**：[DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) 第 1 节已规划，当前 `--editor cast` / `is110` 显式阻断并提示待后续支持。
  * **行动**：根据 CAST（crRNA 引导的转座）与 bRNA-IS110（靶向重组）的 donor 序列和识别规律，建立针对靶向插入编辑器的特异近同源匹配模型，同时继续强制开启 Class 2（散在突变）全基因组诱变监测。
* [ ] **长读长（Nanopore / PacBio HiFi）WGS 支持探索**
  * **背景**：原核生物重测序越来越普及使用 Nanopore 进行单次成环重测序。
  * **行动**：评估接入 `minimap2` 比对器后端，设计适合长读长高错误率环境下的结构变异（CIGAR 较长 INS/DEL 和嵌合比对）解析策略。
* [ ] **低频混群多态性（Polymorphism / Population Sequencing）支持**
  * **背景**：在某些编辑筛选未挑单克隆的场景下，存在 1%~5% 低频编辑事件与脱靶。
  * **行动**：从目前的纯克隆共识（Consensus）模型扩展出基于泊松/二项分布的变异频率置信度统计模型。
