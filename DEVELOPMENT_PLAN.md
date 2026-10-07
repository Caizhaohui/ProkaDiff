# ProkaDiff 开发计划

> 当前执行顺序，更新于 2026-10-07。一次只执行一个里程碑；每个里程碑同时交付实现、测试、文档，并按 [AGENTS.md](AGENTS.md) 的计算节点规则验证。本文是根目录的执行入口；阶段字母只在各自历史任务书内有效。

## 目标与边界

将 **Evidence → Differential Variant → Intended Verification → Mechanistic Association → AuditResult** 连成一条可追溯的数据链。所有产品输出从同一 `AuditResult` 生成，报告中的每一项都能回指差分事件及其实际证据。观察、空间关联、机制假说和实验因果遵循 [科学解释模型](docs/scientific_model.md)，不能相互替代。

继续遵守现有产品约束：必须同时提供出发株与编辑株 WGS；NCBI 参考只是坐标骨架；运行时调用系统 Bowtie2，RA/MC/JC 由 Rust 实现，breseq 只作 oracle；BAM 使用 `noodles`；分类顺序为 structural → near_homolog → scattered_snv；第一期编辑器仅 cas9、cas12a、dsb。I/O 以 [schema](docs/schema.md) 为准，对拍以 [parity](docs/parity.md) 为准。不得在登录节点运行 FASTQ、Bowtie2、breseq 或真实数据测试。

## 文档来源与阶段对照

本计划综合本地 `markdown/` 中的旧版 `DEVELOPMENT_PLAN.md`、`NEXT_STEPS.md`、Breseq-Parity + Reference-based Off-target 任务书、Post-edit Genome Audit 任务书及执行计划、Validation & Scientific-Correctness Sprint、Real-world Validation Sprint。它们记录了不同开发周期，不能将相同的 Phase A/B/C 当作同一个阶段。当前行为以代码、[schema](docs/schema.md)、[parity](docs/parity.md) 和可复核的仓库内结果为准；历史任务书中的完成标记不自动构成验证证据。

| 历史阶段 | 当前状态与用途 |
| --- | --- |
| 基础计划阶段 0–4 | Workspace、GD、引擎、双株差分和旧版输出已有实现；继续维护，不重复执行旧清单。 |
| 审计计划 Phase A–G | 逐编辑状态、多维注释、`AuditResult`、TSV、Markdown、结构注释与批次规格已有不同程度实现；本计划处理跨阶段的数据契约和可追溯性。 |
| 科学修正 FIX-001–020 | 作为回归约束与 oracle 门禁；不能把夹具通过等同于外部验证。 |
| 真实数据验证 Phase A–H | 作为分数据集验证和发布门禁；每个数据集只凭仓库内可复核的输入、比较结果和溯源提升状态。 |

## 当前进度

- 已实现：
  - M0 验证声明与可复核输入收口已完成（2026-09-18）。
  - M1 Evidence → Differential Variant 已闭环完成（2026-09-22）：定义了规范化差分事件（`CanonicalEvent`）、确定性内容截断哈希身份（`DV1 EventId`）及其碰撞防护；实现了基于最小费用最大流的确定性容差亲本扣除算法（先最大基数、再最小坐标距离代价、后规范字节序决胜；严格限定 JC 符号重叠相等、MOB duplication_size 相等、DEL 双方均 >2 bp 预条件、±5 bp 边界）；实现了对称且独立的 MOB 组分 JC 吸收；建立了权威 `EvidenceReferences` 和三态 `EvidenceSummary`（仅从真实 GD 证据记录或亲本引用解析，缺失写 `NA`，禁止由变异类型逆向推断）；移除了所有全 N 伪参考，所有坐标与碱基严格校验。全工作区 257 项单元测试及 `qcpu_18i` 上全部 M1 端到端双株 FASTQ 差分与 breseq 对拍作业均已通过。
  - M2 Differential Variant → Intended Verification 已闭环完成（2026-09-23）：逐编辑评估权威引用 M1 的 `EventId`；建立了 `IntendedEventRelationship` 并严格区分 `ExpectedConstituent`、`PartialObservation` 与 `UnexpectedAtLocus`；实现评估先于掩码，只有唯一归属且无异常判定的预期组成事件才移出非预期视图，靶位异常等位基因、额外变异及异常接合点完整保留；实现了基于生物拟合度与确定性 `EventId` 决胜的删除及 cassette 候选选择，消除输入行序依赖；确立了 MC 纯证据诊断契约（无 `EventId`、缺失写 `Unknown/NA`、亲本 MC 自动过滤）；保持了全部既有公开 TSV 表头、列序及 summary 键不变，`edit_outcomes.tsv` 仅在输出边界将 `EventId` 映射为代表数字 GD ID；全工作区单元测试与 M2 专项回归测试全过；`qcpu_18i` 上 Layer 1 `synth_parent_child` FASTQ 差分流水线与 `edit_outcomes.tsv` 闭环验证通过（Job 2784683）。
  - M3 Intended Verification → Mechanistic Association = **COMPLETE (2026-09-24)**：`EventId` 是内部权威突变身份；候选位点只扫描一次；`CandidateSearchStatus` 显式区分未扫描、零候选、有候选但事件无关联；不存在 `SITE_UNKNOWN`、伪造 mismatch/distance/PAM 或位置回退关联；OnTarget 缺失测量输出 `NA`；未解析的 `EventId` 触发 typed error；全部合格关联均保留；空间与序列机制关联不证明 Cas 切割因果。全工作区 `cargo test --workspace` 在 `qcpu_18i` 通过（Job 2791150）；M3 FASTQ/Bowtie2 E2E 通过（Job 2791151）。
- M4 Mechanistic Association → AuditResult 与全部输出 = **COMPLETE (2026-09-24)**：以 `DV1 EventId` 作为公开差分事件身份，由完整 `AuditResult` 与共享 `AuditSummary` 生成全部产品；集中管理 V1/V2 模式，V2 默认且相对 V1 仅追加列；M3 闭环修订 `v0.3.0-m3` / `0f79495` 的七个 V1 产品文件与当前 V1 输出逐文件字节比较一致；V2 身份/schema、阶段化发布及回滚恢复、跨输出一致性均通过回归验证。最新 `qcpu_18i` FASTQ/Bowtie2 产品验证 Job 2792159 完成，覆盖默认 V2 与显式 V1 产品生成。
- M5-SR1 Repeat-Ambiguous Junction Recovery = **COMPLETE (2026-09-28)**，冻结提交 `f30f93be5439d2f8eb32f79719fe8482713fcf5a`。科学结论：ProkaDiff 保留高拷贝重复相关接头证据；在拷贝身份未独立解析时，不发出拷贝特异结构调用。诊断旁路不直接生成 JC、DEL、MOB、EventId、DifferentialEvent 或关联。验收：合成 Job 2815876 与焦点 B21_3_1 Job 2815888 通过。**溯源限定**：这两份作业不是 `f30f93b` 的干净检出；日志 HEAD 为 `08d41edd7f6ce7f6d287fb19cebdedd4e3bcb09c`，工作树含当时未提交的 SR1 实现；2815876 保存的 `git diff HEAD` SHA256 为 `e6111f92e94cb6c5986d37a2dff07256db4a4f49f08bae223472b1dee10701e7`。`f30f93b` 随后冻结已验证实现，不能写成这两份作业直接跑了干净的 `f30f93b`。
- M5-BL21-P1 单样本参考相对复核 = **COMPLETE (2026-09-28)**。三个异常克隆各自相对参考 `GCF_013167015.1` 运行，不是 MATCHED_PARENT，也不是 PEER_COMPARATOR。作业：B21_3_1 Job 2815888，B21_4_1 Job 2816157，B21_2_1 Job 2816170。科学结论：三克隆均有覆盖 lacZ 的参考相对结构异常；ProkaDiff 保留重复相关接头证据，但未将重复家族解析到单一拷贝；三份均为 `exact_breakpoint_supported=false`。2816157 与 2816170 明确记录 `SR1_COMMIT=f30f93b`。对照与来源分离见 [M5_BL21_P1_RECONCILIATION.md](benchmark/real_world/M5_BL21_P1_RECONCILIATION.md)。`MAX_STAGE2_FAMILY_GEOMETRIES=5000` 在这三份 BL21 数据中未触发；支持表述仅为“当前上限在所测 BL21 数据集中未触发”，不能外推为对重复丰富细菌基因组普遍足够。
- M5 PEER_COMPARATOR 矩阵 = **COMPLETE (2026-10-06 执行 / 2026-10-07 收口)**。验证提交 `b558371a85367bcbb1b97a63647aff3d663f2bd7`。作业：2888235（B21_3_1 vs B21_4_1，G1+G2）、2888236（B21_4_1 vs B21_3_1，G1+G2）、2888237（B21_2_1 vs B21_3_1，G1+G2）；三份均 `COMPLETED` / `0:0`；六个 guide 格子产品验证 PASS。六格均为 M2=`MISSING`、预期相关差分事件 0、DV1 matched/unexpected 列表空、`false_complete=0`、`gate_passed=true`、搜索已执行（G1：1 候选含 0-mismatch；G2：2 候选含 0 与 4 mismatch）、关联数 0、七个 M4 产品一致。方向差分事件数：3→4 为 32，4→3 为 44，2→3 为 54；方向不对称可接受。权威收口记录：[M5_CLOSEOUT.md](benchmark/real_world/M5_CLOSEOUT.md)。
- M5 整体 = **COMPLETE (2026-10-07)**。M6 未开始；M5 不再阻塞启动 M6，但本次收口不启动 M6。
- 待证实（M6+）：`VALIDATION_REGISTRY.tsv` 对 PRJNA1088182 和 PRJNA884016 的完成声明，须与随仓库交付的非 FASTQ 输入、truth、实际 calls、比较结果一致。局部忽略的结果目录不能作为可复核发布证据。

## 顺序与验收

按下表自上而下执行。任何里程碑未达到验收条件时，保持在该里程碑；后续里程碑不提前实现。每个里程碑均包含实现、针对性测试与文档更新。原始 FASTQ、BAM、私有测序数据不入库。

| 顺序 | 里程碑 | 交付与通过条件 |
| --- | --- | --- |
| **M0：验证声明与可复核输入收口** | 已完成（2026-09-18） | PRJNA1088182/PRJNA884016 的 manifest、truth 与 SRA 元数据已进入变更集；缺少 ProkaDiff calls 和比较结果的字段已退回 `PENDING`。`check_validation_assets.sh` 检查 clean-checkout 所需的版本控制资产、truth 行数及登记状态；验证说明已写入 real-world README。 |
| **M1：Evidence → Differential Variant** | 已完成（2026-09-22） | 交付稳定的差分事件身份（`DV1 EventId`）与规范化前像生成；实现权威 `EvidenceReferences` 溯源（无 `GdKind` 逆向重构，缺失字段置 `NA`）；实现基于最小费用最大流的确定性容差亲本扣除算法（JC 符号重叠相等、MOB duplication_size 相等、DEL 双方均 >2 bp 预条件、±5 bp 容差边界、最大基数优先与字节序决胜）；实现对称独立的 MOB 组分 JC 吸收；严格坐标、跨度与 DNA/IUPAC 碱基校验；纯内存层 0 单元测试 257 项全过；`qcpu_18i` 集群验证通过（Job 2769120 端到端双株差分 FASTQ 流程，Job 2769121 IS/MOB 对拍，Job 2769122 SNP/INS/DEL 差分对拍，Job 2769123 结构 DEL 对拍）。 |
| **M2：Differential Variant → Intended Verification** | 已完成（2026-09-23） | 逐编辑评估权威引用 M1 的 EventId，保留 Complete/Partial/Missing/UnexpectedStructure 与断点判定；异常靶位事件不因掩码而丢失；支持一事件匹配多声明、cassette 双接头及异常结构判定；edit_outcomes.tsv 正确回指并在输出边界兼容映射；纯内存单元测试与 qcpu_18i 集群端到端验证通过（Job 2784683）。 |
| **M3：Intended Verification → Mechanistic Association** | **COMPLETE (2026-09-24)** | 内部突变身份权威为 `EventId`；单次候选扫描；搜索状态显式区分；无 `SITE_UNKNOWN` 或伪造 mismatch/distance/PAM；无位置关联回退；OnTarget 缺失测量写 `NA`；未知 EventId 投影为 typed error；保留所有合格关联；关联不等于切割因果。`qcpu_18i` workspace 测试 Job 2791150 与 M3 FASTQ/Bowtie2 E2E Job 2791151 通过。 |
| **M4：Mechanistic Association → AuditResult 与全部输出** | **COMPLETE (2026-09-24)** | 公开差分事件身份为 `DV1 EventId`；单一 `AuditResult` 与 `AuditSummary` 提供完整产品状态和共享计数；V2 为默认并仅相对 V1 追加列，V1 与 M3 闭环版本 `v0.3.0-m3` / `0f79495` 的七个产品文件逐字节相同。V2 身份/schema、阶段化发布及回滚恢复、跨输出一致性均通过回归验证。`qcpu_18i` FASTQ/Bowtie2 Job 2792159 的默认 V2 与显式 V1 产品生成通过。冻结的不变量见下节。 |
| **M5：BL21 真实数据复核** | **COMPLETE (2026-10-07)** | M5-SR1、M5-BL21-P1 与 PEER_COMPARATOR 六格矩阵均已完成。权威收口：[M5_CLOSEOUT.md](benchmark/real_world/M5_CLOSEOUT.md)。矩阵作业 2888235–2888237（提交 `b558371`）；单样本 2815888 / 2816157 / 2816170；SR1 合成 2815876。保留无核实亲本、无 Cas9 因果、无精确重复拷贝解析等限制。 |
| **M6：PRJNA1088182 队列复核** | 未开始（M5 不再阻塞） | 以随仓库交付的 manifest、publication truth 和评估器逐样本验证编辑状态与关联陈述；仅在测量结果完整时更新登记表。 |
| **M7：PRJNA884016 组装 truth 复核** | M6 后 | 以可解析的组装 truth、ProkaDiff calls 和结构比较结果核对 MOB/JC 及断点误差；不得以组装差异直接替代测序检出结果。 |
| **M8：REL606 引擎扩展对拍** | M7 后 | 记录与 breseq 的双向 GD 比较、差异解释及同作业 wall/RSS；性能或 parity 声明须满足 [parity](docs/parity.md) 的门禁。 |
| **M9：报告审阅与发布门禁** | M8 后 | 人工审阅真实报告，核对每个 finding 的事件与证据回指、科学措辞和溯源；汇总 M5–M8 的可复核记录及尚未验证的能力，只发布有证据支持的结论。 |

M0–M5 已闭环通过。M6–M9 尚未开始；下一项是 M6（PRJNA1088182 队列复核）。批次模式、CAST/IS110、长读长、混群以及未经外部 oracle 验证的评分保持在后续路线图。

## 冻结的 M4 不变量

- `DV1 EventId` 是公开差分事件的权威身份；`VAR_####` 仅作显示标签。
- 原始数字 GD ID 仅作来源记录/兼容标识；`mut_<GD ID>` 仅作兼容输出标识。二者均不能替代 EventId 作为关联键。
- `AuditResult` 是唯一科学与产品聚合；`AuditSummary` 是重复公开计数的共享投影。
- `AuditResult.associations` 是权威关联集合；`variant_site_links` 仅用于兼容投影。
- `event_index` 解析每个真实 `DifferentialEvent`，包括正常变异视图中隐藏的 MOB companion JC；仅证据 MC/RA 记录不获得 EventId。
- Schema V2 是默认；Schema V1 与 M3 遗留产品字节兼容；V2 相对 V1 仅追加已批准字段。
- 产品发布先写入 staging；普通发布错误通过已验证的回滚与恢复处理。多文件最终发布期间发生硬进程崩溃（包括 SIGKILL）仍属于已接受的文件系统级限制。

## 每个里程碑启动与结束规则

启动前读本文、相关代码和测试；触及 I/O 或 oracle 时读 [schema](docs/schema.md) 与 [parity](docs/parity.md)。先跑基线 `cargo test --workspace`，但任何会调用 FASTQ/Bowtie2/breseq 的测试必须改在 `qcpu_18i` 执行；登录节点只允许纯内存层 0 或静态检查。结束时记录实际测试、作业号、差异和未解决项；验证状态由可核查记录决定，不由任务书勾选决定。
