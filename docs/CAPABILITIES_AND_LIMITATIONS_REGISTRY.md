# ProkaDiff Capabilities and Limitations Registry
## 跨里程碑能力与科学局限性权威注册表 (M5–M8)

> **版本：** v1.0  
> **生效里程碑：** Milestone M9 (报告审阅与发布门禁)  
> **治理规范：** Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)  
> **核心准则：** 只有经过可复核证据支持的结论才允许声明；所有已知科学局限如实固定，禁止过度外推与虚假宣传。  
> **因果声明边界：** Spatial association is consistent with a possible guide-dependent event but does not establish cleavage causality.（空间邻近性仅支持候选向导依赖事件假说，不建立 Cas 切割因果。）  

---

## 1. 核心定位与能力层级 (Overview)

ProkaDiff 是面向原核工程菌株的配对全基因组测序（Starter WGS + Edited WGS）高精度变异发现与编辑完整性审计系统。
依据 [`docs/scientific_model.md`](docs/scientific_model.md)，系统严格遵守“**观察优先、机制关联次之、因果证实最后**”的证据层级：

| 能力维度 | 当前已验证能力 (Verified Capabilities) | 明确已知局限 (Documented Limitations) | 路线图未验证/禁止声明 (Out of Scope) |
|---|---|---|---|
| **点突变与小片段 (SNV / Small InDels)** | 单样本 REL606 对拍 100% 召回与精度（28 SNP, 3 INS, 3 DEL 完全一致，delta=0）。 | 低频亚克隆变异（ALT fraction < 0.75）在共识阈值下默认保守过滤。 | 严禁在缺乏正交验证下宣称散发 SNV 属于 Cas 脱靶突变。 |
| **移动元件插入 (IS / MOB)** | 5/5 经典原核转座子插入 100% 精确召回（元件类型、插入坐标、正反链、靶位重复 TSD 长度完全一致）。 | 旁系同源位点易引起单侧接头多重映射，需结合双向接头守恒判定。 | 严禁在缺乏配对亲本对照时将所有 IS 位点作为编辑衍生变异报告。 |
| **染色体大尺度结构变异 (Large Deletions / JCs)** | 在非高拷贝重复区域能够稳定捕获断裂接头（Novel Junctions）与删除断点；支持双向接头守恒平衡对账。 | 处于高拷贝重复元件（如 IS150、IS186）内部的接头，短读长数据无法唯一确定具体的物理整合拷贝（`exact_breakpoint_supported=false`）。 | 严禁采用“最小物理距离就近归因”强行指派拷贝。 |
| **预期编辑完整性审计 (Intended Verification)** | 支持多维状态判定（Complete, Partial, Missing, UnexpectedStructure）；异常靶位结构保留。 | 仅对用户通过 `--intended` 明确声明的序列与结构进行比对审计。 | 预期编辑验证不自动代表菌株表型或代谢功能的合规。 |
| **机制关联与假说生成 (Mechanistic Association)** | 全基因组候选向导依赖位点单次扫描（Single-scan）；输出规范空间几何与失配特征。 | 空间邻近性仅代表几何相关，**绝不代表 Cas 切割因果**。 | 严禁将几何候选标称为已证实脱靶；仅允许统一规范标注为 `candidate guide-dependent off-target`。 |
| **计算性能与吞吐 (Performance & Resource)** | Slurm 计算节点实测核心 evidence 生成耗时 500.28s，较传统 breseq 管线（2648.46s）实现 **5.29x 加速**。 | 内存开销约为 breseq 的 1.84 倍（Peak RSS 3.34 GB vs 1.81 GB），需在具有充裕内存的节点运行。 | 5.29x 加速严格属于 command-level evidence 阶段对比，非端到端产品全流程宣称。 |

---

## 2. 真实数据集验证矩阵与收口边界 (Real-World Benchmarks)

### 2.1 M5: BL21 真实克隆复核 (PRJNA1088182 & GCF_013167015.1)
- **数据集角色**: 单样本参考相对复核与同伴克隆对照矩阵 (`PEER_COMPARATOR`)。
- **已证实能力**:
  - 在 B21_3_1, B21_4_1, B21_2_1 中准确捕获覆盖 lacZ 的参考相对结构异常接头。
  - 六格 PEER_COMPARATOR 对照矩阵（作业 2888235–2888237）全部成功完成，六格产品验证 PASS，M2 判定均为 `MISSING`（预期相关差分事件数为 0，关联数为 0）。方向差分事件数实测如实记录：3→4 为 32，4→3 为 44，2→3 为 54，方向不对称性被如实保留。
- **保留限制与未验证项**:
  - BL21 数据集缺少公开发表的未编辑出发株正交 WGS，非严谨 MATCHED_PARENT。
  - 高拷贝重复序列（IS 家族）未解析到具体物理单拷贝（`exact_breakpoint_supported=false`）。
  - 不证明 Cas9 存在序列特异性切割因果。
  - `MAX_STAGE2_FAMILY_GEOMETRIES=5000` 上限仅在所测 BL21 数据集中未触发，不外推为对所有细菌基因组普遍充裕。

### 2.2 M6: PRJNA1088182 队列复核 (MG1655 & BW25113)
- **数据集角色**: 真实科研编辑菌株队列（Widney & Copley 2024）。
- **已证实能力**:
  - Tier A 参考相对发表突变中，2/12 突变达到共识阈值被严格检出。
  - 其余 10 个位点在读段水平具有低频支持（ALT fraction 0.198–0.691），通过读段级诊断得到如实验证。
  - 实施了 1,319 行双向接头守恒对账（M6-P3c）。
- **保留限制与未验证项**:
  - Tier B 因亲本 WGS 缺失（`NO_PUBLIC_RUN_FOUND`）保持 `NOT_EVALUABLE`。
  - 7 个过滤接头经注释依赖过滤去除，其生物学假阳性状态未经验证。
  - 缺乏三代长读长或正交生物学金标准；登记表 `cohort_25` 行严格保持 `PENDING`。

### 2.3 M7: PRJNA884016 组装真值复核 (PRJNA884016 Assembly Truth)
- **数据集角色**: 基于二代短读长与三代组装真值对比的极端质控。
- **已证实能力**:
  - 染色体组装与参考基因组 100% 严格一致。
  - 识别出历史 4 MOB + 7 JC 属于局部参考截断假象，建立了 fail-closed 隔离保护。
- **保留限制与未验证项**:
  - 该数据集无阳性结构变异，阳性 MOB/JC 灵敏度被评定为 `NOT_EVALUABLE`。
  - 登记表 `prjna884016` 行严格保持 `PENDING`。

### 2.4 M8: REL606 进化克隆单样本对拍与性能 (SRR030257)
- **数据集角色**: breseq 官方经典标准克隆数据集。
- **已证实能力**:
  - 所有 39 个非 JC 突变（28 SNP, 3 INS, 3 DEL, 5 MOB）1-to-1 完全一致（0 over, 0 under, delta=0）。
  - 基于 Slurm 作业 2769249 实测证实 5.29x 计算加速（500.28s vs 2648.46s）。
- **保留限制与未验证项**:
  - 结构接头（JC）严格门保持 `FROZEN_HOLD` (`MATCHED=10, OVER=24, UNDER=8, PASS=false`)。
  - 存量工作区源码 `crates/prokadiff-evidence/src/engine/jc_cluster.rs:1281` 中仍存在历史引入的 `nearer_unique_repeat_representative` 就近距离规则调用。作业 3072550 实测证明其破坏既有对拍（从 10 降到 9），已被 Claude Opus 5.5 科学裁决明确驳回；该路径已被严格冻结，严禁作为生产逻辑推广合并。多拷贝移动元件（IS150/IS186/IS1）代表拷贝选取受限于短读长与序列二义性，需待后续结构变异架构重构解决。
  - 5.29x 性能加速范围严格限定为 command-level evidence 阶段对比，非端到端完整管线声明。

---

## 3. 发布门禁防御性检查清单 (Release Gate Checklist)

在任何正式发布、版本 Tag 或对外报告生成前，必须通过 `scripts/audit_release_gates.py` 的自动化检查：
1. **科学用语红线**:
   - `check_scientific_language`: 0 处违禁因果词汇；向导依赖关联必须附带因果否定限定。
2. **数据链回指红线**:
   - `check_evidence_chain_traceability`: TSV 报告中 100% 具备有效 `DV1 EventId` 或规范位点回指；缺失字段必须以 `NA` 表达，严禁凭空构造。
3. **隔离与冻结完整性**:
   - `check_quarantine_and_hold_integrity`: `VALIDATION_REGISTRY.tsv` 中未完成项目严格保持 `PENDING`；G1 严格门保持 `PASS=false`。
