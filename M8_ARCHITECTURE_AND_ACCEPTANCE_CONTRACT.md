# Milestone M8: Architecture and Acceptance Contract
## REL606 引擎扩展对拍与性能基准阶段契约

- **Milestone ID**: `M8_REL606_ENGINE_PARITY_AND_PERFORMANCE`
- **Controller**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **Operating Policy**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **Strategy Document**: [`docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md`](docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md)
- **Status**: `M8_CONTRACT_FORMALIZED`
- **Date**: 2026-11-03

---

## 1. 阶段目标与适用范围 (Objectives & Scope)

### 1.1 核心目标
1. **REL606 单样本引擎对拍基准固化**:
   - 在真实大肠杆菌进化克隆数据集 `REL606` (`SRR030257`) 上，系统性验证 ProkaDiff 对比 breseq 0.40.2 oracle 的突变召回与精度。
   - 证明所有非 JC 共识突变类型（SNP, INS, DEL, MOB）达到 100% 严格一致（0 over, 0 under），保证变异引擎核心逻辑的高度准确。
2. **Wall-clock 耗时与 Peak RSS 资源测定**:
   - 依据 Rule 5 与 Rule 15，在 Slurm 计算节点（`qcpu_18i`）上，独立测量并解耦记录 ProkaDiff 与 breseq 在相同数据集、相同 CPU 资源下的真实 Wall-clock 耗时与内存峰值（Peak RSS），输出机器可读的性能对比表。
3. **G1 结构接头 (JC) 科学现状受控管理**:
   - 严格遵循人类用户关于 G1 严格门（`matched=10, OVER=24, UNDER=8, PASS=false`）的不可变指令，不使用未经科学证明的就近距离归因规则，不改动生产引擎，不提交未授权作业，完整保留归因证据。
4. **全套交付件规范收口**:
   - 消除脆弱的回归测试断言，确保全工作区单元测试与全套 Python 回归测试在计算节点上 100% 绿色通过，形成完整的阶段集成审计报告。

---

## 2. 核心技术不变量与科学红线 (Invariants & Red Lines)

1. **Path B 引擎与 Clean-room 原则 (Rule 4, Rule 5)**:
   - 运行时严禁调用 breseq 二进制；breseq 仅作为外部对照 Oracle。
   - 严禁复制 breseq (GPL-2.0) 源代码，所有实现均基于公开算法原理自研。
2. **非 JC 变异零容忍红线**:
   - 在消除 Contig 版本命名差异（`REL606.6` $\leftrightarrow$ `REL606`）后，`SNP`, `INS`, `DEL`, `MOB` 必须保持 0 over-call 与 0 under-call。
3. **G1 严格门受控冻结**:
   - 严格禁止在生产分支上为了迎合对拍而引入“就近距离归因规则”。
   - 在缺乏充分生物学与序列证据前，G1 状态严格保持 `BLOCKED_AUTHORIZATION`，未解决的 OVER/UNDER 必须在报告中作为已知限制如实记录。
4. **Slurm 执行与登录节点保护 (Rule 11)**:
   - 严禁在集群登录节点运行任何编译、测试或计算密集型作业。单元测试与快速集成在 `qfree` 执行，全基因组对齐与 breseq 对照在 `qcpu_18i` 执行。
5. **Git 与工作树安全原则**:
   - 严格单一写入者，禁止执行 `git clean -fd` 或 `git reset --hard`。
   - 生产分支保持稳定，所有发布和合并授权保留给人类用户。

---

## 3. 验收门禁矩阵 (Mandatory Acceptance Gates)

| 门禁 ID | 门禁名称 | 验收标准 | 验证工具 / 证据 |
| --- | --- | --- | --- |
| **Gate A** | **全工作区 Rust 单元测试** | `cargo test --workspace` 441 项测试 100% 通过（0 failed, 3 ignored）。 | `bnode2` (Slurm) 执行记录 (task-1249) |
| **Gate B** | **Python 回归测试套件全绿** | `pytest scripts/ benchmark/real_world/scripts/` 197 项测试 100% 通过。 | Python 3.12 虚拟环境测试日志 (task-1243) |
| **Gate C** | **REL606 非 JC 突变对拍闭环** | 39/39 非 JC 核心突变 1-to-1 完全匹配，0 over, 0 under。 | `scripts/layer2_gd_compare.py` / TSV 证据 |
| **Gate D** | **独立性能与资源测量台账** | ProkaDiff 与 breseq 耗时与 Peak RSS 独立分离记录，禁止捏造数字。 | `m8_perf_summary.tsv` / Slurm `sacct` / `/usr/bin/time -v` |
| **Gate E** | **G1 现状受控审计与限制保留** | 严格门状态、差异分类、Hold 证据完整记录，无越权修改。 | `WAVE11E_HOLD.md`, `jc_acceptance_gate.py` |
| **Gate F** | **独立同行审查与最终集成报告** | Codex Sol 独立审查通过，形成完整的 M8 交付报告。 | `M8_FINAL_INTEGRATION_REPORT.md` |

---

## 4. 交付件清单 (Deliverables)

1. `STAGE_ROADMAP.md`
2. `M8_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`
3. `M8_TASK_LEDGER.md`
4. `M8_WORKER_DISPATCH_LEDGER.md`
5. `benchmark/results/m8_rel606_non_jc_parity.tsv`
6. `benchmark/results/m8_perf_summary.tsv`
7. `M8_FINAL_INTEGRATION_REPORT.md`
8. `docs/development/HANDOFF.md` 更新
