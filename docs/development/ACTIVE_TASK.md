<!-- LIVE VERSION MOVED: authoritative HANDOFF.md and ACTIVE_TASK.md now live at the repository root (2026-11-05 controller handoff). This file is historical. -->

# Active Task

- **Project**: ProkaDiff
- **Milestone**: `M10` (多克隆工程批次审计模式)
- **Active Task ID**: `M10-T02`
- **Controller**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **Reviewer**: OpenAI Codex CLI / `gpt-6.1-sol`
- **Developer P1**: OpenCode CLI / GLM-5.3
- **Governing Policy**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **Slurm Allocation**: Job `3095377`, compute node `bnode2` (`qcpu_18i`), 8 threads
- **Status**: `M10_T02_READY_FOR_DISPATCH`

## 历史里程碑收口 (Milestone M9 Closeout)

- **M9 状态**: **COMPLETE (SIGNED_OFF)** (2026-11-04)
- **审查裁决**: Codex Sol 于 Round 12 出具正式审查裁决：`[PASS] Milestone M9 Accepted and Signed Off — Round 12`。
- **质量门禁汇总**:
  - `scripts/audit_release_gates.py`: 3 gates, 0 findings, **Overall: PASS** (科学用语因果越界防护、双向证据回指溯源、隔离区完整性校验均通过)。
  - `scripts/test_audit_release_gates.py`: **195 passed in 41.40s** (compute node `bnode2`).
  - Python 全工作区套件: **392 passed in 96.36s** (compute node `bnode2`).
  - Rust 全工作区套件: **441 passed, 0 failed, 3 ignored** (compute node `bnode2`).

---

## 当前里程碑目标与进展 (Milestone M10 Progress)

1. **`M10-T01` (M10 契约制定与多克隆数据模型确立)**: **COMPLETE**.
   - 制定了 `M10_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`，基于 `docs/batch_mode_spec.md` 确立了一祖先对多克隆架构、`RecurrenceStatus` 与 `BatchVariantRecord` 模型，明确了 Taskbook §66 严格非推荐性原则（严禁评选“最佳克隆”）。
   - 初始化了 `M10_TASK_LEDGER.md` 与 `M10_WORKER_DISPATCH_LEDGER.md`。
2. **`M10-T02` (Batch Manifest 解析器与 CLI 参数拓展)**: **READY**.
   - 目标：在 `crates/prokadiff-report` (或对应模块) 中实现 `BatchManifest` 解析与校验器，在 `crates/prokadiff/src/cli.rs` 中增加 `--batch-manifest <path>` 选项。
   - 校验要求：严格检查克隆 ID 唯一性、FASTQ 路径存在性、单端/双端支持、每克隆可选 spacer/intended 覆盖；缺少公共 starter 或 ref 时 fail-closed 报错。
