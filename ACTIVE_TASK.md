# Active Task

- **Project**: ProkaDiff
- **Milestone**: `M11` (亚克隆低频变异与混群审计模式 / Subclonal Pool-seq Variant Audit Mode)
- **Active Task ID**: `M11-T06` (阶段集成与独立审查验收)
- **Controller**: OpenCode CLI / GLM-5.3-Flash (acting Controller per user-authorized handoff 2026-11-05; former Controller Antigravity / Gemini 3.8 Flash Medium weekly quota exhausted)
- **Reviewer**: OpenAI Codex CLI / `gpt-6.1-sol` (Round 3 pending on current snapshot)
- **Developer P1**: OpenCode CLI / GLM-5.3
- **Governing Policy**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **Slurm Allocation**: Job `3095377`, compute node `bnode2` (`qcpu_18i`), reused via `srun --overlap`
- **Status**: `RELEASED v0.11.0-m11` — M11 accepted (Codex Sol R6 `[PASS]`) and released by user authorization 2026-11-05: `main` pushed to `origin` (`1a0c7b5..da58fe0`), annotated tag `v0.11.0-m11` pushed; GitHub Release object to be published from the tag by the user (API credentials unavailable in dev environment; notes prepared at `.prokadiff/RELEASE_NOTES_v0.11.0-m11.md`). Next milestone requires fresh user authorization.

---

## M11 正式验收 (2026-11-05) — Codex Sol Round 6 `[PASS] Accepted`

- **Final verdict**: `[PASS] Milestone M11 Round 6 Accepted`, no open P1/P2 findings (`.prokadiff/codex_m11_r6_review.txt`).
- **Gates**: A–E all PASS (Gate A under user-authorized revision, contract §7).
- **Final evidence**: workspace 494 passed / 0 failed / 3 ignored; clippy 0 warnings; 256 oracle tests green; release gates 3×PASS — qfree vnode1 jobs 3130686/3130717.
- **Baseline commit**: `b30fbad` (user-authorized, local only; push/tag/release retained by user).
- **Acceptance report**: [`M11_FINAL_INTEGRATION_REPORT.md`](M11_FINAL_INTEGRATION_REPORT.md).
- **Review trail**: R1–R6 (REJECTED ×5 with all findings repaired within the three-attempt rule; ACCEPTED ×1). Full trail in `M11_WORKER_DISPATCH_LEDGER.md` DISP-M11-006…017.
- **Next action (user)**: decide on push/tag/release and authorize the next roadmap sprint (`STAGE_ROADMAP.md`). No autonomous dispatch remains for M11 per POLICY §13.4 (terminal state reached).

---

## M11 当前状态快照 (2026-11-05T11:01+08:00)

1. **M11-T01…T05**: COMPLETE (see `M11_TASK_LEDGER.md`).
2. **M11-T06**: IN_PROGRESS.
   - Codex Sol Round 1: REJECTED (9 findings, 09:43) → repaired.
   - Codex Sol Round 2: REJECTED (4×P1 + 2×P2, 10:19, `.prokadiff/codex_m11_r2_review.txt`).
   - Post-R2 repair round completed 10:20–10:27 (6 files; all 6 findings code-verified fixed; snapshot unreviewed at Gemini quota exhaustion).
   - Takeover recovery validation (2026-11-05): `cargo test --workspace` **481 passed / 0 failed / 3 ignored** on bnode2 (alloc 3095377); clippy 0 warnings; release gates 3×PASS / 0 findings.
   - Next: Codex Sol Round 3 review of current snapshot.
3. **Preserved constraints**: default consensus mode byte-for-byte unchanged; M11 low-frequency scope = SNP/SUB/small INS/DEL (≤2 bp) only; no low-frequency MOB/JC/AMP; no 5% sensitivity claims (AF 0.05 is a configurable threshold, not an experimental sensitivity floor); titration tests use pseudorandom **non-repeating** synthetic references — documented as synthetic validation only, not evidence for repetitive genomic regions; G1 `FROZEN_HOLD` and Registry `PENDING` states untouched.

---

## 历史里程碑收口 (Milestone M10 Closeout, 归档)

- **M9 状态**: **COMPLETE (SIGNED_OFF)** (2026-11-04)
- **审查裁决**: Codex Sol 于 Round 12 出具正式审查裁决：`[PASS] Milestone M9 Accepted and Signed Off — Round 12`。

---

## 当前里程碑目标与进展 (Milestone M10 Progress)

1. **`M10-T01` (M10 契约制定与多克隆数据模型确立)**: **COMPLETE**.
   - 制定了 `M10_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`，基于 `docs/batch_mode_spec.md` 确立了一祖先对多克隆架构、`RecurrenceStatus` 与 `BatchVariantRecord` 模型，明确了 Taskbook §66 严格非推荐性原则（严禁评选“最佳克隆”）。
   - 初始化了 `M10_TASK_LEDGER.md` 与 `M10_WORKER_DISPATCH_LEDGER.md`。
2. **`M10-T02` (Batch Manifest 解析器与 CLI 参数拓展)**: **COMPLETE**.
   - 在 `crates/prokadiff-report/src/batch/manifest.rs` 实现 `BatchManifest` 解析校验器（单端/双端、每克隆可选 spacer/intended 覆盖、防路径穿越 `validate_clone_id`、防破行 `RaggedRow`、预检文件存在性扫描）。
   - 在 `crates/prokadiff/src/cli.rs` 增加 `--batch-manifest <path>` 选项，严格执行 Rule 6（缺失 `--starter` 立即非零报错）。
3. **`M10-T03` (多克隆变异对齐与复现性分类器)**: **COMPLETE**.
   - 在 `crates/prokadiff-report/src/batch/recurrence.rs` 实现 `RecurrenceStatus` 与 `align_batch_variants`。
   - 跨克隆对齐以 `CanonicalEvent::to_preimage_bytes()` 规范线协议哈希作为同一性依据，保证端点互换 JC（`chr1:100` -> `chr2:200` 与 `chr2:200` -> `chr1:100`）精确合并，且区分不同 CON 供体区域。
4. **`M10-T04` (批次报告与变异矩阵生成器)**: **COMPLETE**.
   - 在 `crates/prokadiff-report/src/batch/render.rs` 实现 `batch_report.md` 与 `batch_variants.tsv` 渲染器。
   - 严格遵循 Taskbook §66 非推荐性原则，明确输出观测指标与重现性统计，严禁评价或推荐任何克隆为 "best", "safest", "recommended"。
   - 区分预期状态：`Complete (k/N)`、`Partial (k/N)`、`None observed (0/N)`、`Not specified`。
5. **`M10-T05` (多克隆批次模式端到端集成测试)**: **COMPLETE**.
   - 在 `crates/prokadiff/tests/batch_e2e.rs` 构建 10 项端到端测试，覆盖 CLI 参数、Rule 6 约束、单组件 clone_id 校验、预检拦截、Cas9 guide 传递与缺失拦截、完整合成队列（含 Bowtie2 比对）运行。
   - 10/10 测试全部通过 (48.61s, `bnode2`)。
   - 在 `crates/prokadiff-report/src/batch/render.rs:628` 补充 `test_extract_clone_observations_canonical_jc_swap_and_con` 真实 fixture 测试，验证规范化线协议合并与区分。
6. **`M10-T06` (M10 阶段集成与独立审查验收)**: **COMPLETE (ACCEPTED)**.
   - 质量门禁全绿：
     - `scripts/audit_release_gates.py`: **3 gates PASS, 0 findings** (发布门禁全绿)。
     - `cargo clippy --workspace --all-targets -- -D warnings`: **0 warnings, PASS**.
     - `cargo test --workspace`: **463 passed, 0 failed, 3 ignored** (compute node `bnode2`, Job `3095377`).
   - 证据日志归档：`.prokadiff/m10_workspace_test.log`, `.prokadiff/m10_batch_e2e.log`, `.prokadiff/m10_report_test.log`。
   - 输出正式验收报告：[`M10_FINAL_INTEGRATION_REPORT.md`](M10_FINAL_INTEGRATION_REPORT.md)。
   - 更新台账：`M10_TASK_LEDGER.md` (全部 COMPLETE), `M10_WORKER_DISPATCH_LEDGER.md` (ACCEPTED)。
