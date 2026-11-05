<!-- LIVE VERSION MOVED: authoritative HANDOFF.md and ACTIVE_TASK.md now live at the repository root (2026-11-05 controller handoff). This file is historical. -->

# ProkaDiff Development Handoff

## Milestone M9 Formal Sign-Off Closeout & Milestone M10 Activation — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **PREVIOUS_MILESTONE**: `M9` (报告全面审阅与发布门禁) — **COMPLETE (SIGNED_OFF)**
  - Codex Sol completed Round 12 review: **`[PASS] Milestone M9 Accepted and Signed Off — Round 12`**.
  - All findings closed: strict GD mutation type matching, RA/MC exclusion, evidence grammar empty token/duplicate key rejection, documentation test counts synchronized.
  - Final Quality Gates:
    - Release gate test suite: **195 passed in 41.40s** (`test_audit_release_gates.py` on `bnode2`).
    - Full Python test suite: **392 passed in 96.36s** on `bnode2` (58 parity, 195 release gate, 136 real-world, 3 layer-2).
    - Rust workspace test suite: **441 passed, 0 failed, 3 ignored** on `bnode2`.
    - Release gate CLI: `python3 -B scripts/audit_release_gates.py --repo-root .` -> 3 gates, 0 findings, **Overall: PASS**.
- **CURRENT_MILESTONE**: `M10` (多克隆工程批次审计模式 / Multi-Clone Engineering Batch Audit Mode)
- **CURRENT_STATUS**: `M10_T02_READY_FOR_DISPATCH`
  - Completed `M10-T01`:
    - Authored `M10_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md` establishing the one-starter vs many-clones paradigm based on `docs/batch_mode_spec.md` and Taskbook §64–67.
    - Locked the Taskbook §66 Strict Non-Prescriptive Policy (ProkaDiff never ranks or declares a "best" clone).
    - Initialized `M10_TASK_LEDGER.md` and `M10_WORKER_DISPATCH_LEDGER.md`.
    - Updated `STAGE_ROADMAP.md` reflecting M9 completion and M10 activation.
  - Next Action: Dispatch Developer P1 (OpenCode GLM-5.3) for `M10-T02` (`BatchManifest` TSV parser and `--batch-manifest` CLI option).
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`)

---

## M9-T06 Codex Sol Round 11 Review Repair & Round 12 Ready — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M9` (报告全面审阅与发布门禁)
- **CURRENT_STATUS**: `M9_T06_ROUND12_READY`
  - Remediated Codex Sol Round 11 feedback (1 P1 and 1 P2 findings):
    1. **[P1] Strict GD mutation type resolution & evidence-only exclusion**:
       - `GD_MUTATION_TYPES` restricted to standard GenomeDiff mutation types (`SNP`, `SUB`, `DEL`, `INS`, `MOB`, `AMP`, `CON`, `INV`, `JC`), strictly excluding evidence records (`RA`, `MC`).
       - Removed all `("", seq_id, position)` type-free fallback entries. Matching in `_audit_product_table` enforces exact `(decl_type, decl_seq, decl_pos) in gd_mutations`.
       - Intercepts type divergence (e.g. declaring `SNP` at `chr1:100` when sibling `.gd` only has `DEL` at `chr1:100`) and evidence-only records (e.g. sibling `.gd` only has `RA` at `chr1:100`).
    2. **[P1] Strict evidence schema grammar enforcement completed**:
       - Enforced no empty tokens (e.g. consecutive semicolons `RA=1;;MC=NA;JC=NA` or trailing semicolons are rejected).
       - Enforced unique keys: duplicate keys (e.g. `RA=0;RA=1;MC=NA;JC=NA`) are strictly caught as findings.
    3. **[P2] Complete documentation & test suite reconciliation**:
       - Reconciled Section 3.2 in `M9_FINAL_INTEGRATION_REPORT.md` (updated from legacy 379/182 to current 392/195).
       - Synchronized all documents (`M9_FINAL_INTEGRATION_REPORT.md`, `M9_TASK_LEDGER.md`, `M9_WORKER_DISPATCH_LEDGER.md`, `ACTIVE_TASK.md`, `HANDOFF.md`) to exact measured counts:
         - Release gate tests: **195 passed in 41.40s** (`test_audit_release_gates.py` on `bnode2`)
         - Full Python workspace suite: **392 passed in 96.36s** on `bnode2` (58 parity, 195 release gate, 136 real-world cohort & regression, 3 layer-2 comparison)
         - Full Rust workspace suite: **441 passed, 0 failed, 3 ignored** (1 in `cli_e2e`, 2 in `prokadiff_evidence`)
         - Release gate audit CLI: 3 gates, 0 findings, **Overall: PASS**
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`)
- **NEXT_ACTION**:
  - Dispatch Codex Sol (`gpt-6.1-sol`) for M9-T06 Round 12 final independent sign-off review.

---

## M9-T06 Codex Sol Round 10 Review Repair & Round 11 Ready — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M9` (报告全面审阅与发布门禁)
- **CURRENT_STATUS**: `M9_T06_ROUND11_READY`
  - Remediated Codex Sol Round 10 feedback (4 P1 and 1 P2 findings):
    1. **[P1] Evidence grammar vs GD ID confusion resolved**:
       - Replaced ad-hoc integer GD ID extraction from `evidence` with strict schema-compliant evidence grammar parser `_validate_evidence_grammar(evidence) -> tuple[bool, bool, str | None]`.
       - Enforces standard schema `(RA|MC|JC)=(1|0|NA)` per `docs/schema.md:232` and `audit.rs:130`. Intercepts all malformed evidence tokens (`RA=garbage`, `FAKE=yes`, `999`) and all non-positive entries (`RA=NA;MC=NA;JC=NA`).
    2. **[P1] True GD mutation locus resolution**:
       - Implemented `_load_parent_gd_mutations(parent, gd_cache)` resolving actual mutation identities `(gd_type, seq_id, position)` from sibling `.gd` files (`edited.gd` / `output.gd`). If sibling `.gd` files exist, `post_edit_variants.tsv` rows must match actual mutations declared in the GD files.
    3. **[P1] Report finding reconciliation gaps closed**:
       - Removed `#` heading skip; audits table rows and markdown headings (`# Finding DV1_...`) identically.
       - Audits `report.md` unconditionally whenever it exists (including empty/missing declaration tables, flagging orphan DV1 references).
       - Added contig mismatch reporting (`locus_seq != decl["seq_id"]`) and coordinate mismatch reporting (`locus_pos != decl["position"]`).
    4. **[P1] Event content cross-checks expanded**:
       - Added `end` coordinate check.
       - Added `ref` and `alt` checks: for SNPs (and any rows where both are valid DNA letter sequences), `ref` and `alt` must match between `unintended.tsv` and `post_edit_variants.tsv`. Non-DNA representations (e.g. `.` vs `none` on DEL/JC) are respected without false positives.
    5. **[P2] Conflicting test counts in documentation reconciled**:
       - Synchronized all documents (`M9_FINAL_INTEGRATION_REPORT.md`, `M9_TASK_LEDGER.md`, `M9_WORKER_DISPATCH_LEDGER.md`, `ACTIVE_TASK.md`, `HANDOFF.md`) to exact measured counts:
         - Release gate tests: **191 passed in 33.62s** (`test_audit_release_gates.py` on `bnode2`)
         - Full Python workspace suite: **388 passed in 95.54s** on `bnode2` (58 parity, 191 release gate, 136 real-world cohort & regression, 3 layer-2 comparison)
         - Full Rust workspace suite: **441 passed, 0 failed, 3 ignored** (1 in `cli_e2e`, 2 in `prokadiff_evidence`)
         - Release gate audit CLI: 3 gates, 0 findings, **Overall: PASS**
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`)
- **NEXT_ACTION**:
  - Dispatch Codex Sol (`gpt-6.1-sol`) for M9-T06 Round 11 final independent sign-off review.

---

## M9-T06 Codex Sol Round 9 Review Repair & Round 10 Ready — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M9` (报告全面审阅与发布门禁)
- **CURRENT_STATUS**: `M9_T06_ROUND10_READY`
  - Remediated Codex Sol Round 9 feedback:
    1. **[P1] Bidirectional GD-to-report and event content reconciliation**:
       - Replaced loose ID set lookups with rich record dictionaries (`_sibling_event_records`) mapping `event_id` to its declared metadata (`seq_id`, `position`, `end`, `gd_type`, `ref`, `alt`, `evidence`).
       - Implemented coordinate and mutation content consistency verification: sibling tables (`unintended.tsv`, `mutation_offtarget_links.tsv`, `edit_outcomes.tsv`) referencing a declared `event_id` must match `seq_id`, `position`/`mutation_position`, and `gd_type`/`mutation_type`. Mismatched coordinates or types are intercepted as findings.
       - Enforced supporting GD evidence contract: every `event_id` declared in `post_edit_variants.tsv` must contain non-empty, non-NA supporting GD evidence (e.g. `RA=1`, `MC=1`, `JC=1`). Empty, NA, or `RA=NA;MC=NA;JC=NA` without positive evidence are intercepted as findings.
       - Enforced sibling GD record resolution: when sibling `.gd` files are present, GD record IDs referenced in the `evidence` column (e.g. `RA=1`, `JC=2`) must resolve against existing record IDs in the sibling `.gd` files. Unresolved IDs are intercepted as findings.
       - Implemented `report.md` finding reconciliation (`_audit_report_findings`): every DV1 event referenced in `report.md` must resolve against the declared events in `post_edit_variants.tsv` (orphan findings intercepted), and any associated genomic coordinate (e.g. `chr1:100`) must match the declared event coordinates.
    2. **[P2] Handoff and ledger state reconciliation**:
       - Synchronized `docs/development/HANDOFF.md` and root `HANDOFF.md` to document Round 9 blockers and Round 10 remediation.
       - Reconciled all test counts across `M9_FINAL_INTEGRATION_REPORT.md`, `M9_TASK_LEDGER.md`, `M9_WORKER_DISPATCH_LEDGER.md`, and `ACTIVE_TASK.md`:
         - Scanner tests: **182 passed** (`test_audit_release_gates.py` on `bnode2`)
         - Full Python workspace suite: **379 passed** on `bnode2` (58 parity, 182 release gate, 136 real-world cohort & regression, 3 layer-2 comparison)
         - Full Rust workspace suite: **441 passed, 0 failed, 3 ignored**
         - Release gate audit CLI: 3 gates, 0 findings, **Overall: PASS**
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`)
- **NEXT_ACTION**:
  - Dispatch Codex Sol (`gpt-6.1-sol`) for M9-T06 Round 10 final independent sign-off review.

---

## M9-T06 Codex Sol Round 7 Review Repair & Round 8 Ready — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M9` (报告全面审阅与发布门禁)
- **CURRENT_STATUS**: `M9_T06_ROUND8_READY`
  - Remediated Codex Sol Round 7 feedback:
    1. **Comma-separated clause leakage across remainder eliminated**: In `_is_scope_prohibition_governing`, regex for remainder was extended to `r"[,，、.;:!?。！？；：]|\b(?:but|however|yet|nevertheless|nonetheless|although|whereas|though)\b"`, and `claim_items[-1]` is used to extract the nearest governing item. Probes with comma-separated clauses (`Out of scope: any claim that contamination exists, the deletion was-caused-by-Cas9.` and `Out of scope: any claim that contamination exists, this evidence establishes-cleavage-causality.`) as well as line-wrapped variants are strictly intercepted as findings.
    2. **Document & ledger reconciliation**: All test numbers across `M9_FINAL_INTEGRATION_REPORT.md`, `M9_TASK_LEDGER.md`, `M9_WORKER_DISPATCH_LEDGER.md`, `ACTIVE_TASK.md`, and `docs/development/ACTIVE_TASK.md` have been reconciled to the verified results:
       - Scanner tests: **176 passed** (`test_audit_release_gates.py` in 34.92s on `bnode2`)
       - Full Python workspace suite: **373 passed** in 99.34s on `bnode2` (58 parity, 176 release gate, 136 real-world cohort & regression, 3 layer-2 comparison)
       - Full Rust workspace suite: **441 passed, 0 failed, 3 ignored** (1 in `cli_e2e`, 2 in `prokadiff_evidence`)
       - Release gate audit CLI: 3 gates, 0 findings, **Overall: PASS**
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`)
- **NEXT_ACTION**:
  - Dispatch Codex Sol (`gpt-6.1-sol`) for M9-T06 Round 8 final independent sign-off review.

---

## M9-T06 Codex Sol Round 4 Review Repair — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M9` (报告全面审阅与发布门禁)
- **CURRENT_STATUS**: `M9_T06_ROUND4_REPAIRS_APPLIED`
  - Developer P1 (OpenCode / GLM-5.3) 修复 Codex Sol Round 4 审查提出的 2 项科学措辞旁路（反例原文不内嵌于本文件，见 `scripts/test_audit_release_gates.py`）：
    1. **肯定式 claim 框架被当作否定**: `CAUSALITY_FRAMING_NEGATION_RE` 移除裸 `claim that/of/to` 备选（该形式把主语+动词的肯定声明句误判为否定框架）；仅保留真正的否定框架——no claim that / do not (never) claim that / no evidence that / cannot prove that / 禁止类框架 / 中文框架——并显式保留 M5_CLOSEOUT 冻结措辞所用的 `any claim that` 限定枚举框架（新增 `\bany\s+claim\b` 备选）。肯定式声明句现在正确命中 finding。
    2. **无关前缀否定洗白后续肯定子句**: `_is_causality_assertion_negated` 双重收紧——(a) 紧邻动词/形容词否定（`CAUSALITY_IMMEDIATE_NEGATION_RE`）仅在新增 `CLAUSE_DELIMITER_RE`（ASCII/中文逗号与冒号）切分的**紧邻子句**内求值；(b) 声明框架否定必须**无对比屏障**地统辖断言——新增 `CONTRASTIVE_BARRIER_RE`（but/however/yet/nevertheless/nonetheless/although/whereas/though 或逗号/分号），屏障出现在框架匹配与因果断言之间时，框架只统辖前一分句，不再豁免断言（"no evidence of ..., but ..." 与 "cannot prove ..., but ..." 类旁路现在命中 finding）。
  - 措辞对齐（收紧规则的连带整改，语义不变）: `benchmark/real_world/M5_BL21_P1_RECONCILIATION.md` 声明边界表因果行标签由裸 claim 框架措辞改为 `any claim that` 限定枚举措辞（与 M5_CLOSEOUT 冻结措辞同型；三列 NOT SUPPORTED 裁决不变）。收紧后真实仓库全量扫描确认该处是唯一新命中。
  - 测试: `scripts/test_audit_release_gates.py` 新增 10 个用例项（Round 4 四反例 4 项、无逗号对比连词屏障 1 项、跨子句紧邻否定不豁免 1 项、真正的否定框架 no claim that / no evidence that / do not claim that / never claim that 放行 4 项）。
  - 验证证据 (bnode2, Slurm 3095377): pytest **156/156 passed in 17.02s**（原 146 项全数保留通过）；CLI `--repo-root .` **PASS / 0 findings / exit 0**；Round 4 全部 4 个反例经 CLI 端到端复核正确 FAIL（exit 1 且命中对应 finding 类别），干净对照 exit 0。
  - 记录待议（本轮未改动）: `CAUSALITY_ESTABLISHMENT_RES` 中 `establishes?` 字面量仅匹配 establishes（裸不定式 establish 从不进入该模式管控，例如必需 disclaimer 句本身；其放行来自"未匹配"而非紧邻否定）。当前真实仓库零影响（仅有的两处裸 establish 均为紧邻否定的 disclaimer 原句）；是否改写为 `establish(?:es|ed)?` 以覆盖肯定式裸不定式，留待 Controller/Reviewer 裁决。
  - 未提交：改动保留于工作树（scripts/audit_release_gates.py、scripts/test_audit_release_gates.py、benchmark/real_world/M5_BL21_P1_RECONCILIATION.md、docs/development/{ACTIVE_TASK,HANDOFF}.md），由 Controller 决定集成提交。
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`)
- **NEXT_ACTION**:
  - 将 Round 4 修复交付件提交 Codex Sol 进行 M9-T06 Round 5（或最终）独立审查签收。

---

## M9-T06 Codex Sol Round 3 Review Repair — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M9` (报告全面审阅与发布门禁)
- **CURRENT_STATUS**: `M9_T06_ROUND3_REPAIRS_APPLIED`
  - Developer P1 (OpenCode / GLM-5.3) 修复 Codex Sol Round 3 审查提出的 2 项技术缺口（反例原文不内嵌于本文件，见 `scripts/test_audit_release_gates.py`）：
    1. **肯定式因果断言被后缀否定洗白**: `scripts/audit_release_gates.py` 新增 `_is_causality_assertion_negated(line, start, end)` 助手并替换 `DIRECT_CAUSALITY_RES` 与 `CAUSALITY_ESTABLISHMENT_RES` 两处的整句否定搜索。判定收紧为**仅认匹配所在句子内、匹配之前的子句前缀**：(a) 紧邻匹配的动词/形容词否定（`CAUSALITY_IMMEDIATE_NEGATION_RE`：可选助动词 + 否定词 + 有界副词/小品词串，覆盖 was not caused / is not 编辑器名-induced / does not establish / cannot establish / fails to establish / without establishing / never establishes / non- 前缀形容词，及中文"并非由/不建立"类紧邻形式）；(b) 前缀内统辖性声明框架否定（`CAUSALITY_FRAMING_NEGATION_RE`：no claim that / cannot prove (that) / no evidence that / any claim that / 禁止类框架 / 中文框架）。匹配之后的后缀否定（如 "without ambiguity"、"not by contamination"）与其他句子的否定不再豁免；原 `CAUSALITY_SENTENCE_GUARD_RE` 与 `_sentence_of` 移除。
    2. **产品表坐标覆盖缺口**: `COORDINATE_COLUMN_NAMES` 新增 `mutation_position` 与 `side2_position`；`_check_coordinate_cells` 读取行内 `gd_type` 并执行第二侧契约——非 JC 行（`gd_type` 列存在且 != "JC"）允许空白哨兵 `("", ".", "NA")` 放行；任何行上已填充值必须是 `1..MAX_BACTERIAL_COORD` 内正整数（`-1`/`abc`/`10000001` 均触发 finding）；`gd_type == "JC"` 行空白 `side2_position` 触发 finding。
  - 测试: `scripts/test_audit_release_gates.py` 新增 25 个用例项（2 个 Round 3 后缀穿透反例、8 个紧邻前缀否定放行形式、mutation_position 3 反例、side2_position JC 填充 3 反例、非 JC 填充 3 反例、非 JC 空白 3 放行、JC 空白 3 反例）。
  - 措辞对齐（收紧规则的连带整改，语义不变）: `benchmark/real_world/M5_BL21_P1_RECONCILIATION.md` 声明边界表因果行标签改为声明框架措辞（原标签依赖后缀 NOT SUPPORTED 列作为否定，Round 3 规则不再接受）；`docs/development/ACTIVE_TASK.md` M9-T06 条目改写（不再内嵌形容词式触发性英文短语，遵循 Round 2 起的描述惯例）。真实仓库预扫描确认这两处是收紧后仅有的新命中；M5_CLOSEOUT 的声明框架措辞在新规则下保持放行。
  - 验证证据 (bnode2, Slurm 3095377): pytest **146/146 passed in 17.56s**（原 121 项全数保留通过）；CLI `--repo-root .` **PASS / 0 findings / exit 0**；Round 3 全部反例经 CLI 端到端复核正确 FAIL（后缀穿透 2、mutation_position 3、side2 JC/非 JC 填充 6、JC 空白 3），非 JC 空白 3 例与干净对照 exit 0；真实仓库 52 个 unintended.tsv / 43 个 mutation_offtarget_links.tsv 预扫描确认 side2_position（JC 行全部为界内整数、非 JC 行全部空白）与 mutation_position 在新规则下零新命中。
  - 未提交：改动保留于工作树（scripts/audit_release_gates.py、scripts/test_audit_release_gates.py、benchmark/real_world/M5_BL21_P1_RECONCILIATION.md、docs/development/{ACTIVE_TASK,HANDOFF}.md），由 Controller 决定集成提交。
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`)
- **NEXT_ACTION**:
  - 将 Round 3 修复交付件提交 Codex Sol 进行 M9-T06 Round 4（或最终）独立审查签收。

---

## M9-T06 Codex Sol Round 2 Review Repair — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M9` (报告全面审阅与发布门禁)
- **CURRENT_STATUS**: `M9_T06_ROUND2_REPAIRS_APPLIED`
  - Developer P1 (OpenCode / GLM-5.3) 修复 Codex Sol Round 2 对 `scripts/audit_release_gates.py` 的 3 项技术缺口（描述不内嵌触发性英文原句，反例原文见测试用例）：
    1. **因果防护缺口**: 新增形容词式直接因果模式（编辑器名以连字符或空格接 induced 的写法，如 `Cas(?:9|12a)?[- ]induced`，Reviewer 英文反例原句现在命中）；直接因果断言的否定护栏收紧到**匹配所在精确句子**（分隔符 `. ; ! ?` 换行与中文句读），同段落其他句子（反例中前一句“未检出污染”类否定）不再洗白后一句的直接因果断言；因果**确立**类断言不再被文档级 disclaimer 豁免，仅当断言自身句子含否定（句式 "does not establish" 一类）时放行（必需 disclaimer 句子自身因句内否定而通过）。
    2. **坐标与引用校验缺口**: `_check_site_reference` 空白 `site_id` 记为错误；`_check_coordinate_cells` 非数值（负数/字母/空白单元格）记为错误，仅放行干净 `NA` 或 1..10,000,000 内正整数；通用产品表分支（如 `mechanistic_associations.tsv`）现在同样校验坐标列。
    3. **空产品表**: 文件名在 `PRODUCT_FILENAMES` 中但无表头行（0 字节/空）的表记为 EVIDCHAIN 失败，不再静默跳过。
  - 测试: `scripts/test_audit_release_gates.py` 新增/更新 15 个用例（编辑器 induced 形容词的连字符与空格形式、跨句子否定不豁免、disclaimer 加确立断言、否定确立句放行、坐标负数/字母/越界/空白、NA 坐标放行、空白 site_id、0 字节产品表、通用分支坐标校验），并将原“文档级 disclaimer 豁免确立断言”测试反转为必须命中。
  - 验证证据 (bnode2, Slurm 3095377): pytest **121/121 passed in 16.01s**；CLI `--repo-root .` **PASS / 0 findings / exit 0**；Codex Sol Round 2 全部反例经 CLI 端到端复核正确 FAIL（形容词式 1 finding、跨句否定 1、disclaimer 加确立断言 1、坐标与 site_id 组合 2–3、空表 4），干净对照 exit 0。规则收紧前预扫描确认真实仓库 43+9 产品目录与 docs/benchmark 文本在新规则下零新命中。
  - 未提交：改动保留于工作树（scripts/audit_release_gates.py、scripts/test_audit_release_gates.py、docs/development/{ACTIVE_TASK,HANDOFF}.md），由 Controller 决定集成提交。
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`)
- **NEXT_ACTION**:
  - 将 Round 2 修复交付件提交 Codex Sol 进行 M9-T06 Round 3（或最终）独立审查签收。

---

## M9-T06 Codex Sol Review Repair — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M9` (报告全面审阅与发布门禁)
- **CURRENT_STATUS**: `M9_T06_REVIEW_REPAIRS_APPLIED`
  - Developer P1 (OpenCode / GLM-5.3) 全量修复 Codex Sol 对 `scripts/audit_release_gates.py` / `scripts/test_audit_release_gates.py` 的 P1 findings：兄弟表完整性（site/event 孤儿引用、缺表即失败）、坐标合理性边界 (1..10 Mb)、句级直接因果管控（文档级 disclaimer 不再洗白）、策划文档 guide-dependent 讨论强制 disclaimer（生成报告树豁免）、零文件扫描 fail-closed、quarantine 标志非 false 字面量全拒绝、超限 in-scope 文档不再静默跳过。
  - 验证证据 (bnode2, Slurm 3095377): pytest **106/106 passed in 18.71s**；CLI `--repo-root .` **PASS / 0 findings / exit 0**；Codex Sol 四个反例经 CLI 端到端复核全部正确 FAIL，干净对照 exit 0。真实仓库 43+9 个产品目录在新规则下 100% 干净。
  - 文档：`docs/CAPABILITIES_AND_LIMITATIONS_REGISTRY.md` 补充规范因果边界声明。
  - 未提交：改动保留于工作树（scripts/audit_release_gates.py、scripts/test_audit_release_gates.py、docs/CAPABILITIES_AND_LIMITATIONS_REGISTRY.md、docs/development/{ACTIVE_TASK,HANDOFF}.md），由 Controller 决定集成提交。
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`)
- **NEXT_ACTION**:
  - 将修复交付件提交 Codex Sol 进行 M9-T06 Round 2 独立审查签收。

---

## Milestone M8 Formal Closure & Milestone M9 Activation — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M9` (报告全面审阅与发布门禁)
- **CURRENT_STATUS**: `M9_IN_ACTIVE_REPAIR`
  - **Milestone M8 闭环结论**: Codex Sol 独立审查最终签收裁决为 **PASS WITH RECORDED LIMITATIONS**。
    - **Gate A (Rust Engine)**: 441/441 单元与集成测试全绿 (0 failed, 3 ignored on `bnode2`)。
    - **Gate B (Regression Health)**: 197/197 Python 回归与集成测试全绿 (48.01s on `bnode2`)。
    - **Gate C (Non-JC Parity)**: 39/39 非 JC 突变严格 1-to-1 对拍完全一致 (28 SNP, 3 INS, 3 DEL, 5 MOB), 0 over, 0 under, delta=0。
    - **Gate D (Resource Benchmark)**: Slurm 2769249 真实测得 Wall-clock ProkaDiff 500.28s vs breseq 2648.46s (**5.29x 加速**)，Peak RSS 3.34GB vs 1.81GB。
    - **Gate E (G1 Baseline Control)**: 严格门保持 `FROZEN_HOLD` (`MATCHED=10, OVER=24, UNDER=8, PASS=false`)，源码中遗留的距离归因规则严格隔离，严守不修改生产引擎红线。
    - **Gate F (Multi-AI Traceability)**: 独立同行审查全闭环，交付件与台账一致。
  - **Milestone M9 执行状态**: 正在落实 Codex Sol Round 1 审查反馈整改（解决证据链来源回指、行级因果否定校验、fail-closed 完整性保障、能力台账历史数据对齐）。
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **EXECUTION_MODE**: `DELEGATION_FIRST_STRICT`
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **BASE_MAIN_SHA**: `1a0c7b553753bfa1fd6f0c0741c0a8bd43d4cf3c`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`), 8 threads
- **NEXT_ACTION**:
  - 推进 Developer P1 完成发布门禁工具修复与负向测试补充，重新调度 Codex Sol 进行 Round 2 签收。

---

## Milestone M8 Active Integration & Baseline Closure — 2026-11-04

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: `M8` (REL606 引擎扩展对拍、资源测定与基线受控验收)
- **CURRENT_STATUS**: `M8_ACTIVE_INTEGRATION_REVIEW` (Round 2 审查整改已落地)
  - **Gate A (Rust Engine)**: 441/441 单元与集成测试全绿 (0 failed, 3 ignored on `bnode2`).
  - **Gate B (Regression Health)**: 197/197 Python 回归与集成测试全绿 (48.01s on `bnode2`)，`git status --porcelain` 严格守护。
  - **Gate C (Non-JC Parity)**: 39/39 非 JC 突变严格 1-to-1 对拍完全一致 (28 SNP, 3 INS, 3 DEL, 5 MOB), 0 over, 0 under, delta=0，基于 Kuhn 迭代增广路最大二分匹配算法，Codex Sol 审查 PASS。
  - **Gate D (Resource Benchmark)**: Slurm 2769249 真实测得 Wall-clock ProkaDiff 500.28s vs breseq 2648.46s (**5.29x 加速**)，Peak RSS 3.34GB vs 1.81GB，Codex Sol 审查 PASS。
  - **Gate E (G1 Baseline Control)**: 严格门重算一致 (`MATCHED=10, OVER=24, UNDER=8, PASS=false`)。诚实记录为未通过并受控冻结 (`FROZEN_HOLD`)，源码中历史遗留的 `nearer_unique_repeat_representative` 距离归因规则明确披露并隔离，严守不修改生产引擎红线。
  - **Gate F (Multi-AI Traceability)**: Codex Sol Round 1 审查反馈全量闭环整改，台账、路线图与交接文档全面同步。
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **EXECUTION_MODE**: `DELEGATION_FIRST_STRICT`
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **BASE_MAIN_SHA**: `1a0c7b553753bfa1fd6f0c0741c0a8bd43d4cf3c`
- **SLURM_COMPUTE**: Job `3095377`, compute node `bnode2` (`qcpu_18i`), 8 threads
- **NEXT_ACTION**:
  - 推进 Codex Sol 对 M8 整改交付件的 Round 2 签收。
  - 依据 `DEVELOPMENT_PLAN.md` 筹备后续阶段任务。

---

## Multi-AI collaborative strategy migration and fastBact alignment — 2026-11-03

- **CURRENT_PROJECT**: ProkaDiff
- **CURRENT_MILESTONE**: Multi-AI Development Strategy Standardization & Governance Alignment
- **CURRENT_STATUS**: `MULTI_AI_STRATEGY_ALIGNED`
  - Fully aligned with fastBact verified **Multi-AI Collaboration V3.2** (`~/.config/ai-collaboration/POLICY.md`).
  - Cleared 19 legacy root-level orchestration clutter files into `docs/development/archive/legacy_orchestration_202611/`.
  - Terminated legacy Python supervisor loops (`parity_closure_supervisor.py`) and archived stale locks in `benchmark/results/parity_closure/state/`.
  - Refactored `AGENTS.md` to cleanly load `~/.config/ai-collaboration/POLICY.md` without bloat, while strictly preserving ProkaDiff Rules 1–15 and execution constraints.
  - Formulated [`docs/development/MULTI_AI_DEVELOPMENT_STRATEGY.md`](MULTI_AI_DEVELOPMENT_STRATEGY.md) matching fastBact's four-in-one milestone governance standard.
- **CONTROLLER**: Antigravity CLI / Gemini 3.8 Flash Medium (`AUTONOMOUS_EXECUTION_CONTROLLER`)
- **OPERATING_POLICY**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **EXECUTION_MODE**: `DELEGATION_FIRST_STRICT`
- **ACTIVE_WORKTREE**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff`
- **BASE_MAIN_SHA**: `1a0c7b553753bfa1fd6f0c0741c0a8bd43d4cf3c`
- **ROLE_BINDINGS**:
  - Controller: Antigravity CLI / Gemini 3.8 Flash Medium
  - Developer P1: OpenCode CLI (`glm-5.3-flash` routine, `glm-5.3` complex)
  - Developer P2: Cursor Pro (`grok-4.7-high-fast`, `grok-4.7-xhigh-fast`) as backup
  - Independent Reviewer: Codex CLI (`gpt-6.1-terra` / `gpt-6.1-sol`)
  - Validator: Codex CLI (`gpt-6-luna`) + deterministic tools on Slurm compute nodes
  - Architect: Codex Astra (`gpt-6-astra`) / Claude Opus 5.5 (`claude-opus-5-5`)
  - Final Authorization: Human Owner (`FINAL_AUTHORIZATION_AUTHORITY`)
- **SCIENTIFIC_STATUS_AND_INVARIANTS**:
  - M7 status: `COMPLETE_WITH_LIMITATIONS`.
  - PRJNA884016 positive truth remains strictly quarantined.
  - G1 Clonal JC parity closure: Baseline Slurm `3066886` (matched 10, OVER 24, UNDER 8). `WAVE11E_HOLD.md` remains in force. No uncommitted/dirty production overrides allowed without approved contract.
  - Rules 1–15 (Path B bowtie2 alignment, Rust evidence extraction, clean-room vs breseq, noodles BAM I/O, classification order) fully active.
- **IMMEDIATE_ACTION**:
  - Multi-AI continuous development strategy finalized and operational.
  - Ready for active milestone dispatch via OpenCode / GLM under Antigravity Controller.
- **NEXT_ACTION**:
  - Await user guidance on next specific milestone or technical task (e.g. M8 milestone progression or G1 clonal JC parity resolution).

---

## Historical handoff archive (M7 closure & G1 qualification) — 2026-11-03T15:06

- **Date**: 2026-11-03T15:06+08:00
- **Context**: Transition from experimental supervisor loops and local failover back to standard Multi-AI framework.
- **Key Evidence**:
  - Codex Sol qualification verdict: `PASS` on `wave11b2/WAVE11B2_QUALIFICATION.md`.
  - Claude Opus 5.5 scientific architect consultation: `WAVE11B1_VERDICT: INSUFFICIENT_EVIDENCE`, `PRODUCTION_CHANGE: NO`. Verified spend 0.2611 USD; cumulative spend 0.7004 USD of 5.00 USD cap.
  - Construct delta documented in `wave11b2/WAVE11B2_CONSTRUCT_DELTA.md`.
  - Regression verified: 19 tests passed in `test_m7_p0b_scoring.py` on Slurm Job `3095377`.
