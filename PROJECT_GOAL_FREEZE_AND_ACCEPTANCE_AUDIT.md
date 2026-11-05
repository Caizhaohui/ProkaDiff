# ProkaDiff — Project Goal Freeze & Final Acceptance Audit

**Document ID:** `PROJECT_GOAL_FREEZE_AND_ACCEPTANCE_AUDIT.md`  
**Status:** `GOAL_FROZEN` · **Project verdict:** **`PROJECT_ACCEPTANCE_PENDING`**  
**Audit date:** 2026-10-28  
**Controller:** Cursor `composer-2.5` (read-only audit; no silent rewrite of M7 closeout)  
**Oracle pin:** breseq **0.40.2** + bowtie2 **2.5.4** (`testdata/VERSIONS.txt`)

---

## 1. Frozen project goal (binding)

ProkaDiff shall be a high-performance, **breseq-compatible** prokaryotic whole-genome variant caller and annotator, and on that foundation support post-edit genome audit.

**Computational gold standard for this stage:** breseq **0.40.2** (not independent biological truth authority).

**PROJECT_SUCCESS requires all three gates:**

| Gate | Meaning |
| --- | --- |
| `VARIANT_CALL_PARITY_PASS` | Accepted mutation calls: bidirectional non-whitelist over/under = 0 vs breseq on frozen workloads |
| `VARIANT_ANNOTATION_PARITY_PASS` | Core annotation fields with shared semantics match where both tools emit them |
| `PERFORMANCE_ADVANTAGE_PASS` | Same inputs/ref/threads/task scope: ProkaDiff wall-clock faster than breseq, with measured RSS/versions/SHA; no accuracy sacrifice |

**Explicit non-equivalence:** unit-test pass, `M7_COMPLETE_WITH_LIMITATIONS`, Multi-AI review ACCEPT/PASS ≠ `PROJECT_SUCCESS`.

**Off-target / audit dependency order (frozen):**

```text
Whole-genome Variant Calling
→ breseq Parity Validation
→ Starter/Edited Differential Subtraction
→ Intended Edit Assessment
→ Candidate Off-target Association
→ Final Genome Audit Report
```

Until upstream WGS calling passes the corresponding breseq acceptance gate: do not claim formal off-target acceptance; do not reverse-construct observed variants from predicted sites; do not treat candidate associations as proven editing causality; do not expand off-target algorithms to paper over caller gaps. **Keep existing off-target modules; no new algorithms/modes.**

**Allowed development (only):** fix call/annotation parity gaps; real algorithm defects affecting concordance; prevent positive-event undercalls; necessary regressions; evidence-backed performance work; reproducible final benchmarks.

**Paused:** new variant types, new sequencing platforms, GUI, new viz systems, off-target algorithm expansion, unnecessary large refactors, unrelated features.

**M7 historical closeout must not be silently rewritten.** New fixes use a separate bounded **parity-closure** ledger (this document + future `PARITY_CLOSURE_*.md` entries). Do not auto-start M8. Do not auto-push GitHub.

---

## 2. Audit scope & evidence anchors (actual SHA / paths)

| Tree | Path | HEAD (verified) |
| --- | --- | --- |
| Main worktree | `04_ProkaDiff` | `1a0c7b5` (`fix(benchmark): fail closed when PRJNA884016… M7-P0b`) — **does not contain M7 closeout commit** |
| M7 worktree | `04_ProkaDiff-m7-worktree` branch `m7-autonomous-dev` | **`04ee053`** (`docs(m7): close M7 with limitations after P2 oracle completion`) |
| M7 closeout | `…/M7_CLOSEOUT.md` | Status **`M7_COMPLETE_WITH_LIMITATIONS`**; internal §9 still cites reproduction HEAD `e2eb617` (doc drift vs `04ee053` — note only, do not rewrite history) |
| Parity contract | `docs/parity.md` | Retained as verified comparison rules |
| Publishable perf | `benchmark/results/verified/clonal_{2424014,2433013,2673729}.md` | Same-job wall/RSS + comparator exit |
| Latest clonal notes | `benchmark/results/clonal_2673729/parity_notes.txt` | Authoritative numeric JC/red leftovers |

---

## 3. Implemented variant surface (code + measured GD)

### 3.1 Emitted / first-line accepted mutation types

| Type | Engine emit | GD I/O | Measured Clonal `2673729` ProkaDiff | Same-node breseq |
| --- | --- | --- | --- | --- |
| SNP | yes | yes | 28 | 28 |
| INS | yes | yes | 3 | 3 |
| DEL | yes | yes | 3 | 3 |
| SUB | RA path / GD | yes | 0 | 0 (this run) |
| MOB | yes | yes | 5 | 5 |
| JC | yes | yes | 33 | 27 |
| AMP | **not emitted** | I/O only | 0 | 0 |
| CON | **not emitted** | I/O only | 0 | 0 |
| INV | I/O only | I/O only | 0 | — |
| RA / MC | evidence | yes | MC 275 | RA 45 / MC 4 |
| UN | not ProkaDiff product | — | — | 789 |

Evidence-only / diagnostic: `repeat_ambiguous_junctions.tsv` (M5-SR1) — not GD mutations, no EventId.

### 3.2 Product / audit stack (downstream of caller)

M0–M6 closed in plan terms (M6 = `CLOSED_WITH_LIMITATIONS`). Pipeline modules exist: evidence → differential (`EventId`) → intended → association → `AuditResult` / report. **Downstream audit is not a substitute for call/annotation parity.**

---

## 4. Gate-by-gate verdict

### 4.1 `VARIANT_CALL_PARITY_PASS` → **FAIL (not met)**

**What is evidenced as good**

- Clonal same-node job **`2673729`**: red-line `SNP/INS/DEL/MOB/AMP/CON` bidirectional leftover **`over_red=0` / `under_red=0`** (`parity_notes.txt`); MOB 5/5 present in both GDs.
- Synth layer-2 history (DEL / SNP-indel / IS-MOB) has red-line / JC-tol successes recorded in `VALIDATION_REPORT.md` / job notes (not re-litigated here).
- `docs/parity.md` whitelist retained (UN, HTML, BAM names, repeat-ambiguous JC ±tol, frequency rounding; RA short indel 3′ alignment required).

**What blocks PASS**

1. **Clonal JC bidirectional non-empty (same-node breseq 0.40.2, job 2673729):**  
   `over_jc_unmatched=29`, `under_jc_unmatched=23`, of which `under_jc_oracle_rejected=9`, `jc_matched_tol=4`.  
   Comparator `scripts/layer2_gd_compare.py` **exits 0 on red-line only** — JC leftovers do **not** fail the script. README “0 false positives” refers to **`over_red`**, not full GD / JC empty.  
   → Formal accepted-call parity (including JC as accepted mutations when not rejected/whitelisted) is **not** satisfied.
2. **BL21 real WGS vs breseq (registry ACTIVE but JC not validated):** after depth-aware filter, residual JC FP ~29–38 and recall ~0.40–0.60 (`VALIDATION_REGISTRY` / RW-002 RESOLVED only as FP *reduction*, not parity). `mob_validated=false`, `jc_validated=false`.
3. **M7 PRJNA884016:** P1b **`FAILED_STRICT_ZERO_GATE`** (27 unintended JC vs truth `NONE`); ProkaDiff CALL_BURDEN starter JC41 / edited JC32 vs P2 breseq oracle JC9 / JC12. Positive structural truth **`NOT_EVALUABLE`**. M7 ≠ call-parity PASS.
4. **AMP/CON:** no engine emission path; any future oracle AMP/CON would be undercalls (currently absent on Clonal).
5. **Stale `parity_summary.tsv`** (commit `b541121`) shows clonal `jc_over=33` / `jc_under=27` — must not be cited as current PASS; use job `parity_notes.txt` + verified/*.md.

### 4.2 `VARIANT_ANNOTATION_PARITY_PASS` → **FAIL (not met)**

**Evidence**

- breseq Clonal `annotated.gd` SNP lines carry shared-semantic fields (`gene_name`, `locus_tag`, `aa_*`, `codon_*`, `mutation_category`, …).
- ProkaDiff Clonal `output.gd` SNP lines are **identity-only** (no `gene=` / aa / codon attrs on measured sample).
- Product TSV gene/locus annotation (`prokadiff-classify` / report) is **GenBank overlap projection**, not a verified breseq Annotate field-by-field compare harness.
- `layer2_gd_compare.py` explicitly: **“annotations ignored”** for identity.
- No repository artifact named / gated as annotation-parity PASS vs breseq 0.40.2.

**Scope note (frozen):** compare only core fields with corresponding semantics on both sides; do not require breseq HTML/R plots (`docs/parity.md` whitelist). Still: **no verified core-field annotation concordance → gate fails.**

### 4.3 `PERFORMANCE_ADVANTAGE_PASS` → **PARTIAL (Clonal wall-clock only; not project PASS)**

| Job | Node | Threads | ProkaDiff wall_s | breseq wall_s | ProkaDiff RSS_kB | breseq RSS_kB |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 2424014 | bnode29 | 8 | 544.55 | 3098.01 | 3432960 | 1856948 |
| 2433013 | bnode2 | 8 | 500.04 | 2674.81 | 3428164 | 1856884 |
| 2673729 | bnode1 | 8 | 488.61 | 2640.40 | 3420124 | 1856968 |

- **Wall-clock:** ProkaDiff faster on frozen Clonal workload (≈5.1–5.7× on these three same-job records).  
- **RSS:** ProkaDiff ≈1.8× higher peak RSS — recorded, not hidden.  
- **Publishable set:** only Clonal under `benchmark/results/verified/`. Do **not** promote a single-job multiplier as “overall speedup.”  
- **Accuracy coupling:** user freeze forbids claiming performance advantage obtained by skipping necessary concordance. With JC call leftovers and missing annotation parity, **project-level `PERFORMANCE_ADVANTAGE_PASS` is not awarded** even though Clonal wall-clock evidence exists.

---

## 5. Dataset gates (actual, not README)

| Dataset | Registry / closeout | Call vs breseq | Annotation vs breseq | Perf verified |
| --- | --- | --- | --- | --- |
| Synth SNP/INS/DEL, DEL, IS-MOB | layer-2 jobs / VALIDATION_REPORT | Red-line historically OK; JC tol case-by-case | Not gated | Incomplete / not in `verified/` |
| REL606 Clonal (SRR030257) | `rel606` ACTIVE; oracle flags true in registry | Red-line OK @2673729; **JC leftovers** | **Not shown** | **Yes** (verified/*.md) |
| BL21 user clones | ACTIVE; intended OK on B21_3_1; mob/jc **false** | Residual JC FP / FN | Product only | No same-job breseq wall in `verified/` |
| PRJNA1088182 | `cohort_25` **PENDING**; M6 `CLOSED_WITH_LIMITATIONS` | Publication recovery 2/12; Tier B NOT_EVALUABLE | Repeat-filter A/B only | Not project perf gate |
| PRJNA884016 | **PENDING**; M7 `COMPLETE_WITH_LIMITATIONS` | P1b zero-gate **FAILED**; metrics NOT_EVALUABLE | N/A | N/A |

---

## 6. Severity-ranked gap list (finite; blockers first)

### P0 — Blocks `VARIANT_CALL_PARITY_PASS`

| ID | Gap | Evidence | Allowed fix class |
| --- | --- | --- | --- |
| G1 | Clonal accepted JC over/under vs same-node breseq | 2673729: 29 over / 23 under (9 oracle-rejected) | Algorithm / accept-filter / reject semantics; regressions; **no whitelist expansion without contract** |
| G2 | BL21 residual JC FP/FN vs breseq | Registry jc_validated=false; FP≈29–38 post-filter | Same as G1; prevent TP loss |
| G3 | M7 differential JC burden vs NONE / vs breseq oracle scale | P1b 27 unintended JC; PD JC41/32 vs breseq 9/12 | Caller concordance; do not weaken zero-gate |

### P1 — Blocks `VARIANT_ANNOTATION_PARITY_PASS`

| ID | Gap | Evidence | Allowed fix class |
| --- | --- | --- | --- |
| G4 | No breseq-core annotation field parity harness / GD gene·aa·codon attrs | Clonal GD: breseq annotated vs PD bare SNP | Annotate core shared fields; field-limited compare; regressions |

### P2 — Blocks project-level `PERFORMANCE_ADVANTAGE_PASS` / evidence hygiene

| ID | Gap | Evidence | Allowed fix class |
| --- | --- | --- | --- |
| G5 | Perf PASS cannot close while G1–G4 open | Freeze §6 accuracy coupling | After call+annotation gates: refresh same-job verified records |
| G6 | `parity_summary.tsv` / README overclaim risk | Stale TSV; README cites over_red as “0 FP” | Doc/evidence hygiene only after metrics re-verified |

### P3 — Open / deferred (do not expand scope)

| ID | Gap | Note |
| --- | --- | --- |
| G7 | RW-003 MC vs breseq UN metric fairness | OPEN; engine UN-equivalent deferred |
| G8 | AMP/CON/INV emission | Not required until oracle emits them on frozen workloads |
| G9 | Off-target score oracles (FlashFry/CRISPRitz) | BLOCKED_EXT_DEP; **paused** under freeze |
| G10 | M8 “REL606 engine extension” plan row | Superseded by this freeze as *parity closure*, not new milestone marketing |
| G11 | GUI / platforms / CAST-IS110 / long-read | **Paused** |

---

## 7. Unnecessary / out-of-scope development (pause)

- New variant types, platforms, GUI, visualization systems  
- Off-target algorithm expansion or new association modes  
- Expanding scientific truth / registry flags because Slurm completed  
- Restoring quarantined PRJNA884016 4 MOB+7 JC  
- Auto M8 start, push/merge/tag without owner authorization  
- Silent edits to `M7_CLOSEOUT.md` historical facts  

---

## 8. Final project verdict

| Gate | Result |
| --- | --- |
| `VARIANT_CALL_PARITY_PASS` | **FAIL** |
| `VARIANT_ANNOTATION_PARITY_PASS` | **FAIL** |
| `PERFORMANCE_ADVANTAGE_PASS` | **NOT AWARDED** (Clonal wall-clock evidence exists; project gate blocked by accuracy coupling + incomplete multi-workload verified set) |

### **`PROJECT_ACCEPTANCE_PENDING`**

Not `PROJECT_SUCCESS`. Not declared `PROJECT_BLOCKED` (gaps are finite and in-scope for parity closure).

**Next engineering (only):** close G1→G3 (call), then G4 (annotation), then refresh verified same-job benchmarks for G5 — each under a bounded parity-closure record. Stop feature development for this stage upon three-gate PASS; then freeze version + docs only.

---

## 9. Principle (frozen)

> First guarantee whole-genome variant **calling** and **annotation** concordance with breseq; only then treat gene-edit off-target association as formally acceptable; and prove ProkaDiff runs **faster** on frozen comparable workloads — then end this stage’s software development.
