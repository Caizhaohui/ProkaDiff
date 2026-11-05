# M7 Multi-AI dispatch acceptance (V1.4)

**Report date:** 2026-10-28  
**Overall runtime:** `MULTI_AI_RUNTIME_PASS` (Opus 5.5 **`CLAUDE_OPUS55_DISPATCH_PASS`** via R2; gateway **$10 reserve unverified**)  
**M7 development:** `M7_IN_PROGRESS`  
**Historical note:** `DEVELOPER_ROUTING_BYPASSED` for P1 pipeline commits `4d2c810` / `da9d6dc` (Controller direct implementation).

## Dispatch matrix

| Role | CLI | Requested model | Observed model | Task ID | Mode | Exit | Artifact verification | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Controller | cursor-agent | composer-2.5 | **UNVERIFIED** | m7_v1_4_orchestration | orchestration | 0 | session active | **UNVERIFIED** |
| Developer P1 | opencode | zhipuai-coding-plan/glm-5.3-flash | zhipuai-coding-plan/glm-5.3-flash (JSON + read tool) | M7_DISPATCH_GLM_FLASH | read-only | 0 | status=PENDING truth_level=ASSEMBLY_CONFIRMED ✓ | **DISPATCH_PASS** |
| Developer P1+ | opencode | zhipuai-coding-plan/glm-5.3 | zhipuai-coding-plan/glm-5.3 | M7_DISPATCH_GLM_53 | read-only | 0 | constant NOT_EVALUABLE ✓ | **DISPATCH_PASS** |
| Developer P2 | cursor-agent -p | grok-4.7-high-fast | grok-4.7-high-fast (--trust) | M7_DISPATCH_GROK_HIGH | plan | 0 | partition/cpus/mem ✓ | **DISPATCH_PASS** |
| Developer P2+ | cursor-agent -p | grok-4.7-xhigh-fast | grok-4.7-xhigh-fast | M7_DISPATCH_GROK_XHIGH | plan | 0 | job-name=m7p1_stage ✓ | **DISPATCH_PASS** |
| Independent Reviewer | codex | gpt-6.1-sol | gpt-6.1-sol | M7_DISPATCH_SOL | read-only | 0 | gzip_sha256 verification YES ✓ | **DISPATCH_PASS** |
| Independent Reviewer | codex | gpt-6.1-terra | gpt-6.1-terra | M7_DISPATCH_TERRA | read-only | 0 | SLURM_JOB_ID guard YES ✓ | **DISPATCH_PASS** |
| Technical Architect | codex | gpt-6-astra | gpt-6-astra | M7_DISPATCH_ASTRA | read-only | 0 | post_edit_variants.tsv ✓ | **DISPATCH_PASS** |
| Validator | codex | gpt-6-luna | gpt-6-luna | M7_DISPATCH_LUNA | read-only | 0 | EXIT_CODE=0 from evidence file ✓ | **DISPATCH_PASS** |
| Scientific Architect | claude `-p` | claude-opus-5-5 | **claude-opus-5-5** (`modelUsage`) | M7_CLAUDE_OPUS55_SMOKE_R2 | plan | 0 | NO MOB/JC truth ✓ M7-P0b | **DISPATCH_PASS** |
| Scientific Architect | claude `-p` | claude-opus-5-5 | claude-opus-5-5 | M7_CLAUDE_OPUS55_SMOKE | plan | 0 | budget cap $0.15 | **DISPATCH_PARTIAL** |
| Scientific Architect (historical preflight) | claude | claude-opus-5-5 | — | M7 (prior) | — | — | loggedIn false only | **AUTH_BLOCKED** (not erased) |
| Scientific Architect (Opus 5 smoke) | claude | **claude-opus-5** | **claude-opus-5** | M7_CLAUDE_OPUS5_DISPATCH_SMOKE | plan | 0 | max_turns | **DISPATCH_PARTIAL** |

## Developer handoff acceptance

| Test | Task ID | Executor | Controller wrote impl? | Git/artifact | Tests | Review | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| OpenCode impl | M7_OPENCODE_DEV_HANDOFF | OpenCode GLM Flash | **No** | `M7_P1_FASTQ_STAGING.md` (79 lines) | (doc only) | Sol pending/ACCEPT | **ROLE_ACCEPTANCE_PASS** |
| Grok write | M7_GROK_HANDOFF_WRITE | Grok high-fast | **No** | `.prokadiff/grok_handoff_probe.md` | N/A (probe) | N/A | **ROLE_ACCEPTANCE_PASS** |

## Luna validation (separate from dispatch)

| Task ID | Evidence | Verdict |
| --- | --- | --- |
| M7_LUNA_VALIDATION | pytest + sacct blocked in sandbox | **FAIL** (infra) |
| M7_LUNA_VALIDATION_R2 | `.prokadiff/m7_pytest_evidence.txt`, `.prokadiff/sacct_3049255.txt`, `run_provenance.tsv` | **LUNA_VALIDATION_PASS** |

## Slurm evidence (M7-P1)

| JobID | Run | State | App |
| --- | --- | --- | --- |
| 3049255 | SRR21706776 | COMPLETED 0:0 | provenance COMPLETE, 2867192 pairs |
| 3049302 | SRR21706795 | COMPLETED 0:0 | provenance COMPLETE |
| 3049330 | P1b differential | RUNNING | negative concordance smoke |

## Commands (reference)

```bash
# OpenCode dispatch
opencode run --standalone --format json -m zhipuai-coding-plan/glm-5.3-flash -- "TASK_ID: ..."

# Grok dispatch (requires --trust for non-interactive worktree)
agent -p --trust --model grok-4.7-high-fast --mode plan -- "TASK_ID: ..."

# Codex
codex exec -m gpt-6-luna -s read-only --ephemeral -C <worktree> "..."
```

## Overall PASS checklist

1. OpenCode GLM routing — **verified**
2. Both GLM models dispatched — **verified**
3. Grok separate process — **verified** (`--trust` required non-interactive)
4. Grok isolation / no recursive agents — **verified** (handoff probe)
5. Codex review — **verified**
6. Luna dispatch + current validation — **verified** (R2)
7. No unauthorized recursive AI — **verified**
8. Single writer worktree — **verified** (M7 worktree)
9. Execution evidence on disk — **verified**
10. Scientific/budget gates — Opus 5.5 **R2 PASS**; Opus 5 / Opus 5.5 R1 **PARTIAL**; ledger cumulative **$0.43928375**
