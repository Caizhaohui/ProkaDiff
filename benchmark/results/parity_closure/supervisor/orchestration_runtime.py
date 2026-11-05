#!/usr/bin/env python3
"""Deterministic task queue for the parity-closure supervisor.

The supervisor process is the only dispatch owner. This module does not
decide science. It keeps a dependency queue alive across worker and Slurm
completions, including when some tasks are externally blocked.
"""

from __future__ import annotations

import fcntl
import json
import os
import re
import signal
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

READY = "READY"
RUNNING = "RUNNING"
WAITING_WORKER = "WAITING_WORKER"
WAITING_SLURM = "WAITING_SLURM"
BLOCKED_EXTERNAL = "BLOCKED_EXTERNAL"
REVIEW_REQUIRED = "REVIEW_REQUIRED"
COMPLETE = "COMPLETE"
FAILED_RETRYABLE = "FAILED_RETRYABLE"
FAILED_TERMINAL = "FAILED_TERMINAL"

TERMINAL_TASK = {COMPLETE, FAILED_TERMINAL}
ACTIVE_TASK = {RUNNING, WAITING_WORKER, WAITING_SLURM}
ELIGIBLE = {READY, FAILED_RETRYABLE}

CONTROLLER_CLI_DEFAULT = "agent"
CONTROLLER_MODEL_DEFAULT = "grok-4.7-high-fast"


class LockBusy(Exception):
    def __init__(self, path: Path, meta: str):
        super().__init__(f"lock busy: {path} {meta}")
        self.path = path
        self.meta = meta


class GeminiDispatchRefused(Exception):
    pass


def iso_now() -> str:
    return datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")


def atomic_write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def pid_alive(pid: int | None) -> bool:
    if not pid or int(pid) <= 0:
        return False
    try:
        os.kill(int(pid), 0)
    except OSError:
        return False
    # Reap a zombie only when it is our child. A live process we did not
    # spawn raises ChildProcessError here and must still count as alive.
    try:
        waited, _status = os.waitpid(int(pid), os.WNOHANG)
        if waited == int(pid):
            return False
    except ChildProcessError:
        pass
    except OSError:
        pass
    stat_path = Path(f"/proc/{int(pid)}/stat")
    if stat_path.exists():
        try:
            fields = stat_path.read_text().split()
            if len(fields) > 2 and fields[2] == "Z":
                return False
        except OSError:
            return False
    return True


def assert_controller_backend(cli: str, model: str) -> None:
    banned = ("agy", "gemini", "antigravity")
    blob = f"{cli} {model}".lower()
    if any(token in blob for token in banned):
        raise GeminiDispatchRefused(cli, model)


def cursor_launch_argv(prompt: str, cli: str | None = None, model: str | None = None) -> list[str]:
    """Non-interactive Cursor plan invocation. Never selects Gemini."""
    cli = cli or os.environ.get("PARITY_AGENT_CLI", CONTROLLER_CLI_DEFAULT)
    model = model or os.environ.get("PARITY_CONTROLLER_MODEL", CONTROLLER_MODEL_DEFAULT)
    assert_controller_backend(cli, model)
    return [
        cli,
        "-p",
        "--output-format",
        "json",
        "--mode",
        "plan",
        "--model",
        model,
        prompt,
    ]


def parse_agent_stdout(text: str) -> dict[str, Any]:
    """Accept a single JSON document or line-delimited Cursor/agy events."""
    raw = text.strip()
    if not raw:
        return {"parse_error": "empty"}
    try:
        payload = json.loads(raw)
        if isinstance(payload, dict):
            return _normalize_payload(payload)
    except json.JSONDecodeError:
        pass
    chosen: dict[str, Any] | None = None
    for line in raw.splitlines():
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(payload, dict):
            continue
        chosen = _normalize_payload(payload)
        if payload.get("type") == "result" or "result" in payload or "response" in payload:
            break
    return chosen or {"parse_error": "no_json_result"}


def _normalize_payload(payload: dict[str, Any]) -> dict[str, Any]:
    result = payload.get("result", payload.get("response", ""))
    return {
        "result": result if isinstance(result, str) else json.dumps(result)[:2000],
        "session_id": payload.get("session_id") or payload.get("conversation_id"),
        "usage": payload.get("usage"),
        "status": payload.get("status") or payload.get("type"),
    }


@dataclass
class FileLease:
    path: Path
    handle: Any

    def release(self) -> None:
        if self.handle is None:
            return
        try:
            fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
        finally:
            self.handle.close()
            self.handle = None


def try_lease(path: Path) -> FileLease:
    """Kernel flock is the ownership. Stale pid text does not count."""
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = open(path, "a+")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        handle.seek(0)
        meta = handle.read()
        handle.close()
        raise LockBusy(path, meta.strip())
    handle.seek(0)
    handle.truncate()
    handle.write(
        f"pid={os.getpid()}\nholder=supervisor-queue\nstarted={iso_now()}\nhost={os.uname().nodename}\n"
    )
    handle.flush()
    return FileLease(path, handle)


def task_template(
    task_id: str,
    task_type: str,
    worker: str,
    *,
    dependencies: list[str] | None = None,
    status: str = READY,
    priority: int = 100,
    maximum_retries: int = 2,
    command: list[str] | None = None,
    output_paths: list[str] | None = None,
    allowed_while_other_tasks_blocked: bool = True,
    blocked_reason: str = "",
    resume_condition: str = "",
    rate_limit_logs: list[str] | None = None,
    fallback_worker: str = "",
    fallback_command: list[str] | None = None,
    fallback_rate_limit_logs: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "task_id": task_id,
        "task_type": task_type,
        "worker": worker,
        "status": status,
        "dependencies": list(dependencies or []),
        "blocked_reason": blocked_reason,
        "resume_condition": resume_condition,
        "allowed_while_other_tasks_blocked": allowed_while_other_tasks_blocked,
        "retry_count": 0,
        "maximum_retries": maximum_retries,
        "execution_pid": None,
        "output_paths": list(output_paths or []),
        "rate_limit_logs": list(rate_limit_logs or []),
        "fallback_worker": fallback_worker,
        "fallback_command": list(fallback_command or []),
        "fallback_rate_limit_logs": list(fallback_rate_limit_logs or []),
        "escalation_history": [],
        "resume_after": "",
        "source_hash": "",
        "completion_timestamp": "",
        "priority": priority,
        "command": list(command or []),
        "exit_code": None,
        "slurm_job_id": "",
    }


def production_queue(repo: Path, *, claude_ok: bool) -> list[dict[str, Any]]:
    """Owner-specified G1 queue. Production attribution stays untouched."""
    sup = repo / "benchmark/results/parity_closure/supervisor"
    b0 = repo / "benchmark/results/parity_closure/g1_clonal_jc/wave11b0"
    b0r = repo / "benchmark/results/parity_closure/g1_clonal_jc/wave11b0r"
    b2 = repo / "benchmark/results/parity_closure/g1_clonal_jc/wave11b2"
    review = repo / "benchmark/results/parity_closure/g1_clonal_jc/review/codex_sol_wave11b_review.md"
    tasks = [
        task_template(
            "WAVE11B0_SHADOW_AS",
            "analysis",
            "opencode/glm-5.3",
            priority=10,
            command=["bash", str(sup / "dispatch_g1_opencode_wave11b0.sh")],
            rate_limit_logs=[str(b0 / "opencode_wave11b0_stdout.txt")],
            fallback_worker="cursor/grok-4.7-high-fast",
            fallback_command=["bash", str(sup / "dispatch_g1_grok_wave11b0.sh")],
            fallback_rate_limit_logs=[str(b0 / "grok_wave11b0_stdout.log")],
            output_paths=[
                str(b0 / "WAVE11B0_AS_SHADOW_RESULTS.tsv"),
                str(b0 / "WAVE11B0_READ_PROVENANCE.tsv"),
                str(b0 / "WAVE11B0_TIE_AUDIT.tsv"),
                str(b0 / "WAVE11B0_REGRESSION_REPORT.md"),
            ],
        ),
        task_template(
            "WAVE11B0R_ANALYZER_REPAIR",
            "analysis",
            "opencode/glm-5.3",
            priority=15,
            dependencies=["WAVE11B0_SHADOW_AS"],
            command=["bash", str(sup / "dispatch_g1_opencode_wave11b0r.sh")],
            rate_limit_logs=[str(b0r / "opencode_wave11b0r_stdout.txt")],
            fallback_worker="cursor/grok-4.7-high-fast",
            fallback_command=["bash", str(sup / "dispatch_g1_grok_wave11b0r.sh")],
            fallback_rate_limit_logs=[str(b0r / "grok_wave11b0r_stdout.log")],
            output_paths=[
                str(b0r / "WAVE11B0R_AS_SHADOW_RESULTS.tsv"),
                str(b0r / "WAVE11B0R_READ_PROVENANCE.tsv"),
                str(b0r / "WAVE11B0R_TIE_AUDIT.tsv"),
                str(b0r / "WAVE11B0R_REGRESSION_REPORT.md"),
            ],
        ),
        task_template(
            "WAVE11B2_3072550_REGRESSION",
            "analysis",
            "opencode/glm-5.3",
            priority=20,
            command=["bash", str(sup / "dispatch_g1_opencode_wave11b2.sh")],
            rate_limit_logs=[str(b2 / "opencode_wave11b2_stdout.txt")],
            fallback_worker="cursor/grok-4.7-high-fast",
            fallback_command=["bash", str(sup / "dispatch_g1_grok_wave11b2.sh")],
            fallback_rate_limit_logs=[str(b2 / "grok_wave11b2_stdout.log")],
            output_paths=[str(b2 / "WAVE11B2_REGRESSION_REPORT.md")],
        ),
        task_template(
            "WAVE11B1_SCIENTIFIC_REVIEW",
            "scientific_review",
            "claude",
            priority=30,
            status=READY if claude_ok else BLOCKED_EXTERNAL,
            allowed_while_other_tasks_blocked=False,
            blocked_reason="" if claude_ok else "claude_gateway_balance_unverified",
            resume_condition="claude_gateway_balance_verified_with_10usd_reserve_and_project_cap",
            dependencies=["WAVE11B0_SHADOW_AS", "WAVE11B0R_ANALYZER_REPAIR"],
        ),
        task_template(
            "WAVE11D_CODEX_REVIEW",
            "independent_review",
            "codex/gpt-6.1-sol",
            priority=40,
            dependencies=[
                "WAVE11B0_SHADOW_AS",
                "WAVE11B0R_ANALYZER_REPAIR",
                "WAVE11B2_3072550_REGRESSION",
            ],
            command=["bash", str(sup / "dispatch_codex_wave11b_review.sh")],
            output_paths=[str(review)],
        ),
        task_template(
            "WAVE11E_PRODUCTION_AS",
            "production_integration",
            "held",
            priority=50,
            status=BLOCKED_EXTERNAL,
            allowed_while_other_tasks_blocked=False,
            dependencies=[
                "WAVE11B0_SHADOW_AS",
                "WAVE11B0R_ANALYZER_REPAIR",
                "WAVE11B1_SCIENTIFIC_REVIEW",
                "WAVE11D_CODEX_REVIEW",
            ],
            blocked_reason="awaiting_scientific_approval_and_independent_review",
            resume_condition="WAVE11B1_COMPLETE_and_WAVE11D_COMPLETE",
        ),
        task_template(
            "WAVE11F_SLURM_VALIDATION",
            "slurm_validation",
            "slurm",
            priority=60,
            dependencies=["WAVE11E_PRODUCTION_AS"],
            allowed_while_other_tasks_blocked=False,
            resume_condition="production_integration_complete_and_no_duplicate_g1_job",
        ),
    ]
    if not claude_ok:
        for task in tasks:
            if task["task_id"] == "WAVE11B1_SCIENTIFIC_REVIEW":
                task["status"] = BLOCKED_EXTERNAL
    return tasks


def claude_balance_ok(budget_path: Path, *, now: datetime | None = None, max_age_hours: float = 24.0) -> tuple[bool, str]:
    """Live gateway observation only. The local spend ledger is not a balance."""
    if not budget_path.is_file():
        return False, "claude_budget_file_missing"
    try:
        data = json.loads(budget_path.read_text())
    except json.JSONDecodeError:
        return False, "claude_budget_unreadable"
    stamp = str(data.get("last_verified_balance_at_utc") or "")
    balance = data.get("last_verified_balance_usd")
    reserve = float(data.get("desired_reserve_usd") or 10.0)
    if balance is None or not stamp:
        return False, "claude_balance_unverified"
    try:
        verified = datetime.fromisoformat(stamp.replace("Z", "+00:00"))
    except ValueError:
        return False, "claude_balance_timestamp_invalid"
    current = now or datetime.now(timezone.utc)
    if verified.tzinfo is None:
        verified = verified.replace(tzinfo=timezone.utc)
    age_h = (current - verified).total_seconds() / 3600.0
    if age_h > max_age_hours:
        return False, f"claude_balance_stale_age_hours={age_h:.1f}"
    if float(balance) - reserve < 0:
        return False, "claude_reserve_would_be_consumed"
    return True, "claude_balance_verified"


def index_tasks(tasks: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {task["task_id"]: task for task in tasks}


def dependencies_met(task: dict[str, Any], by_id: dict[str, dict[str, Any]]) -> bool:
    for dep in task.get("dependencies") or []:
        other = by_id.get(dep)
        if other is None or other.get("status") != COMPLETE:
            return False
    return True


def eligible_tasks(tasks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_id = index_tasks(tasks)
    ready = []
    for task in tasks:
        if task.get("status") not in ELIGIBLE:
            continue
        if int(task.get("retry_count") or 0) > int(task.get("maximum_retries") or 0):
            continue
        if task.get("status") == FAILED_RETRYABLE and int(task.get("retry_count") or 0) >= int(
            task.get("maximum_retries") or 0
        ):
            continue
        if not dependencies_met(task, by_id):
            continue
        if not task.get("command"):
            continue
        ready.append(task)
    ready.sort(key=lambda item: (int(item.get("priority", 100)), item["task_id"]))
    return ready


def blocked_summary(tasks: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for task in tasks:
        if task.get("status") in {BLOCKED_EXTERNAL, FAILED_TERMINAL, REVIEW_REQUIRED}:
            rows.append(
                {
                    "task_id": task["task_id"],
                    "status": task["status"],
                    "blocked_reason": task.get("blocked_reason") or "",
                    "resume_condition": task.get("resume_condition") or "",
                }
            )
    return rows


def select_decision(ready: list[dict[str, Any]]) -> str | None:
    if not ready:
        return None
    return ready[0]["task_id"]


def classify_exit(task: dict[str, Any], exit_code: int | None) -> None:
    task["exit_code"] = exit_code
    task["execution_pid"] = None
    if exit_code == 0:
        task["status"] = COMPLETE
        task["completion_timestamp"] = iso_now()
        task["blocked_reason"] = ""
        return
    task["retry_count"] = int(task.get("retry_count") or 0) + 1
    if task["retry_count"] >= int(task.get("maximum_retries") or 0):
        task["status"] = FAILED_TERMINAL
        task["blocked_reason"] = f"exit={exit_code}"
        task["completion_timestamp"] = iso_now()
    else:
        task["status"] = FAILED_RETRYABLE
        task["blocked_reason"] = f"exit={exit_code}"


VERDICT_RE = re.compile(
    r"(?:WAVE\d+[A-Z]*_VERDICT\s*:\s*|verdict[^A-Za-z]{0,8})(PASS|APPROVE\w*|ACCEPT\w*|NEEDS_REVISION|FAIL\w*|REJECT\w*|BLOCKED)",
    re.IGNORECASE,
)

REJECTED_VERDICTS = {
    "NEEDS_REVISION",
    "FAIL",
    "FAILED",
    "REJECT",
    "REJECTED",
    "BLOCKED",
    "CHANGES_REQUESTED",
}


def review_verdict(task: dict[str, Any]) -> str:
    """Read the reviewer's own verdict from the artifact it was asked to write.

    Exit code 0 only proves the reviewer process ran. A review that says
    NEEDS_REVISION or FAIL must not be recorded as COMPLETE, otherwise the
    downstream production-integration dependency would open on a rejected
    review.
    """
    for item in task.get("output_paths") or []:
        path = Path(item)
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        match = VERDICT_RE.search(text[:8000])
        if match:
            return match.group(1).upper()
    return ""


def apply_review_verdict(task: dict[str, Any], log_path: Path) -> bool:
    """Downgrade an exit-0 review task whose verdict rejects the work."""
    if task.get("task_type") != "independent_review":
        return False
    verdict = review_verdict(task)
    if not verdict:
        task["review_verdict"] = "UNSPECIFIED"
        log_line(log_path, f"REVIEW_VERDICT_MISSING {task['task_id']} treated as not approved")
        task["status"] = REVIEW_REQUIRED
        task["blocked_reason"] = "review_verdict_missing"
        return True
    task["review_verdict"] = verdict
    if any(token in verdict for token in REJECTED_VERDICTS) or verdict.startswith(
        ("NEEDS_REVISION", "FAIL", "REJECT", "BLOCKED")
    ):
        task["status"] = REVIEW_REQUIRED
        task["blocked_reason"] = f"review_verdict={verdict}"
        log_line(log_path, f"REVIEW_REJECTED {task['task_id']} verdict={verdict} status={task['status']}")
        return True
    log_line(log_path, f"REVIEW_APPROVED {task['task_id']} verdict={verdict}")
    return False


RATE_LIMIT_MARKERS = (
    "provider.rate-limit",
    '"code":"1308"',
    "使用上限",
    "rate limit exceeded",
    "429 Too Many Requests",
)

RESET_RE = re.compile(r"(\d{4}-\d{2}-\d{2})[ T](\d{2}:\d{2}:\d{2})")


def _log_tail(path: Path, limit: int = 65536) -> str:
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - limit))
            return handle.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def detect_rate_limit(task: dict[str, Any]) -> str:
    """Return an ISO resume time when the worker log shows a provider quota block."""
    blob = ""
    for item in task.get("rate_limit_logs") or []:
        blob += _log_tail(Path(item))
    if not blob:
        return ""
    if not any(marker in blob for marker in RATE_LIMIT_MARKERS):
        return ""
    match = None
    for match in RESET_RE.finditer(blob):
        pass
    if match is None:
        return ""
    try:
        when = datetime.strptime(f"{match.group(1)} {match.group(2)}", "%Y-%m-%d %H:%M:%S")
    except ValueError:
        return ""
    return when.astimezone().isoformat(timespec="seconds")


def classify_rate_limited(task: dict[str, Any], exit_code: int | None, resume_at: str) -> None:
    """Quota/rate limiting is an external block, not a task failure."""
    task["exit_code"] = exit_code
    task["execution_pid"] = None
    task["status"] = BLOCKED_EXTERNAL
    task["blocked_reason"] = f"worker_rate_limit exit={exit_code}"
    task["resume_condition"] = f"provider_quota_reset_at={resume_at}"
    task["resume_after"] = resume_at
    task["completion_timestamp"] = ""


def escalate_worker(task: dict[str, Any], log_path: Path, resume_at: str, *, reason: str = "provider_rate_limit") -> bool:
    """GLM quota exhausted -> continue the same task on the fallback developer.

    Mirrors the fastBact policy: OpenCode GLM is primary, Cursor
    grok-4.7-high-fast takes over when the GLM pool is rate limited instead of
    parking the task until the subscription window resets.
    """
    if not task.get("fallback_command"):
        return False
    history = task.setdefault("escalation_history", [])
    if any(item.get("to_worker") == task.get("fallback_worker") for item in history):
        return False
    history.append(
        {
            "from_worker": task.get("worker"),
            "to_worker": task.get("fallback_worker"),
            "reason": reason,
            "quota_reset_at": resume_at,
            "at": iso_now(),
        }
    )
    task["worker"] = task.get("fallback_worker") or task.get("worker")
    task["command"] = list(task["fallback_command"])
    task["rate_limit_logs"] = list(task.get("fallback_rate_limit_logs") or [])
    task["status"] = READY
    task["blocked_reason"] = ""
    task["resume_condition"] = "fallback_worker_active"
    task["resume_after"] = ""
    task["execution_pid"] = None
    task["exit_code"] = None
    log_line(
        log_path,
        f"ESCALATE_WORKER {task['task_id']} from={history[-1]['from_worker']} "
        f"to={task['worker']} reason={reason} quota_reset_at={resume_at}",
    )
    return True


def quota_pool_status(state_path: Path) -> dict[str, Any]:
    """Read the shared cross-project quota state (ai-routing).

    fastBact checks this before choosing a developer. An unreadable file means
    availability is unconfirmed, never "available".
    """
    if not state_path.is_file():
        return {"known": False, "status": "UNKNOWN", "next_probe": ""}
    try:
        data = json.loads(state_path.read_text())
        pool = (data.get("pools") or {}).get("bigmodel-coding-plan") or {}
    except (json.JSONDecodeError, AttributeError):
        return {"known": False, "status": "UNKNOWN", "next_probe": ""}
    return {
        "known": True,
        "status": str(pool.get("status") or "UNKNOWN"),
        "next_probe": str(pool.get("next_probe_at_utc") or ""),
    }


def glm_pool_blocked(state_path: Path, now: datetime | None = None) -> tuple[bool, str]:
    """True when the GLM coding plan pool is rate limited and still in cooldown."""
    pool = quota_pool_status(state_path)
    if pool["status"] != "RATE_LIMITED":
        return False, pool["status"]
    stamp = pool["next_probe"]
    if not stamp:
        return True, "RATE_LIMITED"
    try:
        until = datetime.fromisoformat(stamp.replace("Z", "+00:00"))
    except ValueError:
        return True, "RATE_LIMITED"
    current = now or datetime.now(timezone.utc)
    if until.tzinfo is None:
        until = until.replace(tzinfo=timezone.utc)
    if current >= until:
        return False, f"RATE_LIMITED_COOLDOWN_ELAPSED_UNTIL_{stamp}"
    return True, f"RATE_LIMITED_UNTIL_{stamp}"


def pre_escalate_if_pool_blocked(
    task: dict[str, Any], quota_state: Path, log_path: Path, now: datetime | None = None
) -> bool:
    """Skip a known-limited GLM pool instead of burning a dispatch on it."""
    if not str(task.get("worker", "")).startswith("opencode"):
        return False
    blocked, why = glm_pool_blocked(quota_state, now)
    if not blocked:
        return False
    return escalate_worker(task, log_path, why, reason="provider_quota_state")


def _completed_after(later: str, earlier: str) -> bool:
    if not later or not earlier:
        return False
    return later > earlier


def requeue_review_tasks(tasks: list[dict[str, Any]], log_path: Path) -> list[str]:
    """Requeue a rejected review only after an upstream deliverable is newer than that review.

    Dependencies that merely exist and already completed do not requeue the review.
    That guard stops the loop in which a finished NEEDS_REVISION review is dispatched
    again against the same artifacts.
    """
    by_id = index_tasks(tasks)
    events: list[str] = []
    for task in tasks:
        if task.get("status") != REVIEW_REQUIRED:
            continue
        if not dependencies_met(task, by_id):
            continue
        if not task.get("command"):
            continue
        reviewed_at = str(task.get("completion_timestamp") or "")
        if reviewed_at:
            newer = [
                dep
                for dep in task.get("dependencies") or []
                if _completed_after(str(by_id[dep].get("completion_timestamp") or ""), reviewed_at)
            ]
            if not newer:
                continue
        task["status"] = READY
        task["blocked_reason"] = ""
        event = f"REQUEUE_REVIEW {task['task_id']} dependencies_satisfied"
        log_line(log_path, event)
        events.append(event)
    return events


def refresh_rate_limits(tasks: list[dict[str, Any]], log_path: Path, now: datetime | None = None) -> list[str]:
    """Reactivate tasks whose provider quota window has ended."""
    current = now or datetime.now(timezone.utc).astimezone()
    events = []
    for task in tasks:
        if task.get("status") != BLOCKED_EXTERNAL:
            continue
        if task.get("resume_condition") == "fallback_worker_active":
            continue
        stamp = task.get("resume_after") or ""
        if not stamp or not str(task.get("blocked_reason", "")).startswith("worker_rate_limit"):
            continue
        try:
            when = datetime.fromisoformat(stamp)
        except ValueError:
            continue
        if when.tzinfo is None:
            when = when.astimezone()
        if current < when:
            continue
        task["status"] = READY
        task["blocked_reason"] = ""
        task["resume_condition"] = ""
        task["resume_after"] = ""
        event = f"REACTIVATE_RATE_LIMITED {task['task_id']} quota_reset_passed={stamp}"
        log_line(log_path, event)
        events.append(event)
    return events


def apply_slurm(task: dict[str, Any], probe: dict[str, Any]) -> str | None:
    if task.get("status") != WAITING_SLURM:
        return None
    if not probe.get("terminal"):
        return "WAIT_SLURM"
    code = str(probe.get("exit_code") or "")
    ok = probe.get("state") == "COMPLETED" and code.startswith("0:")
    classify_exit(task, 0 if ok else 1)
    return "COLLECT_SLURM"


def read_exit_file(path: Path) -> int | None:
    if not path.is_file():
        return None
    text = path.read_text().strip()
    if not text:
        return None
    try:
        return int(text.splitlines()[-1].strip())
    except ValueError:
        return None


def file_sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def log_line(log_path: Path, message: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    line = f"{iso_now()} {message}\n"
    with log_path.open("a") as handle:
        handle.write(line)
    print(line, end="")


def slurm_probe(job_id: str, stub_path: Path | None = None) -> dict[str, Any]:
    if stub_path and stub_path.is_file():
        data = json.loads(stub_path.read_text())
        row = data.get(str(job_id)) or {}
        return {
            "known": bool(row),
            "state": row.get("state", ""),
            "exit_code": row.get("exit_code", ""),
            "terminal": row.get("state") in {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL"},
            "running": row.get("state") == "RUNNING",
        }
    if not job_id or not str(job_id).isdigit():
        return {"known": False, "terminal": False, "running": False}
    try:
        out = subprocess.check_output(
            ["sacct", "-j", str(job_id), "--format=State,ExitCode", "-n", "-P"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return {"known": False, "terminal": False, "running": False}
    if not out:
        return {"known": False, "terminal": False, "running": False}
    state, exit_code = out.splitlines()[0].split("|", 1)
    return {
        "known": True,
        "state": state,
        "exit_code": exit_code,
        "terminal": state in {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL"},
        "running": state == "RUNNING",
    }


def ensure_queue_state(state: dict[str, Any], tasks: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    queue = state.setdefault("task_queue", {})
    queue.setdefault("version", 1)
    queue.setdefault("tasks", tasks or [])
    queue.setdefault("controller_cycles", 0)
    queue.setdefault("last_decision", "")
    return state


def reconcile(
    tasks: list[dict[str, Any]],
    *,
    exit_dir: Path,
    children: dict[str, subprocess.Popen[Any]],
    slurm_stub: Path | None,
    log_path: Path,
) -> list[str]:
    events: list[str] = []
    for task in tasks:
        tid = task["task_id"]
        proc = children.get(tid)
        if proc is not None:
            code = proc.poll()
            if code is None:
                task["status"] = WAITING_WORKER
                task["execution_pid"] = proc.pid
                continue
            children.pop(tid, None)
            resume_at = "" if code == 0 else detect_rate_limit(task)
            if resume_at:
                classify_rate_limited(task, code, resume_at)
                events.append(
                    f"WORKER_RATE_LIMITED {tid} exit={code} status={task['status']} resume_at={resume_at}"
                )
                log_line(log_path, events[-1])
                if escalate_worker(task, log_path, resume_at):
                    events.append(f"WORKER_ESCALATED {tid} worker={task['worker']} status={task['status']}")
                    log_line(log_path, events[-1])
            else:
                classify_exit(task, code)
                events.append(f"WORKER_EXIT {tid} exit={code} status={task['status']}")
                log_line(log_path, events[-1])
                if code == 0 and apply_review_verdict(task, log_path):
                    events.append(f"WORKER_VERDICT {tid} status={task['status']} verdict={task.get('review_verdict','-')}")
                    log_line(log_path, events[-1])
            continue
        if task.get("status") in {WAITING_WORKER, RUNNING} and task.get("execution_pid"):
            if pid_alive(task.get("execution_pid")):
                events.append(f"ADOPT_RUNNING_WORKER {tid} pid={task['execution_pid']}")
                log_line(log_path, events[-1])
                task["status"] = WAITING_WORKER
                continue
            code = read_exit_file(exit_dir / f"{tid}.exit")
            code = 1 if code is None else code
            resume_at = detect_rate_limit(task) if code != 0 else ""
            if resume_at:
                classify_rate_limited(task, code, resume_at)
                events.append(
                    f"REAP_WORKER_RATE_LIMITED {tid} exit={code} status={task['status']} resume_at={resume_at}"
                )
                log_line(log_path, events[-1])
                if escalate_worker(task, log_path, resume_at):
                    events.append(f"REAP_WORKER_ESCALATED {tid} worker={task['worker']} status={task['status']}")
                    log_line(log_path, events[-1])
            else:
                classify_exit(task, code)
                events.append(f"REAP_WORKER {tid} exit={task['exit_code']} status={task['status']}")
                log_line(log_path, events[-1])
                if code == 0 and apply_review_verdict(task, log_path):
                    events.append(f"REAP_WORKER_VERDICT {tid} status={task['status']} verdict={task.get('review_verdict','-')}")
                    log_line(log_path, events[-1])
        if task.get("status") == WAITING_SLURM and task.get("slurm_job_id"):
            probe = slurm_probe(str(task["slurm_job_id"]), slurm_stub)
            action = apply_slurm(task, probe)
            if action == "COLLECT_SLURM":
                events.append(f"COLLECT_SLURM {tid} state={probe.get('state')} exit={probe.get('exit_code')}")
                log_line(log_path, events[-1])
            elif action == "WAIT_SLURM":
                events.append(f"WAIT_SLURM {tid}")
    return events


def worktree_busy(tasks: list[dict[str, Any]], children: dict[str, subprocess.Popen[Any]]) -> bool:
    if children:
        return True
    for task in tasks:
        if task.get("status") in {RUNNING, WAITING_WORKER} and pid_alive(task.get("execution_pid")):
            return True
    return False


def spawn_task(
    task: dict[str, Any],
    *,
    exit_dir: Path,
    log_path: Path,
    source_hash: str,
) -> subprocess.Popen[Any]:
    exit_dir.mkdir(parents=True, exist_ok=True)
    exit_file = exit_dir / f"{task['task_id']}.exit"
    if exit_file.exists():
        exit_file.unlink()
    command = task.get("command") or []
    if not command:
        raise RuntimeError(f"{task['task_id']} has no command")
    quoted = " ".join(subprocess.list2cmdline([part]) for part in command)
    wrapper = f"{quoted}; ec=$?; printf '%s\\n' \"$ec\" > {subprocess.list2cmdline([str(exit_file)])}; exit $ec"
    proc = subprocess.Popen(
        ["bash", "-c", wrapper],
        start_new_session=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    task["status"] = WAITING_WORKER
    task["execution_pid"] = proc.pid
    task["source_hash"] = source_hash
    task["completion_timestamp"] = ""
    log_line(log_path, f"SPAWN {task['task_id']} pid={proc.pid} worker={task.get('worker')}")
    return proc


def refresh_external(tasks: list[dict[str, Any]], claude_ok: bool, reason: str) -> None:
    by_id = index_tasks(tasks)
    review = by_id.get("WAVE11B1_SCIENTIFIC_REVIEW")
    if review is None:
        return
    if claude_ok and review.get("status") == BLOCKED_EXTERNAL and dependencies_met(review, by_id):
        review["status"] = READY
        review["blocked_reason"] = ""
    elif not claude_ok and review.get("status") in {READY, FAILED_RETRYABLE}:
        review["status"] = BLOCKED_EXTERNAL
        review["blocked_reason"] = reason


def run_queue_loop(
    *,
    state_path: Path,
    control_path: Path,
    log_path: Path,
    supervisor_lock: Path,
    writer_lock: Path,
    exit_dir: Path,
    sleep_s: float,
    once: bool,
    max_cycles: int | None,
    source_hash: str,
    claude_budget: Path | None,
    slurm_stub: Path | None,
    quota_state: Path | None = None,
    stop_when_idle_cycles: int | None = None,
) -> int:
    """Persistent dispatch loop. A finished worker returns here instead of exiting."""
    lease = try_lease(supervisor_lock)
    children: dict[str, subprocess.Popen[Any]] = {}
    writer: FileLease | None = None
    idle_cycles = 0
    cycles = 0
    stopping = {"flag": False}

    def _stop(_signum: int, _frame: Any) -> None:
        stopping["flag"] = True
        log_line(log_path, f"LOOP_SIGNAL signum={_signum}")

    signal.signal(signal.SIGTERM, _stop)
    signal.signal(signal.SIGINT, _stop)
    try:
        assert_controller_backend(
            os.environ.get("PARITY_AGENT_CLI", CONTROLLER_CLI_DEFAULT),
            os.environ.get("PARITY_CONTROLLER_MODEL", CONTROLLER_MODEL_DEFAULT),
        )
    except GeminiDispatchRefused as exc:
        log_line(log_path, f"GEMINI_DISPATCH_REFUSED {exc}")
        lease.release()
        return 3
    log_line(log_path, f"QUEUE_LOOP_START pid={os.getpid()} once={once}")
    try:
        while not stopping["flag"]:
            cycles += 1
            control = load_json(control_path) if control_path.exists() else {"pause": False, "stop": False}
            if control.get("stop"):
                log_line(log_path, "LOOP_STOP control.stop")
                break
            state = ensure_queue_state(load_json(state_path))
            tasks = state["task_queue"]["tasks"]
            if control.get("pause"):
                log_line(log_path, "PAUSED operator control; loop remains alive")
                atomic_write_json(state_path, state)
                if once:
                    break
                time.sleep(sleep_s)
                continue
            ok, why = (False, "claude_budget_not_configured")
            if claude_budget is not None:
                ok, why = claude_balance_ok(claude_budget)
            refresh_external(tasks, ok, why)
            refresh_rate_limits(tasks, log_path)
            requeue_review_tasks(tasks, log_path)
            reconcile(tasks, exit_dir=exit_dir, children=children, slurm_stub=slurm_stub, log_path=log_path)
            if quota_state is not None:
                for candidate in tasks:
                    if candidate.get("status") == READY:
                        pre_escalate_if_pool_blocked(candidate, quota_state, log_path)
            ready = eligible_tasks(tasks)
            busy = worktree_busy(tasks, children)
            if ready and not busy:
                if writer is None:
                    try:
                        writer = try_lease(writer_lock)
                    except LockBusy as exc:
                        log_line(log_path, f"WRITER_LOCK_BUSY {exc.meta!r}")
                        ready = []
                if ready and writer is not None:
                    chosen = select_decision(ready)
                    state["task_queue"]["controller_cycles"] = int(state["task_queue"].get("controller_cycles") or 0) + 1
                    state["task_queue"]["last_decision"] = chosen
                    log_line(
                        log_path,
                        f"CONTROLLER_CYCLE n={state['task_queue']['controller_cycles']} dispatch={chosen} ready={','.join(t['task_id'] for t in ready)}",
                    )
                    task = index_tasks(tasks)[chosen]
                    if task.get("status") == FAILED_RETRYABLE:
                        pass
                    proc = spawn_task(task, exit_dir=exit_dir, log_path=log_path, source_hash=source_hash)
                    children[chosen] = proc
                    idle_cycles = 0
            else:
                idle_cycles += 1
                held = [task["task_id"] for task in tasks if task.get("status") in ACTIVE_TASK]
                log_line(
                    log_path,
                    "QUEUE_IDLE "
                    + json.dumps(
                        {
                            "held": held,
                            "ready": [task["task_id"] for task in ready],
                            "blocked": blocked_summary(tasks),
                            "claude": why,
                        },
                        ensure_ascii=True,
                    ),
                )
                if not held and writer is not None:
                    writer.release()
                    writer = None
            state["checkpoint"] = iso_now()
            state["status"] = "QUEUE_RUNNING" if any(t.get("status") in ACTIVE_TASK for t in tasks) else "QUEUE_MONITORING"
            state["next_action"] = "QUEUE_LOOP"
            atomic_write_json(state_path, state)
            if once:
                break
            if max_cycles is not None and cycles >= max_cycles:
                log_line(log_path, f"QUEUE_MAX_CYCLES cycles={cycles}")
                break
            if stop_when_idle_cycles is not None and idle_cycles >= stop_when_idle_cycles and not children:
                log_line(log_path, f"QUEUE_IDLE_LIMIT cycles={idle_cycles}")
                break
            time.sleep(sleep_s)
    finally:
        if writer is not None:
            writer.release()
        lease.release()
        log_line(log_path, f"QUEUE_LOOP_END pid={os.getpid()} cycles={cycles}")
    return 0


def install_production_queue(state_path: Path, repo: Path, budget_path: Path) -> list[str]:
    state = load_json(state_path) if state_path.exists() else {}
    ok, why = claude_balance_ok(budget_path)
    incoming = production_queue(repo, claude_ok=ok)
    queue = state.setdefault("task_queue", {"version": 1, "tasks": [], "controller_cycles": 0, "last_decision": ""})
    existing = {task["task_id"] for task in queue.get("tasks") or []}
    added = []
    for task in incoming:
        if task["task_id"] in existing:
            continue
        queue["tasks"].append(task)
        added.append(task["task_id"])
    queue["claude_balance_reason"] = why
    state["controller_platform"] = "cursor-cli"
    state["controller_model"] = CONTROLLER_MODEL_DEFAULT
    state["governance"] = "PROKADIFF_MULTI_AI_V5_QUEUE"
    atomic_write_json(state_path, state)
    return added
