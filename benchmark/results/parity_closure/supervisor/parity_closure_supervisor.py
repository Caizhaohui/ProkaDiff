#!/usr/bin/env python3
"""Parity Closure Supervisor — independent of Cursor chat session lifetime.

Orchestration only on login node: locks, state, agent dispatch, Slurm sacct.
Does not run cargo, bowtie2, breseq, or WGS benchmarks itself.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_DEFAULT = Path("/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff")
AGENT_CLI = os.environ.get("PARITY_AGENT_CLI", "agent")
CONTROLLER_MODEL = os.environ.get("PARITY_CONTROLLER_MODEL", "grok-4.7-high-fast")
MAX_CONTROLLER_ROUNDS_PER_LOOP = int(os.environ.get("PARITY_MAX_CONTROLLER_ROUNDS", "1"))
MAX_STALL_ROUNDS = int(os.environ.get("PARITY_MAX_STALL_ROUNDS", "2"))
MAX_RETRIES = int(os.environ.get("PARITY_MAX_RETRIES", "3"))


def iso_now() -> str:
    return datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")


def log(msg: str, log_path: Path) -> None:
    line = f"{iso_now()} {msg}\n"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as f:
        f.write(line)
    print(line, end="")


@dataclass
class Paths:
    repo: Path
    state_dir: Path
    state_json: Path
    control_json: Path
    supervisor_lock: Path
    controller_lock: Path
    writer_lock: Path
    runs_dir: Path
    log: Path
    baseline_w3: Path

    @classmethod
    def from_repo(cls, repo: Path) -> Paths:
        sd = repo / "benchmark/results/parity_closure/state"
        return cls(
            repo=repo,
            state_dir=sd,
            state_json=sd / "parity_closure_state.json",
            control_json=sd / "supervisor_control.json",
            supervisor_lock=sd / "supervisor.lock",
            controller_lock=sd / "controller.lock",
            writer_lock=sd / "writer.lock",
            runs_dir=sd / "controller_runs",
            log=sd / "supervisor.log",
            baseline_w3=sd / "baseline_w3_floors.json",
        )


class FileLock:
    def __init__(self, path: Path, holder: str, log: Path):
        self.path = path
        self.holder = holder
        self.log = log
        self.fp = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.fp = open(self.path, "a+")
        try:
            fcntl.flock(self.fp.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            self.fp.seek(0)
            meta = self.fp.read()
            log(f"LOCK_BUSY {self.path.name} meta={meta.strip()!r}", self.log)
            self.fp.close()
            self.fp = None
            raise SystemExit(75)
        self.fp.seek(0)
        self.fp.truncate()
        self.fp.write(
            f"pid={os.getpid()}\nholder={self.holder}\nstarted={iso_now()}\nhost={os.uname().nodename}\n"
        )
        self.fp.flush()
        return self

    def __exit__(self, *args):
        if self.fp:
            fcntl.flock(self.fp.fileno(), fcntl.LOCK_UN)
            self.fp.close()


def pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def read_lock_pid(lock_path: Path) -> int | None:
    if not lock_path.exists():
        return None
    text = lock_path.read_text(errors="replace")
    m = re.search(r"^pid=(\d+)", text, re.M)
    return int(m.group(1)) if m else None


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def atomic_write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(data, f, indent=2)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def default_control() -> dict[str, Any]:
    return {
        "pause": False,
        "stop": False,
        "pause_reason": "",
        "updated": iso_now(),
    }


def ensure_state_schema(state: dict[str, Any]) -> dict[str, Any]:
    sup = state.setdefault("supervisor", {})
    sup.setdefault("version", 1)
    sup.setdefault("continuation_test_round", 0)
    sup.setdefault("continuation_test_pass", False)
    sup.setdefault("stall_rounds", 0)
    sup.setdefault("controller_round", 0)
    sup.setdefault("last_progress_signature", "")
    sup.setdefault("engineering_stalled", False)
    sup.setdefault("auto_resume", "SUPERVISOR_MANAGED")
    state.setdefault("writer_owner", None)
    state.setdefault("protected_floor", {
        "slurm_job": "3063624",
        "matched_min": 6,
        "mob": "5/5",
        "breseq_ids": ["89", "90", "91", "94", "96", "110"],
    })
    return state


def git_snapshot(repo: Path) -> dict[str, str]:
    def run(cmd: list[str]) -> str:
        return subprocess.check_output(cmd, cwd=repo, text=True, stderr=subprocess.DEVNULL).strip()

    sha = run(["git", "rev-parse", "HEAD"])
    dirty = "yes" if subprocess.call(["git", "diff", "--quiet"], cwd=repo) else "no"
    diff_stat = ""
    if dirty == "yes":
        diff_stat = run(["git", "diff", "--stat"])
    return {"sha": sha, "dirty": dirty, "diff_stat_head": diff_stat[:4000]}


def slurm_terminal(job_id: str) -> dict[str, Any]:
    if not job_id or not job_id.isdigit():
        return {"known": False}
    try:
        out = subprocess.check_output(
            ["sacct", "-j", job_id, "--format=State,ExitCode", "-n", "-P"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return {"known": False}
    if not out:
        return {"known": False}
    line = out.splitlines()[0]
    state, exit_code = line.split("|", 1)
    return {
        "known": True,
        "state": state,
        "exit_code": exit_code,
        "terminal": state in {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL"},
        "running": state == "RUNNING",
    }


def worker_pids_alive(state: dict[str, Any]) -> bool:
    alive: list[int] = []
    for pid in state.get("worker_pids") or []:
        p = int(pid)
        if pid_alive(p):
            alive.append(p)
    if alive != state.get("worker_pids"):
        state["worker_pids"] = alive
    return bool(alive)


def progress_signature(state: dict[str, Any]) -> str:
    jc = state.get("jc_gate") or {}
    return "|".join(
        str(jc.get(k, ""))
        for k in ("last_matched", "last_over", "last_under", "best_floor_job")
    )


def verify_cli(cli_binary: str | None = None) -> dict[str, Any]:
    cli = cli_binary or AGENT_CLI
    info: dict[str, Any] = {"agent": cli}
    try:
        ver_flag = "--version" if "agy" in cli or "codex" in cli or "opencode" in cli else "-v"
        info["version"] = subprocess.check_output(
            [cli, ver_flag], text=True, stderr=subprocess.STDOUT
        ).strip()
    except Exception as e:
        info["version_error"] = str(e)
    try:
        help_text = subprocess.check_output([cli, "--help"], text=True, stderr=subprocess.STDOUT)
        info["supports_print"] = "-p" in help_text or "--print" in help_text
        info["supports_json"] = "json" in help_text
        info["supports_resume"] = "--resume" in help_text or "--continue" in help_text or "-c" in help_text
    except Exception as e:
        info["help_error"] = str(e)
    return info


def build_controller_prompt(
    paths: Paths,
    state: dict[str, Any],
    mode: str,
    extra: str = "",
) -> str:
    base = f"""TASK_ID: PARITY_CLOSURE_CONTROLLER_ROUND
Role: Cursor CLI grok-4.7-high-fast Controller ONLY — orchestrate; do NOT implement large Rust production changes (delegate OpenCode GLM). Never launch agy or Gemini.
Repo: {paths.repo}
Read: benchmark/results/parity_closure/state/parity_closure_state.json
Read: HANDOFF.md ACTIVE_TASK.md AGENTS.md (orchestration rules)
Rules: Multi-AI Collaboration V3.2 (~/.config/ai-collaboration/POLICY.md), Slurm on compute, no whitelist expansion, no clonal_2673729 overwrite, no github push.
Current status: {state.get('status')} next_action: {state.get('next_action')}
Protected floor: w3 job 3063624 matched=6 over=22 under=12 MOB 5/5 BR89/90/91/94/96/110.
{extra}
At end: update parity_closure_state.json checkpoint atomically via Python or documented fields; update HANDOFF.md briefly.
"""
    if mode == "continuation_test":
        r = state.get("supervisor", {}).get("continuation_test_round", 0)
        return (
            base
            + f"""
CONTINUATION_TEST: Update benchmark/results/parity_closure/state/parity_closure_state.json:
1. Increment supervisor.continuation_test_round from {r} to {r + 1}.
2. Set supervisor.last_controller_note to "CONTINUATION_TEST round {r + 1}: checkpoint verified via cursor grok-4.7-high-fast".
Run python3 or edit_file on benchmark/results/parity_closure/state/parity_closure_state.json now to save this state.
Do not dispatch Slurm or Grok in this test round.
"""
        )
    return base + extra


def run_agent_controller(
    paths: Paths,
    prompt: str,
    run_id: str,
    model: str = CONTROLLER_MODEL,
    plan_only: bool = False,
) -> dict[str, Any]:
    run_dir = paths.runs_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "prompt.txt").write_text(prompt)
    if "agy" in AGENT_CLI or "gemini" in model:
        log(f"GEMINI_DISPATCH_REFUSED cli={AGENT_CLI} model={model}", paths.log)
        return {"run_id": run_id, "exit_code": 3, "model": model, "refused": "gemini"}
    cmd = [
        AGENT_CLI,
        "-p",
        "--output-format",
        "json",
        "--mode",
        "plan" if plan_only else "ask",
        "--model",
        model,
        prompt,
    ]
    env = os.environ.copy()
    log(f"AGENT_START run_id={run_id} model={model} plan={plan_only}", paths.log)
    proc = subprocess.run(
        cmd,
        cwd=paths.repo,
        capture_output=True,
        text=True,
        env=env,
        timeout=int(os.environ.get("PARITY_AGENT_TIMEOUT_SEC", "900")),
    )
    (run_dir / "stdout.raw").write_text(proc.stdout)
    (run_dir / "stderr.log").write_text(proc.stderr)
    meta: dict[str, Any] = {
        "run_id": run_id,
        "exit_code": proc.returncode,
        "model": model,
        "started": iso_now(),
        "git": git_snapshot(paths.repo),
    }
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import orchestration_runtime as runtime

    parsed = runtime.parse_agent_stdout(proc.stdout)
    session_id = parsed.get("session_id")
    if parsed.get("parse_error"):
        meta["parse_error"] = parsed["parse_error"]
    else:
        meta["result"] = str(parsed.get("result") or "")[:2000]
        meta["session_id"] = session_id
        meta["usage"] = parsed.get("usage")
        meta["status"] = parsed.get("status")
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    log(
        f"AGENT_END run_id={run_id} exit={proc.returncode} session_id={session_id}",
        paths.log,
    )
    return meta


def update_stall_counter(state: dict[str, Any]) -> dict[str, Any]:
    sig = progress_signature(state)
    sup = state.setdefault("supervisor", {})
    prev = sup.get("last_progress_signature", "")
    if sig and sig == prev:
        sup["stall_rounds"] = int(sup.get("stall_rounds", 0)) + 1
    else:
        sup["stall_rounds"] = 0
        sup["last_progress_signature"] = sig
    if int(sup.get("stall_rounds", 0)) >= MAX_STALL_ROUNDS:
        state["status"] = "ENGINEERING_STALLED"
        state["block_reason"] = (
            f"No jc_gate progress for {MAX_STALL_ROUNDS} supervisor rounds"
        )
        sup["engineering_stalled"] = True
    return state


def decide_next_action(state: dict[str, Any], control: dict[str, Any]) -> str:
    if control.get("stop"):
        return "STOP"
    if control.get("pause"):
        return "PAUSE"
    if state.get("status") in {"PROJECT_SUCCESS", "HARD_BLOCKED", "ACCEPTED"}:
        return "TERMINAL"
    if state.get("supervisor", {}).get("engineering_stalled"):
        return "STALLED"
    if worker_pids_alive(state):
        return "WAIT_WORKER"
    slurm_id = str(state.get("slurm_id") or "")
    sacct = slurm_terminal(slurm_id)
    if sacct.get("running"):
        return "WAIT_SLURM"
    if state.get("status") == "WAITING_FOR_SLURM" and slurm_id and not sacct.get("terminal"):
        return "WAIT_SLURM"
    return "RUN_CONTROLLER"


def continuation_test(paths: Paths) -> int:
    """Run two controller rounds; verify round 2 sees round 1 checkpoint."""
    state = ensure_state_schema(load_json(paths.state_json))
    state["supervisor"]["continuation_test_round"] = 0
    state["supervisor"]["continuation_test_pass"] = False
    state["supervisor"]["last_controller_note"] = "test_reset"
    atomic_write_json(paths.state_json, state)

    r1 = run_agent_controller(
        paths,
        build_controller_prompt(paths, state, "continuation_test"),
        run_id=f"continuation_test_r1_{int(time.time())}",
        plan_only=False,
    )
    if r1["exit_code"] != 0:
        log(f"CONTINUATION_TEST_FAIL r1 exit={r1['exit_code']}", paths.log)
        return 1
    state = ensure_state_schema(load_json(paths.state_json))
    if int(state["supervisor"].get("continuation_test_round", 0)) != 1:
        log(
            f"CONTINUATION_TEST_FAIL expected round=1 got={state['supervisor'].get('continuation_test_round')}",
            paths.log,
        )
        return 1

    time.sleep(2)
    r2 = run_agent_controller(
        paths,
        build_controller_prompt(paths, state, "continuation_test"),
        run_id=f"continuation_test_r2_{int(time.time())}",
        plan_only=False,
    )
    if r2["exit_code"] != 0:
        log(f"CONTINUATION_TEST_FAIL r2 exit={r2['exit_code']}", paths.log)
        return 1
    state = ensure_state_schema(load_json(paths.state_json))
    if int(state["supervisor"].get("continuation_test_round", 0)) != 2:
        log(
            f"CONTINUATION_TEST_FAIL expected round=2 got={state['supervisor'].get('continuation_test_round')}",
            paths.log,
        )
        return 1
    if read_lock_pid(paths.controller_lock) and pid_alive(read_lock_pid(paths.controller_lock) or -1):
        log("CONTINUATION_TEST_FAIL controller lock still held", paths.log)
        return 1
    state["supervisor"]["continuation_test_pass"] = True
    state["supervisor"]["auto_resume"] = "AUTO_RESUME_PASS"
    atomic_write_json(paths.state_json, state)
    log("AUTO_CONTINUATION_TEST_PASS", paths.log)
    return 0


def pause_test(paths: Paths) -> int:
    ctrl = default_control()
    ctrl["pause"] = True
    ctrl["pause_reason"] = "supervisor_pause_test"
    atomic_write_json(paths.control_json, ctrl)
    state = ensure_state_schema(load_json(paths.state_json))
    if decide_next_action(state, ctrl) != "PAUSE":
        log("PAUSE_TEST_FAIL decide != PAUSE", paths.log)
        return 1
    ctrl["pause"] = False
    ctrl["pause_reason"] = ""
    atomic_write_json(paths.control_json, ctrl)
    if decide_next_action(state, ctrl) == "PAUSE":
        log("PAUSE_TEST_FAIL still paused after resume", paths.log)
        return 1
    log("PAUSE_RESUME_TEST_PASS", paths.log)
    return 0


def dispatch_grok_wave7(paths: Paths, state: dict[str, Any]) -> dict[str, Any]:
    script = paths.repo / "benchmark/results/parity_closure/supervisor/dispatch_g1_grok_wave7.sh"
    if not script.exists():
        return {"ok": False, "error": "missing dispatch script"}
    with FileLock(paths.writer_lock, "grok_wave7", paths.log):
        proc = subprocess.Popen(
            ["bash", str(script)],
            cwd=paths.repo,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    state["worker_pids"] = [proc.pid]
    state["status"] = "WAITING_FOR_WORKER"
    state["worker_classification"] = "DEVELOPER_RUNNING"
    state["developer_model"] = os.environ.get("PARITY_GROK_MODEL", "grok-4.7-xhigh-fast")
    state["worker_wave"] = "wave-7-grok"
    state["next_action"] = "MONITOR_GROK_WAVE7_THEN_SBATCH"
    state["checkpoint"] = iso_now()
    atomic_write_json(paths.state_json, state)
    log(f"GROK_DISPATCH pid={proc.pid}", paths.log)
    return {"ok": True, "pid": proc.pid}


def dispatch_grok_wave8(paths: Paths, state: dict[str, Any]) -> dict[str, Any]:
    # Prefer OpenCode GLM-5.3 per AGENTS.md Developer P1 Advanced routing
    opencode_script = paths.repo / "benchmark/results/parity_closure/supervisor/dispatch_g1_opencode_wave8.sh"
    script = opencode_script if opencode_script.exists() else paths.repo / "benchmark/results/parity_closure/supervisor/dispatch_g1_grok_wave8.sh"
    wave_name = "wave-8-opencode" if opencode_script.exists() else "wave-8-grok"
    dev_model = "zhipuai-coding-plan/glm-5.3" if opencode_script.exists() else os.environ.get("PARITY_GROK_MODEL", "grok-4.7-xhigh-fast")
    if not script.exists():
        return {"ok": False, "error": "missing dispatch script"}
    with FileLock(paths.writer_lock, wave_name, paths.log):
        proc = subprocess.Popen(
            ["bash", str(script)],
            cwd=paths.repo,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    state["worker_pids"] = [proc.pid]
    state["status"] = "WAITING_FOR_WORKER"
    state["worker_classification"] = "DEVELOPER_RUNNING"
    state["developer_model"] = dev_model
    state["worker_wave"] = wave_name
    state["next_action"] = "MONITOR_WAVE8_THEN_SBATCH"
    state["checkpoint"] = iso_now()
    atomic_write_json(paths.state_json, state)
    log(f"WAVE8_DISPATCH wave={wave_name} model={dev_model} pid={proc.pid}", paths.log)
    return {"ok": True, "pid": proc.pid}


def supervisor_once(paths: Paths, *, dispatch_controller: bool) -> int:
    if not paths.control_json.exists():
        atomic_write_json(paths.control_json, default_control())
    control = load_json(paths.control_json)
    state = ensure_state_schema(load_json(paths.state_json))
    state["git_worktree"] = git_snapshot(paths.repo)
    state["checkpoint"] = iso_now()
    action = decide_next_action(state, control)
    log(f"DECIDE action={action} status={state.get('status')}", paths.log)

    if action == "STOP":
        state["status"] = "STOPPED_BY_OPERATOR"
        atomic_write_json(paths.state_json, state)
        return 0
    if action == "PAUSE":
        log("PAUSED operator control", paths.log)
        return 0
    if action in {"WAIT_WORKER", "WAIT_SLURM", "TERMINAL", "STALLED"}:
        atomic_write_json(paths.state_json, state)
        return 0

    na = str(state.get("next_action") or "")
    if (
        action == "RUN_CONTROLLER"
        and "WAVE7" in na
        and not worker_pids_alive(state)
        and "DISPATCH" in na
    ):
        dispatch_grok_wave7(paths, state)
        return 0
    if (
        action == "RUN_CONTROLLER"
        and "WAVE8" in na
        and not worker_pids_alive(state)
        and "DISPATCH" in na
    ):
        dispatch_grok_wave8(paths, state)
        return 0

    if not dispatch_controller:
        atomic_write_json(paths.state_json, state)
        return 0

    with FileLock(paths.controller_lock, "controller_agent", paths.log):
        sup = state.setdefault("supervisor", {})
        sup["controller_round"] = int(sup.get("controller_round", 0)) + 1
        round_n = sup["controller_round"]
        prompt = build_controller_prompt(
            paths,
            state,
            "orchestrate",
            extra=state.get("next_action") or "",
        )
        meta = run_agent_controller(
            paths,
            prompt,
            run_id=f"controller_round_{round_n}_{int(time.time())}",
        )
        state = ensure_state_schema(load_json(paths.state_json))
        state.setdefault("supervisor", {})["last_controller_run"] = meta
        state = update_stall_counter(state)
        atomic_write_json(paths.state_json, state)
    return 0 if meta.get("exit_code") == 0 else 1


def write_baseline_w3(paths: Paths) -> None:
    if paths.baseline_w3.exists():
        return
    data = {
        "slurm_job": "3063624",
        "gate": {"matched": 6, "over": 22, "under": 12, "mob": "5/5"},
        "protected_breseq_jc_ids": ["89", "90", "91", "94", "96", "110"],
        "note": "Floor checkpoint; do not replace with w6b matched-only gains.",
        "frozen_gd": "benchmark/results/parity_closure/g1_clonal_jc/slurm/3063624/prokadiff/output.gd",
    }
    atomic_write_json(paths.baseline_w3, data)


def main() -> int:
    ap = argparse.ArgumentParser(description="Parity Closure Supervisor")
    ap.add_argument("--repo", type=Path, default=REPO_DEFAULT)
    ap.add_argument("--once", action="store_true", help="Single supervisor iteration")
    ap.add_argument("--loop", action="store_true", help="Loop with sleep (orchestration only)")
    ap.add_argument("--sleep", type=float, default=120)
    ap.add_argument("--no-controller", action="store_true", help="Decide only, no agent")
    ap.add_argument("--continuation-test", action="store_true")
    ap.add_argument("--pause-test", action="store_true")
    ap.add_argument("--verify-cli", action="store_true")
    ap.add_argument("--install-queue", action="store_true", help="Add the G1 task queue if missing")
    args = ap.parse_args()
    paths = Paths.from_repo(args.repo.resolve())

    if "agy" in AGENT_CLI or "gemini" in CONTROLLER_MODEL:
        log(f"GEMINI_DISPATCH_REFUSED cli={AGENT_CLI} model={CONTROLLER_MODEL}", paths.log)
        return 3

    if args.verify_cli:
        print(json.dumps(verify_cli(), indent=2))
        return 0

    paths.state_dir.mkdir(parents=True, exist_ok=True)
    write_baseline_w3(paths)

    if args.continuation_test or args.pause_test:
        cli_info = verify_cli()
        atomic_write_json(
            paths.state_dir / "agent_cli_verified.json",
            {"verified_at": iso_now(), **cli_info},
        )

    if args.continuation_test:
        with FileLock(paths.supervisor_lock, "supervisor_test", paths.log):
            return continuation_test(paths)

    if args.pause_test:
        with FileLock(paths.supervisor_lock, "supervisor_pause_test", paths.log):
            return pause_test(paths)

    if args.install_queue or args.loop or args.once:
        return run_persistent_queue(paths, args)

    def run_body() -> int:
        with FileLock(paths.supervisor_lock, "supervisor", paths.log):
            return supervisor_once(paths, dispatch_controller=not args.no_controller)

    log("REFUSING_IMPLICIT_ONESHOT pass --loop for production or --once for a single queue cycle", paths.log)
    return 2


def run_persistent_queue(paths: Paths, args: argparse.Namespace) -> int:
    """Production loop. Dispatch completion returns to the queue instead of exiting."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import orchestration_runtime as runtime

    budget = Path(os.environ.get("PARITY_CLAUDE_BUDGET", str(Path.home() / ".local/state/ai-routing/claude-budget.json")))
    if args.install_queue:
        added = runtime.install_production_queue(paths.state_json, paths.repo, budget)
        log(f"QUEUE_INSTALL added={','.join(added) if added else 'none'}", paths.log)
        if not args.loop and not args.once:
            return 0
    if not paths.control_json.exists():
        atomic_write_json(paths.control_json, default_control())
    if not paths.state_json.exists():
        atomic_write_json(paths.state_json, {"task_queue": {"version": 1, "tasks": [], "controller_cycles": 0}})
    stub = os.environ.get("PARITY_SLURM_STUB")
    max_cycles = os.environ.get("PARITY_MAX_CYCLES")
    source = runtime.file_sha256(Path(__file__).resolve())
    log(
        f"QUEUE_MODE loop={args.loop} once={args.once} sleep={args.sleep} source={source[:12]}",
        paths.log,
    )
    try:
        return runtime.run_queue_loop(
            state_path=paths.state_json,
            control_path=paths.control_json,
            log_path=paths.log,
            supervisor_lock=paths.supervisor_lock,
            writer_lock=paths.writer_lock,
            exit_dir=paths.state_dir / "task_exits",
            sleep_s=float(args.sleep),
            once=bool(args.once and not args.loop),
            max_cycles=int(max_cycles) if max_cycles else None,
            source_hash=source,
            claude_budget=budget,
            slurm_stub=Path(stub) if stub else None,
            quota_state=Path(
                os.environ.get(
                    "PARITY_QUOTA_STATE",
                    str(Path.home() / ".local/state/ai-routing/quota-state.json"),
                )
            ),
        )
    except runtime.LockBusy as exc:
        log(f"SUPERVISOR_LOCK_BUSY {exc}", paths.log)
        return 75


if __name__ == "__main__":
    sys.exit(main())
