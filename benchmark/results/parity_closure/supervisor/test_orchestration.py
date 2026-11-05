#!/usr/bin/env python3
"""Deterministic orchestration tests. They do not call OpenCode, Codex, or Slurm."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import orchestration_runtime as runtime

SUPERVISOR = HERE / "parity_closure_supervisor.py"
BUDGET = Path("/hpcfs/fhome/caizhh/.local/state/ai-routing/claude-budget.json")


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def make_repo(root: Path) -> tuple[Path, Path]:
    state = root / "benchmark/results/parity_closure/state"
    state.mkdir(parents=True)
    write_json(state / "supervisor_control.json", {"pause": False, "stop": False, "pause_reason": ""})
    write_json(
        state / "parity_closure_state.json",
        {"task_queue": {"version": 1, "tasks": [], "controller_cycles": 0, "last_decision": ""}},
    )
    return root, state


def put_tasks(state: Path, tasks: list[dict]) -> None:
    path = state / "parity_closure_state.json"
    data = json.loads(path.read_text())
    data["task_queue"]["tasks"] = tasks
    write_json(path, data)


def run_sup(root: Path, args: list[str], env: dict | None = None, timeout: float = 30) -> subprocess.CompletedProcess[str]:
    merged = os.environ.copy()
    merged["PARITY_AGENT_CLI"] = "agent"
    merged["PARITY_CONTROLLER_MODEL"] = "grok-4.7-high-fast"
    merged["PARITY_CLAUDE_BUDGET"] = str(BUDGET)
    if env:
        merged.update(env)
    return subprocess.run(
        [sys.executable, str(SUPERVISOR), "--repo", str(root), *args],
        text=True,
        capture_output=True,
        timeout=timeout,
        env=merged,
        check=False,
    )


class QueueRules(unittest.TestCase):
    def test_2_and_3_claude_block_does_not_stop_ready_work(self) -> None:
        tasks = runtime.production_queue(Path("/tmp/prokadiff-orch"), claude_ok=False)
        by_id = runtime.index_tasks(tasks)
        self.assertEqual(by_id["WAVE11B1_SCIENTIFIC_REVIEW"]["status"], runtime.BLOCKED_EXTERNAL)
        ready = [task["task_id"] for task in runtime.eligible_tasks(tasks)]
        self.assertIn("WAVE11B0_SHADOW_AS", ready)
        self.assertIn("WAVE11B2_3072550_REGRESSION", ready)
        self.assertNotIn("WAVE11B1_SCIENTIFIC_REVIEW", ready)
        self.assertNotIn("WAVE11E_PRODUCTION_AS", ready)
        self.assertEqual(runtime.select_decision(runtime.eligible_tasks(tasks)), "WAVE11B0_SHADOW_AS")

    def test_6_stale_lock_text_is_not_ownership(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "writer.lock"
            path.write_text("pid=999999\nholder=stale\n")
            lease = runtime.try_lease(path)
            self.assertIn(f"pid={os.getpid()}", path.read_text())
            lease.release()

    def test_10_second_writer_lease_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "writer.lock"
            first = runtime.try_lease(path)
            with self.assertRaises(runtime.LockBusy):
                runtime.try_lease(path)
            first.release()
            second = runtime.try_lease(path)
            second.release()

    def test_9_cursor_argv_never_selects_gemini(self) -> None:
        argv = runtime.cursor_launch_argv("ping", cli="agent", model="grok-4.7-high-fast")
        self.assertNotIn("agy", argv)
        self.assertNotIn("--force", argv)
        self.assertIn("grok-4.7-high-fast", argv)
        with self.assertRaises(runtime.GeminiDispatchRefused):
            runtime.cursor_launch_argv("ping", cli="agy", model="grok-4.7-high-fast")
        with self.assertRaises(runtime.GeminiDispatchRefused):
            runtime.cursor_launch_argv("ping", cli="agent", model="gemini-3.8-flash-medium")

    def test_json_result_objects_parse(self) -> None:
        single = json.dumps({"type": "result", "result": "ok", "session_id": "s"})
        parsed = runtime.parse_agent_stdout(single)
        self.assertEqual(parsed["result"], "ok")
        lines = '\n'.join(
            [
                json.dumps({"type": "system"}),
                json.dumps({"conversation_id": "c", "response": "done", "status": "SUCCESS"}),
            ]
        )
        parsed = runtime.parse_agent_stdout(lines)
        self.assertEqual(parsed["result"], "done")
        self.assertEqual(parsed["session_id"], "c")

    def test_stale_claude_balance_is_not_inferred_from_spend(self) -> None:
        ok, reason = runtime.claude_balance_ok(BUDGET)
        self.assertFalse(ok)
        self.assertIn("stale", reason)

    def test_13_provider_rate_limit_is_external_block_not_terminal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "worker_stdout.txt"
            log.write_text(
                '{"type":"error","error":{"type":"provider.rate-limit","status":429,'
                '"message":"已达到 5 小时的使用上限。您的限额将在 2026-11-02 16:49:43 重置。"}}\n'
            )
            task = runtime.task_template(
                "B2",
                "analysis",
                "opencode/glm-5.3",
                maximum_retries=1,
                command=["/bin/false"],
                rate_limit_logs=[str(log)],
            )
            resume_at = runtime.detect_rate_limit(task)
            self.assertTrue(resume_at.startswith("2026-11-02T16:49:43"), resume_at)
            runtime.classify_rate_limited(task, 1, resume_at)
            self.assertEqual(task["status"], runtime.BLOCKED_EXTERNAL)
            self.assertEqual(task["retry_count"], 0)
            self.assertEqual(runtime.eligible_tasks([task]), [])
            self.assertEqual(runtime.blocked_summary([task])[0]["blocked_reason"], "worker_rate_limit exit=1")
            # A quota block must not consume the retry budget.
            runtime.refresh_rate_limits([task], log, now=datetime.fromisoformat("2026-11-02T16:00:00+08:00"))
            self.assertEqual(task["status"], runtime.BLOCKED_EXTERNAL)
            runtime.refresh_rate_limits([task], log, now=datetime.fromisoformat("2026-11-02T17:00:00+08:00"))
            self.assertEqual(task["status"], runtime.READY)
            self.assertEqual([t["task_id"] for t in runtime.eligible_tasks([task])], ["B2"])
            # Ordinary failures still become terminal after the retry budget.
            plain = runtime.task_template("X", "analysis", "w", maximum_retries=1, command=["/bin/false"])
            runtime.classify_exit(plain, 1)
            self.assertEqual(plain["status"], runtime.FAILED_TERMINAL)

    def test_14_glm_quota_escalates_to_grok_instead_of_waiting(self) -> None:
        """fastBact policy: GLM primary, Cursor grok takes over on rate limit."""
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "sup.log"
            quota = Path(tmp) / "glm.log"
            quota.write_text(
                '{"type":"error","error":{"type":"provider.rate-limit","status":429,'
                '"message":"已达到 5 小时的使用上限。您的限额将在 2026-11-02 16:49:43 重置。"}}\n'
            )
            task = runtime.task_template(
                "B2",
                "analysis",
                "opencode/glm-5.3",
                maximum_retries=2,
                command=["bash", "glm.sh"],
                rate_limit_logs=[str(quota)],
                fallback_worker="cursor/grok-4.7-high-fast",
                fallback_command=["bash", "grok.sh"],
                fallback_rate_limit_logs=[str(Path(tmp) / "grok.log")],
            )
            resume_at = runtime.detect_rate_limit(task)
            runtime.classify_rate_limited(task, 1, resume_at)
            self.assertTrue(runtime.escalate_worker(task, log, resume_at))
            self.assertEqual(task["status"], runtime.READY)
            self.assertEqual(task["worker"], "cursor/grok-4.7-high-fast")
            self.assertEqual(task["command"], ["bash", "grok.sh"])
            self.assertEqual(task["fallback_rate_limit_logs"], [str(Path(tmp) / "grok.log")])
            self.assertEqual(task["retry_count"], 0)
            self.assertEqual(task["escalation_history"][0]["reason"], "provider_rate_limit")
            self.assertIn("ESCALATE_WORKER B2", log.read_text())
            # The escalated task must run immediately, not wait for the GLM reset.
            self.assertEqual([t["task_id"] for t in runtime.eligible_tasks([task])], ["B2"])
            runtime.refresh_rate_limits([task], log, now=datetime.fromisoformat("2026-11-02T17:30:00+08:00"))
            self.assertEqual(task["worker"], "cursor/grok-4.7-high-fast")
            self.assertEqual(task["status"], runtime.READY)
            # A second rate limit on the fallback must not bounce back to GLM.
            self.assertFalse(runtime.escalate_worker(task, log, resume_at))
            self.assertEqual(task["worker"], "cursor/grok-4.7-high-fast")

    def test_15_production_queue_declares_grok_failover(self) -> None:
        tasks = runtime.production_queue(Path("/tmp/prokadiff-orch"), claude_ok=False)
        by_id = runtime.index_tasks(tasks)
        for task_id in (
            "WAVE11B0_SHADOW_AS",
            "WAVE11B0R_ANALYZER_REPAIR",
            "WAVE11B2_3072550_REGRESSION",
        ):
            task = by_id[task_id]
            self.assertEqual(task["fallback_worker"], "cursor/grok-4.7-high-fast")
            self.assertTrue(task["fallback_command"], task_id)
            self.assertNotIn("--force", " ".join(task["fallback_command"]))
        argv = runtime.cursor_launch_argv("ping", cli="cursor-agent", model="grok-4.7-high-fast")
        self.assertIn("grok-4.7-high-fast", argv)
        with self.assertRaises(runtime.GeminiDispatchRefused):
            runtime.cursor_launch_argv("ping", cli="cursor-agent", model="gemini-3.8-flash-medium")

    def test_17_rejected_review_is_not_complete(self) -> None:
        """Exit 0 from a reviewer must not open the production-integration gate."""
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "sup.log"
            review = Path(tmp) / "review.md"
            review.write_text("**WAVE11D verdict: NEEDS_REVISION.** two attribution defects.\n")
            task = runtime.task_template(
                "WAVE11D_CODEX_REVIEW",
                "independent_review",
                "codex/gpt-6.1-sol",
                status=runtime.WAITING_WORKER,
                command=["/bin/true"],
                output_paths=[str(review)],
            )
            runtime.classify_exit(task, 0)
            self.assertEqual(task["status"], runtime.COMPLETE)
            changed = runtime.apply_review_verdict(task, log)
            self.assertTrue(changed)
            self.assertEqual(task["status"], runtime.REVIEW_REQUIRED)
            self.assertIn("NEEDS_REVISION", task["blocked_reason"])
            self.assertEqual(task["review_verdict"], "NEEDS_REVISION")
            by_id = runtime.index_tasks([task])
            self.assertNotEqual(task["status"], runtime.COMPLETE)
            self.assertFalse(runtime.dependencies_met({"dependencies": ["WAVE11D_CODEX_REVIEW"]}, by_id))

            # A PASS verdict keeps COMPLETE and opens the gate.
            review.write_text("WAVE11D_VERDICT: PASS\nall checks reproduced\n")
            approved = runtime.task_template(
                "WAVE11D_CODEX_REVIEW",
                "independent_review",
                "codex/gpt-6.1-sol",
                status=runtime.WAITING_WORKER,
                command=["/bin/true"],
                output_paths=[str(review)],
            )
            runtime.classify_exit(approved, 0)
            self.assertFalse(runtime.apply_review_verdict(approved, log))
            self.assertEqual(approved["status"], runtime.COMPLETE)
            self.assertEqual(approved["review_verdict"], "PASS")

            # No identifiable verdict must not silently count as approval.
            review.write_text("I looked at the files and they seem fine.\n")
            silent = runtime.task_template(
                "WAVE11D_CODEX_REVIEW",
                "independent_review",
                "codex/gpt-6.1-sol",
                status=runtime.WAITING_WORKER,
                command=["/bin/true"],
                output_paths=[str(review)],
            )
            runtime.classify_exit(silent, 0)
            self.assertTrue(runtime.apply_review_verdict(silent, log))
            self.assertEqual(silent["status"], runtime.REVIEW_REQUIRED)
            self.assertEqual(silent["review_verdict"], "UNSPECIFIED")

    def test_19_rejected_review_requeues_after_upstream_repair(self) -> None:
        """This is the stall the live queue hit: REVIEW_REQUIRED was not eligible."""
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "sup.log"
            dep = runtime.task_template("B0R", "analysis", "opencode/glm-5.3", status=runtime.WAITING_WORKER, command=["/bin/true"])
            review = runtime.task_template(
                "D",
                "independent_review",
                "codex/gpt-6.1-sol",
                status=runtime.REVIEW_REQUIRED,
                blocked_reason="review_verdict=NEEDS_REVISION",
                dependencies=["B0R"],
                command=["bash", "codex.sh"],
            )
            tasks = [dep, review]
            # Repair still running -> reviewer must NOT requeue yet.
            self.assertEqual(runtime.requeue_review_tasks(tasks, log), [])
            self.assertEqual(review["status"], runtime.REVIEW_REQUIRED)
            # Repair completes -> reviewer requeues and becomes dispatchable.
            dep["status"] = runtime.COMPLETE
            self.assertEqual(len(runtime.requeue_review_tasks(tasks, log)), 1)
            self.assertEqual(review["status"], runtime.READY)
            self.assertEqual(review["blocked_reason"], "")
            self.assertIn("D", [t["task_id"] for t in runtime.eligible_tasks(tasks)])
            self.assertIn("REQUEUE_REVIEW D", log.read_text())
            # A review that already finished after its inputs must not be dispatched again.
            review["status"] = runtime.REVIEW_REQUIRED
            review["blocked_reason"] = "review_verdict=NEEDS_REVISION"
            review["completion_timestamp"] = "2026-11-03T14:03:44+0800"
            dep["completion_timestamp"] = "2026-11-03T11:54:13+0800"
            self.assertEqual(runtime.requeue_review_tasks(tasks, log), [])
            self.assertEqual(review["status"], runtime.REVIEW_REQUIRED)
            dep["completion_timestamp"] = "2026-11-03T15:00:00+0800"
            self.assertEqual(len(runtime.requeue_review_tasks(tasks, log)), 1)

    def test_18_analyzer_repair_gates_codex_review(self) -> None:
        tasks = runtime.production_queue(Path("/tmp/prokadiff-orch"), claude_ok=False)
        by_id = runtime.index_tasks(tasks)
        self.assertIn("WAVE11B0R_ANALYZER_REPAIR", by_id)
        self.assertIn("WAVE11B0R_ANALYZER_REPAIR", by_id["WAVE11D_CODEX_REVIEW"]["dependencies"])
        self.assertIn("WAVE11B0R_ANALYZER_REPAIR", by_id["WAVE11B1_SCIENTIFIC_REVIEW"]["dependencies"])
        self.assertIn("WAVE11B0R_ANALYZER_REPAIR", by_id["WAVE11E_PRODUCTION_AS"]["dependencies"])
        # The repair must not be eligible until the reviewed baseline completed.
        by_id["WAVE11B0_SHADOW_AS"]["status"] = runtime.WAITING_WORKER
        self.assertNotIn(
            "WAVE11B0R_ANALYZER_REPAIR",
            [t["task_id"] for t in runtime.eligible_tasks(list(by_id.values()))],
        )
        by_id["WAVE11B0_SHADOW_AS"]["status"] = runtime.COMPLETE
        self.assertIn(
            "WAVE11B0R_ANALYZER_REPAIR",
            [t["task_id"] for t in runtime.eligible_tasks(list(by_id.values()))],
        )

    def test_16_known_glm_rate_limit_skips_glm_before_dispatch(self) -> None:
        """fastBact pre-dispatch rule: check the shared pool first, do not burn a turn."""
        with tempfile.TemporaryDirectory() as tmp:
            quota = Path(tmp) / "quota-state.json"
            quota.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "pools": {
                            "bigmodel-coding-plan": {
                                "status": "RATE_LIMITED",
                                "next_probe_at_utc": "2026-11-02T08:49:43Z",
                            }
                        },
                    }
                )
            )
            log = Path(tmp) / "sup.log"
            task = runtime.task_template(
                "B2",
                "analysis",
                "opencode/glm-5.3",
                command=["bash", "glm.sh"],
                fallback_worker="cursor/grok-4.7-high-fast",
                fallback_command=["bash", "grok.sh"],
            )
            # Cooldown still running -> go straight to Grok.
            blocked, why = runtime.glm_pool_blocked(quota, now=datetime.fromisoformat("2026-11-02T16:00:00+08:00"))
            self.assertTrue(blocked)
            self.assertIn("RATE_LIMITED_UNTIL", why)
            self.assertTrue(runtime.pre_escalate_if_pool_blocked(task, quota, log, now=datetime.fromisoformat("2026-11-02T16:00:00+08:00")))
            self.assertEqual(task["worker"], "cursor/grok-4.7-high-fast")
            self.assertEqual(task["escalation_history"][0]["reason"], "provider_quota_state")
            # Cooldown elapsed -> GLM is allowed again for a retry probe.
            blocked, _ = runtime.glm_pool_blocked(quota, now=datetime.fromisoformat("2026-11-02T17:00:00+08:00"))
            self.assertFalse(blocked)
            fresh = runtime.task_template(
                "B0",
                "analysis",
                "opencode/glm-5.3",
                command=["bash", "glm.sh"],
                fallback_worker="cursor/grok-4.7-high-fast",
                fallback_command=["bash", "grok.sh"],
            )
            self.assertFalse(runtime.pre_escalate_if_pool_blocked(fresh, quota, log, now=datetime.fromisoformat("2026-11-02T17:00:00+08:00")))
            self.assertEqual(fresh["worker"], "opencode/glm-5.3")
            # Unreadable quota state is UNKNOWN, never treated as available.
            self.assertEqual(runtime.quota_pool_status(Path(tmp) / "absent.json")["status"], "UNKNOWN")
            available = runtime.task_template(
                "C", "analysis", "opencode/glm-5.3", fallback_command=["bash", "grok.sh"]
            )
            self.assertFalse(runtime.pre_escalate_if_pool_blocked(available, Path(tmp) / "absent.json", log))
            # A Grok worker must never be pre-escalated onto itself.
            grok = runtime.task_template(
                "D",
                "analysis",
                "cursor/grok-4.7-high-fast",
                fallback_command=["bash", "grok.sh"],
            )
            self.assertFalse(runtime.pre_escalate_if_pool_blocked(grok, quota, log))


class SupervisorProcess(unittest.TestCase):
    def test_1_worker_completion_starts_next_cycle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root, state = make_repo(Path(tmp))
            marker = state / "second.txt"
            put_tasks(
                state,
                [
                    runtime.task_template("A", "synthetic", "fake", priority=1, command=["/bin/true"]),
                    runtime.task_template(
                        "B",
                        "synthetic",
                        "fake",
                        priority=2,
                        dependencies=["A"],
                        command=["/bin/sh", "-c", f"echo ok > {marker}"],
                    ),
                ],
            )
            proc = run_sup(root, ["--loop", "--sleep", "0.1"], {"PARITY_MAX_CYCLES": "8"})
            self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
            data = json.loads((state / "parity_closure_state.json").read_text())
            by_id = runtime.index_tasks(data["task_queue"]["tasks"])
            self.assertEqual(by_id["A"]["status"], runtime.COMPLETE)
            self.assertEqual(by_id["B"]["status"], runtime.COMPLETE)
            self.assertTrue(marker.is_file())
            log = (state / "supervisor.log").read_text()
            self.assertGreaterEqual(log.count("CONTROLLER_CYCLE"), 2)
            self.assertLess(log.index("WORKER_EXIT A"), log.index("CONTROLLER_CYCLE n=2"))

    def test_4_slurm_completion_collects_and_dispatches(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root, state = make_repo(Path(tmp))
            stub = state / "slurm_stub.json"
            write_json(stub, {"42": {"state": "COMPLETED", "exit_code": "0:0"}})
            marker = state / "after_slurm.txt"
            slurm = runtime.task_template("S", "slurm", "slurm", priority=1)
            slurm["status"] = runtime.WAITING_SLURM
            slurm["slurm_job_id"] = "42"
            nxt = runtime.task_template(
                "N",
                "synthetic",
                "fake",
                priority=2,
                dependencies=["S"],
                command=["/bin/sh", "-c", f"echo collected > {marker}"],
            )
            put_tasks(state, [slurm, nxt])
            proc = run_sup(
                root,
                ["--loop", "--sleep", "0.1"],
                {"PARITY_MAX_CYCLES": "6", "PARITY_SLURM_STUB": str(stub)},
            )
            self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
            data = json.loads((state / "parity_closure_state.json").read_text())
            by_id = runtime.index_tasks(data["task_queue"]["tasks"])
            self.assertEqual(by_id["S"]["status"], runtime.COMPLETE)
            self.assertEqual(by_id["N"]["status"], runtime.COMPLETE)
            log = (state / "supervisor.log").read_text()
            self.assertIn("COLLECT_SLURM S", log)
            self.assertIn("dispatch=N", log)

    def test_5_crash_does_not_duplicate_worker(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root, state = make_repo(Path(tmp))
            token = f"orch-sleep-{os.getpid()}-{time.time_ns()}"
            put_tasks(
                state,
                [runtime.task_template("A", "synthetic", "fake", command=["/bin/sleep", "30"], priority=1)],
            )
            env = os.environ.copy()
            env.update(
                {
                    "PARITY_AGENT_CLI": "agent",
                    "PARITY_CONTROLLER_MODEL": "grok-4.7-high-fast",
                    "PARITY_CLAUDE_BUDGET": str(BUDGET),
                }
            )
            first = subprocess.Popen(
                [sys.executable, str(SUPERVISOR), "--repo", str(root), "--loop", "--sleep", "0.1"],
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
            pid = None
            try:
                deadline = time.time() + 20
                pid = None
                while time.time() < deadline:
                    if (state / "parity_closure_state.json").exists():
                        data = json.loads((state / "parity_closure_state.json").read_text())
                        task = data["task_queue"]["tasks"][0]
                        if task.get("status") == runtime.WAITING_WORKER and task.get("execution_pid"):
                            pid = int(task["execution_pid"])
                            break
                    time.sleep(0.1)
                else:
                    self.fail("first supervisor did not record a live worker")
                os.kill(first.pid, 9)
                first.wait(timeout=5)
                second = run_sup(root, ["--loop", "--sleep", "0.1"], {"PARITY_MAX_CYCLES": "2"})
                self.assertEqual(second.returncode, 0, second.stderr + second.stdout)
                log = (state / "supervisor.log").read_text()
                self.assertEqual(log.count("SPAWN A"), 1)
                self.assertIn("ADOPT_RUNNING_WORKER A", log)
            finally:
                if pid:
                    try:
                        os.killpg(int(pid), 9)
                    except OSError:
                        try:
                            os.kill(int(pid), 9)
                        except OSError:
                            pass

    def test_7_loop_survives_multiple_dispatches(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root, state = make_repo(Path(tmp))
            tasks = []
            previous = None
            for name in ("A", "B", "C"):
                tasks.append(
                    runtime.task_template(
                        name,
                        "synthetic",
                        "fake",
                        priority=ord(name),
                        dependencies=[previous] if previous else [],
                        command=["/bin/true"],
                    )
                )
                previous = name
            put_tasks(state, tasks)
            env = os.environ.copy()
            env.update({"PARITY_AGENT_CLI": "agent", "PARITY_CONTROLLER_MODEL": "grok-4.7-high-fast", "PARITY_CLAUDE_BUDGET": str(BUDGET)})
            proc = subprocess.Popen(
                [sys.executable, str(SUPERVISOR), "--repo", str(root), "--loop", "--sleep", "0.1"],
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
            try:
                deadline = time.time() + 25
                while time.time() < deadline:
                    data = json.loads((state / "parity_closure_state.json").read_text())
                    statuses = [task["status"] for task in data["task_queue"]["tasks"]]
                    if statuses == [runtime.COMPLETE, runtime.COMPLETE, runtime.COMPLETE]:
                        break
                    time.sleep(0.1)
                else:
                    self.fail(json.dumps(data["task_queue"]["tasks"]))
                os.kill(proc.pid, 0)
                log = (state / "supervisor.log").read_text()
                self.assertGreaterEqual(log.count("CONTROLLER_CYCLE"), 3)
                self.assertNotIn("QUEUE_LOOP_END", log)
            finally:
                os.kill(proc.pid, 15)
                proc.wait(timeout=10)

    def test_8_restart_resumes_checkpoint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root, state = make_repo(Path(tmp))
            marker = state / "resumed.txt"
            put_tasks(
                state,
                [
                    runtime.task_template("A", "synthetic", "fake", priority=1, command=["/bin/true"]),
                    runtime.task_template(
                        "B",
                        "synthetic",
                        "fake",
                        priority=2,
                        dependencies=["A"],
                        command=["/bin/sh", "-c", f"echo resumed > {marker}"],
                    ),
                ],
            )
            first = run_sup(root, ["--loop", "--sleep", "0.05"], {"PARITY_MAX_CYCLES": "1"})
            self.assertEqual(first.returncode, 0, first.stderr + first.stdout)
            mid = json.loads((state / "parity_closure_state.json").read_text())
            self.assertEqual(mid["task_queue"]["tasks"][0]["status"], runtime.WAITING_WORKER)
            second = run_sup(root, ["--loop", "--sleep", "0.05"], {"PARITY_MAX_CYCLES": "6"})
            self.assertEqual(second.returncode, 0, second.stderr + second.stdout)
            data = json.loads((state / "parity_closure_state.json").read_text())
            by_id = runtime.index_tasks(data["task_queue"]["tasks"])
            self.assertEqual(by_id["A"]["status"], runtime.COMPLETE)
            self.assertEqual(by_id["B"]["status"], runtime.COMPLETE)
            self.assertEqual(marker.read_text().strip(), "resumed")

    def test_9_process_refuses_gemini_without_launching_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root, state = make_repo(Path(tmp))
            proc = run_sup(
                root,
                ["--loop", "--sleep", "0.1"],
                {"PARITY_AGENT_CLI": "agy", "PARITY_CONTROLLER_MODEL": "gemini-3.8-flash-medium", "PARITY_MAX_CYCLES": "1"},
            )
            self.assertEqual(proc.returncode, 3, proc.stderr + proc.stdout)
            log = (state / "supervisor.log").read_text()
            self.assertIn("GEMINI_DISPATCH_REFUSED", log)
            self.assertNotIn("AGENT_START", log)

    def test_11_failed_worker_is_terminal_and_does_not_stick(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root, state = make_repo(Path(tmp))
            marker = state / "sibling.txt"
            bad = runtime.task_template("BAD", "synthetic", "fake", priority=1, maximum_retries=1, command=["/bin/false"])
            good = runtime.task_template(
                "GOOD",
                "synthetic",
                "fake",
                priority=2,
                command=["/bin/sh", "-c", f"echo sibling > {marker}"],
            )
            put_tasks(state, [bad, good])
            proc = run_sup(root, ["--loop", "--sleep", "0.05"], {"PARITY_MAX_CYCLES": "8"})
            self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
            data = json.loads((state / "parity_closure_state.json").read_text())
            by_id = runtime.index_tasks(data["task_queue"]["tasks"])
            self.assertEqual(by_id["BAD"]["status"], runtime.FAILED_TERMINAL)
            self.assertEqual(by_id["GOOD"]["status"], runtime.COMPLETE)
            self.assertNotEqual(by_id["BAD"]["status"], runtime.WAITING_WORKER)
            self.assertTrue(marker.is_file())

    def test_12_all_external_blocks_keep_the_loop_alive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root, state = make_repo(Path(tmp))
            blocked = runtime.task_template(
                "ONLY",
                "scientific_review",
                "claude",
                status=runtime.BLOCKED_EXTERNAL,
                blocked_reason="claude_gateway_balance_unverified",
                resume_condition="balance_refresh",
                allowed_while_other_tasks_blocked=False,
            )
            put_tasks(state, [blocked])
            env = os.environ.copy()
            env.update({"PARITY_AGENT_CLI": "agent", "PARITY_CONTROLLER_MODEL": "grok-4.7-high-fast", "PARITY_CLAUDE_BUDGET": str(BUDGET)})
            proc = subprocess.Popen(
                [sys.executable, str(SUPERVISOR), "--repo", str(root), "--loop", "--sleep", "0.2"],
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
            try:
                deadline = time.time() + 8
                log = ""
                while time.time() < deadline:
                    os.kill(proc.pid, 0)
                    if (state / "supervisor.log").is_file():
                        log = (state / "supervisor.log").read_text()
                        if "QUEUE_IDLE" in log:
                            break
                    time.sleep(0.1)
                self.assertTrue(log, "supervisor.log was not written while the process stayed alive")
                self.assertIn("QUEUE_IDLE", log)
                self.assertIn("claude_gateway_balance_unverified", log)
                self.assertNotIn("SPAWN", log)
                self.assertIsNone(proc.poll())
            finally:
                os.kill(proc.pid, 15)
                proc.wait(timeout=10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
