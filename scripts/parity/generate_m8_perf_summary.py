#!/usr/bin/env python3
"""generate_m8_perf_summary.py: Extract measured /usr/bin/time -v logs into a perf TSV.

Parses existing Slurm qcpu_18i measurement logs (e.g.
``benchmark/results/clonal_2769249/{breseq,prokadiff}.time``) and writes
``benchmark/results/m8_perf_summary.tsv``.

Zero-fabrication contract (ProkaDiff Rules 5 & 15):
- every numeric cell is parsed directly from a ``/usr/bin/time -v`` log;
- ``speedup_vs_breseq`` is the only derived column and is computed as
  ``breseq_elapsed_sec / prokadiff_elapsed_sec``;
- missing fields, unparsable values, or a non-zero ``Exit status`` abort the
  run with no output file (a failed measurement can never support a perf
  claim);
- nothing is defaulted, interpolated, or remembered from previous runs.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

COLUMNS = [
    "tool",
    "dataset",
    "job_id",
    "threads",
    "elapsed_sec",
    "elapsed_human",
    "user_sec",
    "sys_sec",
    "cpu_percent",
    "peak_rss_kb",
    "peak_rss_mb",
    "speedup_vs_breseq",
    "notes",
]

# GNU time -v line patterns. Every value is captured verbatim from the log.
TIME_V_PATTERNS = {
    "command": re.compile(r"^\tCommand being timed: (.+)$"),
    "user_sec": re.compile(r"^\tUser time \(seconds\): ([0-9]+(?:\.[0-9]+)?)$"),
    "sys_sec": re.compile(r"^\tSystem time \(seconds\): ([0-9]+(?:\.[0-9]+)?)$"),
    "cpu_percent": re.compile(r"^\tPercent of CPU this job got: ([0-9]+)%$"),
    "elapsed_raw": re.compile(
        r"^\tElapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)$"
    ),
    "peak_rss_kb": re.compile(r"^\tMaximum resident set size \(kbytes\): ([0-9]+)$"),
    "exit_status": re.compile(r"^\tExit status: ([0-9]+)$"),
}

THREADS_PATTERNS = {
    "breseq": re.compile(r"(?:^|\s)-j(?:\s+|=)(\d+)"),
    "prokadiff": re.compile(r"(?:^|\s)--threads(?:\s+|=)(\d+)"),
}


@dataclass(frozen=True)
class TimeVRecord:
    """One parsed /usr/bin/time -v log. All values verbatim from the file."""

    tool: str
    log_path: Path
    command: str
    user_sec: str
    sys_sec: str
    cpu_percent: str
    elapsed_raw: str
    elapsed_sec: float
    peak_rss_kb: str
    exit_status: int


def parse_elapsed(text: str) -> float:
    """Convert a GNU time elapsed token to seconds.

    Accepted forms: ``h:mm:ss.cc``, ``m:ss.cc`` (fraction optional), or a bare
    ``ss.cc`` second count. Anything else raises.
    """
    token = text.strip()
    parts = token.split(":")
    try:
        if len(parts) == 1:
            hours, minutes, seconds = 0, 0, float(parts[0])
        elif len(parts) == 2:
            hours, minutes, seconds = 0, int(parts[0]), float(parts[1])
        elif len(parts) == 3:
            hours, minutes, seconds = int(parts[0]), int(parts[1]), float(parts[2])
        else:
            raise ValueError(f"unparsable elapsed token: {text!r}")
    except ValueError as exc:
        raise ValueError(f"unparsable elapsed token: {text!r}") from exc
    if hours < 0 or minutes < 0 or seconds < 0:
        raise ValueError(f"negative elapsed token: {text!r}")
    return hours * 3600.0 + minutes * 60.0 + seconds


def format_elapsed_human(seconds: float) -> str:
    """Render seconds as ``44m 08.46s`` or ``1h 02m 03.50s``."""
    if seconds < 0:
        raise ValueError(f"negative elapsed: {seconds!r}")
    hours, rem = divmod(seconds, 3600.0)
    minutes, secs = divmod(rem, 60.0)
    if hours >= 1:
        return f"{int(hours)}h {int(minutes):02d}m {secs:05.2f}s"
    return f"{int(minutes)}m {secs:05.2f}s"


def _sanitize(text: str) -> str:
    """Keep TSV cells single-line and tab-free."""
    return text.replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


def parse_time_v(path: Path, tool: str) -> TimeVRecord:
    """Parse one /usr/bin/time -v log, failing closed on any defect."""
    if not path.is_file():
        raise FileNotFoundError(f"missing time -v log: {path}")
    found: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            for key, pattern in TIME_V_PATTERNS.items():
                match = pattern.match(line.rstrip("\n"))
                if match:
                    found[key] = match.group(1)
                    break
    missing = [key for key in TIME_V_PATTERNS if key not in found]
    if missing:
        raise ValueError(f"{path}: missing time -v fields: {', '.join(sorted(missing))}")
    exit_status = int(found["exit_status"])
    if exit_status != 0:
        raise ValueError(
            f"{path}: Exit status {exit_status} != 0; failed measurement "
            "cannot support a performance claim (Rules 5 & 15)"
        )
    elapsed_sec = parse_elapsed(found["elapsed_raw"])
    if elapsed_sec <= 0:
        raise ValueError(f"{path}: non-positive elapsed time {found['elapsed_raw']!r}")
    if int(found["peak_rss_kb"]) <= 0:
        raise ValueError(f"{path}: non-positive peak RSS {found['peak_rss_kb']!r}")
    return TimeVRecord(
        tool=tool,
        log_path=path,
        command=_sanitize(found["command"]),
        user_sec=found["user_sec"],
        sys_sec=found["sys_sec"],
        cpu_percent=found["cpu_percent"],
        elapsed_raw=found["elapsed_raw"],
        elapsed_sec=elapsed_sec,
        peak_rss_kb=found["peak_rss_kb"],
        exit_status=exit_status,
    )


def extract_threads(command: str, tool: str) -> int:
    """Pull the aligned thread count out of the timed command line."""
    pattern = THREADS_PATTERNS.get(tool)
    if pattern is None:
        raise ValueError(f"no thread-flag pattern registered for tool {tool!r}")
    match = pattern.search(command)
    if not match:
        raise ValueError(
            f"{tool}: could not locate thread count ({'-j' if tool == 'breseq' else '--threads'}) "
            f"in timed command: {command!r}"
        )
    threads = int(match.group(1))
    if threads <= 0:
        raise ValueError(f"{tool}: non-positive thread count {threads}")
    return threads


def derive_dataset_and_job(results_dir: Path) -> tuple[str, str]:
    """Split a results dir named ``<dataset>_<jobid>`` (e.g. clonal_2769249)."""
    match = re.match(r"^(?P<dataset>.+)_(?P<job>\d+)$", results_dir.name)
    if not match:
        raise ValueError(
            f"cannot derive dataset/job_id from directory name {results_dir.name!r}; "
            "pass --dataset and --job-id explicitly"
        )
    return match.group("dataset"), match.group("job")


def compute_speedup(breseq_elapsed_sec: float, prokadiff_elapsed_sec: float) -> float:
    """breseq wall-clock divided by prokadiff wall-clock; both strictly > 0."""
    if breseq_elapsed_sec <= 0 or prokadiff_elapsed_sec <= 0:
        raise ValueError("speedup requires strictly positive elapsed times on both sides")
    return breseq_elapsed_sec / prokadiff_elapsed_sec


def build_rows(
    breseq: TimeVRecord,
    prokadiff: TimeVRecord,
    dataset: str,
    job_id: str,
    extra_notes: str = "",
) -> list[dict[str, str]]:
    """Assemble the two summary rows. Order is fixed: breseq baseline first."""
    breseq_threads = extract_threads(breseq.command, "breseq")
    prokadiff_threads = extract_threads(prokadiff.command, "prokadiff")
    speedup = compute_speedup(breseq.elapsed_sec, prokadiff.elapsed_sec)

    def notes(record: TimeVRecord, suffix: str) -> str:
        parts = [
            f"source={record.log_path}",
            "measured_by=/usr/bin/time -v",
            suffix,
        ]
        if extra_notes:
            parts.append(_sanitize(extra_notes))
        return "; ".join(parts)

    breseq_row = {
        "tool": "breseq",
        "dataset": dataset,
        "job_id": job_id,
        "threads": str(breseq_threads),
        "elapsed_sec": f"{breseq.elapsed_sec:.2f}",
        "elapsed_human": format_elapsed_human(breseq.elapsed_sec),
        "user_sec": breseq.user_sec,
        "sys_sec": breseq.sys_sec,
        "cpu_percent": breseq.cpu_percent,
        "peak_rss_kb": breseq.peak_rss_kb,
        "peak_rss_mb": f"{int(breseq.peak_rss_kb) / 1024.0:.2f}",
        "speedup_vs_breseq": "1.00",
        "notes": notes(
            breseq,
            "baseline (speedup numerator: breseq_elapsed_sec / prokadiff_elapsed_sec)",
        ),
    }
    prokadiff_row = {
        "tool": "prokadiff",
        "dataset": dataset,
        "job_id": job_id,
        "threads": str(prokadiff_threads),
        "elapsed_sec": f"{prokadiff.elapsed_sec:.2f}",
        "elapsed_human": format_elapsed_human(prokadiff.elapsed_sec),
        "user_sec": prokadiff.user_sec,
        "sys_sec": prokadiff.sys_sec,
        "cpu_percent": prokadiff.cpu_percent,
        "peak_rss_kb": prokadiff.peak_rss_kb,
        "peak_rss_mb": f"{int(prokadiff.peak_rss_kb) / 1024.0:.2f}",
        "speedup_vs_breseq": f"{speedup:.2f}",
        "notes": notes(
            prokadiff,
            "speedup=breseq_elapsed_sec/prokadiff_elapsed_sec; "
            "scope=prokadiff evidence vs breseq pipeline (command-level resource "
            f"measurement of recorded job {job_id}, not end-to-end pipeline claim)",
        ),
    }
    if breseq_threads != prokadiff_threads:
        raise ValueError(
            f"thread counts not aligned (breseq -j {breseq_threads} vs prokadiff "
            f"--threads {prokadiff_threads}); a speedup claim requires the same "
            "thread setting on both sides"
        )
    return [breseq_row, prokadiff_row]


def write_tsv(rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize measured /usr/bin/time -v logs for breseq vs prokadiff "
            "into benchmark/results/m8_perf_summary.tsv. Numbers are parsed "
            "from the logs only; nothing is fabricated (Rules 5 & 15)."
        )
    )
    parser.add_argument(
        "--breseq-time",
        type=Path,
        default=REPO_ROOT / "benchmark/results/clonal_2769249/breseq.time",
        help="path to the breseq /usr/bin/time -v log",
    )
    parser.add_argument(
        "--prokadiff-time",
        type=Path,
        default=REPO_ROOT / "benchmark/results/clonal_2769249/prokadiff.time",
        help="path to the prokadiff /usr/bin/time -v log",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=REPO_ROOT / "benchmark/results/m8_perf_summary.tsv",
        help="output TSV path",
    )
    parser.add_argument(
        "--dataset",
        default=None,
        help="dataset name (default: derived from the results directory name)",
    )
    parser.add_argument(
        "--job-id",
        default=None,
        help="Slurm job id (default: derived from the results directory name)",
    )
    parser.add_argument(
        "--note",
        default="",
        help="extra provenance note appended to every row (e.g. node, partition)",
    )
    args = parser.parse_args(argv)

    breseq = parse_time_v(args.breseq_time, "breseq")
    prokadiff = parse_time_v(args.prokadiff_time, "prokadiff")

    dataset = args.dataset
    job_id = args.job_id
    if dataset is None or job_id is None:
        derived_dataset, derived_job = derive_dataset_and_job(
            args.breseq_time.resolve().parent
        )
        dataset = dataset if dataset is not None else derived_dataset
        job_id = job_id if job_id is not None else derived_job

    rows = build_rows(breseq, prokadiff, dataset, job_id, extra_notes=args.note)
    write_tsv(rows, args.output)

    speedup = compute_speedup(breseq.elapsed_sec, prokadiff.elapsed_sec)
    print(f"[m8_perf_summary] dataset={dataset} job_id={job_id}")
    print(
        f"[m8_perf_summary] breseq   elapsed={breseq.elapsed_sec:.2f}s "
        f"peak_rss={breseq.peak_rss_kb}KB"
    )
    print(
        f"[m8_perf_summary] prokadiff elapsed={prokadiff.elapsed_sec:.2f}s "
        f"peak_rss={prokadiff.peak_rss_kb}KB"
    )
    print(f"[m8_perf_summary] speedup_vs_breseq={speedup:.2f}x")
    print(f"[m8_perf_summary] wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
