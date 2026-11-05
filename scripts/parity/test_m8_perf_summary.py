#!/usr/bin/env python3
"""Unit tests for generate_m8_perf_summary.py (no FASTQ / no Slurm / no network).

Hermetic tests pin the parser, the speedup math, and the TSV schema against
synthetic fixtures. Optional integration tests against the real measured logs
in benchmark/results/clonal_2769249/ skip when those (gitignored) files are
absent.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from generate_m8_perf_summary import (
    COLUMNS,
    TimeVRecord,
    build_rows,
    compute_speedup,
    derive_dataset_and_job,
    extract_threads,
    format_elapsed_human,
    main,
    parse_elapsed,
    parse_time_v,
    write_tsv,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
REAL_DIR = REPO_ROOT / "benchmark" / "results" / "clonal_2769249"
REAL_BRESEQ = REAL_DIR / "breseq.time"
REAL_PROKADIFF = REAL_DIR / "prokadiff.time"

# Verbatim measured values from Slurm job 2769249 (qcpu_18i), used to pin the
# delivered numbers without depending on the (gitignored) real log files.
BRESEQ_ELAPSED_RAW = "44:08.46"
BRESEQ_ELAPSED_SEC = 2648.46
PROKADIFF_ELAPSED_RAW = "8:20.28"
PROKADIFF_ELAPSED_SEC = 500.28

BRESEQ_COMMAND = (
    'breseq -j 8 -o /tmp/out_breseq -r REL606.gbk '
    'SRR030257_1.fastq SRR030257_2.fastq'
)
PROKADIFF_COMMAND = (
    '/opt/prokadiff evidence --ref REL606.gbk '
    '--fastq SRR030257_1.fastq --fastq SRR030257_2.fastq '
    '--threads 8 --outdir /tmp/out_pd --keep-bam'
)


def _time_v_text(
    command: str,
    user: str = "100.00",
    sys: str = "10.00",
    cpu: str = "500",
    elapsed: str = "1:00.00",
    rss_kb: str = "1000000",
    exit_status: str = "0",
) -> str:
    return (
        f"\tCommand being timed: {command}\n"
        f"\tUser time (seconds): {user}\n"
        f"\tSystem time (seconds): {sys}\n"
        f"\tPercent of CPU this job got: {cpu}%\n"
        f"\tElapsed (wall clock) time (h:mm:ss or m:ss): {elapsed}\n"
        f"\tAverage shared text size (kbytes): 0\n"
        f"\tAverage unshared data size (kbytes): 0\n"
        f"\tAverage stack size (kbytes): 0\n"
        f"\tAverage total size (kbytes): 0\n"
        f"\tMaximum resident set size (kbytes): {rss_kb}\n"
        f"\tAverage resident set size (kbytes): 0\n"
        f"\tMajor (requiring I/O) page faults: 0\n"
        f"\tMinor (reclaiming a frame) page faults: 1000\n"
        f"\tVoluntary context switches: 1\n"
        f"\tInvoluntary context switches: 2\n"
        f"\tSwaps: 0\n"
        f"\tFile system inputs: 0\n"
        f"\tFile system outputs: 0\n"
        f"\tSocket messages sent: 0\n"
        f"\tSocket messages received: 0\n"
        f"\tSignals delivered: 0\n"
        f"\tPage size (bytes): 4096\n"
        f"\tExit status: {exit_status}\n"
    )


def _breseq_log(tmp_path: Path) -> Path:
    path = tmp_path / "breseq.time"
    path.write_text(
        _time_v_text(
            BRESEQ_COMMAND,
            user="13494.64",
            sys="448.26",
            cpu="526",
            elapsed=BRESEQ_ELAPSED_RAW,
            rss_kb="1857580",
        )
    )
    return path


def _prokadiff_log(tmp_path: Path) -> Path:
    path = tmp_path / "prokadiff.time"
    path.write_text(
        _time_v_text(
            PROKADIFF_COMMAND,
            user="2988.79",
            sys="135.09",
            cpu="624",
            elapsed=PROKADIFF_ELAPSED_RAW,
            rss_kb="3420460",
        )
    )
    return path


# ---------------------------------------------------------------- elapsed ----

@pytest.mark.parametrize(
    ("token", "expected"),
    [
        ("44:08.46", 2648.46),
        ("8:20.28", 500.28),
        ("1:02:03", 3723.0),
        ("1:02:03.50", 3723.5),
        ("0:00.12", 0.12),
        ("12.34", 12.34),
        ("59:59.99", 3599.99),
    ],
)
def test_parse_elapsed_formats(token, expected):
    assert parse_elapsed(token) == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("token", ["ab:cd", "1:2:3:4", "-5:00", "", "1:xx.00"])
def test_parse_elapsed_rejects_garbage(token):
    with pytest.raises(ValueError):
        parse_elapsed(token)


@pytest.mark.parametrize(
    ("seconds", "expected"),
    [
        (2648.46, "44m 08.46s"),
        (500.28, "8m 20.28s"),
        (3723.5, "1h 02m 03.50s"),
        (0.12, "0m 00.12s"),
    ],
)
def test_format_elapsed_human(seconds, expected):
    assert format_elapsed_human(seconds) == expected


# ------------------------------------------------------------------ parser ----

def test_parse_time_v_extracts_all_fields(tmp_path):
    log = _breseq_log(tmp_path)
    rec = parse_time_v(log, "breseq")
    assert rec.tool == "breseq"
    assert rec.command == BRESEQ_COMMAND
    assert rec.user_sec == "13494.64"
    assert rec.sys_sec == "448.26"
    assert rec.cpu_percent == "526"
    assert rec.elapsed_raw == BRESEQ_ELAPSED_RAW
    assert rec.elapsed_sec == pytest.approx(BRESEQ_ELAPSED_SEC, abs=1e-9)
    assert rec.peak_rss_kb == "1857580"
    assert rec.exit_status == 0


def test_parse_time_v_missing_field_fails_closed(tmp_path):
    path = tmp_path / "incomplete.time"
    path.write_text(_time_v_text("tool -j 8").replace(
        "\tMaximum resident set size (kbytes): 1000000\n", ""
    ))
    with pytest.raises(ValueError, match="peak_rss_kb"):
        parse_time_v(path, "breseq")


def test_parse_time_v_nonzero_exit_rejected(tmp_path):
    path = tmp_path / "failed.time"
    path.write_text(_time_v_text("tool -j 8", exit_status="1"))
    with pytest.raises(ValueError, match="Exit status 1"):
        parse_time_v(path, "breseq")


def test_parse_time_v_missing_file_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):
        parse_time_v(tmp_path / "absent.time", "breseq")


# ----------------------------------------------------------------- threads ----

def test_extract_threads_breseq_and_prokadiff():
    assert extract_threads(BRESEQ_COMMAND, "breseq") == 8
    assert extract_threads(PROKADIFF_COMMAND, "prokadiff") == 8


def test_extract_threads_equals_form():
    assert extract_threads("prokadiff run --threads=16 --outdir o", "prokadiff") == 16


def test_extract_threads_missing_rejected():
    with pytest.raises(ValueError, match="thread count"):
        extract_threads("breseq -o out reference.gbk reads.fastq", "breseq")
    with pytest.raises(ValueError, match="thread count"):
        extract_threads("prokadiff evidence --ref r.gbk", "prokadiff")


# ------------------------------------------------------- dataset / job id ----

def test_derive_dataset_and_job():
    assert derive_dataset_and_job(Path("/x/benchmark/results/clonal_2769249")) == (
        "clonal",
        "2769249",
    )
    with pytest.raises(ValueError):
        derive_dataset_and_job(Path("/x/benchmark/results/ambiguous_base"))


# ---------------------------------------------------------------- speedup ----

def test_speedup_math_matches_measured_5_29():
    ratio = compute_speedup(BRESEQ_ELAPSED_SEC, PROKADIFF_ELAPSED_SEC)
    assert ratio == pytest.approx(5.293955, abs=1e-5)
    assert f"{ratio:.2f}" == "5.29"


@pytest.mark.parametrize(
    ("breseq_s", "prokadiff_s"),
    [(0.0, 10.0), (10.0, 0.0), (-1.0, 10.0)],
)
def test_speedup_requires_positive_elapsed(breseq_s, prokadiff_s):
    with pytest.raises(ValueError):
        compute_speedup(breseq_s, prokadiff_s)


# ------------------------------------------------------------- rows/schema ----

def _make_records(tmp_path: Path) -> tuple[TimeVRecord, TimeVRecord]:
    return (
        parse_time_v(_breseq_log(tmp_path), "breseq"),
        parse_time_v(_prokadiff_log(tmp_path), "prokadiff"),
    )


def test_build_rows_schema_and_measured_values(tmp_path):
    breseq, prokadiff = _make_records(tmp_path)
    rows = build_rows(breseq, prokadiff, "clonal", "2769249")
    assert [r["tool"] for r in rows] == ["breseq", "prokadiff"]  # fixed order

    br, pd = rows
    for row in rows:
        assert set(row) == set(COLUMNS)  # exact schema
    assert br["dataset"] == pd["dataset"] == "clonal"
    assert br["job_id"] == pd["job_id"] == "2769249"
    assert br["threads"] == pd["threads"] == "8"
    assert br["elapsed_sec"] == "2648.46"
    assert br["elapsed_human"] == "44m 08.46s"
    assert br["user_sec"] == "13494.64"
    assert br["sys_sec"] == "448.26"
    assert br["cpu_percent"] == "526"
    assert br["peak_rss_kb"] == "1857580"
    assert br["peak_rss_mb"] == f"{1857580 / 1024:.2f}"
    assert br["speedup_vs_breseq"] == "1.00"
    assert pd["elapsed_sec"] == "500.28"
    assert pd["elapsed_human"] == "8m 20.28s"
    assert pd["user_sec"] == "2988.79"
    assert pd["sys_sec"] == "135.09"
    assert pd["cpu_percent"] == "624"
    assert pd["peak_rss_kb"] == "3420460"
    assert pd["peak_rss_mb"] == f"{3420460 / 1024:.2f}"
    assert pd["speedup_vs_breseq"] == "5.29"
    # every row cites its source log (provenance, zero fabrication)
    assert "source=" in br["notes"] and str(breseq.log_path) in br["notes"]
    assert "source=" in pd["notes"] and str(prokadiff.log_path) in pd["notes"]
    # FIX-M8-T05 labeling: breseq is the baseline operand of the speedup
    # numerator (breseq_elapsed_sec / prokadiff_elapsed_sec), not a denominator
    assert (
        "baseline (speedup numerator: breseq_elapsed_sec / prokadiff_elapsed_sec)"
        in br["notes"]
    )
    assert "speedup denominator" not in br["notes"]
    # FIX-M8-T05 scope: command-level measurement, not an end-to-end claim
    assert "speedup=breseq_elapsed_sec/prokadiff_elapsed_sec" in pd["notes"]
    assert (
        "scope=prokadiff evidence vs breseq pipeline (command-level resource "
        "measurement of recorded job 2769249, not end-to-end pipeline claim)"
        in pd["notes"]
    )


def test_build_rows_rejects_mismatched_threads(tmp_path):
    breseq, _ = _make_records(tmp_path)
    # breseq -j 8 vs prokadiff --threads 4 must not produce a speedup claim
    mismatch_log = _write_text(
        tmp_path,
        "mismatch.time",
        _time_v_text(
            PROKADIFF_COMMAND.replace("--threads 8", "--threads 4"),
            elapsed=PROKADIFF_ELAPSED_RAW,
        ),
    )
    mismatched = parse_time_v(mismatch_log, "prokadiff")
    with pytest.raises(ValueError, match="thread counts not aligned"):
        build_rows(breseq, mismatched, "clonal", "2769249")


def _write_text(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path


def test_notes_extra_appended_and_sanitized(tmp_path):
    breseq, prokadiff = _make_records(tmp_path)
    rows = build_rows(breseq, prokadiff, "clonal", "2769249", extra_notes="node=g01\tqcpu_18i\nline2")
    for row in rows:
        assert "node=g01 qcpu_18i line2" in row["notes"]
        assert "\t" not in row["notes"] and "\n" not in row["notes"]


# --------------------------------------------------------------- end-to-end ----

def test_main_end_to_end_writes_valid_tsv(tmp_path, monkeypatch):
    job_dir = tmp_path / "clonal_2769249"
    job_dir.mkdir()
    (job_dir / "breseq.time").write_text(
        _time_v_text(
            BRESEQ_COMMAND,
            user="13494.64",
            sys="448.26",
            cpu="526",
            elapsed=BRESEQ_ELAPSED_RAW,
            rss_kb="1857580",
        )
    )
    (job_dir / "prokadiff.time").write_text(
        _time_v_text(
            PROKADIFF_COMMAND,
            user="2988.79",
            sys="135.09",
            cpu="624",
            elapsed=PROKADIFF_ELAPSED_RAW,
            rss_kb="3420460",
        )
    )
    out = tmp_path / "m8_perf_summary.tsv"
    rc = main(
        [
            "--breseq-time", str(job_dir / "breseq.time"),
            "--prokadiff-time", str(job_dir / "prokadiff.time"),
            "--output", str(out),
            "--note", "fixture=hermetic",
        ]
    )
    assert rc == 0
    with out.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == COLUMNS
        rows = list(reader)
    assert len(rows) == 2
    assert rows[0]["tool"] == "breseq" and rows[1]["tool"] == "prokadiff"
    assert rows[0]["dataset"] == rows[1]["dataset"] == "clonal"
    assert rows[0]["job_id"] == rows[1]["job_id"] == "2769249"
    assert rows[0]["speedup_vs_breseq"] == "1.00"
    assert rows[1]["speedup_vs_breseq"] == "5.29"
    assert rows[1]["elapsed_sec"] == "500.28"
    assert all("fixture=hermetic" in r["notes"] for r in rows)
    # FIX-M8-T05: corrected label + explicit command scope survive end-to-end
    assert (
        "baseline (speedup numerator: breseq_elapsed_sec / prokadiff_elapsed_sec)"
        in rows[0]["notes"]
    )
    assert "scope=prokadiff evidence vs breseq pipeline" in rows[1]["notes"]
    assert "not end-to-end pipeline claim" in rows[1]["notes"]


def test_write_tsv_is_deterministic(tmp_path):
    breseq, prokadiff = _make_records(tmp_path)
    rows = build_rows(breseq, prokadiff, "clonal", "2769249")
    out1, out2 = tmp_path / "a.tsv", tmp_path / "b.tsv"
    write_tsv(rows, out1)
    write_tsv(rows, out2)
    assert out1.read_bytes() == out2.read_bytes()


# ------------------------------------------------- real measured logs (opt) ----

@pytest.mark.skipif(
    not (REAL_BRESEQ.is_file() and REAL_PROKADIFF.is_file()),
    reason="real qcpu_18i time -v logs for job 2769249 not present in this checkout",
)
def test_real_logs_match_delivered_values(tmp_path):
    breseq = parse_time_v(REAL_BRESEQ, "breseq")
    prokadiff = parse_time_v(REAL_PROKADIFF, "prokadiff")
    # Pin the delivered measurement: if the logs are replaced by a re-measure,
    # this guard forces m8_perf_summary.tsv and this test to be regenerated.
    assert breseq.elapsed_raw == BRESEQ_ELAPSED_RAW
    assert prokadiff.elapsed_raw == PROKADIFF_ELAPSED_RAW
    assert breseq.peak_rss_kb == "1857580"
    assert prokadiff.peak_rss_kb == "3420460"
    out = tmp_path / "real.tsv"
    rc = main(
        [
            "--breseq-time", str(REAL_BRESEQ),
            "--prokadiff-time", str(REAL_PROKADIFF),
            "--output", str(out),
        ]
    )
    assert rc == 0
    with out.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[1]["speedup_vs_breseq"] == "5.29"
    delivered = REAL_DIR.parent / "m8_perf_summary.tsv"
    if delivered.is_file():
        assert delivered.read_bytes() == out.read_bytes(), (
            "benchmark/results/m8_perf_summary.tsv is stale; rerun "
            "scripts/parity/generate_m8_perf_summary.py"
        )
