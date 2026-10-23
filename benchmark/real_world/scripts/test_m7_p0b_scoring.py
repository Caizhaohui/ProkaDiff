#!/usr/bin/env python3
"""
benchmark/real_world/scripts/test_m7_p0b_scoring.py

Regression suite for Milestone M7-P0b:
1. Quarantined / invalid_legacy truth paths are rejected fail-closed (exit 2)
   by both comparators; --test-gd is never subject to that rejection.
2. Zero-denominator scoring emits the NOT_EVALUABLE sentinel, never 0.0/1.0.
3. Real positive truth with matched calls still emits numeric metrics.
4. aggregate_results.py propagates NOT_EVALUABLE verbatim.
5. PRJNA884016 assets stay negative-control only: manifest repudiation,
   registry booleans false / PENDING, no stale MUMmer sidecars in the active
   truth directory.

Run with: python3 -m pytest -q benchmark/real_world/scripts/test_m7_p0b_scoring.py
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS_DIR = REPO_ROOT / "benchmark" / "real_world" / "scripts"
TRUTH_DIR = REPO_ROOT / "benchmark" / "real_world" / "truth" / "prjna884016"
QUARANTINED_DIR = TRUTH_DIR / "quarantined_legacy"
REGISTRY_PATH = REPO_ROOT / "benchmark" / "real_world" / "VALIDATION_REGISTRY.tsv"
PRJNA884016_MANIFEST = REPO_ROOT / "benchmark" / "real_world" / "manifests" / "prjna884016.tsv"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from aggregate_results import aggregate_dir
from compare_mob_truth import compare_mob
from compare_variant_calls import NOT_EVALUABLE, match_mutations, parse_gd

GD_HEADER = "#=GENOME_DIFF\t1.0\n"
ONE_MOB_RECORD = "MOB\t1\t.\tchr\t1000\tIS1\t1\t9\n"


def _write_gd(tmp_path: Path, name: str, body: str = "") -> Path:
    path = tmp_path / name
    path.write_text(GD_HEADER + body, encoding="utf-8")
    return path


def _run_compare(script: str, test_gd: Path, truth_gd: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / script), "--test-gd", str(test_gd), "--truth-gd", str(truth_gd), *extra],
        capture_output=True,
        text=True,
    )


# ---------------------------------------------------------------------------
# 1. Fail-closed quarantined truth loading
# ---------------------------------------------------------------------------


def test_variant_comparator_rejects_quarantined_truth_paths(tmp_path: Path) -> None:
    plain_test = _write_gd(tmp_path, "test.gd", ONE_MOB_RECORD)
    cases = [
        QUARANTINED_DIR / "SY22A1_truth_invalid_legacy.gd",
        QUARANTINED_DIR / "MANIFEST.tsv",
        _write_gd(tmp_path, "curated_invalid_legacy.gd", ONE_MOB_RECORD),
    ]
    for truth in cases:
        proc = _run_compare("compare_variant_calls.py", plain_test, truth, "--tsv")
        assert proc.returncode == 2, f"expected exit 2 for {truth}, got {proc.returncode}: {proc.stderr}"
        assert "quarantined truth cannot be used as positive truth" in proc.stderr


def test_mob_comparator_rejects_quarantined_truth_paths(tmp_path: Path) -> None:
    plain_test = _write_gd(tmp_path, "test.gd", ONE_MOB_RECORD)
    cases = [
        QUARANTINED_DIR / "SY22A1_truth_invalid_legacy.gd",
        QUARANTINED_DIR / "MANIFEST.tsv",
        _write_gd(tmp_path, "curated_invalid_legacy.gd", ONE_MOB_RECORD),
    ]
    for truth in cases:
        proc = _run_compare("compare_mob_truth.py", plain_test, truth, "--tsv")
        assert proc.returncode == 2, f"expected exit 2 for {truth}, got {proc.returncode}: {proc.stderr}"
        assert "quarantined truth cannot be used as positive truth" in proc.stderr


def test_invalid_legacy_name_is_allowed_as_test_input(tmp_path: Path) -> None:
    """The rejection targets --truth-gd only; a legacy-named --test-gd runs normally."""
    legacy_named_test = _write_gd(tmp_path, "calls_invalid_legacy.gd", ONE_MOB_RECORD)
    plain_truth = _write_gd(tmp_path, "truth.gd", ONE_MOB_RECORD)
    for script in ("compare_variant_calls.py", "compare_mob_truth.py"):
        proc = _run_compare(script, legacy_named_test, plain_truth, "--tsv")
        assert proc.returncode == 0, proc.stderr


def test_ordinary_positive_truth_is_not_rejected(tmp_path: Path) -> None:
    """Ordinary curated truth (BL21-style) must pass the quarantine gate and score numerically."""
    test_gd = _write_gd(tmp_path, "test.gd", ONE_MOB_RECORD)
    bl21_style_truth = _write_gd(tmp_path, "bl21_curated_truth.gd", ONE_MOB_RECORD)
    proc = _run_compare("compare_variant_calls.py", test_gd, bl21_style_truth, "--tsv")
    assert proc.returncode == 0, proc.stderr
    assert "variant_precision\t1.0000" in proc.stdout


# ---------------------------------------------------------------------------
# 2. Zero-denominator scoring sentinels
# ---------------------------------------------------------------------------


def test_empty_truth_and_empty_calls_are_not_evaluable(tmp_path: Path) -> None:
    test_gd = _write_gd(tmp_path, "empty_test.gd")
    truth_gd = _write_gd(tmp_path, "empty_truth.gd")
    metrics = match_mutations(parse_gd(str(test_gd)), parse_gd(str(truth_gd)))
    assert metrics["overall_tp"] == 0
    assert metrics["overall_fp"] == 0
    assert metrics["overall_fn"] == 0
    for key in (
        "variant_precision",
        "variant_recall",
        "variant_f1",
        "breakpoint_exact_fraction",
        "breakpoint_within_1bp",
        "breakpoint_within_5bp",
    ):
        assert metrics[key] == NOT_EVALUABLE, key

    proc = _run_compare("compare_variant_calls.py", test_gd, truth_gd, "--tsv")
    assert proc.returncode == 0, proc.stderr
    for expected in (
        "variant_precision\tNOT_EVALUABLE",
        "variant_recall\tNOT_EVALUABLE",
        "variant_f1\tNOT_EVALUABLE",
        "MOB_precision\tNOT_EVALUABLE",
        "JC_recall\tNOT_EVALUABLE",
        "breakpoint_exact_fraction\tNOT_EVALUABLE",
        "breakpoint_within_1bp\tNOT_EVALUABLE",
        "breakpoint_within_5bp\tNOT_EVALUABLE",
    ):
        assert expected in proc.stdout, expected
    for fabricated in (
        "variant_precision\t1.0000",
        "variant_precision\t0.0000",
        "variant_recall\t1.0000",
        "breakpoint_exact_fraction\t1.0000",
    ):
        assert fabricated not in proc.stdout, fabricated

    human = _run_compare("compare_variant_calls.py", test_gd, truth_gd)
    assert human.returncode == 0, human.stderr
    assert "NOT_EVALUABLE" in human.stdout
    assert "100.00%" not in human.stdout
    assert "0.00%" not in human.stdout


def test_empty_truth_with_one_extra_call(tmp_path: Path) -> None:
    test_gd = _write_gd(tmp_path, "one_call.gd", ONE_MOB_RECORD)
    truth_gd = _write_gd(tmp_path, "empty_truth.gd")
    metrics = match_mutations(parse_gd(str(test_gd)), parse_gd(str(truth_gd)))
    assert metrics["overall_tp"] == 0
    assert metrics["overall_fp"] == 1
    assert metrics["overall_fn"] == 0
    assert metrics["variant_precision"] == 0.0
    assert isinstance(metrics["variant_precision"], float)
    assert metrics["variant_recall"] == NOT_EVALUABLE
    assert metrics["variant_f1"] == NOT_EVALUABLE
    assert metrics["breakpoint_exact_fraction"] == NOT_EVALUABLE
    assert metrics["breakpoint_within_1bp"] == NOT_EVALUABLE
    assert metrics["breakpoint_within_5bp"] == NOT_EVALUABLE

    proc = _run_compare("compare_variant_calls.py", test_gd, truth_gd, "--tsv")
    assert proc.returncode == 0, proc.stderr
    assert "variant_precision\t0.0000" in proc.stdout
    assert "MOB_precision\t0.0000" in proc.stdout
    assert "variant_recall\tNOT_EVALUABLE" in proc.stdout
    assert "MOB_recall\tNOT_EVALUABLE" in proc.stdout
    assert "breakpoint_within_5bp\tNOT_EVALUABLE" in proc.stdout


def test_matched_positive_event_keeps_numeric_scores(tmp_path: Path) -> None:
    test_gd = _write_gd(tmp_path, "matched_test.gd", ONE_MOB_RECORD)
    truth_gd = _write_gd(tmp_path, "matched_truth.gd", ONE_MOB_RECORD)
    metrics = match_mutations(parse_gd(str(test_gd)), parse_gd(str(truth_gd)))
    assert metrics["overall_tp"] == 1
    assert metrics["overall_fp"] == 0
    assert metrics["overall_fn"] == 0
    assert metrics["variant_precision"] == 1.0
    assert metrics["variant_recall"] == 1.0
    assert metrics["variant_f1"] == 1.0
    assert metrics["breakpoint_exact_fraction"] == 1.0
    assert not isinstance(metrics["variant_precision"], str)
    assert not isinstance(metrics["variant_recall"], str)

    proc = _run_compare("compare_variant_calls.py", test_gd, truth_gd, "--tsv")
    assert proc.returncode == 0, proc.stderr
    for numeric in (
        "variant_precision\t1.0000",
        "variant_recall\t1.0000",
        "MOB_recall\t1.0000",
        "breakpoint_exact_fraction\t1.0000",
        "breakpoint_within_5bp\t1.0000",
    ):
        assert numeric in proc.stdout, numeric


def test_mob_comparator_matched_positive_event_is_numeric() -> None:
    call = [{"seq_id": "chr", "pos": 1000, "element": "IS1", "strand": "1", "repeat_size": 9}]
    m = compare_mob(list(call), list(call))
    assert m["precision"] == 1.0
    assert m["recall"] == 1.0
    assert not isinstance(m["precision"], str)
    assert not isinstance(m["recall"], str)


def test_mob_comparator_empty_truth_stays_not_evaluable() -> None:
    empty = compare_mob([], [])
    assert empty["precision"] == NOT_EVALUABLE
    assert empty["recall"] == NOT_EVALUABLE
    one_fp = compare_mob([{"seq_id": "chr", "pos": 100, "element": "IS1", "strand": "1", "repeat_size": 9}], [])
    assert one_fp["precision"] == 0.0
    assert one_fp["recall"] == NOT_EVALUABLE


# ---------------------------------------------------------------------------
# 3. Aggregation pass-through
# ---------------------------------------------------------------------------


def test_aggregate_results_passes_not_evaluable_through(tmp_path: Path, capsys) -> None:
    ds_dir = tmp_path / "results" / "negative_control"
    ds_dir.mkdir(parents=True)
    (ds_dir / "mob_metrics.tsv").write_text(
        "metric\tvalue\n"
        "mob_precision\t0.5000\n"
        "mob_recall\tNOT_EVALUABLE\n"
        "variant_recall\tNOT_EVALUABLE\n",
        encoding="utf-8",
    )
    aggregate_dir(str(tmp_path / "results"))
    lines = [line.split("\t") for line in capsys.readouterr().out.strip().splitlines()]
    header = lines[0]
    row = next(line for line in lines[1:] if line[0] == "negative_control")
    cells = dict(zip(header, row))
    assert cells["mob_recall"] == "NOT_EVALUABLE"
    assert cells["variant_recall"] == "NOT_EVALUABLE"
    assert cells["mob_precision"] == "0.5000"
    for coerced in ("0", "1", "0.0", "1.0", "0.0000", "1.0000"):
        assert cells["mob_recall"] != coerced
        assert cells["variant_recall"] != coerced


def test_aggregate_results_does_not_invent_recall_from_missing_file(tmp_path: Path, capsys) -> None:
    (tmp_path / "results" / "empty_ds").mkdir(parents=True)
    aggregate_dir(str(tmp_path / "results"))
    lines = [line.split("\t") for line in capsys.readouterr().out.strip().splitlines()]
    header = lines[0]
    row = next(line for line in lines[1:] if line[0] == "empty_ds")
    cells = dict(zip(header, row))
    assert cells["variant_recall"] == "N/A"
    assert cells["mob_recall"] == "N/A"


def _aggregate_cells(capsys, dataset: str) -> dict:
    lines = [line.split("\t") for line in capsys.readouterr().out.strip().splitlines()]
    header = lines[0]
    row = next(line for line in lines[1:] if line[0] == dataset)
    return dict(zip(header, row))


def test_aggregate_aliases_variant_comparator_shaped_keys(tmp_path: Path, capsys) -> None:
    """compare_variant_calls.py --tsv emits MOB_/JC_/SNP_-prefixed keys; their
    NOT_EVALUABLE sentinels must reach the canonical lowercase columns."""
    ds_dir = tmp_path / "results" / "variant_shape"
    ds_dir.mkdir(parents=True)
    (ds_dir / "variant_metrics.tsv").write_text(
        "metric\tvalue\n"
        "variant_precision\tNOT_EVALUABLE\n"
        "variant_recall\tNOT_EVALUABLE\n"
        "variant_f1\tNOT_EVALUABLE\n"
        "SNP_tp\t0\n"
        "SNP_fp\t0\n"
        "SNP_fn\t0\n"
        "SNP_precision\tNOT_EVALUABLE\n"
        "SNP_recall\tNOT_EVALUABLE\n"
        "MOB_tp\t0\n"
        "MOB_fp\t0\n"
        "MOB_fn\t0\n"
        "MOB_precision\tNOT_EVALUABLE\n"
        "MOB_recall\tNOT_EVALUABLE\n"
        "JC_tp\t0\n"
        "JC_fp\t0\n"
        "JC_fn\t0\n"
        "JC_precision\tNOT_EVALUABLE\n"
        "JC_recall\tNOT_EVALUABLE\n"
        "breakpoint_exact_fraction\tNOT_EVALUABLE\n"
        "breakpoint_within_1bp\tNOT_EVALUABLE\n"
        "breakpoint_within_5bp\tNOT_EVALUABLE\n"
        "MC_fp_explained_by_un\t0\n"
        "MC_fp_total\t0\n",
        encoding="utf-8",
    )
    aggregate_dir(str(tmp_path / "results"))
    cells = _aggregate_cells(capsys, "variant_shape")
    for column in (
        "mob_precision",
        "mob_recall",
        "jc_precision",
        "jc_recall",
        "snp_precision",
        "snp_recall",
        "variant_precision",
        "variant_recall",
    ):
        assert cells[column] == "NOT_EVALUABLE", column
        for coerced in ("0", "1", "0.0", "1.0", "0.0000", "1.0000"):
            assert cells[column] != coerced, column


def test_aggregate_aliases_mob_comparator_bare_keys(tmp_path: Path, capsys) -> None:
    """compare_mob_truth.py --tsv emits bare precision/recall plus the marker
    metrics exact_fraction / median_error_bp; the bare keys map onto mob_* only
    and are never promoted to variant_recall."""
    ds_dir = tmp_path / "results" / "mob_shape"
    ds_dir.mkdir(parents=True)
    (ds_dir / "mob_metrics.tsv").write_text(
        "metric\tvalue\n"
        "tp\t0\n"
        "fp\t0\n"
        "fn\t0\n"
        "precision\tNOT_EVALUABLE\n"
        "recall\tNOT_EVALUABLE\n"
        "exact_fraction\tNOT_EVALUABLE\n"
        "within_5bp_fraction\tNOT_EVALUABLE\n"
        "median_error_bp\tNOT_EVALUABLE\n",
        encoding="utf-8",
    )
    aggregate_dir(str(tmp_path / "results"))
    cells = _aggregate_cells(capsys, "mob_shape")
    assert cells["mob_recall"] == "NOT_EVALUABLE"
    assert cells["mob_precision"] == "NOT_EVALUABLE"
    assert cells["variant_recall"] == "N/A"
    assert cells["variant_precision"] == "N/A"


def test_aggregate_merges_conflicting_files_fail_closed(tmp_path: Path, capsys) -> None:
    """A metrics file's NOT_EVALUABLE must survive a later summary.tsv 1.0000;
    identical values pass through and concrete disagreements collapse to
    NOT_EVALUABLE instead of keeping the later file's number."""
    ds_dir = tmp_path / "results" / "conflict"
    ds_dir.mkdir(parents=True)
    (ds_dir / "mob_metrics.tsv").write_text(
        "metric\tvalue\n"
        "mob_recall\tNOT_EVALUABLE\n"
        "mob_precision\t0.5000\n"
        "jc_recall\t0.7500\n"
        "snp_recall\t0.9000\n",
        encoding="utf-8",
    )
    (ds_dir / "summary.tsv").write_text(
        "metric\tvalue\n"
        "mob_recall\t1.0000\n"
        "mob_precision\t0.5000\n"
        "jc_recall\t0.7500\n"
        "snp_recall\t0.8000\n",
        encoding="utf-8",
    )
    aggregate_dir(str(tmp_path / "results"))
    cells = _aggregate_cells(capsys, "conflict")
    assert cells["mob_recall"] == "NOT_EVALUABLE"
    for fabricated in ("1.0000", "1.0", "1", "0.0000"):
        assert cells["mob_recall"] != fabricated
    assert cells["mob_precision"] == "0.5000"
    assert cells["jc_recall"] == "0.7500"
    assert cells["snp_recall"] == "NOT_EVALUABLE"


# ---------------------------------------------------------------------------
# 4. PRJNA884016 negative-control asset hygiene
# ---------------------------------------------------------------------------


def test_manifest_still_repudiates_historical_events() -> None:
    with PRJNA884016_MANIFEST.open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sy22a1 = next(r for r in rows if r["sample_name"] == "SY22A1")
    assert "Negative concordance control" in sy22a1["notes"]
    assert "Historical 4 MOB + 7 JC invalidated" in sy22a1["notes"]

    # The active truth assets carry zero positive structural events.
    gd_text = (TRUTH_DIR / "SY22A1_truth.gd").read_text(encoding="utf-8")
    data_lines = [l for l in gd_text.splitlines() if l.strip() and not l.startswith("#")]
    assert data_lines == []
    for tsv_name in ("mob_truth.tsv", "jc_truth.tsv"):
        lines = [
            l
            for l in (TRUTH_DIR / tsv_name).read_text(encoding="utf-8").splitlines()
            if l.strip() and not l.startswith("#")
        ]
        assert len(lines) == 1, f"{tsv_name} must contain only the header line"


def test_quarantine_manifest_rows_are_all_quarantined_artifacts() -> None:
    with (QUARANTINED_DIR / "MANIFEST.tsv").open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    quarantined_names = {p.name for p in QUARANTINED_DIR.iterdir()}
    for row in rows:
        assert row["quarantine_classification"].startswith("QUARANTINED"), row["file_name"]
        assert row["file_name"] in quarantined_names, row["file_name"]


def test_registry_prjna884016_row_remains_false_and_pending() -> None:
    with REGISTRY_PATH.open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    row = next(r for r in rows if r["dataset"] == "prjna884016")
    for column in (
        "variant_validated",
        "intended_validated",
        "mob_validated",
        "jc_validated",
        "report_reviewed",
    ):
        assert row[column] == "false", column
    assert row["status"] == "PENDING"


def test_active_truth_directory_has_no_stale_mummer_sidecars() -> None:
    for extension in ("delta", "report", "1delta", "mdelta", "1coords", "mcoords", "qdiff", "snps"):
        stale = TRUTH_DIR / f"SY22A1_vs_SY221.{extension}"
        assert not stale.exists(), f"stale MUMmer sidecar must not live in the active truth dir: {stale}"
    # The zeroed rdiff and the negative-control truth assets stay in place.
    assert (TRUTH_DIR / "SY22A1_vs_SY221.rdiff").exists()
    assert (TRUTH_DIR / "SY22A1_truth.gd").exists()
    assert (TRUTH_DIR / "mob_truth.tsv").exists()
    assert (TRUTH_DIR / "jc_truth.tsv").exists()


# ---------------------------------------------------------------------------
# 5. Scope guard
# ---------------------------------------------------------------------------


def test_this_suite_never_references_frozen_m6_assets() -> None:
    """Frozen M6 dataset files are neither imported nor modified by these tests."""
    source = Path(__file__).read_text(encoding="utf-8").lower()
    forbidden = [
        "prjna" + "1088182",
        "m6" + "_p1a",
        "m6" + "_p2a",
        "m6" + "_p3",
        "stage" + "_m6",
    ]
    for marker in forbidden:
        assert marker not in source, f"forbidden M6 reference in this suite: {marker}"
