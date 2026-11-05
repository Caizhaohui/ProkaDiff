#!/usr/bin/env python3
"""Unit tests for jc_acceptance_gate.py (no FASTQ / no Slurm)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from jc_acceptance_gate import parse_jc_gd, run_gate

FIX = Path(__file__).resolve().parent / "_fixtures_jc_gate"


@pytest.fixture(scope="module", autouse=True)
def _fixtures(tmp_path_factory):
    # Build tiny GDs in a temp dir shared via FIX symlink pattern: write beside tests.
    FIX.mkdir(exist_ok=True)
    (FIX / "pd_match.gd").write_text(
        "JC\t1\t.\tchr\t100\t1\tchr\t200\t-1\t0\n"
        "JC\t2\t.\tchr\t50\t1\tchr\t60\t-1\t0\n"
    )
    (FIX / "br_match.gd").write_text(
        "JC\t10\t.\tchr\t100\t1\tchr\t200\t-1\t0\n"
        "JC\t11\t.\tchr\t50\t1\tchr\t60\t-1\t0\treject=COVERAGE_EVENNESS_SKEW\n"
        "JC\t12\t.\tchr\t300\t1\tchr\t400\t-1\t0\n"
    )
    (FIX / "pd_tol.gd").write_text(
        "JC\t1\t.\tchr\t105\t1\tchr\t200\t-1\t0\n"
    )
    (FIX / "br_tol.gd").write_text(
        "JC\t10\t.\tchr\t100\t1\tchr\t202\t-1\t0\n"
    )
    yield


def test_rejected_excluded_from_accepted_truth():
    r = run_gate(FIX / "pd_match.gd", FIX / "br_match.gd", tol=5, allow_repeat_copy_equiv=False)
    assert r["counts"]["breseq_jc_rejected"] == 1
    assert r["counts"]["matched_exact_tol"] == 1  # only 100/200
    # PD id2 matches rejected BR11 → classified separately, not accepted under
    assert r["counts"]["pd_matches_breseq_rejected"] == 1
    assert r["counts"]["accepted_under_unmatched"] == 1  # BR12
    assert r["pass"] is False


def test_tolerance_and_overlap_required():
    r = run_gate(FIX / "pd_tol.gd", FIX / "br_tol.gd", tol=5, allow_repeat_copy_equiv=False)
    assert r["counts"]["matched_exact_tol"] == 1
    assert r["pass"] is True


def test_overlap_mismatch_is_not_match(tmp_path):
    pd = tmp_path / "pd.gd"
    br = tmp_path / "br.gd"
    pd.write_text("JC\t1\t.\tchr\t100\t1\tchr\t200\t-1\t9\n")
    br.write_text("JC\t10\t.\tchr\t100\t1\tchr\t200\t-1\t0\n")
    r = run_gate(pd, br, tol=5, allow_repeat_copy_equiv=False)
    assert r["counts"]["matched_exact_tol"] == 0
    assert r["counts"]["accepted_over_unmatched"] == 1
    assert r["counts"]["accepted_under_unmatched"] == 1


def test_normalize_contig_version_matches_locus_to_version(tmp_path):
    pd = tmp_path / "pd.gd"
    br = tmp_path / "br.gd"
    pd.write_text("JC\t1\t.\tREL606.6\t100\t1\tREL606.6\t200\t-1\t0\n")
    br.write_text("JC\t10\t.\tREL606\t100\t1\tREL606\t200\t-1\t0\n")
    raw = run_gate(pd, br, tol=5, allow_repeat_copy_equiv=False)
    assert raw["counts"]["matched_exact_tol"] == 0
    norm = run_gate(
        pd, br, tol=5, allow_repeat_copy_equiv=False, strip_contig_version=True
    )
    assert norm["counts"]["matched_exact_tol"] == 1
    assert norm["pass"] is True


def test_parse_reject_attr():
    evs = parse_jc_gd(FIX / "br_match.gd", "breseq")
    rej = [e for e in evs if e.rejected]
    assert len(rej) == 1
    assert "COVERAGE_EVENNESS_SKEW" in rej[0].reject_reasons
