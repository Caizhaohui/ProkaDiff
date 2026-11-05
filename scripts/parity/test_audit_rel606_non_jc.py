#!/usr/bin/env python3
"""Regression tests for audit_rel606_non_jc.py (strict identity contract).

Covers FIX-M8-T03-AUDIT-SCRIPT-RIGOR (Codex Sol review): tampering with an
INS new_seq, DEL size, or MOB strand / duplication_size / repeat_name (and
SNP new_seq, position, or 1-to-1 multiplicity) must break the match and make
the audit exit 1. No FASTQ / no Slurm / no real benchmark data required.

Also covers FIX-M8-T03-RECORD-ACCOUNTING-GAP (Codex Sol re-review): record
accounting must be by positional index, never by GD id, so a duplicate-id
ProkaDiff record that does not match breseq still counts as OVER (and the
breseq mirror as UNDER), and parse_gd() warns about duplicate record ids.

Also covers FIX-M8-T03-MAXIMUM-BIPARTITE-MATCHING (Codex Sol review):
match_kind() must be a true maximum-cardinality bipartite matching; greedy
closest-first pairing under-counted matches when tol > 0.

Also covers FIX-M8-T03-ITERATIVE-AUGMENTING-PATH (Codex Sol re-review):
the augmenting-path search must be iterative (BFS + deque), never
recursive — the recursive DFS hit CPython's recursion limit on a chain of
1050 SNPs with tol=1.
"""
from __future__ import annotations

import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from audit_rel606_non_jc import candidate_ok, main, match_kind, parse_gd

SCRIPT = Path(__file__).resolve().parent / "audit_rel606_non_jc.py"

TSV_HEADER = (
    "mutation_type", "contig", "pd_pos", "breseq_pos", "pd_allele_or_elem",
    "breseq_allele_or_elem", "delta_bp", "status", "notes",
)

# One consensus record per redline kind, mirroring the REL606 clonal GDs
# (trailing evidence annotations are ignored by the parser).
BASE_PD = (
    "SNP\t1\t.\tREL606\t161041\tG\tpd_ra_depth=58\n"
    "INS\t2\t.\tREL606\t475292\tG\tpd_ra_depth=37\n"
    "DEL\t3\t.\tREL606\t4126706\t1\tpd_ra_depth=18\n"
    "MOB\t4\t10,11\tREL606\t4524522\tIS186\t1\t6\tpd_note=x\n"
)
BASE_BR = (
    "SNP\t101\t43\tREL606\t161041\tG\n"
    "INS\t102\t48\tREL606\t475292\tG\n"
    "DEL\t103\t82\tREL606\t4126706\t1\trepeat_length=1\n"
    "MOB\t104\t106,107\tREL606\t4524522\tIS186\t1\t6\n"
)


def run_audit(monkeypatch, tmp_path, pd_text, br_text, extra=()):
    """Write fixture GDs, run main(), return (exit_code, tsv_rows)."""
    pd_gd = tmp_path / "pd.gd"
    br_gd = tmp_path / "br.gd"
    out = tmp_path / "out.tsv"
    pd_gd.write_text(pd_text)
    br_gd.write_text(br_text)
    argv = [
        "audit_rel606_non_jc.py",
        "--pd-gd", str(pd_gd),
        "--breseq-gd", str(br_gd),
        "--out", str(out),
        "--fallback-gd", str(tmp_path / "missing_fallback.gd"),
        *extra,
    ]
    monkeypatch.setattr(sys, "argv", argv)
    rc = main()
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == list(TSV_HEADER)
    rows = [dict(zip(TSV_HEADER, line.split("\t"))) for line in lines[1:]]
    return rc, rows


def statuses(rows, kind):
    return Counter(r["status"] for r in rows if r["mutation_type"] == kind)


def test_exact_identity_passes(monkeypatch, tmp_path):
    rc, rows = run_audit(monkeypatch, tmp_path, BASE_PD, BASE_BR)
    assert rc == 0
    assert statuses(rows, "SNP") == {"MATCHED": 1}
    assert statuses(rows, "INS") == {"MATCHED": 1}
    assert statuses(rows, "DEL") == {"MATCHED": 1}
    assert statuses(rows, "MOB") == {"MATCHED": 1}
    assert all(r["delta_bp"] == "0" for r in rows)


def test_contig_version_suffix_normalized(monkeypatch, tmp_path):
    pd = BASE_PD.replace("REL606", "REL606.6")
    rc, rows = run_audit(monkeypatch, tmp_path, pd, BASE_BR)
    assert rc == 0
    assert len(rows) == 4 and all(r["status"] == "MATCHED" for r in rows)


def test_snp_allele_tamper_fails(monkeypatch, tmp_path):
    rc, rows = run_audit(
        monkeypatch, tmp_path,
        BASE_PD.replace("161041\tG", "161041\tT"), BASE_BR,
    )
    assert rc == 1
    assert statuses(rows, "SNP") == {"OVER": 1, "UNDER": 1}


def test_ins_allele_tamper_fails(monkeypatch, tmp_path):
    rc, rows = run_audit(
        monkeypatch, tmp_path,
        BASE_PD.replace("475292\tG", "475292\tA"), BASE_BR,
    )
    assert rc == 1
    assert statuses(rows, "INS") == {"OVER": 1, "UNDER": 1}
    over = next(r for r in rows if r["status"] == "OVER")
    assert "allele pd=A breseq=G" in over["notes"]


def test_del_size_tamper_fails(monkeypatch, tmp_path):
    rc, rows = run_audit(
        monkeypatch, tmp_path,
        BASE_PD.replace("4126706\t1\t", "4126706\t2\t"), BASE_BR,
    )
    assert rc == 1
    assert statuses(rows, "DEL") == {"OVER": 1, "UNDER": 1}
    over = next(r for r in rows if r["status"] == "OVER")
    assert "size pd=2 breseq=1" in over["notes"]


def test_mob_strand_tamper_fails(monkeypatch, tmp_path):
    rc, rows = run_audit(
        monkeypatch, tmp_path,
        BASE_PD.replace("IS186\t1\t6", "IS186\t-1\t6"), BASE_BR,
    )
    assert rc == 1
    assert statuses(rows, "MOB") == {"OVER": 1, "UNDER": 1}
    over = next(r for r in rows if r["status"] == "OVER")
    assert "strand pd=-1 breseq=1" in over["notes"]


def test_mob_dup_size_tamper_fails(monkeypatch, tmp_path):
    rc, rows = run_audit(
        monkeypatch, tmp_path,
        BASE_PD.replace("IS186\t1\t6", "IS186\t1\t5"), BASE_BR,
    )
    assert rc == 1
    assert statuses(rows, "MOB") == {"OVER": 1, "UNDER": 1}
    over = next(r for r in rows if r["status"] == "OVER")
    assert "dup_size pd=5 breseq=6" in over["notes"]


def test_mob_repeat_name_tamper_fails(monkeypatch, tmp_path):
    rc, rows = run_audit(
        monkeypatch, tmp_path,
        BASE_PD.replace("\tIS186\t", "\tIS150\t"), BASE_BR,
    )
    assert rc == 1
    assert statuses(rows, "MOB") == {"OVER": 1, "UNDER": 1}
    over = next(r for r in rows if r["status"] == "OVER")
    assert "repeat_name pd=IS150 breseq=IS186" in over["notes"]


def test_position_shift_not_matched_at_default_tol(monkeypatch, tmp_path):
    rc, rows = run_audit(
        monkeypatch, tmp_path,
        BASE_PD.replace("475292", "475297"), BASE_BR,
    )
    assert rc == 1
    assert statuses(rows, "INS") == {"OVER": 1, "UNDER": 1}


def test_nonzero_tol_is_diagnostic_only_and_never_relaxes_attributes(
    monkeypatch, tmp_path, capsys
):
    # Position-only shift within the diagnostic tolerance may match...
    rc, rows = run_audit(
        monkeypatch, tmp_path,
        BASE_PD.replace("475292", "475297"), BASE_BR, extra=("--tol-bp", "5"),
    )
    assert rc == 0
    assert statuses(rows, "INS") == {"MATCHED": 1}
    assert "diagnostic override" in capsys.readouterr().err
    # ...but an attribute difference never matches, whatever tol is.
    rc, rows = run_audit(
        monkeypatch, tmp_path,
        BASE_PD.replace("IS186\t1\t6", "IS186\t1\t5"), BASE_BR,
        extra=("--tol-bp", "10"),
    )
    assert rc == 1
    assert statuses(rows, "MOB") == {"OVER": 1, "UNDER": 1}


def test_matching_is_strictly_one_to_one(monkeypatch, tmp_path):
    # Two identical ProkaDiff INS identities vs one breseq INS: only one may
    # be consumed; the duplicate must stay OVER.
    pd = BASE_PD + "INS\t9\t.\tREL606\t475292\tG\n"
    rc, rows = run_audit(monkeypatch, tmp_path, pd, BASE_BR)
    assert rc == 1
    assert statuses(rows, "INS") == {"MATCHED": 1, "OVER": 1}


def test_maximum_matching_greedy_counterexample(monkeypatch, tmp_path):
    """FIX-M8-T03-MAXIMUM-BIPARTITE-MATCHING: Codex Sol's counterexample.

    SNPs at ProkaDiff positions 2, 3 vs breseq positions 1, 2 with tol=1.
    Greedy closest-first selection took (2, 2) [delta 0], leaving pd 3
    unmatched (breseq 1 is 2bp away) and reporting matched=1. True
    maximum bipartite matching is {(2, 1), (3, 2)}: matched=2, over=0,
    under=0.
    """
    pd = (
        "SNP\t1\t.\tREL606\t2\tA\n"
        "SNP\t2\t.\tREL606\t3\tA\n"
    )
    br = (
        "SNP\t101\t43\tREL606\t1\tA\n"
        "SNP\t102\t43\tREL606\t2\tA\n"
    )
    pd_gd = tmp_path / "pd_counter.gd"
    br_gd = tmp_path / "br_counter.gd"
    pd_gd.write_text(pd)
    br_gd.write_text(br)
    pd_records, _ = parse_gd(pd_gd)
    br_records, _ = parse_gd(br_gd)
    res = match_kind(pd_records["SNP"], br_records["SNP"], tol=1)
    assert len(res["matched"]) == 2
    assert res["over"] == []
    assert res["under"] == []
    # End-to-end with --tol-bp 1: both SNPs MATCHED, exit 0.
    rc, rows = run_audit(monkeypatch, tmp_path, pd, br, extra=("--tol-bp", "1"))
    assert rc == 0
    assert statuses(rows, "SNP") == {"MATCHED": 2}


def test_stress_1050_snps_no_recursion_error(tmp_path):
    """FIX-M8-T03-ITERATIVE-AUGMENTING-PATH: Codex Sol's exact stress test.

    n=1050 SNPs: ProkaDiff positions 2*i, breseq positions 2*i+1, tol=1.
    Every pd SNP edges only to its two neighbouring br SNPs, so each
    augmenting path for pd_i must chain through pd_{i-1}, pd_{i-2}, ...
    The recursive DFS try_augment() followed that chain on the call stack
    (depth ~1050 > CPython's default limit 1000) and raised
    RecursionError. The iterative BFS augmenting path must complete with
    matched=1050, over=0, under=0.

    The tightened-recursion-limit rerun proves the search is genuinely
    iterative: even with sys.setrecursionlimit(150) it still returns the
    full maximum matching (the old recursive version died at ~150 frames
    already), while pytest's own frame depth stays well under 150.
    """
    n = 1050
    pd_lines = [f"SNP\t{i}\t.\tREL606\t{2 * i}\tA" for i in range(n)]
    br_lines = [f"SNP\t{1000 + i}\t43\tREL606\t{2 * i + 1}\tA" for i in range(n)]
    pd_gd = tmp_path / "pd_stress.gd"
    br_gd = tmp_path / "br_stress.gd"
    pd_gd.write_text("\n".join(pd_lines) + "\n")
    br_gd.write_text("\n".join(br_lines) + "\n")
    pd_records, _ = parse_gd(pd_gd)
    br_records, _ = parse_gd(br_gd)
    assert len(pd_records["SNP"]) == n
    assert len(br_records["SNP"]) == n

    res = match_kind(pd_records["SNP"], br_records["SNP"], tol=1)
    assert len(res["matched"]) == 1050
    assert res["over"] == []
    assert res["under"] == []

    # Same computation under a deliberately tiny recursion limit: an
    # iterative implementation is unaffected, a recursive one cannot run.
    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(150)
    try:
        res = match_kind(pd_records["SNP"], br_records["SNP"], tol=1)
    finally:
        sys.setrecursionlimit(old_limit)
    assert len(res["matched"]) == 1050
    assert res["over"] == []
    assert res["under"] == []


def test_duplicate_pd_id_unmatched_record_still_accounted(
    monkeypatch, tmp_path, capsys
):
    """FIX-M8-T03-RECORD-ACCOUNTING-GAP: Codex Sol's counterexample.

    Two ProkaDiff SNP records sharing GD id '1'; only the first matches the
    single breseq SNP. Id-keyed accounting used to consume BOTH pd records
    with one id and report pd=2, breseq=1, matched=1, over=0 -> exit 0.
    Index-keyed accounting must report pd=2, breseq=1, matched=1, over=1
    and exit 1, and parse_gd() must warn about the duplicate id.
    """
    pd = (
        "SNP\t1\t.\tREL606\t161041\tG\n"
        "SNP\t1\t.\tREL606\t999999\tA\n"
    )
    br = "SNP\t101\t43\tREL606\t161041\tG\n"
    rc, rows = run_audit(monkeypatch, tmp_path, pd, br)
    assert rc == 1
    assert statuses(rows, "SNP") == {"MATCHED": 1, "OVER": 1}
    over_row = next(r for r in rows if r["status"] == "OVER")
    assert over_row["pd_pos"] == "999999"
    captured = capsys.readouterr()
    assert "Per-type SNP: pd=2 breseq=1 matched=1" in captured.out
    assert "duplicate REDLINE record ids" in captured.err
    assert "'1' x2" in captured.err


def test_two_pd_records_distinct_ids_one_match_still_accounted(
    monkeypatch, tmp_path
):
    """Distinct ids, same shape: exactly one match -> over=1, exit 1.

    Locks the strict accounting invariant for the non-duplicate variant of
    the counterexample: every pd record appears in exactly one status row.
    """
    pd = (
        "SNP\t1\t.\tREL606\t161041\tG\n"
        "SNP\t2\t.\tREL606\t999999\tA\n"
    )
    br = "SNP\t101\t43\tREL606\t161041\tG\n"
    rc, rows = run_audit(monkeypatch, tmp_path, pd, br)
    assert rc == 1
    assert statuses(rows, "SNP") == {"MATCHED": 1, "OVER": 1}
    over_row = next(r for r in rows if r["status"] == "OVER")
    assert over_row["pd_pos"] == "999999"


def test_duplicate_breseq_id_unmatched_record_still_accounted(
    monkeypatch, tmp_path
):
    """Mirror of the accounting gap on the breseq (UNDER) side.

    Two breseq SNP records sharing id '101'; only one matches. The
    unmatched duplicate must surface as UNDER (breseq=2, pd=1, matched=1,
    under=1 -> exit 1), not silently vanish via id-keyed consumption.
    """
    pd = "SNP\t1\t.\tREL606\t161041\tG\n"
    br = (
        "SNP\t101\t43\tREL606\t161041\tG\n"
        "SNP\t101\t43\tREL606\t888888\tT\n"
    )
    rc, rows = run_audit(monkeypatch, tmp_path, pd, br)
    assert rc == 1
    assert statuses(rows, "SNP") == {"MATCHED": 1, "UNDER": 1}
    under_row = next(r for r in rows if r["status"] == "UNDER")
    assert under_row["breseq_pos"] == "888888"


def test_parse_gd_warns_on_duplicate_ids_no_warning_when_unique(tmp_path, capsys):
    pd_gd = tmp_path / "pd.gd"
    pd_gd.write_text(BASE_PD)  # ids 1..4, unique
    parse_gd(pd_gd)
    assert "duplicate REDLINE record ids" not in capsys.readouterr().err
    dup_gd = tmp_path / "dup.gd"
    dup_gd.write_text(BASE_PD + "SNP\t4\t.\tREL606\t200001\tT\n")  # id '4' twice
    parse_gd(dup_gd)
    err = capsys.readouterr().err
    assert "duplicate REDLINE record ids" in err
    assert "'4' x2" in err


def test_candidate_ok_rejects_attribute_diff_even_with_tolerance(tmp_path):
    pd_gd = tmp_path / "pd.gd"
    br_gd = tmp_path / "br.gd"
    pd_gd.write_text(BASE_PD)
    br_gd.write_text(BASE_BR)
    pd_records, _ = parse_gd(pd_gd)
    br_records, _ = parse_gd(br_gd)
    strict = candidate_ok(pd_records["MOB"][0], br_records["MOB"][0], tol=0)
    assert strict == (True, 0)
    tampered = dict(pd_records["MOB"][0], dup=5)
    assert candidate_ok(tampered, br_records["MOB"][0], tol=100) == (False, 0)
    shifted = dict(pd_records["INS"][0], pos=475297)
    assert candidate_ok(shifted, br_records["INS"][0], tol=0) == (False, 5)
    assert candidate_ok(shifted, br_records["INS"][0], tol=5) == (True, 5)


def test_cli_exit_codes_end_to_end(tmp_path):
    """The real process must exit 0 on identity and 1 on a tampered allele."""
    pd_gd = tmp_path / "pd.gd"
    br_gd = tmp_path / "br.gd"
    out = tmp_path / "out.tsv"
    pd_gd.write_text(BASE_PD)
    br_gd.write_text(BASE_BR)
    common = [
        sys.executable, str(SCRIPT),
        "--pd-gd", str(pd_gd), "--breseq-gd", str(br_gd),
        "--out", str(out),
        "--fallback-gd", str(tmp_path / "missing.gd"),
    ]
    assert subprocess.run(common, capture_output=True, text=True).returncode == 0
    pd_gd.write_text(BASE_PD.replace("475292\tG", "475292\tA"))
    proc = subprocess.run(common, capture_output=True, text=True)
    assert proc.returncode == 1
    assert "Result: FAIL" in proc.stdout
    # Missing input file exits 2 (input error), distinct from parity FAIL.
    proc = subprocess.run(
        common[:4] + [str(tmp_path / "nope.gd")] + common[5:],
        capture_output=True, text=True,
    )
    assert proc.returncode == 2
