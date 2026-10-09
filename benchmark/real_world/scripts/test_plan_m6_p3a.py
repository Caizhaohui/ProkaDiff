import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plan_m6_p3a import OUT, build_targets, classify_diagnostic


def _measure(**overrides):
    row = {
        "reference_base": "G",
        "published_ref": "G",
        "caller_shaped_depth": 40,
        "caller_shaped_alt": 0,
        "raw_alt": 0,
        "q20_alt": 0,
        "q20_depth": 40,
        "low_mapq_alt": 0,
    }
    row.update(overrides)
    return row


def test_ten_unrecovered_alleles_stay_unmeasured():
    rows = build_targets()
    assert len(rows) == 10
    assert len({(row["sample"], row["run"], row["coordinate"], row["ref"], row["alt"]) for row in rows}) == 10
    assert all(row["original_recovery_status"] == "NOT_RECOVERED" for row in rows)
    assert all(row["p2a_diagnostic"] == "UNKNOWN" for row in rows)
    assert all(row["p2b_diagnosis"] == "UNRESOLVED" for row in rows)
    assert all(row["bam_artifact"] == "ABSENT" for row in rows)
    assert all(row["depth_status"] == "NOT_MEASURED" for row in rows)
    assert all(row["alt"] != "+AAT" for row in rows)


def test_priority_subset_matches_authoritative_rows():
    rows = build_targets()
    by_run = {}
    for row in rows:
        by_run.setdefault(row["run"], []).append(row)
    first = by_run["SRR28349204"]
    assert len(first) == 1
    assert first[0]["sample"] == "Cas9_3 _and2others"
    assert first[0]["reference"] == "NC_000913.3"
    assert (first[0]["coordinate"], first[0]["ref"], first[0]["alt"]) == ("1571257", "G", "T")
    second = by_run["SRR28349176"]
    assert {(row["coordinate"], row["ref"], row["alt"]) for row in second} == {
        ("2776013", "T", "C"),
        ("4035824", "T", "C"),
    }
    assert {row["reference"] for row in second} == {"NC_000913.3"}
    third = by_run["SRR28349230"]
    assert len(third) == 1
    assert third[0]["sample"] == "I-SceI_7_8_9"
    assert third[0]["reference"] == "CP009273.1"
    assert (third[0]["coordinate"], third[0]["ref"], third[0]["alt"]) == ("2872949", "C", "A")
    assert all(row["diagnostic_priority"] == "1" for row in first + second + third)
    assert sum(row["diagnostic_priority"] == "2" for row in rows) == 6


def test_classification_uses_counts_not_mixture_percentages():
    assert classify_diagnostic(_measure()) == "ALT_ABSENT_WITH_ADEQUATE_DEPTH"
    assert classify_diagnostic(_measure(caller_shaped_depth=4, q20_depth=4)) == "INSUFFICIENT_COVERAGE"
    assert classify_diagnostic(_measure(caller_shaped_alt=8, raw_alt=8, q20_alt=8)) == "LOW_ALT_FREQUENCY"
    assert classify_diagnostic(_measure(caller_shaped_alt=32, raw_alt=32, q20_alt=32)) == "ALT_PRESENT_ABOVE_CONSENSUS_THRESHOLD"
    assert classify_diagnostic(_measure(raw_alt=6, q20_alt=0, caller_shaped_alt=0)) == "LOW_QUALITY_SUPPORT"
    assert classify_diagnostic(_measure(reference_base="A")) == "REFERENCE_DISAGREEMENT"
    assert classify_diagnostic(_measure(caller_shaped_depth=10, caller_shaped_alt=0, q20_depth=10)) == "UNRESOLVED"


def test_versioned_manifest_matches_builder_and_contracts_are_unchanged():
    stored = list(csv.DictReader(OUT.open(encoding="utf-8"), delimiter="\t"))
    assert stored == build_targets()
    repo = Path(__file__).resolve().parents[3]
    diff = subprocess.check_output(
        ["git", "diff", "--stat", "--", "crates", "benchmark/real_world/VALIDATION_REGISTRY.tsv"],
        cwd=repo,
        text=True,
    )
    assert diff.strip() == ""
