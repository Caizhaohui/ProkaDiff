import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plan_m6_p3a import classify_diagnostic

RESULTS = Path(__file__).resolve().parents[1] / "manifests/prjna1088182_m6_p3a_read_level.tsv"


def test_ten_unrecovered_alleles_keep_p2a_status_and_balance():
    rows = list(csv.DictReader(RESULTS.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 10
    keys = {(row["run"], row["coordinate"], row["ref"], row["alt"], row["reference"]) for row in rows}
    assert keys == {
        ("SRR28349204", "1571257", "G", "T", "NC_000913.3"),
        ("SRR28349176", "2776013", "T", "C", "NC_000913.3"),
        ("SRR28349176", "4035824", "T", "C", "NC_000913.3"),
        ("SRR28349230", "2872949", "C", "A", "CP009273.1"),
        ("SRR28349225", "1635887", "A", "G", "CP009273.1"),
        ("SRR28349221", "3108123", "G", "A", "CP009273.1"),
        ("SRR28349165", "161068", "A", "C", "CP009273.1"),
        ("SRR28349232", "2653047", "T", "A", "CP009273.1"),
        ("SRR28349232", "4115617", "C", "A", "CP009273.1"),
        ("SRR28349228", "3073571", "G", "A", "CP009273.1"),
    }
    for row in rows:
        assert row["p2a_recovery"] == "NOT_RECOVERED"
        assert row["p2a_diagnostic"] == "UNKNOWN"
        assert row["p2b_diagnosis"] == "UNRESOLVED"
        assert row["read_level_diagnosis"] == "LOW_ALT_FREQUENCY"
        for tier in ("raw", "caller_shaped"):
            total = int(row[f"{tier}_ref"]) + int(row[f"{tier}_alt"]) + int(row[f"{tier}_other"])
            assert total == int(row[f"{tier}_depth"])
        assert int(row["caller_shaped_depth"]) > 0
        assert float(row["caller_shaped_alt_fraction"]) < 0.75
        assert int(row["plus_caller_shaped_alt"]) > 0
        assert int(row["minus_caller_shaped_alt"]) > 0
        assert "precision" not in row["read_level_diagnosis"].lower()
        again = classify_diagnostic(
            {
                "reference_base": row["ref"],
                "published_ref": row["ref"],
                "caller_shaped_depth": int(row["caller_shaped_depth"]),
                "caller_shaped_alt": int(row["caller_shaped_alt"]),
                "raw_alt": int(row["raw_alt"]),
                "q20_alt": int(row["q20_alt"]),
                "q20_depth": int(row["q20_depth"]),
                "low_mapq_alt": 0,
            }
        )
        assert again == "LOW_ALT_FREQUENCY"


def test_contracts_unchanged():
    repo = Path(__file__).resolve().parents[3]
    diff = subprocess.check_output(
        ["git", "diff", "--stat", "--", "crates", "benchmark/real_world/VALIDATION_REGISTRY.tsv"],
        cwd=repo,
        text=True,
    )
    assert diff.strip() == ""


def test_read_level_tsv_schema_and_numeric_integrity():
    rows = list(csv.DictReader(RESULTS.open(encoding="utf-8"), delimiter="\t"))
    expected_header = [
        "sample",
        "run",
        "reference",
        "coordinate",
        "ref",
        "alt",
        "p2a_recovery",
        "p2a_diagnostic",
        "p2b_diagnosis",
        "raw_depth",
        "raw_ref",
        "raw_alt",
        "raw_other",
        "mapq10_depth",
        "mapq10_alt",
        "caller_shaped_depth",
        "caller_shaped_ref",
        "caller_shaped_alt",
        "caller_shaped_other",
        "caller_shaped_alt_fraction",
        "q20_depth",
        "q20_alt",
        "plus_caller_shaped_alt",
        "minus_caller_shaped_alt",
        "nearby_cigar_reads",
        "duplicate_primary_skipped",
        "read_level_diagnosis",
        "confidence",
        "slurm_array",
    ]
    with RESULTS.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    assert header == expected_header

    int_cols = [
        "raw_depth", "raw_ref", "raw_alt", "raw_other",
        "mapq10_depth", "mapq10_alt",
        "caller_shaped_depth", "caller_shaped_ref", "caller_shaped_alt", "caller_shaped_other",
        "q20_depth", "q20_alt",
        "plus_caller_shaped_alt", "minus_caller_shaped_alt",
        "nearby_cigar_reads", "duplicate_primary_skipped",
    ]
    for row in rows:
        for col in int_cols:
            val = row[col]
            assert val != "", f"Empty string in numeric column {col}"
            assert val.isdigit(), f"Non-integer value {val} in column {col}"
        # Strand conservation
        assert int(row["plus_caller_shaped_alt"]) + int(row["minus_caller_shaped_alt"]) == int(row["caller_shaped_alt"])
        # Alt fraction float precision
        af = float(row["caller_shaped_alt_fraction"])
        expected_af = int(row["caller_shaped_alt"]) / int(row["caller_shaped_depth"])
        assert abs(af - expected_af) < 1e-5


def test_classify_diagnostic_exact_boundaries():
    base_template = {
        "reference_base": "G",
        "published_ref": "G",
        "caller_shaped_depth": 100,
        "caller_shaped_alt": 0,
        "raw_alt": 0,
        "q20_alt": 0,
        "q20_depth": 100,
        "low_mapq_alt": 0,
    }

    # Boundary 1: AF = 0.7499 vs 0.7500
    m_below = dict(base_template, caller_shaped_depth=10000, caller_shaped_alt=7499, raw_alt=7499, q20_alt=7499, q20_depth=10000)
    assert classify_diagnostic(m_below) == "LOW_ALT_FREQUENCY"

    m_at = dict(base_template, caller_shaped_depth=10000, caller_shaped_alt=7500, raw_alt=7500, q20_alt=7500, q20_depth=10000)
    assert classify_diagnostic(m_at) == "ALT_PRESENT_ABOVE_CONSENSUS_THRESHOLD"

    # Boundary 2: Depth = 4 vs 5
    # Depth 4 with raw_alt > 0 and q20_alt > 0 -> UNRESOLVED (usable < 5)
    m_depth4 = dict(base_template, caller_shaped_depth=4, caller_shaped_alt=3, raw_alt=3, q20_alt=3, q20_depth=4)
    assert classify_diagnostic(m_depth4) == "UNRESOLVED"

    # Depth 5 with alt 4 (AF = 0.80 >= 0.75) -> ALT_PRESENT_ABOVE_CONSENSUS_THRESHOLD
    m_depth5_above = dict(base_template, caller_shaped_depth=5, caller_shaped_alt=4, raw_alt=4, q20_alt=4, q20_depth=5)
    assert classify_diagnostic(m_depth5_above) == "ALT_PRESENT_ABOVE_CONSENSUS_THRESHOLD"

    # Depth 5 with alt 3 (AF = 0.60 < 0.75) -> LOW_ALT_FREQUENCY
    m_depth5_below = dict(base_template, caller_shaped_depth=5, caller_shaped_alt=3, raw_alt=3, q20_alt=3, q20_depth=5)
    assert classify_diagnostic(m_depth5_below) == "LOW_ALT_FREQUENCY"

    # Boundary 3: Depth = 19 vs 20 with alt = 0
    m_depth19_zero = dict(base_template, caller_shaped_depth=19, caller_shaped_alt=0, raw_alt=0, q20_alt=0, q20_depth=19)
    assert classify_diagnostic(m_depth19_zero) == "UNRESOLVED"

    m_depth20_zero = dict(base_template, caller_shaped_depth=20, caller_shaped_alt=0, raw_alt=0, q20_alt=0, q20_depth=20)
    assert classify_diagnostic(m_depth20_zero) == "ALT_ABSENT_WITH_ADEQUATE_DEPTH"
