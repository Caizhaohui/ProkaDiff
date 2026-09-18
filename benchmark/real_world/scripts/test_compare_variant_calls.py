#!/usr/bin/env python3
"""
Regression fixture for the RW-003 metric-side fix in compare_variant_calls.py:
- UN records must be parsed but must NOT be scored in match_mutations() (see
  docs/schema.md: "UN 不对拍失败").
- count_mc_fp_explained_by_un() must correctly attribute MC false positives
  that overlap a truth UN region.

Run with: python3 -m pytest benchmark/real_world/scripts/test_compare_variant_calls.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compare_variant_calls import parse_gd, match_mutations, count_mc_fp_explained_by_un


TRUTH_GD = """\
#=GENOME_DIFF	1.0
SNP	1	.	chr	100	A
MC	2	.	chr	1000	1200	200
UN	3	.	chr	5000	5300
"""

# test.gd: one matching SNP, one matching MC, and TWO MC false positives —
# one overlapping the truth UN region (5000-5300), one that overlaps nothing
# (a genuine, unexplained false positive).
TEST_GD = """\
#=GENOME_DIFF	1.0
SNP	1	.	chr	100	A
MC	2	.	chr	1000	1210	210
MC	3	.	chr	5050	5250	200
MC	4	.	chr	9000	9100	100
"""


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content)
    return str(p)


def test_un_is_parsed_but_not_scored(tmp_path):
    truth_path = _write(tmp_path, "truth.gd", TRUTH_GD)
    truth_records = parse_gd(truth_path)

    un_records = [r for r in truth_records if r["type"] == "UN"]
    assert len(un_records) == 1
    assert un_records[0]["start"] == 5000
    assert un_records[0]["end"] == 5300

    # Test set has no UN records at all -> if UN were scored via the generic
    # min(len(tests), len(truths)) branch, it would add 1 phantom FN and drop
    # overall recall. It must not.
    test_records = [r for r in truth_records if r["type"] != "UN"]  # pretend test == truth minus UN
    metrics = match_mutations(test_records, truth_records)
    assert "UN" not in metrics["by_type"]
    # SNP + MC are the only truth types actually scored.
    assert metrics["overall_fn"] == 0


def test_mc_fp_explained_by_un(tmp_path):
    truth_path = _write(tmp_path, "truth.gd", TRUTH_GD)
    test_path = _write(tmp_path, "test.gd", TEST_GD)

    truth_records = parse_gd(truth_path)
    test_records = parse_gd(test_path)

    metrics = match_mutations(test_records, truth_records)
    mc = metrics["by_type"]["MC"]
    assert mc["tp"] == 1
    assert mc["fp"] == 2
    assert mc["fn"] == 0

    test_mc = [r for r in test_records if r["type"] == "MC"]
    truth_mc = [r for r in truth_records if r["type"] == "MC"]
    truth_un = [r for r in truth_records if r["type"] == "UN"]

    explained, total = count_mc_fp_explained_by_un(test_mc, truth_mc, truth_un)
    assert total == 2
    # Exactly one of the two FPs (chr:5050-5250) overlaps the UN region
    # (chr:5000-5300); the other (chr:9000-9100) does not.
    assert explained == 1


def test_rejected_records_skipped_by_default(tmp_path):
    # A truth file containing one accepted JC and one rejected candidate
    # (e.g. breseq reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF)
    gd_content = """\
#=GENOME_DIFF\t1.0
JC\t1\t.\tchr\t100\t1\tchr\t200\t-1\t0
JC\t2\t.\tchr\t300\t1\tchr\t400\t-1\t0\treject=COVERAGE_EVENNESS_SKEW\tfrequency=0.01
"""
    p = _write(tmp_path, "sample.gd", gd_content)

    # By default, reject= records must be skipped
    default_records = parse_gd(p)
    assert len(default_records) == 1
    assert default_records[0]["id"] == "1"

    # With include_rejected=True, all records are returned
    all_records = parse_gd(p, include_rejected=True)
    assert len(all_records) == 2
    assert all_records[1]["id"] == "2"


if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([__file__, "-v"]))
