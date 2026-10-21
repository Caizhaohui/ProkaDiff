import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from measure_m6_p3a_locus import analyze_sam


def _sam(body: str) -> list[str]:
    header = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:NC_000913.3\tLN:100"]
    return header + [line for line in body.strip().splitlines() if line]


def test_primary_counts_partition_depth_and_skip_secondary():
    lines = _sam(
        "\n".join(
            [
                "r1\t0\tNC_000913.3\t10\t40\t5M\t*\t0\t0\tGGGGT\tIIIII\tAS:i:5",
                "r1\t256\tNC_000913.3\t10\t0\t5M\t*\t0\t0\tTTTTA\tIIIII",
                "r2\t16\tNC_000913.3\t10\t40\t5M\t*\t0\t0\tAGGGA\tIIIII",
                "r3\t0\tNC_000913.3\t10\t40\t5M\t*\t0\t0\tGGCGG\t#####",
            ]
        )
    )
    # Position 12 is the middle base of a 5M alignment starting at 10: G,G,C for r1,r2,r3.
    counts = analyze_sam(lines, 12, "G", "T")
    assert counts["duplicate_primary_skipped"] == 0
    assert counts["raw_depth"] == 3
    assert counts["raw_ref"] + counts["raw_alt"] + counts["raw_other"] == counts["raw_depth"]
    assert counts["caller_shaped_ref"] + counts["caller_shaped_alt"] + counts["caller_shaped_other"] == counts["caller_shaped_depth"]
    assert counts["mapq10_ref"] + counts["mapq10_alt"] + counts["mapq10_other"] == counts["mapq10_depth"]
    assert counts["q20_ref"] + counts["q20_alt"] + counts["q20_other"] == counts["q20_depth"]
    assert counts["plus_caller_shaped"] + counts["minus_caller_shaped"] == counts["caller_shaped_depth"]
    assert counts["raw_ref"] == 2
    assert counts["raw_other"] == 1
    assert counts["minus_caller_shaped"] == 1
    assert counts["plus_caller_shaped"] >= 1


def test_end_trim_and_low_quality_are_separate_filters():
    lines = _sam("edge\t0\tNC_000913.3\t10\t40\t3M\t*\t0\t0\tTGG\tI#I")
    counts = analyze_sam(lines, 10, "G", "T")
    assert counts["raw_alt"] == 1
    assert counts["caller_shaped_depth"] == 0
    assert counts["q20_alt"] == 1
    noisy = _sam("low\t0\tNC_000913.3\t10\t40\t5M\t*\t0\t0\tGGTTT\tII#II")
    low = analyze_sam(noisy, 12, "G", "T")
    assert low["raw_alt"] == 1
    assert low["q20_alt"] == 0
    assert low["q20_depth"] == 0


def test_xs_score_removes_mapq_and_duplicate_qname_is_not_counted_twice():
    lines = _sam(
        "\n".join(
            [
                "amb\t0\tNC_000913.3\t10\t40\t5M\t*\t0\t0\tGGTTT\tIIIII\tAS:i:4\tXS:i:4",
                "amb\t0\tNC_000913.3\t10\t40\t5M\t*\t0\t0\tGGTTT\tIIIII",
            ]
        )
    )
    counts = analyze_sam(lines, 12, "G", "T")
    assert counts["raw_depth"] == 1
    assert counts["duplicate_primary_skipped"] == 1
    assert counts["mapq10_depth"] == 0
    assert counts["low_mapq_alt"] == 1
    assert counts["raw_ref"] + counts["raw_alt"] + counts["raw_other"] == counts["raw_depth"]


def test_paired_mates_with_one_qname_are_both_counted():
    lines = _sam(
        "\n".join(
            [
                "pair\t99\tNC_000913.3\t10\t40\t5M\t=\t10\t5\tGGGGT\tIIIII",
                "pair\t147\tNC_000913.3\t10\t40\t5M\t=\t10\t-5\tGGAGA\tIIIII",
            ]
        )
    )
    counts = analyze_sam(lines, 12, "G", "T")
    assert counts["raw_depth"] == 2
    assert counts["duplicate_primary_skipped"] == 0
    assert counts["plus_caller_shaped"] == 1
    assert counts["minus_caller_shaped"] == 1
    assert counts["raw_ref"] == 1
    assert counts["raw_other"] == 1


def test_reverse_flag_is_strand_not_a_second_read():
    lines = _sam("rev\t16\tNC_000913.3\t10\t30\t5M\t*\t0\t0\tGGGGT\tIIIII")
    counts = analyze_sam(lines, 12, "G", "T")
    assert counts["raw_depth"] == 1
    assert counts["minus_caller_shaped"] == 1
    assert counts["plus_caller_shaped"] == 0
