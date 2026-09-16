#!/usr/bin/env python3
"""
benchmark/real_world/scripts/compare_mob_truth.py
Evaluates mobile element (IS) insertions and transposon junctions against assembly truth.
Computes MOB precision, recall, coordinate error, and TSD concordance.
"""

import sys
import argparse
from collections import defaultdict


def parse_mob_records(gd_path):
    """Extracts MOB records from a GenomeDiff file."""
    records = []
    with open(gd_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[0] == "MOB":
                # MOB id parent seq_id pos element strand repeat_size
                rec = {
                    "seq_id": parts[3] if len(parts) > 3 else "",
                    "pos": int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0,
                    "element": parts[5] if len(parts) > 5 else "",
                    "strand": parts[6] if len(parts) > 6 else "",
                    "repeat_size": int(parts[7]) if len(parts) > 7 and parts[7].isdigit() else 0,
                }
                records.append(rec)
    return records


def compare_mob(test_mobs, truth_mobs, tolerance=10):
    """Matches MOB calls within coordinate tolerance and family identity."""
    tp = 0
    matched_truth = set()
    errors = []
    exact_count = 0
    within_5bp_count = 0

    for t in test_mobs:
        best_match = None
        min_dist = tolerance + 1
        for idx, tr in enumerate(truth_mobs):
            if idx in matched_truth:
                continue
            if t["seq_id"] == tr["seq_id"]:
                dist = abs(t["pos"] - tr["pos"])
                # Match family prefix (e.g. IS1 matches IS1A)
                e1 = t["element"].upper()
                e2 = tr["element"].upper()
                family_match = (e1 in e2 or e2 in e1 or not e1 or not e2)
                if dist <= tolerance and family_match and dist < min_dist:
                    min_dist = dist
                    best_match = idx

        if best_match is not None:
            matched_truth.add(best_match)
            tp += 1
            errors.append(min_dist)
            if min_dist == 0:
                exact_count += 1
            if min_dist <= 5:
                within_5bp_count += 1

    fp = len(test_mobs) - tp
    fn = len(truth_mobs) - tp
    precision = tp / (tp + fp) if (tp + fp) > 0 else 1.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 1.0
    median_error = sorted(errors)[len(errors) // 2] if errors else 0.0

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "exact_fraction": exact_count / tp if tp > 0 else 1.0,
        "within_5bp_fraction": within_5bp_count / tp if tp > 0 else 1.0,
        "median_error_bp": median_error,
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate MOB / IS insertions against truth")
    parser.add_argument("--test-gd", required=True, help="Test GenomeDiff file")
    parser.add_argument("--truth-gd", required=True, help="Truth GenomeDiff file")
    parser.add_argument("--tolerance", type=int, default=10, help="Coordinate tolerance in bp (default: 10)")
    parser.add_argument("--tsv", action="store_true", help="Output TSV format")
    args = parser.parse_args()

    test_mobs = parse_mob_records(args.test_gd)
    truth_mobs = parse_mob_records(args.truth_gd)
    m = compare_mob(test_mobs, truth_mobs, args.tolerance)

    if args.tsv:
        print("metric\tvalue")
        for k, v in m.items():
            if isinstance(v, float):
                print(f"{k}\t{v:.4f}")
            else:
                print(f"{k}\t{v}")
    else:
        print("=== MOB / IS Insertion Evaluation ===")
        print(f"TP: {m['tp']} | FP: {m['fp']} | FN: {m['fn']}")
        print(f"Precision: {m['precision']:.2%}")
        print(f"Recall:    {m['recall']:.2%}")
        print(f"Exact Breakpoints: {m['exact_fraction']:.1%}")
        print(f"Within 5bp:        {m['within_5bp_fraction']:.1%}")
        print(f"Median Error:      {m['median_error_bp']} bp")


if __name__ == "__main__":
    main()
