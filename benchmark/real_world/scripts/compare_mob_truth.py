#!/usr/bin/env python3
"""
benchmark/real_world/scripts/compare_mob_truth.py
Evaluates mobile element (IS) insertions and transposon junctions against assembly truth.
Computes MOB precision, recall, coordinate error, and TSD concordance.
"""

import sys
import argparse
from collections import defaultdict
from pathlib import Path

_QUARANTINED_TRUTH_DIR = (
    Path(__file__).resolve().parent.parent / "truth" / "prjna884016" / "quarantined_legacy"
)


def reject_quarantined_truth(truth_gd):
    """Fail closed (exit 2): quarantined PRJNA884016 legacy truth must never be
    loaded as positive truth (M7-P0b). Applies to --truth-gd only; --test-gd is
    never rejected by this check. Ordinary positive truth (e.g. BL21 curated)
    passes through untouched."""
    resolved = Path(truth_gd).resolve()
    quarantined = "invalid_legacy" in resolved.name or _QUARANTINED_TRUTH_DIR.resolve() in resolved.parents
    if not quarantined:
        return
    print(
        f"ERROR: quarantined truth cannot be used as positive truth: {truth_gd} "
        "(PRJNA884016 quarantined_legacy / invalid_legacy, M7-P0b)",
        file=sys.stderr,
    )
    sys.exit(2)


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

    if (tp + fp) > 0:
        precision = tp / (tp + fp)
    else:
        precision = "NOT_EVALUABLE"

    if (tp + fn) > 0:
        recall = tp / (tp + fn)
    else:
        # Zero verified positive events in truth: recall is NOT_EVALUABLE (§6, M7-P0a)
        recall = "NOT_EVALUABLE"

    if tp > 0:
        exact_fraction = exact_count / tp
        within_5bp_fraction = within_5bp_count / tp
        median_error = sorted(errors)[len(errors) // 2] if errors else 0.0
    else:
        exact_fraction = "NOT_EVALUABLE"
        within_5bp_fraction = "NOT_EVALUABLE"
        median_error = "NOT_EVALUABLE"

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "exact_fraction": exact_fraction,
        "within_5bp_fraction": within_5bp_fraction,
        "median_error_bp": median_error,
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate MOB / IS insertions against truth")
    parser.add_argument("--test-gd", required=True, help="Test GenomeDiff file")
    parser.add_argument("--truth-gd", required=True, help="Truth GenomeDiff file")
    parser.add_argument("--tolerance", type=int, default=10, help="Coordinate tolerance in bp (default: 10)")
    parser.add_argument("--tsv", action="store_true", help="Output TSV format")
    args = parser.parse_args()

    reject_quarantined_truth(args.truth_gd)

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
        p_str = f"{m['precision']:.2%}" if isinstance(m['precision'], (int, float)) else str(m['precision'])
        r_str = f"{m['recall']:.2%}" if isinstance(m['recall'], (int, float)) else str(m['recall'])
        ef_str = f"{m['exact_fraction']:.1%}" if isinstance(m['exact_fraction'], (int, float)) else str(m['exact_fraction'])
        w5_str = f"{m['within_5bp_fraction']:.1%}" if isinstance(m['within_5bp_fraction'], (int, float)) else str(m['within_5bp_fraction'])
        me_str = f"{m['median_error_bp']} bp" if isinstance(m['median_error_bp'], (int, float)) else str(m['median_error_bp'])
        print(f"Precision: {p_str}")
        print(f"Recall:    {r_str}")
        print(f"Exact Breakpoints: {ef_str}")
        print(f"Within 5bp:        {w5_str}")
        print(f"Median Error:      {me_str}")


if __name__ == "__main__":
    main()
