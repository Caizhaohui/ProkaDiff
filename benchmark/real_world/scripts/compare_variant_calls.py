#!/usr/bin/env python3
"""
benchmark/real_world/scripts/compare_variant_calls.py
Compares ProkaDiff variant calls against oracle (breseq) or curated/assembly truth.
Computes precision, recall, F1, and breakpoint tolerance metrics.
"""

import sys
import argparse
import re
from collections import defaultdict
from pathlib import Path

NOT_EVALUABLE = "NOT_EVALUABLE"

_QUARANTINED_TRUTH_DIR = (
    Path(__file__).resolve().parent.parent / "truth" / "prjna884016" / "quarantined_legacy"
)


def safe_ratio(numerator, denominator):
    """Zero-denominator metrics are NOT_EVALUABLE, never a fabricated 0.0/1.0 (M7-P0b)."""
    if denominator == 0:
        return NOT_EVALUABLE
    return numerator / denominator


def format_metric(value, spec):
    """Print the NOT_EVALUABLE sentinel verbatim; never float-format it (M7-P0b)."""
    if isinstance(value, str):
        return value
    return format(value, spec)


def reject_quarantined_truth(truth_gd):
    """Fail closed (exit 2): quarantined PRJNA884016 legacy truth must never be
    scored as positive truth (M7-P0b). Applies to --truth-gd only; --test-gd is
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


def parse_gd(filepath, include_rejected=False):
    """Parses a GenomeDiff (.gd) file into a list of mutation dicts."""
    records = []
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if not include_rejected and any(p.startswith("reject=") for p in parts):
                continue
            record_type = parts[0]
            if record_type in ("SNP", "SUB", "INS", "DEL", "MOB", "AMP", "CON", "INV"):
                # Standard mutation: TYPE id parent seq_id pos ...
                rec = {
                    "type": record_type,
                    "id": parts[1] if len(parts) > 1 else "",
                    "parent_ids": parts[2] if len(parts) > 2 else "",
                    "seq_id": parts[3].split(".")[0] if len(parts) > 3 else "",
                    "pos": int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0,
                    "raw": parts,
                }
                if record_type == "SNP":
                    rec["new_seq"] = parts[5] if len(parts) > 5 else ""
                elif record_type == "INS":
                    rec["new_seq"] = parts[5] if len(parts) > 5 else ""
                elif record_type == "DEL":
                    rec["size"] = int(parts[5]) if len(parts) > 5 and parts[5].isdigit() else 0
                elif record_type == "MOB":
                    rec["element"] = parts[5] if len(parts) > 5 else ""
                    rec["strand"] = parts[6] if len(parts) > 6 else ""
                    rec["repeat_size"] = int(parts[7]) if len(parts) > 7 and parts[7].isdigit() else 0
                records.append(rec)
            elif record_type == "JC":
                # JC id parent seq1 pos1 strand1 seq2 pos2 strand2 overlap
                rec = {
                    "type": "JC",
                    "id": parts[1] if len(parts) > 1 else "",
                    "parent_ids": parts[2] if len(parts) > 2 else "",
                    "seq_id_1": parts[3].split(".")[0] if len(parts) > 3 else "",
                    "pos_1": int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0,
                    "strand_1": int(parts[5]) if len(parts) > 5 and parts[5] in ("1", "-1") else 1,
                    "seq_id_2": parts[6].split(".")[0] if len(parts) > 6 else "",
                    "pos_2": int(parts[7]) if len(parts) > 7 and parts[7].isdigit() else 0,
                    "strand_2": int(parts[8]) if len(parts) > 8 and parts[8] in ("1", "-1") else 1,
                    "overlap": int(parts[9]) if len(parts) > 9 and parts[9].lstrip("-").isdigit() else 0,
                    "raw": parts,
                }
                records.append(rec)
            elif record_type == "MC":
                # Missing coverage: MC id parent seq_id start end size
                rec = {
                    "type": "MC",
                    "id": parts[1] if len(parts) > 1 else "",
                    "seq_id": parts[3].split(".")[0] if len(parts) > 3 else "",
                    "start": int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0,
                    "end": int(parts[5]) if len(parts) > 5 and parts[5].isdigit() else 0,
                    "raw": parts,
                }
                records.append(rec)
            elif record_type == "UN":
                # Unknown/ambiguous region: UN id parent seq_id start end
                rec = {
                    "type": "UN",
                    "id": parts[1] if len(parts) > 1 else "",
                    "seq_id": parts[3].split(".")[0] if len(parts) > 3 else "",
                    "start": int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0,
                    "end": int(parts[5]) if len(parts) > 5 and parts[5].isdigit() else 0,
                    "raw": parts,
                }
                records.append(rec)
    return records


def count_mc_fp_explained_by_un(test_mc, truth_mc, truth_un, tol_bp=50):
    """
    Of the test MC records that do not match any truth MC record (within
    tol_bp on both start/end), count how many overlap a truth UN (unknown/
    ambiguous, not-confidently-called) region. See RW-003: breseq declines
    to call MC in some ambiguous-coverage regions and reports UN instead,
    so an "MC false positive" there is not necessarily a ProkaDiff defect.
    This is purely descriptive metadata; it does not change tp/fp/fn.
    """
    matched_test_ids = set()
    for tr in truth_mc:
        for t in test_mc:
            if t["id"] in matched_test_ids:
                continue
            if t["seq_id"] == tr["seq_id"] and abs(t["start"] - tr["start"]) <= tol_bp and abs(t["end"] - tr["end"]) <= tol_bp:
                matched_test_ids.add(t["id"])
                break

    fp_records = [t for t in test_mc if t["id"] not in matched_test_ids]
    explained = 0
    for t in fp_records:
        if any(u["seq_id"] == t["seq_id"] and t["start"] <= u["end"] and t["end"] >= u["start"] for u in truth_un):
            explained += 1
    return explained, len(fp_records)


def canonicalize_jc(jc):
    """Canonicalizes a JC record so (side1, side2) ordering is deterministic."""
    s1 = (jc["seq_id_1"], jc["pos_1"], jc["strand_1"])
    s2 = (jc["seq_id_2"], jc["pos_2"], jc["strand_2"])
    if s1 <= s2:
        return (jc["seq_id_1"], jc["pos_1"], jc["strand_1"], jc["seq_id_2"], jc["pos_2"], jc["strand_2"])
    else:
        return (jc["seq_id_2"], jc["pos_2"], -jc["strand_2"], jc["seq_id_1"], jc["pos_1"], -jc["strand_1"])


def match_mutations(test_records, truth_records, jc_tol_bp=5):
    """
    Matches test calls against truth records.
    Returns categorized counts (TP, FP, FN) and breakpoint accuracy stats.
    """
    by_type_test = defaultdict(list)
    by_type_truth = defaultdict(list)

    for r in test_records:
        by_type_test[r["type"]].append(r)
    for r in truth_records:
        by_type_truth[r["type"]].append(r)

    results = {}
    total_exact_bp = 0
    total_1bp = 0
    total_5bp = 0
    total_tested_bp = 0

    all_types = sorted(list(set(by_type_test.keys()) | set(by_type_truth.keys())))

    for mtype in all_types:
        # UN (breseq "unknown/ambiguous, not confidently called") is not a
        # mutation call and is never compared for match/mismatch purposes
        # (see docs/schema.md: "UN 不对拍失败"). It is parsed only so
        # count_mc_fp_explained_by_un() can use it as context for MC FPs.
        if mtype == "UN":
            continue

        tests = by_type_test[mtype]
        truths = by_type_truth[mtype]

        if mtype == "SNP":
            matched_truth = set()
            tp = 0
            for t in tests:
                for idx, tr in enumerate(truths):
                    if idx in matched_truth:
                        continue
                    if t["seq_id"] == tr["seq_id"] and t["pos"] == tr["pos"] and t.get("new_seq") == tr.get("new_seq"):
                        matched_truth.add(idx)
                        tp += 1
                        total_exact_bp += 1
                        total_1bp += 1
                        total_5bp += 1
                        total_tested_bp += 1
                        break
            fp = len(tests) - tp
            fn = len(truths) - tp
            results[mtype] = {"tp": tp, "fp": fp, "fn": fn}

        elif mtype in ("INS", "DEL"):
            matched_truth = set()
            tp = 0
            for t in tests:
                for idx, tr in enumerate(truths):
                    if idx in matched_truth:
                        continue
                    if t["seq_id"] == tr["seq_id"]:
                        diff = abs(t["pos"] - tr["pos"])
                        if diff <= jc_tol_bp:
                            matched_truth.add(idx)
                            tp += 1
                            total_tested_bp += 1
                            if diff == 0:
                                total_exact_bp += 1
                            if diff <= 1:
                                total_1bp += 1
                            if diff <= 5:
                                total_5bp += 1
                            break
            fp = len(tests) - tp
            fn = len(truths) - tp
            results[mtype] = {"tp": tp, "fp": fp, "fn": fn}

        elif mtype == "MOB":
            matched_truth = set()
            tp = 0
            for t in tests:
                for idx, tr in enumerate(truths):
                    if idx in matched_truth:
                        continue
                    if t["seq_id"] == tr["seq_id"]:
                        diff = abs(t["pos"] - tr["pos"])
                        # Check element name match if available
                        elem_match = True
                        if t.get("element") and tr.get("element"):
                            elem_match = (t["element"].lower() == tr["element"].lower())
                        if diff <= jc_tol_bp and elem_match:
                            matched_truth.add(idx)
                            tp += 1
                            total_tested_bp += 1
                            if diff == 0:
                                total_exact_bp += 1
                            if diff <= 1:
                                total_1bp += 1
                            if diff <= 5:
                                total_5bp += 1
                            break
            fp = len(tests) - tp
            fn = len(truths) - tp
            results[mtype] = {"tp": tp, "fp": fp, "fn": fn}

        elif mtype == "JC":
            matched_truth = set()
            tp = 0
            for t in tests:
                c_test = canonicalize_jc(t)
                for idx, tr in enumerate(truths):
                    if idx in matched_truth:
                        continue
                    c_tr = canonicalize_jc(tr)
                    if c_test[0] == c_tr[0] and c_test[3] == c_tr[3] and c_test[2] == c_tr[2] and c_test[5] == c_tr[5]:
                        d1 = abs(c_test[1] - c_tr[1])
                        d2 = abs(c_test[4] - c_tr[4])
                        if d1 <= jc_tol_bp and d2 <= jc_tol_bp:
                            matched_truth.add(idx)
                            tp += 1
                            total_tested_bp += 1
                            max_d = max(d1, d2)
                            if max_d == 0:
                                total_exact_bp += 1
                            if max_d <= 1:
                                total_1bp += 1
                            if max_d <= 5:
                                total_5bp += 1
                            break
            fp = len(tests) - tp
            fn = len(truths) - tp
            results[mtype] = {"tp": tp, "fp": fp, "fn": fn}

        elif mtype == "MC":
            matched_truth = set()
            tp = 0
            for t in tests:
                for idx, tr in enumerate(truths):
                    if idx in matched_truth:
                        continue
                    if t["seq_id"] == tr["seq_id"]:
                        d_start = abs(t["start"] - tr["start"])
                        d_end = abs(t["end"] - tr["end"])
                        if d_start <= 50 and d_end <= 50:
                            matched_truth.add(idx)
                            tp += 1
                            break
            fp = len(tests) - tp
            fn = len(truths) - tp
            results[mtype] = {"tp": tp, "fp": fp, "fn": fn}

        else:
            # Other types
            tp = min(len(tests), len(truths))
            fp = len(tests) - tp
            fn = len(truths) - tp
            results[mtype] = {"tp": tp, "fp": fp, "fn": fn}

    # Summary aggregations
    overall_tp = sum(v["tp"] for v in results.values())
    overall_fp = sum(v["fp"] for v in results.values())
    overall_fn = sum(v["fn"] for v in results.values())

    precision = safe_ratio(overall_tp, overall_tp + overall_fp)
    recall = safe_ratio(overall_tp, overall_tp + overall_fn)
    if isinstance(precision, str) or isinstance(recall, str):
        # F1 is undefined when either parent metric has a zero denominator (M7-P0b).
        f1 = NOT_EVALUABLE
    else:
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    bp_exact = safe_ratio(total_exact_bp, total_tested_bp)
    bp_1bp = safe_ratio(total_1bp, total_tested_bp)
    bp_5bp = safe_ratio(total_5bp, total_tested_bp)

    return {
        "by_type": results,
        "overall_tp": overall_tp,
        "overall_fp": overall_fp,
        "overall_fn": overall_fn,
        "variant_precision": precision,
        "variant_recall": recall,
        "variant_f1": f1,
        "breakpoint_exact_fraction": bp_exact,
        "breakpoint_within_1bp": bp_1bp,
        "breakpoint_within_5bp": bp_5bp,
    }


def main():
    parser = argparse.ArgumentParser(description="Compare ProkaDiff variant calls against truth or oracle")
    parser.add_argument("--test-gd", required=True, help="Test GenomeDiff file (ProkaDiff output)")
    parser.add_argument("--truth-gd", required=True, help="Truth GenomeDiff file (oracle breseq or curated)")
    parser.add_argument("--tolerance", type=int, default=5, help="Breakpoint tolerance in bp (default: 5)")
    parser.add_argument("--tsv", action="store_true", help="Output tab-separated metrics format")
    parser.add_argument("--include-rejected", action="store_true", help="Include records with reject= attributes (default: False)")
    args = parser.parse_args()

    reject_quarantined_truth(args.truth_gd)

    test_records = parse_gd(args.test_gd, include_rejected=args.include_rejected)
    truth_records = parse_gd(args.truth_gd, include_rejected=args.include_rejected)

    metrics = match_mutations(test_records, truth_records, jc_tol_bp=args.tolerance)

    test_mc = [r for r in test_records if r["type"] == "MC"]
    truth_mc = [r for r in truth_records if r["type"] == "MC"]
    truth_un = [r for r in truth_records if r["type"] == "UN"]
    mc_fp_explained_by_un, mc_fp_total = count_mc_fp_explained_by_un(test_mc, truth_mc, truth_un)

    if args.tsv:
        print("metric\tvalue")
        print(f"variant_precision\t{format_metric(metrics['variant_precision'], '.4f')}")
        print(f"variant_recall\t{format_metric(metrics['variant_recall'], '.4f')}")
        print(f"variant_f1\t{format_metric(metrics['variant_f1'], '.4f')}")
        for mtype in ("SNP", "INS", "DEL", "MOB", "JC", "MC"):
            d = metrics["by_type"].get(mtype, {"tp": 0, "fp": 0, "fn": 0})
            tp, fp, fn = d["tp"], d["fp"], d["fn"]
            prec = safe_ratio(tp, tp + fp)
            rec = safe_ratio(tp, tp + fn)
            print(f"{mtype}_tp\t{tp}")
            print(f"{mtype}_fp\t{fp}")
            print(f"{mtype}_fn\t{fn}")
            print(f"{mtype}_precision\t{format_metric(prec, '.4f')}")
            print(f"{mtype}_recall\t{format_metric(rec, '.4f')}")
        print(f"breakpoint_exact_fraction\t{format_metric(metrics['breakpoint_exact_fraction'], '.4f')}")
        print(f"breakpoint_within_1bp\t{format_metric(metrics['breakpoint_within_1bp'], '.4f')}")
        print(f"breakpoint_within_5bp\t{format_metric(metrics['breakpoint_within_5bp'], '.4f')}")
        print(f"MC_fp_explained_by_un\t{mc_fp_explained_by_un}")
        print(f"MC_fp_total\t{mc_fp_total}")
    else:
        p_str = format_metric(metrics["variant_precision"], ".2%")
        r_str = format_metric(metrics["variant_recall"], ".2%")
        print("=== Variant Comparison Summary ===")
        print(f"Precision: {p_str} ({metrics['overall_tp']}/{(metrics['overall_tp'] + metrics['overall_fp'])})")
        print(f"Recall:    {r_str} ({metrics['overall_tp']}/{(metrics['overall_tp'] + metrics['overall_fn'])})")
        print(f"F1 Score:  {format_metric(metrics['variant_f1'], '.4f')}")
        print(
            f"Breakpoint accuracy: Exact={format_metric(metrics['breakpoint_exact_fraction'], '.1%')}, "
            f"<=1bp={format_metric(metrics['breakpoint_within_1bp'], '.1%')}, "
            f"<=5bp={format_metric(metrics['breakpoint_within_5bp'], '.1%')}"
        )
        print("\nBreakdown by Mutation Type:")
        for mtype, d in sorted(metrics["by_type"].items()):
            tp, fp, fn = d["tp"], d["fp"], d["fn"]
            prec = safe_ratio(tp, tp + fp)
            rec = safe_ratio(tp, tp + fn)
            print(f"  {mtype:5s} | TP: {tp:3d} | FP: {fp:3d} | FN: {fn:3d} | Prec: {format_metric(prec, '.1%')} | Rec: {format_metric(rec, '.1%')}")
        if mc_fp_total > 0:
            print(f"\nMC false positives explained by breseq UN (ambiguous) regions: {mc_fp_explained_by_un}/{mc_fp_total}")


if __name__ == "__main__":
    main()
