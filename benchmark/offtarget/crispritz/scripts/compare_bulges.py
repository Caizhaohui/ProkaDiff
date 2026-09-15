#!/usr/bin/env python3
"""compare_bulges.py — ProkaDiff vs CRISPRitz site-level parity for bulge search.

STATUS: BLOCKED_EXTERNAL_DEPENDENCY — requires golden/normalized_crispritz.tsv

Acceptance criteria:
  precision = 1.000 (no FP)
  recall    = 1.000 (no FN)

Any discrepancy is written to benchmark/offtarget/BULGE_ISSUES.md format.
"""
import sys
import csv
import argparse
from pathlib import Path


def site_key(row):
    """Canonical key for site matching."""
    return (
        row.get("seq_id", ""),
        row.get("start", ""),
        row.get("end", ""),
        row.get("strand", ""),
        row.get("bulge_type", ""),
        row.get("bulge_size", ""),
    )


def compare(golden_path: str, prokadiff_path: str) -> int:
    golden = Path(golden_path)
    prok = Path(prokadiff_path)

    if not golden.is_file():
        print(f"BLOCKED: golden file not found: {golden_path}", file=sys.stderr)
        print("Run scripts/run_crispritz.sh then scripts/normalize_crispritz.py first.", file=sys.stderr)
        return 2

    if not prok.is_file():
        print(f"Error: ProkaDiff output not found: {prokadiff_path}", file=sys.stderr)
        return 1

    with open(golden, newline="") as f:
        golden_sites = set(site_key(r) for r in csv.DictReader(f, delimiter="\t"))
    with open(prok, newline="") as f:
        prok_sites = set(site_key(r) for r in csv.DictReader(f, delimiter="\t"))

    tp = golden_sites & prok_sites
    fp = prok_sites - golden_sites
    fn = golden_sites - prok_sites

    n_golden = len(golden_sites)
    n_prok = len(prok_sites)
    precision = len(tp) / n_prok if n_prok > 0 else 0.0
    recall = len(tp) / n_golden if n_golden > 0 else 0.0

    print(f"[compare_bulges] CRISPRitz sites: {n_golden}, ProkaDiff sites: {n_prok}")
    print(f"  TP: {len(tp)}, FP: {len(fp)}, FN: {len(fn)}")
    print(f"  Precision: {precision:.6f}, Recall: {recall:.6f}")

    if fp:
        print("\nFalse Positives (ProkaDiff-only):")
        for s in sorted(fp):
            print(f"  {s}")

    if fn:
        print("\nFalse Negatives (CRISPRitz-only):")
        for s in sorted(fn):
            print(f"  {s}")

    if precision == 1.0 and recall == 1.0:
        print("\nBULGE GATE: PASS")
        return 0
    else:
        print("\nBULGE GATE: NOT MET")
        return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ProkaDiff vs CRISPRitz bulge site parity")
    parser.add_argument("--golden", required=True, help="normalized_crispritz.tsv")
    parser.add_argument("--prokadiff", required=True, help="ProkaDiff bulge search TSV")
    args = parser.parse_args()
    sys.exit(compare(args.golden, args.prokadiff))
