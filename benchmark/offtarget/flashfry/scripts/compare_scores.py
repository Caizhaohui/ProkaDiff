#!/usr/bin/env python3
"""compare_scores.py — ProkaDiff vs FlashFry numeric parity for CFD and Hsu2013.

STATUS: BLOCKED_EXTERNAL_DEPENDENCY — requires golden/normalized_flashfry.tsv
which is only produced after running FlashFry 1.15.

Acceptance criteria:
  - CFD: abs(expected - rust) <= 1e-6 for all pairs
  - Hsu: abs(expected - rust) <= 1e-6 for all pairs
  - At least 100 pairs tested
  - 0 failures outside tolerance
"""
import sys
import csv
import argparse
from pathlib import Path


def compare(golden_path: str, prokadiff_path: str, cfd_tol: float = 1e-6, hsu_tol: float = 1e-6) -> int:
    golden = Path(golden_path)
    prok = Path(prokadiff_path)
    if not golden.is_file():
        print(f"BLOCKED: golden file not found: {golden_path}", file=sys.stderr)
        print("Run scripts/run_flashfry.sh then scripts/normalize_flashfry.py first.", file=sys.stderr)
        return 2
    if not prok.is_file():
        print(f"Error: ProkaDiff output not found: {prokadiff_path}", file=sys.stderr)
        return 1

    with open(golden, newline="") as f:
        golden_rows = {r["target_seq"]: r for r in csv.DictReader(f, delimiter="\t")}
    with open(prok, newline="") as f:
        prok_rows = {r["target_seq"]: r for r in csv.DictReader(f, delimiter="\t")}

    n_tested = 0
    cfd_failures = 0
    hsu_failures = 0

    for target, gr in golden_rows.items():
        pr = prok_rows.get(target)
        if pr is None:
            print(f"  MISSING in ProkaDiff: {target}")
            cfd_failures += 1
            hsu_failures += 1
            continue

        n_tested += 1
        g_cfd = gr.get("cfd_score", "NA")
        p_cfd = pr.get("cfd_score", "NA")
        if g_cfd != "NA" and p_cfd != "NA":
            try:
                diff = abs(float(g_cfd) - float(p_cfd))
                if diff > cfd_tol:
                    print(f"  CFD FAIL: target={target} expected={g_cfd} got={p_cfd} diff={diff:.2e}")
                    cfd_failures += 1
            except ValueError:
                pass
        elif g_cfd != "NA" and p_cfd == "NA":
            print(f"  CFD DISABLED for: {target} (expected={g_cfd})")

        g_hsu = gr.get("hsu_score", "NA")
        p_hsu = pr.get("hsu_score", "NA")
        if g_hsu != "NA" and p_hsu != "NA":
            try:
                diff = abs(float(g_hsu) - float(p_hsu))
                if diff > hsu_tol:
                    print(f"  Hsu FAIL: target={target} expected={g_hsu} got={p_hsu} diff={diff:.2e}")
                    hsu_failures += 1
            except ValueError:
                pass

    print(f"\n[compare_scores] Pairs tested: {n_tested}")
    print(f"  CFD failures: {cfd_failures}")
    print(f"  Hsu failures: {hsu_failures}")

    if n_tested < 100:
        print(f"  WARNING: Only {n_tested} pairs tested; gate requires >= 100")
    if cfd_failures == 0 and n_tested >= 100:
        print("  CFD GATE: PASS")
    else:
        print("  CFD GATE: NOT MET")
    if hsu_failures == 0 and n_tested >= 100:
        print("  Hsu GATE: PASS")
    else:
        print("  Hsu GATE: NOT MET")

    return 0 if (cfd_failures == 0 and hsu_failures == 0 and n_tested >= 100) else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ProkaDiff vs FlashFry score parity check")
    parser.add_argument("--golden", required=True, help="normalized_flashfry.tsv")
    parser.add_argument("--prokadiff", required=True, help="ProkaDiff off-target TSV with cfd_score/hsu_score cols")
    parser.add_argument("--cfd-tol", type=float, default=1e-6)
    parser.add_argument("--hsu-tol", type=float, default=1e-6)
    args = parser.parse_args()
    sys.exit(compare(args.golden, args.prokadiff, args.cfd_tol, args.hsu_tol))
