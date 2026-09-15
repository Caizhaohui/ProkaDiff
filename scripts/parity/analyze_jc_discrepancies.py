#!/usr/bin/env python3
"""analyze_jc_discrepancies.py — Classify JC mismatches between ProkaDiff and breseq.

Categories (FIX-016):
  coordinate_jitter     — same sides, coords within N bp
  strand_mismatch       — same coords, different orientation
  duplicate             — ProkaDiff emits same JC twice
  repeat_placement      — JC assigned to wrong copy in repeat region
  low_support_extra     — ProkaDiff-only JC with suspiciously low read support
  breseq_only           — present in breseq, absent in ProkaDiff
  prokadiff_only        — present in ProkaDiff, absent in breseq

Usage:
  python3 analyze_jc_discrepancies.py \
    --prokadiff-gd output/output.gd \
    --breseq-gd    breseq_output/output.gd \
    --outdir       analysis/jc_disc/ \
    [--jitter-bp   10]
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


EVIDENCE_KINDS = {"RA", "MC", "UN"}


def parse_gd_jc(path: Path) -> list[dict]:
    """Extract only JC records from a GenomeDiff file."""
    records = []
    if not path.is_file():
        return records
    for line in path.read_text(errors="replace").splitlines():
        ls = line.strip()
        if not ls or ls.startswith("#") or ls.startswith("="):
            continue
        parts = ls.split("\t")
        if len(parts) < 7 or parts[0] != "JC":
            continue
        records.append({
            "kind": "JC",
            "id": parts[1],
            "fields": parts[3:],
            "raw_line": line,
        })
    return records


def parse_jc_sides(fields: list[str]) -> tuple | None:
    """Parse (seq1, pos1, strand1, seq2, pos2, strand2) from JC fields."""
    if len(fields) < 6:
        return None
    try:
        s1, p1, o1 = fields[0], int(fields[1]), fields[2]
        s2, p2, o2 = fields[3], int(fields[4]), fields[5]
        return s1, p1, o1, s2, p2, o2
    except (ValueError, IndexError):
        return None


def canonical_sides(s1, p1, o1, s2, p2, o2):
    """Return sides in canonical order for comparison."""
    side1 = (s1, o1, p1)
    side2 = (s2, o2, p2)
    if side1 > side2:
        side1, side2 = side2, side1
    return side1, side2


def classify_discrepancies(prok_jc: list[dict], breseq_jc: list[dict], jitter_bp: int = 10):
    results = {
        "coordinate_jitter": [],
        "strand_mismatch": [],
        "duplicate": [],
        "prokadiff_only": [],
        "breseq_only": [],
    }

    # Check for duplicates in ProkaDiff output
    prok_keys = defaultdict(list)
    for r in prok_jc:
        parsed = parse_jc_sides(r["fields"])
        if parsed:
            key = canonical_sides(*parsed)
            prok_keys[key].append(r)
    for key, recs in prok_keys.items():
        if len(recs) > 1:
            results["duplicate"].extend(recs[1:])

    # Match ProkaDiff ↔ breseq
    prok_parsed = [(r, parse_jc_sides(r["fields"])) for r in prok_jc]
    breseq_parsed = [(r, parse_jc_sides(r["fields"])) for r in breseq_jc]

    matched_prok = set()
    matched_breseq = set()

    # Exact match
    breseq_exact = {}
    for j, (r, p) in enumerate(breseq_parsed):
        if p:
            breseq_exact[canonical_sides(*p)] = j

    for i, (pr, pp) in enumerate(prok_parsed):
        if pp is None:
            continue
        key = canonical_sides(*pp)
        if key in breseq_exact:
            matched_prok.add(i)
            matched_breseq.add(breseq_exact[key])

    # Classify unmatched ProkaDiff JC
    for i, (pr, pp) in enumerate(prok_parsed):
        if i in matched_prok or pp is None:
            continue
        ps1, ps2 = canonical_sides(*pp)

        found_jitter = False
        found_strand = False
        for j, (br, bp) in enumerate(breseq_parsed):
            if j in matched_breseq or bp is None:
                continue
            bs1, bs2 = canonical_sides(*bp)

            # Check strand mismatch (same coords, different strand)
            if (ps1[0] == bs1[0] and ps2[0] == bs2[0] and
                    abs(ps1[2] - bs1[2]) <= jitter_bp and abs(ps2[2] - bs2[2]) <= jitter_bp):
                if ps1[1] != bs1[1] or ps2[1] != bs2[1]:
                    results["strand_mismatch"].append({
                        "prok": pr["raw_line"],
                        "breseq": br["raw_line"],
                    })
                    found_strand = True
                    matched_prok.add(i)
                    matched_breseq.add(j)
                    break
                # Same strand, just coord jitter
                elif abs(ps1[2] - bs1[2]) <= jitter_bp and abs(ps2[2] - bs2[2]) <= jitter_bp:
                    results["coordinate_jitter"].append({
                        "prok": pr["raw_line"],
                        "breseq": br["raw_line"],
                        "delta1": abs(ps1[2] - bs1[2]),
                        "delta2": abs(ps2[2] - bs2[2]),
                    })
                    found_jitter = True
                    matched_prok.add(i)
                    matched_breseq.add(j)
                    break

        if not found_jitter and not found_strand and i not in matched_prok:
            results["prokadiff_only"].append(pr)

    # Unmatched breseq JC
    for j, (br, bp) in enumerate(breseq_parsed):
        if j not in matched_breseq:
            results["breseq_only"].append(br)

    return results


def write_report(results: dict, outdir: Path, jitter_bp: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    summary = []
    for cat, items in results.items():
        summary.append({"category": cat, "count": len(items)})
        if items:
            cat_file = outdir / f"{cat}.txt"
            with open(cat_file, "w") as f:
                for item in items:
                    if isinstance(item, dict) and "prok" in item:
                        f.write(f"ProkaDiff: {item['prok']}\n")
                        f.write(f"breseq:    {item.get('breseq', '')}\n")
                        if "delta1" in item:
                            f.write(f"delta:     side1={item['delta1']}bp, side2={item['delta2']}bp\n")
                        f.write("\n")
                    elif isinstance(item, dict) and "raw_line" in item:
                        f.write(item["raw_line"] + "\n")

    summary_file = outdir / "jc_discrepancy_summary.tsv"
    with open(summary_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["category", "count"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    print(f"[analyze_jc_discrepancies] jitter_bp={jitter_bp}")
    for row in summary:
        print(f"  {row['category']:30s}: {row['count']}")
    print(f"  Written to {outdir}/")


def main():
    parser = argparse.ArgumentParser(description="Classify JC discrepancies between ProkaDiff and breseq")
    parser.add_argument("--prokadiff-gd", required=True, type=Path)
    parser.add_argument("--breseq-gd", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--jitter-bp", type=int, default=10)
    args = parser.parse_args()

    prok_jc = parse_gd_jc(args.prokadiff_gd)
    breseq_jc = parse_gd_jc(args.breseq_gd)

    print(f"[analyze_jc_discrepancies] ProkaDiff JC: {len(prok_jc)}, breseq JC: {len(breseq_jc)}")

    results = classify_discrepancies(prok_jc, breseq_jc, args.jitter_bp)
    write_report(results, args.outdir, args.jitter_bp)
    return 0


if __name__ == "__main__":
    sys.exit(main())
