#!/usr/bin/env python3
"""summarize_parity.py: Aggregates parity.tsv rows into benchmark/results/parity_summary.tsv.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

COLUMNS = [
    "dataset",
    "breseq_version",
    "bowtie2_version",
    "prokadiff_commit",
    "breseq_mutations",
    "prokadiff_mutations",
    "strict_over",
    "strict_under",
    "snp_over",
    "snp_under",
    "sub_over",
    "sub_under",
    "ins_over",
    "ins_under",
    "del_over",
    "del_under",
    "mob_over",
    "mob_under",
    "amp_over",
    "amp_under",
    "jc_over",
    "jc_under",
    "wall_breseq",
    "wall_prokadiff",
    "rss_breseq",
    "rss_prokadiff",
]


def main():
    parser = argparse.ArgumentParser(description="Aggregate parity.tsv files into parity_summary.tsv.")
    parser.add_argument("inputs", nargs="*", type=Path, help="Paths to parity.tsv files or directories containing parity.tsv")
    parser.add_argument("--output", "-o", type=Path, default=Path("benchmark/results/parity_summary.tsv"))
    args = parser.parse_args()

    tsv_files = []
    for inp in args.inputs:
        if inp.is_file() and inp.name.endswith(".tsv"):
            tsv_files.append(inp)
        elif inp.is_dir():
            for p in inp.rglob("parity.tsv"):
                tsv_files.append(p)

    if not tsv_files and not args.inputs:
        # Default scan benchmark/results/
        default_dir = Path("benchmark/results")
        if default_dir.is_dir():
            tsv_files = list(default_dir.rglob("parity.tsv"))

    print(f"[summarize_parity] Found {len(tsv_files)} parity TSV files.")

    rows = []
    seen = set()
    for tsv in sorted(tsv_files):
        try:
            with open(tsv, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for r in reader:
                    # Filter/validate columns
                    cleaned = {col: r.get(col, "NA") for col in COLUMNS}
                    key = (cleaned["dataset"], cleaned["prokadiff_commit"], cleaned["breseq_version"])
                    if key not in seen:
                        seen.add(key)
                        rows.append(cleaned)
        except Exception as e:
            sys.stderr.write(f"Warning: could not read {tsv}: {e}\n")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"[summarize_parity] Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
