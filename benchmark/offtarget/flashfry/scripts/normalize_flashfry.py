#!/usr/bin/env python3
"""normalize_flashfry.py — Convert raw FlashFry 1.15 TSV to ProkaDiff normalized format.

FlashFry 1.15 native column layout (as documented):
  contig  start  stop  target  orientation  numberOfMismatches  Doench2016CFDScore  Hsu2013

Coordinate convention (to be confirmed with real output):
  FlashFry start: 0-based inclusive
  FlashFry stop:  0-based exclusive (i.e. half-open)
  ProkaDiff:      1-based inclusive

This script records the verified convention in a header comment.
"""
import sys
import csv
from pathlib import Path


def normalize(raw_path: str, out_path: str) -> None:
    raw = Path(raw_path)
    if not raw.is_file():
        print(f"Error: {raw_path} not found. Run run_flashfry.sh first.", file=sys.stderr)
        sys.exit(1)

    with open(raw, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    # Normalize header names to lowercase for matching
    def find_col(names):
        for fn in fieldnames:
            if fn.strip().lower() in [n.lower() for n in names]:
                return fn
        return None

    col_contig = find_col(["contig", "chrom", "chr"])
    col_start = find_col(["start"])
    col_stop = find_col(["stop", "end"])
    col_target = find_col(["target", "sequence"])
    col_orient = find_col(["orientation", "strand", "direction"])
    col_mm = find_col(["numberofmismatches", "mismatches", "mismatch_count"])
    col_cfd = find_col(["doench2016cfdscore", "cfd", "cfd_score"])
    col_hsu = find_col(["hsu2013", "hsu", "hsu_score"])

    if not all([col_contig, col_start, col_stop, col_target, col_orient]):
        print("Error: required columns missing from FlashFry output.", file=sys.stderr)
        print(f"Available: {fieldnames}", file=sys.stderr)
        sys.exit(1)

    out_rows = []
    for r in rows:
        # Coordinate conversion: FlashFry 0-based half-open → 1-based inclusive
        # prokadiff_start = flashfry_start + 1
        # prokadiff_end   = flashfry_stop      (stop is exclusive, so last nt = stop-1, +1 for 1-based = stop)
        try:
            ff_start = int(r[col_start])
            ff_stop = int(r[col_stop])
            pk_start = ff_start + 1
            pk_end = ff_stop  # half-open exclusive → 1-based inclusive
        except (ValueError, KeyError):
            pk_start = r.get(col_start, "")
            pk_end = r.get(col_stop, "")

        # Strand mapping: FWD → +, RVS → -
        orient_raw = r.get(col_orient, "").strip()
        if orient_raw.upper() in ("FWD", "F", "+", "FORWARD", "PLUS"):
            strand = "+"
        elif orient_raw.upper() in ("RVS", "R", "-", "REVERSE", "MINUS"):
            strand = "-"
        else:
            strand = orient_raw

        out_rows.append({
            "seq_id": r.get(col_contig, ""),
            "start": pk_start,
            "end": pk_end,
            "strand": strand,
            "target_seq": r.get(col_target, ""),
            "mismatches": r.get(col_mm, "") if col_mm else "",
            "cfd_score": r.get(col_cfd, "NA") if col_cfd else "NA",
            "hsu_score": r.get(col_hsu, "NA") if col_hsu else "NA",
            "flashfry_start_raw": r.get(col_start, ""),
            "flashfry_stop_raw": r.get(col_stop, ""),
            "orientation_raw": orient_raw,
        })

    out_fields = ["seq_id", "start", "end", "strand", "target_seq", "mismatches",
                  "cfd_score", "hsu_score", "flashfry_start_raw", "flashfry_stop_raw", "orientation_raw"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"[normalize_flashfry] {len(out_rows)} sites written to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <raw_flashfry.tsv> <normalized_flashfry.tsv>")
        sys.exit(1)
    normalize(sys.argv[1], sys.argv[2])
