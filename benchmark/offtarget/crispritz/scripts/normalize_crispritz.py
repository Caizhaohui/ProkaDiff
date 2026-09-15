#!/usr/bin/env python3
"""normalize_crispritz.py — Convert raw CRISPRitz 2.8.1 output to ProkaDiff normalized format.

CRISPRitz output format (tab-separated, typical columns):
  crRNA  DNA_bulge  RNA_bulge  Chromosome  Start  End  Strand  Mismatches  Bulge_type  Bulge_size  ...

This script produces a TSV with columns matching the ProkaDiff off-target schema.
"""
import sys
import csv
from pathlib import Path


def normalize(raw_path: str, out_path: str) -> None:
    raw = Path(raw_path)
    if not raw.is_file():
        print(f"BLOCKED: {raw_path} not found.", file=sys.stderr)
        print("Run scripts/run_crispritz.sh first.", file=sys.stderr)
        sys.exit(2)

    with open(raw, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    def find(names):
        for fn in fieldnames:
            if fn.strip().lower() in [n.lower() for n in names]:
                return fn
        return None

    # CRISPRitz 2.8.1 typical columns (lowercase matching)
    col_crrna = find(["crrna", "guide", "spacer"])
    col_chrom = find(["chromosome", "chr", "chrom", "seq_id", "contig"])
    col_start = find(["start", "pos"])
    col_end = find(["end", "stop"])
    col_strand = find(["strand", "direction"])
    col_mm = find(["mismatches", "mismatch_count", "mm"])
    col_btype = find(["bulge_type", "bulge type", "btype"])
    col_bsize = find(["bulge_size", "bulge size", "bsize", "dna_bulge", "rna_bulge"])

    out_rows = []
    for r in rows:
        # Strand normalization: crispritz uses +/-
        strand_raw = r.get(col_strand, "").strip() if col_strand else ""
        if strand_raw in ("+", "1", "fwd", "forward"):
            strand = "+"
        elif strand_raw in ("-", "-1", "rvs", "reverse"):
            strand = "-"
        else:
            strand = strand_raw

        # Bulge type normalization
        btype_raw = r.get(col_btype, "none").strip().lower() if col_btype else "none"
        if "dna" in btype_raw:
            bulge_type = "DNA"
        elif "rna" in btype_raw:
            bulge_type = "RNA"
        else:
            bulge_type = "none"

        out_rows.append({
            "seq_id": r.get(col_chrom, "") if col_chrom else "",
            "start": r.get(col_start, "") if col_start else "",
            "end": r.get(col_end, "") if col_end else "",
            "strand": strand,
            "mismatches": r.get(col_mm, "") if col_mm else "",
            "bulge_type": bulge_type,
            "bulge_size": r.get(col_bsize, "0") if col_bsize else "0",
            "guide": r.get(col_crrna, "") if col_crrna else "",
        })

    out_fields = ["seq_id", "start", "end", "strand", "mismatches", "bulge_type", "bulge_size", "guide"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"[normalize_crispritz] {len(out_rows)} sites written to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <raw_crispritz.txt> <normalized_crispritz.tsv>")
        sys.exit(1)
    normalize(sys.argv[1], sys.argv[2])
