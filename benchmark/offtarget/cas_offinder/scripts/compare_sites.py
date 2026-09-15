#!/usr/bin/env python3
"""compare_sites.py: Compare candidate off-target sites between Cas-OFFinder and ProkaDiff.

Computes TP, FP, FN, precision, recall for exact site matching.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

def parse_prokadiff_tsv(path: Path) -> set[tuple]:
    """Returns set of (seq_id, start, end, strand, mismatches, target_seq)"""
    sites = set()
    if not path.is_file():
        return sites
    lines = path.read_text(errors="replace").splitlines()
    for line in lines:
        l = line.strip()
        if not l or l.startswith("site_id") or l.startswith("#"):
            continue
        parts = l.split("\t")
        if len(parts) < 8:
            continue
        # header: site_id, seq_id, start, end, strand, target_seq, pam, mismatches, ...
        seq_id = parts[1]
        start = int(parts[2])
        end = int(parts[3])
        strand = parts[4]
        target_seq = parts[5].upper()
        mismatches = int(parts[7])
        sites.add((seq_id, start, end, strand, mismatches, target_seq))
    return sites

def parse_cas_offinder_tsv(path: Path) -> set[tuple]:
    """Returns set of (seq_id, start, end, strand, mismatches, target_seq)"""
    sites = set()
    if not path.is_file():
        return sites
    lines = path.read_text(errors="replace").splitlines()
    for line in lines:
        l = line.strip()
        if not l or l.startswith("#"):
            continue
        parts = l.split("\t")
        if len(parts) < 6:
            continue
        # Cas-OFFinder v2: CrRNA, Chromosome, Location(0-based), Target_Sequence, Direction(+/-), Mismatches
        seq_id = parts[1]
        start = int(parts[2]) + 1
        target_seq = parts[3].upper()
        end = start + len(target_seq) - 1
        strand = parts[4]
        mismatches = int(parts[5])
        sites.add((seq_id, start, end, strand, mismatches, target_seq))
    return sites

def main():
    parser = argparse.ArgumentParser(description="Compare off-target sites between Cas-OFFinder and ProkaDiff")
    parser.add_argument("--oracle", required=True, help="Cas-OFFinder output TSV")
    parser.add_argument("--prokadiff", required=True, help="ProkaDiff off-target output TSV")
    parser.add_argument("--output", default="comparison_summary.tsv", help="Summary TSV output path")
    args = parser.parse_args()

    oracle_sites = parse_cas_offinder_tsv(Path(args.oracle))
    prok_sites = parse_prokadiff_tsv(Path(args.prokadiff))

    tp_sites = oracle_sites.intersection(prok_sites)
    fp_sites = prok_sites.difference(oracle_sites)
    fn_sites = oracle_sites.difference(prok_sites)

    tp = len(tp_sites)
    fp = len(fp_sites)
    fn = len(fn_sites)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 1.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 1.0

    print(f"=== Cas-OFFinder Parity Comparison ===")
    print(f"Oracle sites: {len(oracle_sites)}")
    print(f"ProkaDiff sites: {len(prok_sites)}")
    print(f"True Positives (TP): {tp}")
    print(f"False Positives (FP): {fp}")
    print(f"False Negatives (FN): {fn}")
    print(f"Precision: {precision:.6f}")
    print(f"Recall: {recall:.6f}")

    with open(args.output, "w") as f:
        f.write("oracle_sites\tprokadiff_sites\ttp\tfp\tfn\tprecision\trecall\n")
        f.write(f"{len(oracle_sites)}\t{len(prok_sites)}\t{tp}\t{fp}\t{fn}\t{precision:.6f}\t{recall:.6f}\n")

if __name__ == "__main__":
    main()
