#!/usr/bin/env python3
"""compare_gd.py: Two-tier parity comparison between ProkaDiff and breseq GenomeDiff outputs.

Level A (Strict): Exact coordinate, type, and allele/feature match.
Level B (Tolerant): Diagnostic coordinate tolerance for JC, MOB, and structural events.

Outputs:
  - over_call.gd (mutations in ProkaDiff not in breseq)
  - under_call.gd (mutations in breseq not in ProkaDiff)
  - parity.tsv (single-dataset parity summary row)
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

MUT_KINDS = ["SNP", "SUB", "INS", "DEL", "MOB", "AMP", "CON", "INV"]
JC_KIND = "JC"
EVIDENCE_KINDS = {"RA", "MC", "UN"}


def parse_time_file(path: Path | None):
    wall_s = ""
    rss_kb = ""
    if not path or not path.is_file():
        return wall_s, rss_kb
    text = path.read_text(errors="replace")
    for line in text.splitlines():
        if "Elapsed (wall clock) time" in line:
            m = re.search(r":\s*(\d+):(\d+):(\d+(?:\.\d+)?)$", line)
            if m:
                wall_s = f"{int(m.group(1)) * 3600 + int(m.group(2)) * 60 + float(m.group(3)):.2f}"
            else:
                m = re.search(r":\s*(\d+):(\d+(?:\.\d+)?)$", line)
                if m:
                    wall_s = f"{int(m.group(1)) * 60 + float(m.group(2)):.2f}"
        if "Maximum resident set size" in line:
            m = re.search(r":\s*(\d+)", line)
            if m:
                rss_kb = m.group(1)
    return wall_s, rss_kb


def get_cmd_version(cmd: list[str]) -> str:
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=10)
        lines = res.stdout.strip().splitlines()
        return lines[0] if lines else "unknown"
    except Exception:
        return "unavailable"


def get_git_commit(repo_root: Path) -> str:
    try:
        res = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=repo_root, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=5)
        return res.stdout.strip() if res.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def parse_gd(path: Path):
    """Parses GenomeDiff file into records."""
    records = []
    headers = []
    if not path.is_file():
        return headers, records
    for line in path.read_text(errors="replace").splitlines():
        line_s = line.strip()
        if not line_s or line_s.startswith("#") or line_s.startswith("="):
            headers.append(line)
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        kind = parts[0]
        # Ignore raw evidence lines for strict mutation parity, but keep JC
        if kind in EVIDENCE_KINDS:
            continue
        records.append({
            "kind": kind,
            "id": parts[1],
            "parent_ids": parts[2],
            "fields": parts[3:],
            "raw_line": line,
        })
    return headers, records


def strict_key(rec):
    return (rec["kind"], tuple(rec["fields"]))


def parse_jc_coords(fields):
    """(seq1, pos1, strand1, seq2, pos2, strand2)"""
    if len(fields) < 6:
        return None
    try:
        s1, p1, o1 = fields[0], int(fields[1]), fields[2]
        s2, p2, o2 = fields[3], int(fields[4]), fields[5]
        side1 = (s1, o1, p1)
        side2 = (s2, o2, p2)
        if side1 > side2:
            side1, side2 = side2, side1
        return side1, side2
    except (ValueError, IndexError):
        return None


def jc_close(a, b, tol_bp):
    ca = parse_jc_coords(a["fields"])
    cb = parse_jc_coords(b["fields"])
    if not ca or not cb:
        return False
    (a1, a2), (b1, b2) = ca, cb
    if (a1[0], a1[1], a2[0], a2[1]) != (b1[0], b1[1], b2[0], b2[1]):
        return False
    return abs(a1[2] - b1[2]) <= tol_bp and abs(a2[2] - b2[2]) <= tol_bp


def compare_gd_records(prok_records, breseq_records, tol_bp=5):
    """Computes strict and tolerant differences."""
    breseq_keys = {strict_key(r): r for r in breseq_records}
    prok_keys = {strict_key(r): r for r in prok_records}

    strict_over = [r for r in prok_records if strict_key(r) not in breseq_keys]
    strict_under = [r for r in breseq_records if strict_key(r) not in prok_keys]

    # Per-type counts
    over_by_type = {}
    under_by_type = {}
    for k in MUT_KINDS + [JC_KIND]:
        over_by_type[k] = sum(1 for r in strict_over if r["kind"] == k)
        under_by_type[k] = sum(1 for r in strict_under if r["kind"] == k)

    # Diagnostic tolerant matching for JC
    over_jc = [r for r in strict_over if r["kind"] == JC_KIND]
    under_jc = [r for r in strict_under if r["kind"] == JC_KIND]
    used_under_jc = set()
    matched_jc_tol = 0
    for oj in over_jc:
        for idx, uj in enumerate(under_jc):
            if idx in used_under_jc:
                continue
            if jc_close(oj, uj, tol_bp):
                used_under_jc.add(idx)
                matched_jc_tol += 1
                break

    return {
        "strict_over": strict_over,
        "strict_under": strict_under,
        "over_by_type": over_by_type,
        "under_by_type": under_by_type,
        "matched_jc_tol": matched_jc_tol,
        "unmatched_over_jc_tol": len(over_jc) - matched_jc_tol,
        "unmatched_under_jc_tol": len(under_jc) - matched_jc_tol,
    }


def run_gdtools_subtract(file1: Path, file2: Path, outpath: Path, gdtools_bin="gdtools"):
    """Uses gdtools SUBTRACT file1 file2 > outpath if available."""
    if shutil.which(gdtools_bin):
        try:
            with open(outpath, "w") as out:
                subprocess.run([gdtools_bin, "SUBTRACT", str(file1), str(file2)], stdout=out, check=True)
            return True
        except Exception as e:
            sys.stderr.write(f"[compare_gd] gdtools SUBTRACT failed: {e}\n")
    return False


def main():
    parser = argparse.ArgumentParser(description="Compare ProkaDiff and breseq GenomeDiff outputs.")
    parser.add_argument("--prokadiff-gd", required=True, type=Path)
    parser.add_argument("--breseq-gd", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--jc-tol-bp", type=int, default=5)
    parser.add_argument("--prokadiff-time", type=Path, default=None)
    parser.add_argument("--breseq-time", type=Path, default=None)
    parser.add_argument("--prokadiff-commit", default=None)
    parser.add_argument("--breseq-version", default=None)
    parser.add_argument("--bowtie2-version", default=None)
    parser.add_argument("--gdtools-bin", default="gdtools")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    if not args.prokadiff_gd.is_file():
        sys.stderr.write(f"Error: missing prokadiff GD {args.prokadiff_gd}\n")
        return 1
    if not args.breseq_gd.is_file():
        sys.stderr.write(f"Error: missing breseq GD {args.breseq_gd}\n")
        return 1

    _, prok_records = parse_gd(args.prokadiff_gd)
    _, breseq_records = parse_gd(args.breseq_gd)

    comp = compare_gd_records(prok_records, breseq_records, args.jc_tol_bp)

    # Write over_call.gd and under_call.gd
    over_gd_path = args.outdir / "over_call.gd"
    under_gd_path = args.outdir / "under_call.gd"

    # Prefer gdtools SUBTRACT if possible, fallback to Python-generated GD
    used_gdtools = run_gdtools_subtract(args.prokadiff_gd, args.breseq_gd, over_gd_path, args.gdtools_bin)
    run_gdtools_subtract(args.breseq_gd, args.prokadiff_gd, under_gd_path, args.gdtools_bin)

    if not used_gdtools:
        with open(over_gd_path, "w") as f:
            f.write("#=GENOMEDIFF\n")
            for r in comp["strict_over"]:
                f.write(r["raw_line"] + "\n")
        with open(under_gd_path, "w") as f:
            f.write("#=GENOMEDIFF\n")
            for r in comp["strict_under"]:
                f.write(r["raw_line"] + "\n")

    # Time and environment
    wall_prok, rss_prok = parse_time_file(args.prokadiff_time or (args.outdir / "prokadiff.time"))
    wall_breseq, rss_breseq = parse_time_file(args.breseq_time or (args.outdir / "breseq.time"))

    repo_root = Path(__file__).resolve().parent.parent.parent
    prok_commit = args.prokadiff_commit or get_git_commit(repo_root)
    breseq_v = args.breseq_version or get_cmd_version(["breseq", "--version"])
    bt2_v = args.bowtie2_version or get_cmd_version(["bowtie2", "--version"])

    # Build parity.tsv dictionary
    tsv_row = {
        "dataset": args.dataset,
        "breseq_version": breseq_v,
        "bowtie2_version": bt2_v,
        "prokadiff_commit": prok_commit,
        "breseq_mutations": len(breseq_records),
        "prokadiff_mutations": len(prok_records),
        "strict_over": len(comp["strict_over"]),
        "strict_under": len(comp["strict_under"]),
        "snp_over": comp["over_by_type"].get("SNP", 0),
        "snp_under": comp["under_by_type"].get("SNP", 0),
        "sub_over": comp["over_by_type"].get("SUB", 0),
        "sub_under": comp["under_by_type"].get("SUB", 0),
        "ins_over": comp["over_by_type"].get("INS", 0),
        "ins_under": comp["under_by_type"].get("INS", 0),
        "del_over": comp["over_by_type"].get("DEL", 0),
        "del_under": comp["under_by_type"].get("DEL", 0),
        "mob_over": comp["over_by_type"].get("MOB", 0),
        "mob_under": comp["under_by_type"].get("MOB", 0),
        "amp_over": comp["over_by_type"].get("AMP", 0),
        "amp_under": comp["under_by_type"].get("AMP", 0),
        "jc_over": comp["over_by_type"].get("JC", 0),
        "jc_under": comp["under_by_type"].get("JC", 0),
        "wall_breseq": wall_breseq if wall_breseq else "NA",
        "wall_prokadiff": wall_prok if wall_prok else "NA",
        "rss_breseq": rss_breseq if rss_breseq else "NA",
        "rss_prokadiff": rss_prok if rss_prok else "NA",
    }

    parity_tsv_path = args.outdir / "parity.tsv"
    with open(parity_tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(tsv_row.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(tsv_row)

    # Also diagnostic json
    diag_path = args.outdir / "parity_diagnostic.json"
    with open(diag_path, "w") as f:
        json.dump({
            "summary": tsv_row,
            "matched_jc_tol": comp["matched_jc_tol"],
            "unmatched_over_jc_tol": comp["unmatched_over_jc_tol"],
            "unmatched_under_jc_tol": comp["unmatched_under_jc_tol"],
        }, f, indent=2)

    print(f"[compare_gd] Parity evaluation for dataset '{args.dataset}':")
    print(f"  breseq mutations: {len(breseq_records)}, prokadiff mutations: {len(prok_records)}")
    print(f"  strict_over: {len(comp['strict_over'])}, strict_under: {len(comp['strict_under'])}")
    print(f"  JC matched with tol={args.jc_tol_bp}bp: {comp['matched_jc_tol']}")
    print(f"  Wrote {over_gd_path}, {under_gd_path}, {parity_tsv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
