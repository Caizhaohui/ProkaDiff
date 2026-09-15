#!/usr/bin/env python3
"""compare_gd.py: Two-tier parity comparison between ProkaDiff and breseq GenomeDiff outputs.

Level A (Strict): Exact coordinate, type, and allele/feature match.
  - Primary gate: gdtools SUBTRACT output (when gdtools is available)
  - Diagnostic: Python exact record comparison (record_exact_over/under)
Level B (Tolerant): Diagnostic coordinate tolerance for JC events.
Level C (Normalized): INS/DEL coordinate normalization for homopolymer/repeat regions.

Outputs:
  - over_call.gd (mutations in ProkaDiff not in breseq — via gdtools SUBTRACT when available)
  - under_call.gd (mutations in breseq not in ProkaDiff — via gdtools SUBTRACT when available)
  - parity.tsv (single-dataset parity summary row)
  - parity_diagnostic.json (detailed breakdown)
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


def normalize_indel_coord(kind: str, seq_id: str, pos: int, alt: str) -> tuple:
    """Normalise INS/DEL to right-aligned canonical form for repeat regions.
    Currently a stub — full normalization requires reference sequence.
    Returns (kind, seq_id, pos, alt) unchanged; replace with reference-aware
    logic once prokadiff normalize-gd is available.
    """
    return (kind, seq_id, pos, alt)


def normalized_key(rec):
    """Attempt normalized comparison key for INS/DEL (right-align stub)."""
    kind = rec["kind"]
    if kind in ("INS", "DEL") and len(rec["fields"]) >= 3:
        seq_id = rec["fields"][0]
        try:
            pos = int(rec["fields"][1])
        except ValueError:
            return strict_key(rec)
        alt = rec["fields"][2] if len(rec["fields"]) > 2 else ""
        return normalize_indel_coord(kind, seq_id, pos, alt)
    return strict_key(rec)


def compare_gd_records(prok_records, breseq_records, tol_bp=5):
    """Computes strict, record-exact, tolerant, and normalized differences."""
    breseq_keys = {strict_key(r): r for r in breseq_records}
    prok_keys = {strict_key(r): r for r in prok_records}

    # ── Tier 1: Python record-exact comparison (diagnostic) ─────────────────
    record_exact_over = [r for r in prok_records if strict_key(r) not in breseq_keys]
    record_exact_under = [r for r in breseq_records if strict_key(r) not in prok_keys]

    # ── Tier 2: Normalized parity (INS/DEL right-align stub) ────────────────
    breseq_norm_keys = {normalized_key(r): r for r in breseq_records}
    prok_norm_keys = {normalized_key(r): r for r in prok_records}
    normalized_over = [r for r in prok_records if normalized_key(r) not in breseq_norm_keys]
    normalized_under = [r for r in breseq_records if normalized_key(r) not in prok_norm_keys]

    # Per-type counts (using record_exact for per-type breakdown)
    over_by_type = {}
    under_by_type = {}
    for k in MUT_KINDS + [JC_KIND]:
        over_by_type[k] = sum(1 for r in record_exact_over if r["kind"] == k)
        under_by_type[k] = sum(1 for r in record_exact_under if r["kind"] == k)

    # ── Diagnostic tolerant matching for JC ─────────────────────────────────
    over_jc = [r for r in record_exact_over if r["kind"] == JC_KIND]
    under_jc = [r for r in record_exact_under if r["kind"] == JC_KIND]
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
        "record_exact_over": record_exact_over,
        "record_exact_under": record_exact_under,
        "normalized_over": normalized_over,
        "normalized_under": normalized_under,
        "over_by_type": over_by_type,
        "under_by_type": under_by_type,
        "matched_jc_tol": matched_jc_tol,
        "unmatched_over_jc_tol": len(over_jc) - matched_jc_tol,
        "unmatched_under_jc_tol": len(under_jc) - matched_jc_tol,
    }


def run_gdtools_subtract(file1: Path, file2: Path, outpath: Path, gdtools_bin="gdtools") -> tuple[bool, int]:
    """Uses gdtools SUBTRACT file1 file2 > outpath if available.
    Returns (success, count_records).
    """
    gdtools_path = shutil.which(gdtools_bin)
    if gdtools_path:
        try:
            with open(outpath, "w") as out:
                result = subprocess.run(
                    [gdtools_path, "SUBTRACT", str(file1), str(file2)],
                    stdout=out, check=True, stderr=subprocess.PIPE
                )
            # Count non-header, non-comment, non-empty lines (mutations)
            count = 0
            for line in outpath.read_text(errors="replace").splitlines():
                ls = line.strip()
                if ls and not ls.startswith("#") and not ls.startswith("="):
                    parts = ls.split("\t")
                    if len(parts) >= 4 and parts[0] not in EVIDENCE_KINDS:
                        count += 1
            return True, count
        except Exception as e:
            sys.stderr.write(f"[compare_gd] gdtools SUBTRACT failed: {e}\n")
    return False, 0


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

    # ── Write over_call.gd and under_call.gd via gdtools SUBTRACT ───────────
    over_gd_path = args.outdir / "over_call.gd"
    under_gd_path = args.outdir / "under_call.gd"

    # Try gdtools SUBTRACT (Tier 1 primary gate)
    gdtools_used, gdtools_over_count = run_gdtools_subtract(
        args.prokadiff_gd, args.breseq_gd, over_gd_path, args.gdtools_bin
    )
    _, gdtools_under_count = run_gdtools_subtract(
        args.breseq_gd, args.prokadiff_gd, under_gd_path, args.gdtools_bin
    )

    gdtools_status = "available" if gdtools_used else "unavailable"

    if not gdtools_used:
        # Fallback: write Python-computed over/under (for diagnostic only)
        with open(over_gd_path, "w") as f:
            f.write("#=GENOMEDIFF\n")
            for r in comp["record_exact_over"]:
                f.write(r["raw_line"] + "\n")
        with open(under_gd_path, "w") as f:
            f.write("#=GENOMEDIFF\n")
            for r in comp["record_exact_under"]:
                f.write(r["raw_line"] + "\n")
        # When gdtools unavailable, count from Python comparison
        gdtools_over_count = len(comp["record_exact_over"])
        gdtools_under_count = len(comp["record_exact_under"])

    # ── Time and environment ─────────────────────────────────────────────────
    wall_prok, rss_prok = parse_time_file(args.prokadiff_time or (args.outdir / "prokadiff.time"))
    wall_breseq, rss_breseq = parse_time_file(args.breseq_time or (args.outdir / "breseq.time"))

    repo_root = Path(__file__).resolve().parent.parent.parent
    prok_commit = args.prokadiff_commit or get_git_commit(repo_root)
    breseq_v = args.breseq_version or get_cmd_version(["breseq", "--version"])
    bt2_v = args.bowtie2_version or get_cmd_version(["bowtie2", "--version"])

    # Resolve gdtools version
    gdtools_bin_path = shutil.which(args.gdtools_bin)
    gdtools_v = get_cmd_version([args.gdtools_bin, "--version"]) if gdtools_bin_path else "unavailable"

    # ── Build parity.tsv dictionary ──────────────────────────────────────────
    tsv_row = {
        "dataset": args.dataset,
        "breseq_version": breseq_v,
        "bowtie2_version": bt2_v,
        "gdtools_version": gdtools_v,
        "gdtools_status": gdtools_status,
        "prokadiff_commit": prok_commit,
        # FIX-012: gdtools SUBTRACT counts are the primary strict parity gate
        "breseq_mutations": len(breseq_records),
        "prokadiff_mutations": len(prok_records),
        "gdtools_strict_over": gdtools_over_count,
        "gdtools_strict_under": gdtools_under_count,
        # FIX-013: normalized parity (right-align stub — will improve with real normalization)
        "normalized_over": len(comp["normalized_over"]),
        "normalized_under": len(comp["normalized_under"]),
        # Diagnostic: Python exact comparison (not the primary gate)
        "record_exact_over": len(comp["record_exact_over"]),
        "record_exact_under": len(comp["record_exact_under"]),
        # Per-type breakdown (based on record_exact for diagnostics)
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
        "con_over": comp["over_by_type"].get("CON", 0),
        "con_under": comp["under_by_type"].get("CON", 0),
        "inv_over": comp["over_by_type"].get("INV", 0),
        "inv_under": comp["under_by_type"].get("INV", 0),
        "jc_record_over": comp["over_by_type"].get("JC", 0),
        "jc_record_under": comp["under_by_type"].get("JC", 0),
        "jc_tolerant_match": comp["matched_jc_tol"],
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
            "gdtools_status": gdtools_status,
            "matched_jc_tol": comp["matched_jc_tol"],
            "unmatched_over_jc_tol": comp["unmatched_over_jc_tol"],
            "unmatched_under_jc_tol": comp["unmatched_under_jc_tol"],
        }, f, indent=2)

    print(f"[compare_gd] Parity evaluation for dataset '{args.dataset}':")
    print(f"  breseq mutations: {len(breseq_records)}, prokadiff mutations: {len(prok_records)}")
    print(f"  gdtools_status: {gdtools_status}")
    if gdtools_used:
        print(f"  gdtools_strict_over: {gdtools_over_count}, gdtools_strict_under: {gdtools_under_count}  [PRIMARY GATE]")
    else:
        print(f"  gdtools unavailable — using Python exact comparison (diagnostic only):")
        print(f"  record_exact_over: {len(comp['record_exact_over'])}, record_exact_under: {len(comp['record_exact_under'])}")
        print("  WARNING: without gdtools, strict parity counts are approximate.")
    print(f"  normalized_over: {len(comp['normalized_over'])}, normalized_under: {len(comp['normalized_under'])}  [FIX-013]")
    print(f"  JC matched with tol={args.jc_tol_bp}bp: {comp['matched_jc_tol']}")
    print(f"  Wrote {over_gd_path}, {under_gd_path}, {parity_tsv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
