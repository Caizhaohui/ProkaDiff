#!/usr/bin/env python3
"""Layer-2 Genome Diff leftover + wall/RSS CSV (compute node only).

Red-line types: SNP / INS / DEL / MOB / AMP / CON (exact key).
JC: greedy match after canonicalizing sides, ±JC_TOL_BP (default 5).
UN / RA / MC are ignored. Does not invent numbers; missing inputs → exit 1.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import socket
import subprocess
import sys
from pathlib import Path

RED = {"SNP", "INS", "DEL", "MOB", "AMP", "CON"}
JC = "JC"
IGNORE = {"UN", "RA", "MC"}
JC_TOL_BP = 5


def parse_gd_muts(path: Path):
    rows = []
    if not path.is_file():
        return rows
    for line in path.read_text(errors="replace").splitlines():
        if not line or line.startswith("#") or line.startswith("="):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        kind = parts[0]
        if kind in IGNORE:
            continue
        if kind not in RED and kind != JC:
            continue
        fields = parts[3:]
        rows.append({"kind": kind, "fields": fields, "line": line})
    return rows


def exact_key(row):
    return (row["kind"], tuple(row["fields"]))


def leftover_exact(left, right):
    keys = {exact_key(r) for r in right}
    return [r for r in left if exact_key(r) not in keys]


def jc_canonical(fields):
    if len(fields) < 6:
        return None
    try:
        a = (fields[0], fields[2], int(fields[1]))
        b = (fields[3], fields[5], int(fields[4]))
    except (ValueError, IndexError):
        return None
    if (a[0], a[1], a[2]) > (b[0], b[1], b[2]):
        a, b = b, a
    return a, b


def jc_close(a, b, tol):
    ca = jc_canonical(a["fields"])
    cb = jc_canonical(b["fields"])
    if ca is None or cb is None:
        return False
    (a1, a2), (b1, b2) = ca, cb
    if (a1[0], a1[1], a2[0], a2[1]) != (b1[0], b1[1], b2[0], b2[1]):
        return False
    return abs(a1[2] - b1[2]) <= tol and abs(a2[2] - b2[2]) <= tol


def pair_jc(over, under, tol):
    over_jc = [r for r in over if r["kind"] == JC]
    under_jc = [r for r in under if r["kind"] == JC]
    used = set()
    matched = 0
    for j1 in over_jc:
        for i, j2 in enumerate(under_jc):
            if i in used:
                continue
            if jc_close(j1, j2, tol):
                used.add(i)
                matched += 1
                break
    return (
        matched,
        len(over_jc) - matched,
        len(under_jc) - len(used),
    )


def compare(left, right, tol):
    over = leftover_exact(left, right)
    under = leftover_exact(right, left)
    over_red = [r for r in over if r["kind"] in RED]
    under_red = [r for r in under if r["kind"] in RED]
    jc_matched, over_jc, under_jc = pair_jc(over, under, tol)
    return {
        "over_red": over_red,
        "under_red": under_red,
        "jc_matched": jc_matched,
        "over_jc_unmatched": over_jc,
        "under_jc_unmatched": under_jc,
        "over_red_n": len(over_red),
        "under_red_n": len(under_red),
    }


def fmt_preview(rows, n=6):
    keys = []
    for r in rows[:n]:
        keys.append((r["kind"], tuple(r["fields"][:6])))
    return repr(keys)


def parse_time(path: Path):
    wall_s = ""
    rss_kb = ""
    if not path.is_file():
        return wall_s, rss_kb
    text = path.read_text(errors="replace")
    for line in text.splitlines():
        if "Elapsed (wall clock) time" in line:
            m = re.search(r":\s*(\d+):(\d+):(\d+(?:\.\d+)?)$", line)
            if m:
                wall_s = str(int(m.group(1)) * 3600 + int(m.group(2)) * 60 + float(m.group(3)))
            else:
                m = re.search(r":\s*(\d+):(\d+(?:\.\d+)?)$", line)
                if m:
                    wall_s = str(int(m.group(1)) * 60 + float(m.group(2)))
        if "Maximum resident set size" in line:
            m = re.search(r":\s*(\d+)", line)
            if m:
                rss_kb = m.group(1)
    return wall_s, rss_kb


def ver(cmd):
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True, timeout=30)
        return out.strip().splitlines()[0][:120]
    except Exception as e:
        return f"unavailable:{e}"


def summarize(tag, stats):
    return (
        f"{tag}_over_red={stats['over_red_n']};{tag}_under_red={stats['under_red_n']};"
        f"{tag}_over_jc_unmatched={stats['over_jc_unmatched']};"
        f"{tag}_under_jc_unmatched={stats['under_jc_unmatched']};"
        f"{tag}_jc_matched_tol={stats['jc_matched']};"
        f"{tag}_over_red_preview={fmt_preview(stats['over_red'])};"
        f"{tag}_under_red_preview={fmt_preview(stats['under_red'])}"
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--jobout", required=True, type=Path)
    p.add_argument("--workload", required=True)
    p.add_argument("--threads", required=True)
    p.add_argument("--prokdiff-gd", required=True, type=Path)
    p.add_argument("--breseq-gd", required=True, type=Path)
    p.add_argument("--oracle-gd", type=Path, default=None)
    p.add_argument("--csv", type=Path, default=None)
    p.add_argument("--jc-tol-bp", type=int, default=JC_TOL_BP)
    args = p.parse_args()

    rust_gd = args.prokdiff_gd
    breseq_gd = args.breseq_gd
    if not rust_gd.is_file():
        print(f"error: missing ProkDiff GD {rust_gd}", file=sys.stderr)
        return 1
    if not breseq_gd.is_file():
        print(f"error: missing breseq GD {breseq_gd}", file=sys.stderr)
        return 1

    rust = parse_gd_muts(rust_gd)
    breseq = parse_gd_muts(breseq_gd)
    vs_node = compare(rust, breseq, args.jc_tol_bp)
    notes = [
        f"jc_tol_bp={args.jc_tol_bp}",
        "red=SNP,INS,DEL,MOB,AMP,CON exact; JC greedy ±tol after side canonicalization; UN/RA/MC ignored",
        "oracle_breseq_env=0.39.0_not_0.40.x",
        summarize("vs_samenode", vs_node),
    ]
    red_fail = vs_node["over_red_n"] or vs_node["under_red_n"]

    if args.oracle_gd is not None:
        if not args.oracle_gd.is_file():
            print(f"error: missing oracle GD {args.oracle_gd}", file=sys.stderr)
            return 1
        oracle = parse_gd_muts(args.oracle_gd)
        vs_off = compare(rust, oracle, args.jc_tol_bp)
        notes.append(summarize("vs_official", vs_off))
        red_fail = red_fail or vs_off["over_red_n"] or vs_off["under_red_n"]

    notes_s = ";".join(notes)
    csv_path = args.csv or (args.jobout.parent / f"{args.workload}.csv")
    node = socket.gethostname()
    breseq_v = ver(["breseq", "--version"])
    bt2_v = ver(["bowtie2", "--version"])
    new_file = not csv_path.exists()
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("a", newline="") as fh:
        w = csv.writer(fh)
        if new_file:
            w.writerow(
                [
                    "workload",
                    "tool",
                    "threads",
                    "wall_s",
                    "max_rss_kb",
                    "node",
                    "job_id",
                    "breseq_version",
                    "bowtie2_version",
                    "notes",
                ]
            )
        for tool, tfile in (
            ("prokdiff", args.jobout / "prokdiff.time"),
            ("breseq", args.jobout / "breseq.time"),
        ):
            wall, rss = parse_time(tfile)
            w.writerow(
                [
                    args.workload,
                    tool,
                    args.threads,
                    wall,
                    rss,
                    node,
                    os.environ.get("SLURM_JOB_ID", ""),
                    breseq_v,
                    bt2_v,
                    notes_s,
                ]
            )

    summary_path = args.jobout / "parity_notes.txt"
    summary_path.write_text(notes_s + "\n")
    print(f"wrote {csv_path}")
    print(notes_s)
    if red_fail:
        print("WARNING: SNP/INS/DEL/MOB/AMP/CON leftover is non-empty after exact match.")
        return 2
    print("Red-line types bidirectional leftover empty (JC ±tol applied separately).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
