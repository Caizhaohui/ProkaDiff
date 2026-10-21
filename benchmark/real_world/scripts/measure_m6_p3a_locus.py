#!/usr/bin/env python3
"""Count one publication coordinate on stage-1 SAM/BAM text.

Primary alignments only. Secondary and supplementary records are ignored.
Phred >= 20 is a diagnostic filter, not an RA rule. End-trim of one aligned
base approximates production END_TRIM and is not a complete RA reproduction.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _tags(fields: list[str]) -> dict[str, str]:
    out = {}
    for field in fields[11:]:
        if ":" in field:
            key, _typ, value = field.split(":", 2)
            out[key] = value
    return out


def adjusted_mapq(mapq: int, tags: dict[str, str]) -> int:
    if "AS" in tags and "XS" in tags:
        try:
            if int(tags["XS"]) >= int(tags["AS"]):
                return 0
        except ValueError:
            return mapq
    return mapq


def _walk(pos: int, cigar: str, seq: str, qual: str, target: int):
    """Return coverage facts for one primary alignment."""
    ref = pos
    qpos = 0
    hits: list[tuple[int, str | None, int | None]] = []
    covers = False
    base = None
    phred = None
    nearby = False
    for length, op in CIGAR_RE.findall(cigar):
        length = int(length)
        if op in {"M", "=", "X"}:
            for _ in range(length):
                called = seq[qpos].upper() if qpos < len(seq) else "N"
                score = (ord(qual[qpos]) - 33) if qpos < len(qual) and qual not in {"*", ""} else None
                hits.append((ref, called, score))
                if ref == target:
                    covers = True
                    base = called
                    phred = score
                ref += 1
                qpos += 1
        elif op in {"D", "N"}:
            for _ in range(length):
                if op == "D" and length <= 2:
                    hits.append((ref, None, None))
                if ref == target:
                    covers = True
                    base = None
                    phred = None
                if abs(ref - target) <= 5:
                    nearby = True
                ref += 1
            if op == "D" and length > 2 and abs(ref - length - target) <= 5:
                nearby = True
        elif op == "I":
            if abs(ref - target) <= 5:
                nearby = True
            qpos += length
        elif op == "S":
            boundary = ref
            if abs(boundary - target) <= 5:
                nearby = True
            qpos += length
        elif op == "H":
            continue
        else:
            qpos += length
    if len(hits) > 2:
        kept = {item[0] for item in hits[1:-1]}
    else:
        kept = set()
    return covers, base, phred, nearby, target in kept


def empty_counts() -> dict[str, object]:
    return {
        "raw_depth": 0,
        "raw_ref": 0,
        "raw_alt": 0,
        "raw_other": 0,
        "mapq10_depth": 0,
        "mapq10_ref": 0,
        "mapq10_alt": 0,
        "mapq10_other": 0,
        "caller_shaped_depth": 0,
        "caller_shaped_ref": 0,
        "caller_shaped_alt": 0,
        "caller_shaped_other": 0,
        "q20_depth": 0,
        "q20_ref": 0,
        "q20_alt": 0,
        "q20_other": 0,
        "low_mapq_alt": 0,
        "plus_caller_shaped": 0,
        "minus_caller_shaped": 0,
        "plus_caller_shaped_alt": 0,
        "minus_caller_shaped_alt": 0,
        "nearby_cigar_reads": 0,
        "duplicate_primary_skipped": 0,
        "base_quality_hist": Counter(),
        "mapping_quality_hist": Counter(),
    }


def _add(bucket: str, counts: dict[str, object], base: str | None, ref: str, alt: str) -> None:
    counts[f"{bucket}_depth"] = int(counts[f"{bucket}_depth"]) + 1
    if base == ref:
        counts[f"{bucket}_ref"] = int(counts[f"{bucket}_ref"]) + 1
    elif base == alt:
        counts[f"{bucket}_alt"] = int(counts[f"{bucket}_alt"]) + 1
    else:
        counts[f"{bucket}_other"] = int(counts[f"{bucket}_other"]) + 1


def analyze_sam(lines: list[str], position: int, published_ref: str, published_alt: str) -> dict[str, object]:
    counts = empty_counts()
    seen: set[str] = set()
    for line in lines:
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        qname, flag_s, _rname, pos_s, mapq_s, cigar, _rnext, _pnext, _tlen, seq, qual = fields[:11]
        flag = int(flag_s)
        if flag & 0x4 or flag & 0x100 or flag & 0x800:
            continue
        mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0
        identity = (qname, mate)
        if identity in seen:
            counts["duplicate_primary_skipped"] = int(counts["duplicate_primary_skipped"]) + 1
            continue
        seen.add(identity)
        mapq = int(mapq_s)
        tags = _tags(fields)
        usable_mapq = adjusted_mapq(mapq, tags)
        covers, base, phred, nearby, in_trimmed = _walk(int(pos_s), cigar, seq, qual, position)
        if not covers:
            continue
        if nearby:
            counts["nearby_cigar_reads"] = int(counts["nearby_cigar_reads"]) + 1
        hist_q = counts["mapping_quality_hist"]
        assert isinstance(hist_q, Counter)
        hist_q[mapq] += 1
        if phred is not None:
            hist_b = counts["base_quality_hist"]
            assert isinstance(hist_b, Counter)
            hist_b[phred] += 1
        _add("raw", counts, base, published_ref, published_alt)
        if usable_mapq >= 10:
            _add("mapq10", counts, base, published_ref, published_alt)
        else:
            if base == published_alt:
                counts["low_mapq_alt"] = int(counts["low_mapq_alt"]) + 1
        if usable_mapq >= 10 and in_trimmed:
            _add("caller_shaped", counts, base, published_ref, published_alt)
            if flag & 0x10:
                counts["minus_caller_shaped"] = int(counts["minus_caller_shaped"]) + 1
                if base == published_alt:
                    counts["minus_caller_shaped_alt"] = int(counts["minus_caller_shaped_alt"]) + 1
            else:
                counts["plus_caller_shaped"] = int(counts["plus_caller_shaped"]) + 1
                if base == published_alt:
                    counts["plus_caller_shaped_alt"] = int(counts["plus_caller_shaped_alt"]) + 1
        if phred is not None and phred >= 20:
            _add("q20", counts, base, published_ref, published_alt)
    return counts


def fraction(numer: int, denom: int) -> str:
    if denom == 0:
        return "NA"
    return f"{numer / denom:.6f}"


def hist_text(counter: Counter) -> str:
    if not counter:
        return "."
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main() -> int:
    parser = argparse.ArgumentParser(description="Measure one locus from samtools view output.")
    parser.add_argument("--sam", required=True)
    parser.add_argument("--position", type=int, required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--alt", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    lines = open(args.sam, encoding="utf-8", errors="replace").read().splitlines()
    counts = analyze_sam(lines, args.position, args.ref, args.alt)
    row = {
        "sample": args.sample,
        "run": args.run,
        "reference": args.reference,
        "coordinate": str(args.position),
        "published_ref": args.ref,
        "published_alt": args.alt,
        "raw_depth": counts["raw_depth"],
        "raw_ref": counts["raw_ref"],
        "raw_alt": counts["raw_alt"],
        "raw_other": counts["raw_other"],
        "mapq10_depth": counts["mapq10_depth"],
        "mapq10_ref": counts["mapq10_ref"],
        "mapq10_alt": counts["mapq10_alt"],
        "mapq10_other": counts["mapq10_other"],
        "caller_shaped_depth": counts["caller_shaped_depth"],
        "caller_shaped_ref": counts["caller_shaped_ref"],
        "caller_shaped_alt": counts["caller_shaped_alt"],
        "caller_shaped_other": counts["caller_shaped_other"],
        "caller_shaped_alt_fraction": fraction(int(counts["caller_shaped_alt"]), int(counts["caller_shaped_depth"])),
        "q20_depth": counts["q20_depth"],
        "q20_ref": counts["q20_ref"],
        "q20_alt": counts["q20_alt"],
        "q20_other": counts["q20_other"],
        "low_mapq_alt": counts["low_mapq_alt"],
        "plus_caller_shaped": counts["plus_caller_shaped"],
        "minus_caller_shaped": counts["minus_caller_shaped"],
        "plus_caller_shaped_alt": counts["plus_caller_shaped_alt"],
        "minus_caller_shaped_alt": counts["minus_caller_shaped_alt"],
        "nearby_cigar_reads": counts["nearby_cigar_reads"],
        "duplicate_primary_skipped": counts["duplicate_primary_skipped"],
        "base_quality_hist": hist_text(counts["base_quality_hist"]),
        "mapping_quality_hist": hist_text(counts["mapping_quality_hist"]),
        "count_policy": "primary_only_one_record_per_qname_and_mate_flag",
        "caller_shaped_note": "MAPQ>=10 after XS>=AS, drop one aligned base at each end; not full RA homopolymer masking",
    }
    fields = list(row)
    with open(args.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({key: str(value) for key, value in row.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
