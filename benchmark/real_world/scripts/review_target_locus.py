#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path


FIELDS = [
    "sample_id", "target_id", "state", "evidence_level", "evidence_source",
    "observed_records", "exact_breakpoint_supported", "notes",
]
ALLOWED_TYPES = {"MC", "RA", "DEL", "JC", "MOB", "UN", "SNP"}


class ReviewError(Exception):
    pass


@dataclass(frozen=True)
class Record:
    kind: str
    line: str
    fields: tuple[str, ...]


def read_gd(path: Path) -> list[Record]:
    records: list[Record] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = tuple(line.split("\t"))
        if fields[0] not in ALLOWED_TYPES:
            continue
        if len(fields) < 5:
            raise ReviewError(f"malformed {fields[0]} GD record in {path}: {line}")
        records.append(Record(fields[0], line, fields))
    return records


def record_interval(record: Record) -> tuple[str, int, int] | None:
    f = record.fields
    try:
        if record.kind == "MC":
            return f[3], int(f[4]), int(f[5])
        if record.kind == "DEL":
            start, size = int(f[4]), int(f[5])
            return f[3], start, start + max(size, 1)
        if record.kind in {"RA", "SNP", "INS", "MOB"}:
            return f[3], int(f[4]), int(f[4])
        if record.kind == "UN":
            return f[3], int(f[4]), int(f[5])
        if record.kind == "JC":
            seq1, pos1 = f[3], int(f[4])
            seq2, pos2 = f[6], int(f[7])
            if seq1 == seq2:
                return seq1, min(pos1, pos2), max(pos1, pos2)
            return None
    except (IndexError, ValueError) as error:
        raise ReviewError(f"cannot parse coordinates from GD record: {record.line}") from error
    return None


def canonical_seq_id(value: str) -> str:
    return value.rsplit(".", 1)[0] if value.rsplit(".", 1)[-1].isdigit() else value


def relevant(record: Record, seq_id: str, start: int, end: int) -> bool:
    f = record.fields
    if record.kind == "JC":
        if "ignore=CIRCULAR_CHROMOSOME" in record.line:
            return False
        try:
            sides = ((f[3], int(f[4]), int(f[5])), (f[6], int(f[7]), int(f[8])))
        except (IndexError, ValueError) as error:
            raise ReviewError(f"cannot parse JC coordinates: {record.line}") from error
        if any(canonical_seq_id(s) == canonical_seq_id(seq_id) and start <= p <= end for s, p, _ in sides):
            return True
        if canonical_seq_id(sides[0][0]) == canonical_seq_id(sides[1][0]) == canonical_seq_id(seq_id):
            p1, s1 = sides[0][1], sides[0][2]
            p2, s2 = sides[1][1], sides[1][2]
            if p1 > p2:
                (p1, s1), (p2, s2) = (p2, s2), (p1, s1)
            if s1 == -1 and s2 == 1 and p1 < start <= end < p2:
                if (p2 - p1) <= 100_000 and (start - p1) <= 50_000 and (p2 - end) <= 50_000:
                    return True
        return False
    interval = record_interval(record)
    return interval is not None and canonical_seq_id(interval[0]) == canonical_seq_id(seq_id) and interval[1] <= end and interval[2] >= start


def classify(records: list[Record], seq_id: str, start: int, end: int, source: str = "PROKADIFF") -> tuple[str, str, bool, list[Record]]:
    observed = [record for record in records if relevant(record, seq_id, start, end)]
    supported_jc = [record for record in observed if record.kind == "JC" and (
        source != "PROKADIFF" or any(
            item.startswith("pd_support_reads=") and item.partition("=")[2].isdigit() and int(item.partition("=")[2]) > 0
            for item in record.fields[9:]
        )
    )]
    if supported_jc:
        return "EXACT_STRUCTURE_SUPPORTED", "BREAKPOINT", True, observed
    if any(record.kind in {"DEL", "MOB"} for record in observed):
        return "ABNORMAL_STRUCTURE_SUPPORTED", "PRODUCT_MUTATION", False, observed
    if any(record.kind in {"MC", "RA", "SNP", "INS"} for record in observed):
        return "ABNORMAL_TARGET_EVIDENCE", "COVERAGE" if any(r.kind == "MC" for r in observed) else "READ_EVIDENCE", False, observed
    if any(record.kind in {"UN", "JC"} for record in observed):
        return "INSUFFICIENT_EVIDENCE", "UNRESOLVED", False, observed
    return "NO_ABNORMALITY_DETECTED", "NO_QUALIFYING_RECORD", False, observed


def row_for_gd(sample: str, target: str, source: str, gd_path: Path, seq: str, start: int, end: int) -> dict[str, str]:
    records = read_gd(gd_path)
    state, level, exact, observed = classify(records, seq, start, end, source)
    if source == "PROKADIFF" and state == "ABNORMAL_TARGET_EVIDENCE" and any(r.kind == "MC" for r in observed):
        notes = "Target-region MC supports coverage loss; no supported DEL or JC breakpoint record. No DifferentialEvent or EventId is created."
    else:
        notes = "No relevant record in this source call set does not establish biological normality. No DifferentialEvent, EventId, or IntendedEditAssessment is created."
    return {
        "sample_id": sample, "target_id": target, "state": state,
        "evidence_level": level, "evidence_source": source,
        "observed_records": ";".join(record.line for record in observed) or "NONE",
        "exact_breakpoint_supported": str(exact).lower(), "notes": notes,
    }


def curated_row(sample: str, target: str, truth_path: Path, seq: str, start: int, end: int) -> dict[str, str]:
    with truth_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    truth_id = f"TRUTH_{sample.replace('B21_', 'BL21_')}_DEL"
    matched = [row for row in rows if row.get("truth_id") == truth_id and row.get("seq_id") == seq and int(row["start"]) <= end and int(row["end"]) >= start]
    exact = any(re.search(r"JC:[^;]*:\d+reads", row.get("evidence", "")) for row in matched)
    observed = ";".join(f"{r['truth_id']}:{r.get('evidence', '')}" for r in matched) or "NONE"
    state = "EXACT_STRUCTURE_SUPPORTED" if exact else "ABNORMAL_STRUCTURE_SUPPORTED" if matched else "NO_ABNORMALITY_DETECTED"
    return {
        "sample_id": sample, "target_id": target, "state": state,
        "evidence_level": "CURATED_BREAKPOINT" if exact else "CURATED_STRUCTURE" if matched else "NO_CURATED_RECORD",
        "evidence_source": "CURATED_TRUTH", "observed_records": observed,
        "exact_breakpoint_supported": str(exact).lower(),
        "notes": "External curated interpretation; not represented as a ProkaDiff observation.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--target-id", required=True)
    parser.add_argument("--seq-id", required=True)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--prokadiff-gd", type=Path)
    parser.add_argument("--breseq-gd", type=Path)
    parser.add_argument("--curated-truth", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        rows = []
        for source, path in (("PROKADIFF", args.prokadiff_gd), ("BRESEQ", args.breseq_gd)):
            if path:
                rows.append(row_for_gd(args.sample, args.target_id, source, path, args.seq_id, args.start, args.end))
        if args.curated_truth:
            rows.append(curated_row(args.sample, args.target_id, args.curated_truth, args.seq_id, args.start, args.end))
        if not rows:
            raise ReviewError("provide at least one evidence source")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    except (ReviewError, OSError, ValueError) as error:
        print(f"target-locus review failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
