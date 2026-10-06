#!/usr/bin/env python3
import argparse
import csv
import re
import sys
from pathlib import Path

from review_target_locus import Record, read_gd, record_interval, relevant
from validate_bl21_registry import read_tsv


FIELDS = [
    "source", "sample", "comparison_role", "record_type", "coordinates",
    "evidence", "event_id_if_applicable", "breakpoint_supported", "interpretation",
]


def supported_prokadiff_jc(record: Record) -> bool:
    return any(
        item.startswith("pd_support_reads=")
        and item.partition("=")[2].isdigit()
        and int(item.partition("=")[2]) > 0
        for item in record.fields[9:]
    )


def record_coordinates(record: Record) -> str:
    f = record.fields
    if record.kind == "JC":
        return f"{f[3]}:{f[4]}/{f[5]}->{f[6]}:{f[7]}/{f[8]}"
    interval = record_interval(record)
    if interval is None:
        return "NA"
    seq_id, start, end = interval
    return f"{seq_id}:{start}" if start == end else f"{seq_id}:{start}-{end}"


def gd_rows(source: str, sample: str, path: Path, role: str, seq: str, start: int, end: int) -> list[dict[str, str]]:
    rows = []
    for record in read_gd(path):
        if not relevant(record, seq, start, end):
            continue
        supported = record.kind == "JC" and (source != "PROKADIFF" or supported_prokadiff_jc(record))
        if record.kind == "MC":
            interpretation = "Target-region coverage observation; not a differential mutation call."
        elif record.kind == "JC" and not supported:
            interpretation = "Junction record lacks source-specific support for a resolved breakpoint."
        elif record.kind == "JC":
            interpretation = "Source-reported junction breakpoint; source identity remains explicit."
        elif record.kind in {"DEL", "MOB"}:
            interpretation = "Source-reported structural product mutation; not independently merged with other sources."
        else:
            interpretation = "Source-specific target-region observation."
        rows.append({
            "source": source, "sample": sample, "comparison_role": role,
            "record_type": record.kind, "coordinates": record_coordinates(record),
            "evidence": record.line, "event_id_if_applicable": "NA",
            "breakpoint_supported": str(supported).lower(), "interpretation": interpretation,
        })
    return rows


def curated_rows(sample: str, truth_path: Path, seq: str, start: int, end: int) -> list[dict[str, str]]:
    _, records = read_tsv(truth_path)
    truth_id = f"TRUTH_{sample.replace('B21_', 'BL21_')}_DEL"
    matched = [
        row for row in records
        if row.get("truth_id") == truth_id and row.get("seq_id") == seq
        and int(row["start"]) <= end and int(row["end"]) >= start
    ]
    result = []
    for row in matched:
        evidence = row.get("evidence", "")
        supported = bool(re.search(r"JC:[^;]*:\d+reads", evidence))
        result.append({
            "source": "CURATED", "sample": sample, "comparison_role": "REFERENCE_ONLY",
            "record_type": row.get("variant_type", "STRUCTURE"),
            "coordinates": f"{row['seq_id']}:{row['start']}-{row['end']}",
            "evidence": evidence, "event_id_if_applicable": "NA",
            "breakpoint_supported": str(supported).lower(),
            "interpretation": "External curated interpretation; not a ProkaDiff or breseq observation.",
        })
    return result


def assessment_rows(sample: str, product_dirs: list[Path]) -> list[dict[str, str]]:
    result = []
    for product_dir in product_dirs:
        guide = product_dir.name
        _, rows = read_tsv(product_dir / "edit_outcomes.tsv")
        for row in rows:
            event_ids = [
                value for field in ("matched_event_ids_dv1", "unexpected_event_ids_dv1")
                for value in row.get(field, "").split(",") if value not in {"", "NA", "NONE"}
            ]
            result.append({
                "source": "PROKADIFF", "sample": sample, "comparison_role": "PEER_COMPARATOR",
                "record_type": "INTENDED_ASSESSMENT",
                "coordinates": f"{row.get('seq_id', 'NA')}:{row.get('expected_start', 'NA')}-{row.get('expected_end', 'NA')}",
                "evidence": f"guide={guide};status={row.get('status', 'NA')}",
                "event_id_if_applicable": ",".join(event_ids) or "NA",
                "breakpoint_supported": "false",
                "interpretation": "M4 differential intended assessment relative to the declared peer comparator.",
            })
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a source-separated target-locus reconciliation TSV.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--seq-id", required=True)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--prokadiff-gd", required=True, type=Path)
    parser.add_argument("--breseq-gd", required=True, type=Path)
    parser.add_argument("--curated-truth", required=True, type=Path)
    parser.add_argument("--products", required=True, nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        rows = gd_rows("PROKADIFF", args.sample, args.prokadiff_gd, "REFERENCE_ONLY", args.seq_id, args.start, args.end)
        rows += gd_rows("BRESEQ", args.sample, args.breseq_gd, "REFERENCE_ONLY", args.seq_id, args.start, args.end)
        rows += curated_rows(args.sample, args.curated_truth, args.seq_id, args.start, args.end)
        rows += assessment_rows(args.sample, args.products)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    except (OSError, ValueError, KeyError) as error:
        print(f"target-locus reconciliation failed: {error}", file=sys.stderr)
        return 1
    print(f"reconciliation_rows\t{len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
