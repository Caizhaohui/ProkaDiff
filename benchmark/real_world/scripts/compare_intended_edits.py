#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path
from typing import Final

from validation_context import peer_complete_has_differential_support


STATUSES: Final = {"COMPLETE", "PARTIAL", "MISSING", "UNEXPECTED_STRUCTURE"}


class ProductError(Exception):
    pass


def read_edit_outcomes(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ProductError(f"edit_outcomes.tsv not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"edit_id", "seq_id", "expected_start", "expected_end", "status"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ProductError("edit_outcomes.tsv lacks required M4 columns")
        rows = [dict(row) for row in reader]
    for row in rows:
        row["status"] = row["status"].upper()
        if row["status"] not in STATUSES:
            raise ProductError(f"unknown intended edit status {row['status']!r} for {row['edit_id']}")
    if not rows:
        raise ProductError("edit_outcomes.tsv has no intended edit rows")
    return rows


def evaluate_outcomes(rows: list[dict[str, str]], expected: str) -> dict[str, int]:
    metrics = {
        "complete": 0,
        "partial": 0,
        "missing": 0,
        "unexpected_structure": 0,
        "false_complete": 0,
        "status_mismatch": 0,
        "total_evaluated": len(rows),
    }
    for row in rows:
        observed = row["status"]
        metrics[observed.lower()] += 1
        if observed != expected:
            metrics["status_mismatch"] += 1
        if expected != "COMPLETE" and observed == "COMPLETE":
            metrics["false_complete"] += 1
    return metrics


def evaluate_peer_outcomes(rows: list[dict[str, str]]) -> dict[str, int]:
    metrics = evaluate_outcomes(rows, rows[0]["status"])
    metrics["status_mismatch"] = 0
    metrics["false_complete"] = sum(
        row["status"] == "COMPLETE" and not peer_complete_has_differential_support(row)
        for row in rows
    )
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description="Check M4 intended-edit rows. Peer comparisons do not preset a status.")
    parser.add_argument("--edit-outcomes", required=True, type=Path, help="Current product edit_outcomes.tsv")
    parser.add_argument("--expected-status", choices=sorted(STATUSES))
    parser.add_argument("--comparison-role", choices=("REFERENCE_ONLY", "MATCHED_PARENT", "PEER_COMPARATOR"))
    parser.add_argument("--tsv", action="store_true")
    args = parser.parse_args()
    try:
        rows = read_edit_outcomes(args.edit_outcomes)
    except (ProductError, OSError) as error:
        print(f"intended-edit validation failed: {error}", file=sys.stderr)
        return 2
    if args.comparison_role == "PEER_COMPARATOR":
        if args.expected_status:
            parser.error("PEER_COMPARATOR decides status from differential evidence and does not accept a preset status")
        try:
            metrics = evaluate_peer_outcomes(rows)
        except ValueError as error:
            print(f"intended-edit validation failed: {error}", file=sys.stderr)
            return 2
        gate_passed = metrics["false_complete"] == 0 and metrics["status_mismatch"] == 0
    elif args.expected_status:
        metrics = evaluate_outcomes(rows, args.expected_status)
        gate_passed = metrics["false_complete"] == 0 and metrics["status_mismatch"] == 0
    else:
        parser.error("--expected-status is required unless --comparison-role PEER_COMPARATOR")
    if args.tsv:
        print("metric\tvalue")
        for key, value in metrics.items():
            print(f"{key}\t{value}")
        print(f"gate_passed\t{str(gate_passed).lower()}")
    else:
        for key, value in metrics.items():
            print(f"{key}\t{value}")
        print(f"gate_passed\t{str(gate_passed).lower()}")
    return 0 if gate_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
