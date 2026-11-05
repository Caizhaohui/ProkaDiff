#!/usr/bin/env python3
"""Reconcile M6-P2a evidence GD files with frozen publication truth.

Evidence GD numeric IDs are record identifiers. They are not DV1 EventIds.
Publication integers match only at the same coordinate and allele. This script
does not change caller thresholds and does not treat NONE as a false-positive set.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass

MUTATION_TYPES = ("SNP", "SUB", "INS", "DEL", "JC", "MOB")
EVIDENCE_TYPES = ("MC",)
RECOVERED = {"RECOVERED_EXACT", "RECOVERED_EQUIVALENT"}
DIAGNOSTIC_REASONS = {
    "NO_READ_COVERAGE",
    "LOW_ALLELE_SUPPORT",
    "QUALITY_FILTER",
    "MAPPING_FILTER",
    "CALLER_THRESHOLD",
    "REPRESENTATION_DIFFERENCE",
    "PUBLICATION_PROVENANCE_ISSUE",
    "UNKNOWN",
    ".",
}


class MalformedGdError(ValueError):
    pass


@dataclass(frozen=True)
class GdRecord:
    record_type: str
    gd_record_id: str
    parent_ids: str
    seq_id: str
    position: int
    new_seq: str
    size: int | None
    raw: str

    @property
    def dv1_event_id(self) -> str:
        return "NOT_EVALUABLE"


def _require_int(value: str, line_no: int, label: str) -> int:
    if not value.isdigit():
        raise MalformedGdError(f"line {line_no}: {label} is not an integer")
    return int(value)


def parse_gd(text: str) -> list[GdRecord]:
    records: list[GdRecord] = []
    for line_no, line in enumerate(text.splitlines(), start=1):
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            raise MalformedGdError(f"line {line_no}: fewer than 5 columns")
        record_type, gd_id, parent_ids, seq_id = parts[0], parts[1], parts[2], parts[3]
        if record_type == "" or seq_id == "":
            raise MalformedGdError(f"line {line_no}: missing type or sequence id")
        if not gd_id.isdigit():
            raise MalformedGdError(f"line {line_no}: GD record id is not an integer")
        position = _require_int(parts[4], line_no, "position")
        new_seq = ""
        size: int | None = None
        if record_type == "SNP":
            if len(parts) < 6 or parts[5] == "":
                raise MalformedGdError(f"line {line_no}: SNP is missing the new base")
            new_seq = parts[5]
        elif record_type == "INS":
            if len(parts) < 6 or parts[5] == "":
                raise MalformedGdError(f"line {line_no}: INS is missing the inserted sequence")
            new_seq = parts[5]
        elif record_type == "SUB":
            if len(parts) < 7:
                raise MalformedGdError(f"line {line_no}: SUB is missing size or new sequence")
            size = _require_int(parts[5], line_no, "SUB size")
            new_seq = parts[6]
        elif record_type == "DEL":
            if len(parts) < 6:
                raise MalformedGdError(f"line {line_no}: DEL is missing size")
            size = _require_int(parts[5], line_no, "DEL size")
        elif record_type == "MC":
            if len(parts) < 6:
                raise MalformedGdError(f"line {line_no}: MC is missing the end coordinate")
            _require_int(parts[5], line_no, "MC end")
        elif record_type == "JC":
            if len(parts) < 10:
                raise MalformedGdError(f"line {line_no}: JC is missing the second side")
            _require_int(parts[7], line_no, "JC side-2 position")
        elif record_type == "MOB":
            if len(parts) < 8:
                raise MalformedGdError(f"line {line_no}: MOB is missing repeat fields")
        elif record_type not in MUTATION_TYPES and record_type not in EVIDENCE_TYPES:
            if not parts[4].isdigit():
                raise MalformedGdError(f"line {line_no}: unrecognized record is not positioned")
        else:
            raise MalformedGdError(f"line {line_no}: unhandled record type {record_type}")
        records.append(
            GdRecord(
                record_type=record_type,
                gd_record_id=gd_id,
                parent_ids=parent_ids,
                seq_id=seq_id,
                position=position,
                new_seq=new_seq,
                size=size,
                raw=line,
            )
        )
    return records


def inventory(records: list[GdRecord]) -> dict[str, int]:
    counts = {name: 0 for name in (*MUTATION_TYPES, *EVIDENCE_TYPES, "other")}
    for record in records:
        if record.record_type in counts:
            counts[record.record_type] += 1
        else:
            counts["other"] += 1
    counts["total_mutation_calls"] = sum(counts[name] for name in MUTATION_TYPES)
    counts["MC_count"] = counts["MC"]
    return counts


def repeat_diagnostic_count(text: str) -> int:
    rows = [line for line in text.splitlines() if line]
    if not rows:
        return 0
    if rows[0].startswith("sample_id\t"):
        rows = rows[1:]
    return len(rows)


def allele_key(sample_name: str, position: str, ref_allele: str, alt_allele: str) -> tuple[str, str, str, str]:
    return (sample_name, position, ref_allele, alt_allele)


def eligible_unique_alleles(audit_rows: list[dict[str, str]]) -> set[tuple[str, str, str, str]]:
    alleles: set[tuple[str, str, str, str]] = set()
    for row in audit_rows:
        if row.get("tier_a_eligibility") == "ELIGIBLE_EXACT_RECOVERY" and row.get("denominator") == "include":
            alleles.add(
                allele_key(
                    row["sample_name"],
                    row["reported_position"],
                    row["ref_allele"],
                    row["alt_allele"],
                )
            )
    return alleles


def _describe(record: GdRecord) -> str:
    if record.record_type == "SNP":
        return f"SNP {record.seq_id} {record.position} {record.new_seq} gd_record_id={record.gd_record_id}"
    if record.record_type in {"INS", "SUB"}:
        return f"{record.record_type} {record.seq_id} {record.position} {record.new_seq} gd_record_id={record.gd_record_id}"
    return f"{record.record_type} {record.seq_id} {record.position} gd_record_id={record.gd_record_id}"


def classify_publication_event(
    *,
    coordinate_status: str,
    denominator: str,
    position: int,
    ref_allele: str,
    alt_allele: str,
    execution_reference: str,
    reference_base: str | None,
    records: list[GdRecord],
) -> dict[str, str]:
    """Classify one publication row. Unresolved coordinates stay out of the denominator."""
    blank = {
        "recovery": "NOT_EVALUABLE",
        "diagnostic_reason": ".",
        "observed_event": "NA",
        "gd_record_id": "NA",
        "dv1_event_id": "NOT_EVALUABLE",
        "eligible_denominator": "exclude",
        "m4_products": "NOT_EVALUABLE_BY_TIER_A",
        "review_stop": ".",
    }
    if (
        coordinate_status != "COORDINATE_VERIFIED"
        or denominator != "include"
        or ref_allele in {"", "."}
        or alt_allele in {"", "."}
        or alt_allele.startswith("+")
    ):
        return blank
    if reference_base is None or reference_base != ref_allele:
        return {
            "recovery": "NOT_RECOVERED",
            "diagnostic_reason": "PUBLICATION_PROVENANCE_ISSUE",
            "observed_event": "NA",
            "gd_record_id": "NA",
            "dv1_event_id": "NOT_EVALUABLE",
            "eligible_denominator": "include",
            "m4_products": "NOT_EVALUABLE_BY_TIER_A",
            "review_stop": "execution_reference_base_disagrees_with_publication_allele",
        }

    def on_site(record: GdRecord) -> bool:
        return record.seq_id == execution_reference and record.position == position

    exact = [
        record
        for record in records
        if record.record_type == "SNP" and on_site(record) and record.new_seq == alt_allele
    ]
    if exact:
        return {
            "recovery": "RECOVERED_EXACT",
            "diagnostic_reason": ".",
            "observed_event": _describe(exact[0]),
            "gd_record_id": exact[0].gd_record_id,
            "dv1_event_id": "NOT_EVALUABLE",
            "eligible_denominator": "include",
            "m4_products": "NOT_EVALUABLE_BY_TIER_A",
            "review_stop": ".",
        }
    equivalent = [
        record
        for record in records
        if record.record_type == "SUB"
        and on_site(record)
        and record.size == 1
        and record.new_seq == alt_allele
    ]
    if equivalent:
        return {
            "recovery": "RECOVERED_EQUIVALENT",
            "diagnostic_reason": ".",
            "observed_event": _describe(equivalent[0]),
            "gd_record_id": equivalent[0].gd_record_id,
            "dv1_event_id": "NOT_EVALUABLE",
            "eligible_denominator": "include",
            "m4_products": "NOT_EVALUABLE_BY_TIER_A",
            "review_stop": ".",
        }
    representation = [
        record
        for record in records
        if record.record_type in {"INS", "DEL"} and on_site(record)
    ]
    if representation:
        return {
            "recovery": "REPRESENTATION_DIFFERENCE",
            "diagnostic_reason": "REPRESENTATION_DIFFERENCE",
            "observed_event": _describe(representation[0]),
            "gd_record_id": representation[0].gd_record_id,
            "dv1_event_id": "NOT_EVALUABLE",
            "eligible_denominator": "include",
            "m4_products": "NOT_EVALUABLE_BY_TIER_A",
            "review_stop": ".",
        }
    conflict = [
        record
        for record in records
        if record.record_type in {"SNP", "SUB"} and on_site(record)
    ]
    return {
        "recovery": "NOT_RECOVERED",
        "diagnostic_reason": "UNKNOWN",
        "observed_event": _describe(conflict[0]) if conflict else "NA",
        "gd_record_id": conflict[0].gd_record_id if conflict else "NA",
        "dv1_event_id": "NOT_EVALUABLE",
        "eligible_denominator": "include",
        "m4_products": "NOT_EVALUABLE_BY_TIER_A",
        "review_stop": ".",
    }


def deduplicated_evaluations(
    audit_rows: list[dict[str, str]],
    records: list[GdRecord],
    execution_reference: str,
    reference_base_at: dict[int, str],
) -> list[dict[str, str]]:
    """One evaluation per unique publication allele. Duplicate rows are not extra truth."""
    seen: set[tuple[str, str, str, str]] = set()
    evaluations: list[dict[str, str]] = []
    for row in audit_rows:
        key = allele_key(row["sample_name"], row["reported_position"], row["ref_allele"], row["alt_allele"])
        if key in seen:
            continue
        seen.add(key)
        position_text = row["reported_position"]
        position = int(position_text) if position_text.isdigit() else -1
        classified = classify_publication_event(
            coordinate_status=row["coordinate_status"],
            denominator=row["denominator"],
            position=position,
            ref_allele=row["ref_allele"],
            alt_allele=row["alt_allele"],
            execution_reference=execution_reference,
            reference_base=reference_base_at.get(position),
            records=records,
        )
        classified["sample_name"] = row["sample_name"]
        classified["publication_event"] = f"{position_text} {row['ref_allele']}>{row['alt_allele']}"
        classified["duplicate_row"] = row.get("duplicate_row", "false")
        evaluations.append(classified)
    return evaluations


def write_tsv(path: str, rows: list[dict[str, str]], columns: list[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
