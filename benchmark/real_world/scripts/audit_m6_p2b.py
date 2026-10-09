#!/usr/bin/env python3
"""Diagnose M6-P2a publication-event discrepancies from retained evidence.

Aligned BAM files were deleted after each evidence job. This audit does not
invent depth or allele fractions for sites that have no GD record.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
AUDIT = REPO / "benchmark/real_world/manifests/prjna1088182_truth_coordinate_audit.tsv"
EVENTS = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p2a_event_recovery.tsv"
RELATIONSHIPS = REPO / "benchmark/real_world/manifests/prjna1088182_relationships.tsv"
OUT = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p2b_event_diagnostics.tsv"

COLUMNS = [
    "sample",
    "run",
    "strain",
    "reference",
    "published_event",
    "publication_source",
    "mixture_status",
    "reference_base",
    "reference_base_matches",
    "p2a_recovery",
    "p2a_diagnostic",
    "gd_record_at_locus",
    "mc_interval_covers_locus",
    "bam_artifact",
    "total_read_depth",
    "alt_supporting_reads",
    "alt_allele_fraction",
    "ref_supporting_reads",
    "base_quality_distribution",
    "mapping_quality_distribution",
    "strand_balance",
    "diagnosis",
    "confidence",
    "observation",
    "interpretation",
]


def mixture_status(mixture_evidence: str) -> str:
    """Classify mixture from the spreadsheet field, not from SampleName."""
    values = []
    for part in mixture_evidence.split("||"):
        for field in part.split("|"):
            field = field.strip()
            if field.startswith("mixed="):
                values.append(field.split("=", 1)[1].strip())
    if any(value and value != "N/A" for value in values):
        return "DOCUMENTED_POOLED"
    if values and all(value == "N/A" for value in values):
        return "UNRESOLVED_MIXTURE"
    return "UNRESOLVED_MIXTURE"


def load_sequence(accession: str) -> str:
    if accession == "NC_000913.3":
        path = REPO / "data/references/MG1655.gbk"
        chunks = []
        origin = False
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("ORIGIN"):
                origin = True
                continue
            if not origin:
                continue
            if line.startswith("//"):
                break
            chunks.append("".join(ch for ch in line if ch.isalpha()))
        return "".join(chunks).upper()
    if accession == "CP009273.1":
        path = REPO / "data/references/m6_identity/CP009273.1.fasta"
        return "".join(
            line.strip()
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
            if line and not line.startswith(">")
        ).upper()
    raise ValueError(accession)


def parse_gd(path: Path) -> list[list[str]]:
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        rows.append(line.split("\t"))
    return rows


def locus_records(records: list[list[str]], accession: str, position: int) -> tuple[str, str]:
    exact = []
    covered = False
    for parts in records:
        if len(parts) < 6 or parts[3] != accession:
            continue
        if parts[0] == "MC":
            start, end = int(parts[4]), int(parts[5])
            if start <= position <= end:
                covered = True
            continue
        if parts[0] in {"SNP", "INS", "DEL", "SUB"} and int(parts[4]) == position:
            exact.append(parts[0] + ":" + parts[1])
    return (",".join(exact) if exact else ".", "true" if covered else "false")


def unique_audit_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    seen = set()
    kept = []
    for row in rows:
        if row["denominator"] == "exclude" or row["tier_a_eligibility"] != "ELIGIBLE_EXACT_RECOVERY":
            if row["sample_name"] == "lambdaRed_3_4" and row["alt_allele"] == "+AAT":
                key = (row["sample_name"], row["reported_position"], row["ref_allele"], row["alt_allele"])
                if key not in seen:
                    seen.add(key)
                    kept.append(row)
            continue
        key = (row["sample_name"], row["reported_position"], row["ref_allele"], row["alt_allele"])
        if key in seen:
            continue
        seen.add(key)
        kept.append(row)
    return kept


def build_rows() -> list[dict[str, str]]:
    audit = list(csv.DictReader(AUDIT.open(encoding="utf-8"), delimiter="\t"))
    events = list(csv.DictReader(EVENTS.open(encoding="utf-8"), delimiter="\t"))
    relationships = {
        row["sample_name"]: row
        for row in csv.DictReader(RELATIONSHIPS.open(encoding="utf-8"), delimiter="\t")
    }
    by_event = {
        (row["sample"], row["publication_event"]): row
        for row in events
        if row["eligible_denominator_flag"] == "include" or row["publication_event"].endswith("+AAT") or "+AAT" in row["publication_event"]
    }
    sequences = {
        "NC_000913.3": load_sequence("NC_000913.3"),
        "CP009273.1": load_sequence("CP009273.1"),
    }
    gd_cache: dict[str, list[list[str]]] = {}
    built = []
    for row in unique_audit_rows(audit):
        position = int(row["reported_position"])
        accession = row["execution_reference"]
        label = f"{row['reported_position']} {row['ref_allele']}>{row['alt_allele']}"
        event = by_event[(row["sample_name"], label)]
        gd_path = REPO / event["evidence_location"]
        if event["evidence_location"] not in gd_cache:
            gd_cache[event["evidence_location"]] = parse_gd(gd_path)
        exact, mc_cover = locus_records(gd_cache[event["evidence_location"]], accession, position)
        base = "."
        matches = "not_applicable"
        if row["ref_allele"] not in {"", "."} and not row["alt_allele"].startswith("+"):
            base = sequences[accession][position - 1]
            matches = "true" if base == row["ref_allele"] else "false"
        unrecovered = event["recovery_classification"] == "NOT_RECOVERED"
        excluded = row["denominator"] == "exclude"
        if excluded:
            diagnosis = "NOT_EVALUABLE"
            confidence = "HIGH"
            observation = "Frozen coordinate audit keeps 3112908 +AAT COORDINATE_UNRESOLVED. It is excluded from the recovery denominator."
            interpretation = "No coordinate was assigned, so absence of an INS is not a negative call."
            depth = alt_n = alt_f = ref_n = bq = mq = strand = "NOT_APPLICABLE"
        elif unrecovered:
            diagnosis = "UNRESOLVED"
            confidence = "LOW"
            observation = (
                "No SNP, SUB, INS, or DEL GD record is at this coordinate on the execution reference. "
                "No MC interval covers the coordinate. aligned.bam was removed because keep_bam is false."
            )
            interpretation = (
                "The missing call is observed. Read depth, allele fraction, base quality, and mapping quality were not measured. "
                "A pooled component below the RA frequency 0.75 rule is a hypothesis, not a demonstrated filter failure."
            )
            depth = alt_n = alt_f = ref_n = bq = mq = strand = "NOT_MEASURED"
        else:
            diagnosis = event["recovery_classification"]
            confidence = "HIGH"
            observation = f"GD record {exact} matches the publication position and alternate allele. DV1 EventId remains NOT_EVALUABLE."
            interpretation = "Exact recovery does not establish editing causality."
            depth = alt_n = alt_f = ref_n = "SEE_GD_RA_FIELDS"
            bq = mq = strand = "NOT_IN_GD"
        built.append(
            {
                "sample": row["sample_name"],
                "run": event["run"],
                "strain": relationships[row["sample_name"]]["sra_organism"],
                "reference": accession,
                "published_event": label,
                "publication_source": "PUBLICATION_REPORTED",
                "mixture_status": mixture_status(relationships[row["sample_name"]]["mixture_evidence"]),
                "reference_base": base,
                "reference_base_matches": matches,
                "p2a_recovery": event["recovery_classification"],
                "p2a_diagnostic": event["diagnostic_reason"],
                "gd_record_at_locus": exact,
                "mc_interval_covers_locus": mc_cover,
                "bam_artifact": "ABSENT",
                "total_read_depth": depth,
                "alt_supporting_reads": alt_n,
                "alt_allele_fraction": alt_f,
                "ref_supporting_reads": ref_n,
                "base_quality_distribution": bq,
                "mapping_quality_distribution": mq,
                "strand_balance": strand,
                "diagnosis": diagnosis,
                "confidence": confidence,
                "observation": observation,
                "interpretation": interpretation,
            }
        )
    return built


def main() -> int:
    parser = argparse.ArgumentParser(description="Write the M6-P2b per-event diagnostic table.")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    rows = build_rows()
    if not args.write:
        print(f"rows={len(rows)} unresolved={sum(row['diagnosis']=='UNRESOLVED' for row in rows)}")
        return 0
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUT} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
