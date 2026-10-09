#!/usr/bin/env python3
"""Freeze the M6-P3a targeted diagnostic inventory from the P2b table.

This does not create publication truth and does not measure reads.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
P2B = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p2b_event_diagnostics.tsv"
OUT = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3a_targets.tsv"

# First execution subset. Identities are checked against the P2b table.
PRIORITY_RUNS = ("SRR28349204", "SRR28349176", "SRR28349230")

COLUMNS = [
    "sample",
    "run",
    "reference",
    "coordinate",
    "ref",
    "alt",
    "original_recovery_status",
    "p2a_diagnostic",
    "p2b_diagnosis",
    "bam_artifact",
    "depth_status",
    "diagnostic_priority",
]


def parse_event(label: str) -> tuple[str, str, str]:
    position, alleles = label.split()
    ref, alt = alleles.split(">")
    return position, ref, alt


def build_targets(rows: list[dict[str, str]] | None = None) -> list[dict[str, str]]:
    if rows is None:
        with P2B.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    unresolved = [row for row in rows if row["diagnosis"] == "UNRESOLVED"]
    keys = [(row["sample"], row["run"], row["published_event"]) for row in unresolved]
    if len(keys) != 10 or len(set(keys)) != 10:
        raise ValueError(f"expected 10 unique unresolved alleles, found {len(keys)}")
    targets = []
    for row in unresolved:
        if row["p2a_recovery"] != "NOT_RECOVERED" or row["p2a_diagnostic"] != "UNKNOWN":
            raise ValueError(f"{row['published_event']} lost its P2a classification")
        if row["total_read_depth"] != "NOT_MEASURED" or row["bam_artifact"] != "ABSENT":
            raise ValueError(f"{row['published_event']} no longer records missing BAM evidence")
        if row["alt_allele_fraction"] in {"0", "0.0"}:
            raise ValueError("NOT_MEASURED was replaced with zero")
        position, ref, alt = parse_event(row["published_event"])
        if row["reference_base"] != ref or row["reference_base_matches"] != "true":
            raise ValueError(f"{row['published_event']} reference base disagrees with P2b")
        targets.append(
            {
                "sample": row["sample"],
                "run": row["run"],
                "reference": row["reference"],
                "coordinate": position,
                "ref": ref,
                "alt": alt,
                "original_recovery_status": row["p2a_recovery"],
                "p2a_diagnostic": row["p2a_diagnostic"],
                "p2b_diagnosis": row["diagnosis"],
                "bam_artifact": row["bam_artifact"],
                "depth_status": row["total_read_depth"],
                "diagnostic_priority": "1" if row["run"] in PRIORITY_RUNS else "2",
            }
        )
    priority = {row["run"] for row in targets if row["diagnostic_priority"] == "1"}
    if priority != set(PRIORITY_RUNS):
        raise ValueError(f"priority runs {priority} do not match the planned subset")
    return targets


def classify_diagnostic(measurement: dict[str, int | str]) -> str:
    """Classify one locus from measured counts. Mixture percentages are not inputs."""
    if measurement["reference_base"] != measurement["published_ref"]:
        return "REFERENCE_DISAGREEMENT"
    usable = int(measurement["caller_shaped_depth"])
    alt = int(measurement["caller_shaped_alt"])
    raw_alt = int(measurement["raw_alt"])
    q20_alt = int(measurement["q20_alt"])
    low_mapq_alt = int(measurement["low_mapq_alt"])
    if usable < 5 and raw_alt == 0:
        return "INSUFFICIENT_COVERAGE"
    if raw_alt > 0 and q20_alt == 0 and int(measurement["q20_depth"]) >= 5:
        return "LOW_QUALITY_SUPPORT"
    if alt == 0 and low_mapq_alt > 0 and usable < 5:
        return "MAPPING_AMBIGUITY"
    if usable < 5:
        return "UNRESOLVED"
    if usable < 20 and alt == 0:
        return "UNRESOLVED"
    fraction = alt / usable
    if alt == 0 and usable >= 20:
        return "ALT_ABSENT_WITH_ADEQUATE_DEPTH"
    if fraction < 0.75:
        return "LOW_ALT_FREQUENCY"
    if fraction >= 0.75:
        return "ALT_PRESENT_ABOVE_CONSENSUS_THRESHOLD"
    return "UNRESOLVED"


def main() -> int:
    parser = argparse.ArgumentParser(description="Write the M6-P3a target manifest.")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    rows = build_targets()
    if not args.write:
        print(f"targets={len(rows)} priority1={sum(row['diagnostic_priority']=='1' for row in rows)}")
        return 0
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUT} rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
