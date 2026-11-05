#!/usr/bin/env python3
"""Freeze the M6-P2a execution cohort from the versioned relationship table.

Eligibility uses each run's own runinfo identity and spots_with_mates.
avgLength is not a selection rule. The two audited SRA1825115 structures are
not copied onto the other 23 runs from that submission.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from m6_reference_policy import (  # noqa: E402
    BW25113_EXECUTION,
    MG1655_EXECUTION,
    RELATIONSHIPS,
    ReferencePolicyError,
    expected_accessions,
    read_tsv,
    require_execution_reference,
)

REPO = Path(__file__).resolve().parents[3]
RUNINFO = REPO / "benchmark/real_world/manifests/prjna1088182_sra_runinfo.csv"
COORDINATES = REPO / "benchmark/real_world/manifests/prjna1088182_truth_coordinate_audit.tsv"
TRUTH = REPO / "benchmark/real_world/truth/prjna1088182_publication_truth.tsv"
OUT = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p2a_execution_cohort.tsv"

SHA256 = {
    MG1655_EXECUTION: "e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0",
    BW25113_EXECUTION: "3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce",
}

SRA_STAT_VERIFIED = {
    "SRR28349244": "TRUE_PAIRED_READS_AVAILABLE",
    "SRR28349243": "TRUE_PAIRED_READS_AVAILABLE",
    "SRR28357593": "SINGLE_BIOLOGICAL_READ_PER_SPOT",
    "SRR28357592": "SINGLE_BIOLOGICAL_READ_PER_SPOT",
}

REUSED_VERIFIED_P1B = ("SRR28349244", "SRR28349243")

COLUMNS = [
    "sample_name",
    "biosample",
    "experiment",
    "submission",
    "run",
    "strain",
    "library_layout_reported",
    "observed_read_structure",
    "execution_reference",
    "reference_sequence_sha256",
    "truth_class",
    "eligible_publication_event_count",
    "execution_eligibility",
    "exclusion_reason",
]


def classify_run(
    *,
    submission: str,
    spots: int,
    spots_with_mates: int,
    library_layout: str,
    identity_ok: bool,
    reference_ok: bool,
    sra_stat_structure: str | None,
) -> tuple[str, str, str]:
    """Return eligibility, observed structure, and exclusion reason.

    avgLength is intentionally absent. A long average read length cannot
    make a run eligible.
    """
    if not identity_ok:
        return "BLOCKED_RUN_IDENTITY", "RUNINFO_IDENTITY_MISMATCH", "SampleName, BioSample, or Experiment disagrees with the relationship table"
    if not reference_ok:
        return "BLOCKED_REFERENCE", "REFERENCE_POLICY_FAILED", "execution reference disagrees with the frozen organism policy"
    if library_layout != "PAIRED":
        return "BLOCKED_METADATA", f"REPORTED_LAYOUT_{library_layout or 'MISSING'}", "reported LibraryLayout is not PAIRED"
    if spots < 0 or spots_with_mates < 0:
        return "BLOCKED_METADATA", "RUNINFO_SPOT_COUNTS_MISSING", "spots or spots_with_mates is missing"

    if submission == "SRA1825115" and spots_with_mates == 0:
        if sra_stat_structure == "SINGLE_BIOLOGICAL_READ_PER_SPOT":
            return (
                "BLOCKED_READ_STRUCTURE",
                "SINGLE_BIOLOGICAL_READ_PER_SPOT",
                "sra-stat spot_count_mates=0; one biological read per spot; reported LibraryLayout PAIRED is not paired biological reads",
            )
        return (
            "BLOCKED_READ_STRUCTURE",
            "RUNINFO_SPOTS_WITH_MATES_ZERO",
            "this run's runinfo spots_with_mates=0; paired biological reads are not supported; sra-stat structure was not measured and is not inferred from the two audited SRA1825115 runs",
        )

    if (
        submission == "SRA1824580"
        and spots > 0
        and spots_with_mates == spots
        and library_layout == "PAIRED"
    ):
        if sra_stat_structure == "TRUE_PAIRED_READS_AVAILABLE":
            structure = "TRUE_PAIRED_READS_AVAILABLE"
            reason = "."
        elif sra_stat_structure in (None, ""):
            structure = "RUNINFO_SPOTS_WITH_MATES_EQUAL_SPOTS"
            reason = "."
        else:
            return (
                "BLOCKED_READ_STRUCTURE",
                sra_stat_structure,
                "sra-stat structure contradicts runinfo spots_with_mates equal to spots",
            )
        return "ELIGIBLE_PAIRED", structure, reason

    return (
        "BLOCKED_METADATA",
        "RUNINFO_MATE_COUNT_NOT_PAIRED",
        f"submission {submission} spots={spots} spots_with_mates={spots_with_mates} does not meet the paired-run rule",
    )


def truth_by_sample(coordinates: list[dict[str, str]], truth_rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    by_name: dict[str, dict[str, object]] = {}
    for row in truth_rows:
        name = row["sample_name"]
        slot = by_name.setdefault(name, {"none": False, "rows": 0})
        slot["rows"] = int(slot["rows"]) + 1
        if row["reported_mutation"] == "NONE":
            slot["none"] = True
    for row in coordinates:
        name = row["sample_name"]
        slot = by_name.setdefault(name, {"none": False, "rows": 0})
        alleles = slot.setdefault("eligible", set())
        unresolved = slot.setdefault("unresolved", [])
        assert isinstance(alleles, set)
        assert isinstance(unresolved, list)
        if row["tier_a_eligibility"] == "ELIGIBLE_EXACT_RECOVERY" and row["denominator"] == "include":
            alleles.add((row["reported_position"], row["ref_allele"], row["alt_allele"]))
        elif row["denominator"] == "exclude" or row["coordinate_status"] == "COORDINATE_UNRESOLVED":
            unresolved.append((row["reported_position"], row["ref_allele"], row["alt_allele"], row["coordinate_status"]))
    summary: dict[str, dict[str, object]] = {}
    for name, slot in by_name.items():
        eligible = slot.get("eligible", set())
        unresolved = slot.get("unresolved", [])
        assert isinstance(eligible, set)
        assert isinstance(unresolved, list)
        if unresolved and not eligible:
            truth_class = "COORDINATE_UNRESOLVED"
        elif eligible:
            truth_class = "COORDINATE_VERIFIED"
        elif slot.get("none"):
            truth_class = "PUBLICATION_NONE"
        else:
            truth_class = "TRUTH_MISSING"
        summary[name] = {
            "truth_class": truth_class,
            "eligible_publication_event_count": len(eligible),
        }
    return summary


def build_rows() -> list[dict[str, str]]:
    relationships = read_tsv(RELATIONSHIPS)
    if len(relationships) != 25:
        raise ReferencePolicyError(f"expected 25 relationship rows, found {len(relationships)}")
    runinfo = {row["Run"]: row for row in csv.DictReader(RUNINFO.open(encoding="utf-8", newline=""))}
    coordinates = read_tsv(COORDINATES)
    truth_rows = read_tsv(TRUTH)
    truth = truth_by_sample(coordinates, truth_rows)
    built: list[dict[str, str]] = []
    for rel in relationships:
        sample = rel["sample_name"]
        if sample not in truth:
            raise ReferencePolicyError(f"{sample} has no publication truth row")
        publication, execution = expected_accessions(rel["sra_organism"])
        if rel["publication_reference"] != publication or rel["execution_reference"] != execution:
            raise ReferencePolicyError(f"{sample} relationship reference drifted")
        info = truth[sample]
        pairs = (
            ("run_long", "long_experiment", "long_submission"),
            ("run_short", "short_experiment", "short_submission"),
        )
        for run_key, exp_key, sub_key in pairs:
            run = rel[run_key]
            meta = runinfo.get(run)
            identity_ok = meta is not None and meta["SampleName"] == sample and meta["BioSample"] == rel["biosample"] and meta["Experiment"] == rel[exp_key] and meta["Submission"] == rel[sub_key] and meta["BioProject"] == "PRJNA1088182"
            reference_ok = False
            if identity_ok and meta is not None:
                try:
                    require_execution_reference(run, execution, relationships)
                    reference_ok = meta["ScientificName"] == rel["sra_organism"]
                except ReferencePolicyError:
                    reference_ok = False
            spots = int(meta["spots"]) if meta is not None else -1
            mates = int(meta["spots_with_mates"]) if meta is not None else -1
            layout = meta["LibraryLayout"] if meta is not None else ""
            submission = meta["Submission"] if meta is not None else rel[sub_key]
            eligibility, structure, reason = classify_run(
                submission=submission,
                spots=spots,
                spots_with_mates=mates,
                library_layout=layout,
                identity_ok=identity_ok,
                reference_ok=reference_ok,
                sra_stat_structure=SRA_STAT_VERIFIED.get(run),
            )
            built.append(
                {
                    "sample_name": sample,
                    "biosample": rel["biosample"],
                    "experiment": rel[exp_key],
                    "submission": rel[sub_key],
                    "run": run,
                    "strain": rel["sra_organism"],
                    "library_layout_reported": layout or ".",
                    "observed_read_structure": structure,
                    "execution_reference": execution,
                    "reference_sequence_sha256": SHA256[execution],
                    "truth_class": str(info["truth_class"]),
                    "eligible_publication_event_count": str(info["eligible_publication_event_count"]),
                    "execution_eligibility": eligibility,
                    "exclusion_reason": reason,
                }
            )
    return built


def new_execution_runs(rows: list[dict[str, str]] | None = None) -> list[str]:
    """Eligible paired runs that are not the two verified P1b smoke executions."""
    rows = rows if rows is not None else build_rows()
    reused = set(REUSED_VERIFIED_P1B)
    return [
        row["run"]
        for row in rows
        if row["execution_eligibility"] == "ELIGIBLE_PAIRED" and row["run"] not in reused
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description="Freeze or list the M6-P2a cohort.")
    parser.add_argument("--list-new-runs", action="store_true", help="print eligible runs that still need evidence")
    args = parser.parse_args()
    rows = build_rows()
    if args.list_new_runs:
        for run in new_execution_runs(rows):
            print(run)
        return 0
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    eligible = sum(row["execution_eligibility"] == "ELIGIBLE_PAIRED" for row in rows)
    print(f"wrote {OUT} rows={len(rows)} eligible_paired={eligible}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
