#!/usr/bin/env python3
"""Fail closed when an M6 run is paired with the wrong execution reference.

Authoritative selector: prjna1088182_relationships.tsv column execution_reference.
The cohort manifest must carry the same accession. It is not a second policy.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

MG1655_EXECUTION = "NC_000913.3"
MG1655_PUBLICATION = "U00096.3"
BW25113_EXECUTION = "CP009273.1"
BW25113_PUBLICATION = "CP009273.1"

REPO_ROOT = Path(__file__).resolve().parents[3]
RELATIONSHIPS = REPO_ROOT / "benchmark/real_world/manifests/prjna1088182_relationships.tsv"
COHORT = REPO_ROOT / "benchmark/real_world/manifests/prjna1088182_cohort.tsv"
COORDINATES = REPO_ROOT / "benchmark/real_world/manifests/prjna1088182_truth_coordinate_audit.tsv"

SMOKE_RUNS = {
    "SRR28349244": "Cas9_1",
    "SRR28357593": "Cas9_1",
    "SRR28349243": "Cas9_2",
    "SRR28357592": "Cas9_2",
}


class ReferencePolicyError(Exception):
    pass


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_mg1655(organism: str) -> bool:
    return "MG1655" in organism


def is_bw25113(organism: str) -> bool:
    return "BW25113" in organism


def expected_accessions(organism: str) -> tuple[str, str]:
    if is_mg1655(organism) and not is_bw25113(organism):
        return MG1655_PUBLICATION, MG1655_EXECUTION
    if is_bw25113(organism):
        return BW25113_PUBLICATION, BW25113_EXECUTION
    raise ReferencePolicyError(f"unclassified SRA organism: {organism}")


def index_runs(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        for key in ("cohort_run", "run_short", "run_long"):
            run = row[key]
            if run in indexed and indexed[run]["sample_name"] != row["sample_name"]:
                raise ReferencePolicyError(f"run {run} maps to more than one sample")
            indexed[run] = row
    return indexed


def require_execution_reference(run_id: str, requested: str, rows: list[dict[str, str]]) -> str:
    """Return the frozen execution accession or raise if the request contradicts it."""
    row = index_runs(rows).get(run_id)
    if row is None:
        raise ReferencePolicyError(f"run {run_id} is not in the M6 relationship table")
    publication, execution = expected_accessions(row["sra_organism"])
    if row["publication_reference"] != publication or row["execution_reference"] != execution:
        raise ReferencePolicyError(
            f"{row['sample_name']} {run_id} relationship reference "
            f"{row['publication_reference']}/{row['execution_reference']} "
            f"disagrees with organism policy {publication}/{execution}"
        )
    if requested != execution:
        raise ReferencePolicyError(
            f"{row['sample_name']} {run_id} requested {requested}; execution reference is {execution}"
        )
    if is_bw25113(row["sra_organism"]) and requested == MG1655_EXECUTION:
        raise ReferencePolicyError("BW25113 must not execute against NC_000913.3")
    if is_mg1655(row["sra_organism"]) and requested == BW25113_EXECUTION:
        raise ReferencePolicyError("MG1655 must not execute against CP009273.1")
    return execution


def check_manifests(
    relationships: list[dict[str, str]] | None = None,
    cohort: list[dict[str, str]] | None = None,
    coordinates: list[dict[str, str]] | None = None,
) -> None:
    relationships = relationships if relationships is not None else read_tsv(RELATIONSHIPS)
    cohort = cohort if cohort is not None else read_tsv(COHORT)
    coordinates = coordinates if coordinates is not None else read_tsv(COORDINATES)
    if len(relationships) != 25:
        raise ReferencePolicyError(f"expected 25 relationship rows, found {len(relationships)}")
    by_name = {row["sample_name"]: row for row in relationships}
    if len(by_name) != 25:
        raise ReferencePolicyError("relationship sample names are not unique")
    for row in relationships:
        publication, execution = expected_accessions(row["sra_organism"])
        if row["publication_reference"] != publication or row["execution_reference"] != execution:
            raise ReferencePolicyError(f"{row['sample_name']} failed organism/reference policy")
        if row["execution_unit"] != "SRA_Run":
            raise ReferencePolicyError(f"{row['sample_name']} execution_unit is not SRA_Run")
        if row["publication_truth_unit"] != "SampleName/BioSample":
            raise ReferencePolicyError(f"{row['sample_name']} publication truth unit changed")
        if row["short_class"] != "INDEPENDENT_WGS_RUN" or row["long_class"] != "INDEPENDENT_WGS_RUN":
            raise ReferencePolicyError(f"{row['sample_name']} run class is not an independent WGS run")
        if row["parent_class"] != "NO_PUBLIC_RUN_FOUND":
            raise ReferencePolicyError(f"{row['sample_name']} parent class is not NO_PUBLIC_RUN_FOUND")
        if row["starter_role"] != "computational_comparator_not_verified_parent":
            raise ReferencePolicyError(f"{row['sample_name']} starter_id was promoted to a parent")
        if row["cohort_role"] != "DEVELOPMENT_SUBSET":
            raise ReferencePolicyError(f"{row['sample_name']} is no longer marked DEVELOPMENT_SUBSET")
    if len(cohort) != 25:
        raise ReferencePolicyError(f"expected 25 cohort rows, found {len(cohort)}")
    for row in cohort:
        rel = by_name.get(row["sample_name"])
        if rel is None:
            raise ReferencePolicyError(f"cohort sample {row['sample_name']} has no relationship row")
        if row.get("execution_reference") != rel["execution_reference"]:
            raise ReferencePolicyError(
                f"{row['sample_name']} cohort execution_reference disagrees with relationships"
            )
        if row.get("publication_reference") != rel["publication_reference"]:
            raise ReferencePolicyError(
                f"{row['sample_name']} cohort publication_reference disagrees with relationships"
            )
        if row.get("reference") != rel["execution_reference"]:
            raise ReferencePolicyError(
                f"{row['sample_name']} cohort reference is not the execution accession"
            )
        if "MG1655.gbk" in row.get("reference", "") or "MG1655.gbk" in row.get("execution_reference", ""):
            raise ReferencePolicyError(f"{row['sample_name']} still uses a legacy MG1655 path as the reference selector")
    runs = index_runs(relationships)
    for run_id in SMOKE_RUNS:
        if require_execution_reference(run_id, MG1655_EXECUTION, relationships) != MG1655_EXECUTION:
            raise ReferencePolicyError(f"smoke run {run_id} is not NC_000913.3")
        if runs[run_id]["sample_name"] != SMOKE_RUNS[run_id]:
            raise ReferencePolicyError(f"smoke run {run_id} changed sample identity")
    eligible = [row for row in coordinates if row["tier_a_eligibility"] == "ELIGIBLE_EXACT_RECOVERY"]
    excluded = [row for row in coordinates if row["denominator"] == "exclude"]
    if len(eligible) != 13:
        raise ReferencePolicyError(f"expected 13 eligible truth rows, found {len(eligible)}")
    alleles = {
        (row["sample_id"], row["reported_position"], row["ref_allele"], row["alt_allele"])
        for row in eligible
    }
    if len(alleles) != 12:
        raise ReferencePolicyError(f"expected 12 independent eligible alleles, found {len(alleles)}")
    if len(excluded) != 1 or excluded[0]["sample_name"] != "lambdaRed_3_4":
        raise ReferencePolicyError("the unresolved +AAT row is not the sole excluded truth row")
    if excluded[0]["reported_position"] != "3112908" or excluded[0]["alt_allele"] != "+AAT":
        raise ReferencePolicyError("excluded truth row is not 3112908 +AAT")
    if excluded[0]["coordinate_status"] != "COORDINATE_UNRESOLVED":
        raise ReferencePolicyError("+AAT was promoted out of COORDINATE_UNRESOLVED")
    for row in eligible:
        if row["coordinate_status"] != "COORDINATE_VERIFIED" or row["denominator"] != "include":
            raise ReferencePolicyError("an eligible row is not a verified included coordinate")


def main() -> int:
    parser = argparse.ArgumentParser(description="Check the frozen M6 Tier-A reference policy.")
    parser.add_argument("--check", action="store_true", help="validate the versioned manifests")
    args = parser.parse_args()
    if not args.check:
        parser.error("--check is required")
    try:
        check_manifests()
    except (ReferencePolicyError, OSError, KeyError) as error:
        print(f"M6 reference policy failed: {error}", file=sys.stderr)
        return 1
    print("M6 Tier-A reference policy passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
