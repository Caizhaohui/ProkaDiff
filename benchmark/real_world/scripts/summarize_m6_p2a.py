#!/usr/bin/env python3
"""Collect M6-P2a staging and evidence records into small TSV summaries.

Missing outputs stay in the tables. They are not dropped from the denominator.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_m6_p2a_cohort import REUSED_VERIFIED_P1B  # noqa: E402
from m6_reference_policy import COORDINATES, read_tsv  # noqa: E402
from reconcile_m6_p2a import (  # noqa: E402
    MalformedGdError,
    deduplicated_evaluations,
    inventory,
    parse_gd,
    repeat_diagnostic_count,
)

REPO = Path(__file__).resolve().parents[3]
COHORT = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p2a_execution_cohort.tsv"
P1B = REPO / "benchmark/real_world/results/prjna1088182/m6_p1b"
P2A = REPO / "benchmark/real_world/results/prjna1088182/m6_p2a"
RUN_OUT = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p2a_run_results.tsv"
EVENT_OUT = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p2a_event_recovery.tsv"

RUN_COLUMNS = [
    "sample",
    "run",
    "strain",
    "reference",
    "input_status",
    "execution_status",
    "job_id",
    "exit_code",
    "wall_time",
    "max_rss",
    "total_mutation_calls",
    "MC_count",
    "repeat_diagnostic_count",
    "publication_eligible_event_count",
    "recovered_exact",
    "recovered_equivalent",
    "representation_difference",
    "not_recovered",
    "not_evaluable",
    "snp_calls",
    "sub_calls",
    "ins_calls",
    "del_calls",
    "jc_calls",
    "mob_calls",
    "evidence_gd",
    "binary_sha256",
    "prokadiff_commit",
]

EVENT_COLUMNS = [
    "sample",
    "run",
    "publication_event",
    "truth_source",
    "truth_coordinate_status",
    "reference",
    "eligible_denominator_flag",
    "prokadiff_observed_event",
    "recovery_classification",
    "diagnostic_reason",
    "evidence_location",
    "gd_record_id",
    "dv1_event_id",
    "review_stop",
]


def parse_telemetry(text: str) -> dict[str, str]:
    wall, rss, exit_code = ".", ".", "."
    for line in text.splitlines():
        if "Elapsed (wall clock) time" in line:
            match = re.search(r"(\d+:\d{2}:\d{2}|\d+:\d{2}(?:\.\d+)?)\s*$", line)
            if match:
                wall = match.group(1)
        elif "Maximum resident set size" in line and ":" in line:
            rss = line.split(":", 1)[1].strip()
            if rss.isdigit():
                rss = f"{rss}K"
        elif line.strip().startswith("Exit status:"):
            exit_code = line.split(":", 1)[1].strip()
    return {"wall_time": wall, "max_rss": rss, "exit_code": exit_code}


def _blank_counts() -> dict[str, str]:
    return {
        "total_mutation_calls": ".",
        "MC_count": ".",
        "repeat_diagnostic_count": ".",
        "recovered_exact": ".",
        "recovered_equivalent": ".",
        "representation_difference": ".",
        "not_recovered": ".",
        "not_evaluable": ".",
        "snp_calls": ".",
        "sub_calls": ".",
        "ins_calls": ".",
        "del_calls": ".",
        "jc_calls": ".",
        "mob_calls": ".",
    }


def reference_bases(accession: str) -> dict[int, str]:
    from m6_p2a_stage_gate import REFERENCE_FILES, chromosome_bases

    bases = chromosome_bases(REFERENCE_FILES[accession])
    return {index + 1: base for index, base in enumerate(bases)}


def load_audit_by_sample() -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in read_tsv(COORDINATES):
        grouped.setdefault(row["sample_name"], []).append(row)
    return grouped


def command_is_p1b_compatible(telemetry_text: str, run: str) -> bool:
    command = ""
    for line in telemetry_text.splitlines():
        if line.strip().startswith("Command being timed:"):
            command = line
            break
    if " evidence " not in command and " evidence\n" not in command and "prokadiff evidence" not in command:
        return False
    if "--starter" in command:
        return False
    if "MG1655.gbk" not in command:
        return False
    if f"{run}_1.fastq.gz" not in command or f"{run}_2.fastq.gz" not in command:
        return False
    return True


def summarize(cohort_rows: list[dict[str, str]] | None = None) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    cohort_rows = cohort_rows if cohort_rows is not None else list(csv.DictReader(COHORT.open(encoding="utf-8"), delimiter="\t"))
    audit = load_audit_by_sample()
    base_cache: dict[str, dict[int, str]] = {}
    run_rows: list[dict[str, str]] = []
    event_rows: list[dict[str, str]] = []
    for row in cohort_rows:
        sample = row["sample_name"]
        run = row["run"]
        accession = row["execution_reference"]
        result = {
            "sample": sample,
            "run": run,
            "strain": row["strain"],
            "reference": accession,
            "input_status": ".",
            "execution_status": row["execution_eligibility"],
            "job_id": ".",
            "exit_code": ".",
            "wall_time": ".",
            "max_rss": ".",
            "publication_eligible_event_count": row["eligible_publication_event_count"],
            "evidence_gd": ".",
            "binary_sha256": ".",
            "prokadiff_commit": ".",
            **_blank_counts(),
        }
        if row["execution_eligibility"] != "ELIGIBLE_PAIRED":
            result["input_status"] = "NOT_STAGED_FOR_P2A"
            run_rows.append(result)
            continue
        evidence_dir = P1B / run if run in REUSED_VERIFIED_P1B else P2A / run
        gd_path = evidence_dir / "output.gd"
        telemetry_path = evidence_dir / "telemetry.txt"
        if run in REUSED_VERIFIED_P1B:
            result["execution_status"] = "REUSED_VERIFIED_P1B"
            result["job_id"] = "2897505"
            result["prokadiff_commit"] = "702505105196ade5eed449d553fb6656d31dead8"
            result["binary_sha256"] = "02077c761b1c9bd83cf1b626a3e2e955d8d2a1c63abb3e35852421e0dbc20bfd"
        else:
            result["execution_status"] = "NOT_RUN"
        prov = REPO / "data/prjna1088182/fastq" / run / "run_provenance.tsv"
        if prov.is_file():
            staged = list(csv.DictReader(prov.open(encoding="utf-8"), delimiter="\t"))
            if len(staged) == 1 and staged[0].get("staging_status") == "READY_FOR_EXECUTION":
                result["input_status"] = "READY_FOR_EXECUTION"
            else:
                result["input_status"] = staged[0].get("staging_status", "STAGING_NOT_READY") if staged else "STAGING_NOT_READY"
        else:
            result["input_status"] = "STAGING_PROVENANCE_MISSING"
        if not gd_path.is_file() or not telemetry_path.is_file():
            run_rows.append(result)
            continue
        telemetry = telemetry_path.read_text(encoding="utf-8", errors="replace")
        if run in REUSED_VERIFIED_P1B and not command_is_p1b_compatible(telemetry, run):
            result["execution_status"] = "P1B_NOT_COMPATIBLE"
            run_rows.append(result)
            continue
        parsed = parse_telemetry(telemetry)
        result.update(parsed)
        if run not in REUSED_VERIFIED_P1B:
            result["execution_status"] = "COMPLETED" if result["exit_code"] == "0" else "FAILED"
            command = (evidence_dir / "command.txt").read_text(encoding="utf-8", errors="replace") if (evidence_dir / "command.txt").is_file() else ""
            for line in command.splitlines():
                if line.startswith("slurm_job_id "):
                    result["job_id"] = line.split(" ", 1)[1]
                elif line.startswith("commit "):
                    result["prokadiff_commit"] = line.split(" ", 1)[1]
                elif line.startswith("binary_sha256 "):
                    result["binary_sha256"] = line.split(" ", 1)[1]
        try:
            records = parse_gd(gd_path.read_text(encoding="utf-8", errors="replace"))
        except MalformedGdError as error:
            result["execution_status"] = f"MALFORMED_GD:{error}"
            run_rows.append(result)
            continue
        counts = inventory(records)
        repeat_path = evidence_dir / "repeat_ambiguous_junctions.tsv"
        repeat_count = repeat_diagnostic_count(repeat_path.read_text(encoding="utf-8", errors="replace")) if repeat_path.is_file() else 0
        result.update(
            {
                "total_mutation_calls": str(counts["total_mutation_calls"]),
                "MC_count": str(counts["MC_count"]),
                "repeat_diagnostic_count": str(repeat_count),
                "snp_calls": str(counts["SNP"]),
                "sub_calls": str(counts["SUB"]),
                "ins_calls": str(counts["INS"]),
                "del_calls": str(counts["DEL"]),
                "jc_calls": str(counts["JC"]),
                "mob_calls": str(counts["MOB"]),
                "evidence_gd": str(gd_path.relative_to(REPO)),
            }
        )
        if accession not in base_cache:
            base_cache[accession] = reference_bases(accession)
        sample_audit = audit.get(sample, [])
        evaluations = deduplicated_evaluations(sample_audit, records, accession, base_cache[accession])
        tallies = {"RECOVERED_EXACT": 0, "RECOVERED_EQUIVALENT": 0, "REPRESENTATION_DIFFERENCE": 0, "NOT_RECOVERED": 0, "NOT_EVALUABLE": 0}
        if not sample_audit and row["truth_class"] == "PUBLICATION_NONE":
            event_rows.append(
                {
                    "sample": sample,
                    "run": run,
                    "publication_event": "NONE",
                    "truth_source": "PUBLICATION_REPORTED",
                    "truth_coordinate_status": "PUBLICATION_NONE",
                    "reference": accession,
                    "eligible_denominator_flag": "exclude",
                    "prokadiff_observed_event": "NA",
                    "recovery_classification": "NOT_EVALUABLE",
                    "diagnostic_reason": ".",
                    "evidence_location": str(gd_path.relative_to(REPO)),
                    "gd_record_id": "NA",
                    "dv1_event_id": "NOT_EVALUABLE",
                    "review_stop": ".",
                }
            )
        for item in evaluations:
            tallies[item["recovery"]] += 1
            event_rows.append(
                {
                    "sample": sample,
                    "run": run,
                    "publication_event": item["publication_event"],
                    "truth_source": "PUBLICATION_REPORTED",
                    "truth_coordinate_status": next(
                        audit_row["coordinate_status"]
                        for audit_row in sample_audit
                        if f"{audit_row['reported_position']} {audit_row['ref_allele']}>{audit_row['alt_allele']}" == item["publication_event"]
                    ),
                    "reference": accession,
                    "eligible_denominator_flag": item["eligible_denominator"],
                    "prokadiff_observed_event": item["observed_event"],
                    "recovery_classification": item["recovery"],
                    "diagnostic_reason": item["diagnostic_reason"],
                    "evidence_location": str(gd_path.relative_to(REPO)),
                    "gd_record_id": item["gd_record_id"],
                    "dv1_event_id": item["dv1_event_id"],
                    "review_stop": item["review_stop"],
                }
            )
        result["recovered_exact"] = str(tallies["RECOVERED_EXACT"])
        result["recovered_equivalent"] = str(tallies["RECOVERED_EQUIVALENT"])
        result["representation_difference"] = str(tallies["REPRESENTATION_DIFFERENCE"])
        result["not_recovered"] = str(tallies["NOT_RECOVERED"])
        result["not_evaluable"] = str(tallies["NOT_EVALUABLE"])
        run_rows.append(result)
    return run_rows, event_rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Write M6-P2a run and event summaries.")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    runs, events = summarize()
    if not args.write:
        print(f"runs={len(runs)} events={len(events)}")
        return 0
    with RUN_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RUN_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(runs)
    with EVENT_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EVENT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(events)
    print(f"wrote {RUN_OUT}")
    print(f"wrote {EVENT_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
