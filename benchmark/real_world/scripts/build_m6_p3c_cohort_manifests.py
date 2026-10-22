#!/usr/bin/env python3
"""Build M6-P3c Benchmark Quintet Cohort Manifests.

Aggregates execution telemetry, 4D structural reconciliation, decoupled MOB concordance,
and bidirectional JC event accounting across all 5 benchmark samples in PRJNA1088182:
  1. SRR28349173 (Cas9_110_111, MG1655 / NC_000913.3)
  2. SRR28349176 (Cas9_104_105, MG1655 / NC_000913.3)
  3. SRR28349244 (Cas9_1,       MG1655 / NC_000913.3) [Pilot]
  4. SRR28349231 (I-SceI_3..6,  BW25113 / CP009273.1)
  5. SRR28349226 (control_1..4, BW25113 / CP009273.1)
"""

import csv
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
RESULTS_DIR = REPO / "benchmark/real_world/results/prjna1088182/m6_p3c"
MANIFESTS_DIR = REPO / "benchmark/real_world/manifests"

RUN_METADATA = [
    {
        "run": "SRR28349173",
        "sample_name": "Cas9_110_111",
        "slurm_job_id": "2999176",
        "node": "bnode1",
        "reference": "NC_000913.3",
    },
    {
        "run": "SRR28349176",
        "sample_name": "Cas9_104_105",
        "slurm_job_id": "2999177",
        "node": "bnode2",
        "reference": "NC_000913.3",
    },
    {
        "run": "SRR28349244",
        "sample_name": "Cas9_1",
        "slurm_job_id": "2997962",
        "node": "bnode2",
        "reference": "NC_000913.3",
    },
    {
        "run": "SRR28349231",
        "sample_name": "I-SceI_3_4_5_6",
        "slurm_job_id": "3000148",
        "node": "bnode1",
        "reference": "CP009273.1",
    },
    {
        "run": "SRR28349226",
        "sample_name": "control_1_2_3_4",
        "slurm_job_id": "2999175",
        "node": "bnode2",
        "reference": "CP009273.1",
    },
]


def parse_time_v(telemetry_path: Path) -> dict:
    text = telemetry_path.read_text(encoding="utf-8")
    data = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if "User time (seconds):" in line:
            data["user_cpu_seconds"] = line.split(":", 1)[1].strip()
        elif "System time (seconds):" in line:
            data["system_cpu_seconds"] = line.split(":", 1)[1].strip()
        elif "Percent of CPU this job got:" in line:
            data["cpu_percent"] = line.split(":", 1)[1].strip().rstrip("%")
        elif "Elapsed (wall clock) time" in line:
            # line is e.g. "Elapsed (wall clock) time (h:mm:ss or m:ss): 29:32.48"
            v = line.split("):", 1)[1].strip()
            parts = v.split(":")
            if len(parts) == 2:
                mins, secs = parts
                total_s = float(mins) * 60 + float(secs)
            elif len(parts) == 3:
                hrs, mins, secs = parts
                total_s = float(hrs) * 3600 + float(mins) * 60 + float(secs)
            else:
                total_s = float(v)
            data["wall_clock_seconds"] = f"{total_s:.2f}"
        elif "Maximum resident set size (kbytes):" in line:
            data["max_rss_kb"] = line.split(":", 1)[1].strip()
        elif "Exit status:" in line:
            data["exit_code"] = line.split(":", 1)[1].strip()
    return data


def parse_muts(path: Path) -> dict[str, set[int]]:
    muts = {"SNP": set(), "INS": set(), "DEL": set(), "MOB": set(), "SUB": set()}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        toks = line.split("\t")
        t = toks[0]
        if t in muts:
            pos = int(toks[4])
            muts[t].add(pos)
    return muts


def build_cohort_manifests():
    MANIFESTS_DIR.mkdir(parents=True, exist_ok=True)

    telemetry_rows = []
    reconciliation_rows = []
    mob_rows = []
    matrix_rows = []

    for meta in RUN_METADATA:
        run = meta["run"]
        sample_name = meta["sample_name"]
        run_dir = RESULTS_DIR / run

        # 1. Telemetry
        time_data = parse_time_v(run_dir / "breseq/telemetry.txt")
        summary_data = json.loads((run_dir / "reconciliation_summary.json").read_text(encoding="utf-8"))

        pk_gd = REPO / f"benchmark/real_world/results/prjna1088182/m6_p2a/{run}/output.gd"
        if not pk_gd.is_file():
            pk_gd = REPO / f"benchmark/real_world/results/prjna1088182/m6_p1b/{run}/output.gd"
        breseq_gd = run_dir / "breseq/output/output.gd"

        pk_muts = parse_muts(pk_gd)
        b_muts = parse_muts(breseq_gd)

        snps_conc = f"{len(pk_muts['SNP'] & b_muts['SNP'])}/{len(pk_muts['SNP'])}"
        ins_conc = f"{len(pk_muts['INS'] & b_muts['INS'])}/{len(pk_muts['INS'])}"
        del_conc = f"{len(pk_muts['DEL'] & b_muts['DEL'])}/{len(pk_muts['DEL'])}"
        mob_conc = f"{len(pk_muts['MOB'] & b_muts['MOB'])}/{len(pk_muts['MOB'])}"

        # Counts from summary and files
        pk_total_jcs = summary_data["total_records"]
        comp_total_jcs = summary_data.get("matching_comparison", {}).get("comparator_total")
        if comp_total_jcs is None or comp_total_jcs == 0:
            comp_total_jcs = sum(1 for line in breseq_gd.read_text(encoding="utf-8").splitlines() if line.startswith("JC\t"))
        exact_matches = summary_data["dim3_corroboration"].get("CORROBORATED", 0)
        family_matches = summary_data["dim3_corroboration"].get("CORROBORATED_AT_FAMILY_LEVEL_UNRESOLVED_AT_LOCUS", 0)
        uncalled_jcs = summary_data["dim3_corroboration"].get("NOT_CORROBORATED", 0)

        t_row = {
            "run": run,
            "sample_name": sample_name,
            "slurm_job_id": meta["slurm_job_id"],
            "node": meta["node"],
            "partition": "qcpu_18i",
            "exit_code": time_data.get("exit_code", "0"),
            "wall_clock_seconds": time_data.get("wall_clock_seconds", "."),
            "user_cpu_seconds": time_data.get("user_cpu_seconds", "."),
            "system_cpu_seconds": time_data.get("system_cpu_seconds", "."),
            "cpu_percent": time_data.get("cpu_percent", "."),
            "max_rss_kb": time_data.get("max_rss_kb", "."),
            "prokadiff_jcs": pk_total_jcs,
            "breseq_jcs": comp_total_jcs,
            "exact_matches": exact_matches,
            "family_matches": family_matches,
            "uncalled_jcs": uncalled_jcs,
            "snps_concordant": snps_conc,
            "ins_concordant": ins_conc,
            "del_concordant": del_conc,
            "mob_concordant": mob_conc,
            "dim1_unique_anchor": summary_data["dim1_placement"].get("UNIQUE_ANCHOR", 0),
            "dim1_repeat_anchored": summary_data["dim1_placement"].get("REPEAT_ANCHORED", 0),
            "dim1_repeat_ambiguous": summary_data["dim1_placement"].get("REPEAT_AMBIGUOUS", 0),
            "dim2_strong": summary_data["dim2_read_evidence"].get("STRONG", 0),
            "dim2_limited": summary_data["dim2_read_evidence"].get("LIMITED", 0),
            "dim3_corroborated": exact_matches,
            "dim3_family_level": family_matches,
            "dim3_not_corroborated": uncalled_jcs,
            "dim4_biological_validation": summary_data["dim4_biological_validation"],
        }
        telemetry_rows.append(t_row)

        # 2. Structural Reconciliation rows
        reconciled_p = run_dir / "reconciled_structural_events.tsv"
        with reconciled_p.open(encoding="utf-8") as f:
            r_reader = csv.DictReader(f, delimiter="\t")
            for row in r_reader:
                row_copy = {"run": run, "sample_name": sample_name}
                row_copy.update(row)
                reconciliation_rows.append(row_copy)

        # 3. MOB concordance rows
        mob_p = run_dir / "mob_concordance.tsv"
        with mob_p.open(encoding="utf-8") as f:
            m_reader = csv.DictReader(f, delimiter="\t")
            for row in m_reader:
                row_copy = {"run": run, "sample_name": sample_name}
                row_copy.update(row)
                mob_rows.append(row_copy)

        # 4. Bidirectional JC matrix rows
        matrix_p = run_dir / "bidirectional_jc_matrix.tsv"
        with matrix_p.open(encoding="utf-8") as f:
            mat_reader = csv.DictReader(f, delimiter="\t")
            for row in mat_reader:
                row_copy = {"run": run, "sample_name": sample_name}
                row_copy.update(row)
                matrix_rows.append(row_copy)

    # Write cohort telemetry
    telemetry_out = MANIFESTS_DIR / "prjna1088182_m6_p3c_cohort_telemetry.tsv"
    with telemetry_out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(telemetry_rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(telemetry_rows)
    print(f"Wrote {len(telemetry_rows)} cohort telemetry records to {telemetry_out}")

    # Write cohort reconciliation
    reconciliation_out = MANIFESTS_DIR / "prjna1088182_m6_p3c_cohort_reconciliation.tsv"
    with reconciliation_out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(reconciliation_rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(reconciliation_rows)
    print(f"Wrote {len(reconciliation_rows)} cohort reconciliation records to {reconciliation_out}")

    # Write cohort MOB concordance
    mob_out = MANIFESTS_DIR / "prjna1088182_m6_p3c_cohort_mob_concordance.tsv"
    mob_fields = ["run", "sample_name", "prokadiff_id", "comparator_id", "seq_id", "prokadiff_pos", "comparator_pos",
                  "site_concordance", "prokadiff_family", "comparator_family", "family_concordance",
                  "prokadiff_strand", "comparator_strand", "strand_concordance", "prokadiff_tsd", "comparator_tsd",
                  "tsd_concordance", "composite_classification", "notes"]
    with mob_out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=mob_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(mob_rows)
    print(f"Wrote {len(mob_rows)} cohort MOB concordance records to {mob_out}")

    # Write cohort bidirectional matrix
    matrix_out = MANIFESTS_DIR / "prjna1088182_m6_p3c_cohort_bidirectional_matrix.tsv"
    mat_fields = ["run", "sample_name", "category", "prokadiff_id", "breseq_id", "prokadiff_coords", "breseq_coords",
                  "prokadiff_support", "breseq_support", "delta_breakend1", "delta_breakend2", "total_distance",
                  "reconciliation_status", "notes"]
    with matrix_out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=mat_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(matrix_rows)
    print(f"Wrote {len(matrix_rows)} cohort bidirectional matrix records to {matrix_out}")


if __name__ == "__main__":
    build_cohort_manifests()
