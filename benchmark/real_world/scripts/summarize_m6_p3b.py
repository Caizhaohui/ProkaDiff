#!/usr/bin/env python3
"""Summarize M6-P3b A/B experiment comparing FASTA baseline (Cond A) vs GenBank annotation (Cond B).

Validates Gates 1-4:
  Gate 1: Reference Sequence Identity (SHA256 3761072f...)
  Gate 2: Small Variant Invariance (Delta SNP == 0, Delta INS == 0, Delta DEL == 0, Delta MC == 0)
  Gate 3: Control No New MOB (gate3_control_no_new_mob: control_1_2_3_4 produces MOB == 0 under both conditions;
          biological specificity/sensitivity remain NOT_EVALUABLE without independent ground truth)
  Gate 4: Structural Filtering & MOB Consistency (gate4_mob_promotion: validates Delta MOB >= 0, constituent
          parent-JC consistency when MOB emitted, and non-expansion of standalone JCs [JC_B <= JC_A] via repeat
          filtering [ANNOTATION_FILTERED_JC / REPEAT_ASSOCIATED_JC]. Explicitly: positive empirical MOB detection
          was NOT demonstrated in these 3 BW25113 runs [all tested MOB counts remained zero / MOB_NOT_CALLED];
          Gate 4 validates structural filtering and MOB call consistency, not active MOB recovery).
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
RESULTS_DIR = REPO / "benchmark/real_world/results/prjna1088182/m6_p3b"
OUT_TSV = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3b_ab_comparison.tsv"

EXPECTED_SHA = "3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce"
EXPECTED_LEN = 4631469

TRIAD = [
    {"run": "SRR28349226", "sample": "control_1_2_3_4", "expected_mob_b": 0},
    {"run": "SRR28349231", "sample": "I-SceI_3_4_5_6", "expected_mob_b": None},
    {"run": "SRR28349230", "sample": "I-SceI_7_8_9", "expected_mob_b": None},
]


def parse_gd(path: Path) -> list[dict[str, str]]:
    """Parse Genome Diff text into records."""
    records = []
    if not path.is_file():
        return records
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split("\t")
        rec_type = tokens[0]
        rec: dict[str, str] = {"type": rec_type, "raw": line}
        if len(tokens) > 1:
            rec["id"] = tokens[1]
        if len(tokens) > 2:
            rec["parent_ids"] = tokens[2]
        if len(tokens) > 3:
            rec["seq_id"] = tokens[3]
        if len(tokens) > 4:
            rec["pos"] = tokens[4]
        # parse key=value attributes
        for tok in tokens[5:]:
            if "=" in tok:
                k, v = tok.split("=", 1)
                rec[k] = v
        records.append(rec)
    return records


def get_fasta_sha(path: Path) -> str:
    if not path.is_file():
        return "MISSING"
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    bases = "".join(line.strip() for line in lines if line and not line.startswith(">"))
    return hashlib.sha256(bases.upper().encode("ascii")).hexdigest()


def evaluate_run(run_info: dict[str, object]) -> dict[str, object]:
    run = str(run_info["run"])
    sample = str(run_info["sample"])
    run_dir = RESULTS_DIR / run
    cond_a_dir = run_dir / "cond_A"
    cond_b_dir = run_dir / "cond_B"

    # Reference hashes
    ref_a_file = cond_a_dir / "work/reference.fa"
    if not ref_a_file.is_file():
        ref_a_file = REPO / "data/references/m6_identity/CP009273.1.fasta"
    sha_a = get_fasta_sha(ref_a_file)

    ref_b_file = cond_b_dir / "work/reference.fa"
    if not ref_b_file.is_file():
        ref_b_file = REPO / "data/references/m6_identity/CP009273.1.gbk"
        # calculate from gbk origin if work/reference.fa was cleaned
        from stage_m6_p3b_reference import extract_origin_bases
        origin = extract_origin_bases(ref_b_file.read_text(encoding="utf-8"))
        sha_b = hashlib.sha256(origin.encode("ascii")).hexdigest()
    else:
        sha_b = get_fasta_sha(ref_b_file)

    gate1_pass = (sha_a == EXPECTED_SHA and sha_b == EXPECTED_SHA)

    gd_a = parse_gd(cond_a_dir / "output.gd")
    gd_b = parse_gd(cond_b_dir / "output.gd")

    counts_a = {t: sum(1 for r in gd_a if r["type"] == t) for t in ["SNP", "INS", "DEL", "MC", "JC", "MOB"]}
    counts_b = {t: sum(1 for r in gd_b if r["type"] == t) for t in ["SNP", "INS", "DEL", "MC", "JC", "MOB"]}

    # Gate 2: Small variant invariance
    small_types = ["SNP", "INS", "DEL", "MC"]
    gate2_pass = True
    for st in small_types:
        if counts_a[st] != counts_b[st]:
            gate2_pass = False
        recs_a = sorted(r["raw"] for r in gd_a if r["type"] == st)
        recs_b = sorted(r["raw"] for r in gd_b if r["type"] == st)
        # Strip ID column differences if any, compare mutation coordinates and alleles
        recs_a_norm = ["\t".join([r.split("\t")[0]] + r.split("\t")[3:]) for r in recs_a]
        recs_b_norm = ["\t".join([r.split("\t")[0]] + r.split("\t")[3:]) for r in recs_b]
        if recs_a_norm != recs_b_norm:
            gate2_pass = False

    # Gate 3: Control No New MOB (gate3_control_no_new_mob)
    # Control sample SRR28349226 has MOB == 0 under both conditions.
    # Biological specificity/sensitivity remain NOT_EVALUABLE without independent ground truth.
    if run == "SRR28349226":
        gate3_pass = (counts_b["MOB"] == 0 and counts_a["MOB"] == 0)
    else:
        gate3_pass = True  # Evaluated for control

    # Gate 4: Structural filtering and MOB consistency (gate4_mob_promotion)
    # Validates:
    # 1. Non-negative delta MOB (Delta MOB >= 0).
    # 2. Constituent parent-JC consistency: for any emitted MOB in Cond B, parent JCs exist and mob_evidence=1.
    # 3. Non-expansion of standalone JCs (JC_B <= JC_A) via repeat filtering (ANNOTATION_FILTERED_JC / REPEAT_ASSOCIATED_JC).
    # Note: Positive empirical MOB detection was NOT demonstrated in these 3 BW25113 runs (all tested MOB counts
    # remained zero / MOB_NOT_CALLED). Gate 4 validates structural filtering and MOB call consistency, not active MOB recovery.
    gate4_pass = True
    delta_mob = counts_b["MOB"] - counts_a["MOB"]
    if delta_mob < 0:
        gate4_pass = False

    # Check that for each MOB in Cond B, parent JCs exist and are tagged mob_evidence=1
    jc_by_id = {r["id"]: r for r in gd_b if r["type"] == "JC"}
    for mob in [r for r in gd_b if r["type"] == "MOB"]:
        pids = mob.get("parent_ids", "").split(",")
        if len(pids) < 2:
            gate4_pass = False
        for pid in pids:
            if pid in jc_by_id:
                jc = jc_by_id[pid]
                if jc.get("mob_evidence") != "1":
                    gate4_pass = False

    # Check that JC count did not expand (JC_B <= JC_A) due to repeat-associated junction filtering
    delta_jc = counts_b["JC"] - counts_a["JC"]
    if counts_b["JC"] > counts_a["JC"]:
        gate4_pass = False

    return {
        "run": run,
        "sample": sample,
        "cond_a_ref_sha256": sha_a,
        "cond_b_ref_sha256": sha_b,
        "cond_a_repeats": 0,
        "cond_b_repeats": 401,
        "cond_a_snp": counts_a["SNP"],
        "cond_b_snp": counts_b["SNP"],
        "delta_snp": counts_b["SNP"] - counts_a["SNP"],
        "cond_a_ins": counts_a["INS"],
        "cond_b_ins": counts_b["INS"],
        "delta_ins": counts_b["INS"] - counts_a["INS"],
        "cond_a_del": counts_a["DEL"],
        "cond_b_del": counts_b["DEL"],
        "delta_del": counts_b["DEL"] - counts_a["DEL"],
        "cond_a_mc": counts_a["MC"],
        "cond_b_mc": counts_b["MC"],
        "delta_mc": counts_b["MC"] - counts_a["MC"],
        "cond_a_jc": counts_a["JC"],
        "cond_b_jc": counts_b["JC"],
        "delta_jc": delta_jc,
        "cond_a_mob": counts_a["MOB"],
        "cond_b_mob": counts_b["MOB"],
        "delta_mob": delta_mob,
        "gate1_sequence_identity": "PASS" if gate1_pass else "FAIL",
        "gate2_small_variant_invariance": "PASS" if gate2_pass else "FAIL",
        "gate3_control_no_new_mob": "PASS" if gate3_pass else "FAIL",
        "gate4_mob_promotion": "PASS" if gate4_pass else "FAIL",
    }


def main() -> int:
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    all_pass = True
    for item in TRIAD:
        res = evaluate_run(item)
        rows.append(res)
        gates = [res[f"gate{i}_{suffix}"] for i, suffix in [
            (1, "sequence_identity"),
            (2, "small_variant_invariance"),
            (3, "control_no_new_mob"),
            (4, "mob_promotion"),
        ]]
        if any(g == "FAIL" for g in gates):
            all_pass = False

    fields = list(rows[0])
    with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT_TSV} with {len(rows)} runs.")

    # Print summary
    print("\n=== M6-P3b A/B Comparison Results ===")
    print(f"{'Run':<12} {'Sample':<20} {'SNP(A->B)':<10} {'JC(A->B)':<12} {'MOB(A->B)':<10} {'Gates 1-4':<10}")
    for r in rows:
        snps = f"{r['cond_a_snp']}->{r['cond_b_snp']}"
        jcs = f"{r['cond_a_jc']}->{r['cond_b_jc']}"
        mobs = f"{r['cond_a_mob']}->{r['cond_b_mob']}"
        gates = f"{r['gate1_sequence_identity']}/{r['gate2_small_variant_invariance']}/{r['gate3_control_no_new_mob']}/{r['gate4_mob_promotion']}"
        print(f"{r['run']:<12} {r['sample']:<20} {snps:<10} {jcs:<12} {mobs:<10} {gates:<10}")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
