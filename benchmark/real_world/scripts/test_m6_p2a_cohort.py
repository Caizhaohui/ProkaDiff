import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_m6_p2a_cohort import (
    OUT,
    SHA256,
    build_rows,
    classify_run,
    new_execution_runs,
)
from m6_reference_policy import BW25113_EXECUTION, MG1655_EXECUTION, check_manifests


def test_frozen_reference_policy_still_passes():
    check_manifests()


def test_long_average_length_without_mates_is_not_eligible():
    eligibility, structure, reason = classify_run(
        submission="SRA1825115",
        spots=100,
        spots_with_mates=0,
        library_layout="PAIRED",
        identity_ok=True,
        reference_ok=True,
        sra_stat_structure=None,
    )
    assert eligibility == "BLOCKED_READ_STRUCTURE"
    assert structure == "RUNINFO_SPOTS_WITH_MATES_ZERO"
    assert "not inferred" in reason


def test_sra_stat_single_read_stays_blocked_even_if_avg_length_would_have_been_long():
    eligibility, structure, _reason = classify_run(
        submission="SRA1825115",
        spots=1899109,
        spots_with_mates=0,
        library_layout="PAIRED",
        identity_ok=True,
        reference_ok=True,
        sra_stat_structure="SINGLE_BIOLOGICAL_READ_PER_SPOT",
    )
    assert eligibility == "BLOCKED_READ_STRUCTURE"
    assert structure == "SINGLE_BIOLOGICAL_READ_PER_SPOT"


def test_runinfo_mate_match_is_eligible_without_using_length():
    eligibility, structure, reason = classify_run(
        submission="SRA1824580",
        spots=2252350,
        spots_with_mates=2252350,
        library_layout="PAIRED",
        identity_ok=True,
        reference_ok=True,
        sra_stat_structure=None,
    )
    assert eligibility == "ELIGIBLE_PAIRED"
    assert structure == "RUNINFO_SPOTS_WITH_MATES_EQUAL_SPOTS"
    assert reason == "."


def test_identity_failure_blocks_before_mate_counts():
    eligibility, _structure, _reason = classify_run(
        submission="SRA1824580",
        spots=10,
        spots_with_mates=10,
        library_layout="PAIRED",
        identity_ok=False,
        reference_ok=True,
        sra_stat_structure="TRUE_PAIRED_READS_AVAILABLE",
    )
    assert eligibility == "BLOCKED_RUN_IDENTITY"


def test_cohort_keeps_every_sample_and_both_runs():
    rows = build_rows()
    assert len(rows) == 50
    names = {row["sample_name"] for row in rows}
    assert len(names) == 25
    by_elig = {}
    for row in rows:
        by_elig.setdefault(row["execution_eligibility"], []).append(row)
    assert len(by_elig["ELIGIBLE_PAIRED"]) == 25
    assert len(by_elig["BLOCKED_READ_STRUCTURE"]) == 25
    assert set(by_elig) == {"ELIGIBLE_PAIRED", "BLOCKED_READ_STRUCTURE"}
    for row in by_elig["ELIGIBLE_PAIRED"]:
        assert row["submission"] == "SRA1824580"
        assert row["exclusion_reason"] == "."
        assert row["library_layout_reported"] == "PAIRED"
    for row in by_elig["BLOCKED_READ_STRUCTURE"]:
        assert row["submission"] == "SRA1825115"
        assert row["exclusion_reason"] != "."
    verified_single = [row for row in rows if row["observed_read_structure"] == "SINGLE_BIOLOGICAL_READ_PER_SPOT"]
    assert {row["run"] for row in verified_single} == {"SRR28357593", "SRR28357592"}
    verified_paired = [row for row in rows if row["observed_read_structure"] == "TRUE_PAIRED_READS_AVAILABLE"]
    assert {row["run"] for row in verified_paired} == {"SRR28349244", "SRR28349243"}
    unverified_zero = [row for row in rows if row["observed_read_structure"] == "RUNINFO_SPOTS_WITH_MATES_ZERO"]
    assert len(unverified_zero) == 23


def test_strain_specific_reference_and_truth_denominators():
    rows = build_rows()
    cas9 = [row for row in rows if row["sample_name"] == "Cas9_1"]
    assert {row["execution_reference"] for row in cas9} == {MG1655_EXECUTION}
    assert {row["reference_sequence_sha256"] for row in cas9} == {SHA256[MG1655_EXECUTION]}
    assert {row["eligible_publication_event_count"] for row in cas9} == {"1"}
    bw = [row for row in rows if row["sample_name"] == "control_5_6_7_8"]
    assert {row["execution_reference"] for row in bw} == {BW25113_EXECUTION}
    assert {row["reference_sequence_sha256"] for row in bw} == {SHA256[BW25113_EXECUTION]}
    assert "NC_000913.3" not in {row["execution_reference"] for row in bw}
    unresolved = [row for row in rows if row["sample_name"] == "lambdaRed_3_4"]
    assert {row["truth_class"] for row in unresolved} == {"COORDINATE_UNRESOLVED"}
    assert {row["eligible_publication_event_count"] for row in unresolved} == {"0"}
    duplicated = [row for row in rows if row["sample_name"] == "I-SceI_7_8_9"]
    assert {row["eligible_publication_event_count"] for row in duplicated} == {"1"}
    none = [row for row in rows if row["sample_name"] == "Cas9_2"]
    assert {row["truth_class"] for row in none} == {"PUBLICATION_NONE"}
    assert {row["eligible_publication_event_count"] for row in none} == {"0"}
    unique = {}
    for row in rows:
        if row["execution_eligibility"] != "ELIGIBLE_PAIRED":
            continue
        unique[row["sample_name"]] = int(row["eligible_publication_event_count"])
    assert sum(unique.values()) == 12


def test_local_reference_hashes_match_the_frozen_accessions():
    from m6_p2a_stage_gate import gate_run, verify_local_reference
    from m6_reference_policy import BW25113_EXECUTION, MG1655_EXECUTION

    assert verify_local_reference(MG1655_EXECUTION) == SHA256[MG1655_EXECUTION]
    assert verify_local_reference(BW25113_EXECUTION) == SHA256[BW25113_EXECUTION]
    paired = gate_run("SRR28349177")
    assert paired["execution_reference"] == MG1655_EXECUTION
    assert paired["sample_name"] == "Cas9_103"


def test_blocked_run_cannot_pass_the_execution_gate():
    import pytest
    from m6_p2a_stage_gate import gate_run
    from m6_reference_policy import ReferencePolicyError

    with pytest.raises(ReferencePolicyError, match="BLOCKED_READ_STRUCTURE"):
        gate_run("SRR28357554")


def test_new_execution_list_omits_verified_p1b_runs():
    runs = new_execution_runs()
    assert len(runs) == 23
    assert "SRR28349244" not in runs
    assert "SRR28349243" not in runs
    assert len(set(runs)) == 23


def test_versioned_cohort_file_matches_builder():
    generated = build_rows()
    if not OUT.exists():
        pytest.skip("cohort manifest is written by build_m6_p2a_cohort.py")
    with OUT.open(encoding="utf-8", newline="") as handle:
        stored = list(csv.DictReader(handle, delimiter="\t"))
    assert stored == generated
