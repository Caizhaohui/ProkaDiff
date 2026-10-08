import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from m6_reference_policy import (
    BW25113_EXECUTION,
    MG1655_EXECUTION,
    ReferencePolicyError,
    check_manifests,
    read_tsv,
    RELATIONSHIPS,
    require_execution_reference,
)


def test_versioned_manifests_match_frozen_reference_policy():
    check_manifests()


def test_bw25113_request_for_mg1655_reference_fails_closed():
    rows = read_tsv(RELATIONSHIPS)
    bw = next(row for row in rows if "BW25113" in row["sra_organism"])
    with pytest.raises(ReferencePolicyError, match="NC_000913.3"):
        require_execution_reference(bw["run_long"], MG1655_EXECUTION, rows)


def test_mg1655_request_for_bw25113_reference_fails_closed():
    rows = read_tsv(RELATIONSHIPS)
    mg = next(row for row in rows if row["sample_name"] == "Cas9_1")
    with pytest.raises(ReferencePolicyError, match="CP009273.1"):
        require_execution_reference(mg["run_long"], BW25113_EXECUTION, rows)


def test_sra1825115_is_single_biological_read_and_not_initial_paired_smoke():
    import csv
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "manifests/prjna1088182_m6_p1a_fastq_provenance.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    by = {row["run"]: row for row in rows}
    for run in ("SRR28349244", "SRR28349243"):
        assert by[run]["reported_library_layout"] == "PAIRED"
        assert by[run]["read_structure_class"] == "TRUE_PAIRED_READS_AVAILABLE"
        assert by[run]["p1b_initial_paired_smoke"] == "ELIGIBLE_INITIAL_PAIRED_SMOKE"
        assert by[run]["r1_reads"] == by[run]["r2_reads"]
    for run in ("SRR28357593", "SRR28357592"):
        assert by[run]["reported_library_layout"] == "PAIRED"
        assert by[run]["spots_with_mates"] == "0"
        assert by[run]["read_structure_class"] == "SINGLE_BIOLOGICAL_READ_PER_SPOT"
        assert by[run]["p1b_initial_paired_smoke"] == "NOT_ELIGIBLE_INITIAL_PAIRED_SMOKE"
        assert by[run]["r2_reads"] == ""
        assert by[run]["read_structure_class"] != "SECOND_READ_TECHNICAL_ONLY"


def test_unauthorized_run_is_rejected():
    from m6_reference_policy import assert_authorized_smoke_run
    with pytest.raises(ReferencePolicyError, match="not one of the four"):
        assert_authorized_smoke_run("SRR28349176")


def test_smoke_runs_are_independent_mg1655_units():
    rows = read_tsv(RELATIONSHIPS)
    for run_id in ("SRR28349244", "SRR28357593", "SRR28349243", "SRR28357592"):
        assert require_execution_reference(run_id, MG1655_EXECUTION, rows) == MG1655_EXECUTION
