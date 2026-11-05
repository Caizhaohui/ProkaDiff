import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from m6_reference_policy import COORDINATES, read_tsv
from reconcile_m6_p2a import (
    MalformedGdError,
    deduplicated_evaluations,
    eligible_unique_alleles,
    inventory,
    parse_gd,
    classify_publication_event,
    repeat_diagnostic_count,
)

P1B = Path(__file__).resolve().parents[1] / "results/prjna1088182/m6_p1b"


def _event(records, **overrides):
    payload = {
        "coordinate_status": "COORDINATE_VERIFIED",
        "denominator": "include",
        "position": 301690,
        "ref_allele": "C",
        "alt_allele": "T",
        "execution_reference": "NC_000913.3",
        "reference_base": "C",
        "records": records,
    }
    payload.update(overrides)
    return classify_publication_event(**payload)


def test_exact_snp_recovers_and_does_not_use_the_gd_id_as_dv1():
    records = parse_gd("SNP\t3\t.\tNC_000913.3\t301690\tT\tpd_ra_depth=86\n")
    result = _event(records)
    assert result["recovery"] == "RECOVERED_EXACT"
    assert result["gd_record_id"] == "3"
    assert result["dv1_event_id"] == "NOT_EVALUABLE"
    assert result["m4_products"] == "NOT_EVALUABLE_BY_TIER_A"
    assert "301690" in result["observed_event"]
    assert result["observed_event"].split()[0] == "SNP"


def test_wrong_alternate_allele_is_not_recovered():
    records = parse_gd("SNP\t3\t.\tNC_000913.3\t301690\tA\n")
    result = _event(records)
    assert result["recovery"] == "NOT_RECOVERED"
    assert result["diagnostic_reason"] == "UNKNOWN"
    assert "301690 A" in result["observed_event"]


def test_wrong_reference_sequence_is_not_exact_recovery():
    records = parse_gd("SNP\t9\t.\tCP009273.1\t301690\tT\n")
    result = _event(records)
    assert result["recovery"] == "NOT_RECOVERED"
    assert result["dv1_event_id"] == "NOT_EVALUABLE"


def test_reference_base_disagreement_is_a_provenance_stop():
    records = parse_gd("SNP\t3\t.\tNC_000913.3\t301690\tT\n")
    result = _event(records, reference_base="G")
    assert result["recovery"] == "NOT_RECOVERED"
    assert result["diagnostic_reason"] == "PUBLICATION_PROVENANCE_ISSUE"
    assert result["review_stop"] != "."
    assert result["eligible_denominator"] == "include"


def test_duplicated_truth_counts_once():
    row = {
        "sample_name": "I-SceI_7_8_9",
        "reported_position": "2872949",
        "ref_allele": "C",
        "alt_allele": "A",
        "coordinate_status": "COORDINATE_VERIFIED",
        "denominator": "include",
        "duplicate_row": "true",
        "tier_a_eligibility": "ELIGIBLE_EXACT_RECOVERY",
    }
    records = parse_gd("SNP\t4\t.\tCP009273.1\t2872949\tA\n")
    evaluations = deduplicated_evaluations([row, dict(row)], records, "CP009273.1", {2872949: "C"})
    assert len(evaluations) == 1
    assert evaluations[0]["recovery"] == "RECOVERED_EXACT"
    assert len(eligible_unique_alleles([row, dict(row)])) == 1


def test_publication_none_has_no_false_positive_denominator():
    counts = inventory(parse_gd("SNP\t1\t.\tNC_000913.3\t100\tA\nJC\t2\t.\tNC_000913.3\t1\t1\tNC_000913.3\t10\t-1\t0\n"))
    assert "false_positive" not in counts
    assert counts["total_mutation_calls"] == 2
    assert counts["MC_count"] == 0
    assert eligible_unique_alleles([]) == set()


def test_unresolved_insertion_stays_not_evaluable_even_if_an_ins_is_present():
    records = parse_gd("INS\t8\t.\tCP009273.1\t3112908\tAAT\n")
    result = _event(
        records,
        coordinate_status="COORDINATE_UNRESOLVED",
        denominator="exclude",
        position=3112908,
        ref_allele=".",
        alt_allele="+AAT",
        execution_reference="CP009273.1",
        reference_base=None,
    )
    assert result["recovery"] == "NOT_EVALUABLE"
    assert result["eligible_denominator"] == "exclude"


def test_empty_gd_is_not_recovered_without_inventing_a_coverage_reason():
    result = _event(parse_gd("#=GENOME_DIFF\t1.0\n"))
    assert result["recovery"] == "NOT_RECOVERED"
    assert result["diagnostic_reason"] == "UNKNOWN"
    assert result["observed_event"] == "NA"


def test_malformed_gd_is_rejected():
    with pytest.raises(MalformedGdError):
        parse_gd("SNP\tnot-an-id\t.\tNC_000913.3\t301690\tT\n")
    with pytest.raises(MalformedGdError):
        parse_gd("SNP\t1\t.\tNC_000913.3\n")


def test_structural_record_at_the_same_coordinate_is_a_representation_difference():
    records = parse_gd("INS\t4\t.\tNC_000913.3\t301690\tT\n")
    result = _event(records)
    assert result["recovery"] == "REPRESENTATION_DIFFERENCE"
    assert result["diagnostic_reason"] == "REPRESENTATION_DIFFERENCE"
    assert result["recovery"] != "RECOVERED_EXACT"


def test_one_base_sub_is_equivalent_not_a_new_coordinate_rule():
    records = parse_gd("SUB\t5\t.\tNC_000913.3\t301690\t1\tT\n")
    result = _event(records)
    assert result["recovery"] == "RECOVERED_EQUIVALENT"
    assert result["gd_record_id"] == "5"
    assert result["dv1_event_id"] == "NOT_EVALUABLE"


def test_gd_record_id_is_not_read_as_the_publication_coordinate():
    records = parse_gd("SNP\t301690\t.\tNC_000913.3\t1\tT\n")
    result = _event(records)
    assert result["recovery"] == "NOT_RECOVERED"
    assert records[0].gd_record_id == "301690"
    assert records[0].position == 1
    assert records[0].dv1_event_id == "NOT_EVALUABLE"


def test_mc_and_repeat_diagnostics_are_not_mutation_identities():
    text = "\n".join(
        [
            "#=GENOME_DIFF\t1.0",
            "SNP\t1\t.\tNC_000913.3\t10\tA",
            "MC\t2\t.\tNC_000913.3\t15\t40\t0\t0",
            "JC\t3\t.\tNC_000913.3\t1\t1\tNC_000913.3\t4641652\t-1\t0",
            "MOB\t4\t3\tNC_000913.3\t20\trepeat\t-1\t4",
        ]
    )
    counts = inventory(parse_gd(text))
    assert counts["SNP"] == 1
    assert counts["MC_count"] == 1
    assert counts["JC"] == 1
    assert counts["MOB"] == 1
    assert counts["total_mutation_calls"] == 3
    assert repeat_diagnostic_count("sample_id\tx\nrow\n") == 1


def test_telemetry_parser_and_blocked_runs_stay_in_the_summary():
    from summarize_m6_p2a import parse_telemetry, summarize

    parsed = parse_telemetry(
        "\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:35:20\n"
        "\tMaximum resident set size (kbytes): 2918432\n"
        "\tExit status: 0\n"
    )
    assert parsed == {"wall_time": "1:35:20", "max_rss": "2918432K", "exit_code": "0"}
    short = parse_telemetry("\tElapsed (wall clock) time (h:mm:ss or m:ss): 40:02.05\n")
    assert short["wall_time"] == "40:02.05"
    runs, events = summarize()
    assert len(runs) == 50
    blocked = [row for row in runs if row["execution_status"] == "BLOCKED_READ_STRUCTURE"]
    assert len(blocked) == 25
    assert all(row["input_status"] == "NOT_STAGED_FOR_P2A" for row in blocked)
    reused = [row for row in runs if row["run"] == "SRR28349244"]
    assert reused[0]["execution_status"] == "REUSED_VERIFIED_P1B"
    assert reused[0]["recovered_exact"] == "1"
    assert reused[0]["jc_calls"] == "49"
    cas9_2 = [row for row in events if row["run"] == "SRR28349243" and row["publication_event"] == "NONE"]
    assert cas9_2[0]["recovery_classification"] == "NOT_EVALUABLE"
    assert cas9_2[0]["eligible_denominator_flag"] == "exclude"


def test_frozen_audit_has_twelve_unique_eligible_alleles_and_excludes_the_insertion():
    rows = read_tsv(COORDINATES)
    alleles = eligible_unique_alleles(rows)
    assert len(alleles) == 12
    excluded = [row for row in rows if row["denominator"] == "exclude"]
    assert len(excluded) == 1
    assert excluded[0]["sample_name"] == "lambdaRed_3_4"
    assert excluded[0]["alt_allele"] == "+AAT"
    evaluations = deduplicated_evaluations(excluded, [], "CP009273.1", {})
    assert evaluations[0]["recovery"] == "NOT_EVALUABLE"
    assert evaluations[0]["eligible_denominator"] == "exclude"


@pytest.mark.skipif(not (P1B / "SRR28349244/output.gd").exists(), reason="local P1b GD is gitignored")
def test_p1b_cas9_1_gd_still_recovers_301690():
    gd = (P1B / "SRR28349244/output.gd").read_text(encoding="utf-8")
    records = parse_gd(gd)
    counts = inventory(records)
    assert counts["SNP"] == 7
    assert counts["SUB"] == 0
    assert counts["INS"] == 2
    assert counts["DEL"] == 1
    assert counts["JC"] == 49
    assert counts["MOB"] == 3
    assert counts["MC_count"] == 86
    assert counts["total_mutation_calls"] == 62
    result = _event(records)
    assert result["recovery"] == "RECOVERED_EXACT"
    assert result["gd_record_id"] == "3"
    assert result["dv1_event_id"] == "NOT_EVALUABLE"
    rerun = (P1B / "SRR28349244_rerun/output.gd").read_bytes()
    assert rerun == gd.encode()
    neg = parse_gd((P1B / "SRR28349243/output.gd").read_text(encoding="utf-8"))
    neg_counts = inventory(neg)
    assert neg_counts["SNP"] == 6
    assert neg_counts["JC"] == 48
    assert neg_counts["MOB"] == 2
    assert neg_counts["MC_count"] == 88
    assert neg_counts["total_mutation_calls"] == 58
    absent = _event(neg)
    assert absent["recovery"] == "NOT_RECOVERED"
    repeat_rows = repeat_diagnostic_count((P1B / "SRR28349244/repeat_ambiguous_junctions.tsv").read_text(encoding="utf-8"))
    assert repeat_rows > counts["JC"]
