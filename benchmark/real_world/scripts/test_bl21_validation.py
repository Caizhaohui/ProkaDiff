import csv
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from validate_bl21_products import ProductValidationError, validate_products
from validate_bl21_registry import REPO_ROOT, validate_registry
from compare_intended_edits import evaluate_outcomes, evaluate_peer_outcomes, read_edit_outcomes
from review_target_locus import FIELDS as REVIEW_FIELDS, row_for_gd
from validation_context import ComparisonRole, parse_declaration


EVENT_ID = "DV1_0123456789abcdef0123456789abcdef"
SPACER = "ATTTCGCTGGTGGTCAGATG"


def put(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def output_dir(
    tmp_path: Path,
    status: str = "UNEXPECTED_STRUCTURE",
    *,
    relation: str | None = None,
    matched: str | None = None,
    unexpected: str | None = None,
    gd_type: str = "DEL",
    variant_event_id: str = EVENT_ID,
    link_event_id: str | None = None,
    link_type: str = "JC",
) -> Path:
    if relation is None:
        relation = "NONE" if status == "MISSING" else "UNEXPECTED_AT_TARGET"
    if matched is None:
        matched = "NONE"
    if unexpected is None:
        unexpected = "NONE" if status == "MISSING" else EVENT_ID
    if link_event_id is None:
        link_event_id = variant_event_id
    path = tmp_path / "B21_3_1_BL21_G1"
    path.mkdir(parents=True)
    put(
        path / "summary.txt",
        "editor\tcas9\nintended_provided\tyes\nstructural\t1\nnear_homolog\t0\nscattered_snv\t0\n",
    )
    put(
        path / "report.md",
        "# ProkaDiff Genome Audit Report\n"
        "B21_3_1_1.clean.fq.gz\n"
        f"| Guide spacer | `{SPACER}` |\n"
        "- **Candidate guide-dependent off-target events:** 1\n"
        "Candidate guide sites were found, but no post-edit differential variant was associated within the configured window.\n"
        f"| **edit_1** | `{status}` | `{variant_event_id}` |\n"
        f"**Total post-edit differential variants:** 1\n{variant_event_id}\n",
    )
    put(
        path / "edit_outcomes.tsv",
        "edit_id\tkind\tseq_id\texpected_start\texpected_end\tstatus\tmatched_event_ids\tleft_boundary_status\tright_boundary_status\texpected_size\tobserved_size\tunexpected_events\tnotes\tmatched_event_ids_dv1\tunexpected_event_ids_dv1\n"
        f"edit_1\tdel\tNZ_CP053602.1\t334876\t335735\t{status}\t17\tNA\tNA\t860\t19585\t18\tnone\t{matched}\t{unexpected}\n",
    )
    put(
        path / "post_edit_variants.tsv",
        "variant_id\tseq_id\tposition\tend\tgd_type\tref\talt\torigin_status\tintended_relation\tsize_class\treview_priority\tguide_relation\tguide_mismatches\tpam\tdist_to_site\tmobile_element\tevidence\tgene\tlocus_tag\tfeature_type\tlegacy_class\tevent_id\n"
        f"VAR_0001\tNZ_CP053602.1\t331956\t351540\t{gd_type}\t.\t.\tPOST_EDIT_DIFF\t{relation}\tSTRUCTURAL\tHIGH_ATTENTION\tNONE\tNA\tNA\tNA\tNONE\tRA=NA;MC=0x;JC=385reads\tNA\tNA\tNA\tstructural\t{variant_event_id}\n",
    )
    put(
        path / "unintended.tsv",
        "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\tside2_seq_id\tside2_position\tevent_id\n"
        f"NZ_CP053602.1\t331956\t351540\t{gd_type}\t.\t.\tstructural\tcas9\t\t\t\t\t\t{variant_event_id}\n",
    )
    put(
        path / "offtarget_sites.tsv",
        "site_id\tseq_id\tstart\tend\tstrand\ttarget_seq\tpam\tmismatches\tbulge_type\tbulge_size\tsearch_backend\tcfd_score\thsu_score\n"
        "SITE_000001\tNZ_CP053602.1\t331953\t331957\t+\tACGT\tAGG\t0\tNone\t0\trust_exact\tNA\tNA\n",
    )
    put(
        path / "mutation_offtarget_links.tsv",
        "mutation_id\tsite_id\tmutation_type\tmutation_position\tsite_start\tsite_end\tdistance_to_site\tmismatches\tpam\tcfd_score\tassociation_window\tevent_id\tjunction_side\tsearch_backend\tbulge_type\tbulge_size\n"
        f"mut_17\tSITE_000001\t{link_type}\t331955\t331953\t331957\t0\t0\tAGG\tNA\t50\t{link_event_id}\t{'side_1' if link_type == 'JC' else 'NA'}\trust_exact\tNone\t0\n",
    )
    return path


def test_versioned_registry_validates_private_pairs_and_reference_identity() -> None:
    rows = validate_registry(Path("/hpcfs/fhome/caizhh/19_BL21_edited"))
    assert rows == 3


def test_compact_truth_records_reference_context_without_peer_status() -> None:
    fixture = REPO_ROOT / "benchmark/real_world/fixtures/bl21_compact_truth.tsv"
    with fixture.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert "intended_status" not in rows[0]
    peer_statuses = {"COMPLETE", "PARTIAL", "MISSING", "UNEXPECTED_STRUCTURE"}
    structural = [row for row in rows if row["record_type"] == "DEL"]
    assert {row["reference_context"] for row in structural} == {"CURATED_REFERENCE_STRUCTURE"}
    assert peer_statuses.isdisjoint(row["reference_context"] for row in rows)
    assert all(row["exact_jc"] != "NA" for row in structural)
    with (REPO_ROOT / "benchmark/real_world/truth/bl21_curated_truth.tsv").open(encoding="utf-8", newline="") as handle:
        curated = list(csv.DictReader(handle, delimiter="\t"))
    curated_coordinates = {
        (row["seq_id"], row["start"], row["end"])
        for row in curated
        if row["is_structural"] == "true"
    }
    assert {(row["seq_id"], row["start"], row["end"]) for row in structural} == curated_coordinates
    mc_only = next(row for row in rows if row["case_id"] == "mc_only_no_exact_jc")
    assert mc_only["record_type"] == "MC"
    assert mc_only["event_id"] == "NA"
    assert mc_only["exact_jc"] == "NA"
    assert "JC:" not in mc_only["evidence"]


def test_current_m4_products_accept_aberrant_target_structure(tmp_path: Path) -> None:
    counts = validate_products(output_dir(tmp_path), "B21_3_1", "BL21_G1")
    assert counts["unexpected_event_ids"] == 1
    assert counts["associations"] == 1


def test_current_m4_products_reject_false_complete(tmp_path: Path) -> None:
    with pytest.raises(ProductValidationError, match="unexpected"):
        validate_products(
            output_dir(tmp_path, "COMPLETE", relation="UNEXPECTED_AT_TARGET", matched=EVENT_ID),
            "B21_3_1", "BL21_G1",
        )


def test_peer_comparator_allows_missing_when_no_intended_event_relation_exists(tmp_path: Path) -> None:
    counts = validate_products(output_dir(tmp_path, "MISSING"), "B21_3_1", "BL21_G1", "PEER_COMPARATOR")
    assert counts["matched_event_ids"] == 0
    assert counts["unexpected_event_ids"] == 0


def test_peer_comparator_allows_empty_association_output_for_missing(tmp_path: Path) -> None:
    path = output_dir(tmp_path, "MISSING")
    put(
        path / "mutation_offtarget_links.tsv",
        "mutation_id\tsite_id\tmutation_type\tmutation_position\tsite_start\tsite_end\tdistance_to_site\tmismatches\tpam\tcfd_score\tassociation_window\tevent_id\tjunction_side\tsearch_backend\tbulge_type\tbulge_size\n",
    )
    report = path / "report.md"
    put(report, report.read_text(encoding="utf-8").replace("Candidate guide-dependent off-target events:** 1", "Candidate guide-dependent off-target events:** 0"))
    counts = validate_products(path, "B21_3_1", "BL21_G1", "PEER_COMPARATOR")
    assert counts["associations"] == 0


def test_peer_comparator_missing_rejects_unexpected_event_reference(tmp_path: Path) -> None:
    path = output_dir(tmp_path, "MISSING")
    outcomes = path / "edit_outcomes.tsv"
    text = outcomes.read_text(encoding="utf-8").replace("\tNONE\tNONE\n", f"\tNONE\t{EVENT_ID}\n")
    put(outcomes, text)
    with pytest.raises(ProductValidationError, match="MISSING status must not claim"):
        validate_products(path, "B21_3_1", "BL21_G1", "PEER_COMPARATOR")


def test_intended_comparator_reads_current_edit_outcomes_product(tmp_path: Path) -> None:
    path = tmp_path / "edit_outcomes.tsv"
    put(path, "edit_id\tseq_id\texpected_start\texpected_end\tstatus\nedit_1\tNZ_CP053602.1\t334876\t335735\tUNEXPECTED_STRUCTURE\n")
    rows = read_edit_outcomes(path)
    metrics = evaluate_outcomes(rows, "UNEXPECTED_STRUCTURE")
    assert metrics["false_complete"] == 0
    assert metrics["status_mismatch"] == 0
    assert metrics["unexpected_structure"] == 1


def test_peer_comparator_metrics_accept_supported_complete() -> None:
    missing = evaluate_peer_outcomes([{"status": "MISSING", "matched_event_ids_dv1": "NONE", "unexpected_event_ids_dv1": "NONE"}])
    unsupported = evaluate_peer_outcomes([{"status": "COMPLETE", "matched_event_ids_dv1": "NONE", "unexpected_event_ids_dv1": "NONE"}])
    supported = evaluate_peer_outcomes([
        {"status": "COMPLETE", "matched_event_ids_dv1": EVENT_ID, "unexpected_event_ids_dv1": "NONE"},
    ])
    assert missing["false_complete"] == 0 and missing["status_mismatch"] == 0
    assert unsupported["false_complete"] == 1
    assert supported["false_complete"] == 0 and supported["status_mismatch"] == 0


def test_evidence_only_mc_cannot_become_a_variant_or_event(tmp_path: Path) -> None:
    path = output_dir(tmp_path)
    variants = path / "post_edit_variants.tsv"
    text = variants.read_text(encoding="utf-8").replace("\tDEL\t", "\tMC\t").replace(f"\t{EVENT_ID}\n", "\tNA\n")
    put(variants, text)
    with pytest.raises(ProductValidationError, match="evidence-only MC"):
        validate_products(path, "B21_3_1", "BL21_G1")


def test_legacy_mutation_id_cannot_be_used_as_public_event_id(tmp_path: Path) -> None:
    path = output_dir(tmp_path)
    variants = path / "post_edit_variants.tsv"
    put(variants, variants.read_text(encoding="utf-8").replace(EVENT_ID, "mut_17"))
    with pytest.raises(ProductValidationError, match="invalid DV1 EventId"):
        validate_products(path, "B21_3_1", "BL21_G1")


def test_association_measurements_must_match_authoritative_site(tmp_path: Path) -> None:
    path = output_dir(tmp_path)
    links = path / "mutation_offtarget_links.tsv"
    put(links, links.read_text(encoding="utf-8").replace("\t0\tAGG\tNA", "\t1\tAGG\tNA"))
    with pytest.raises(ProductValidationError, match="mismatches does not match"):
        validate_products(path, "B21_3_1", "BL21_G1")


def test_report_cannot_claim_unsupported_exact_breakpoint(tmp_path: Path) -> None:
    path = output_dir(tmp_path)
    report = path / "report.md"
    put(report, report.read_text(encoding="utf-8") + "\nAn exact breakpoint was identified.\n")
    with pytest.raises(ProductValidationError, match="exact breakpoint"):
        validate_products(path, "B21_3_1", "BL21_G1")


def test_one_guide_output_cannot_be_labeled_as_other_declared_guide(tmp_path: Path) -> None:
    with pytest.raises(ProductValidationError, match="guide identity mismatch"):
        validate_products(output_dir(tmp_path), "B21_3_1", "BL21_G2")


def test_comparison_contexts_are_explicit_and_do_not_infer_parenthood() -> None:
    reference = parse_declaration({"comparison_role": "REFERENCE_ONLY", "reference": "GCF_X", "biological_parent_verified": "false"})
    parent = parse_declaration({"comparison_role": "MATCHED_PARENT", "starter_id": "parent", "biological_parent_verified": "true"})
    peer = parse_declaration({"comparison_role": "PEER_COMPARATOR", "starter_id": "clone", "biological_parent_verified": "false"})
    assert reference.role is ComparisonRole.REFERENCE_ONLY and not reference.has_differential_product
    assert parent.role is ComparisonRole.MATCHED_PARENT and parent.has_differential_product
    assert peer.role is ComparisonRole.PEER_COMPARATOR and not peer.biological_parent_verified
    with pytest.raises(ValueError, match="verified biological parent"):
        parse_declaration({"comparison_role": "MATCHED_PARENT", "starter_id": "parent"})


def test_single_sample_mc_only_review_does_not_create_event_identity(tmp_path: Path) -> None:
    gd = tmp_path / "evidence.gd"
    put(gd, "#=GENOME_DIFF\t1.0\nMC\t903\t.\tNZ_CP053602.1\t331956\t352202\t0\t0\n")
    row = row_for_gd("B21_3_1", "edit_1", "PROKADIFF", gd, "NZ_CP053602.1", 334876, 335735)
    assert row["state"] == "ABNORMAL_TARGET_EVIDENCE"
    assert row["evidence_level"] == "COVERAGE"
    assert row["exact_breakpoint_supported"] == "false"
    assert row["evidence_source"] == "PROKADIFF"
    assert set(row) == set(REVIEW_FIELDS)
    assert "EventId" not in row and "event_id" not in row


def test_single_sample_review_requires_supported_jc_for_exact_breakpoint(tmp_path: Path) -> None:
    gd = tmp_path / "evidence.gd"
    put(gd, "MC\t1\t.\tseq\t100\t200\t0\t0\nJC\t2\t.\tseq\t99\t-1\tseq\t201\t1\t0\n")
    unsupported = row_for_gd("sample", "edit", "PROKADIFF", gd, "seq", 120, 180)
    assert unsupported["state"] == "ABNORMAL_TARGET_EVIDENCE"
    assert unsupported["exact_breakpoint_supported"] == "false"
    put(gd, "JC\t2\t.\tseq\t99\t-1\tseq\t201\t1\t0\n")
    unresolved = row_for_gd("sample", "edit", "PROKADIFF", gd, "seq", 120, 180)
    assert unresolved["state"] == "INSUFFICIENT_EVIDENCE"
    put(gd, "JC\t2\t.\tseq\t99\t-1\tseq\t201\t1\t0\tpd_support_reads=18\n")
    supported = row_for_gd("sample", "edit", "PROKADIFF", gd, "seq", 120, 180)
    assert supported["state"] == "EXACT_STRUCTURE_SUPPORTED"
    assert supported["evidence_level"] == "BREAKPOINT"
    assert supported["exact_breakpoint_supported"] == "true"


def test_single_sample_review_ignores_distal_and_circular_junctions(tmp_path: Path) -> None:
    gd = tmp_path / "evidence.gd"
    put(
        gd,
        "#=GENOME_DIFF\t1.0\n"
        "MC\t903\t.\tNZ_CP053602.1\t331956\t352202\t0\t0\n"
        "JC\t970\t.\tNZ_CP053602.1\t1\t1\tNZ_CP053602.1\t4560251\t-1\t0\tpd_support_reads=400\n"
        "JC\t978\t.\tNZ_CP053602.1\t175961\t-1\tNZ_CP053602.1\t441429\t1\t0\tpd_support_reads=115\n",
    )
    row = row_for_gd("B21_3_1", "edit_1", "PROKADIFF", gd, "NZ_CP053602.1", 334876, 335735)
    assert row["state"] == "ABNORMAL_TARGET_EVIDENCE"
    assert row["evidence_level"] == "COVERAGE"
    assert row["exact_breakpoint_supported"] == "false"
    assert "970" not in row["observed_records"]
    assert "978" not in row["observed_records"]
    assert "903" in row["observed_records"]


FAKE_EVENT = "DV1_ffffffffffffffffffffffffffffffff"
COMPANION_EVENT = "DV1_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"


def abnormal_review(path: Path) -> Path:
    path.write_text(
        "sample_id\ttarget_id\tstate\tevidence_level\tevidence_source\tobserved_records\texact_breakpoint_supported\tnotes\n"
        "B21_3_1\tedit_1\tABNORMAL_TARGET_EVIDENCE\tCOVERAGE\tPROKADIFF\tMC\tfalse\tcoverage loss only\n",
        encoding="utf-8",
    )
    return path


def test_complete_passes_when_differential_evidence_supports_the_expected_edit(tmp_path: Path) -> None:
    products = output_dir(tmp_path, "COMPLETE", relation="EXPECTED", matched=EVENT_ID, unexpected="NONE")
    counts = validate_products(products, "B21_3_1", "BL21_G1", "PEER_COMPARATOR", abnormal_review(tmp_path / "review.tsv"))
    assert counts["matched_event_ids"] == 1
    assert counts["unexpected_event_ids"] == 0
    part = output_dir(tmp_path / "part", "COMPLETE", relation="PART_OF_EXPECTED", matched=EVENT_ID, unexpected="NONE")
    assert validate_products(part, "B21_3_1", "BL21_G1")["matched_event_ids"] == 1


def test_complete_fails_without_matched_event_or_with_contradictory_relationship(tmp_path: Path) -> None:
    with pytest.raises(ProductValidationError, match="lacks a matched"):
        validate_products(
            output_dir(tmp_path / "none", "COMPLETE", relation="NONE", matched="NONE", unexpected="NONE"),
            "B21_3_1", "BL21_G1",
        )
    with pytest.raises(ProductValidationError, match="incompatible relationship"):
        validate_products(
            output_dir(tmp_path / "relation", "COMPLETE", relation="UNEXPECTED_AT_TARGET", matched=EVENT_ID, unexpected="NONE"),
            "B21_3_1", "BL21_G1",
        )


def test_missing_rejects_expected_and_part_of_expected_relationships(tmp_path: Path) -> None:
    validate_products(output_dir(tmp_path / "clear", "MISSING"), "B21_3_1", "BL21_G1", "PEER_COMPARATOR")
    with pytest.raises(ProductValidationError, match="intended-related"):
        validate_products(output_dir(tmp_path / "expected", "MISSING", relation="EXPECTED"), "B21_3_1", "BL21_G1")
    with pytest.raises(ProductValidationError, match="intended-related"):
        validate_products(output_dir(tmp_path / "part", "MISSING", relation="PART_OF_EXPECTED"), "B21_3_1", "BL21_G1")


def test_event_ids_resolve_through_the_differential_index_not_the_report(tmp_path: Path) -> None:
    resolved = validate_products(output_dir(tmp_path / "matched", "COMPLETE", relation="EXPECTED", matched=EVENT_ID, unexpected="NONE"), "B21_3_1", "BL21_G1")
    assert resolved["matched_event_ids"] == 1
    unexpected = validate_products(output_dir(tmp_path / "unexpected"), "B21_3_1", "BL21_G1")
    assert unexpected["unexpected_event_ids"] == 1
    missing_id = output_dir(tmp_path / "missing_id", "UNEXPECTED_STRUCTURE", matched=FAKE_EVENT)
    report = missing_id / "report.md"
    report.write_text(report.read_text(encoding="utf-8") + FAKE_EVENT + "\n", encoding="utf-8")
    with pytest.raises(ProductValidationError, match="does not resolve"):
        validate_products(missing_id, "B21_3_1", "BL21_G1")
    bad_link = output_dir(tmp_path / "bad_link", link_event_id=FAKE_EVENT)
    report = bad_link / "report.md"
    report.write_text(report.read_text(encoding="utf-8") + FAKE_EVENT + "\n", encoding="utf-8")
    with pytest.raises(ProductValidationError, match="does not resolve"):
        validate_products(bad_link, "B21_3_1", "BL21_G1")


def test_hidden_mob_companion_event_resolves_without_a_variant_row(tmp_path: Path) -> None:
    products = output_dir(tmp_path, "MISSING", gd_type="MOB", link_event_id=COMPANION_EVENT)
    counts = validate_products(products, "B21_3_1", "BL21_G1", "PEER_COMPARATOR")
    assert counts["hidden_mob_companions"] == 1
    assert COMPANION_EVENT not in (products / "post_edit_variants.tsv").read_text(encoding="utf-8")
    assert COMPANION_EVENT not in (products / "report.md").read_text(encoding="utf-8")


def test_evidence_only_ra_cannot_become_a_variant_or_event(tmp_path: Path) -> None:
    path = output_dir(tmp_path)
    variants = path / "post_edit_variants.tsv"
    text = variants.read_text(encoding="utf-8").replace("\tDEL\t", "\tRA\t").replace(f"\t{EVENT_ID}\n", "\tNA\n")
    put(variants, text)
    with pytest.raises(ProductValidationError, match="evidence-only RA"):
        validate_products(path, "B21_3_1", "BL21_G1")


def test_repeat_ambiguous_sidecar_cannot_supply_event_identity(tmp_path: Path) -> None:
    products = output_dir(tmp_path, "MISSING")
    sidecar = products / "repeat_ambiguous_junctions.tsv"
    put(sidecar, "sample_id\tstate\texact_breakpoint_supported\tanchor_position_1\nB21_3_1\tSUPPORTED_AMBIGUOUS\tfalse\t331955\n")
    counts = validate_products(products, "B21_3_1", "BL21_G1", "PEER_COMPARATOR")
    assert counts["hidden_mob_companions"] == 0
    put(sidecar, f"sample_id\tstate\tevent_id\nB21_3_1\tSUPPORTED_AMBIGUOUS\t{FAKE_EVENT}\n")
    with pytest.raises(ProductValidationError, match="EventId identity"):
        validate_products(products, "B21_3_1", "BL21_G1")


def test_external_breakpoint_wording_is_not_attributed_to_prokadiff(tmp_path: Path) -> None:
    products = output_dir(tmp_path, "MISSING")
    report = products / "report.md"
    text = report.read_text(encoding="utf-8")
    text += "\nbreseq reports an exact breakpoint.\n"
    text += "curated interpretation suggests IS1A involvement.\n"
    put(report, text)
    validate_products(products, "B21_3_1", "BL21_G1", "PEER_COMPARATOR")
    for sentence in (
        "ProkaDiff called the exact repeat-copy breakpoint.",
        "The exact IS1A copy was resolved.",
        "This is an editing-induced structural event.",
        "Cas9 caused this event.",
        "The call is a verified parent-derived event.",
    ):
        put(report, text + sentence + "\n")
        with pytest.raises(ProductValidationError, match="unsupported claim"):
            validate_products(products, "B21_3_1", "BL21_G1")


def test_one_guide_failure_still_attempts_the_other_guide(tmp_path: Path) -> None:
    log = tmp_path / "attempts"
    runner = tmp_path / "runner.sh"
    runner.write_text("#!/bin/bash\nprintf '%s\\n' \"$1\" >> \"$ATTEMPT_LOG\"\n[[ \"$1\" != \"$FAIL_GUIDE\" ]]\n", encoding="utf-8")
    runner.chmod(0o755)
    helper = REPO_ROOT / "benchmark/real_world/scripts/attempt_peer_guides.sh"
    for failed in ("BL21_G1", "BL21_G2"):
        log.write_text("", encoding="utf-8")
        completed = subprocess.run(
            ["bash", "-c", f'source "{helper}"; export ATTEMPT_LOG="{log}"; export FAIL_GUIDE="{failed}"; attempt_peer_guides "{runner}"; printf "%s" "$ATTEMPTED_GUIDES"'],
            check=False, capture_output=True, text=True,
        )
        assert completed.returncode == 0, completed.stderr
        assert completed.stdout == "BL21_G1 BL21_G2"
        assert log.read_text(encoding="utf-8").split() == ["BL21_G1", "BL21_G2"]


def test_job_scripts_name_columns_and_require_the_manifest_comparator() -> None:
    matrix = (REPO_ROOT / "benchmark/real_world/scripts/run_bl21_validation.sbatch").read_text(encoding="utf-8")
    focused = (REPO_ROOT / "benchmark/real_world/scripts/run_p1_3_diff.sbatch").read_text(encoding="utf-8")
    assert 'manifest_value "$sample" spacer_1' in matrix
    assert 'manifest_value "$sample" spacer_2' in matrix
    assert 'manifest_value "$sample" 11' not in matrix
    assert 'manifest_value "$sample" 13' not in matrix
    assert 'manifest_value "$EDITED" starter_id' in focused
    assert 'manifest starter_id' in focused
    assert "attempt_peer_guides run_peer_guide" in matrix
    assert "attempt_peer_guides run_peer_guide" in focused
