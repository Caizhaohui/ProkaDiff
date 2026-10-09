import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_m6_p2b import OUT, build_rows, mixture_status


def test_mixture_uses_the_spreadsheet_field_not_the_sample_name():
    assert mixture_status("DY1164|n=3|mixed=2 others not in file|pct=32") == "DOCUMENTED_POOLED"
    assert mixture_status("DY1134_fbaB236_C0017492|n=1|mixed=N/A|pct=N/A") == "UNRESOLVED_MIXTURE"
    assert mixture_status("DY1316-1|n=103|mixed=N/A|pct=N/A") == "UNRESOLVED_MIXTURE"


def test_twelve_alleles_and_ten_unresolved_are_accounted_for():
    rows = build_rows()
    eligible = [row for row in rows if row["p2a_recovery"] != "NOT_EVALUABLE"]
    assert len(eligible) == 12
    unresolved = [row for row in rows if row["diagnosis"] == "UNRESOLVED"]
    assert len(unresolved) == 10
    assert all(row["p2a_diagnostic"] == "UNKNOWN" for row in unresolved)
    assert all(row["bam_artifact"] == "ABSENT" for row in unresolved)
    assert all(row["total_read_depth"] == "NOT_MEASURED" for row in unresolved)
    assert all(row["alt_allele_fraction"] == "NOT_MEASURED" for row in unresolved)
    assert all(row["mc_interval_covers_locus"] == "false" for row in unresolved)
    assert all(row["reference_base_matches"] == "true" for row in unresolved)
    recovered = [row for row in rows if row["p2a_recovery"] == "RECOVERED_EXACT"]
    assert {row["published_event"] for row in recovered} == {"301690 C>T", "3536100 G>T"}
    assert all(row["diagnosis"] == "RECOVERED_EXACT" for row in recovered)
    excluded = [row for row in rows if row["published_event"] == "3112908 .>+AAT"]
    assert len(excluded) == 1
    assert excluded[0]["diagnosis"] == "NOT_EVALUABLE"
    assert excluded[0]["total_read_depth"] != "0"


def test_unrecovered_events_keep_sample_and_reference_mapping():
    rows = {row["published_event"]: row for row in build_rows()}
    assert rows["1571257 G>T"]["reference"] == "NC_000913.3"
    assert rows["1571257 G>T"]["sample"] == "Cas9_3 _and2others"
    assert rows["1635887 A>G"]["reference"] == "CP009273.1"
    assert "MG1655" in rows["2776013 T>C"]["strain"]
    assert "BW25113" in rows["2872949 C>A"]["strain"]
    assert rows["2872949 C>A"]["sample"] == "I-SceI_7_8_9"
    pooled = [row for row in rows.values() if row["diagnosis"] == "UNRESOLVED"]
    assert all(row["mixture_status"] == "DOCUMENTED_POOLED" for row in pooled)


def test_gd_ids_are_not_written_as_event_ids_and_registry_is_untouched():
    rows = build_rows()
    for row in rows:
        assert row["diagnosis"] != "3"
        assert "EventId=3" not in row["observation"]
    repo = Path(__file__).resolve().parents[3]
    diff = subprocess.check_output(
        ["git", "diff", "--stat", "--", "crates", "benchmark/real_world/VALIDATION_REGISTRY.tsv"],
        cwd=repo,
        text=True,
    )
    assert diff.strip() == ""


def test_versioned_diagnostic_table_matches_the_builder():
    if not OUT.exists():
        return
    stored = list(csv.DictReader(OUT.open(encoding="utf-8"), delimiter="\t"))
    assert stored == build_rows()
