import hashlib
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from stage_m6_p3b_reference import (
    DEST as GBK_PATH,
    EXPECTED_LEN,
    EXPECTED_SHA,
    FASTA_PATH,
    count_repeat_features,
    extract_fasta_bases,
    extract_origin_bases,
)


def test_cp009273_gbk_exists_and_matches_fasta_sequence():
    assert GBK_PATH.exists(), f"Missing {GBK_PATH}"
    text = GBK_PATH.read_text(encoding="utf-8")
    origin_bases = extract_origin_bases(text)
    assert len(origin_bases) == EXPECTED_LEN
    digest = hashlib.sha256(origin_bases.encode("ascii")).hexdigest()
    assert digest == EXPECTED_SHA

    assert FASTA_PATH.exists(), f"Missing {FASTA_PATH}"
    fasta_bases = extract_fasta_bases(FASTA_PATH)
    assert len(fasta_bases) == EXPECTED_LEN
    fasta_digest = hashlib.sha256(fasta_bases.encode("ascii")).hexdigest()
    assert fasta_digest == EXPECTED_SHA

    assert origin_bases == fasta_bases, "ORIGIN sequence does not match FASTA sequence bit-identically"


def test_cp009273_gbk_repeat_annotations():
    assert GBK_PATH.exists()
    text = GBK_PATH.read_text(encoding="utf-8")
    n_rep, n_mob = count_repeat_features(text)
    assert n_rep == 353, f"Expected 353 repeat_region features, found {n_rep}"
    assert n_mob == 48, f"Expected 48 mobile_element features, found {n_mob}"
    assert n_rep + n_mob == 401, f"Expected 401 total repeat features, found {n_rep + n_mob}"

    # Verify version and locus header
    assert "VERSION     CP009273.1" in text
    assert "LOCUS       CP009273" in text


def test_contracts_and_crates_unchanged():
    diff = subprocess.check_output(
        ["git", "diff", "--stat", "--", "crates", "benchmark/real_world/VALIDATION_REGISTRY.tsv"],
        cwd=REPO,
        text=True,
    )
    assert diff.strip() == ""


def test_reference_provenance_manifest():
    prov_file = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3b_reference_provenance.tsv"
    assert prov_file.exists(), f"Missing {prov_file}"
    lines = prov_file.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 2
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    assert header[0] == "accession"
    data = dict(zip(header, row))
    assert data["accession"] == "CP009273.1"
    assert data["length_bp"] == "4631469"
    assert data["sequence_sha256"] == EXPECTED_SHA
    assert data["total_repeats_parsed"] == "401"


def test_ab_comparison_manifest():
    cmp_file = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3b_ab_comparison.tsv"
    assert cmp_file.exists(), f"Missing {cmp_file}"
    lines = cmp_file.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 4
    header = lines[0].split("\t")
    assert "gate3_control_no_new_mob" in header
    assert "gate3_control_specificity" not in header
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))
        assert row["gate1_sequence_identity"] == "PASS"
        assert row["gate2_small_variant_invariance"] == "PASS"
        assert row["gate3_control_no_new_mob"] == "PASS"
        assert row["gate4_mob_promotion"] == "PASS"
