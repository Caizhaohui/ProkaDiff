#!/usr/bin/env python3
"""
benchmark/real_world/scripts/test_prjna884016_quarantine.py

Regression test suite for Milestone M7-P0a:
PRJNA884016 Assembly Truth Quarantine and Fail-Closed Invariants.

Proves:
1. Truncated assemblies cannot generate accepted truth (AssemblyTruncationError).
2. Checksum mismatches fail closed (SequenceChecksumMismatchError).
3. Replicon mismatch fails closed (RepliconMismatchError).
4. MOB family and TSD length are never established from approximate segment length alone (UnsupportedHeuristicError).
5. Empty truth sets yield NOT_EVALUABLE for recall, precision, and breakpoint tolerances.
6. Original erroneous 4 MOB + 7 JC are rejected in active truth assets.
7. Historical invalid assets are preserved with verified checksums in quarantined_legacy/.
8. Registry invariants hold: Line 7 booleans are false, status PENDING; Line 6 M6 untouched.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS_DIR = REPO_ROOT / "benchmark" / "real_world" / "scripts"
TRUTH_DIR = REPO_ROOT / "benchmark" / "real_world" / "truth" / "prjna884016"
QUARANTINED_DIR = TRUTH_DIR / "quarantined_legacy"
REGISTRY_PATH = REPO_ROOT / "benchmark" / "real_world" / "VALIDATION_REGISTRY.tsv"
MANIFEST_PATH = REPO_ROOT / "benchmark" / "real_world" / "manifests" / "prjna884016.tsv"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from build_prjna884016_truth import (
    CANONICAL_LENGTH,
    CANONICAL_SEQ_SHA256,
    EXPECTED_MUT_ACC,
    EXPECTED_REF_ACC,
    AssemblyTruncationError,
    RepliconMismatchError,
    SequenceChecksumMismatchError,
    UnsupportedHeuristicError,
    forbid_heuristic_mob_assignment,
    validate_assembly_fasta,
)
from compare_mob_truth import compare_mob, parse_mob_records


def test_truncated_assembly_fails_closed(tmp_path: Path) -> None:
    """Proves that a truncated FASTA contig raises AssemblyTruncationError."""
    trunc_fasta = tmp_path / "truncated.fna"
    trunc_fasta.write_text(f">{EXPECTED_MUT_ACC} truncated contig\n" + "A" * 4492785 + "\n")
    with pytest.raises(AssemblyTruncationError) as exc_info:
        validate_assembly_fasta(trunc_fasta, EXPECTED_MUT_ACC)
    assert "4492785 bp" in str(exc_info.value)
    assert "cannot generate truth" in str(exc_info.value)


def test_sequence_checksum_mismatch_fails_closed(tmp_path: Path) -> None:
    """Proves that a full-length contig with non-canonical sequence raises SequenceChecksumMismatchError."""
    mutated_fasta = tmp_path / "mutated_full.fna"
    mutated_fasta.write_text(f">{EXPECTED_REF_ACC} non-canonical\n" + "N" * CANONICAL_LENGTH + "\n")
    with pytest.raises(SequenceChecksumMismatchError) as exc_info:
        validate_assembly_fasta(mutated_fasta, EXPECTED_REF_ACC)
    assert "does not match canonical checksum" in str(exc_info.value)


def test_replicon_mismatch_fails_closed(tmp_path: Path) -> None:
    """Proves that a FASTA contig with wrong accession raises RepliconMismatchError."""
    wrong_acc_fasta = tmp_path / "wrong_acc.fna"
    wrong_acc_fasta.write_text(">NC_000913.3 MG1655\n" + "A" * CANONICAL_LENGTH + "\n")
    with pytest.raises(RepliconMismatchError) as exc_info:
        validate_assembly_fasta(wrong_acc_fasta, EXPECTED_REF_ACC)
    assert "does not contain expected accession" in str(exc_info.value)


def test_heuristic_mob_family_and_tsd_rejected() -> None:
    """Proves that heuristic MOB family assignment based on approximate duplication length is rejected."""
    # IS1-like length (~769 bp)
    with pytest.raises(UnsupportedHeuristicError) as exc_is1:
        forbid_heuristic_mob_assignment("DUP", 769)
    assert "strictly prohibited" in str(exc_is1.value)

    # IS10-like length (~1333 bp)
    with pytest.raises(UnsupportedHeuristicError) as exc_is10:
        forbid_heuristic_mob_assignment("DUP", 1333)
    assert "strictly prohibited" in str(exc_is10.value)


def test_empty_truth_scoring_returns_not_evaluable() -> None:
    """Proves that zero positive truth events result in NOT_EVALUABLE rather than 100% or 0%."""
    # Test comparing 0 calls against 0 truth
    metrics_empty = compare_mob([], [])
    assert metrics_empty["tp"] == 0
    assert metrics_empty["fp"] == 0
    assert metrics_empty["fn"] == 0
    assert metrics_empty["precision"] == "NOT_EVALUABLE"
    assert metrics_empty["recall"] == "NOT_EVALUABLE"
    assert metrics_empty["exact_fraction"] == "NOT_EVALUABLE"
    assert metrics_empty["within_5bp_fraction"] == "NOT_EVALUABLE"
    assert metrics_empty["median_error_bp"] == "NOT_EVALUABLE"

    # Test comparing 1 unverified call against 0 truth (1 FP)
    dummy_call = [{"seq_id": "CP106998.1", "pos": 100, "element": "IS1", "strand": "+", "repeat_size": 9}]
    metrics_fp = compare_mob(dummy_call, [])
    assert metrics_fp["tp"] == 0
    assert metrics_fp["fp"] == 1
    assert metrics_fp["fn"] == 0
    assert metrics_fp["precision"] == 0.0
    assert metrics_fp["recall"] == "NOT_EVALUABLE"
    assert metrics_fp["exact_fraction"] == "NOT_EVALUABLE"


def test_legacy_corrupted_truth_records_rejected() -> None:
    """Proves that active truth files reject the 4 MOB + 7 JC erroneous calls."""
    active_gd = TRUTH_DIR / "SY22A1_truth.gd"
    assert active_gd.is_file()
    mob_count = 0
    jc_count = 0
    with active_gd.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if parts[0] == "MOB":
                mob_count += 1
            elif parts[0] == "JC":
                jc_count += 1
    assert mob_count == 0, f"Active GD must have 0 MOB records, found {mob_count}"
    assert jc_count == 0, f"Active GD must have 0 JC records, found {jc_count}"

    # Verify mob_truth.tsv has 0 data rows
    mob_tsv = TRUTH_DIR / "mob_truth.tsv"
    with mob_tsv.open("r", encoding="utf-8") as handle:
        mob_rows = [line for line in handle if line.strip() and not line.startswith("#")]
    assert len(mob_rows) == 1, "mob_truth.tsv must contain only the header line"

    # Verify jc_truth.tsv has 0 data rows
    jc_tsv = TRUTH_DIR / "jc_truth.tsv"
    with jc_tsv.open("r", encoding="utf-8") as handle:
        jc_rows = [line for line in handle if line.strip() and not line.startswith("#")]
    assert len(jc_rows) == 1, "jc_truth.tsv must contain only the header line"


def test_quarantined_legacy_assets_preserved_and_hashes_intact() -> None:
    """Proves that historical invalid assets are preserved with verified checksums."""
    manifest_tsv = QUARANTINED_DIR / "MANIFEST.tsv"
    assert manifest_tsv.is_file(), "quarantined_legacy/MANIFEST.tsv must exist"

    with manifest_tsv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            fname = row["file_name"]
            expected_hash = row["sha256"]
            fpath = QUARANTINED_DIR / fname
            assert fpath.is_file(), f"Quarantined asset {fname} must exist"
            with fpath.open("rb") as fhandle:
                actual_hash = hashlib.sha256(fhandle.read()).hexdigest()
            assert actual_hash == expected_hash, f"Hash mismatch for {fname}"


def test_validation_registry_invariants() -> None:
    """Proves that VALIDATION_REGISTRY.tsv adheres to frozen schema invariants."""
    with REGISTRY_PATH.open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    # Line 6: M6 cohort_25 remains untouched
    m6_row = next(r for r in rows if r["dataset"] == "prjna1088182")
    assert m6_row["sample_id"] == "cohort_25"
    assert m6_row["truth_level"] == "PUBLICATION_REPORTED"
    assert m6_row["variant_validated"] == "false"
    assert m6_row["intended_validated"] == "false"
    assert m6_row["mob_validated"] == "false"
    assert m6_row["jc_validated"] == "false"
    assert m6_row["report_reviewed"] == "false"
    assert m6_row["status"] == "PENDING"

    # Line 7: M7 SY22A1 is quarantined
    m7_row = next(r for r in rows if r["dataset"] == "prjna884016")
    assert m7_row["sample_id"] == "SY22A1"
    assert m7_row["truth_level"] == "ASSEMBLY_CONFIRMED"
    assert m7_row["variant_validated"] == "false"
    assert m7_row["intended_validated"] == "false"
    assert m7_row["mob_validated"] == "false"
    assert m7_row["jc_validated"] == "false"
    assert m7_row["report_reviewed"] == "false"
    assert m7_row["status"] == "PENDING"
    assert "0 chromosomal variants" in m7_row["notes"]
    assert "Historical 4 MOB + 7 JC invalidated" in m7_row["notes"]
    assert "NOT_EVALUABLE" in m7_row["notes"]


def test_manifest_prjna884016_notes() -> None:
    """Proves that manifests/prjna884016.tsv repudiates chromosomal positive transpositions."""
    with MANIFEST_PATH.open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sy22a1 = next(r for r in rows if r["sample_name"] == "SY22A1")
    assert "Negative concordance control" in sy22a1["notes"]
    assert "Historical 4 MOB + 7 JC invalidated" in sy22a1["notes"]
