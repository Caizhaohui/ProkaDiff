#!/usr/bin/env python3
"""Fail closed before an M6-P2a download or evidence job.

The cohort manifest selects the run. The relationship table selects the
execution accession. This script checks both, then hashes the local chromosome.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_m6_p2a_cohort import OUT, SHA256, build_rows  # noqa: E402
from m6_reference_policy import (  # noqa: E402
    BW25113_EXECUTION,
    MG1655_EXECUTION,
    RELATIONSHIPS,
    ReferencePolicyError,
    read_tsv,
    require_execution_reference,
)

REPO = Path(__file__).resolve().parents[3]
REFERENCE_FILES = {
    MG1655_EXECUTION: REPO / "data/references/MG1655.gbk",
    BW25113_EXECUTION: REPO / "data/references/m6_identity/CP009273.1.fasta",
}
REFERENCE_LENGTH = {
    MG1655_EXECUTION: 4_641_652,
    BW25113_EXECUTION: 4_631_469,
}


def chromosome_bases(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    if path.suffix.lower() in {".gbk", ".gb", ".genbank"} or text.startswith("LOCUS"):
        origin = False
        chunks: list[str] = []
        for line in text.splitlines():
            if line.startswith("ORIGIN"):
                origin = True
                continue
            if not origin:
                continue
            if line.startswith("//"):
                break
            chunks.append("".join(ch for ch in line if ch.isalpha()))
        return "".join(chunks).upper()
    chunks = []
    for line in text.splitlines():
        if not line or line.startswith(">"):
            continue
        chunks.append("".join(ch for ch in line if ch.isalpha()))
    return "".join(chunks).upper()


def verify_local_reference(accession: str) -> str:
    path = REFERENCE_FILES.get(accession)
    if path is None or not path.is_file():
        raise ReferencePolicyError(f"local reference for {accession} is missing")
    bases = chromosome_bases(path)
    digest = hashlib.sha256(bases.encode()).hexdigest()
    if len(bases) != REFERENCE_LENGTH[accession] or digest != SHA256[accession]:
        raise ReferencePolicyError(
            f"{accession} sequence mismatch len={len(bases)} sha256={digest} path={path}"
        )
    return digest


def gate_run(run_id: str, rows: list[dict[str, str]] | None = None) -> dict[str, str]:
    rows = rows if rows is not None else build_rows()
    match = [row for row in rows if row["run"] == run_id]
    if len(match) != 1:
        raise ReferencePolicyError(f"{run_id} is not a single row of the P2a cohort")
    row = match[0]
    if row["execution_eligibility"] != "ELIGIBLE_PAIRED":
        raise ReferencePolicyError(
            f"{run_id} is {row['execution_eligibility']}: {row['exclusion_reason']}"
        )
    relationships = read_tsv(RELATIONSHIPS)
    execution = require_execution_reference(run_id, row["execution_reference"], relationships)
    digest = verify_local_reference(execution)
    if digest != row["reference_sequence_sha256"]:
        raise ReferencePolicyError(f"{run_id} cohort hash disagrees with the local sequence")
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description="Check one M6-P2a run before staging or evidence.")
    parser.add_argument("--run", required=True)
    parser.add_argument("--check-references", action="store_true")
    args = parser.parse_args()
    try:
        if args.check_references:
            for accession in (MG1655_EXECUTION, BW25113_EXECUTION):
                digest = verify_local_reference(accession)
                print(f"reference_ok {accession} {digest}")
        row = gate_run(args.run)
    except (ReferencePolicyError, OSError, KeyError) as error:
        print(f"M6-P2a gate failed: {error}", file=sys.stderr)
        return 1
    print(
        "\t".join(
            [
                row["run"],
                row["sample_name"],
                row["execution_reference"],
                row["reference_sequence_sha256"],
                str(REFERENCE_LENGTH[row["execution_reference"]]),
                str(REFERENCE_FILES[row["execution_reference"]]),
            ]
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
