#!/usr/bin/env python3
"""
benchmark/real_world/scripts/build_prjna884016_truth.py

QUARANTINED AND DISABLED FOR AUTHORITATIVE VALIDATION (M7-P0a).

This generator was previously used to construct positive MOB/JC truth for SY22A1.
Forensic audit in M7-P0 proved that:
1. Complete circular assemblies CP107005.1 and CP106998.1 are 100.000% bit-identical
   across all 4,687,476 bp (0 chromosomal mutations).
2. The historical 4 MOB + 7 JC calls were artifacts of a 194.7 kb local file truncation
   (4,492,785 bp vs expected 4,687,476 bp).
3. Transposition events in Sheng et al. (2023) occurred exclusively on the cas-expression
   plasmid, which was never deposited in NCBI.

This script now enforces fail-closed integrity checks and refuses to generate
unsupported positive structural truth from approximate duplication lengths or truncated contigs.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Final

CANONICAL_LENGTH: Final[int] = 4687476
CANONICAL_SEQ_SHA256: Final[str] = "524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a"
EXPECTED_REF_ACC: Final[str] = "CP106998.1"
EXPECTED_MUT_ACC: Final[str] = "CP107005.1"


class TruthGeneratorError(Exception):
    """Base error for truth generator integrity violations."""


class AssemblyTruncationError(TruthGeneratorError):
    """Raised when an assembly contig is truncated or missing terminal sequences."""


class SequenceChecksumMismatchError(TruthGeneratorError):
    """Raised when sequence checksum does not match canonical NCBI record."""


class RepliconMismatchError(TruthGeneratorError):
    """Raised when replicon header does not match expected accession."""


class UnsupportedHeuristicError(TruthGeneratorError):
    """Raised when an algorithm attempts to infer IS family or TSD from segment length."""


class PositiveTruthProhibitedError(TruthGeneratorError):
    """Raised when attempting to emit positive variants for 0-difference assemblies."""


def validate_assembly_fasta(path: Path, expected_acc: str) -> tuple[str, str, int]:
    """
    Validates assembly FASTA completeness, replicon accession, and normalized checksum.
    Returns (header, sequence, length).
    """
    if not path.is_file():
        raise FileNotFoundError(f"Assembly FASTA not found: {path}")

    with path.open("r", encoding="ascii", errors="replace") as handle:
        header = handle.readline().strip()
        if not header.startswith(">"):
            raise RepliconMismatchError(f"Missing FASTA header in {path}")
        if expected_acc not in header:
            raise RepliconMismatchError(f"Header '{header}' does not contain expected accession '{expected_acc}'")

        seq = "".join(line.strip() for line in handle if not line.startswith(">"))

    length = len(seq)
    if length != CANONICAL_LENGTH:
        raise AssemblyTruncationError(
            f"Contig length {length} bp in {path} does not match canonical length {CANONICAL_LENGTH} bp. "
            f"Missing {CANONICAL_LENGTH - length} bp terminal sequence; truncated assemblies cannot generate truth."
        )

    seq_hash = hashlib.sha256(seq.upper().encode("ascii")).hexdigest()
    if seq_hash != CANONICAL_SEQ_SHA256:
        raise SequenceChecksumMismatchError(
            f"Sequence SHA256 {seq_hash} in {path} does not match canonical checksum {CANONICAL_SEQ_SHA256}."
        )

    return header, seq, length


def forbid_heuristic_mob_assignment(feature: str, length: int) -> None:
    """
    Fail-closed guard: No inferred event may be promoted to validated truth solely
    from approximate segment length.
    """
    if feature == "DUP":
        if 700 <= length <= 850 or 1200 <= length <= 1450:
            raise UnsupportedHeuristicError(
                f"Attempted to assign MOB family from approximate segment length ({length} bp). "
                "Inferring IS family and fixed TSD length solely from approximate duplication length is strictly prohibited."
            )


def verify_negative_truth(ref_seq: str, mut_seq: str) -> None:
    """
    Verifies that the assemblies are bit-identical and forbids emitting positive truth.
    """
    if ref_seq == mut_seq:
        # Assemblies are 100.000% identical: 0 chromosomal variants
        return
    raise TruthGeneratorError("Unexpected sequence mismatch between canonical reference and mutant assemblies.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fail-closed truth integrity validator for PRJNA884016 (M7-P0a quarantine)"
    )
    parser.add_argument(
        "--check-integrity",
        action="store_true",
        help="Run fail-closed integrity checks on assembly files",
    )
    parser.add_argument(
        "--ref-fasta",
        type=Path,
        default=Path("data/prjna884016/assemblies/SY221_ref.fna"),
        help="Path to reference FASTA (CP106998.1)",
    )
    parser.add_argument(
        "--mut-fasta",
        type=Path,
        default=Path("data/prjna884016/assemblies/SY22A1.fna"),
        help="Path to mutant FASTA (CP107005.1)",
    )
    args = parser.parse_args()

    # Fail closed by default if invoked without explicit check
    if not args.check_integrity:
        sys.stderr.write(
            "ERROR: build_prjna884016_truth.py is QUARANTINED and DISABLED for authoritative validation.\n"
            "Positive structural truth (4 MOB + 7 JC) is invalidated as a local file truncation artifact.\n"
            "Official chromosome assemblies CP107005.1 and CP106998.1 are 100.000% bit-identical.\n"
            "Use --check-integrity to run fail-closed verification.\n"
        )
        sys.exit(1)

    print("Running fail-closed assembly integrity checks...")
    try:
        ref_header, ref_seq, ref_len = validate_assembly_fasta(args.ref_fasta, EXPECTED_REF_ACC)
        print(f"Reference: {ref_header[:50]}... ({ref_len} bp) [PASS]")
    except Exception as exc:
        print(f"Reference validation failed: {exc}")
        sys.exit(2)

    try:
        mut_header, mut_seq, mut_len = validate_assembly_fasta(args.mut_fasta, EXPECTED_MUT_ACC)
        print(f"Mutant: {mut_header[:50]}... ({mut_len} bp) [PASS]")
    except AssemblyTruncationError as exc:
        print(f"FAIL-CLOSED: Truncated mutant assembly rejected: {exc}")
        sys.exit(3)
    except Exception as exc:
        print(f"Mutant validation failed: {exc}")
        sys.exit(4)

    verify_negative_truth(ref_seq, mut_seq)
    print("Zero-difference chromosomal identity verified. Positive truth generation is prohibited.")


if __name__ == "__main__":
    main()
