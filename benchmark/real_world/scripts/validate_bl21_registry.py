#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path
from typing import Final

from validation_context import ComparisonRole, parse_declaration


REPO_ROOT: Final = Path(__file__).resolve().parents[3]
MANIFEST: Final = REPO_ROOT / "benchmark/real_world/manifests/bl21_user.tsv"
SAMPLES: Final = {"B21_2_1", "B21_3_1", "B21_4_1"}
SPACERS: Final = {"BL21_G1": "ATTTCGCTGGTGGTCAGATG", "BL21_G2": "AATATCGGTGGCCGTGGTGT"}
INTENDED: Final = ("NZ_CP053602.1", "334876", "335735", "del")
STARTERS: Final = {"B21_2_1": "B21_3_1", "B21_3_1": "B21_4_1", "B21_4_1": "B21_3_1"}
TRUTH_COORDS: Final = {
    "TRUTH_BL21_3_1_DEL": ("NZ_CP053602.1", "331956", "351540"),
    "TRUTH_BL21_4_1_DEL": ("NZ_CP053602.1", "307992", "351540"),
    "TRUTH_BL21_2_1_DEL": ("NZ_CP053602.1", "307062", "352239"),
}


class RegistryError(Exception):
    pass


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RegistryError(f"registry has no header: {path}")
        return reader.fieldnames, [dict(row) for row in reader]


def fasta_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current: str | None = None
    with path.open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in sequences:
                    raise RegistryError(f"duplicate FASTA contig {current} in {path}")
                sequences[current] = ""
            elif current is not None:
                sequences[current] += "".join(line.split()).upper()
            elif line.strip():
                raise RegistryError(f"sequence before FASTA header in {path}")
    if not sequences:
        raise RegistryError(f"FASTA has no contigs: {path}")
    return sequences


def genbank_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current: str | None = None
    in_origin = False
    with path.open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith("VERSION"):
                fields = line.split()
                if len(fields) < 2:
                    raise RegistryError(f"malformed GenBank VERSION line in {path}")
                current = fields[1]
            elif line.startswith("ORIGIN"):
                if current is None or current in sequences:
                    raise RegistryError(f"missing or duplicate GenBank VERSION before ORIGIN in {path}")
                sequences[current] = ""
                in_origin = True
            elif line.startswith("//"):
                in_origin = False
                current = None
            elif in_origin and current is not None:
                sequences[current] += "".join(re.findall(r"[A-Za-z]", line)).upper()
    if not sequences:
        raise RegistryError(f"GenBank has no VERSION/ORIGIN records: {path}")
    return sequences


def sequence_digests(sequences: dict[str, str]) -> dict[str, tuple[int, str]]:
    return {
        name: (len(sequence), hashlib.sha256(sequence.encode("ascii")).hexdigest())
        for name, sequence in sequences.items()
    }


def validate_registry(data_root: Path) -> int:
    fields, rows = read_tsv(MANIFEST)
    required = {
        "registry_version", "dataset_id", "sample_id", "starter_id", "r1", "r2",
        "reference", "intended_tsv", "editor", "spacer_1_id", "spacer_1",
        "spacer_2_id", "spacer_2", "pam", "truth_source", "sample_type", "notes",
        "comparison_role", "biological_parent_verified",
    }
    missing = required.difference(fields)
    if missing:
        raise RegistryError(f"BL21 registry schema lacks columns: {', '.join(sorted(missing))}")
    if len(rows) != 3:
        raise RegistryError(f"BL21 registry must contain three sample rows, found {len(rows)}")
    samples = {row["sample_id"] for row in rows}
    if samples != SAMPLES:
        raise RegistryError(f"sample name mismatch: expected {sorted(SAMPLES)}, found {sorted(samples)}")

    fastqs: list[Path] = []
    references: set[Path] = set()
    intended_paths: set[Path] = set()
    for row in rows:
        sample = row["sample_id"]
        if row["registry_version"] != "BL21_INPUTS_V2":
            raise RegistryError(f"{sample}: unsupported registry_version {row['registry_version']!r}")
        try:
            declaration = parse_declaration(row)
        except (KeyError, ValueError) as error:
            raise RegistryError(f"{sample}: invalid comparison declaration: {error}") from error
        if declaration.role is not ComparisonRole.PEER_COMPARATOR or declaration.biological_parent_verified:
            raise RegistryError(f"{sample}: BL21 clone comparison must be PEER_COMPARATOR with unverified parent relationship")
        if row["dataset_id"] != "bl21_user" or row["editor"] != "cas9" or row["pam"] != "NGG":
            raise RegistryError(f"{sample}: dataset/editor/PAM declaration mismatch")
        if row["starter_id"] != STARTERS[sample]:
            raise RegistryError(f"{sample}: starter sample mismatch; expected {STARTERS[sample]}, found {row['starter_id']!r}")
        if row["sample_type"] != "single_clone":
            raise RegistryError(f"{sample}: expected sample_type=single_clone")
        if row["truth_source"] != "MANUALLY_CURATED":
            raise RegistryError(f"{sample}: expected manually curated truth source")
        for mate, key, number in (("R1", "r1", 1), ("R2", "r2", 2)):
            path = Path(row[key])
            suffix = f"/{sample}/{sample}_{number}.clean.fq.gz"
            if not str(path).endswith(suffix):
                raise RegistryError(f"{sample}: {mate} path does not match sample/mate declaration: {path}")
            if not path.is_file():
                raise RegistryError(f"{sample}: missing FASTQ {mate} pair member: {path}")
            fastqs.append(path)

        expected_ref = data_root / "BL21/ncbi_dataset/data/GCF_013167015.1/GCF_013167015.1_ASM1316701v1_genomic.fna"
        reference = Path(row["reference"])
        if "GCF_013167015.1" not in str(reference) or reference.resolve() != expected_ref.resolve():
            raise RegistryError(f"{sample}: reference mismatch; expected {expected_ref}, found {reference}")
        if not reference.is_file():
            raise RegistryError(f"{sample}: missing reference FASTA: {reference}")
        references.add(reference)
        intended_paths.add((REPO_ROOT / row["intended_tsv"]).resolve())

        for id_key, sequence_key, expected_id in (
            ("spacer_1_id", "spacer_1", "BL21_G1"),
            ("spacer_2_id", "spacer_2", "BL21_G2"),
        ):
            sequence = row[sequence_key].upper()
            if row[id_key] != expected_id or sequence != SPACERS[expected_id]:
                raise RegistryError(f"{sample}: incomplete or mismatched {expected_id} spacer declaration")
            if len(sequence) != 20 or re.fullmatch(r"[ACGT]{20}", sequence) is None:
                raise RegistryError(f"{sample}: {expected_id} must be an explicit 20-nt DNA spacer")
    if len(set(fastqs)) != 6:
        raise RegistryError("BL21 FASTQ registry must resolve to six distinct files (three complete R1/R2 pairs)")
    if len(references) != 1 or len(intended_paths) != 1:
        raise RegistryError("BL21 samples must share one reference and one intended-edit table")

    intended = next(iter(intended_paths))
    if not intended.is_file():
        raise RegistryError(f"missing intended-edit table: {intended}")
    _, edits = read_tsv(intended)
    expected = {"seq_id": INTENDED[0], "start": INTENDED[1], "end": INTENDED[2], "kind": INTENDED[3]}
    if len(edits) != 1 or any(edits[0].get(key, "").lower() != value.lower() for key, value in expected.items()):
        raise RegistryError(f"intended coordinate mismatch: expected one del at {INTENDED[0]}:{INTENDED[1]}-{INTENDED[2]}")

    truth = REPO_ROOT / "benchmark/real_world/truth/bl21_curated_truth.tsv"
    truth_fields, truth_rows = read_tsv(truth)
    if not {"truth_id", "seq_id", "start", "end", "is_structural", "evidence"}.issubset(truth_fields):
        raise RegistryError("curated BL21 truth table is missing structural/evidence fields")
    truth_by_id = {row["truth_id"]: row for row in truth_rows}
    if set(truth_by_id) != set(TRUTH_COORDS) or any(row["is_structural"] != "true" for row in truth_rows):
        raise RegistryError("curated truth sample names or structural declarations do not match the BL21 registry")
    for truth_id, expected_coordinates in TRUTH_COORDS.items():
        row = truth_by_id[truth_id]
        observed_coordinates = (row["seq_id"], row["start"], row["end"])
        if observed_coordinates != expected_coordinates:
            raise RegistryError(f"curated truth coordinate mismatch for {truth_id}: expected {expected_coordinates}, found {observed_coordinates}")

    fasta_path = next(iter(references))
    genbank_path = fasta_path.parent / "genomic.gbff"
    if not genbank_path.is_file():
        raise RegistryError(f"missing matching GenBank reference: {genbank_path}")
    fasta = sequence_digests(fasta_sequences(fasta_path))
    genbank = sequence_digests(genbank_sequences(genbank_path))
    if fasta != genbank:
        raise RegistryError("FASTA/GenBank contig mismatch: accession, sequence length, and SHA-256 must match")
    if INTENDED[0] not in fasta:
        raise RegistryError(f"intended contig {INTENDED[0]} is absent from GCF_013167015.1")
    return len(rows)


def write_contig_identity(data_root: Path, output: Path) -> None:
    _, rows = read_tsv(MANIFEST)
    if not rows:
        raise RegistryError("BL21 registry has no samples for reference identity output")
    fasta_path = Path(rows[0]["reference"])
    if fasta_path.resolve() != (data_root / "BL21/ncbi_dataset/data/GCF_013167015.1/GCF_013167015.1_ASM1316701v1_genomic.fna").resolve():
        raise RegistryError(f"reference identity output path mismatch: {fasta_path}")
    fasta = sequence_digests(fasta_sequences(fasta_path))
    genbank = sequence_digests(genbank_sequences(fasta_path.parent / "genomic.gbff"))
    if fasta != genbank:
        raise RegistryError("FASTA/GenBank contig identity changed before identity report was written")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig_accession", "length", "sha256", "fasta_genbank_match"])
        for name, (length, digest) in sorted(fasta.items()):
            writer.writerow([name, length, digest, "true"])


def main() -> int:
    parser = argparse.ArgumentParser(description="Metadata-only BL21 registry preflight; FASTQ contents are never read.")
    parser.add_argument("--data-root", type=Path, default=Path("/hpcfs/fhome/caizhh/19_BL21_edited"))
    parser.add_argument("--identity-output", type=Path)
    args = parser.parse_args()
    try:
        samples = validate_registry(args.data_root)
        if args.identity_output:
            write_contig_identity(args.data_root, args.identity_output)
    except (RegistryError, OSError) as error:
        print(f"BL21 registry validation failed: {error}", file=sys.stderr)
        return 1
    print(f"BL21_INPUTS_V2 valid: {samples} samples, six FASTQ files, two explicit spacers per sample, peer-comparator role, GCF_013167015.1 FASTA/GenBank identity verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
