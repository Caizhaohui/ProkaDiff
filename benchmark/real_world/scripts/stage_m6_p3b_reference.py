#!/usr/bin/env python3
"""Stage and verify the CP009273.1 GenBank reference for M6-P3b.

Downloads NCBI Entrez CP009273.1 (E. coli BW25113) complete genome GenBank record,
verifies bit-identical nucleotide sequence with CP009273.1.fasta:
  Length: 4,631,469 bp
  SHA256: 3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce
"""

from __future__ import annotations

import hashlib
import sys
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
DEST = REPO / "data/references/m6_identity/CP009273.1.gbk"
FASTA_PATH = REPO / "data/references/m6_identity/CP009273.1.fasta"

EXPECTED_LEN = 4631469
EXPECTED_SHA = "3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce"
NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP009273.1&rettype=gbwithparts&retmode=text"


def extract_origin_bases(gbk_text: str) -> str:
    """Extract uppercase nucleotide sequence from GenBank ORIGIN block."""
    in_origin = False
    chunks: list[str] = []
    for line in gbk_text.splitlines():
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            if line.startswith("//"):
                break
            parts = line.split()
            for token in parts[1:]:
                chunks.append("".join(ch for ch in token if ch.isalpha()).upper())
    return "".join(chunks)


def extract_fasta_bases(fasta_path: Path) -> str:
    lines = fasta_path.read_text(encoding="utf-8").splitlines()
    chunks = ["".join(ch for ch in line if ch.isalpha()).upper() for line in lines if line and not line.startswith(">")]
    return "".join(chunks)


def count_repeat_features(gbk_text: str) -> tuple[int, int]:
    """Count repeat_region and mobile_element features."""
    n_repeat = 0
    n_mob = 0
    for line in gbk_text.splitlines():
        trimmed = line.strip()
        if trimmed.startswith("repeat_region"):
            n_repeat += 1
        elif trimmed.startswith("mobile_element"):
            n_mob += 1
    return n_repeat, n_mob


def sanitize_gbk_wrapped_lines(text: str) -> str:
    """Ensure feature qualifier wrapped lines do not accidentally begin with 'locus '."""
    bad = (
        '                     /note="restriction of DNA at 5-methylcytosine residues; at\n'
        '                     locus of e14 element"'
    )
    good = (
        '                     /note="restriction of DNA at 5-methylcytosine residues;\n'
        '                     at locus of e14 element"'
    )
    return text.replace(bad, good)


def stage_cp009273_gbk(dest: Path = DEST, force: bool = False) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        text = dest.read_text(encoding="utf-8", errors="replace")
        text = sanitize_gbk_wrapped_lines(text)
        dest.write_text(text, encoding="utf-8")
        bases = extract_origin_bases(text)
        digest = hashlib.sha256(bases.encode("ascii")).hexdigest()
        if len(bases) == EXPECTED_LEN and digest == EXPECTED_SHA:
            print(f"Existing {dest} is valid (length={len(bases)}, sha256={digest}).")
            return dest
        print(f"Existing {dest} failed validation; redownloading...")

    print(f"Downloading CP009273.1 from NCBI Entrez...")
    req = urllib.request.Request(NCBI_URL, headers={"User-Agent": "ProkaDiff-M6P3b/1.0"})
    max_retries = 3
    for attempt in range(1, max_retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                content = resp.read()
            text = content.decode("utf-8", errors="replace")
            text = sanitize_gbk_wrapped_lines(text)
            bases = extract_origin_bases(text)
            digest = hashlib.sha256(bases.encode("ascii")).hexdigest()
            if len(bases) != EXPECTED_LEN:
                raise ValueError(f"Length mismatch: got {len(bases)}, expected {EXPECTED_LEN}")
            if digest != EXPECTED_SHA:
                raise ValueError(f"SHA256 mismatch: got {digest}, expected {EXPECTED_SHA}")
            dest.write_text(text, encoding="utf-8")
            print(f"Successfully staged {dest} ({len(text)} bytes, {len(bases)} bp).")
            return dest
        except Exception as err:
            print(f"Attempt {attempt}/{max_retries} failed: {err}")
            if attempt == max_retries:
                raise
            time.sleep(5)
    raise RuntimeError("Failed to stage CP009273.1.gbk")


def main() -> int:
    dest = stage_cp009273_gbk()
    text = dest.read_text(encoding="utf-8")
    bases = extract_origin_bases(text)
    digest = hashlib.sha256(bases.encode("ascii")).hexdigest()
    fasta_bases = extract_fasta_bases(FASTA_PATH)
    if bases != fasta_bases:
        print("ERROR: GenBank ORIGIN does not exactly match CP009273.1.fasta!", file=sys.stderr)
        return 1
    n_rep, n_mob = count_repeat_features(text)
    print(f"Verified CP009273.1.gbk:")
    print(f"  Length: {len(bases)} bp")
    print(f"  SHA256: {digest}")
    print(f"  Match with FASTA: 100% exact")
    print(f"  repeat_region features: {n_rep}")
    print(f"  mobile_element features: {n_mob}")
    print(f"  total repeat annotations: {n_rep + n_mob}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
