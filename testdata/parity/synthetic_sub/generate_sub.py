#!/usr/bin/env python3
"""generate_sub.py: Generate synthetic reference and reads containing adjacent substitutions.

Target:
  At position 201..202 (1-based), reference is 'AC'.
  All reads covering 201..202 have 'TT' (adjacent 2-bp substitution).
  Background: 30x coverage across 600 bp.
"""

from pathlib import Path
import random

def make_ref(length=600):
    random.seed(123)
    bases = ["A", "C", "G", "T"]
    seq = [random.choice(bases) for _ in range(length)]
    # Pin position 201..203 (1-based, 0-based 200..202) to ACG
    seq[200] = "A"
    seq[201] = "C"
    seq[202] = "G"
    return "".join(seq)

def make_fastq_record(name, seq):
    qual = "I" * len(seq)  # Q40
    return f"@{name}\n{seq}\n+\n{qual}\n"

def main():
    outdir = Path(__file__).resolve().parent
    outdir.mkdir(parents=True, exist_ok=True)

    ref_seq = make_ref(600)
    ref_file = outdir / "reference.fa"
    with open(ref_file, "w") as f:
        f.write(f">ref_sub\n{ref_seq}\n")

    # Mutated genome has TT instead of AC at 201..202 (0-based 200..201)
    mut_seq = list(ref_seq)
    mut_seq[200] = "T"
    mut_seq[201] = "T"
    mut_seq = "".join(mut_seq)

    read_len = 50
    reads = []
    r_idx = 0

    # 30x coverage: stride of 2 bp across 600 bp with 50 bp reads ~ 25x-30x
    for start in range(0, len(mut_seq) - read_len + 1, 2):
        r_idx += 1
        subseq = mut_seq[start : start + read_len]
        reads.append(make_fastq_record(f"read_{r_idx}", subseq))

    fq_file = outdir / "reads.fq"
    with open(fq_file, "w") as f:
        f.writelines(reads)

    print(f"Generated {ref_file} and {fq_file} ({len(reads)} reads)")

if __name__ == "__main__":
    main()
