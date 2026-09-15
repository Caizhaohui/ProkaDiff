#!/usr/bin/env python3
"""generate_ambiguous.py: Generate synthetic reference and reads to test ambiguous N handling.

Cases:
  Case 1 (pos 100, Ref=A): 10x A reads, 5x N reads.
  Case 2 (pos 200, Ref=G): 10x G reads, 20x N reads.
  Case 3 (pos 300, Ref=C): Reads with >50% N across read length.
  Case 4 (pos 400, Ref=T): Reads with low base quality and N-rich tiles.
  Background: 25x normal coverage across a 600 bp synthetic reference.
"""

from pathlib import Path
import random

def make_ref(length=600):
    random.seed(42)
    bases = ["A", "C", "G", "T"]
    seq = [random.choice(bases) for _ in range(length)]
    # Pin specific positions (1-based: 100, 200, 300, 400)
    seq[99] = "A"
    seq[199] = "G"
    seq[299] = "C"
    seq[399] = "T"
    return "".join(seq)

def make_fastq_record(name, seq, qual=None):
    if qual is None:
        qual = "I" * len(seq)  # Q40
    return f"@{name}\n{seq}\n+\n{qual}\n"

def main():
    outdir = Path(__file__).resolve().parent
    outdir.mkdir(parents=True, exist_ok=True)

    ref_seq = make_ref(600)
    ref_file = outdir / "reference.fa"
    with open(ref_file, "w") as f:
        f.write(f">ref_ambiguous\n{ref_seq}\n")

    read_len = 50
    reads = []
    r_idx = 0

    # 1. Background coverage: ~20x spanning across 1..600
    for start in range(0, len(ref_seq) - read_len + 1, 2):
        r_idx += 1
        subseq = ref_seq[start : start + read_len]
        reads.append(make_fastq_record(f"bg_{r_idx}", subseq))

    # Case 1: pos 100 (0-based 99, Ref=A). 10x A, 5x N spanning pos 100
    # Read starts at 75..125, pos 99 is at read offset 24
    for i in range(10):
        r_idx += 1
        subseq = list(ref_seq[75 : 75 + read_len])
        subseq[24] = "A"
        reads.append(make_fastq_record(f"case1_A_{i}", "".join(subseq)))
    for i in range(5):
        r_idx += 1
        subseq = list(ref_seq[75 : 75 + read_len])
        subseq[24] = "N"
        reads.append(make_fastq_record(f"case1_N_{i}", "".join(subseq)))

    # Case 2: pos 200 (0-based 199, Ref=G). 10x G, 20x N spanning pos 200
    # Read starts at 175..225, pos 199 is at read offset 24
    for i in range(10):
        r_idx += 1
        subseq = list(ref_seq[175 : 175 + read_len])
        subseq[199 - 175] = "G"
        reads.append(make_fastq_record(f"case2_G_{i}", "".join(subseq)))
    for i in range(20):
        r_idx += 1
        subseq = list(ref_seq[175 : 175 + read_len])
        subseq[199 - 175] = "N"
        reads.append(make_fastq_record(f"case2_N_{i}", "".join(subseq)))

    # Case 3: pos 300 (0-based 299). Reads with >50% N (e.g. 30 Ns out of 50 bp)
    for i in range(8):
        r_idx += 1
        subseq = list(ref_seq[275 : 275 + read_len])
        for n_pos in range(10, 40):
            subseq[n_pos] = "N"
        reads.append(make_fastq_record(f"case3_highN_{i}", "".join(subseq)))

    # Case 4: pos 400 (0-based 399). Low quality N-rich reads (qual = '!' or '#', Q0-Q2)
    for i in range(8):
        r_idx += 1
        subseq = list(ref_seq[375 : 375 + read_len])
        for n_pos in range(20, 30):
            subseq[n_pos] = "N"
        quals = ["I"] * read_len
        for n_pos in range(20, 30):
            quals[n_pos] = "!"  # Q0
        reads.append(make_fastq_record(f"case4_lowqualN_{i}", "".join(subseq), "".join(quals)))

    fq_file = outdir / "reads.fq"
    with open(fq_file, "w") as f:
        f.writelines(reads)

    print(f"Generated {ref_file} and {fq_file} ({len(reads)} reads)")

if __name__ == "__main__":
    main()
