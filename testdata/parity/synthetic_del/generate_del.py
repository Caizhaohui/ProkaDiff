#!/usr/bin/env python3
"""generate_del.py: Generate synthetic reference and reads containing a 50-bp deletion.

Target:
  Reference: 1500 bp.
  Deletion: 50 bp from position 501 to 550 (1-based, size 50).
  Reads: 50 bp reads, realistic shotgun sampling on both strands.
  Split reads spanning the breakpoint connect 500 to 551.
"""

from pathlib import Path
import random

def make_ref(length=1500):
    random.seed(42)
    bases = ["A", "C", "G", "T"]
    seq = [random.choice(bases) for _ in range(length)]
    # Ensure no repetitive homopolymers across deletion boundaries
    seq[495:505] = list("ACGTACGTAC")
    seq[545:555] = list("TGCATGCATG")
    return "".join(seq)

def revcomp(s):
    tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(tbl)[::-1]

def make_fastq_record(name, seq):
    qual = "I" * len(seq)  # Q40
    return f"@{name}\n{seq}\n+\n{qual}\n"

def main():
    outdir = Path(__file__).resolve().parent
    outdir.mkdir(parents=True, exist_ok=True)

    ref_seq = make_ref(1500)
    ref_file = outdir / "reference.fa"
    with open(ref_file, "w") as f:
        f.write(f">ref_del\n{ref_seq}\n")

    # Mutated genome deletes 50 bp: indices 500..550 (exclusive)
    del_start = 500
    del_len = 50
    mut_seq = ref_seq[:del_start] + ref_seq[del_start + del_len:]

    read_len = 50
    reads = []

    # Realistic shotgun sampling across both strands (~50x coverage)
    random.seed(12345)
    num_reads = 1500
    for r_idx in range(1, num_reads + 1):
        start = random.randint(0, len(mut_seq) - read_len)
        subseq = mut_seq[start : start + read_len]
        if random.random() < 0.5:
            subseq = revcomp(subseq)
        reads.append(make_fastq_record(f"read_{r_idx}", subseq))

    fq_file = outdir / "reads.fq"
    with open(fq_file, "w") as f:
        f.writelines(reads)

    # Expected oracle GD
    oracle_file = outdir / "oracle.gd"
    with open(oracle_file, "w") as f:
        f.write(f"#=GENOME_DIFF\t1.0\nDEL\t1\t.\tref_del\t501\t50\n")

    print(f"Generated {ref_file}, {fq_file} ({len(reads)} shotgun reads), and {oracle_file}")

if __name__ == "__main__":
    main()
