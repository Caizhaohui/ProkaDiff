#!/usr/bin/env python3
"""generate_ref.py: Generates deterministic test reference genome for Cas-OFFinder parity.

Contains two contigs:
  - `chr`: 1000 bp
  - `plasmid`: 500 bp

Cases included:
  1. cas9_exact_plus (chr: 51..73)
  2. cas9_exact_minus (chr: 151..173)
  3. cas9_1mm_plus (chr: 251..273)
  4. cas9_1mm_minus (chr: 351..373)
  5. cas9_2mm (chr: 451..473)
  6. cas9_3mm (chr: 551..573)
  7. cas9_4mm (chr: 651..673)
  8. cas9_5mm (chr: 751..773)
  9. near_contig_start (chr: 1..23)
  10. near_contig_end (chr: 978..1000)
  11. plasmid_exact_plus (plasmid: 101..123)
  12. plasmid_2mm_minus (plasmid: 251..273)
"""

from pathlib import Path

SPACER = "GAGTCCGAGCAGAAGAAGAA" # 20 bp

def revcomp(s):
    tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(tbl)[::-1]

def make_background(length):
    # Repeat a safe 4-mer without GG or CC to avoid spurious random PAMs
    pattern = "ATATAGATAT"
    s = (pattern * ((length // len(pattern)) + 2))[:length]
    return list(s)

def place(seq_list, pos_0based, target_str):
    for i, c in enumerate(target_str):
        seq_list[pos_0based + i] = c

def main():
    outdir = Path(__file__).resolve().parent
    outdir.mkdir(parents=True, exist_ok=True)

    chr_seq = make_background(1000)
    plasmid_seq = make_background(500)

    # 1. near_contig_start (0..23) -> 1..23 (0 mm, Plus)
    place(chr_seq, 0, f"{SPACER}TGG")

    # 2. cas9_exact_plus (50..73) -> 51..73 (0 mm, Plus)
    place(chr_seq, 50, f"{SPACER}CGG")

    # 3. cas9_exact_minus (150..173) -> 151..173 (0 mm, Minus)
    # top strand has revcomp(target + PAM)
    place(chr_seq, 150, revcomp(f"{SPACER}TGG"))

    # 4. cas9_1mm_plus (250..273) -> 251..273 (1 mm, Plus)
    # pos 1 (5' distal) G -> A
    s_1mm = "A" + SPACER[1:]
    place(chr_seq, 250, f"{s_1mm}CGG")

    # 5. cas9_1mm_minus (350..373) -> 351..373 (1 mm, Minus)
    # pos 20 (proximal) A -> C
    s_1mm_prox = SPACER[:19] + "C"
    place(chr_seq, 350, revcomp(f"{s_1mm_prox}CGG"))

    # 6. cas9_2mm (450..473) -> 451..473 (2 mm, Plus)
    # pos 1 G->A, pos 18 G->C
    s_2mm = "A" + SPACER[1:17] + "C" + SPACER[18:]
    place(chr_seq, 450, f"{s_2mm}AGG")

    # 7. cas9_3mm (550..573) -> 551..573 (3 mm, Plus)
    # pos 1 G->A, pos 3 G->T, pos 18 G->C
    s_3mm = "A" + SPACER[1] + "T" + SPACER[3:17] + "C" + SPACER[18:]
    place(chr_seq, 550, f"{s_3mm}CGG")

    # 8. cas9_4mm (650..673) -> 651..673 (4 mm, Plus)
    # pos 1, 3, 15 (A->T), 18
    s_4mm = "A" + SPACER[1] + "T" + SPACER[3:14] + "T" + SPACER[15:17] + "C" + SPACER[18:]
    place(chr_seq, 650, f"{s_4mm}CGG")

    # 9. cas9_5mm (750..773) -> 751..773 (5 mm, Plus)
    # pos 1, 3, 13 (A->C), 15, 18
    s_5mm = "A" + SPACER[1] + "T" + SPACER[3:12] + "C" + SPACER[13:14] + "T" + SPACER[15:17] + "C" + SPACER[18:]
    place(chr_seq, 750, f"{s_5mm}CGG")

    # 10. near_contig_end (977..1000) -> 978..1000 (0 mm, Plus)
    place(chr_seq, 977, f"{SPACER}CGG")

    # Plasmid:
    # 11. plasmid_exact_plus (100..123) -> 101..123 (0 mm, Plus)
    place(plasmid_seq, 100, f"{SPACER}TGG")

    # 12. plasmid_2mm_minus (250..273) -> 251..273 (2 mm, Minus)
    place(plasmid_seq, 250, revcomp(f"{s_2mm}TGG"))

    ref_file = outdir / "reference.fa"
    with open(ref_file, "w") as f:
        f.write(f">chr\n{''.join(chr_seq)}\n")
        f.write(f">plasmid\n{''.join(plasmid_seq)}\n")

    print(f"Generated {ref_file} with chr (1000 bp) and plasmid (500 bp)")

if __name__ == "__main__":
    main()
