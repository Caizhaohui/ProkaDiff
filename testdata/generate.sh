#!/usr/bin/env bash
# Layer-1 synthetic FASTQ generators. MUST run on a compute node (qcpu_18i).
#
#   testdata/generate.sh [outdir]                              # synth_snp_indel
#   testdata/generate.sh synth_snp_indel [outdir]
#   testdata/generate.sh synth_parent_child [outdir]
#   testdata/generate.sh synth_cas9_near [outdir]
set -euo pipefail

host="$(hostname -s 2>/dev/null || hostname)"
if [[ "${host}" == *login* ]]; then
  echo "error: testdata/generate.sh must not run on a login node (${host}). Submit via scripts/submit_*.sh to qcpu_18i." >&2
  exit 1
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE="synth_snp_indel"
OUT=""
if [[ $# -ge 1 ]]; then
  case "$1" in
    synth_snp_indel|synth_parent_child|synth_cas9_near)
      FIXTURE="$1"
      OUT="${2:-}"
      ;;
    *)
      OUT="$1"
      ;;
  esac
fi
OUT="${OUT:-${ROOT}/testdata/generated/${FIXTURE}}"
mkdir -p "${OUT}"

if ! command -v wgsim >/dev/null 2>&1; then
  echo "error: wgsim not found on PATH. Load the BactGenome conda env on a compute node." >&2
  exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "error: python3 is required to build the synthetic reference." >&2
  exit 1
fi

python3 - "${FIXTURE}" "${OUT}" <<'PY'
import random
import sys
from pathlib import Path

fixture = sys.argv[1]
out = Path(sys.argv[2])
bases = "ACGT"
flip = {"A": "G", "G": "A", "C": "T", "T": "C"}


def write_fa(path: Path, name: str, s: str) -> None:
    with path.open("w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(s), 80):
            fh.write(s[i : i + 80] + "\n")


def apply_snps(ref, snps):
    mut = list(ref)
    for p, _r, a in snps:
        mut[p - 1] = a
    return "".join(mut)


if fixture == "synth_snp_indel":
    rng = random.Random(42)
    ref = "".join(rng.choice(bases) for _ in range(10_000))
    snp_pos = [100, 500, 1200, 3400, 7800]
    snps = [(p, ref[p - 1], flip[ref[p - 1]]) for p in snp_pos]
    ins_pos, ins_seq = 2000, "AT"
    del_pos, del_size = 5600, 2
    del_ref = ref[del_pos - 1 : del_pos - 1 + del_size]
    mut = list(ref)
    for p, _r, a in snps:
        mut[p - 1] = a
    mut[del_pos - 1 : del_pos - 1 + del_size] = []
    mut[ins_pos:ins_pos] = list(ins_seq)
    edited = "".join(mut)
    write_fa(out / "ref.fa", "synth", ref)
    write_fa(out / "mutated.fa", "synth", edited)
    gd_lines = ["#=GENOME_DIFF 1.0", "#=REFSEQ file:ref.fa"]
    eid = 1
    for p, _r, a in snps:
        gd_lines.append(f"SNP\t{eid}\t.\tsynth\t{p}\t{a}")
        eid += 1
    gd_lines.append(f"INS\t{eid}\t.\tsynth\t{ins_pos}\t{ins_seq}")
    eid += 1
    gd_lines.append(f"DEL\t{eid}\t.\tsynth\t{del_pos}\t{del_size}")
    (out / "truth.gd").write_text("\n".join(gd_lines) + "\n")
    (out / "mutations.txt").write_text(
        "pos\tkind\tref\talt\n"
        + "".join(f"{p}\tSNP\t{r}\t{a}\n" for p, r, a in snps)
        + f"{ins_pos}\tINS\t.\t{ins_seq}\n"
        + f"{del_pos}\tDEL\t{del_ref}\t.\n"
    )
    print(f"wrote {out / 'ref.fa'} ({len(ref)} bp) and truth.gd")

elif fixture == "synth_parent_child":
    rng = random.Random(43)
    ref = "".join(rng.choice(bases) for _ in range(10_000))
    intended_pos, extra_pos = 2500, 8000
    hist_pos = []
    p = 150
    while len(hist_pos) < 50:
        if p not in (intended_pos, extra_pos) and 80 < p < 9900:
            hist_pos.append(p)
        p += 90
    hist = [(p, ref[p - 1], flip[ref[p - 1]]) for p in hist_pos]
    intended = (intended_pos, ref[intended_pos - 1], flip[ref[intended_pos - 1]])
    extra = (extra_pos, ref[extra_pos - 1], flip[ref[extra_pos - 1]])
    starter = apply_snps(ref, hist)
    edited = apply_snps(ref, hist + [intended, extra])
    write_fa(out / "ref.fa", "synth", ref)
    write_fa(out / "starter.fa", "synth", starter)
    write_fa(out / "edited.fa", "synth", edited)
    (out / "intended.tsv").write_text(
        "seq_id\tstart\tend\tref\talt\tkind\n"
        f"synth\t{intended[0]}\t{intended[0]}\t{intended[1]}\t{intended[2]}\tsnp\n"
    )
    (out / "historical_snps.txt").write_text(
        "pos\tref\talt\n" + "".join(f"{p}\t{r}\t{a}\n" for p, r, a in hist)
    )
    (out / "extra_snp.txt").write_text(f"{extra[0]}\t{extra[1]}\t{extra[2]}\n")
    print(f"wrote parent_child {out} hist={len(hist)} intended={intended[0]} extra={extra[0]}")

elif fixture == "synth_cas9_near":
    rng = random.Random(44)
    ref_l = [rng.choice(bases) for _ in range(10_000)]
    spacer = "GACTGACTGACTGACTGACT"
    plant = 3000  # 1-based start of spacer
    # spacer (20) + NGG (TGG) occupies 3000..3022
    chunk = spacer + "TGG"
    ref_l[plant - 1 : plant - 1 + len(chunk)] = list(chunk)
    snp_pos = 3030
    ref = "".join(ref_l)
    r = ref[snp_pos - 1]
    a = flip[r]
    mut = list(ref)
    mut[snp_pos - 1] = a
    write_fa(out / "ref.fa", "synth", ref)
    write_fa(out / "mutated.fa", "synth", "".join(mut))
    (out / "spacer.txt").write_text(spacer + "\n")
    (out / "snp.txt").write_text(f"{snp_pos}\t{r}\t{a}\n")
    print(f"wrote cas9_near {out} spacer@{plant} snp@{snp_pos}")

else:
    raise SystemExit(f"unknown fixture {fixture}")
PY

wgsim_n=(wgsim -N 4000 -1 100 -2 100 -e 0.001 -r 0 -R 0 -X 0)

case "${FIXTURE}" in
  synth_snp_indel)
    "${wgsim_n[@]}" -S 42 \
      "${OUT}/mutated.fa" \
      "${OUT}/edited_R1.fastq" \
      "${OUT}/edited_R2.fastq"
    echo "generated ${OUT}/edited_R1.fastq ${OUT}/edited_R2.fastq"
    echo "truth: ${OUT}/truth.gd"
    ;;
  synth_parent_child)
    "${wgsim_n[@]}" -S 42 \
      "${OUT}/starter.fa" \
      "${OUT}/starter_R1.fastq" \
      "${OUT}/starter_R2.fastq"
    "${wgsim_n[@]}" -S 43 \
      "${OUT}/edited.fa" \
      "${OUT}/edited_R1.fastq" \
      "${OUT}/edited_R2.fastq"
    echo "generated starter+edited FASTQ under ${OUT}"
    ;;
  synth_cas9_near)
    "${wgsim_n[@]}" -S 42 \
      "${OUT}/ref.fa" \
      "${OUT}/starter_R1.fastq" \
      "${OUT}/starter_R2.fastq"
    "${wgsim_n[@]}" -S 43 \
      "${OUT}/mutated.fa" \
      "${OUT}/edited_R1.fastq" \
      "${OUT}/edited_R2.fastq"
    echo "generated cas9_near FASTQ under ${OUT}"
    ;;
esac
