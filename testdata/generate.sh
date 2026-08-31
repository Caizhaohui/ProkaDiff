#!/usr/bin/env bash
# Layer-1 synthetic FASTQ: ~10 kb FASTA, 5 SNPs + 2 short indels, wgsim with
# no extra mutations (-r 0 -R 0 -X 0).
#
# MUST run on a compute node (Slurm partition qcpu_18i), never on a login node.
set -euo pipefail

host="$(hostname -s 2>/dev/null || hostname)"
if [[ "${host}" == *login* ]]; then
  echo "error: testdata/generate.sh must not run on a login node (${host}). Submit via scripts/submit_parity_synth_snp_indel.sh to qcpu_18i." >&2
  exit 1
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-${ROOT}/testdata/generated/synth_snp_indel}"
mkdir -p "${OUT}"

if ! command -v wgsim >/dev/null 2>&1; then
  echo "error: wgsim not found on PATH. Load the BactGenome conda env on a compute node." >&2
  exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "error: python3 is required to build the synthetic reference." >&2
  exit 1
fi

python3 - "${OUT}" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(42)
bases = "ACGT"
seq = [rng.choice(bases) for _ in range(10_000)]
ref = "".join(seq)

# 1-based positions on the original reference.
snp_pos = [100, 500, 1200, 3400, 7800]
flip = {"A": "G", "G": "A", "C": "T", "T": "C"}
snps = []
for p in snp_pos:
    r = ref[p - 1]
    a = flip[r]
    snps.append((p, r, a))

ins_pos, ins_seq = 2000, "AT"
del_pos, del_size = 5600, 2
del_ref = ref[del_pos - 1 : del_pos - 1 + del_size]

mut = list(ref)
# SNPs first (no length change), then 3' indel, then 5' insert.
for p, _r, a in snps:
    mut[p - 1] = a
mut[del_pos - 1 : del_pos - 1 + del_size] = []
mut[ins_pos:ins_pos] = list(ins_seq)  # insert after 1-based position ins_pos
edited = "".join(mut)

def write_fa(path: Path, name: str, s: str) -> None:
    with path.open("w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(s), 80):
            fh.write(s[i : i + 80] + "\n")

write_fa(out / "ref.fa", "synth", ref)
write_fa(out / "mutated.fa", "synth", edited)

gd_lines = [
    "#=GENOME_DIFF 1.0",
    "#=REFSEQ file:ref.fa",
]
eid = 1
for p, _r, a in snps:
    gd_lines.append(f"SNP\t{eid}\t.\tsynth\t{p}\t{a}")
    eid += 1
gd_lines.append(f"INS\t{eid}\t.\tsynth\t{ins_pos}\t{ins_seq}")
eid += 1
gd_lines.append(f"DEL\t{eid}\t.\tsynth\t{del_pos}\t{del_size}")
(out / "truth.gd").write_text("\n".join(gd_lines) + "\n")

meta = out / "mutations.txt"
meta.write_text(
    "pos\tkind\tref\talt\n"
    + "".join(f"{p}\tSNP\t{r}\t{a}\n" for p, r, a in snps)
    + f"{ins_pos}\tINS\t.\t{ins_seq}\n"
    + f"{del_pos}\tDEL\t{del_ref}\t.\n"
)
print(f"wrote {out / 'ref.fa'} ({len(ref)} bp) and truth.gd")
PY

# Coverage ~80× PE 100 bp on 10 kb: N * 200 / 10000 ≈ 80 → N = 4000.
wgsim -N 4000 -1 100 -2 100 -e 0.001 -r 0 -R 0 -X 0 -S 42 \
  "${OUT}/mutated.fa" \
  "${OUT}/edited_R1.fastq" \
  "${OUT}/edited_R2.fastq"

echo "generated ${OUT}/edited_R1.fastq ${OUT}/edited_R2.fastq"
echo "truth: ${OUT}/truth.gd"
