#!/usr/bin/env python3
"""Read-only diagnosis of the candidate-junction second pass (job 2387259).

Answers, with measured counts only (no files written):
  1. Do the 946 second-pass constructs cover the 27 oracle JC flanks?
  2. For each oracle JC, how many junction-SAM reads span the breakpoint,
     on which strands, at what MAPQ, and with how many bases on each side?

Inputs:
  --gd        breseq output.gd (oracle JC; same parser as inspect_clonal_bam_jc.py)
  --jfa       work/jc/junctions.fa (946 constructs, names jc0..jcN)
  --ref       reference FASTA/GBK (REL606) to rebuild each construct's flanks
  --jsam      work/jc/aligned.sam (second-pass alignments to the constructs)

Construct model (crates/prokadiff-evidence/src/jc_seq.rs):
  flank = max read len (36 for Clonal); breakpoint at offset flank.
  side1 minus -> left flank = rc(seq[pos-flank .. pos])   (1-based pos)
  side1 plus  -> left flank = seq[pos-1 .. pos-1+flank]
  side2 likewise; right flank drops max(overlap,0) leading bases.
"""

import argparse
import re
import sys
from collections import Counter

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")
WINDOW = 20
MIN_COVER = 28          # engine JC_MIN_COVER_BASES
FLANK = 36              # Clonal read length; engine flank = max read len

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(_RC)[::-1]


def parse_cigar(cigar):
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def read_ref(path):
    """FASTA or GenBank ORIGIN -> {name: seq} (first word of header)."""
    seqs = {}
    with open(path) as fh:
        first = fh.readline()
        rest = fh.read()
    if first.startswith(">"):
        name = None
        buf = []
        for line in (first + rest).splitlines():
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append("".join(c for c in line.strip() if c.isalpha()))
        if name is not None:
            seqs[name] = "".join(buf).upper()
    else:
        name = "chr"
        m = re.search(r"LOCUS\s+(\S+)", first)
        if m:
            name = m.group(1)
        in_origin = False
        buf = []
        for line in (first + rest).splitlines():
            t = line.strip()
            if not in_origin:
                if t.upper().startswith("ORIGIN"):
                    in_origin = True
                continue
            if t == "//":
                break
            buf.append("".join(c for c in t if c.isalpha()))
        seqs[name] = "".join(buf).upper()
    return seqs


def read_jfa(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def parse_gd_jc(path):
    jcs = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("JC\t"):
                continue
            f = line.rstrip("\n").split("\t")
            kv = {}
            for t in f[10:]:
                if "=" in t:
                    k, v = t.split("=", 1)
                    kv[k] = v
            jcs.append({
                "id": f[1],
                "seq1": f[3], "pos1": int(f[4]), "str1": int(f[5]),
                "seq2": f[6], "pos2": int(f[7]), "str2": int(f[8]),
                "overlap": f[9],
                "reject": kv.get("reject", ""),
                "freq": kv.get("frequency", ""),
            })
    return jcs


def flank(seq, pos1, minus, n):
    """Reference bases of length n extending away from the junction."""
    if pos1 < 1 or pos1 > len(seq):
        return ""
    p0 = pos1 - 1
    if minus:
        s = max(0, p0 + 1 - n)
        return rc(seq[s:p0 + 1])
    return seq[p0:p0 + n].upper()


def construct_for(jc, refseqs):
    s1 = refseqs.get(jc["seq1"], "")
    s2 = refseqs.get(jc["seq2"], "")
    left = flank(s1, jc["pos1"], jc["str1"] < 0, FLANK)
    right = flank(s2, jc["pos2"], jc["str2"] < 0, FLANK)
    trim = max(0, int(jc["overlap"]))
    right = right[trim:] if trim < len(right) else ""
    return left + right


def norm(s):
    return s.upper()


def find_all(hay, needle):
    out = []
    start = 0
    while True:
        i = hay.find(needle, start)
        if i < 0:
            return out
        out.append(i)
        start = i + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gd", required=True)
    ap.add_argument("--jfa", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--jsam", required=True)
    args = ap.parse_args()

    jcs = parse_gd_jc(args.gd)
    refseqs = read_ref(args.ref)
    jfa = read_jfa(args.jfa)
    print(f"oracle_jc_n\t{len(jcs)}")
    print(f"oracle_jc_reject_n\t{sum(1 for j in jcs if j['reject'])}")
    print(f"constructs_n\t{len(jfa)}")
    print(f"flank\t{FLANK}\tmin_cover\t{MIN_COVER}\twindow\t+/-{WINDOW}")

    # Map each oracle JC to constructs whose sequence matches the oracle
    # junction construct (exact, or the construct appears inside the oracle
    # construct, or 12-mer anchors from each flank co-occur in order).
    jc_constructs = {i: set() for i in range(len(jcs))}
    for i, jc in enumerate(jcs):
        oc = construct_for(jc, refseqs)
        if not oc:
            continue
        anchor_l = oc[:12]
        anchor_r = oc[-12:]
        for name, cs in jfa.items():
            if not cs:
                continue
            if oc in cs or cs in oc:
                jc_constructs[i].add(name)
                continue
            li = cs.find(anchor_l)
            ri = cs.rfind(anchor_r)
            if li >= 0 and ri > li:
                jc_constructs[i].add(name)

    n_with = sum(1 for i in jc_constructs if jc_constructs[i])
    print(f"oracle_jc_with_matching_construct\t{n_with}/{len(jcs)}")
    for i, jc in enumerate(jcs):
        names = sorted(jc_constructs[i], key=lambda x: int(x[2:]))
        print(f"jc_construct_map\tid={jc['id']}\t{jc['seq1']}:{jc['pos1']}({jc['str1']})"
              f"->{jc['seq2']}:{jc['pos2']}({jc['str2']})\treject={jc['reject'] or '-'}"
              f"\tconstructs={','.join(names) if names else 'NONE'}")

    # Parse the junction SAM.
    per = {i: Counter() for i in range(len(jcs))}     # spanning reads
    per_side = {i: Counter() for i in range(len(jcs))}  # min-side buckets
    per_mapq = {i: Counter() for i in range(len(jcs))}
    construct_hits = Counter()
    total_aln = 0
    mapq0 = 0
    span_any_construct = 0
    span_ge28 = 0
    span_both9 = 0
    span_both14 = 0

    with open(args.jsam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            flag = int(f[1])
            if flag & 0x4:
                continue
            rname = f[2]
            if not rname.startswith("jc"):
                continue
            total_aln += 1
            mapq = int(f[4])
            if mapq == 0:
                mapq0 += 1
            pos = int(f[3])
            ops = parse_cigar(f[5])
            span_len = sum(n for n, op in ops if op in REF_CONSUMING)
            cstart = pos - 1
            cend = cstart + span_len
            left = FLANK - cstart
            right = cend - FLANK
            minus = bool(flag & 0x10)
            construct_hits[rname] += 1
            if left <= 0 or right <= 0:
                continue
            span_any_construct += 1
            q_cover = sum(n for n, op in ops if op in "MIS=X")
            if q_cover >= MIN_COVER:
                span_ge28 += 1
            mn = min(left, right)
            if mn >= 9:
                span_both9 += 1
            if mn >= 14:
                span_both14 += 1
            for i in range(len(jcs)):
                if rname in jc_constructs[i]:
                    per[i]["n"] += 1
                    per[i]["minus" if minus else "plus"] += 1
                    if q_cover >= MIN_COVER:
                        per[i]["cover28"] += 1
                    if mn >= 9:
                        per[i]["both9"] += 1
                    if mn >= 14:
                        per[i]["both14"] += 1
                    per_side[i][f"min{min(mn, 20)}"] += 1
                    per_mapq[i][f"mq{mapq}"] += 1

    print()
    print("== junction SAM global ==")
    print(f"mapped_alignments\t{total_aln}")
    print(f"mapped_mapq0\t{mapq0}")
    print(f"span_breakpoint_any\t{span_any_construct}")
    print(f"span_breakpoint_cover>={MIN_COVER}\t{span_ge28}")
    print(f"span_breakpoint_min_side>=9\t{span_both9}")
    print(f"span_breakpoint_min_side>=14\t{span_both14}")
    print("top_constructs_by_hits\t" + ",".join(
        f"{k}:{c}" for k, c in construct_hits.most_common(10)))

    print()
    print("== per oracle JC: spanning reads on its constructs ==")
    print("jc_id\treject\tpos1\tpos2\tfreq\tconstructs\tspan_n\tplus\tminus"
          "\tcover>=28\tmin>=9\tmin>=14\tmapq_hist\tmin_side_hist")
    n_any_span = 0
    n_both_strand = 0
    n_accept_like = 0
    for i, jc in enumerate(jcs):
        c = per[i]
        names = sorted(jc_constructs[i], key=lambda x: int(x[2:]))
        if c["n"]:
            n_any_span += 1
        if c["plus"] and c["minus"]:
            n_both_strand += 1
        if c["both14"] and c["plus"] and c["minus"]:
            n_accept_like += 1
        print("\t".join(str(x) for x in [
            jc["id"], jc["reject"] or "-", jc["pos1"], jc["pos2"],
            jc["freq"] or "-", len(names), c["n"], c["plus"], c["minus"],
            c["cover28"], c["both9"], c["both14"],
            ",".join(f"{k}:{v}" for k, v in sorted(per_mapq[i].items())) or "-",
            ",".join(f"{k}:{v}" for k, v in sorted(per_side[i].items(),
                     key=lambda kv: int(kv[0][3:]))) or "-",
        ]))
    print()
    print("== summary ==")
    print(f"jc_with_any_spanning_read\t{n_any_span}/{len(jcs)}")
    print(f"jc_with_spanning_both_strands\t{n_both_strand}/{len(jcs)}")
    print(f"jc_with_spanning_min14_both_strands\t{n_accept_like}/{len(jcs)}")


if __name__ == "__main__":
    main()
