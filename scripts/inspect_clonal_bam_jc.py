#!/usr/bin/env python3
"""Read-only Clonal BAM inspection: S>=14 clips near oracle JC.

Reads SAM text on stdin (samtools view aligned.bam). Parses 27 oracle JC
from a Genome Diff. Reports measured counts only. Does not write files.

Window: aligned block start or end within +/-20 bp of either JC side.
Engine-adjacent: leading S>=14 with start in window, or trailing S>=14
with end in window (what apply_read actually records).
"""

import argparse
import re
import sys
from collections import Counter

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")
MIN_CLIP = 14
UNIQUE_MAPQ = 10
WINDOW = 20
EXAMPLE_JC = 5
EXAMPLE_READS = 3


def parse_cigar(cigar):
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def lead_trail_s(ops):
    lead = 0
    trail = 0
    if ops:
        if ops[0][1] == "S":
            lead = ops[0][0]
        elif len(ops) > 1 and ops[0][1] == "H" and ops[1][1] == "S":
            lead = ops[1][0]
        if ops[-1][1] == "S":
            trail = ops[-1][0]
        elif len(ops) > 1 and ops[-1][1] == "H" and ops[-2][1] == "S":
            trail = ops[-2][0]
    return lead, trail


def parse_gd_jc(path):
    jcs = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("JC\t"):
                continue
            f = line.rstrip("\n").split("\t")
            kv = {}
            extra_start = 10
            for t in f[extra_start:]:
                if "=" in t:
                    k, v = t.split("=", 1)
                    kv[k] = v
            reject = kv.get("reject", "")
            jcs.append({
                "id": f[1],
                "seq1": f[3],
                "pos1": int(f[4]),
                "str1": int(f[5]),
                "seq2": f[6],
                "pos2": int(f[7]),
                "str2": int(f[8]),
                "overlap": f[9],
                "reject": reject,
                "freq": kv.get("frequency", ""),
                "nj_reads": kv.get("new_junction_read_count", ""),
            })
    return jcs


def windows_for(pos, chrom_len):
    """Inclusive 1-based positions within +/- WINDOW, with circular wrap."""
    lo = pos - WINDOW
    hi = pos + WINDOW
    ranges = []
    if chrom_len and lo < 1:
        ranges.append((chrom_len + lo, chrom_len))
        ranges.append((1, hi))
    elif chrom_len and hi > chrom_len:
        ranges.append((lo, chrom_len))
        ranges.append((1, hi - chrom_len))
    else:
        ranges.append((max(1, lo), hi if not chrom_len else min(chrom_len, hi)))
    return ranges


def in_windows(pos, ranges):
    return any(a <= pos <= b for a, b in ranges)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gd", required=True, help="breseq output.gd (oracle JC)")
    ap.add_argument("--label", default="stdin")
    args = ap.parse_args()

    jcs = parse_gd_jc(args.gd)
    print(f"== {args.label} ==")
    print(f"gd\t{args.gd}")
    print(f"oracle_jc_n\t{len(jcs)}")
    print(f"oracle_jc_reject_n\t{sum(1 for j in jcs if j['reject'])}")
    print(f"window_bp\t{WINDOW}")
    print(f"min_clip\t{MIN_CLIP}")
    print(f"unique_mapq\t{UNIQUE_MAPQ}")

    chrom_len = {}
    # per JC: side -> counters
    # keys: n, plus, minus, mapq_ge10, mapq_lt10,
    #       eng, eng_plus, eng_minus, eng_mapq_ge10, eng_mapq_lt10
    def fresh():
        return Counter()

    side_stats = {(i, side): fresh() for i in range(len(jcs)) for side in (1, 2)}
    examples = {i: [] for i in range(min(EXAMPLE_JC, len(jcs)))}

    total_mapped = 0
    any_s = 0
    s_ge14 = 0
    s_ge14_mapq_ge10 = 0
    s_ge14_mapq_lt10 = 0
    h_ge14 = 0
    split_s_sum_ge14_max_lt14 = 0
    clip_len_ge14_among_s14_mq10 = 0
    max_clip_hist = Counter()
    any_s_hist = Counter()
    cigar_pat = Counter()
    rnames = Counter()
    plus_mapped = 0
    minus_mapped = 0
    n_supp = 0
    n_sa_tag = 0
    global_s_examples = []
    any_s_near = {i: [] for i in range(min(EXAMPLE_JC, len(jcs)))}
    any_read_near = {i: [] for i in range(min(EXAMPLE_JC, len(jcs)))}
    jc_any_s = {(i, side): Counter() for i in range(len(jcs)) for side in (1, 2)}
    jc_max_s = {(i, side): 0 for i in range(len(jcs)) for side in (1, 2)}

    # Precompute windows once chrom_len known; start with linear, wrap later
    # if header provides LN. Header lines come first.
    win = None  # filled after header or lazily

    def build_win():
        out = []
        for j in jcs:
            s1 = windows_for(j["pos1"], chrom_len.get(j["seq1"]))
            s2 = windows_for(j["pos2"], chrom_len.get(j["seq2"]))
            out.append((s1, s2, j["seq1"], j["seq2"]))
        return out

    header_done = False
    for line in sys.stdin:
        if line.startswith("@"):
            if line.startswith("@SQ"):
                tags = dict(
                    t.split(":", 1) for t in line[3:].split() if ":" in t
                )
                if "SN" in tags and "LN" in tags:
                    chrom_len[tags["SN"]] = int(tags["LN"])
            continue
        if not header_done:
            win = build_win()
            header_done = True
            print("bam_sq\t" + ",".join(
                f"{k}:{v}" for k, v in sorted(chrom_len.items())))

        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        qname, flag_s, rname, pos_s, mapq_s, cigar = (
            f[0], f[1], f[2], f[3], f[4], f[5]
        )
        flag = int(flag_s)
        if flag & 0x4:
            continue
        pos = int(pos_s)
        mapq = int(mapq_s)
        ops = parse_cigar(cigar)
        ref_len = sum(n for n, op in ops if op in REF_CONSUMING)
        ref_end = pos + ref_len - 1 if ref_len else pos
        lead_s, trail_s = lead_trail_s(ops)
        s_ops = [n for n, op in ops if op == "S"]
        h_ops = [n for n, op in ops if op == "H"]
        max_s = max(s_ops) if s_ops else 0
        max_h = max(h_ops) if h_ops else 0
        sum_s = sum(s_ops)
        minus = bool(flag & 0x10)
        total_mapped += 1
        if minus:
            minus_mapped += 1
        else:
            plus_mapped += 1
        rnames[rname] += 1
        if flag & 0x800:
            n_supp += 1
        if any(t.startswith("SA:Z:") for t in f[11:]):
            n_sa_tag += 1

        if s_ops:
            any_s += 1
            bucket = (
                "1-5" if max_s <= 5
                else "6-13" if max_s <= 13
                else "14-19" if max_s <= 19
                else "20-39" if max_s <= 39
                else "40+"
            )
            any_s_hist[bucket] += 1
            if len(global_s_examples) < 5:
                global_s_examples.append((qname, flag, pos, mapq, cigar, max_s))
        if max_h >= MIN_CLIP:
            h_ge14 += 1
        if max_s >= MIN_CLIP:
            s_ge14 += 1
            if mapq >= UNIQUE_MAPQ:
                s_ge14_mapq_ge10 += 1
                clip_len_ge14_among_s14_mq10 += 1
            else:
                s_ge14_mapq_lt10 += 1
            bucket = (
                "14-19" if max_s <= 19
                else "20-39" if max_s <= 39
                else "40-59" if max_s <= 59
                else "60-79" if max_s <= 79
                else "80+"
            )
            max_clip_hist[bucket] += 1
            cigar_pat[re.sub(r"\d+", "N", cigar)] += 1
        elif sum_s >= MIN_CLIP:
            split_s_sum_ge14_max_lt14 += 1

        strand = "minus" if minus else "plus"
        mq_bin = "ge10" if mapq >= UNIQUE_MAPQ else "lt10"
        for i, (w1, w2, seq1, seq2) in enumerate(win):
            name_ok1 = rname == seq1 or seq1 in rname or rname in seq1
            name_ok2 = rname == seq2 or seq2 in rname or rname in seq2
            near1 = name_ok1 and (in_windows(pos, w1) or in_windows(ref_end, w1))
            near2 = name_ok2 and (in_windows(pos, w2) or in_windows(ref_end, w2))
            if near1 or near2:
                if i < EXAMPLE_JC and len(any_read_near[i]) < EXAMPLE_READS:
                    any_read_near[i].append(
                        (qname, flag, pos, mapq, cigar, ref_end, strand,
                         max_s, lead_s, trail_s)
                    )
                if max_s > 0 and i < EXAMPLE_JC and len(any_s_near[i]) < EXAMPLE_READS:
                    any_s_near[i].append(
                        (qname, flag, pos, mapq, cigar, ref_end, strand,
                         max_s, lead_s, trail_s)
                    )
                for side, near in ((1, near1), (2, near2)):
                    if not near:
                        continue
                    jc_any_s[(i, side)]["n"] += 1
                    if max_s > 0:
                        jc_any_s[(i, side)]["any_s"] += 1
                        jc_max_s[(i, side)] = max(jc_max_s[(i, side)], max_s)

        if max_s < MIN_CLIP:
            continue

        for i, (w1, w2, seq1, seq2) in enumerate(win):
            near1 = (rname == seq1 or seq1 in rname or rname in seq1) and (
                in_windows(pos, w1) or in_windows(ref_end, w1)
            )
            near2 = (rname == seq2 or seq2 in rname or rname in seq2) and (
                in_windows(pos, w2) or in_windows(ref_end, w2)
            )
            if not (near1 or near2):
                continue
            eng1 = (rname == seq1 or seq1 in rname or rname in seq1) and (
                (lead_s >= MIN_CLIP and in_windows(pos, w1))
                or (trail_s >= MIN_CLIP and in_windows(ref_end, w1))
            )
            eng2 = (rname == seq2 or seq2 in rname or rname in seq2) and (
                (lead_s >= MIN_CLIP and in_windows(pos, w2))
                or (trail_s >= MIN_CLIP and in_windows(ref_end, w2))
            )
            for side, near, eng in ((1, near1, eng1), (2, near2, eng2)):
                if not near:
                    continue
                st = side_stats[(i, side)]
                st["n"] += 1
                st[strand] += 1
                st[f"mapq_{mq_bin}"] += 1
                if eng:
                    st["eng"] += 1
                    st[f"eng_{strand}"] += 1
                    st[f"eng_mapq_{mq_bin}"] += 1
            if i < EXAMPLE_JC and len(examples[i]) < EXAMPLE_READS:
                examples[i].append(
                    (qname, flag, pos, mapq, cigar, ref_end, strand,
                     "s1" if near1 else "", "s2" if near2 else "",
                     lead_s, trail_s)
                )

    if not header_done:
        print("bam_sq\t(none; empty SAM)")

    print()
    print("== global ==")
    print(f"total_mapped\t{total_mapped}")
    print(f"mapped_plus\t{plus_mapped}")
    print(f"mapped_minus\t{minus_mapped}")
    print(f"reads_with_any_S\t{any_s}")
    print("any_S_max_len_bucket\t" + ",".join(
        f"{k}:{any_s_hist[k]}" for k in ["1-5", "6-13", "14-19", "20-39", "40+"]))
    print(f"reads_S_ge{MIN_CLIP}\t{s_ge14}")
    print(f"reads_S_ge{MIN_CLIP}_mapq_ge{UNIQUE_MAPQ}\t{s_ge14_mapq_ge10}")
    print(f"reads_S_ge{MIN_CLIP}_mapq_lt{UNIQUE_MAPQ}\t{s_ge14_mapq_lt10}")
    print(f"of_S_ge{MIN_CLIP}_mapq_ge{UNIQUE_MAPQ}_clip_len_ge{MIN_CLIP}\t"
          f"{clip_len_ge14_among_s14_mq10}")
    print(f"reads_H_ge{MIN_CLIP}\t{h_ge14}")
    print(f"reads_sumS_ge{MIN_CLIP}_maxS_lt{MIN_CLIP}\t{split_s_sum_ge14_max_lt14}")
    print(f"supplementary_flag_0x800\t{n_supp}")
    print(f"reads_with_SA_tag\t{n_sa_tag}")
    print("max_clip_len_bucket_S_ge14\t" + ",".join(
        f"{k}:{max_clip_hist[k]}" for k in
        ["14-19", "20-39", "40-59", "60-79", "80+"]))
    print("top_rnames\t" + ",".join(
        f"{k}:{c}" for k, c in rnames.most_common(5)))
    print("top_cigar_S_ge14\t" + ",".join(
        f"{p}:{c}" for p, c in cigar_pat.most_common(12)))
    print("global_any_S_examples (QNAME FLAG POS MAPQ CIGAR maxS):")
    for q, fl, p, m, c, ms in global_s_examples:
        print(f"  {q}\t{fl}\t{p}\t{m}\t{c}\tmaxS={ms}")
    if not global_s_examples:
        print("  (none)")
    print()

    n_any = 0
    n_mq10 = 0
    n_eng = 0
    n_eng_mq10 = 0
    n_both_sides_eng_mq10 = 0
    n_reject_any = 0
    n_keep_any = 0

    print("== per_jc (user window: S>=14 and start/end +/-20 of side) ==")
    print("jc_id\treject\tseq1\tpos1\tstr1\tseq2\tpos2\tstr2\toverlap\tnj_reads"
          "\ts1_n\ts1_plus\ts1_minus\ts1_mq>=10\ts1_mq<10"
          "\ts1_eng\ts1_eng_plus\ts1_eng_minus\ts1_eng_mq>=10\ts1_eng_mq<10"
          "\ts2_n\ts2_plus\ts2_minus\ts2_mq>=10\ts2_mq<10"
          "\ts2_eng\ts2_eng_plus\ts2_eng_minus\ts2_eng_mq>=10\ts2_eng_mq<10")
    for i, j in enumerate(jcs):
        s1 = side_stats[(i, 1)]
        s2 = side_stats[(i, 2)]
        row = [
            j["id"], j["reject"] or "-", j["seq1"], j["pos1"], j["str1"],
            j["seq2"], j["pos2"], j["str2"], j["overlap"], j["nj_reads"] or "-",
            s1["n"], s1["plus"], s1["minus"], s1["mapq_ge10"], s1["mapq_lt10"],
            s1["eng"], s1["eng_plus"], s1["eng_minus"], s1["eng_mapq_ge10"],
            s1["eng_mapq_lt10"],
            s2["n"], s2["plus"], s2["minus"], s2["mapq_ge10"], s2["mapq_lt10"],
            s2["eng"], s2["eng_plus"], s2["eng_minus"], s2["eng_mapq_ge10"],
            s2["eng_mapq_lt10"],
        ]
        print("\t".join(str(x) for x in row))
        any_n = s1["n"] + s2["n"]
        any_mq = s1["mapq_ge10"] + s2["mapq_ge10"]
        any_eng = s1["eng"] + s2["eng"]
        any_eng_mq = s1["eng_mapq_ge10"] + s2["eng_mapq_ge10"]
        if any_n:
            n_any += 1
            if j["reject"]:
                n_reject_any += 1
            else:
                n_keep_any += 1
        if any_mq:
            n_mq10 += 1
        if any_eng:
            n_eng += 1
        if any_eng_mq:
            n_eng_mq10 += 1
        if s1["eng_mapq_ge10"] and s2["eng_mapq_ge10"]:
            n_both_sides_eng_mq10 += 1

    print()
    print("== summary_across_27 ==")
    print(f"jc_with_any_Sge14_in_window\t{n_any}/{len(jcs)}")
    print(f"  of_which_reject\t{n_reject_any}")
    print(f"  of_which_kept\t{n_keep_any}")
    print(f"jc_with_Sge14_MAPQge10_in_window\t{n_mq10}/{len(jcs)}")
    print(f"jc_with_engine_adjacent_Sge14\t{n_eng}/{len(jcs)}")
    print(f"jc_with_engine_adjacent_Sge14_MAPQge10\t{n_eng_mq10}/{len(jcs)}")
    print(f"jc_both_sides_engine_adjacent_MAPQge10\t{n_both_sides_eng_mq10}/{len(jcs)}")
    print()
    print("== per_jc any-S (not ge14) in +/-20 start/end window ==")
    print("jc_id\treject\tpos1\tpos2\ts1_reads\ts1_anyS\ts1_maxS\ts2_reads\ts2_anyS\ts2_maxS")
    n_win_cov = 0
    n_win_anys = 0
    n_win_maxs_ge14 = 0
    for i, j in enumerate(jcs):
        a1 = jc_any_s[(i, 1)]
        a2 = jc_any_s[(i, 2)]
        m1 = jc_max_s[(i, 1)]
        m2 = jc_max_s[(i, 2)]
        print("\t".join(str(x) for x in [
            j["id"], j["reject"] or "-", j["pos1"], j["pos2"],
            a1["n"], a1["any_s"], m1, a2["n"], a2["any_s"], m2,
        ]))
        if a1["n"] + a2["n"]:
            n_win_cov += 1
        if a1["any_s"] + a2["any_s"]:
            n_win_anys += 1
        if m1 >= MIN_CLIP or m2 >= MIN_CLIP:
            n_win_maxs_ge14 += 1
    print(f"jc_with_any_mapped_in_window\t{n_win_cov}/{len(jcs)}")
    print(f"jc_with_any_S_in_window\t{n_win_anys}/{len(jcs)}")
    print(f"jc_with_maxS_ge14_in_window\t{n_win_maxs_ge14}/{len(jcs)}")
    print()
    print(f"== examples first {EXAMPLE_JC} JC x {EXAMPLE_READS} "
          "S>=14 (none expected if global S_ge14=0) ==")
    for i in range(min(EXAMPLE_JC, len(jcs))):
        j = jcs[i]
        print(f"JC id={j['id']} {j['seq1']}:{j['pos1']}({j['str1']}) -> "
              f"{j['seq2']}:{j['pos2']}({j['str2']}) reject={j['reject'] or '-'}")
        if not examples[i]:
            print("  (none)")
            continue
        for ex in examples[i]:
            q, fl, p, m, c, rend, strand, n1, n2, ls, ts = ex
            print(f"  {q}\t{fl}\t{p}\t{m}\t{c}\tref_end={rend}\t{strand}\t"
                  f"{n1}{n2}\tleadS={ls}\ttrailS={ts}")
    print()
    print(f"== examples first {EXAMPLE_JC} JC: any-S then any-mapped "
          "(QNAME FLAG POS MAPQ CIGAR) ==")
    for i in range(min(EXAMPLE_JC, len(jcs))):
        j = jcs[i]
        print(f"JC id={j['id']} {j['seq1']}:{j['pos1']}({j['str1']}) -> "
              f"{j['seq2']}:{j['pos2']}({j['str2']}) reject={j['reject'] or '-'}")
        print("  any_S_near:")
        if not any_s_near[i]:
            print("    (none)")
        for q, fl, p, m, c, rend, strand, ms, ls, ts in any_s_near[i]:
            print(f"    {q}\t{fl}\t{p}\t{m}\t{c}\tmaxS={ms}\tleadS={ls}\t"
                  f"trailS={ts}\tref_end={rend}\t{strand}")
        print("  any_mapped_near:")
        if not any_read_near[i]:
            print("    (none)")
        for q, fl, p, m, c, rend, strand, ms, ls, ts in any_read_near[i]:
            print(f"    {q}\t{fl}\t{p}\t{m}\t{c}\tmaxS={ms}\tleadS={ls}\t"
                  f"trailS={ts}\tref_end={rend}\t{strand}")


if __name__ == "__main__":
    main()
