#!/usr/bin/env python3
"""Read-only BAM/SAM inspection for the synth_is_mob zero-JC diagnosis.

Reads SAM text on stdin (samtools view aligned.bam, or cat work/aligned.sam)
and reports measured counts only. No files are written.

Geometry below matches the PRE-FIX fixture used by job 2371383 (1440 bp,
insert at 201). testdata/generate.sh has since moved the insert to 601 in a
~4000 bp genome; this script is kept as the historical diagnosis of 2371383.

Geometry (testdata/generate.sh synth_is_mob, pre-fix):
  ref = left(1-400) + motif(401-480) + mid(481-880) + motif(881-960)
        + right(961-1360)
  edited genome inserts a third motif copy between ref positions 200 and 201.
"""

import re
import sys
from collections import Counter

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")

MIN_CLIP = 14  # engine MIN_CLIP_FOR_JC
UNIQUE_MAPQ = 10  # engine UNIQUE_MAPQ

J_END = 200       # aligned block ending here = left side of insert junction
J_START = 201     # aligned block starting here = right side of insert junction
NEAR = 5
MOTIF_COPIES = [(401, 480), (881, 960)]
COPY_EDGE_ENDS = {400, 480, 481, 880, 960, 961}
COPY_EDGE_STARTS = {401, 481, 881, 960, 961}
INDEL_WINDOW = (150, 550)


def parse_cigar(cigar):
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def main():
    label = "stdin"
    if "--label" in sys.argv:
        label = sys.argv[sys.argv.index("--label") + 1]

    total = 0
    minus_total = 0
    with_s = 0
    minus_with_s = 0
    max_clip_hist = Counter()
    clip_ge14 = 0
    clip_ge14_minus = 0
    clip_ge14_mapq_lt10 = 0
    flag_ge14 = Counter()

    clip_pos = Counter()          # (side, pos), strand-unaware (v1 compat)
    clip_pos_strand = Counter()   # (strand, side, pos)
    end200 = []
    start201 = []
    near_end = Counter()
    near_start = Counter()
    j_mapq = Counter()
    j_mapq_lt10 = 0
    j_strand_ge10 = Counter()

    long_indel = 0
    long_indel_examples = []

    copy_edge = Counter()
    copy_edge_strand = Counter()
    motif_intersect_s = 0
    motif_intersect_s_lt10 = 0
    motif_full_nm = Counter()
    motif_full_lt10 = 0
    motif_full_total = 0

    cigar_pattern = Counter()
    examples_pool = []
    minus_examples = []

    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        qname, flag_s, _rname, pos_s, mapq_s, cigar = f[0], f[1], f[2], f[3], f[4], f[5]
        flag = int(flag_s)
        pos = int(pos_s)
        mapq = int(mapq_s)
        if flag & 0x4:
            continue
        total += 1
        minus = bool(flag & 0x10)
        if minus:
            minus_total += 1
        ops = parse_cigar(cigar)
        ref_len = sum(n for n, op in ops if op in REF_CONSUMING)
        ref_end = pos + ref_len - 1
        s_lens = [n for n, op in ops if op == "S"]
        nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), None)
        strand = "minus" if minus else "plus"

        intersects_motif = any(pos <= e and ref_end >= s for s, e in MOTIF_COPIES)

        if not s_lens:
            if intersects_motif:
                motif_full_total += 1
                if mapq < UNIQUE_MAPQ:
                    motif_full_lt10 += 1
                if nm is not None:
                    motif_full_nm[0 if nm == 0 else (1 if nm <= 2 else 3)] += 1
        else:
            with_s += 1
            if minus:
                minus_with_s += 1
            mx = max(s_lens)
            bucket = (
                "0" if mx == 0
                else "1-5" if mx <= 5
                else "6-13" if mx <= 13
                else "14-19" if mx <= 19
                else "20-39" if mx <= 39
                else "40-59" if mx <= 59
                else "60-79" if mx <= 79
                else "80+"
            )
            max_clip_hist[bucket] += 1
            if intersects_motif:
                motif_intersect_s += 1
                if mapq < UNIQUE_MAPQ:
                    motif_intersect_s_lt10 += 1

            if mx >= MIN_CLIP:
                clip_ge14 += 1
                flag_ge14[flag] += 1
                if minus:
                    clip_ge14_minus += 1
                    if len(minus_examples) < 5:
                        minus_examples.append((qname, flag, pos, mapq, cigar))
                if mapq < UNIQUE_MAPQ:
                    clip_ge14_mapq_lt10 += 1
                cigar_pattern[re.sub(r"\d+", "N", cigar)] += 1
                lead_s = ops[0][1] == "S"
                trail_s = ops[-1][1] == "S"
                if trail_s:
                    clip_pos[("end", ref_end)] += 1
                    clip_pos_strand[(strand, "end", ref_end)] += 1
                if lead_s:
                    clip_pos[("start", pos)] += 1
                    clip_pos_strand[(strand, "start", pos)] += 1

                is_j = (trail_s and ref_end == J_END) or (lead_s and pos == J_START)
                if trail_s and ref_end == J_END:
                    end200.append((qname, flag, pos, mapq, cigar, ops[-1][0]))
                if lead_s and pos == J_START:
                    start201.append((qname, flag, pos, mapq, cigar, ops[0][0]))
                if trail_s and ref_end != J_END and abs(ref_end - J_END) <= NEAR:
                    near_end[ref_end] += 1
                if lead_s and pos != J_START and abs(pos - J_START) <= NEAR:
                    near_start[pos] += 1
                if is_j:
                    j_mapq[mapq] += 1
                    if mapq < UNIQUE_MAPQ:
                        j_mapq_lt10 += 1
                    else:
                        j_strand_ge10[strand] += 1
                    examples_pool.append((qname, flag, pos, mapq, cigar))

                if trail_s and ref_end in COPY_EDGE_ENDS:
                    copy_edge[("end", ref_end)] += 1
                    copy_edge_strand[(strand, "end", ref_end)] += 1
                if lead_s and pos in COPY_EDGE_STARTS:
                    copy_edge[("start", pos)] += 1
                    copy_edge_strand[(strand, "start", pos)] += 1

        if pos <= INDEL_WINDOW[1] and ref_end >= INDEL_WINDOW[0]:
            long_ops = [(n, op) for n, op in ops if op in "ID" and n > 2]
            if long_ops:
                long_indel += 1
                if len(long_indel_examples) < 5:
                    long_indel_examples.append((qname, flag, pos, mapq, cigar))

    print(f"== {label} ==")
    print(f"total_mapped_reads\t{total}")
    print(f"total_minus_strand\t{minus_total}")
    print(f"reads_with_S_any\t{with_s}")
    print(f"reads_with_S_any_minus\t{minus_with_s}")
    print("max_clip_len_bucket\t" + ",".join(f"{k}:{max_clip_hist[k]}" for k in
          ["0", "1-5", "6-13", "14-19", "20-39", "40-59", "60-79", "80+"]))
    print(f"reads_with_S_ge{MIN_CLIP}\t{clip_ge14}")
    print(f"reads_with_S_ge{MIN_CLIP}_minus\t{clip_ge14_minus}")
    print(f"reads_with_S_ge{MIN_CLIP}_mapq_lt{UNIQUE_MAPQ}\t{clip_ge14_mapq_lt10}")
    print("flag_counts_S_ge14\t" + ",".join(f"{fl}:{c}" for fl, c in flag_ge14.most_common(10)))
    print()
    print(f"junction_end_at_{J_END}\t{len(end200)}")
    print(f"junction_start_at_{J_START}\t{len(start201)}")
    print("near_miss_end_196-204_excl200\t" + ",".join(f"{p}:{c}" for p, c in sorted(near_end.items())))
    print("near_miss_start_196-206_excl201\t" + ",".join(f"{p}:{c}" for p, c in sorted(near_start.items())))
    print(f"junction_mapq_lt{UNIQUE_MAPQ}\t{j_mapq_lt10}")
    print("junction_mapq_values\t" + ",".join(f"{m}:{c}" for m, c in sorted(j_mapq.items())))
    print(f"junction_strand_mapq_ge{UNIQUE_MAPQ}\t" + ",".join(f"{s}:{c}" for s, c in sorted(j_strand_ge10.items())))
    print()
    print("top_clip_adjacent_positions\t" + ",".join(f"{k[0]}@{k[1]}:{c}" for k, c in clip_pos.most_common(20)))
    print("top_clip_adjacent_by_strand\t" + ",".join(
        f"{k[0]}|{k[1]}@{k[2]}:{c}" for k, c in clip_pos_strand.most_common(30)))
    print("top_cigar_patterns_S_ge14\t" + ",".join(f"{p}:{c}" for p, c in cigar_pattern.most_common(10)))
    print()
    print(f"long_ID_gt2_near_{INDEL_WINDOW[0]}_{INDEL_WINDOW[1]}\t{long_indel}")
    for q, fl, p, m, c in long_indel_examples:
        print(f"  indel_example\t{q}\tflag={fl}\tpos={p}\tmapq={m}\tcigar={c}")
    print()
    print("motif_copy_edge_softclips\t" + ",".join(f"{k[0]}@{k[1]}:{c}" for k, c in sorted(copy_edge.items())))
    print("motif_copy_edge_by_strand\t" + ",".join(
        f"{k[0]}|{k[1]}@{k[2]}:{c}" for k, c in sorted(copy_edge_strand.items())))
    print(f"motif_intersect_with_S\t{motif_intersect_s}")
    print(f"motif_intersect_with_S_mapq_lt{UNIQUE_MAPQ}\t{motif_intersect_s_lt10}")
    print(f"motif_intersect_no_S\t{motif_full_total}")
    print(f"motif_intersect_no_S_mapq_lt{UNIQUE_MAPQ}\t{motif_full_lt10}")
    print("motif_intersect_no_S_NM\t" + ",".join(
        f"{k}:{motif_full_nm[k]}" for k in [0, 1, 3]))
    print()
    print("example_junction_softclip_reads (QNAME FLAG POS MAPQ CIGAR):")
    for q, fl, p, m, c in examples_pool[:5]:
        print(f"  {q}\t{fl}\t{p}\t{m}\t{c}")
    if not examples_pool:
        print("  (none)")
    print("example_minus_strand_S_ge14_reads (QNAME FLAG POS MAPQ CIGAR):")
    for q, fl, p, m, c in minus_examples:
        print(f"  {q}\t{fl}\t{p}\t{m}\t{c}")
    if not minus_examples:
        print("  (none)")


if __name__ == "__main__":
    main()
