#!/usr/bin/env python3
"""Read-only BAM/SAM inspection for the synth_is_mob zero-JC diagnosis.

Reads SAM text on stdin (samtools view aligned.bam, or cat work/aligned.sam)
and reports measured counts only. No files are written.

Geometry below matches the FIXED fixture used by job 2371412 (4000 bp,
insert at 601, motif copies 1201-1280 / 2481-2560). The pre-fix fixture
(1440 bp, insert at 201) was diagnosed by an earlier revision of this
script; see git history.

Geometry (testdata/generate.sh synth_is_mob, post-fix):
  ref = left(1-1200) + motif(1201-1280) + mid(1281-2480) + motif(2481-2560)
        + right(2561-4000)
  edited genome inserts a third motif copy between ref positions 600 and 601
  (insert_at_1based = 601).
  Measured register coincidences (ref vs motif ends):
    motif[:2] == ref[601:602]  ("GT") -> L|M junction may key at 600 or 602
    motif[-2:] == ref[599:600] ("GC") -> M|R junction may key at 601 or 599
  breseq 0.40.2 (job 2371412) reports exactly these two registers:
    JC 599 (+1) -> 2558 (-1)  and  JC 602 (-1) -> 2483 (+1), both freq ~0.6.

With --ref, softclips of reads near the junction are additionally placed on
the reference with the same rules as the engine (exact match, both strands,
trivial continuation skipped) so the side2 landing sites can be compared
against the motif copies.
"""

import argparse
import re
import sys
from collections import Counter

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")

MIN_CLIP = 14  # engine MIN_CLIP_FOR_JC
UNIQUE_MAPQ = 10  # engine UNIQUE_MAPQ

J_END = 600       # nominal 1-based left flank end (insert between 600 and 601)
J_START = 601     # nominal 1-based right flank start
J_WIN = (595, 610)  # diagnosis window around the junction (task-specified)
NEAR = 5
MOTIF_COPIES = [(1201, 1280), (2481, 2560)]
COPY_EDGE_ENDS = {1200, 1280, 1281, 2480, 2560, 2561}
COPY_EDGE_STARTS = {1201, 1280, 1281, 2481, 2560, 2561}
INDEL_WINDOW = (501, 700)


def parse_cigar(cigar):
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def lead_trail_s(ops):
    """(leading S len, trailing S len), tolerating hard clips."""
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


def rc_dna(seq):
    return seq.translate(_RC_TABLE)[::-1]


_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def find_exact(hay, needle):
    out = []
    start = 0
    while True:
        i = hay.find(needle, start)
        if i < 0:
            return out
        out.append(i)
        start = i + 1


def copy_label(pos1):
    for i, (s, e) in enumerate(MOTIF_COPIES, 1):
        if s <= pos1 <= e:
            return f"copy{i}[{s}-{e}]"
    return "outside_copies"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="stdin")
    ap.add_argument("--ref", default=None,
                    help="reference FASTA for engine-style clip placement")
    args = ap.parse_args()
    label = args.label

    ref_seq = None
    if args.ref:
        with open(args.ref) as fh:
            ref_seq = "".join(
                line.strip() for line in fh if not line.startswith(">")
            ).upper()
        print(f"loaded_ref\t{args.ref}\t{len(ref_seq)} bp")

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
    end_j = []
    start_j = []
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

    # Junction-window diagnosis (clip >= MIN_CLIP, aligned block end/start in
    # J_WIN). Engine reference-oriented side1 key per strand, plus clip
    # placement onto the reference when --ref is given.
    win_raw = Counter()         # (strand, cigar_side, pos)
    win_key = Counter()         # (strand, engine side1 key)
    win_mapq = Counter()        # (strand, mapq)
    win_side2 = Counter()       # (strand, side2_pos_1)
    win_side2_copy = Counter()  # (strand, copy label)
    win_unplaced = Counter()    # strand
    win_examples = []

    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        qname, flag_s, _rname, pos_s, mapq_s, cigar = f[0], f[1], f[2], f[3], f[4], f[5]
        seq = f[9]
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
        lead_s, trail_s = lead_trail_s(ops)
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
                if trail_s:
                    clip_pos[("end", ref_end)] += 1
                    clip_pos_strand[(strand, "end", ref_end)] += 1
                if lead_s:
                    clip_pos[("start", pos)] += 1
                    clip_pos_strand[(strand, "start", pos)] += 1

                is_j = (trail_s and ref_end == J_END) or (lead_s and pos == J_START)
                if trail_s and ref_end == J_END:
                    end_j.append((qname, flag, pos, mapq, cigar, ops[-1][0]))
                if lead_s and pos == J_START:
                    start_j.append((qname, flag, pos, mapq, cigar, ops[0][0]))
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

                # Junction-window sides: aligned block end (trailing S) or
                # start (leading S) inside J_WIN.
                win_sides = []
                if trail_s and J_WIN[0] <= ref_end <= J_WIN[1]:
                    win_sides.append(("end", ref_end, trail_s))
                if lead_s and J_WIN[0] <= pos <= J_WIN[1]:
                    win_sides.append(("start", pos, lead_s))
                for side, apos, slen in win_sides:
                    win_raw[(strand, side, apos)] += 1
                    win_mapq[(strand, mapq)] += 1
                    # Engine reference-oriented side1 key (pileup.rs):
                    # trailing S: last match if plus else first match;
                    # leading S: first match if plus else last match.
                    if side == "end":
                        key = ref_end if not minus else pos
                        clip = seq[len(seq) - slen:] if seq != "*" else ""
                    else:
                        key = pos if not minus else ref_end
                        clip = seq[:slen] if seq != "*" else ""
                    clip_is_left = (side == "start") != minus  # lead_s XOR minus
                    win_key[(strand, key)] += 1
                    if len(win_examples) < 12:
                        win_examples.append(
                            (qname, flag, pos, mapq, cigar, side, key))
                    if ref_seq and len(clip) >= MIN_CLIP:
                        m_len = sum(n for n, op in ops if op in "M=X")
                        if m_len < MIN_CLIP:
                            continue
                        placed = []
                        for q, is_rc in ((clip.upper(), False),
                                         (rc_dna(clip.upper()), True)):
                            for hit0 in find_exact(ref_seq, q):
                                hlen = len(q)
                                # engine is_continuation: skip the trivial
                                # continuation of the same alignment
                                if clip_is_left:
                                    hend1 = hit0 + hlen
                                    if hend1 == key or hend1 + 1 == key:
                                        continue
                                else:
                                    hstart1 = hit0 + 1
                                    if hstart1 == key + 1 or hstart1 == key:
                                        continue
                                content_rc_hit = is_rc != minus
                                if clip_is_left == content_rc_hit:
                                    side2 = hit0 + 1
                                else:
                                    side2 = hit0 + hlen
                                placed.append(side2)
                        if placed:
                            for s2 in placed:
                                win_side2[(strand, s2)] += 1
                                win_side2_copy[(strand, copy_label(s2))] += 1
                        else:
                            win_unplaced[strand] += 1

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
    print(f"junction_end_at_{J_END}\t{len(end_j)}")
    print(f"junction_start_at_{J_START}\t{len(start_j)}")
    print(f"near_miss_end_{J_END - NEAR}-{J_END + NEAR}_excl{J_END}\t"
          + ",".join(f"{p}:{c}" for p, c in sorted(near_end.items())))
    print(f"near_miss_start_{J_START - NEAR}-{J_START + NEAR}_excl{J_START}\t"
          + ",".join(f"{p}:{c}" for p, c in sorted(near_start.items())))
    print(f"junction_mapq_lt{UNIQUE_MAPQ}\t{j_mapq_lt10}")
    print("junction_mapq_values\t" + ",".join(f"{m}:{c}" for m, c in sorted(j_mapq.items())))
    print(f"junction_strand_mapq_ge{UNIQUE_MAPQ}\t" + ",".join(f"{s}:{c}" for s, c in sorted(j_strand_ge10.items())))
    print()
    print(f"== junction window {J_WIN[0]}-{J_WIN[1]} (clip>={MIN_CLIP}) ==")
    print("win_raw_by_strand_side_pos\t" + ",".join(
        f"{k[0]}|{k[1]}@{k[2]}:{c}" for k, c in sorted(win_raw.items())))
    print("win_engine_side1_key_by_strand\t" + ",".join(
        f"{k[0]}@{k[1]}:{c}" for k, c in sorted(win_key.items())))
    print("win_mapq_by_strand\t" + ",".join(
        f"{k[0]}|mq{k[1]}:{c}" for k, c in sorted(win_mapq.items())))
    if ref_seq:
        print("win_clip_side2_by_strand\t" + ",".join(
            f"{k[0]}@{k[1]}:{c}" for k, c in sorted(win_side2.items())))
        print("win_clip_side2_copy_by_strand\t" + ",".join(
            f"{k[0]}|{k[1]}:{c}" for k, c in sorted(win_side2_copy.items())))
        print("win_clip_unplaced_by_strand\t" + ",".join(
            f"{s}:{c}" for s, c in sorted(win_unplaced.items())))
    print("win_examples (QNAME FLAG POS MAPQ CIGAR side engine_key):")
    for q, fl, p, m, c, side, key in win_examples:
        print(f"  {q}\t{fl}\t{p}\t{m}\t{c}\t{side}\tkey={key}")
    if not win_examples:
        print("  (none)")
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
