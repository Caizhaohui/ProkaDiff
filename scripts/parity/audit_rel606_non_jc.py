#!/usr/bin/env python3
"""audit_rel606_non_jc.py — M8-T03 REL606 non-JC (SNP/INS/DEL/MOB) redline parity audit.

Compares ProkaDiff consensus (non-JC) REDLINE mutation calls against the
breseq 0.40.2 oracle GD on the REL606 clonal sample, per docs/parity.md.

Strict identity contract (FIX-M8-T03-AUDIT-SCRIPT-RIGOR, Codex Sol review):

  * A ProkaDiff call matches a breseq call iff kind, normalized seq_id
    (trailing GenBank version suffix stripped: 'REL606.6' <-> 'REL606'),
    position (tolerance --tol-bp, default 0 = exact) and ALL kind-specific
    core attributes are equal:
      - SNP: new_seq (alt allele)
      - INS: new_seq (inserted allele)
      - DEL: size (deleted length)
      - MOB: repeat_name AND strand AND duplication_size
  * Any attribute difference is NOT a match: the ProkaDiff call counts as
    OVER and the breseq call as UNDER, and the audit exits non-zero.
  * Matching is strict one-to-one: one ProkaDiff record consumes at most one
    breseq record (duplicate identities on either side remain OVER/UNDER).
    The pairing is a true maximum-cardinality bipartite matching (augmenting
    paths / Kuhn's algorithm; FIX-M8-T03-MAXIMUM-BIPARTITE-MATCHING, Codex
    Sol review) — plain greedy closest-first selection is not maximum when
    tol > 0. The augmenting-path search is fully ITERATIVE (BFS with a
    collections.deque and predecessor maps; FIX-M8-T03-ITERATIVE-AUGMENTING-
    PATH, Codex Sol re-review) so long chains (e.g. 1050 SNPs with tol=1)
    cannot hit CPython's recursion limit.
  * Record accounting is by positional index, never by GD id
    (FIX-M8-T03-RECORD-ACCOUNTING-GAP, Codex Sol re-review): GD ids are not
    guaranteed unique, and id-keyed accounting silently dropped duplicate-id
    records. Every parsed record is strictly accounted for as MATCHED,
    OVER (unmatched ProkaDiff) or UNDER (unmatched breseq). parse_gd()
    additionally warns on stderr when duplicate record ids are present.
  * --tol-bp > 0 is a diagnostic escape hatch only: it may relax position by
    at most N bp but never relaxes attribute identity. The official non-JC
    parity artifact (benchmark/results/m8_rel606_non_jc_parity.tsv) is
    generated with the default tol=0; the non-JC parity contract requires
    exact identity.

REDLINE classes evaluated (consensus records only): SNP, INS, DEL, MOB.
Evidence-only MC and RA records are ignored for redline matching; their
record counts are reported for comparison (breseq UN records are whitelisted
evidence and likewise not matched). JC is out of scope (separate task).
OVER/UNDER rows carry a nearest-counterpart note explaining the first
differing attribute (diagnostic text only; it never affects matching).

Usage:
  python3 scripts/parity/audit_rel606_non_jc.py [--tol-bp 0]

Writes benchmark/results/m8_rel606_non_jc_parity.tsv with columns:
  mutation_type, contig, pd_pos, breseq_pos, pd_allele_or_elem,
  breseq_allele_or_elem, delta_bp, status (MATCHED/OVER/UNDER), notes

Exit codes:
  0  PASS (over == 0 and under == 0)
  1  parity FAIL (over > 0 or under > 0)
  2  input error (missing or malformed GD)
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, deque
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_PD_GD = REPO_ROOT / "benchmark/results/clonal_2769249/prokadiff/output.gd"
FALLBACK_PD_GD = REPO_ROOT / (
    "benchmark/results/parity_closure/g1_clonal_jc/slurm/3066886/prokadiff/output.gd"
)
DEFAULT_BRESEQ_GD = REPO_ROOT / "benchmark/results/clonal_2769249/breseq/output/output.gd"
DEFAULT_OUT = REPO_ROOT / "benchmark/results/m8_rel606_non_jc_parity.tsv"

REDLINE_KINDS = ("SNP", "INS", "DEL", "MOB")
KIND_RANK = {k: i for i, k in enumerate(REDLINE_KINDS)}
EVIDENCE_KINDS = ("MC", "RA", "UN")

CONTIG_VERSION_SUFFIX = re.compile(r"\.\d+$")


def normalize_contig(name: str) -> str:
    """Strip a trailing GenBank version suffix: 'REL606.6' -> 'REL606'."""
    return CONTIG_VERSION_SUFFIX.sub("", name.strip())


def parse_gd(path: Path) -> tuple[dict[str, list[dict]], Counter]:
    """Parse a GenomeDiff 1.0 file.

    Returns (records_by_kind, counts) where records_by_kind holds parsed
    core fields for REDLINE kinds and counts tallies every record kind
    (including evidence-only MC/RA/UN and JC).

    Duplicate record ids among REDLINE records do not abort the parse
    (matching accounts for records by positional index, not id), but they
    violate the GD 1.0 uniqueness expectation and are reported as a
    WARNING on stderr (FIX-M8-T03-RECORD-ACCOUNTING-GAP).
    """
    records: dict[str, list[dict]] = {k: [] for k in REDLINE_KINDS}
    counts: Counter = Counter()
    if not path.is_file():
        raise FileNotFoundError(f"GD file not found: {path}")

    seen_ids: set[str] = set()
    dup_id_counts: Counter = Counter()
    for line in path.read_text(errors="replace").splitlines():
        ls = line.strip()
        if not ls or ls.startswith("#") or ls.startswith("="):
            continue
        parts = ls.split("\t")
        kind = parts[0]
        counts[kind] += 1
        if kind not in REDLINE_KINDS:
            continue
        try:
            rec = {
                "kind": kind,
                "id": parts[1],
                "contig": normalize_contig(parts[3]),
                "pos": int(parts[4]),
                "raw": ls,
            }
            if kind in ("SNP", "INS"):
                rec["allele"] = parts[5].upper()
            elif kind == "DEL":
                rec["size"] = int(parts[5])
            elif kind == "MOB":
                rec["allele"] = parts[5]
                rec["strand"] = int(parts[6])
                rec["dup"] = int(parts[7])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"Malformed {kind} record in {path}: {ls!r}") from exc
        if rec["id"] in seen_ids:
            dup_id_counts[rec["id"]] += 1
        else:
            seen_ids.add(rec["id"])
        records[kind].append(rec)

    if dup_id_counts:
        detail = ", ".join(
            f"'{rid}' x{1 + n}" for rid, n in sorted(dup_id_counts.items())
        )
        print(
            f"WARNING: {path}: duplicate REDLINE record ids (GD 1.0 ids should "
            f"be unique): {detail}; each record is still accounted for "
            "individually by positional index",
            file=sys.stderr,
        )

    for kind in REDLINE_KINDS:
        records[kind].sort(key=lambda r: (r["pos"], r["id"]))
    return records, counts


def elem_label(rec: dict) -> str:
    """Human-readable allele/element label for the TSV."""
    if rec["kind"] in ("SNP", "INS"):
        return rec["allele"]
    if rec["kind"] == "DEL":
        return f"del_size={rec['size']}"
    return f"{rec['allele']}(strand={rec['strand']},dup={rec['dup']})"


def canonical_signature(rec: dict) -> tuple:
    """Content signature (no ids/parents) for fallback cross-check."""
    return (rec["kind"], rec["contig"], rec["pos"], elem_label(rec))


def candidate_ok(pd_rec: dict, br_rec: dict, tol: int) -> tuple[bool, int]:
    """Strict identity test for one candidate pair; returns (ok, delta_bp).

    Match requires: same kind (implicit), same normalized contig, position
    within `tol` (default 0 = exact), and equal kind-specific core
    attributes (SNP/INS new_seq, DEL size, MOB repeat_name + strand +
    duplication_size). Attribute differences are never tolerated, whatever
    `tol` is.
    """
    delta = pd_rec["pos"] - br_rec["pos"]
    if pd_rec["contig"] != br_rec["contig"]:
        return False, delta
    if abs(delta) > tol:
        return False, delta
    kind = pd_rec["kind"]
    if kind in ("SNP", "INS"):
        return pd_rec["allele"] == br_rec["allele"], delta
    if kind == "DEL":
        return pd_rec["size"] == br_rec["size"], delta
    # MOB: repeat_name, strand and duplication_size must all match exactly.
    return (
        pd_rec["allele"] == br_rec["allele"]
        and pd_rec["strand"] == br_rec["strand"]
        and pd_rec["dup"] == br_rec["dup"]
    ), delta


def match_kind(pd_recs: list[dict], br_recs: list[dict], tol: int) -> dict:
    """Deterministic maximum one-to-one matching for one mutation kind.

    True MAXIMUM-CARDINALITY bipartite matching via augmenting paths
    (Kuhn's algorithm), fixing FIX-M8-T03-MAXIMUM-BIPARTITE-MATCHING
    (Codex Sol review): greedy closest-first selection was not maximum
    when tol > 0 (counterexample: SNPs at pd positions 2, 3 vs breseq
    positions 1, 2 with tol=1 greedily took (2, 2), stranding pd 3 and
    reporting matched=1 instead of the feasible maximum matched=2).

    Left vertices are pd indices, right vertices br indices; an edge
    (i, j) exists iff candidate_ok(pd_recs[i], br_recs[j], tol)[0].
    Each left vertex's adjacency is sorted by (abs(delta), breseq pos,
    br index) so ties prefer the closer counterpart, while the
    augmenting-path search still guarantees maximum cardinality.

    Consumed records are tracked by their positional index in pd_recs /
    br_recs, never by GD id: ids are not guaranteed unique, and id-keyed
    accounting silently dropped duplicate-id records from OVER/UNDER
    (FIX-M8-T03-RECORD-ACCOUNTING-GAP, Codex Sol re-review). Every input
    record is strictly accounted for: over = pd indices never matched,
    under = br indices never matched.

    The augmenting-path search is ITERATIVE (FIX-M8-T03-ITERATIVE-
    AUGMENTING-PATH, Codex Sol re-review): a breadth-first search over an
    explicit collections.deque with predecessor maps replaces the previous
    recursive try_augment() DFS, whose call chain followed the augmenting
    path and exceeded CPython's recursion limit on long chains (n=1050
    SNPs with tol=1 raised RecursionError). For each unmatched pd vertex
    the BFS explores alternating pd -> br -> pd layers; when it first
    reaches a free br vertex the path is flipped by walking prev_pd /
    match_pd backwards. No recursion occurs, so arbitrarily large inputs
    are safe.
    """
    adj: list[list[int]] = []
    for pd_rec in pd_recs:
        edges = []
        for br_i, br_rec in enumerate(br_recs):
            ok, delta = candidate_ok(pd_rec, br_rec, tol)
            if ok:
                edges.append((abs(delta), br_rec["pos"], br_i))
        edges.sort()
        adj.append([br_i for _, _, br_i in edges])

    match_br: list[int] = [-1] * len(br_recs)  # br index -> matched pd index
    match_pd: dict[int, int] = {}  # pd index -> matched br index

    for root_pd in range(len(pd_recs)):
        if root_pd in match_pd:
            continue  # already matched by an earlier augmenting path
        # BFS for an augmenting path starting at the free pd vertex
        # root_pd. prev_pd maps each discovered br index to the pd vertex
        # that discovered it; prev_br maps each enqueued (non-root) pd
        # vertex to the br index it was reached through.
        prev_pd: dict[int, int] = {}
        prev_br: dict[int, int] = {}
        free_br = None
        queue = deque([root_pd])
        while queue and free_br is None:
            pd_u = queue.popleft()
            for br_j in adj[pd_u]:
                if br_j in prev_pd:
                    continue  # br vertex already discovered in this BFS
                prev_pd[br_j] = pd_u
                owner = match_br[br_j]
                if owner == -1:
                    free_br = br_j  # free right vertex: augmenting path found
                    break
                prev_br[owner] = br_j
                queue.append(owner)
        if free_br is None:
            continue  # no augmenting path from root_pd; it stays OVER
        # Flip the augmenting path: walk predecessors from the free br
        # vertex back to root_pd, re-matching every pd/br pair on the way.
        br_j = free_br
        while br_j is not None:
            pd_u = prev_pd[br_j]
            next_br = match_pd.get(pd_u)
            match_pd[pd_u] = br_j
            match_br[br_j] = pd_u
            br_j = next_br

    inv = {pd_i: br_i for br_i, pd_i in enumerate(match_br) if pd_i != -1}
    used_pd_idx = set(inv)
    used_br_idx = set(inv.values())
    matched = [(pd_recs[pd_i], br_recs[inv[pd_i]]) for pd_i in sorted(inv)]
    over = [rec for i, rec in enumerate(pd_recs) if i not in used_pd_idx]
    under = [rec for i, rec in enumerate(br_recs) if i not in used_br_idx]
    return {"matched": matched, "over": over, "under": under}


def match_notes(pd_rec: dict, br_rec: dict, tol: int) -> str:
    """Notes for a MATCHED pair (attributes are identical by construction)."""
    delta = pd_rec["pos"] - br_rec["pos"]
    if delta == 0:
        return "identity_exact (pos and all core attributes equal)"
    return (
        f"attributes_exact; pos delta {delta:+d}bp within diagnostic tol {tol}"
    )


def attribute_diffs(pd_rec: dict, br_rec: dict) -> list[str]:
    """List the kind-specific core attribute differences (diagnostic only)."""
    diffs = []
    kind = pd_rec["kind"]
    if kind in ("SNP", "INS"):
        if pd_rec["allele"] != br_rec["allele"]:
            diffs.append(
                f"allele pd={pd_rec['allele']} breseq={br_rec['allele']}"
            )
    elif kind == "DEL":
        if pd_rec["size"] != br_rec["size"]:
            diffs.append(f"size pd={pd_rec['size']} breseq={br_rec['size']}")
    else:  # MOB
        if pd_rec["allele"] != br_rec["allele"]:
            diffs.append(
                f"repeat_name pd={pd_rec['allele']} breseq={br_rec['allele']}"
            )
        if pd_rec["strand"] != br_rec["strand"]:
            diffs.append(
                f"strand pd={pd_rec['strand']} breseq={br_rec['strand']}"
            )
        if pd_rec["dup"] != br_rec["dup"]:
            diffs.append(f"dup_size pd={pd_rec['dup']} breseq={br_rec['dup']}")
    return diffs


def nearest_counterpart_note(rec: dict, opposing: list[dict]) -> str:
    """Describe the closest same-kind, same-contig opposing record.

    Diagnostic text for OVER/UNDER rows only; never affects matching.
    """
    best = None
    for other in opposing:
        if other["contig"] != rec["contig"]:
            continue
        key = (abs(rec["pos"] - other["pos"]), other["pos"])
        if best is None or key < best[0]:
            best = (key, other)
    if best is None:
        return "no same-contig same-kind call on the other side"
    other = best[1]
    bits = [f"nearest Δ={rec['pos'] - other['pos']:+d}bp (pos {other['pos']})"]
    bits.extend(attribute_diffs(rec, other))
    return "; ".join(bits)


def build_rows(pd_records: dict, br_records: dict, tol: int) -> tuple[list, dict]:
    rows, per_type = [], {}
    for kind in REDLINE_KINDS:
        res = match_kind(pd_records[kind], br_records[kind], tol)
        for pd_rec, br_rec in res["matched"]:
            rows.append({
                "mutation_type": kind,
                "contig": pd_rec["contig"],
                "pd_pos": str(pd_rec["pos"]),
                "breseq_pos": str(br_rec["pos"]),
                "pd_allele_or_elem": elem_label(pd_rec),
                "breseq_allele_or_elem": elem_label(br_rec),
                "delta_bp": str(pd_rec["pos"] - br_rec["pos"]),
                "status": "MATCHED",
                "notes": match_notes(pd_rec, br_rec, tol),
            })
        for rec in res["over"]:
            rows.append({
                "mutation_type": kind,
                "contig": rec["contig"],
                "pd_pos": str(rec["pos"]),
                "breseq_pos": ".",
                "pd_allele_or_elem": elem_label(rec),
                "breseq_allele_or_elem": ".",
                "delta_bp": ".",
                "status": "OVER",
                "notes": "OVER: no exact breseq oracle counterpart ("
                         + nearest_counterpart_note(rec, br_records[kind]) + ")",
            })
        for rec in res["under"]:
            rows.append({
                "mutation_type": kind,
                "contig": rec["contig"],
                "pd_pos": ".",
                "breseq_pos": str(rec["pos"]),
                "pd_allele_or_elem": ".",
                "breseq_allele_or_elem": elem_label(rec),
                "delta_bp": ".",
                "status": "UNDER",
                "notes": "UNDER: no exact ProkaDiff call ("
                         + nearest_counterpart_note(rec, pd_records[kind]) + ")",
            })
        per_type[kind] = (len(pd_records[kind]), len(br_records[kind]), len(res["matched"]))
    rows.sort(key=lambda r: (KIND_RANK[r["mutation_type"]], r["contig"], r["pd_pos"]))
    return rows, per_type


def write_tsv(rows: list, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "mutation_type", "contig", "pd_pos", "breseq_pos", "pd_allele_or_elem",
        "breseq_allele_or_elem", "delta_bp", "status", "notes",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row[col] for col in header) + "\n")


def fallback_check(fallback_path: Path, pd_records: dict) -> list[str]:
    """Compare the primary ProkaDiff GD against the G1 baseline fallback."""
    if not fallback_path.is_file():
        return [f"Fallback check: {fallback_path} missing (skipped)"]
    fb_records, _ = parse_gd(fallback_path)
    primary = Counter(
        canonical_signature(r) for kind in REDLINE_KINDS for r in pd_records[kind]
    )
    fallback = Counter(
        canonical_signature(r) for kind in REDLINE_KINDS for r in fb_records[kind]
    )
    fb_only = fallback - primary
    pd_only = primary - fallback
    total_fb = sum(len(fb_records[k]) for k in REDLINE_KINDS)
    lines = [f"Fallback check (G1 baseline {fallback_path.parent.name}): "
             f"non-JC REDLINE calls={total_fb}"]
    if not fb_only and not pd_only:
        lines.append("  Primary vs fallback: IDENTICAL non-JC content")
    else:
        lines.append("  Primary vs fallback: DIFFER")
        for sig, n in sorted(fb_only.items()):
            lines.append(f"  fallback-only x{n}: {sig[0]} {sig[1]}:{sig[2]} {sig[3]}")
        for sig, n in sorted(pd_only.items()):
            lines.append(f"  primary-only x{n}: {sig[0]} {sig[1]}:{sig[2]} {sig[3]}")
    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pd-gd", type=Path, default=DEFAULT_PD_GD,
                    help="ProkaDiff output.gd (primary)")
    ap.add_argument("--fallback-gd", type=Path, default=FALLBACK_PD_GD,
                    help="Fallback ProkaDiff output.gd (G1 baseline 3066886)")
    ap.add_argument("--breseq-gd", type=Path, default=DEFAULT_BRESEQ_GD,
                    help="breseq oracle output.gd")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="Output TSV path")
    ap.add_argument("--tol-bp", type=int, default=0,
                    help="Position tolerance in bp (default 0 = exact, as the "
                         "non-JC parity contract requires; >0 is a diagnostic "
                         "override that never relaxes attribute identity)")
    args = ap.parse_args()

    if args.tol_bp > 0:
        print(
            f"WARNING: --tol-bp {args.tol_bp} > 0 is a diagnostic override; "
            "the non-JC parity contract requires exact positions (tol=0).",
            file=sys.stderr,
        )

    try:
        pd_records, pd_counts = parse_gd(args.pd_gd)
        br_records, br_counts = parse_gd(args.breseq_gd)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    rows, per_type = build_rows(pd_records, br_records, args.tol_bp)
    write_tsv(rows, args.out)

    n_pd = sum(len(pd_records[k]) for k in REDLINE_KINDS)
    n_br = sum(len(br_records[k]) for k in REDLINE_KINDS)
    n_matched = sum(1 for r in rows if r["status"] == "MATCHED")
    n_over = sum(1 for r in rows if r["status"] == "OVER")
    n_under = sum(1 for r in rows if r["status"] == "UNDER")

    print("M8-T03 REL606 non-JC (SNP/INS/DEL/MOB) parity audit")
    print(f"ProkaDiff GD : {args.pd_gd}")
    print(f"breseq GD    : {args.breseq_gd}")
    print(f"Output TSV   : {args.out} ({len(rows)} rows)")
    for kind in REDLINE_KINDS:
        p, b, m = per_type[kind]
        print(f"Per-type {kind}: pd={p} breseq={b} matched={m}")
    print(f"Total ProkaDiff non-JC REDLINE calls: {n_pd}")
    print(f"Total breseq non-JC REDLINE calls: {n_br}")
    print(f"Matched: {n_matched}")
    print(f"Over: {n_over}")
    print(f"Under: {n_under}")
    print(
        "MC (evidence-only, not redline-matched): "
        f"prokadiff={pd_counts.get('MC', 0)} breseq={br_counts.get('MC', 0)}"
    )
    print(
        "Other evidence records ignored (docs/parity.md whitelist): "
        f"prokadiff RA={pd_counts.get('RA', 0)}; "
        f"breseq RA={br_counts.get('RA', 0)} UN={br_counts.get('UN', 0)}; "
        f"JC pd={pd_counts.get('JC', 0)} breseq={br_counts.get('JC', 0)} (JC parity is a separate task)"
    )
    for line in fallback_check(args.fallback_gd, pd_records):
        print(line)

    if n_over == 0 and n_under == 0:
        print(
            f"Result: PASS (matched={n_matched}, over=0, under=0, "
            "strict identity: pos + full core attributes)"
        )
        return 0
    print(f"Result: FAIL (over={n_over}, under={n_under})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
