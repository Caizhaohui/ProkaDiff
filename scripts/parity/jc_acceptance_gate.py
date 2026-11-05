#!/usr/bin/env python3
"""Formal accepted-JC bidirectional gate for ProkaDiff ↔ breseq parity closure.

This is the project acceptance gate for JC. It is independent of
``scripts/layer2_gd_compare.py`` (diagnostic; exit 0 on red-line only).

Contract (docs/parity.md + PROJECT_GOAL_FREEZE):
  - Compare accepted JC only (breseq rows with nonempty reject= are NOT truth).
  - Deterministic one-to-one greedy match after side canonicalization.
  - Require contig + strand identity on both sides; |Δpos| ≤ jc_tol_bp (default 5).
  - Require signed overlap equality (core identity).
  - Optional --normalize-contig-version: strip a trailing GenBank VERSION suffix
    (``\\.\\d+$``) so LOCUS ``REL606`` matches ProkaDiff VERSION ``REL606.6``
    when comparing the same molecule (RW-005 naming vs breseq LOCUS). Documented;
    does not invent coordinate matches.
  - Optional --allow-repeat-copy-equivalence: if both sides annotate one breakend
    as repeat-redundant (breseq side_*_redundant=1 / annotate_key=repeat) OR the
    unmatched pair shares exactly one unique (non-repeat) flank within tol while
    the other flank differs beyond tol, emit EQUIV_REPEAT_COPY (pre-approved
    whitelist class) — never auto-promotes unknown leftovers to PASS.
  - Machine-readable JSON + TSV failure reasons.
  - Exit 0 only when accepted over/under unmatched == 0 after whitelist classes.
  - Exit 2 on accepted JC leftover; exit 1 on I/O / usage errors.

Does not modify historical layer2_gd_compare semantics.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

JC_TOL_DEFAULT = 5
_VERSION_SUFFIX = re.compile(r"\.\d+$")


def normalize_contig_id(seq_id: str, strip_version: bool) -> str:
    """Optionally strip GenBank VERSION numeric suffix (REL606.6 → REL606)."""
    if not strip_version:
        return seq_id
    return _VERSION_SUFFIX.sub("", seq_id)


@dataclass
class JcEvent:
    source: str
    record_id: str
    seq1: str
    pos1: int
    strand1: str
    seq2: str
    pos2: int
    strand2: str
    overlap: int
    rejected: bool
    reject_reasons: str
    attrs: dict[str, str] = field(default_factory=dict)
    raw_line: str = ""
    index: int = -1

    def side_a(self) -> tuple[str, str, int]:
        return (self.seq1, self.strand1, self.pos1)

    def side_b(self) -> tuple[str, str, int]:
        return (self.seq2, self.strand2, self.pos2)

    def canonical(self) -> tuple[tuple[str, str, int], tuple[str, str, int], int]:
        a, b = self.side_a(), self.side_b()
        if a > b:
            a, b = b, a
        return a, b, self.overlap


def _parse_attrs(parts: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for p in parts:
        if "=" in p:
            k, _, v = p.partition("=")
            out[k] = v
    return out


def parse_jc_gd(
    path: Path, source: str, *, strip_contig_version: bool = False
) -> list[JcEvent]:
    events: list[JcEvent] = []
    if not path.is_file():
        raise FileNotFoundError(path)
    for line in path.read_text(errors="replace").splitlines():
        if not line or line.startswith("#") or line.startswith("="):
            continue
        parts = line.split("\t")
        if len(parts) < 10 or parts[0] != "JC":
            continue
        try:
            overlap = int(parts[9])
            ev = JcEvent(
                source=source,
                record_id=parts[1],
                seq1=normalize_contig_id(parts[3], strip_contig_version),
                pos1=int(parts[4]),
                strand1=parts[5],
                seq2=normalize_contig_id(parts[6], strip_contig_version),
                pos2=int(parts[7]),
                strand2=parts[8],
                overlap=overlap,
                rejected=False,
                reject_reasons="",
                attrs=_parse_attrs(parts[10:]),
                raw_line=line,
                index=len(events),
            )
        except (ValueError, IndexError):
            continue
        reject = ev.attrs.get("reject", "")
        if reject:
            ev.rejected = True
            ev.reject_reasons = reject
        events.append(ev)
    return events


def sides_close(
    a: tuple[str, str, int], b: tuple[str, str, int], tol: int
) -> bool:
    return a[0] == b[0] and a[1] == b[1] and abs(a[2] - b[2]) <= tol


def events_match(left: JcEvent, right: JcEvent, tol: int) -> bool:
    la, lb, lo = left.canonical()
    ra, rb, ro = right.canonical()
    if lo != ro:
        return False
    return sides_close(la, ra, tol) and sides_close(lb, rb, tol)


def match_cost(left: JcEvent, right: JcEvent) -> int:
    """Lower is better; used for deterministic min-cost among tol matches."""
    la, lb, _ = left.canonical()
    ra, rb, _ = right.canonical()
    return abs(la[2] - ra[2]) + abs(lb[2] - rb[2])


def one_to_one_match(
    left: list[JcEvent], right: list[JcEvent], tol: int
) -> tuple[list[tuple[JcEvent, JcEvent]], list[JcEvent], list[JcEvent]]:
    """Deterministic 1:1: sort candidate pairs by (cost, left_id, right_id)."""
    candidates: list[tuple[int, str, str, int, int]] = []
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if events_match(a, b, tol):
                candidates.append(
                    (match_cost(a, b), a.record_id, b.record_id, i, j)
                )
    candidates.sort()
    used_l: set[int] = set()
    used_r: set[int] = set()
    pairs: list[tuple[JcEvent, JcEvent]] = []
    for _cost, _lid, _rid, i, j in candidates:
        if i in used_l or j in used_r:
            continue
        used_l.add(i)
        used_r.add(j)
        pairs.append((left[i], right[j]))
    over = [e for i, e in enumerate(left) if i not in used_l]
    under = [e for j, e in enumerate(right) if j not in used_r]
    return pairs, over, under


def unique_flanks(ev: JcEvent) -> list[tuple[str, str, int]]:
    """Prefer breseq annotate_key=gene as unique; else both sides."""
    a_key = ev.attrs.get("side_1_annotate_key", "")
    b_key = ev.attrs.get("side_2_annotate_key", "")
    sides = [ev.side_a(), ev.side_b()]
    if a_key == "gene" and b_key == "repeat":
        return [ev.side_a()]
    if b_key == "gene" and a_key == "repeat":
        return [ev.side_b()]
    if a_key == "repeat" and b_key == "repeat":
        return []
    return sides


def try_repeat_copy_equiv(
    over: list[JcEvent], under: list[JcEvent], tol: int
) -> tuple[list[dict[str, Any]], list[JcEvent], list[JcEvent]]:
    """Pre-approved whitelist: same unique flank + overlap, other flank differs.

    Only applied when oracle (under) marks at least one side as repeat/redundant.
    """
    equivs: list[dict[str, Any]] = []
    used_o: set[int] = set()
    used_u: set[int] = set()
    candidates: list[tuple[int, str, str, int, int]] = []
    for i, o in enumerate(over):
        for j, u in enumerate(under):
            if o.overlap != u.overlap:
                continue
            u_rep = (
                u.attrs.get("side_1_annotate_key") == "repeat"
                or u.attrs.get("side_2_annotate_key") == "repeat"
                or u.attrs.get("side_1_redundant") == "1"
                or u.attrs.get("side_2_redundant") == "1"
            )
            if not u_rep:
                continue
            o_fl = unique_flanks(u)  # unique flank identity from oracle annotation
            if len(o_fl) != 1:
                continue
            uniq = o_fl[0]
            # Require PD to share that unique flank within tol (either PD side)
            pd_hit = any(sides_close(uniq, s, tol) for s in (o.side_a(), o.side_b()))
            if not pd_hit:
                continue
            # And NOT a full event_match (already unmatched)
            if events_match(o, u, tol):
                continue
            cost = min(
                abs(uniq[2] - o.side_a()[2]) if sides_close(uniq, o.side_a(), tol) else 10**9,
                abs(uniq[2] - o.side_b()[2]) if sides_close(uniq, o.side_b(), tol) else 10**9,
            )
            candidates.append((cost, o.record_id, u.record_id, i, j))
    candidates.sort()
    remain_o = list(over)
    remain_u = list(under)
    for _c, _oid, _uid, i, j in candidates:
        if i in used_o or j in used_u:
            continue
        used_o.add(i)
        used_u.add(j)
        equivs.append(
            {
                "class": "EQUIV_REPEAT_COPY",
                "prokadiff_id": over[i].record_id,
                "breseq_id": under[j].record_id,
                "unique_flank": list(unique_flanks(under[j])[0]),
                "prokadiff_sides": [list(over[i].side_a()), list(over[i].side_b())],
                "breseq_sides": [list(under[j].side_a()), list(under[j].side_b())],
                "overlap": over[i].overlap,
                "reason": "shared unique flank within tol; oracle other flank is repeat-annotated; copy placement may differ",
            }
        )
    remain_o = [e for i, e in enumerate(over) if i not in used_o]
    remain_u = [e for j, e in enumerate(under) if j not in used_u]
    return equivs, remain_o, remain_u


def event_brief(e: JcEvent) -> dict[str, Any]:
    return {
        "id": e.record_id,
        "source": e.source,
        "seq1": e.seq1,
        "pos1": e.pos1,
        "strand1": e.strand1,
        "seq2": e.seq2,
        "pos2": e.pos2,
        "strand2": e.strand2,
        "overlap": e.overlap,
        "rejected": e.rejected,
        "reject": e.reject_reasons,
        "mob_evidence": e.attrs.get("mob_evidence", ""),
        "side_1_annotate_key": e.attrs.get("side_1_annotate_key", ""),
        "side_2_annotate_key": e.attrs.get("side_2_annotate_key", ""),
        "canonical": [list(x) if not isinstance(x, int) else x for x in e.canonical()[:2]]
        + [e.overlap],
    }


def run_gate(
    prokadiff_gd: Path,
    breseq_gd: Path,
    tol: int,
    allow_repeat_copy_equiv: bool,
    *,
    strip_contig_version: bool = False,
) -> dict[str, Any]:
    pd_all = parse_jc_gd(
        prokadiff_gd, "prokadiff", strip_contig_version=strip_contig_version
    )
    br_all = parse_jc_gd(
        breseq_gd, "breseq", strip_contig_version=strip_contig_version
    )
    pd_acc = [e for e in pd_all if not e.rejected]
    br_acc = [e for e in br_all if not e.rejected]
    br_rej = [e for e in br_all if e.rejected]

    pairs, over, under = one_to_one_match(pd_acc, br_acc, tol)
    equivs: list[dict[str, Any]] = []
    if allow_repeat_copy_equiv and (over or under):
        equivs, over, under = try_repeat_copy_equiv(over, under, tol)

    # Classify residual PD-only: near a rejected oracle JC?
    near_rejected: list[dict[str, Any]] = []
    still_over: list[JcEvent] = []
    for o in over:
        hit = None
        for r in br_rej:
            if events_match(o, r, tol):
                hit = r
                break
        if hit is not None:
            near_rejected.append(
                {
                    "class": "PD_MATCHES_BRESEQ_REJECTED",
                    "prokadiff_id": o.record_id,
                    "breseq_rejected_id": hit.record_id,
                    "reject": hit.reject_reasons,
                    "note": "ProkaDiff accepted a JC that breseq rejected; not oracle truth",
                }
            )
        else:
            still_over.append(o)

    result = {
        "gate": "JC_ACCEPTANCE_GATE",
        "jc_tol_bp": tol,
        "allow_repeat_copy_equivalence": allow_repeat_copy_equiv,
        "normalize_contig_version": strip_contig_version,
        "prokadiff_gd": str(prokadiff_gd),
        "breseq_gd": str(breseq_gd),
        "counts": {
            "prokadiff_jc_total": len(pd_all),
            "prokadiff_jc_accepted": len(pd_acc),
            "breseq_jc_total": len(br_all),
            "breseq_jc_accepted": len(br_acc),
            "breseq_jc_rejected": len(br_rej),
            "matched_exact_tol": len(pairs),
            "equiv_repeat_copy": len(equivs),
            "pd_matches_breseq_rejected": len(near_rejected),
            "accepted_over_unmatched": len(still_over),
            "accepted_under_unmatched": len(under),
        },
        "matched_pairs": [
            {
                "prokadiff_id": a.record_id,
                "breseq_id": b.record_id,
                "cost": match_cost(a, b),
                "prokadiff": event_brief(a),
                "breseq": event_brief(b),
            }
            for a, b in pairs
        ],
        "equiv_repeat_copy": equivs,
        "pd_matches_breseq_rejected": near_rejected,
        "accepted_over_unmatched": [event_brief(e) for e in still_over],
        "accepted_under_unmatched": [event_brief(e) for e in under],
        "breseq_rejected_inventory": [event_brief(e) for e in br_rej],
        "pass": len(still_over) == 0 and len(under) == 0,
    }
    return result


def write_tsv(result: dict[str, Any], path: Path) -> None:
    rows: list[str] = []
    rows.append(
        "row_class\tprokadiff_id\tbreseq_id\tstatus\tpos_summary\toverlap\treject\treason"
    )

    def pos_sum(brief: dict[str, Any] | None) -> str:
        if not brief:
            return ""
        return (
            f"{brief['seq1']}:{brief['pos1']}({brief['strand1']})"
            f"->{brief['seq2']}:{brief['pos2']}({brief['strand2']})"
        )

    for p in result["matched_pairs"]:
        rows.append(
            "\t".join(
                [
                    "MATCHED",
                    p["prokadiff_id"],
                    p["breseq_id"],
                    "accepted_match",
                    pos_sum(p["prokadiff"]),
                    str(p["prokadiff"]["overlap"]),
                    "",
                    f"cost={p['cost']}",
                ]
            )
        )
    for e in result["equiv_repeat_copy"]:
        rows.append(
            "\t".join(
                [
                    "EQUIV_REPEAT_COPY",
                    e["prokadiff_id"],
                    e["breseq_id"],
                    "whitelist",
                    "",
                    str(e["overlap"]),
                    "",
                    e["reason"],
                ]
            )
        )
    for e in result["pd_matches_breseq_rejected"]:
        rows.append(
            "\t".join(
                [
                    "PD_MATCHES_BRESEQ_REJECTED",
                    e["prokadiff_id"],
                    e["breseq_rejected_id"],
                    "not_truth",
                    "",
                    "",
                    e["reject"],
                    e["note"],
                ]
            )
        )
    for e in result["accepted_over_unmatched"]:
        rows.append(
            "\t".join(
                [
                    "ACCEPTED_OVER",
                    e["id"],
                    "",
                    "FAIL",
                    pos_sum(e),
                    str(e["overlap"]),
                    "",
                    "prokadiff accepted JC with no accepted breseq match",
                ]
            )
        )
    for e in result["accepted_under_unmatched"]:
        rows.append(
            "\t".join(
                [
                    "ACCEPTED_UNDER",
                    "",
                    e["id"],
                    "FAIL",
                    pos_sum(e),
                    str(e["overlap"]),
                    e.get("reject", ""),
                    "breseq accepted JC missing from prokadiff",
                ]
            )
        )
    for e in result["breseq_rejected_inventory"]:
        rows.append(
            "\t".join(
                [
                    "BRESEQ_REJECTED",
                    "",
                    e["id"],
                    "excluded_from_truth",
                    pos_sum(e),
                    str(e["overlap"]),
                    e["reject"],
                    "not mixed into accepted truth",
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prokadiff-gd", type=Path, required=True)
    ap.add_argument("--breseq-gd", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--jc-tol-bp", type=int, default=JC_TOL_DEFAULT)
    ap.add_argument(
        "--allow-repeat-copy-equivalence",
        action="store_true",
        help="Apply pre-approved EQUIV_REPEAT_COPY whitelist (docs/parity.md)",
    )
    ap.add_argument(
        "--strict-no-whitelist",
        action="store_true",
        help="Disable repeat-copy equivalence even if flag set",
    )
    ap.add_argument(
        "--normalize-contig-version",
        action="store_true",
        help="Strip trailing .N VERSION suffix so LOCUS/VERSION seq_ids compare equal",
    )
    args = ap.parse_args()
    allow = args.allow_repeat_copy_equivalence and not args.strict_no_whitelist
    args.outdir.mkdir(parents=True, exist_ok=True)
    try:
        result = run_gate(
            args.prokadiff_gd,
            args.breseq_gd,
            args.jc_tol_bp,
            allow,
            strip_contig_version=args.normalize_contig_version,
        )
    except FileNotFoundError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    json_path = args.outdir / "jc_acceptance_gate.json"
    tsv_path = args.outdir / "jc_acceptance_reconcile.tsv"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    write_tsv(result, tsv_path)

    c = result["counts"]
    print(
        f"JC_ACCEPTANCE_GATE pass={result['pass']} "
        f"matched={c['matched_exact_tol']} "
        f"equiv_repeat_copy={c['equiv_repeat_copy']} "
        f"pd_vs_rejected={c['pd_matches_breseq_rejected']} "
        f"accepted_over={c['accepted_over_unmatched']} "
        f"accepted_under={c['accepted_under_unmatched']} "
        f"breseq_rejected={c['breseq_jc_rejected']}"
    )
    print(f"wrote {json_path}")
    print(f"wrote {tsv_path}")
    if result["pass"]:
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
