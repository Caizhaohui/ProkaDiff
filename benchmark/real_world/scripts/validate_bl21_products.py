#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from validate_bl21_registry import MANIFEST, SPACERS, read_tsv
from validation_context import (
    ComparisonRole,
    IntendedRelation,
    is_intended_related_relationship,
    parse_declaration,
    parse_event_id_list,
    parse_intended_relation,
    relationship_supports_complete,
)


PRODUCT_FILES: Final = (
    "summary.txt", "report.md", "edit_outcomes.tsv", "post_edit_variants.tsv",
    "unintended.tsv", "offtarget_sites.tsv", "mutation_offtarget_links.tsv",
)
CLASS_KEYS: Final = {"structural", "near_homolog", "scattered_snv"}
EVIDENCE_ONLY_TYPES: Final = {"MC", "RA"}
IDENTITY_COLUMNS: Final = {"event_id", "matched_event_ids_dv1", "unexpected_event_ids_dv1"}
CAUSAL_CLAIM: Final = re.compile(
    r"(?:caused by|caused|driven by|attributed to)\s+(?:cas9|cas12a|editing)\b"
    r"|(?:cas9|cas12a)\s+(?:caused|created|produced)\b",
    re.IGNORECASE,
)
EXACT_BREAKPOINT_CLAIM: Final = re.compile(r"\bexact\s+(?:breakpoint|breakpoints|structure|junction)\b", re.IGNORECASE)
UNSUPPORTED_WHEN_UNQUALIFIED: Final = (
    ("exact repeat-copy breakpoint", re.compile(r"exact\s+repeat-copy\s+breakpoint", re.IGNORECASE)),
    ("exact IS1A copy", re.compile(r"exact\s+IS1A\s+copy", re.IGNORECASE)),
    ("editing-induced structural event", re.compile(r"editing-induced\s+structural(?:\s+event)?", re.IGNORECASE)),
    ("Cas9 caused this event", re.compile(r"\bCas9\s+caused\b", re.IGNORECASE)),
    ("verified parent-derived event", re.compile(r"verified\s+parent-derived\s+event", re.IGNORECASE)),
)
EXTERNAL_SOURCE: Final = re.compile(r"\b(?:breseq|curated\s+interpretation|curated\s+truth)\b", re.IGNORECASE)
PROKADIFF_ATTRIBUTION: Final = re.compile(r"\bProkaDiff\b")
NEGATION: Final = re.compile(r"\b(?:not|never|cannot|can't|without|no)\b|\bdoes not\b|\bdo not\b", re.IGNORECASE)


class ProductValidationError(Exception):
    pass


@dataclass(frozen=True)
class IndexedEvent:
    event_id: str
    gd_type: str
    intended_relation: str
    hidden_mob_companion: bool


def require_columns(path: Path, rows: list[dict[str, str]], required: set[str]) -> None:
    fields, _ = read_tsv(path)
    missing = required.difference(fields)
    if missing:
        raise ProductValidationError(f"{path.name} lacks current M4 columns: {', '.join(sorted(missing))}")
    if not rows and path.name not in {"offtarget_sites.tsv", "mutation_offtarget_links.tsv"}:
        raise ProductValidationError(f"{path.name} unexpectedly has no product rows")


def summary_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        values[key] = value
    return values


def one_event_id(value: str, label: str) -> str:
    try:
        items = parse_event_id_list(value)
    except ValueError as error:
        raise ProductValidationError(str(error)) from error
    if len(items) != 1:
        raise ProductValidationError(f"{label} must carry exactly one DV1 EventId")
    return next(iter(items))


def event_id_set(value: str) -> set[str]:
    try:
        return parse_event_id_list(value)
    except ValueError as error:
        raise ProductValidationError(str(error)) from error


def claim_sentences(report: str) -> list[str]:
    return [part.strip() for part in re.split(r"[\n.]+", report) if part.strip()]


def claim_is_negated(sentence: str, start: int) -> bool:
    return NEGATION.search(sentence[max(0, start - 80):start]) is not None


def external_source_only(sentence: str) -> bool:
    return EXTERNAL_SOURCE.search(sentence) is not None and PROKADIFF_ATTRIBUTION.search(sentence) is None


def unsupported_report_claims(report: str, exact_breakpoint_supported: bool) -> str | None:
    for sentence in claim_sentences(report):
        if external_source_only(sentence):
            continue
        for label, pattern in UNSUPPORTED_WHEN_UNQUALIFIED:
            match = pattern.search(sentence)
            if match and not claim_is_negated(sentence, match.start()):
                return label
        causal = CAUSAL_CLAIM.search(sentence)
        if causal and not claim_is_negated(sentence, causal.start()):
            return "nuclease or editing causality"
        exact = EXACT_BREAKPOINT_CLAIM.search(sentence)
        if exact and not exact_breakpoint_supported and not claim_is_negated(sentence, exact.start()):
            return "exact breakpoint"
    return None


def assert_diagnostic_sidecar_isolated(path: Path) -> None:
    """Reject EventId identity in a repeat-ambiguous diagnostic sidecar.

    Rows are not returned and cannot enter the differential-event index.
    Coordinate overlap with an accepted JC is not a failure.
    """
    if not path.is_file() or path.stat().st_size == 0:
        raise ProductValidationError(f"diagnostic sidecar is missing or empty: {path}")
    fields, rows = read_tsv(path)
    identity_columns = sorted(column for column in fields if column.casefold() in IDENTITY_COLUMNS)
    if identity_columns:
        raise ProductValidationError(
            "repeat_ambiguous_junctions.tsv introduces EventId identity columns: " + ", ".join(identity_columns)
        )
    for row in rows:
        for column, value in row.items():
            if value and re.search(r"DV1_[0-9a-f]{32}", value):
                raise ProductValidationError(
                    f"repeat_ambiguous_junctions.tsv cell {column} carries a DV1 EventId"
                )


def authoritative_event_index(
    variants: list[dict[str, str]],
    unintended: list[dict[str, str]],
    links: list[dict[str, str]],
) -> dict[str, IndexedEvent]:
    """Product-level stand-in for AuditResult.event_index.

    Visible differential rows come from the variant and unintended tables.
    A JC association whose EventId is absent from those tables is admitted
    only as a hidden MOB companion when the same products publish a MOB
    differential event. This is the frozen M4 view exception. It is not a
    coordinate join, and diagnostic sidecar rows are not an input.
    """
    index: dict[str, IndexedEvent] = {}
    for row in variants:
        if row["gd_type"] in EVIDENCE_ONLY_TYPES:
            raise ProductValidationError(f"evidence-only {row['gd_type']} record appeared as a differential variant")
        event_id = one_event_id(row["event_id"], f"differential variant {row['variant_id']}")
        try:
            parse_intended_relation(row["intended_relation"])
        except ValueError as error:
            raise ProductValidationError(str(error)) from error
        index[event_id] = IndexedEvent(event_id, row["gd_type"], row["intended_relation"].upper(), False)
    for row in unintended:
        if row["gd_type"] in EVIDENCE_ONLY_TYPES:
            raise ProductValidationError(f"evidence-only {row['gd_type']} record appeared in unintended.tsv")
        event_id = one_event_id(row["event_id"], "unintended.tsv row")
        current = index.get(event_id)
        if current is None:
            index[event_id] = IndexedEvent(event_id, row["gd_type"], IntendedRelation.NONE.value, False)
        elif current.gd_type != row["gd_type"]:
            raise ProductValidationError(f"EventId {event_id} has conflicting public gd_type values")
    mob_published = any(event.gd_type == "MOB" and not event.hidden_mob_companion for event in index.values())
    for row in links:
        event_id = one_event_id(row["event_id"], f"association {row.get('mutation_id', '')}")
        if event_id in index:
            continue
        if row["mutation_type"] == "JC" and mob_published:
            index[event_id] = IndexedEvent(event_id, "JC", IntendedRelation.NONE.value, True)
            continue
        raise ProductValidationError(f"association EventId {event_id} does not resolve to a differential event")
    return index


def require_resolved(event_ids: set[str], index: dict[str, IndexedEvent], label: str) -> None:
    missing = sorted(event_id for event_id in event_ids if event_id not in index)
    if missing:
        raise ProductValidationError(f"{label} EventId does not resolve: {', '.join(missing)}")


def site_search_status(report: str, candidates: list[dict[str, str]]) -> str:
    if "Candidate guide-site search was not performed" in report:
        status = "NOT_PERFORMED"
    elif "Candidate guide-site search was performed and found no candidate sites" in report:
        status = "PERFORMED_NO_CANDIDATES"
    elif candidates or "Candidate guide sites were found" in report:
        status = "PERFORMED_WITH_CANDIDATES"
    else:
        raise ProductValidationError("report.md does not distinguish candidate search status")
    if status == "PERFORMED_WITH_CANDIDATES" and not candidates:
        raise ProductValidationError("report claims guide candidates that are absent from offtarget_sites.tsv")
    if status == "PERFORMED_NO_CANDIDATES" and candidates:
        raise ProductValidationError("report says no guide candidates but offtarget_sites.tsv contains candidates")
    return status


def manifest_role(sample_id: str) -> ComparisonRole:
    _, rows = read_tsv(MANIFEST)
    matches = [row for row in rows if row.get("sample_id") == sample_id]
    if len(matches) != 1:
        raise ProductValidationError(f"manifest must contain exactly one row for {sample_id}")
    try:
        return parse_declaration(matches[0]).role
    except (KeyError, ValueError) as error:
        raise ProductValidationError(f"invalid comparison declaration for {sample_id}: {error}") from error


def validate_complete(row: dict[str, str], index: dict[str, IndexedEvent], matched: set[str], unexpected: set[str]) -> None:
    if not matched:
        raise ProductValidationError("COMPLETE lacks a matched intended EventId")
    if unexpected:
        raise ProductValidationError("COMPLETE retains unexpected intended EventIds")
    require_resolved(matched, index, "COMPLETE matched")
    for event_id in matched:
        relation = index[event_id].intended_relation
        if not relationship_supports_complete(relation):
            raise ProductValidationError(f"COMPLETE EventId {event_id} has incompatible relationship {relation}")
    for event in index.values():
        if event.hidden_mob_companion:
            continue
        if event.intended_relation == IntendedRelation.UNEXPECTED_AT_TARGET.value:
            raise ProductValidationError("COMPLETE contradicts an unexpected-at-target differential event")
        if relationship_supports_complete(event.intended_relation) and event.event_id not in matched:
            raise ProductValidationError(f"COMPLETE omits expected-edit EventId {event.event_id}")


def validate_missing(matched: set[str], unexpected: set[str], index: dict[str, IndexedEvent]) -> None:
    if matched or unexpected:
        raise ProductValidationError("MISSING status must not claim matched or unexpected DV1 events")
    related = sorted(
        event.event_id for event in index.values()
        if not event.hidden_mob_companion and is_intended_related_relationship(event.intended_relation)
    )
    if related:
        raise ProductValidationError(
            "MISSING status conflicts with intended-related differential event(s): " + ", ".join(related)
        )


def validate_products(
    product_dir: Path,
    sample_id: str,
    guide_id: str,
    comparison_role: str = "PEER_COMPARATOR",
    target_review: Path | None = None,
    diagnostic_sidecar: Path | None = None,
) -> dict[str, int | str]:
    if guide_id not in SPACERS:
        raise ProductValidationError(f"unknown guide id {guide_id!r}; both registry guides must be tested")
    try:
        role = ComparisonRole(comparison_role)
    except ValueError as error:
        raise ProductValidationError(f"unknown comparison role {comparison_role!r}") from error
    if role is ComparisonRole.REFERENCE_ONLY:
        raise ProductValidationError("REFERENCE_ONLY is a validation evidence mode, not an M4 differential product mode")
    if manifest_role(sample_id) is not role:
        raise ProductValidationError(f"manifest comparison role does not match requested {role.value}")
    paths = {name: product_dir / name for name in PRODUCT_FILES}
    missing_files = [str(path) for path in paths.values() if not path.is_file()]
    if missing_files:
        raise ProductValidationError(f"missing M4 product file(s): {', '.join(missing_files)}")
    sidecar = diagnostic_sidecar
    if sidecar is None:
        candidate = product_dir / "repeat_ambiguous_junctions.tsv"
        if candidate.is_file():
            sidecar = candidate
    if sidecar is not None:
        assert_diagnostic_sidecar_isolated(sidecar)

    summary = summary_values(paths["summary.txt"])
    report = paths["report.md"].read_text(encoding="utf-8")
    exact_breakpoint_supported = target_review_has_exact_prokadiff_breakpoint(target_review)
    unsupported = unsupported_report_claims(report, exact_breakpoint_supported)
    if unsupported:
        raise ProductValidationError(f"report attributes unsupported claim to ProkaDiff: {unsupported}")
    spacer = SPACERS[guide_id]
    if sample_id not in report or f"| Guide spacer | `{spacer}` |" not in report:
        raise ProductValidationError(f"report sample/guide identity mismatch for {sample_id}/{guide_id}")

    outcome_fields, outcomes = read_tsv(paths["edit_outcomes.tsv"])
    require_columns(paths["edit_outcomes.tsv"], outcomes, {
        "edit_id", "status", "matched_event_ids_dv1", "unexpected_event_ids_dv1",
    })
    if "event_id" in outcome_fields:
        raise ProductValidationError("edit_outcomes.tsv must retain legacy ID cells and append the two DV1 edit columns")

    variant_fields, variants = read_tsv(paths["post_edit_variants.tsv"])
    require_columns(paths["post_edit_variants.tsv"], variants, {"variant_id", "gd_type", "legacy_class", "event_id", "intended_relation"})
    if variant_fields[-1] != "event_id":
        raise ProductValidationError("post_edit_variants.tsv V2 event_id must be appended")

    unintended_fields, unintended = read_tsv(paths["unintended.tsv"])
    require_columns(paths["unintended.tsv"], unintended, {"class", "gd_type", "event_id"})
    if unintended_fields[-1] != "event_id":
        raise ProductValidationError("unintended.tsv V2 event_id must be appended")

    _, sites = read_tsv(paths["offtarget_sites.tsv"])
    require_columns(paths["offtarget_sites.tsv"], sites, {
        "site_id", "seq_id", "start", "end", "pam", "mismatches",
        "bulge_type", "bulge_size", "search_backend",
    })
    site_map = {row["site_id"]: row for row in sites}
    if len(site_map) != len(sites):
        raise ProductValidationError("offtarget_sites.tsv contains duplicate site IDs")
    search_status = site_search_status(report, sites)

    link_fields, links = read_tsv(paths["mutation_offtarget_links.tsv"])
    require_columns(paths["mutation_offtarget_links.tsv"], links, {
        "mutation_id", "site_id", "mutation_type", "mutation_position", "site_start", "site_end",
        "distance_to_site", "mismatches", "pam", "event_id", "junction_side",
        "search_backend", "bulge_type", "bulge_size",
    })
    if link_fields[-5:] != ["event_id", "junction_side", "search_backend", "bulge_type", "bulge_size"]:
        raise ProductValidationError("mutation_offtarget_links.tsv M4 evidence columns are not appended in the approved order")
    if "SITE_UNKNOWN" in {row.get("site_id", "") for row in links}:
        raise ProductValidationError("association output contains fabricated SITE_UNKNOWN identity")

    index = authoritative_event_index(variants, unintended, links)
    variant_by_event_id = {
        event.event_id: row for row, event in (
            (row, index[one_event_id(row["event_id"], row["variant_id"])]) for row in variants
        )
    }
    matched: set[str] = set()
    unexpected: set[str] = set()
    for row in outcomes:
        row_matched = event_id_set(row["matched_event_ids_dv1"])
        row_unexpected = event_id_set(row["unexpected_event_ids_dv1"])
        require_resolved(row_matched | row_unexpected, index, row["edit_id"])
        matched.update(row_matched)
        unexpected.update(row_unexpected)
        status = row["status"].upper()
        if status == "COMPLETE":
            validate_complete(row, index, row_matched, row_unexpected)
        elif status == "MISSING":
            validate_missing(row_matched, row_unexpected, index)
        elif status == "PARTIAL" and not (row_matched or row_unexpected):
            raise ProductValidationError("PARTIAL status has no qualifying DV1 event relationship")
        elif status == "UNEXPECTED_STRUCTURE" and not row_unexpected:
            raise ProductValidationError("UNEXPECTED_STRUCTURE has no authoritative unexpected DV1 EventId")
        elif status not in {"COMPLETE", "MISSING", "PARTIAL", "UNEXPECTED_STRUCTURE"}:
            raise ProductValidationError(f"unknown intended status {status}")

    for row in links:
        site = site_map.get(row["site_id"])
        if site is None:
            raise ProductValidationError(f"association references unknown candidate site {row['site_id']!r}")
        event_id = one_event_id(row["event_id"], f"association {row['mutation_id']}")
        if event_id not in index:
            raise ProductValidationError(f"association EventId {event_id} does not resolve to a differential event")
        for link_field, site_field in (
            ("site_start", "start"), ("site_end", "end"), ("pam", "pam"),
            ("mismatches", "mismatches"), ("search_backend", "search_backend"),
            ("bulge_type", "bulge_type"), ("bulge_size", "bulge_size"),
        ):
            if row.get(link_field, "").casefold() != site.get(site_field, "").casefold():
                raise ProductValidationError(f"association {link_field} does not match its candidate-site evidence")
        if row["mutation_type"] == "JC":
            if row["junction_side"] not in {"side_1", "side_2"}:
                raise ProductValidationError("JC association has invalid junction_side")
        elif row["junction_side"] != "NA":
            raise ProductValidationError("non-JC association must use junction_side=NA")
        distance = row.get("distance_to_site", "")
        if distance not in {"", "NA"}:
            try:
                observed_distance = int(distance)
            except (OverflowError, ValueError) as error:
                raise ProductValidationError("association distance_to_site is not a measured numeric value") from error
            variant = variant_by_event_id.get(event_id)
            if variant is None:
                continue
            mutation_seq = variant["seq_id"]
            try:
                mutation_start = int(variant["position"])
                mutation_end = int(variant["end"])
            except (KeyError, ValueError) as error:
                raise ProductValidationError("associated variant has invalid public coordinates") from error
            if row["mutation_type"] == "JC" and row["junction_side"] == "side_2":
                second_side = re.fullmatch(r"(.+):(\d+)", variant.get("alt", ""))
                if second_side is not None:
                    mutation_seq = second_side.group(1)
                    mutation_start = int(second_side.group(2))
                    mutation_end = mutation_start
            if mutation_seq != site["seq_id"]:
                continue
            site_start = int(site["start"])
            site_end = int(site["end"])
            if mutation_start == mutation_end or row["mutation_type"] == "JC":
                expected_distance = max(site_start - mutation_start, mutation_start - site_end, 0)
            else:
                expected_distance = max(site_start - mutation_end, mutation_start - site_end, 0)
            if observed_distance != expected_distance:
                raise ProductValidationError("association distance_to_site disagrees with EventId-linked coordinates")

    visible_ids = [event.event_id for event in index.values() if not event.hidden_mob_companion]
    omitted = sorted(event_id for event_id in visible_ids if event_id not in report)
    if omitted:
        raise ProductValidationError(f"report.md omits visible differential EventId(s): {', '.join(omitted)}")

    class_counts = {key: 0 for key in CLASS_KEYS}
    for row in unintended:
        if row["class"] not in class_counts:
            raise ProductValidationError(f"unexpected M4 mutation class {row['class']!r}")
        class_counts[row["class"]] += 1
    for key, value in class_counts.items():
        if summary.get(key) != str(value):
            raise ProductValidationError(f"summary/report source mismatch for class {key}: summary={summary.get(key)!r}, TSV={value}")

    total_match = re.search(r"\*\*Total post-edit differential variants:\*\* (\d+)", report)
    if total_match is None or int(total_match.group(1)) != len(variants):
        raise ProductValidationError("report.md post-edit variant count disagrees with post_edit_variants.tsv")
    if any(row["status"].upper() not in report for row in outcomes):
        raise ProductValidationError("report.md does not expose every intended assessment status")
    candidate_event_match = re.search(r"\*\*Candidate guide-dependent off-target events:\*\* (\d+)", report)
    if candidate_event_match is None or int(candidate_event_match.group(1)) != len({row["mutation_id"] for row in links}):
        raise ProductValidationError("report.md candidate-associated event count disagrees with association TSV")
    if summary.get("intended_provided") != "yes":
        raise ProductValidationError("summary.txt does not record the supplied intended table")

    return {
        "variants": len(variants),
        "unintended": len(unintended),
        "intended_edits": len(outcomes),
        "matched_event_ids": len(matched),
        "unexpected_event_ids": len(unexpected),
        "associations": len(links),
        "association_event_ids": len({one_event_id(row["event_id"], "association") for row in links}) if links else 0,
        "candidate_search_status": search_status,
        "candidate_sites": len(sites),
        "candidate_pams": ",".join(sorted({row["pam"] for row in sites})) or "NONE",
        "candidate_mismatches": ",".join(sorted({row["mismatches"] for row in sites})) or "NONE",
        "candidate_bulges": ",".join(sorted({f"{row['bulge_type']}:{row['bulge_size']}" for row in sites})) or "NONE",
        "hidden_mob_companions": sum(event.hidden_mob_companion for event in index.values()),
    }


def target_review_has_exact_prokadiff_breakpoint(path: Path | None) -> bool:
    if path is None:
        return False
    fields, rows = read_tsv(path)
    if not {"evidence_source", "exact_breakpoint_supported"}.issubset(fields):
        raise ProductValidationError("target-locus review lacks breakpoint support fields")
    return any(
        row["evidence_source"] == "PROKADIFF" and row["exact_breakpoint_supported"].lower() == "true"
        for row in rows
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Cross-check one current M4 BL21 product directory.")
    parser.add_argument("--products", type=Path)
    parser.add_argument("--sample")
    parser.add_argument("--guide-id", choices=sorted(SPACERS))
    parser.add_argument("--comparison-role", choices=[role.value for role in ComparisonRole])
    parser.add_argument("--target-review", type=Path)
    parser.add_argument("--diagnostic-sidecar", type=Path)
    args = parser.parse_args()
    try:
        if args.products is None:
            if args.diagnostic_sidecar is None:
                parser.error("--products is required unless --diagnostic-sidecar is set")
            assert_diagnostic_sidecar_isolated(args.diagnostic_sidecar)
            print("diagnostic_sidecar\tisolated")
            return 0
        if not args.sample or not args.guide_id or not args.comparison_role:
            parser.error("--sample, --guide-id, and --comparison-role are required with --products")
        counts = validate_products(
            args.products, args.sample, args.guide_id, args.comparison_role,
            args.target_review, args.diagnostic_sidecar,
        )
    except (ProductValidationError, OSError) as error:
        print(f"BL21 product validation failed: {error}", file=sys.stderr)
        return 1
    for key, value in counts.items():
        print(f"{key}\t{value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
