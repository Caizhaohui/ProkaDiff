import re
from dataclasses import dataclass
from enum import Enum


class ComparisonRole(str, Enum):
    REFERENCE_ONLY = "REFERENCE_ONLY"
    MATCHED_PARENT = "MATCHED_PARENT"
    PEER_COMPARATOR = "PEER_COMPARATOR"


class IntendedRelation(str, Enum):
    """Public projection of the frozen M2 intended-event roles.

    ExpectedConstituent is EXPECTED, or PART_OF_EXPECTED when one declared
    edit has more than one matched constituent. PartialObservation and
    UnexpectedAtLocus are both published as UNEXPECTED_AT_TARGET. NONE is the
    absence of an intended-locus relationship.
    """

    NONE = "NONE"
    EXPECTED = "EXPECTED"
    PART_OF_EXPECTED = "PART_OF_EXPECTED"
    UNEXPECTED_AT_TARGET = "UNEXPECTED_AT_TARGET"


EVENT_ID = re.compile(r"DV1_[0-9a-f]{32}\Z")
EMPTY_EVENT_LIST = {"", "NA", "NONE"}
COMPLETE_SUPPORTING_RELATIONS = {
    IntendedRelation.EXPECTED,
    IntendedRelation.PART_OF_EXPECTED,
}


@dataclass(frozen=True)
class ComparisonDeclaration:
    role: ComparisonRole
    baseline_id: str
    biological_parent_verified: bool

    @property
    def has_differential_product(self) -> bool:
        return self.role is not ComparisonRole.REFERENCE_ONLY


def parse_declaration(row: dict[str, str]) -> ComparisonDeclaration:
    role = ComparisonRole(row["comparison_role"])
    baseline = row.get("starter_id", "").strip()
    if role is ComparisonRole.REFERENCE_ONLY:
        baseline = row.get("reference", "").strip()
    if not baseline:
        raise ValueError(f"{role.value} requires an explicit baseline identifier")
    raw_verified = row.get("biological_parent_verified", "false").strip().lower()
    if raw_verified not in {"true", "false"}:
        raise ValueError("biological_parent_verified must be true or false")
    verified = raw_verified == "true"
    if role is ComparisonRole.MATCHED_PARENT and not verified:
        raise ValueError("MATCHED_PARENT requires verified biological parent metadata")
    if role is not ComparisonRole.MATCHED_PARENT and verified:
        raise ValueError(f"{role.value} cannot declare a verified biological parent")
    return ComparisonDeclaration(role, baseline, verified)


def parse_intended_relation(value: str) -> IntendedRelation:
    try:
        return IntendedRelation(value.strip().upper())
    except ValueError as error:
        raise ValueError(f"unknown intended relationship {value!r}") from error


def is_intended_related_relationship(value: str) -> bool:
    """True for every frozen relationship other than NONE."""
    return parse_intended_relation(value) is not IntendedRelation.NONE


def relationship_supports_complete(value: str) -> bool:
    """Expected-edit evidence compatible with comparator-relative COMPLETE."""
    return parse_intended_relation(value) in COMPLETE_SUPPORTING_RELATIONS


def parse_event_id_list(value: str | None) -> set[str]:
    if value is None or value.strip() in EMPTY_EVENT_LIST:
        return set()
    items = {item.strip() for item in value.split(",") if item.strip() and item.strip() not in EMPTY_EVENT_LIST}
    invalid = sorted(item for item in items if EVENT_ID.fullmatch(item) is None)
    if invalid:
        raise ValueError(f"invalid DV1 EventId value(s): {', '.join(invalid)}")
    return items


def peer_complete_has_differential_support(row: dict[str, str]) -> bool:
    """Structural COMPLETE support from edit-outcome EventId lists alone.

    Relationship compatibility is checked against the variant table by the
    product validator. A single-sample target phenotype is not an input.
    """
    if row.get("status", "").upper() != "COMPLETE":
        return True
    matched = parse_event_id_list(row.get("matched_event_ids_dv1"))
    unexpected = parse_event_id_list(row.get("unexpected_event_ids_dv1"))
    return bool(matched) and not unexpected
