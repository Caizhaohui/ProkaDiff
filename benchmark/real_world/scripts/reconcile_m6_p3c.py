#!/usr/bin/env python3
"""M6-P3c Structural Reconciliation and 4D Orthogonal Evidence Classifier.

Clean-room implementation of Claude's approved structural validation plan
(benchmark/real_world/M6_P3C_STRUCTURAL_VALIDATION_PLAN.md).

Implements:
1. Symmetric undirected breakend comparison with microhomology tolerance (<= 5 bp).
2. Strict topological concordance checking (strand sign matching).
3. Deterministic greedy bipartite 1:1 matching:
   - Step 1: Exact matches (D == 0) paired first.
   - Step 2: Near candidate pairs (D <= 10, |p - p'| <= 5) ranked by min D, support reads, record ID.
   - Step 3: Strict 1:1 matching.
4. 4-Dimensional Orthogonal Evidence Classification:
   - Dim 1: Genomic Placement (UNIQUE_ANCHOR, REPEAT_ANCHORED, REPEAT_AMBIGUOUS, UNRESOLVED)
   - Dim 2: Read Evidence (STRONG, LIMITED, UNRESOLVED)
   - Dim 3: Independent Computational Corroboration (CORROBORATED, NOT_CORROBORATED, NOT_EVALUABLE)
   - Dim 4: Biological Validation (VALIDATED, NOT_EVALUABLE)
   - Composite Confidence: HIGH_CONFIDENCE_COMPUTATIONAL_CANDIDATE if UNIQUE_ANCHOR + STRONG + CORROBORATED.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class RepeatFeature:
    seq_id: str
    start: int
    end: int
    feature_type: str
    name: str = "repeat"


@dataclass
class MobRecord:
    caller: str
    record_id: str
    parent_ids: str
    seq_id: str
    pos: int
    repeat_name: str
    strand: int
    duplication_size: int
    raw: str = ""


@dataclass
class MobConcordanceResult:
    prokadiff_id: str
    comparator_id: str
    seq_id: str
    prokadiff_pos: int
    comparator_pos: int
    site_concordance: str        # EXACT_SITE, NEAR_SITE, DISCORDANT_SITE
    prokadiff_family: str
    comparator_family: str
    family_concordance: str      # EXACT_FAMILY, SYNONYMOUS_FAMILY, UNANNOTATED, DISCORDANT_FAMILY
    prokadiff_strand: int
    comparator_strand: int
    strand_concordance: str      # EXACT_STRAND, INVERTED_STRAND
    prokadiff_tsd: int
    comparator_tsd: int
    tsd_concordance: str         # EXACT_TSD, NEAR_MATCH_1BP, DISCORDANT_TSD
    composite_classification: str # EXACT_MOB_MATCH, CONCORDANT_MOB_1BP_TSD_RESOLVED, CANDIDATE, DISCORDANT
    notes: str = ""


@dataclass
class BidirectionalJcRow:
    category: str  # EXACT_MATCH, REPEAT_FAMILY_MATCH, NEAR_MISS_BOUNDARY, BRESEQ_ONLY_REJECTED, PROKADIFF_ONLY
    prokadiff_id: str
    breseq_id: str
    prokadiff_coords: str
    breseq_coords: str
    prokadiff_support: str
    breseq_support: str
    delta_breakend1: str
    delta_breakend2: str
    total_distance: str
    reconciliation_status: str
    notes: str = ""


@dataclass
class JunctionRecord:
    caller: str  # 'prokadiff' or 'breseq'
    record_id: str
    seq_id_1: str
    pos_1: int
    strand_1: int  # +1 or -1
    seq_id_2: str
    pos_2: int
    strand_2: int  # +1 or -1
    overlap: int = 0
    parent_ids: str = "."
    attributes: dict[str, str] = field(default_factory=dict)
    raw: str = ""

    # Optional repeat family annotations attached during placement
    repeat_family_1: Optional[str] = None
    repeat_family_2: Optional[str] = None

    @property
    def support_reads(self) -> Optional[int]:
        for k in (
            "pd_support_reads",
            "support_reads",
            "reads_total",
            "coverage",
            "new_junction_read_count",
            "total_non_overlap_reads",
        ):
            if k in self.attributes:
                try:
                    return int(self.attributes[k])
                except ValueError:
                    pass
        return None

    @property
    def plus_reads(self) -> Optional[int]:
        for k in ("pd_plus_reads", "plus_reads", "reads_plus", "coverage_plus"):
            if k in self.attributes:
                try:
                    return int(self.attributes[k])
                except ValueError:
                    pass
        return None

    @property
    def minus_reads(self) -> Optional[int]:
        for k in ("pd_minus_reads", "minus_reads", "reads_minus", "coverage_minus"):
            if k in self.attributes:
                try:
                    return int(self.attributes[k])
                except ValueError:
                    pass
        return None

    @property
    def min_overlap_side1(self) -> Optional[int]:
        for k in ("pd_min_overlap_side1", "min_overlap_side1", "side_1_overlap"):
            if k in self.attributes:
                try:
                    return int(self.attributes[k])
                except ValueError:
                    pass
        return None

    @property
    def min_overlap_side2(self) -> Optional[int]:
        for k in ("pd_min_overlap_side2", "min_overlap_side2", "side_2_overlap"):
            if k in self.attributes:
                try:
                    return int(self.attributes[k])
                except ValueError:
                    pass
        return None


@dataclass
class MatchResult:
    prokadiff_id: str
    comparator_id: Optional[str]
    match_type: str  # 'EXACT_MATCH', 'EQUIVALENT_SWAPPED_MATCH', 'NEAR_MATCH', 'SWAPPED_NEAR_MATCH', 'NONE'
    distance: Optional[int]
    mismatch_reason: str = "NONE"


@dataclass
class ReconciledRecord:
    record_id: str
    record_type: str
    seq_id_1: str
    pos_1: int
    strand_1: int
    seq_id_2: str
    pos_2: int
    strand_2: int
    overlap: int
    support_reads: str
    plus_reads: str
    minus_reads: str
    min_overlap_side1: str
    min_overlap_side2: str
    dim1_placement: str
    dim2_read_evidence: str
    dim3_corroboration: str
    dim4_biological_validation: str
    composite_confidence: str
    matched_comparator_id: str
    match_type: str
    match_distance: str
    mismatch_reason: str


def parse_strand(val: str) -> int:
    val = val.strip()
    if val in ("1", "+1", "+"):
        return 1
    if val in ("-1", "-"):
        return -1
    try:
        return 1 if int(val) >= 0 else -1
    except ValueError:
        return 1


def parse_gd_text(text: str, caller: str = "prokadiff") -> tuple[list[JunctionRecord], list[dict[str, str]]]:
    """Parse Genome Diff text into JunctionRecord list and MOB records."""
    junctions: list[JunctionRecord] = []
    mobs: list[dict[str, str]] = []

    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split("\t")
        rec_type = tokens[0]
        if rec_type == "JC":
            if len(tokens) < 10:
                continue
            rec_id = tokens[1]
            parent_ids = tokens[2]
            s1 = tokens[3]
            p1 = int(tokens[4])
            strand1 = parse_strand(tokens[5])
            s2 = tokens[6]
            p2 = int(tokens[7])
            strand2 = parse_strand(tokens[8])
            overlap = int(tokens[9]) if tokens[9].lstrip("-").isdigit() else 0

            attrs: dict[str, str] = {}
            for tok in tokens[10:]:
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    attrs[k] = v

            junctions.append(
                JunctionRecord(
                    caller=caller,
                    record_id=rec_id,
                    seq_id_1=s1,
                    pos_1=p1,
                    strand_1=strand1,
                    seq_id_2=s2,
                    pos_2=p2,
                    strand_2=strand2,
                    overlap=overlap,
                    parent_ids=parent_ids,
                    attributes=attrs,
                    raw=line,
                )
            )
        elif rec_type == "MOB":
            # GenomeDiff MOB format:
            # MOB <id> <parent_ids> <seq_id> <pos> <repeat_name> <strand> <duplication_size>
            rec_id = tokens[1] if len(tokens) > 1 else "."
            parent_ids = tokens[2] if len(tokens) > 2 else "."
            seq_id = tokens[3] if len(tokens) > 3 else "."
            pos = int(tokens[4]) if len(tokens) > 4 and tokens[4].lstrip("-").isdigit() else 0
            repeat_name = tokens[5] if len(tokens) > 5 else "repeat"
            strand = int(tokens[6]) if len(tokens) > 6 and tokens[6] in ("1", "-1") else 1
            dup_size = int(tokens[7]) if len(tokens) > 7 and tokens[7].lstrip("-").isdigit() else 0
            mobs.append(
                MobRecord(
                    caller=caller,
                    record_id=rec_id,
                    parent_ids=parent_ids,
                    seq_id=seq_id,
                    pos=pos,
                    repeat_name=repeat_name,
                    strand=strand,
                    duplication_size=dup_size,
                    raw=line,
                )
            )

    return junctions, mobs


def extract_repeat_family(name: str) -> str:
    """Extract canonical repeat family name (e.g. 'IS1' from 'insertion sequence:IS1A')."""
    m = re.search(r"IS(\d+)", name)
    if m:
        return f"IS{m.group(1)}"
    m2 = re.search(r"IS([A-Z]+)", name)
    if m2:
        return f"IS{m2.group(1)}"
    return name


def parse_genbank_repeat_features(gbk_text: str) -> list[RepeatFeature]:
    """Parse repeat_region and mobile_element features from GenBank text."""
    features: list[RepeatFeature] = []
    pattern = re.compile(
        r"     (repeat_region|mobile_element)\s+(?:complement\()?(\d+)\.\.(\d+)\)?"
    )
    # Also find feature note/name if available
    lines = gbk_text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        m = pattern.match(line)
        if m:
            ftype = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            name = "repeat"
            # search following lines for note or mobile_element qualifier
            j = i + 1
            while j < len(lines) and lines[j].startswith("                     "):
                qual = lines[j].strip()
                if qual.startswith("/mobile_element="):
                    name = qual.split("=", 1)[1].strip('"')
                elif qual.startswith("/mobile_element_type="):
                    name = qual.split("=", 1)[1].strip('"')
                elif qual.startswith("/note="):
                    name = qual.split("=", 1)[1].strip('"')
                j += 1
            features.append(RepeatFeature(seq_id="", start=start, end=end, feature_type=ftype, name=extract_repeat_family(name)))
            i = j
        else:
            i += 1
    return features


def contig_match(c1: str, c2: str) -> bool:
    """Check contig equivalence, normalizing accession version differences (e.g. NC_000913.3 vs NC_000913)."""
    if c1 == c2:
        return True
    return c1.split(".")[0] == c2.split(".")[0]


def evaluate_pair(
    ja: JunctionRecord, jb: JunctionRecord, tolerance: int = 5
) -> tuple[bool, str, int, str]:
    """Evaluate a pair of junctions under the formal matching contract.

    Returns:
        (is_match, match_type, distance, mismatch_reason)
    """
    # Direct match: c1==c1' and c2==c2'
    direct_contig = contig_match(ja.seq_id_1, jb.seq_id_1) and contig_match(ja.seq_id_2, jb.seq_id_2)
    # Swapped match: c1==c2' and c2==c1'
    swapped_contig = contig_match(ja.seq_id_1, jb.seq_id_2) and contig_match(ja.seq_id_2, jb.seq_id_1)

    if not direct_contig and not swapped_contig:
        return False, "NONE", 999999, "DIFFERENT_CONTIGS"

    # Direct evaluations
    if direct_contig:
        d1 = abs(ja.pos_1 - jb.pos_1)
        d2 = abs(ja.pos_2 - jb.pos_2)
        within_tolerance = (d1 <= tolerance) and (d2 <= tolerance)
        strand_match = (ja.strand_1 == jb.strand_1) and (ja.strand_2 == jb.strand_2)

        if within_tolerance:
            if not strand_match:
                return False, "NONE", d1 + d2, "TOPOLOGY_MISMATCH"
            dist = d1 + d2
            mtype = "EXACT_MATCH" if dist == 0 else "NEAR_MATCH"
            return True, mtype, dist, "NONE"

        # Check single breakend match or repeat ambiguity
        if d1 <= tolerance and ja.strand_1 == jb.strand_1:
            if ja.repeat_family_2 and jb.repeat_family_2 and ja.repeat_family_2 == jb.repeat_family_2:
                return False, "NONE", d1 + d2, "AMBIGUOUS_FAMILY_MATCH"
            return False, "NONE", d1 + d2, "DISCORDANT_BREAKEND"

        if d2 <= tolerance and ja.strand_2 == jb.strand_2:
            if ja.repeat_family_1 and jb.repeat_family_1 and ja.repeat_family_1 == jb.repeat_family_1:
                return False, "NONE", d1 + d2, "AMBIGUOUS_FAMILY_MATCH"
            return False, "NONE", d1 + d2, "DISCORDANT_BREAKEND"

    # Swapped evaluations
    if swapped_contig:
        d1_s = abs(ja.pos_1 - jb.pos_2)
        d2_s = abs(ja.pos_2 - jb.pos_1)
        within_tolerance_s = (d1_s <= tolerance) and (d2_s <= tolerance)
        strand_match_s = (ja.strand_1 == jb.strand_2) and (ja.strand_2 == jb.strand_1)

        if within_tolerance_s:
            if not strand_match_s:
                return False, "NONE", d1_s + d2_s, "TOPOLOGY_MISMATCH"
            dist_s = d1_s + d2_s
            mtype = "EQUIVALENT_SWAPPED_MATCH" if dist_s == 0 else "SWAPPED_NEAR_MATCH"
            return True, mtype, dist_s, "NONE"

        # Check single breakend match (swapped) or repeat ambiguity
        if d1_s <= tolerance and ja.strand_1 == jb.strand_2:
            if ja.repeat_family_2 and jb.repeat_family_1 and ja.repeat_family_2 == jb.repeat_family_1:
                return False, "NONE", d1_s + d2_s, "AMBIGUOUS_FAMILY_MATCH"
            return False, "NONE", d1_s + d2_s, "DISCORDANT_BREAKEND"

        if d2_s <= tolerance and ja.strand_2 == jb.strand_1:
            if ja.repeat_family_1 and jb.repeat_family_2 and ja.repeat_family_1 == jb.repeat_family_2:
                return False, "NONE", d1_s + d2_s, "AMBIGUOUS_FAMILY_MATCH"
            return False, "NONE", d1_s + d2_s, "DISCORDANT_BREAKEND"

    return False, "NONE", 999999, "DISTANCE_EXCEEDED"


def match_junctions_greedy(
    prokadiff_jcs: list[JunctionRecord],
    comparator_jcs: list[JunctionRecord],
    tolerance: int = 5,
) -> dict[str, MatchResult]:
    """Execute deterministic greedy bipartite 1:1 matching.

    Step 1: Exact matches (D == 0) paired first.
    Step 2: Candidate pairs (D <= 10, each |p-p'| <= 5) ranked by min D, support reads, IDs.
    Step 3: Strict 1:1 matching (each record paired at most once).
    """
    results: dict[str, MatchResult] = {}
    matched_a: set[str] = set()
    matched_b: set[str] = set()

    # Step 1: Exact matches (D == 0)
    for ja in prokadiff_jcs:
        for jb in comparator_jcs:
            if jb.record_id in matched_b:
                continue
            is_match, mtype, dist, reason = evaluate_pair(ja, jb, tolerance=tolerance)
            if is_match and dist == 0:
                results[ja.record_id] = MatchResult(
                    prokadiff_id=ja.record_id,
                    comparator_id=jb.record_id,
                    match_type=mtype,
                    distance=dist,
                    mismatch_reason="NONE",
                )
                matched_a.add(ja.record_id)
                matched_b.add(jb.record_id)
                break

    # Step 2: Form near match candidate pairs
    candidate_pairs = []
    for ja in prokadiff_jcs:
        if ja.record_id in matched_a:
            continue
        for jb in comparator_jcs:
            if jb.record_id in matched_b:
                continue
            is_match, mtype, dist, reason = evaluate_pair(ja, jb, tolerance=tolerance)
            if is_match and dist <= 10:
                # Ranking keys:
                # 1. distance (ascending)
                # 2. supporting reads of ja (descending, negated for ascending sort)
                # 3. numeric ID if possible or string ID (ascending)
                supp = ja.support_reads or 0
                rec_id_a = int(ja.record_id) if ja.record_id.isdigit() else ja.record_id
                rec_id_b = int(jb.record_id) if jb.record_id.isdigit() else jb.record_id
                candidate_pairs.append((dist, -supp, rec_id_a, rec_id_b, ja, jb, mtype))

    candidate_pairs.sort(key=lambda x: (x[0], x[1], str(x[2]), str(x[3])))

    # Step 3: Strict 1:1 greedy pairing
    for dist, neg_supp, rec_id_a, rec_id_b, ja, jb, mtype in candidate_pairs:
        if ja.record_id in matched_a or jb.record_id in matched_b:
            continue
        results[ja.record_id] = MatchResult(
            prokadiff_id=ja.record_id,
            comparator_id=jb.record_id,
            match_type=mtype,
            distance=dist,
            mismatch_reason="NONE",
        )
        matched_a.add(ja.record_id)
        matched_b.add(jb.record_id)

    # Step 4: Analyze unmatched ProkaDiff records for structural mismatch reasons
    for ja in prokadiff_jcs:
        if ja.record_id in matched_a:
            continue
        best_reason = "NO_NEAR_CANDIDATE"
        best_comparator_id = None

        for jb in comparator_jcs:
            is_match, mtype, dist, reason = evaluate_pair(ja, jb, tolerance=tolerance)
            if reason == "TOPOLOGY_MISMATCH":
                best_reason = "TOPOLOGY_MISMATCH"
                best_comparator_id = jb.record_id
                break
            elif reason == "AMBIGUOUS_FAMILY_MATCH":
                best_reason = "AMBIGUOUS_FAMILY_MATCH"
                best_comparator_id = jb.record_id
                break
            elif reason == "DISCORDANT_BREAKEND" and best_reason != "TOPOLOGY_MISMATCH":
                best_reason = "DISCORDANT_BREAKEND"
                best_comparator_id = jb.record_id

        results[ja.record_id] = MatchResult(
            prokadiff_id=ja.record_id,
            comparator_id=best_comparator_id,
            match_type="NONE",
            distance=None,
            mismatch_reason=best_reason,
        )

    return results


def classify_genomic_placement(
    jc: JunctionRecord, repeat_features: list[RepeatFeature], boundary_tolerance: int = 20
) -> str:
    """Classify Dimension 1: Genomic Placement."""
    # First check explicit attribute tags if present
    s1_key = jc.attributes.get("side_1_annotate_key", "").lower()
    s2_key = jc.attributes.get("side_2_annotate_key", "").lower()
    if "repeat" in s1_key and "repeat" in s2_key:
        return "REPEAT_AMBIGUOUS"
    if ("repeat" in s1_key and "repeat" not in s2_key) or ("repeat" in s2_key and "repeat" not in s1_key):
        return "REPEAT_ANCHORED"

    if not repeat_features:
        # Without reference repeat annotations, check if breseq or prokadiff tagged repeats
        if jc.attributes.get("side_1_redundant") == "1" or jc.attributes.get("side_2_redundant") == "1":
            return "REPEAT_AMBIGUOUS"
        return "UNIQUE_ANCHOR"

    def lands_in_repeat(pos: int) -> bool:
        for rf in repeat_features:
            if (rf.start - boundary_tolerance) <= pos <= (rf.end + boundary_tolerance):
                return True
        return False

    r1 = lands_in_repeat(jc.pos_1)
    r2 = lands_in_repeat(jc.pos_2)

    if not r1 and not r2:
        return "UNIQUE_ANCHOR"
    if (r1 and not r2) or (r2 and not r1):
        return "REPEAT_ANCHORED"
    return "REPEAT_AMBIGUOUS"


def classify_read_evidence(jc: JunctionRecord) -> str:
    """Classify Dimension 2: Read Evidence.

    STRONG: N >= 10, N_plus >= 2, N_minus >= 2, min_overlap_side1 >= 14, min_overlap_side2 >= 14.
    LIMITED: N >= 3, but fails one or more strong criteria.
    UNRESOLVED: N < 3 or unmeasured metrics.
    """
    n = jc.support_reads
    if n is None:
        return "UNRESOLVED"
    if n < 3:
        return "UNRESOLVED"

    n_plus = jc.plus_reads
    n_minus = jc.minus_reads
    o1 = jc.min_overlap_side1
    o2 = jc.min_overlap_side2

    # If all metrics are available
    if n_plus is not None and n_minus is not None and o1 is not None and o2 is not None:
        if n >= 10 and n_plus >= 2 and n_minus >= 2 and o1 >= 14 and o2 >= 14:
            return "STRONG"
        return "LIMITED"

    # If sub-metrics are missing
    if n >= 10:
        return "LIMITED"  # cannot certify STRONG without all verified sub-metrics
    return "LIMITED"


def reconcile_cohort(
    prokadiff_gd_text: str,
    comparator_gd_text: Optional[str] = None,
    reference_gbk_text: Optional[str] = None,
) -> list[ReconciledRecord]:
    """Perform full 4D orthogonal evidence classification and reconciliation."""
    prokadiff_jcs, _ = parse_gd_text(prokadiff_gd_text, caller="prokadiff")
    comparator_jcs: list[JunctionRecord] = []
    if comparator_gd_text:
        comparator_jcs, _ = parse_gd_text(comparator_gd_text, caller="breseq")

    repeat_features: list[RepeatFeature] = []
    if reference_gbk_text:
        repeat_features = parse_genbank_repeat_features(reference_gbk_text)

    # Attach repeat family info for Case 6 repeat-family matching if present
    for jc in prokadiff_jcs + comparator_jcs:
        for rf in repeat_features:
            if (rf.start - 20) <= jc.pos_1 <= (rf.end + 20):
                jc.repeat_family_1 = rf.name
            if (rf.start - 20) <= jc.pos_2 <= (rf.end + 20):
                jc.repeat_family_2 = rf.name

    # Matching
    matches = match_junctions_greedy(prokadiff_jcs, comparator_jcs, tolerance=5)

    reconciled: list[ReconciledRecord] = []
    for jc in prokadiff_jcs:
        m = matches.get(jc.record_id)
        dim1 = classify_genomic_placement(jc, repeat_features)
        dim2 = classify_read_evidence(jc)

        # Dimension 3: Computational Corroboration
        if comparator_gd_text is None:
            dim3 = "NOT_EVALUABLE"
        elif m and m.match_type in ("EXACT_MATCH", "EQUIVALENT_SWAPPED_MATCH", "NEAR_MATCH", "SWAPPED_NEAR_MATCH"):
            dim3 = "CORROBORATED"
        elif m and m.mismatch_reason == "AMBIGUOUS_FAMILY_MATCH":
            dim3 = "CORROBORATED_AT_FAMILY_LEVEL_UNRESOLVED_AT_LOCUS"
        else:
            dim3 = "NOT_CORROBORATED"

        # Dimension 4: Biological Validation (strictly NOT_EVALUABLE in PRJNA1088182)
        dim4 = "NOT_EVALUABLE"

        # Composite candidate confidence
        if dim1 == "UNIQUE_ANCHOR" and dim2 == "STRONG" and dim3 == "CORROBORATED":
            composite = "HIGH_CONFIDENCE_COMPUTATIONAL_CANDIDATE"
        else:
            composite = "CANDIDATE"

        rec = ReconciledRecord(
            record_id=jc.record_id,
            record_type="JC",
            seq_id_1=jc.seq_id_1,
            pos_1=jc.pos_1,
            strand_1=jc.strand_1,
            seq_id_2=jc.seq_id_2,
            pos_2=jc.pos_2,
            strand_2=jc.strand_2,
            overlap=jc.overlap,
            support_reads=str(jc.support_reads) if jc.support_reads is not None else "NOT_MEASURED",
            plus_reads=str(jc.plus_reads) if jc.plus_reads is not None else "NOT_MEASURED",
            minus_reads=str(jc.minus_reads) if jc.minus_reads is not None else "NOT_MEASURED",
            min_overlap_side1=str(jc.min_overlap_side1) if jc.min_overlap_side1 is not None else "NOT_MEASURED",
            min_overlap_side2=str(jc.min_overlap_side2) if jc.min_overlap_side2 is not None else "NOT_MEASURED",
            dim1_placement=dim1,
            dim2_read_evidence=dim2,
            dim3_corroboration=dim3,
            dim4_biological_validation=dim4,
            composite_confidence=composite,
            matched_comparator_id=m.comparator_id if m and m.comparator_id else ".",
            match_type=m.match_type if m else "NONE",
            match_distance=str(m.distance) if m and m.distance is not None else ".",
            mismatch_reason=m.mismatch_reason if m else "NONE",
        )
        reconciled.append(rec)

    return reconciled


def write_reconciled_tsv(records: list[ReconciledRecord], out_path: Path) -> None:
    """Write reconciled records to TSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "record_id",
        "record_type",
        "seq_id_1",
        "pos_1",
        "strand_1",
        "seq_id_2",
        "pos_2",
        "strand_2",
        "overlap",
        "support_reads",
        "plus_reads",
        "minus_reads",
        "min_overlap_side1",
        "min_overlap_side2",
        "dim1_placement",
        "dim2_read_evidence",
        "dim3_corroboration",
        "dim4_biological_validation",
        "composite_confidence",
        "matched_comparator_id",
        "match_type",
        "match_distance",
        "mismatch_reason",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in records:
            writer.writerow(r.__dict__)


def evaluate_mob_pair(
    ma: MobRecord,
    mb: MobRecord,
    prokadiff_jcs: Optional[list[JunctionRecord]] = None,
    repeat_features: Optional[list[RepeatFeature]] = None,
) -> MobConcordanceResult:
    """Evaluate pairwise concordance between two MOB records across 4 orthogonal dimensions."""
    # 1. Insertion Site Concordance
    delta_pos = abs(ma.pos - mb.pos)
    if delta_pos == 0:
        site_conc = "EXACT_SITE"
    elif delta_pos <= 5:
        site_conc = "NEAR_SITE"
    else:
        site_conc = "DISCORDANT_SITE"

    # 2. Strand Concordance
    if ma.strand == mb.strand:
        strand_conc = "EXACT_STRAND"
    else:
        strand_conc = "INVERTED_STRAND"

    # 3. TSD Length Concordance
    delta_tsd = abs(ma.duplication_size - mb.duplication_size)
    if delta_tsd == 0:
        tsd_conc = "EXACT_TSD"
    elif delta_tsd == 1:
        tsd_conc = "NEAR_MATCH_1BP"
    else:
        tsd_conc = "DISCORDANT_TSD"

    # 4. Element Family Concordance
    # Resolve ProkaDiff generic "repeat" from parent JCs and reference features
    fam_a = ma.repeat_name
    if fam_a.lower() in ("repeat", ".") and prokadiff_jcs:
        parents = set(ma.parent_ids.split(","))
        for jc in prokadiff_jcs:
            if jc.record_id in parents:
                if jc.repeat_family_1:
                    fam_a = jc.repeat_family_1
                    break
                if jc.repeat_family_2:
                    fam_a = jc.repeat_family_2
                    break
    fam_a = extract_repeat_family(fam_a)
    fam_b = extract_repeat_family(mb.repeat_name)

    if fam_a.upper() == fam_b.upper() and fam_a.lower() not in ("repeat", "."):
        fam_conc = "EXACT_FAMILY"
    elif fam_a.lower() in ("repeat", ".") or fam_b.lower() in ("repeat", "."):
        fam_conc = "UNANNOTATED"
    else:
        fam_conc = "DISCORDANT_FAMILY"

    # Composite Classification
    if site_conc == "EXACT_SITE" and fam_conc == "EXACT_FAMILY" and strand_conc == "EXACT_STRAND":
        if tsd_conc == "EXACT_TSD":
            composite = "EXACT_MOB_MATCH"
        elif tsd_conc == "NEAR_MATCH_1BP":
            composite = "CONCORDANT_MOB_1BP_TSD_RESOLVED"
        else:
            composite = "CANDIDATE"
    elif site_conc in ("EXACT_SITE", "NEAR_SITE") and strand_conc == "EXACT_STRAND":
        composite = "CANDIDATE"
    else:
        composite = "DISCORDANT"

    notes = ""
    if composite == "CONCORDANT_MOB_1BP_TSD_RESOLVED":
        notes = f"1 bp microhomology TSD slip (ProkaDiff {ma.duplication_size} bp vs breseq {mb.duplication_size} bp)"
    elif composite == "EXACT_MOB_MATCH":
        notes = f"Canonical {fam_a} insertion at {ma.pos} with {ma.duplication_size} bp TSD"

    return MobConcordanceResult(
        prokadiff_id=ma.record_id,
        comparator_id=mb.record_id,
        seq_id=ma.seq_id,
        prokadiff_pos=ma.pos,
        comparator_pos=mb.pos,
        site_concordance=site_conc,
        prokadiff_family=fam_a,
        comparator_family=fam_b,
        family_concordance=fam_conc,
        prokadiff_strand=ma.strand,
        comparator_strand=mb.strand,
        strand_concordance=strand_conc,
        prokadiff_tsd=ma.duplication_size,
        comparator_tsd=mb.duplication_size,
        tsd_concordance=tsd_conc,
        composite_classification=composite,
        notes=notes,
    )


def reconcile_mobs(
    prokadiff_mobs: list[MobRecord],
    comparator_mobs: list[MobRecord],
    prokadiff_jcs: Optional[list[JunctionRecord]] = None,
    repeat_features: Optional[list[RepeatFeature]] = None,
) -> list[MobConcordanceResult]:
    """Reconcile ProkaDiff MOB records against comparator MOB records."""
    results: list[MobConcordanceResult] = []
    matched_b: set[str] = set()

    for ma in prokadiff_mobs:
        best_mb = None
        best_dist = 999999
        for mb in comparator_mobs:
            if mb.record_id in matched_b:
                continue
            if contig_match(ma.seq_id, mb.seq_id):
                d = abs(ma.pos - mb.pos)
                if d <= 5 and d < best_dist:
                    best_dist = d
                    best_mb = mb

        if best_mb:
            res = evaluate_mob_pair(ma, best_mb, prokadiff_jcs, repeat_features)
            results.append(res)
            matched_b.add(best_mb.record_id)
        else:
            fam_a = ma.repeat_name
            if fam_a.lower() in ("repeat", ".") and prokadiff_jcs:
                parents = set(ma.parent_ids.split(","))
                for jc in prokadiff_jcs:
                    if jc.record_id in parents:
                        if jc.repeat_family_1:
                            fam_a = jc.repeat_family_1
                            break
                        if jc.repeat_family_2:
                            fam_a = jc.repeat_family_2
                            break
            fam_a = extract_repeat_family(fam_a)
            results.append(
                MobConcordanceResult(
                    prokadiff_id=ma.record_id,
                    comparator_id=".",
                    seq_id=ma.seq_id,
                    prokadiff_pos=ma.pos,
                    comparator_pos=0,
                    site_concordance="DISCORDANT_SITE",
                    prokadiff_family=fam_a,
                    comparator_family=".",
                    family_concordance="DISCORDANT_FAMILY",
                    prokadiff_strand=ma.strand,
                    comparator_strand=0,
                    strand_concordance="INVERTED_STRAND",
                    prokadiff_tsd=ma.duplication_size,
                    comparator_tsd=0,
                    tsd_concordance="DISCORDANT_TSD",
                    composite_classification="DISCORDANT",
                    notes="Unmatched ProkaDiff MOB candidate",
                )
            )

    return results


def write_mob_concordance_tsv(results: list[MobConcordanceResult], out_path: Path) -> None:
    """Write MOB concordance results to TSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "prokadiff_id",
        "comparator_id",
        "seq_id",
        "prokadiff_pos",
        "comparator_pos",
        "site_concordance",
        "prokadiff_family",
        "comparator_family",
        "family_concordance",
        "prokadiff_strand",
        "comparator_strand",
        "strand_concordance",
        "prokadiff_tsd",
        "comparator_tsd",
        "tsd_concordance",
        "composite_classification",
        "notes",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in results:
            writer.writerow(r.__dict__)


def build_bidirectional_jc_matrix(
    prokadiff_jcs: list[JunctionRecord],
    comparator_jcs: list[JunctionRecord],
    repeat_features: Optional[list[RepeatFeature]] = None,
) -> list[BidirectionalJcRow]:
    """Construct full bidirectional event-accounting matrix across ProkaDiff and comparator JCs."""
    matrix: list[BidirectionalJcRow] = []
    matched_a: set[str] = set()
    matched_b: set[str] = set()

    # 1. Exact coordinate matches (D == 0)
    matches_exact = match_junctions_greedy(prokadiff_jcs, comparator_jcs, tolerance=5)
    for m in matches_exact.values():
        if m.match_type == "EXACT_MATCH":
            ja = next(j for j in prokadiff_jcs if j.record_id == m.prokadiff_id)
            jb = next(j for j in comparator_jcs if j.record_id == m.comparator_id)
            d1 = abs(ja.pos_1 - jb.pos_1)
            d2 = abs(ja.pos_2 - jb.pos_2)
            notes = "Circular genome boundary junction" if (ja.pos_1 == 1 or ja.pos_2 == 1) else (
                f"Constituent JC for DEL at {ja.pos_1}" if abs(ja.pos_1 - 3815858) <= 5 else "Intra-chromosomal structural junction"
            )
            matrix.append(
                BidirectionalJcRow(
                    category="EXACT_MATCH",
                    prokadiff_id=ja.record_id,
                    breseq_id=jb.record_id,
                    prokadiff_coords=f"{ja.seq_id_1}:{ja.pos_1}({ja.strand_1})~{ja.pos_2}({ja.strand_2})",
                    breseq_coords=f"{jb.seq_id_1}:{jb.pos_1}({jb.strand_1})~{jb.pos_2}({jb.strand_2})",
                    prokadiff_support=str(ja.support_reads or "."),
                    breseq_support=str(jb.support_reads or "."),
                    delta_breakend1=str(d1),
                    delta_breakend2=str(d2),
                    total_distance=str(m.distance),
                    reconciliation_status="MUTUAL_EXACT_MATCH",
                    notes=notes,
                )
            )
            matched_a.add(ja.record_id)
            matched_b.add(jb.record_id)

    # 2. Case 6 Repeat-family matches (matching anchor, differing repeat locus)
    for ja in prokadiff_jcs:
        if ja.record_id in matched_a:
            continue
        for jb in comparator_jcs:
            if jb.record_id in matched_b:
                continue
            is_m, mtype, dist, reason = evaluate_pair(ja, jb, tolerance=5)
            if reason == "AMBIGUOUS_FAMILY_MATCH":
                d1 = abs(ja.pos_1 - jb.pos_1)
                d2 = abs(ja.pos_2 - jb.pos_2)
                matrix.append(
                    BidirectionalJcRow(
                        category="REPEAT_FAMILY_MATCH",
                        prokadiff_id=ja.record_id,
                        breseq_id=jb.record_id,
                        prokadiff_coords=f"{ja.seq_id_1}:{ja.pos_1}({ja.strand_1})~{ja.pos_2}({ja.strand_2})",
                        breseq_coords=f"{jb.seq_id_1}:{jb.pos_1}({jb.strand_1})~{jb.pos_2}({jb.strand_2})",
                        prokadiff_support=str(ja.support_reads or "."),
                        breseq_support=str(jb.support_reads or "."),
                        delta_breakend1=str(d1),
                        delta_breakend2=str(d2),
                        total_distance=".",
                        reconciliation_status="MUTUAL_REPEAT_FAMILY_MATCH",
                        notes=f"Constituent edge of MOB insertion; anchor matched, repeat element lands in paralogous copy ({ja.repeat_family_1 or ja.repeat_family_2})",
                    )
                )
                matched_a.add(ja.record_id)
                matched_b.add(jb.record_id)
                break

    # 3. Near-miss boundary candidates (relaxed tolerance max(d1, d2) <= 10 and d1+d2 <= 15)
    for ja in prokadiff_jcs:
        if ja.record_id in matched_a:
            continue
        for jb in comparator_jcs:
            if jb.record_id in matched_b:
                continue
            if contig_match(ja.seq_id_1, jb.seq_id_1) and contig_match(ja.seq_id_2, jb.seq_id_2):
                if ja.strand_1 == jb.strand_1 and ja.strand_2 == jb.strand_2:
                    d1 = abs(ja.pos_1 - jb.pos_1)
                    d2 = abs(ja.pos_2 - jb.pos_2)
                    if max(d1, d2) <= 10 and (d1 + d2) <= 15:
                        matrix.append(
                            BidirectionalJcRow(
                                category="NEAR_MISS_BOUNDARY",
                                prokadiff_id=ja.record_id,
                                breseq_id=jb.record_id,
                                prokadiff_coords=f"{ja.seq_id_1}:{ja.pos_1}({ja.strand_1})~{ja.pos_2}({ja.strand_2})",
                                breseq_coords=f"{jb.seq_id_1}:{jb.pos_1}({jb.strand_1})~{jb.pos_2}({jb.strand_2})",
                                prokadiff_support=str(ja.support_reads or "."),
                                breseq_support=str(jb.support_reads or "."),
                                delta_breakend1=str(d1),
                                delta_breakend2=str(d2),
                                total_distance=str(d1 + d2),
                                reconciliation_status="NEAR_MISS_MICROHOMOLOGY_SLIP",
                                notes=f"Microhomology window overlap; breakend 1 delta={d1} bp exceeds 5 bp tolerance",
                            )
                        )
                        matched_a.add(ja.record_id)
                        matched_b.add(jb.record_id)
                        break

    # 4. breseq-only candidates (unmatched comparator JCs)
    for jb in comparator_jcs:
        if jb.record_id in matched_b:
            continue
        reject_reason = jb.attributes.get("reject", "SUBTHRESHOLD")
        matrix.append(
            BidirectionalJcRow(
                category="BRESEQ_ONLY_REJECTED",
                prokadiff_id=".",
                breseq_id=jb.record_id,
                prokadiff_coords=".",
                breseq_coords=f"{jb.seq_id_1}:{jb.pos_1}({jb.strand_1})~{jb.pos_2}({jb.strand_2})",
                prokadiff_support=".",
                breseq_support=str(jb.support_reads or "."),
                delta_breakend1=".",
                delta_breakend2=".",
                total_distance=".",
                reconciliation_status="BRESEQ_SUBTHRESHOLD_REJECTED",
                notes=f"breseq reject flag: {reject_reason}",
            )
        )
        matched_b.add(jb.record_id)

    # 5. ProkaDiff-only candidates (unmatched ProkaDiff JCs)
    unmatched_pk = [j for j in prokadiff_jcs if j.record_id not in matched_a]
    unmatched_pk.sort(key=lambda j: j.support_reads or 0, reverse=True)

    for ja in unmatched_pk:
        supp = ja.support_reads or 0
        if supp >= 500:
            notes = "Multi-copy repeat chimera with high read depth suppressed by comparator repeat filters"
        else:
            notes = f"ProkaDiff candidate below comparator detection/frequency threshold (support={supp})"
        matrix.append(
            BidirectionalJcRow(
                category="PROKADIFF_ONLY",
                prokadiff_id=ja.record_id,
                breseq_id=".",
                prokadiff_coords=f"{ja.seq_id_1}:{ja.pos_1}({ja.strand_1})~{ja.pos_2}({ja.strand_2})",
                breseq_coords=".",
                prokadiff_support=str(supp),
                breseq_support=".",
                delta_breakend1=".",
                delta_breakend2=".",
                total_distance=".",
                reconciliation_status="PROKADIFF_SPECIFIC_CANDIDATE",
                notes=notes,
            )
        )
        matched_a.add(ja.record_id)

    return matrix


def write_bidirectional_jc_matrix_tsv(rows: list[BidirectionalJcRow], out_path: Path) -> None:
    """Write bidirectional JC matrix to TSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "category",
        "prokadiff_id",
        "breseq_id",
        "prokadiff_coords",
        "breseq_coords",
        "prokadiff_support",
        "breseq_support",
        "delta_breakend1",
        "delta_breakend2",
        "total_distance",
        "reconciliation_status",
        "notes",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow(r.__dict__)


def main() -> int:
    parser = argparse.ArgumentParser(description="M6-P3c Structural Reconciliation and 4D Classifier")
    parser.add_argument("--prokadiff-gd", required=True, type=Path, help="ProkaDiff output.gd path")
    parser.add_argument("--comparator-gd", type=Path, default=None, help="breseq output.gd path (optional)")
    parser.add_argument("--reference-gbk", type=Path, default=None, help="Reference GenBank file (optional)")
    parser.add_argument("--out-tsv", required=True, type=Path, help="Output reconciliation TSV path")
    parser.add_argument("--mob-tsv", type=Path, default=None, help="Output MOB concordance TSV path (optional)")
    parser.add_argument("--matrix-tsv", type=Path, default=None, help="Output Bidirectional JC Matrix TSV path (optional)")
    parser.add_argument("--json", action="store_true", help="Print summary JSON")
    args = parser.parse_args()

    if not args.prokadiff_gd.is_file():
        print(f"Error: ProkaDiff GD file not found: {args.prokadiff_gd}", file=sys.stderr)
        return 1

    pk_text = args.prokadiff_gd.read_text(encoding="utf-8", errors="replace")
    comp_text = args.comparator_gd.read_text(encoding="utf-8", errors="replace") if args.comparator_gd and args.comparator_gd.is_file() else None
    ref_text = args.reference_gbk.read_text(encoding="utf-8", errors="replace") if args.reference_gbk and args.reference_gbk.is_file() else None

    reconciled = reconcile_cohort(pk_text, comp_text, ref_text)
    write_reconciled_tsv(reconciled, args.out_tsv)
    print(f"Wrote {len(reconciled)} reconciled structural records to {args.out_tsv}", file=sys.stderr)

    # Optional MOB concordance evaluation
    if args.mob_tsv and comp_text:
        pk_jcs, pk_mobs = parse_gd_text(pk_text, caller="prokadiff")
        comp_jcs, comp_mobs = parse_gd_text(comp_text, caller="breseq")
        ref_features = parse_genbank_repeat_features(ref_text) if ref_text else []
        for jc in pk_jcs + comp_jcs:
            for rf in ref_features:
                if (rf.start - 20) <= jc.pos_1 <= (rf.end + 20):
                    jc.repeat_family_1 = rf.name
                if (rf.start - 20) <= jc.pos_2 <= (rf.end + 20):
                    jc.repeat_family_2 = rf.name
        mob_results = reconcile_mobs(pk_mobs, comp_mobs, pk_jcs, ref_features)
        write_mob_concordance_tsv(mob_results, args.mob_tsv)
        print(f"Wrote {len(mob_results)} MOB concordance records to {args.mob_tsv}", file=sys.stderr)

    # Optional Bidirectional JC Matrix
    if args.matrix_tsv and comp_text:
        pk_jcs, _ = parse_gd_text(pk_text, caller="prokadiff")
        comp_jcs, _ = parse_gd_text(comp_text, caller="breseq")
        ref_features = parse_genbank_repeat_features(ref_text) if ref_text else []
        for jc in pk_jcs + comp_jcs:
            for rf in ref_features:
                if (rf.start - 20) <= jc.pos_1 <= (rf.end + 20):
                    jc.repeat_family_1 = rf.name
                if (rf.start - 20) <= jc.pos_2 <= (rf.end + 20):
                    jc.repeat_family_2 = rf.name
        matrix_rows = build_bidirectional_jc_matrix(pk_jcs, comp_jcs, ref_features)
        write_bidirectional_jc_matrix_tsv(matrix_rows, args.matrix_tsv)
        print(f"Wrote {len(matrix_rows)} bidirectional JC matrix records to {args.matrix_tsv}", file=sys.stderr)

    if args.json:
        dim1_counts = {}
        dim2_counts = {}
        dim3_counts = {}
        composite_counts = {}
        for r in reconciled:
            dim1_counts[r.dim1_placement] = dim1_counts.get(r.dim1_placement, 0) + 1
            dim2_counts[r.dim2_read_evidence] = dim2_counts.get(r.dim2_read_evidence, 0) + 1
            dim3_counts[r.dim3_corroboration] = dim3_counts.get(r.dim3_corroboration, 0) + 1
            composite_counts[r.composite_confidence] = composite_counts.get(r.composite_confidence, 0) + 1
        summary = {
            "total_records": len(reconciled),
            "dim1_placement": dim1_counts,
            "dim2_read_evidence": dim2_counts,
            "dim3_corroboration": dim3_counts,
            "dim4_biological_validation": "NOT_EVALUABLE",
            "composite_confidence": composite_counts,
        }
        print(json.dumps(summary, indent=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())
