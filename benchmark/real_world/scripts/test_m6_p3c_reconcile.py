import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from reconcile_m6_p3c import (
    BidirectionalJcRow,
    CompatibilityEdge,
    JunctionRecord,
    MatchingComparisonSummary,
    MobConcordanceResult,
    MobRecord,
    RepeatFeature,
    build_bidirectional_jc_matrix,
    build_junction_compatibility_graph,
    classify_genomic_placement,
    classify_read_evidence,
    compare_matching_algorithms,
    compute_candidate_multiplicity,
    evaluate_mob_pair,
    evaluate_pair,
    match_junctions_greedy,
    match_junctions_maximum_cardinality,
    parse_gd_text,
    parse_genbank_repeat_features,
    reconcile_cohort,
    reconcile_mobs,
    write_bidirectional_jc_matrix_tsv,
    write_mob_concordance_tsv,
    write_reconciled_tsv,
)


def test_case_1_exact_match():
    """Case 1: Exact Match (D=0)."""
    ja = JunctionRecord(
        caller="prokadiff",
        record_id="101",
        seq_id_1="NC_000913.3",
        pos_1=1000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=2500,
        strand_2=-1,
        attributes={"pd_support_reads": "20", "pd_plus_reads": "10", "pd_minus_reads": "10", "pd_min_overlap_side1": "25", "pd_min_overlap_side2": "25"},
    )
    jb = JunctionRecord(
        caller="breseq",
        record_id="B1",
        seq_id_1="NC_000913.3",
        pos_1=1000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=2500,
        strand_2=-1,
    )
    is_match, mtype, dist, reason = evaluate_pair(ja, jb, tolerance=5)
    assert is_match is True
    assert mtype == "EXACT_MATCH"
    assert dist == 0
    assert reason == "NONE"

    matches = match_junctions_greedy([ja], [jb], tolerance=5)
    assert matches["101"].match_type == "EXACT_MATCH"
    assert matches["101"].comparator_id == "B1"
    assert matches["101"].distance == 0


def test_case_2_swapped_breakend_order():
    """Case 2: Swapped Breakend Order (D=0, symmetric equivalence)."""
    ja = JunctionRecord(
        caller="prokadiff",
        record_id="102",
        seq_id_1="NC_000913.3",
        pos_1=1500,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=3000,
        strand_2=-1,
    )
    jb = JunctionRecord(
        caller="breseq",
        record_id="B2",
        seq_id_1="NC_000913.3",
        pos_1=3000,
        strand_1=-1,
        seq_id_2="NC_000913.3",
        pos_2=1500,
        strand_2=1,
    )
    is_match, mtype, dist, reason = evaluate_pair(ja, jb, tolerance=5)
    assert is_match is True
    assert mtype == "EQUIVALENT_SWAPPED_MATCH"
    assert dist == 0
    assert reason == "NONE"

    matches = match_junctions_greedy([ja], [jb], tolerance=5)
    assert matches["102"].match_type == "EQUIVALENT_SWAPPED_MATCH"
    assert matches["102"].comparator_id == "B2"
    assert matches["102"].distance == 0


def test_case_3_orientation_strand_mismatch():
    """Case 3: Orientation / Strand Mismatch (TOPOLOGY_MISMATCH, not a match)."""
    ja = JunctionRecord(
        caller="prokadiff",
        record_id="103",
        seq_id_1="NC_000913.3",
        pos_1=1000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=2500,
        strand_2=-1,  # Deletion topology
    )
    jb = JunctionRecord(
        caller="breseq",
        record_id="B3",
        seq_id_1="NC_000913.3",
        pos_1=1000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=2500,
        strand_2=1,   # Inversion topology
    )
    is_match, mtype, dist, reason = evaluate_pair(ja, jb, tolerance=5)
    assert is_match is False
    assert reason == "TOPOLOGY_MISMATCH"

    matches = match_junctions_greedy([ja], [jb], tolerance=5)
    assert matches["103"].match_type == "NONE"
    assert matches["103"].comparator_id == "B3"
    assert matches["103"].mismatch_reason == "TOPOLOGY_MISMATCH"


def test_case_4_single_breakpoint_match():
    """Case 4: Single Breakpoint Match (DISCORDANT_BREAKEND, not a match)."""
    ja = JunctionRecord(
        caller="prokadiff",
        record_id="104",
        seq_id_1="NC_000913.3",
        pos_1=1000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=2500,
        strand_2=-1,
    )
    jb = JunctionRecord(
        caller="breseq",
        record_id="B4",
        seq_id_1="NC_000913.3",
        pos_1=1000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=8500,
        strand_2=-1,  # Discordant second breakend
    )
    is_match, mtype, dist, reason = evaluate_pair(ja, jb, tolerance=5)
    assert is_match is False
    assert reason == "DISCORDANT_BREAKEND"

    matches = match_junctions_greedy([ja], [jb], tolerance=5)
    assert matches["104"].match_type == "NONE"
    assert matches["104"].comparator_id == "B4"
    assert matches["104"].mismatch_reason == "DISCORDANT_BREAKEND"


def test_case_5_competing_candidates_greedy_tiebreaker():
    """Case 5: Competing Candidates within microhomology window (minimal distance tiebreaker)."""
    ja = JunctionRecord(
        caller="prokadiff",
        record_id="105",
        seq_id_1="NC_000913.3",
        pos_1=1000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=2000,
        strand_2=-1,
        attributes={"pd_support_reads": "25"},
    )
    jb1 = JunctionRecord(
        caller="breseq",
        record_id="B5_1",
        seq_id_1="NC_000913.3",
        pos_1=1002,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=2001,
        strand_2=-1,  # D = |1000-1002| + |2000-2001| = 3
    )
    jb2 = JunctionRecord(
        caller="breseq",
        record_id="B5_2",
        seq_id_1="NC_000913.3",
        pos_1=1004,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=2000,
        strand_2=-1,  # D = |1000-1004| + |2000-2000| = 4
    )

    matches = match_junctions_greedy([ja], [jb1, jb2], tolerance=5)
    assert matches["105"].comparator_id == "B5_1"
    assert matches["105"].distance == 3
    assert matches["105"].match_type == "NEAR_MATCH"


def test_case_6_repeat_ambiguity():
    """Case 6: Repeat Ambiguity (matching anchor, differing repeat copy locus of same family)."""
    ja = JunctionRecord(
        caller="prokadiff",
        record_id="106",
        seq_id_1="NC_000913.3",
        pos_1=5000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=250000,
        strand_2=-1,
        repeat_family_1=None,
        repeat_family_2="IS1",
    )
    jb = JunctionRecord(
        caller="breseq",
        record_id="B6",
        seq_id_1="NC_000913.3",
        pos_1=5000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=1400000,
        strand_2=-1,
        repeat_family_1=None,
        repeat_family_2="IS1",
    )
    # Pos 5000 is unique anchor; side 2 is IS1 repeat family at different loci
    is_match, mtype, dist, reason = evaluate_pair(ja, jb, tolerance=5)
    assert is_match is False
    assert reason == "AMBIGUOUS_FAMILY_MATCH"

    matches = match_junctions_greedy([ja], [jb], tolerance=5)
    assert matches["106"].mismatch_reason == "AMBIGUOUS_FAMILY_MATCH"


def test_case_7_no_comparator_match():
    """Case 7: No Comparator Match (NOT_CORROBORATED)."""
    ja = JunctionRecord(
        caller="prokadiff",
        record_id="107",
        seq_id_1="NC_000913.3",
        pos_1=345000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=345800,
        strand_2=-1,
    )
    matches = match_junctions_greedy([ja], [], tolerance=5)
    assert matches["107"].match_type == "NONE"
    assert matches["107"].comparator_id is None
    assert matches["107"].mismatch_reason == "NO_NEAR_CANDIDATE"


def test_dimension_1_genomic_placement():
    """Test Dimension 1: UNIQUE_ANCHOR, REPEAT_ANCHORED, REPEAT_AMBIGUOUS."""
    repeats = [
        RepeatFeature(seq_id="NC_000913.3", start=1000, end=2000, feature_type="repeat_region", name="REP1"),
        RepeatFeature(seq_id="NC_000913.3", start=5000, end=6000, feature_type="mobile_element", name="IS1"),
    ]
    # Both breakends non-repeat
    j_unique = JunctionRecord(
        caller="prokadiff", record_id="1", seq_id_1="NC_000913.3", pos_1=3000, strand_1=1,
        seq_id_2="NC_000913.3", pos_2=4000, strand_2=-1
    )
    assert classify_genomic_placement(j_unique, repeats) == "UNIQUE_ANCHOR"

    # Exactly one breakend in repeat (within 20 bp)
    j_anchored = JunctionRecord(
        caller="prokadiff", record_id="2", seq_id_1="NC_000913.3", pos_1=1010, strand_1=1,
        seq_id_2="NC_000913.3", pos_2=4000, strand_2=-1
    )
    assert classify_genomic_placement(j_anchored, repeats) == "REPEAT_ANCHORED"

    # Both breakends in repeat
    j_ambiguous = JunctionRecord(
        caller="prokadiff", record_id="3", seq_id_1="NC_000913.3", pos_1=1500, strand_1=1,
        seq_id_2="NC_000913.3", pos_2=5500, strand_2=-1
    )
    assert classify_genomic_placement(j_ambiguous, repeats) == "REPEAT_AMBIGUOUS"


def test_dimension_2_read_evidence_classification():
    """Test Dimension 2: STRONG, LIMITED, UNRESOLVED."""
    # STRONG: N>=10, plus>=2, minus>=2, overlap1>=14, overlap2>=14
    j_strong = JunctionRecord(
        caller="prokadiff", record_id="1", seq_id_1="c", pos_1=1, strand_1=1, seq_id_2="c", pos_2=2, strand_2=-1,
        attributes={"pd_support_reads": "15", "pd_plus_reads": "8", "pd_minus_reads": "7", "pd_min_overlap_side1": "20", "pd_min_overlap_side2": "20"}
    )
    assert classify_read_evidence(j_strong) == "STRONG"

    # LIMITED: N=15, but unilateral strand support (plus=15, minus=0)
    j_limited_strand = JunctionRecord(
        caller="prokadiff", record_id="2", seq_id_1="c", pos_1=1, strand_1=1, seq_id_2="c", pos_2=2, strand_2=-1,
        attributes={"pd_support_reads": "15", "pd_plus_reads": "15", "pd_minus_reads": "0", "pd_min_overlap_side1": "20", "pd_min_overlap_side2": "20"}
    )
    assert classify_read_evidence(j_limited_strand) == "LIMITED"

    # LIMITED: N=5 (3 <= N < 10)
    j_limited_count = JunctionRecord(
        caller="prokadiff", record_id="3", seq_id_1="c", pos_1=1, strand_1=1, seq_id_2="c", pos_2=2, strand_2=-1,
        attributes={"pd_support_reads": "5", "pd_plus_reads": "3", "pd_minus_reads": "2", "pd_min_overlap_side1": "20", "pd_min_overlap_side2": "20"}
    )
    assert classify_read_evidence(j_limited_count) == "LIMITED"

    # UNRESOLVED: N=2 (< 3)
    j_unresolved = JunctionRecord(
        caller="prokadiff", record_id="4", seq_id_1="c", pos_1=1, strand_1=1, seq_id_2="c", pos_2=2, strand_2=-1,
        attributes={"pd_support_reads": "2", "pd_plus_reads": "1", "pd_minus_reads": "1", "pd_min_overlap_side1": "20", "pd_min_overlap_side2": "20"}
    )
    assert classify_read_evidence(j_unresolved) == "UNRESOLVED"


def test_full_reconcile_and_tsv_writing():
    """Test full cohort reconciliation and TSV generation."""
    pk_gd = """# Genome Diff
JC\t1\t.\tNC_000913.3\t1000\t1\tNC_000913.3\t2500\t-1\t0\tpd_support_reads=25\tpd_plus_reads=12\tpd_minus_reads=13\tpd_min_overlap_side1=30\tpd_min_overlap_side2=30
JC\t2\t.\tNC_000913.3\t4000\t1\tNC_000913.3\t5000\t-1\t0\tpd_support_reads=15\tpd_plus_reads=8\tpd_minus_reads=7\tpd_min_overlap_side1=20\tpd_min_overlap_side2=20
"""
    comp_gd = """# breseq Genome Diff
JC\tB1\t.\tNC_000913.3\t1000\t1\tNC_000913.3\t2500\t-1\t0\treads_total=30
JC\tB2\t.\tNC_000913.3\t4000\t1\tNC_000913.3\t5000\t1\t0\treads_total=20
"""
    reconciled = reconcile_cohort(pk_gd, comp_gd, None)
    assert len(reconciled) == 2

    # Record 1: exact match with B1
    assert reconciled[0].record_id == "1"
    assert reconciled[0].matched_comparator_id == "B1"
    assert reconciled[0].match_type == "EXACT_MATCH"
    assert reconciled[0].dim1_placement == "UNIQUE_ANCHOR"
    assert reconciled[0].dim2_read_evidence == "STRONG"
    assert reconciled[0].dim3_corroboration == "CORROBORATED"
    assert reconciled[0].dim4_biological_validation == "NOT_EVALUABLE"
    assert reconciled[0].composite_confidence == "HIGH_CONFIDENCE_COMPUTATIONAL_CANDIDATE"

    # Record 2: topology mismatch with B2 (differing strands: -1 vs 1)
    assert reconciled[1].record_id == "2"
    assert reconciled[1].mismatch_reason == "TOPOLOGY_MISMATCH"
    assert reconciled[1].dim3_corroboration == "NOT_CORROBORATED"
    assert reconciled[1].composite_confidence == "CANDIDATE"

    # Test writing TSV
    with tempfile.NamedTemporaryFile("w+", suffix=".tsv") as tf:
        out_p = Path(tf.name)
        write_reconciled_tsv(reconciled, out_p)
        rows = list(csv.DictReader(out_p.open(), delimiter="\t"))
        assert len(rows) == 2
        assert rows[0]["record_id"] == "1"
        assert rows[0]["dim3_corroboration"] == "CORROBORATED"
        assert rows[1]["dim3_corroboration"] == "NOT_CORROBORATED"


def test_cli_execution():
    """Test full CLI entrypoint via subprocess."""
    pk_gd = "JC\t1\t.\tNC_000913.3\t1000\t1\tNC_000913.3\t2500\t-1\t0\tpd_support_reads=25\n"
    comp_gd = "JC\tB1\t.\tNC_000913.3\t1000\t1\tNC_000913.3\t2500\t-1\t0\treads_total=30\n"

    with tempfile.TemporaryDirectory() as tmpdir:
        pk_p = Path(tmpdir) / "prokadiff.gd"
        comp_p = Path(tmpdir) / "breseq.gd"
        out_tsv = Path(tmpdir) / "reconciled.tsv"

        pk_p.write_text(pk_gd)
        comp_p.write_text(comp_gd)

        cmd = [
            sys.executable,
            str(Path(__file__).parent / "reconcile_m6_p3c.py"),
            "--prokadiff-gd",
            str(pk_p),
            "--comparator-gd",
            str(comp_p),
            "--out-tsv",
            str(out_tsv),
            "--json",
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        assert res.returncode == 0, f"CLI failed: {res.stderr}"
        assert out_tsv.is_file()
        summary = json.loads(res.stdout)
        assert summary["total_records"] == 1
        assert summary["dim3_corroboration"]["CORROBORATED"] == 1


def test_greedy_vs_maximum_cardinality_counterexample():
    """Verify synthetic counterexample where greedy matching yields strictly smaller cardinality.

    Theoretical scenario under tolerance <= 5 bp:
      JA1 at pos 1000, JA2 at pos 1004.
      JB1 at pos 1001, JB2 at pos 997.
      Breakend 2 is identical for all at pos 5000 (-1).

    Pairwise distances:
      - (JA1, JB1): D = |1000-1001| = 1 <= 5 (valid)
      - (JA1, JB2): D = |1000-997|  = 3 <= 5 (valid)
      - (JA2, JB1): D = |1004-1001| = 3 <= 5 (valid)
      - (JA2, JB2): D = |1004-997|  = 7 > 5  (INVALID: exceeds tolerance)

    1. Maximum-Cardinality Matching:
       - Match (JA1, JB2) [D=3] and (JA2, JB1) [D=3].
       - Both JA1 and JA2 are matched -> Cardinality = 2.

    2. Greedy Bipartite Matching:
       - Greedy prioritizes candidate pairs by minimal distance D.
       - Candidate (JA1, JB1) has minimal D=1, so it is paired first.
       - JB1 is consumed by JA1.
       - JA2 cannot match JB2 (D=7 > 5), so JA2 remains unmatched.
       - Greedy Cardinality = 1 < Maximum Cardinality = 2.
    """
    ja1 = JunctionRecord(
        caller="prokadiff",
        record_id="JA1",
        seq_id_1="NC_000913.3",
        pos_1=1000,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=5000,
        strand_2=-1,
        attributes={"pd_support_reads": "20"},
    )
    ja2 = JunctionRecord(
        caller="prokadiff",
        record_id="JA2",
        seq_id_1="NC_000913.3",
        pos_1=1004,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=5000,
        strand_2=-1,
        attributes={"pd_support_reads": "20"},
    )
    jb1 = JunctionRecord(
        caller="breseq",
        record_id="JB1",
        seq_id_1="NC_000913.3",
        pos_1=1001,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=5000,
        strand_2=-1,
        attributes={"reads_total": "20"},
    )
    jb2 = JunctionRecord(
        caller="breseq",
        record_id="JB2",
        seq_id_1="NC_000913.3",
        pos_1=997,
        strand_1=1,
        seq_id_2="NC_000913.3",
        pos_2=5000,
        strand_2=-1,
        attributes={"reads_total": "20"},
    )

    # 1. Verify individual pair compatibility
    m11, _, d11, _ = evaluate_pair(ja1, jb1, tolerance=5)
    m12, _, d12, _ = evaluate_pair(ja1, jb2, tolerance=5)
    m21, _, d21, _ = evaluate_pair(ja2, jb1, tolerance=5)
    m22, _, d22, _ = evaluate_pair(ja2, jb2, tolerance=5)

    assert m11 is True and d11 == 1
    assert m12 is True and d12 == 3
    assert m21 is True and d21 == 3
    assert m22 is False  # D=7 > 5 tolerance

    # 2. Maximum cardinality matching exists and has size 2
    # Assignment: JA1 -> JB2 (valid), JA2 -> JB1 (valid)
    assert m12 is True and m21 is True
    max_cardinality = 2

    # 3. Execute greedy bipartite matching
    greedy_results = match_junctions_greedy([ja1, ja2], [jb1, jb2], tolerance=5)

    greedy_matched = [r for r in greedy_results.values() if r.match_type != "NONE"]
    greedy_cardinality = len(greedy_matched)

    # Greedy paired JA1 with JB1 (min distance 1)
    assert greedy_results["JA1"].match_type in ("NEAR_MATCH", "EXACT_MATCH")
    assert greedy_results["JA1"].comparator_id == "JB1"
    assert greedy_results["JA1"].distance == 1

    # JA2 was starved and left unmatched
    assert greedy_results["JA2"].match_type == "NONE"
    assert greedy_results["JA2"].comparator_id is None or greedy_results["JA2"].mismatch_reason != "NONE"

    # Strict cardinality difference demonstrated
    assert greedy_cardinality == 1
    assert greedy_cardinality < max_cardinality

    # 4. Execute maximum-cardinality bipartite matching (Kuhn's augmenting path algorithm)
    max_results = match_junctions_maximum_cardinality([ja1, ja2], [jb1, jb2], tolerance=5)
    max_matched = [r for r in max_results.values() if r.match_type != "NONE"]
    assert len(max_matched) == 2

    # Verify assignments: JA1 -> JB2 (dist=3), JA2 -> JB1 (dist=3)
    assert max_results["JA1"].comparator_id == "JB2"
    assert max_results["JA1"].distance == 3
    assert max_results["JA2"].comparator_id == "JB1"
    assert max_results["JA2"].distance == 3

    # 5. Verify compare_matching_algorithms summary and candidate multiplicity
    pk_mult, comp_mult = compute_candidate_multiplicity([ja1, ja2], [jb1, jb2], tolerance=5)
    assert pk_mult["JA1"] == 2
    assert pk_mult["JA2"] == 1
    assert comp_mult["JB1"] == 2
    assert comp_mult["JB2"] == 1

    _, _, summary = compare_matching_algorithms([ja1, ja2], [jb1, jb2], tolerance=5)
    assert summary.greedy_cardinality == 1
    assert summary.max_cardinality == 2
    assert summary.cardinality_difference == 1
    assert summary.is_identical is False
    assert summary.prokadiff_multiplicity_gt1 == 1
    assert summary.comparator_multiplicity_gt1 == 1


def test_maximum_cardinality_and_multiplicity():
    """Verify compatibility graph builder, multiplicity metrics, and CLI matching mode."""
    ja1 = JunctionRecord(
        caller="prokadiff", record_id="P1", seq_id_1="NC_000913", pos_1=1000, strand_1=1,
        seq_id_2="NC_000913", pos_2=2000, strand_2=-1, attributes={"pd_support_reads": "30"}
    )
    ja2 = JunctionRecord(
        caller="prokadiff", record_id="P2", seq_id_1="NC_000913", pos_1=3000, strand_1=1,
        seq_id_2="NC_000913", pos_2=4000, strand_2=-1, attributes={"pd_support_reads": "25"}
    )
    jb1 = JunctionRecord(
        caller="breseq", record_id="B1", seq_id_1="NC_000913", pos_1=1000, strand_1=1,
        seq_id_2="NC_000913", pos_2=2000, strand_2=-1, attributes={"reads_total": "30"}
    )
    jb2 = JunctionRecord(
        caller="breseq", record_id="B2", seq_id_1="NC_000913", pos_1=3001, strand_1=1,
        seq_id_2="NC_000913", pos_2=4001, strand_2=-1, attributes={"reads_total": "25"}
    )

    graph = build_junction_compatibility_graph([ja1, ja2], [jb1, jb2], tolerance=5)
    assert len(graph["P1"]) == 1
    assert graph["P1"][0].comparator_id == "B1"
    assert graph["P1"][0].distance == 0
    assert len(graph["P2"]) == 1
    assert graph["P2"][0].comparator_id == "B2"
    assert graph["P2"][0].distance == 2

    pk_mult, comp_mult = compute_candidate_multiplicity([ja1, ja2], [jb1, jb2], tolerance=5)
    assert pk_mult["P1"] == 1 and pk_mult["P2"] == 1
    assert comp_mult["B1"] == 1 and comp_mult["B2"] == 1

    _, _, summary = compare_matching_algorithms([ja1, ja2], [jb1, jb2], tolerance=5)
    assert summary.greedy_cardinality == 2
    assert summary.max_cardinality == 2
    assert summary.cardinality_difference == 0
    assert summary.is_identical is True
    assert summary.prokadiff_multiplicity_gt1 == 0
    assert summary.comparator_multiplicity_gt1 == 0


def test_decoupled_mob_concordance_evaluator():
    """Verify decoupled 4-attribute MOB concordance evaluation."""
    # 1. Exact MOB match (IS5 at 889,647, dup=4)
    m_exact_a = MobRecord(
        caller="prokadiff",
        record_id="148",
        parent_ids="120,121",
        seq_id="NC_000913.3",
        pos=889647,
        repeat_name="IS5",
        strand=-1,
        duplication_size=4,
    )
    m_exact_b = MobRecord(
        caller="breseq",
        record_id="4",
        parent_ids="35,38",
        seq_id="NC_000913",
        pos=889647,
        repeat_name="IS5",
        strand=-1,
        duplication_size=4,
    )
    res_exact = evaluate_mob_pair(m_exact_a, m_exact_b)
    assert res_exact.site_concordance == "EXACT_SITE"
    assert res_exact.family_concordance == "EXACT_FAMILY"
    assert res_exact.strand_concordance == "EXACT_STRAND"
    assert res_exact.tsd_concordance == "EXACT_TSD"
    assert res_exact.composite_classification == "EXACT_MOB_MATCH"

    # 2. 1 bp TSD slip resolution (IS1 at 1,873,031: ProkaDiff dup=8 vs breseq dup=9)
    # Must NOT be classified as EXACT_MOB_MATCH!
    m_slip_a = MobRecord(
        caller="prokadiff",
        record_id="147",
        parent_ids="106,132",
        seq_id="NC_000913.3",
        pos=1873031,
        repeat_name="IS1",
        strand=-1,
        duplication_size=8,
    )
    m_slip_b = MobRecord(
        caller="breseq",
        record_id="5",
        parent_ids="34,42",
        seq_id="NC_000913",
        pos=1873031,
        repeat_name="IS1",
        strand=-1,
        duplication_size=9,
    )
    res_slip = evaluate_mob_pair(m_slip_a, m_slip_b)
    assert res_slip.site_concordance == "EXACT_SITE"
    assert res_slip.family_concordance == "EXACT_FAMILY"
    assert res_slip.strand_concordance == "EXACT_STRAND"
    assert res_slip.tsd_concordance == "NEAR_MATCH_1BP"
    assert res_slip.composite_classification == "CONCORDANT_MOB_1BP_TSD_RESOLVED"

    # 3. Discordant site / strand
    m_disc_b = MobRecord(
        caller="breseq",
        record_id="99",
        parent_ids=".",
        seq_id="NC_000913",
        pos=1873031,
        repeat_name="IS1",
        strand=1,  # inverted strand
        duplication_size=8,
    )
    res_disc = evaluate_mob_pair(m_slip_a, m_disc_b)
    assert res_disc.strand_concordance == "INVERTED_STRAND"
    assert res_disc.composite_classification == "DISCORDANT"

    # 4. Cohort reconciliation and TSV writing
    mob_results = reconcile_mobs([m_exact_a, m_slip_a], [m_exact_b, m_slip_b])
    assert len(mob_results) == 2
    assert mob_results[0].composite_classification == "EXACT_MOB_MATCH"
    assert mob_results[1].composite_classification == "CONCORDANT_MOB_1BP_TSD_RESOLVED"

    with tempfile.NamedTemporaryFile("w+", suffix=".tsv") as tf:
        out_p = Path(tf.name)
        write_mob_concordance_tsv(mob_results, out_p)
        rows = list(csv.DictReader(out_p.open(), delimiter="\t"))
        assert len(rows) == 2
        assert rows[0]["composite_classification"] == "EXACT_MOB_MATCH"
        assert rows[1]["composite_classification"] == "CONCORDANT_MOB_1BP_TSD_RESOLVED"


def test_bidirectional_jc_matrix_partitioning():
    """Verify bidirectional JC accounting matrix partitions all JCs without duplicates."""
    # Synthetic mini-cohort covering all 5 categories:
    # 1. Exact match (D=0)
    ja_exact = JunctionRecord(
        caller="prokadiff", record_id="P1", seq_id_1="NC_000913", pos_1=1000, strand_1=1,
        seq_id_2="NC_000913", pos_2=2000, strand_2=-1, attributes={"pd_support_reads": "30"}
    )
    jb_exact = JunctionRecord(
        caller="breseq", record_id="B1", seq_id_1="NC_000913", pos_1=1000, strand_1=1,
        seq_id_2="NC_000913", pos_2=2000, strand_2=-1, attributes={"reads_total": "30"}
    )

    # 2. Case 6 Repeat-family match
    ja_rep = JunctionRecord(
        caller="prokadiff", record_id="P2", seq_id_1="NC_000913", pos_1=5000, strand_1=1,
        seq_id_2="NC_000913", pos_2=250000, strand_2=-1, repeat_family_2="IS1", attributes={"pd_support_reads": "20"}
    )
    jb_rep = JunctionRecord(
        caller="breseq", record_id="B2", seq_id_1="NC_000913", pos_1=5000, strand_1=1,
        seq_id_2="NC_000913", pos_2=1400000, strand_2=-1, repeat_family_2="IS1", attributes={"reads_total": "20"}
    )

    # 3. Near-miss boundary match (breakend 1 delta=6 bp > 5 bp, D=10)
    ja_near = JunctionRecord(
        caller="prokadiff", record_id="P3", seq_id_1="NC_000913", pos_1=1207784, strand_1=1,
        seq_id_2="NC_000913", pos_2=1209623, strand_2=1, attributes={"pd_support_reads": "82"}
    )
    jb_near = JunctionRecord(
        caller="breseq", record_id="B3", seq_id_1="NC_000913", pos_1=1207790, strand_1=1,
        seq_id_2="NC_000913", pos_2=1209619, strand_2=1, attributes={"reads_total": "36"}
    )

    # 4. breseq-only rejected
    jb_rej = JunctionRecord(
        caller="breseq", record_id="B4", seq_id_1="NC_000913", pos_1=800000, strand_1=1,
        seq_id_2="NC_000913", pos_2=805000, strand_2=-1, attributes={"reject": "COVERAGE_EVENNESS_SKEW"}
    )

    # 5. ProkaDiff-only candidate
    ja_only = JunctionRecord(
        caller="prokadiff", record_id="P4", seq_id_1="NC_000913", pos_1=900000, strand_1=1,
        seq_id_2="NC_000913", pos_2=901000, strand_2=-1, attributes={"pd_support_reads": "15"}
    )

    pjcs = [ja_exact, ja_rep, ja_near, ja_only]
    bjcs = [jb_exact, jb_rep, jb_near, jb_rej]

    matrix = build_bidirectional_jc_matrix(pjcs, bjcs)
    assert len(matrix) == 5

    categories = [r.category for r in matrix]
    assert "EXACT_MATCH" in categories
    assert "REPEAT_FAMILY_MATCH" in categories
    assert "NEAR_MISS_BOUNDARY" in categories
    assert "BRESEQ_ONLY_REJECTED" in categories
    assert "PROKADIFF_ONLY" in categories

    # Verify no duplicate assignments
    pk_assigned = [r.prokadiff_id for r in matrix if r.prokadiff_id != "."]
    b_assigned = [r.breseq_id for r in matrix if r.breseq_id != "."]
    assert len(pk_assigned) == len(set(pk_assigned)) == 4
    assert len(b_assigned) == len(set(b_assigned)) == 4
