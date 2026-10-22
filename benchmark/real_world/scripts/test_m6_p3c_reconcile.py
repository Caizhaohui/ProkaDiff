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


def test_reference_hashes_raw_and_normalized():
    """Verify bit-identical CP009273.1 reference files and reject transcription artifacts."""
    import hashlib

    fa_p = REPO / "data/references/m6_identity/CP009273.1.fasta"
    gbk_p = REPO / "data/references/m6_identity/CP009273.1.gbk"
    assert fa_p.is_file(), f"Missing {fa_p}"
    assert gbk_p.is_file(), f"Missing {gbk_p}"

    raw_fa_hash = hashlib.sha256(fa_p.read_bytes()).hexdigest()
    raw_gbk_hash = hashlib.sha256(gbk_p.read_bytes()).hexdigest()

    assert raw_fa_hash == "37eade01311f1f81659fad53bf907a93139075a75d115fd3aec2ef09beac2e90"
    assert raw_gbk_hash == "3440cc8d28d643dda30e8ed9054656bc093abf61895077678ef4ddad4343df3a"

    # Normalized sequence SHA-256 (length 4,631,469 bp)
    seq = "".join(line.strip() for line in fa_p.read_text(encoding="utf-8").splitlines() if not line.startswith(">"))
    assert len(seq) == 4631469
    norm_seq_hash = hashlib.sha256(seq.upper().encode("ascii")).hexdigest()
    assert norm_seq_hash == "3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce"

    # Ensure deprecated / corrupted transcription artifact hash 94ee6a51... is rejected
    corrupted_hash = "94ee6a512c1d0202720d2d3a95c9a1738743a68ea107b3c2e176ebc24f5a3e14"
    assert raw_fa_hash != corrupted_hash
    assert raw_gbk_hash != corrupted_hash
    assert norm_seq_hash != corrupted_hash


def test_bidirectional_matrix_split_corroborated_vs_near_miss():
    """Verify 6-category bidirectional matrix split: 101 corroborated (91 exact + 10 near) vs 10 near-miss."""
    matrix_p = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3c_cohort_bidirectional_matrix.tsv"
    assert matrix_p.is_file(), f"Missing manifest: {matrix_p}"

    with matrix_p.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    assert len(rows) == 1319

    from collections import Counter
    counts = Counter(r["category"] for r in rows)

    assert counts["EXACT_MATCH"] == 91
    assert counts["NEAR_MATCH_WITHIN_TOLERANCE"] == 10
    assert counts["EXACT_MATCH"] + counts["NEAR_MATCH_WITHIN_TOLERANCE"] == 101
    assert counts["REPEAT_FAMILY_MATCH"] == 15
    assert counts["NEAR_MISS_BOUNDARY"] == 10
    assert counts["BRESEQ_ONLY_REJECTED"] == 513
    assert counts["PROKADIFF_ONLY"] == 680


def test_cohort_bijective_accounting_conservation():
    """Verify exact edge-level accounting conservation across ProkaDiff, breseq, and manifests."""
    matrix_p = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3c_cohort_bidirectional_matrix.tsv"
    telemetry_p = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3c_cohort_telemetry.tsv"
    reconciliation_p = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3c_cohort_reconciliation.tsv"

    with matrix_p.open(encoding="utf-8") as f:
        mat_rows = list(csv.DictReader(f, delimiter="\t"))

    with telemetry_p.open(encoding="utf-8") as f:
        tel_rows = list(csv.DictReader(f, delimiter="\t"))

    with reconciliation_p.open(encoding="utf-8") as f:
        rec_rows = list(csv.DictReader(f, delimiter="\t"))

    # Matrix sums
    from collections import Counter
    cat_counts = Counter(r["category"] for r in mat_rows)
    pk_sum = (cat_counts["EXACT_MATCH"] + cat_counts["NEAR_MATCH_WITHIN_TOLERANCE"] +
              cat_counts["REPEAT_FAMILY_MATCH"] + cat_counts["NEAR_MISS_BOUNDARY"] +
              cat_counts["PROKADIFF_ONLY"])
    b_sum = (cat_counts["EXACT_MATCH"] + cat_counts["NEAR_MATCH_WITHIN_TOLERANCE"] +
             cat_counts["REPEAT_FAMILY_MATCH"] + cat_counts["NEAR_MISS_BOUNDARY"] +
             cat_counts["BRESEQ_ONLY_REJECTED"])

    assert pk_sum == 806
    assert b_sum == 639
    assert len(mat_rows) == 1319
    assert len(rec_rows) == 806

    # Verify 1:1 bijectivity within each run in matrix
    runs = set(r["run"] for r in mat_rows)
    assert len(runs) == 5
    for r in runs:
        run_rows = [x for x in mat_rows if x["run"] == r]
        pk_ids = [x["prokadiff_id"] for x in run_rows if x["prokadiff_id"] != "."]
        b_ids = [x["breseq_id"] for x in run_rows if x["breseq_id"] != "."]
        assert len(pk_ids) == len(set(pk_ids)), f"Duplicate ProkaDiff ID in run {r}"
        assert len(b_ids) == len(set(b_ids)), f"Duplicate breseq ID in run {r}"

    # Telemetry column totals
    assert sum(int(x["exact_matches"]) for x in tel_rows) == 91
    assert sum(int(x["near_matches_within_tolerance"]) for x in tel_rows) == 10
    assert sum(int(x["dim3_corroborated"]) for x in tel_rows) == 101
    assert sum(int(x["family_matches"]) for x in tel_rows) == 15
    assert sum(int(x["dim3_family_level"]) for x in tel_rows) == 15
    assert sum(int(x["uncalled_jcs"]) for x in tel_rows) == 690
    assert sum(int(x["dim3_not_corroborated"]) for x in tel_rows) == 690
    assert sum(int(x["prokadiff_jcs"]) for x in tel_rows) == 806
    assert sum(int(x["breseq_jcs"]) for x in tel_rows) == 639


def test_mob_comparator_metrics_no_biological_claims():
    """Verify decoupled MOB comparator concordance metrics and complete absence of biological claims."""
    mob_p = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3c_cohort_mob_concordance.tsv"
    telemetry_p = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3c_cohort_telemetry.tsv"

    with mob_p.open(encoding="utf-8") as f:
        mob_rows = list(csv.DictReader(f, delimiter="\t"))

    with telemetry_p.open(encoding="utf-8") as f:
        tel_rows = list(csv.DictReader(f, delimiter="\t"))

    assert len(mob_rows) == 8

    from collections import Counter
    mc = Counter(r["composite_classification"] for r in mob_rows)
    assert mc["EXACT_MOB_MATCH"] == 6
    assert mc["CONCORDANT_MOB_1BP_TSD_RESOLVED"] == 1
    assert mc["DISCORDANT"] == 1

    # Comparator concordance: 7/8 agreement, 7/7 breseq coverage
    agreed = mc["EXACT_MOB_MATCH"] + mc["CONCORDANT_MOB_1BP_TSD_RESOLVED"]
    assert agreed == 7
    breseq_covered = sum(1 for r in mob_rows if r["comparator_id"] != ".")
    assert breseq_covered == 7

    # BW25113 samples (SRR28349231, SRR28349226) reported zero MOBs by both callers
    bw_mob_rows = [r for r in mob_rows if r["run"] in ("SRR28349231", "SRR28349226")]
    assert len(bw_mob_rows) == 0

    for r in tel_rows:
        if r["run"] in ("SRR28349231", "SRR28349226"):
            assert r["mob_concordant"] == "0/0"
        # Dim 4 biological validation must be NOT_EVALUABLE across cohort
        assert r["dim4_biological_validation"] == "NOT_EVALUABLE"

    # Strict audit: verify no biological sensitivity/specificity claims in manifests
    forbidden_terms = [
        "BIOLOGICAL_SENSITIVITY",
        "BIOLOGICAL_SPECIFICITY",
        "BIOLOGICAL_PRECISION",
        "TRUE_NEGATIVE",
    ]
    mob_text = mob_p.read_text(encoding="utf-8")
    for term in forbidden_terms:
        assert term not in mob_text, f"Forbidden biological term '{term}' found in MOB concordance manifest"


def test_no_contradictory_notes_in_matrix():
    """Verify notes field in bidirectional matrix contains no contradictory tolerance claims."""
    matrix_p = REPO / "benchmark/real_world/manifests/prjna1088182_m6_p3c_cohort_bidirectional_matrix.tsv"
    with matrix_p.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    for r in rows:
        cat = r["category"]
        notes = r["notes"]
        if cat == "NEAR_MATCH_WITHIN_TOLERANCE":
            d1 = int(r["delta_breakend1"])
            d2 = int(r["delta_breakend2"])
            assert d1 <= 5 and d2 <= 5
            assert "exceeds 5 bp tolerance" not in notes
            assert "Within 5 bp tolerance" in notes
            assert r["reconciliation_status"] == "MUTUAL_NEAR_MATCH_WITHIN_TOLERANCE"
        elif cat == "NEAR_MISS_BOUNDARY":
            d1 = int(r["delta_breakend1"])
            d2 = int(r["delta_breakend2"])
            assert max(d1, d2) >= 6
            assert "exceeds 5 bp tolerance" in notes
            assert r["reconciliation_status"] == "NEAR_MISS_MICROHOMOLOGY_SLIP"
