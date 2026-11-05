#!/usr/bin/env python3
"""
scripts/test_audit_release_gates.py

Test suite for Milestone M9-T05 (Gate E): the automated defensive release
gate tool ``scripts/audit_release_gates.py`` (hardened after the
independent Codex Sol review, 2026-11-04, and its Round 2, Round 3, and
Round 4 findings).

Coverage:
1. Repo-level pass: every gate returns zero findings on the real ProkaDiff
   checkout (scientific language, evidence chain, quarantine/G1 hold).
2. Negative unit tests on ``tmp_path`` miniature repositories for each
   audited rule: forbidden off-target phrases, unhedged cleavage causality,
   direct causality claims that a document-wide disclaimer cannot excuse,
   exact-sentence negation guarding (a negation in another sentence of
   the paragraph does not excuse the match — Codex Sol Round 2),
   adjective-form "Cas9-induced" claims, affirmative causality-
   establishment claims that a document-wide disclaimer cannot excuse,
   suffix negations that no longer excuse affirmative claims
   ("... without ambiguity" / "..., not by contamination" — Codex Sol
   Round 3) together with the immediate-prefix negation forms that
   still do, affirmative "We claim that ..." framing that is itself a
   causality assertion and must be flagged, claim-framing negations in
   an earlier clause separated from the assertion by a contrastive
   barrier ("..., but ...") that no longer excuse it, and the governing
   negative framings ("no claim that" / "no evidence that" / "do not
   claim that" / "never claim that") that still do (Codex Sol Round 4),
   Rule 6 / Rule 8 wording, tampered VALIDATION_REGISTRY.tsv
   (including non-literal truthy flag values such as "yes"), tampered
   WAVE11E_HOLD.md, and broken evidence-chain TSVs (fabricated event
   ids, SITE_UNKNOWN, blank site_id cells, unresolved / orphan sites
   and events, missing sibling tables, empty 0-byte product tables,
   non-numeric / negative / out-of-range / blank coordinates in
   ``mutation_position`` and ``side2_position`` (Codex Sol Round 3),
   the side2 second-side contract (blank permitted on non-JC rows,
   required on JC rows), absurd coordinates, placeholder PAM,
   unvalidated scores, ragged rows, dirty manifest sentinels).
3. Fail-closed guarantees: zero-document and zero-table scans are
   themselves findings; oversize in-scope markdown/metadata files are
   flagged rather than silently skipped.
4. CLI execution via subprocess: exit code 0 on a clean repository, exit
   code 1 on a tampered one, markdown table rendering, ``--repo-root``
   default, ``--output`` report writing, and exit code 2 on a bad root.

Run with: python3 -m pytest -v scripts/test_audit_release_gates.py
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).with_name("audit_release_gates.py")
REPO_ROOT = Path(__file__).resolve().parents[1]

_SPEC = importlib.util.spec_from_file_location("audit_release_gates", SCRIPT_PATH)
if _SPEC is None or _SPEC.loader is None:  # pragma: no cover - load guard
    raise RuntimeError(f"cannot load module from {SCRIPT_PATH}")
MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(MODULE)

VALID_DV1 = "DV1_0123456789abcdef0123456789abcdef"
OTHER_DV1 = "DV1_fedcba9876543210fedcba9876543210"
GHOST_DV1 = "DV1_deadbeefdeadbeefdeadbeefdeadbeef"
REQUIRED_DISCLAIMER = MODULE.REQUIRED_DISCLAIMER

REGISTRY_HEADER = (
    "dataset\tsample_id\ttruth_level\tvariant_validated\tintended_validated"
    "\tmob_validated\tjc_validated\treport_reviewed\tstatus\tnotes"
)
REGISTRY_BL21_ROW = (
    "bl21_user\tB21_3_1\tMANUALLY_CURATED\ttrue\ttrue\tfalse\tfalse\ttrue"
    "\tACTIVE\tsynthetic fixture row"
)
REGISTRY_PRJNA1088182_ROW = (
    "prjna1088182\tcohort_25\tPUBLICATION_REPORTED\tfalse\tfalse\tfalse"
    "\tfalse\tfalse\tPENDING\tsynthetic fixture row"
)
REGISTRY_PRJNA884016_ROW = (
    "prjna884016\tSY22A1\tASSEMBLY_CONFIRMED\tfalse\tfalse\tfalse\tfalse"
    "\tfalse\tPENDING\tsynthetic fixture row"
)

HOLD_BODY = """# WAVE11E hold (fixture)

**Decision:** do not start `WAVE11E_PRODUCTION_AS`. Do not submit a G1 Slurm job.

## Preconditions before any production edit

6. A gate comparison against job 3066886. Current baseline is job 3066886:
   matched 10, OVER 24, UNDER 8.

G1 remains `PROJECT_ACCEPTANCE_PENDING`.
"""

V2_POST_EDIT_VARIANTS = "\n".join(
    [
        "variant_id\tseq_id\tposition\tend\tgd_type\tref\talt\torigin_status"
        "\tintended_relation\tsize_class\treview_priority\tguide_relation"
        "\tguide_mismatches\tpam\tdist_to_site\tmobile_element\tevidence"
        "\tgene\tlocus_tag\tfeature_type\tlegacy_class\tevent_id",
        f"VAR_0001\tchr1\t100\t100\tSNP\tG\tA\tPOST_EDIT_DIFF\tNONE\tSMALL"
        f"\tREVIEW\tNONE\tNA\tNA\tNA\tNONE\tRA=1;MC=NA;JC=NA\tdnaA\tTAG_01"
        f"\tCDS\tscattered_snv\t{VALID_DV1}",
        f"VAR_0002\tchr1\t300\t300\tJC\tNA\tNA\tPOST_EDIT_DIFF\tNONE\tLARGE"
        f"\tREVIEW\tNONE\tNA\tNA\tNA\tNONE\tRA=NA;MC=NA;JC=1\tgeneX\tTAG_02"
        f"\tCDS\tstructural\t{OTHER_DV1}",
    ]
) + "\n"

V2_UNINTENDED = "\n".join(
    [
        "\t".join(
            [
                "seq_id",
                "position",
                "end",
                "gd_type",
                "ref",
                "alt",
                "class",
                "editor",
                "pam_profile",
                "offtarget_mismatch",
                "distance_to_site",
                "side2_seq_id",
                "side2_position",
                "event_id",
            ]
        ),
        "\t".join(
            [
                "chr1",
                "100",
                "100",
                "SNP",
                "G",
                "A",
                "scattered_snv",
                "cas9",
                "",
                "",
                "",
                "",
                "",
                VALID_DV1,
            ]
        ),
    ]
) + "\n"

V2_LINKS = "\n".join(
    [
        "mutation_id\tsite_id\tmutation_type\tmutation_position\tsite_start"
        "\tsite_end\tdistance_to_site\tmismatches\tpam\tcfd_score"
        "\tassociation_window\tevent_id\tjunction_side\tsearch_backend"
        "\tbulge_type\tbulge_size",
        f"mut_1\tSITE_000001\tSNP\t100\t3000\t3022\t8\t0\tCGG\tNA\t50"
        f"\t{VALID_DV1}\tNA\trust_exact\tnone\t0",
    ]
) + "\n"

V2_SITES = "\n".join(
    [
        "site_id\tseq_id\tstart\tend\tstrand\ttarget_seq\tpam\tmismatches"
        "\tbulge_type\tbulge_size\tsearch_backend\tcfd_score\thsu_score",
        "SITE_000001\tchr1\t3000\t3022\t+\tGACTGACTGACTGACTGACTTGG\tCGG\t0"
        "\tnone\t0\trust_exact\tNA\tNA",
    ]
) + "\n"

V2_EDIT_OUTCOMES = "\n".join(
    [
        "edit_id\tkind\tseq_id\texpected_start\texpected_end\tstatus"
        "\tmatched_event_ids\tleft_boundary_status\tright_boundary_status"
        "\texpected_size\tobserved_size\tunexpected_events\tnotes"
        "\tmatched_event_ids_dv1\tunexpected_event_ids_dv1",
        f"edit_1\tdel\tchr1\t500\t900\tCOMPLETE\t1\tMATCH\tMATCH\t401\t401"
        f"\tNONE\tNONE\t{VALID_DV1}\tNONE",
    ]
) + "\n"

V1_POST_EDIT_VARIANTS = "\n".join(
    [
        "variant_id\tseq_id\tposition\tend\tgd_type\tref\talt\torigin_status"
        "\tintended_relation\tsize_class\treview_priority\tguide_relation"
        "\tguide_mismatches\tpam\tdist_to_site\tmobile_element\tevidence"
        "\tgene\tlocus_tag\tfeature_type\tlegacy_class",
        "VAR_0001\tchr1\t100\t100\tSNP\tG\tA\tPOST_EDIT_DIFF\tNONE\tSMALL"
        "\tREVIEW\tNONE\tNA\tNA\tNA\tNONE\tRA=1;MC=NA;JC=NA\tdnaA\tTAG_01"
        "\tCDS\tscattered_snv",
    ]
) + "\n"


def _unintended_with_side2(side2_snp: str, side2_jc: str) -> str:
    """``unintended.tsv`` carrying one SNP row and one JC row so the
    ``side2_position`` second-side contract can be exercised (Codex Sol
    Round 3).  Column order mirrors ``V2_UNINTENDED``."""
    header = V2_UNINTENDED.splitlines()[0]
    row_snp = "\t".join(
        [
            "chr1", "100", "100", "SNP", "G", "A", "scattered_snv", "cas9",
            "", "", "", "", side2_snp, VALID_DV1,
        ]
    )
    row_jc = "\t".join(
        [
            "chr1", "300", "300", "JC", "NA", "NA", "structural", "cas9",
            "", "", "", "chr1", side2_jc, OTHER_DV1,
        ]
    )
    return header + "\n" + row_snp + "\n" + row_jc + "\n"


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _write(root: Path, relpath: str, content: str) -> Path:
    path = root / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def _make_registry(
    root: Path,
    *,
    prjna1088182_status: str = "PENDING",
    prjna884016_status: str = "PENDING",
    prjna1088182_flags: tuple[str, ...] = ("false", "false", "false", "false"),
    prjna884016_flags: tuple[str, ...] = ("false", "false", "false", "false"),
    include_prjna1088182: bool = True,
    include_prjna884016: bool = True,
) -> Path:
    lines = [REGISTRY_HEADER, REGISTRY_BL21_ROW]
    if include_prjna1088182:
        lines.append(
            "\t".join(
                [
                    "prjna1088182",
                    "cohort_25",
                    "PUBLICATION_REPORTED",
                    *prjna1088182_flags,
                    "false",
                    prjna1088182_status,
                    "synthetic fixture row",
                ]
            )
        )
    if include_prjna884016:
        lines.append(
            "\t".join(
                [
                    "prjna884016",
                    "SY22A1",
                    "ASSEMBLY_CONFIRMED",
                    *prjna884016_flags,
                    "false",
                    prjna884016_status,
                    "synthetic fixture row",
                ]
            )
        )
    return _write(
        root, "benchmark/real_world/VALIDATION_REGISTRY.tsv", "\n".join(lines) + "\n"
    )


def _make_hold(root: Path, body: str = HOLD_BODY) -> Path:
    return _write(
        root,
        "benchmark/results/parity_closure/g1_clonal_jc/review/WAVE11E_HOLD.md",
        body,
    )


def _make_v2_product_set(root: Path, subdir: str = "validation_1/BL21_G1") -> Path:
    base = f"benchmark/real_world/results/bl21_user/m5_p3/{subdir}"
    _write(root, f"{base}/post_edit_variants.tsv", V2_POST_EDIT_VARIANTS)
    _write(root, f"{base}/unintended.tsv", V2_UNINTENDED)
    _write(root, f"{base}/mutation_offtarget_links.tsv", V2_LINKS)
    _write(root, f"{base}/offtarget_sites.tsv", V2_SITES)
    _write(root, f"{base}/edit_outcomes.tsv", V2_EDIT_OUTCOMES)
    return root / base


def _make_clean_repo(tmp_path: Path) -> Path:
    """A miniature repository that passes every release gate."""
    root = tmp_path / "repo"
    root.mkdir()
    _write(
        root,
        "docs/notes.md",
        "All findings are candidate guide-dependent off-target events.\n"
        f"{REQUIRED_DISCLAIMER}\n",
    )
    _make_registry(root)
    _make_hold(root)
    _make_v2_product_set(root)
    _write(
        root,
        "benchmark/real_world/manifests/recovery.tsv",
        "sample\tdv1_event_id\tevent_id_if_applicable\n"
        f"Cas9_1\tNOT_EVALUABLE\tNA\n",
    )
    return root


# ---------------------------------------------------------------------------
# Gate 1: scientific language
# ---------------------------------------------------------------------------


class TestScientificLanguage:
    def test_real_repo_passes(self) -> None:
        assert MODULE.check_scientific_language(REPO_ROOT) == []

    @pytest.mark.parametrize(
        "phrase",
        [
            "This locus is a confirmed off-target event.",
            "This locus is a confirmed offtarget event.",
            "This locus is a confirmed_off_target event.",
            "This locus is a CONFIRMED OFF-TARGET event.",
        ],
    )
    def test_confirmed_offtarget_flagged(self, tmp_path: Path, phrase: str) -> None:
        root = tmp_path / "repo"
        _write(root, "docs/report.md", f"# Report\n\n{phrase}\n")
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "forbidden phrase" in findings[0]
        assert "docs/report.md" in findings[0]

    def test_confirmed_offtarget_under_benchmark_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "benchmark/real_world/NOTES.md",
            "The cohort shows one confirmed off-target deletion.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "benchmark/real_world/NOTES.md" in findings[0]

    def test_archives_are_ignored(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(root, "docs/keep.md", "Benign curated note.\n")
        _write(
            root,
            "docs/development/archive/legacy.md",
            "This is a confirmed off-target event.\n",
        )
        _write(
            root,
            "benchmark/results/parity_closure/state/lock_archive_20261103/held.md",
            "This is a confirmed off-target event.\n",
        )
        assert MODULE.check_scientific_language(root) == []

    def test_negative_occurrence_in_scientific_model_allowed(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "docs/scientific_model.md",
            "术语规范：统一使用 candidate guide-dependent off-target，"
            "严禁使用 `confirmed off-target`。\n"
            "English note: never write 'confirmed off-target' in a report.\n"
            f"Boundary: {REQUIRED_DISCLAIMER}\n",
        )
        assert MODULE.check_scientific_language(root) == []

    def test_assertive_occurrence_in_scientific_model_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "docs/scientific_model.md",
            "Variant 1 is a confirmed off-target event.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1

    def test_unhedged_causality_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "The 19 kb deletion was caused by Cas9 cleavage at the locus.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "cleavage-causality assertion" in findings[0]
        assert "does not establish cleavage causality" in findings[0]

    def test_direct_causality_despite_document_disclaimer_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol counterexample: a disclaimer elsewhere in the document
        must not whitewash a direct causality claim in another sentence."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "The 19 kb deletion was caused by Cas9.\n\n"
            f"## Limitations\n\n{REQUIRED_DISCLAIMER}\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "direct causality claim" in findings[0]
        assert "does not excuse" in findings[0]

    @pytest.mark.parametrize(
        "sentence",
        [
            "The deletion was not caused by Cas9; the evidence does not "
            "establish cleavage causality.",
            "No claim is made that Cas9 caused this deletion.",
            "We cannot prove the deletion was induced by Cas9.",
        ],
    )
    def test_negated_direct_causality_allowed(
        self, tmp_path: Path, sentence: str
    ) -> None:
        root = tmp_path / "repo"
        _write(root, "docs/report.md", f"{sentence}\n")
        assert MODULE.check_scientific_language(root) == []

    def test_prohibition_paragraph_context_allowed(self, tmp_path: Path) -> None:
        """Claim-framing guards in the enclosing paragraph (e.g. the real
        M5_CLOSEOUT wording) are exempt from direct-causality policing."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/closeout.md",
            "Out of scope for this closure: CAST/IS110, long reads,\n"
            "and any claim that Cas9 caused the BL21 structural events.\n",
        )
        assert MODULE.check_scientific_language(root) == []

    def test_cas_induced_adjective_form_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 2: 'The deletion is Cas9-induced.' is a direct
        causality claim and must be flagged."""
        root = tmp_path / "repo"
        _write(root, "docs/report.md", "The deletion is Cas9-induced.\n")
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "Cas9-induced" in findings[0]
        assert "direct causality claim" in findings[0]

    @pytest.mark.parametrize("editor", ["Cas9", "Cas12a"])
    def test_cas_induced_with_space_or_hyphen_flagged(
        self, tmp_path: Path, editor: str
    ) -> None:
        """The '-induced' pattern tolerates both hyphen and space forms."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            f"The rearrangement is {editor} induced at the locus.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert f"{editor} induced" in findings[0]

    def test_negation_in_other_sentence_does_not_excuse(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 2: 'No contamination was detected. The deletion
        was caused by Cas9.' must flag — the guard must sit in the match's
        own sentence, not elsewhere in the paragraph."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "No contamination was detected. The deletion was caused by "
            "Cas9.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "was caused by Cas9" in findings[0]
        assert "direct causality claim" in findings[0]

    def test_cas_induced_despite_document_disclaimer_flagged(
        self, tmp_path: Path
    ) -> None:
        """A document-wide disclaimer does not whitewash an adjective-form
        direct causality claim either."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "The deletion is Cas9-induced.\n\n"
            f"## Limitations\n\n{REQUIRED_DISCLAIMER}\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "Cas9-induced" in findings[0]

    def test_causality_establishment_despite_document_disclaimer_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 2: affirmative establishment claims are no
        longer excusable by a document-wide disclaimer; only a negation in
        the claim's own sentence ('does not establish') excuses them."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "Section 5 establishes cleavage causality for the proximal event "
            "in the hypothesis tier.\n\n"
            f"Limitation: {REQUIRED_DISCLAIMER}\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "establishes cleavage causality" in findings[0]
        assert "document-wide disclaimer does not excuse" in findings[0]

    def test_establishment_claim_with_document_disclaimer_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 2 counterexample: disclaimer + 'The evidence
        establishes cleavage causality.' must flag."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            f"Limitation: {REQUIRED_DISCLAIMER}\n\n"
            "The evidence establishes cleavage causality.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "establishes cleavage causality" in findings[0]
        assert "establishment claim" in findings[0]

    def test_establishment_suffix_without_phrase_does_not_excuse(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 3 counterexample: 'The evidence establishes
        cleavage causality without ambiguity.' must flag — a 'without'
        phrase in the suffix never negates the affirmative claim it
        follows."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "The evidence establishes cleavage causality without "
            "ambiguity.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "establishes cleavage causality" in findings[0]
        assert "establishment claim" in findings[0]
        assert "suffix" in findings[0]

    def test_direct_causality_suffix_not_by_does_not_excuse(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 3 counterexample: 'The deletion was caused by
        Cas9, not by contamination.' must flag — a negation in the
        suffix after the claim never excuses it."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "The deletion was caused by Cas9, not by contamination.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "was caused by Cas9" in findings[0]
        assert "direct causality claim" in findings[0]
        assert "suffix" in findings[0]

    @pytest.mark.parametrize(
        "sentence",
        [
            "The deletion was not caused by Cas9.",
            "The rearrangement is not Cas9-induced.",
            "Section 5 does not establish cleavage causality.",
            "We cannot establish cleavage causality at this locus.",
            "This analysis fails to establish cleavage causality.",
            "The report never establishes cleavage causality.",
            "The event is non-Cas9-induced by our criteria.",
            "The junction was not directly caused by Cas9.",
        ],
    )
    def test_immediate_prefix_negation_forms_allowed(
        self, tmp_path: Path, sentence: str
    ) -> None:
        """Codex Sol Round 3: verbal/adjectival negation directly
        preceding the claim inside its own clause still exempts the
        match (auxiliary + negator, bare negator, 'fails to' /
        'never' / 'non-' prefixed adjective, adverb gaps)."""
        root = tmp_path / "repo"
        _write(root, "docs/report.md", f"{sentence}\n")
        assert MODULE.check_scientific_language(root) == []

    def test_negated_establishment_sentence_allowed(
        self, tmp_path: Path
    ) -> None:
        """The required disclaimer sentence itself carries the negation
        ('does not establish cleavage causality') and must pass."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            f"{REQUIRED_DISCLAIMER}\n",
        )
        assert MODULE.check_scientific_language(root) == []

    # --- Codex Sol Round 4 ---------------------------------------------
    # Affirmative "claim that" framing is not negation, and a framing
    # negation separated from the assertion by a contrastive barrier
    # governed the prior clause, not the assertion.

    @pytest.mark.parametrize(
        ("sentence", "expected_fragment"),
        [
            (
                "We claim that the deletion was caused by Cas9.",
                "direct causality claim",
            ),
            (
                "We claim that this evidence establishes cleavage "
                "causality.",
                "establishment claim",
            ),
        ],
    )
    def test_affirmative_claim_that_framing_flagged(
        self, tmp_path: Path, sentence: str, expected_fragment: str
    ) -> None:
        """Codex Sol Round 4 counterexamples: an affirmative claim frame
        ('We claim that ...') asserts the causality rather than negating
        it — it must be flagged, not excused as framing negation."""
        root = tmp_path / "repo"
        _write(root, "docs/report.md", f"{sentence}\n")
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "docs/report.md" in findings[0]
        assert expected_fragment in findings[0]

    @pytest.mark.parametrize(
        ("sentence", "expected_fragment"),
        [
            (
                "There is no evidence of contamination, but the deletion "
                "was caused by Cas9.",
                "direct causality claim",
            ),
            (
                "We cannot prove contamination, but this evidence "
                "establishes cleavage causality.",
                "establishment claim",
            ),
        ],
    )
    def test_contrast_barrier_blocks_framing_excuse_flagged(
        self, tmp_path: Path, sentence: str, expected_fragment: str
    ) -> None:
        """Codex Sol Round 4 counterexamples: a framing negation in an
        earlier clause ('There is no evidence of contamination') does
        NOT govern the causal assertion after a contrastive barrier
        (', but ...') and must not excuse it."""
        root = tmp_path / "repo"
        _write(root, "docs/report.md", f"{sentence}\n")
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "docs/report.md" in findings[0]
        assert expected_fragment in findings[0]

    def test_contrast_word_barrier_without_comma_flagged(
        self, tmp_path: Path
    ) -> None:
        """The contrastive barrier does not need a comma: a bare
        contrastive conjunction between the framing and the assertion
        blocks the excuse just the same."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "There is no evidence that contamination occurred but the "
            "deletion was caused by Cas9.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "was caused by Cas9" in findings[0]
        assert "direct causality claim" in findings[0]

    @pytest.mark.parametrize(
        "sentence",
        [
            "We make no claim that the deletion was caused by Cas9.",
            "There is no evidence that the deletion was caused by Cas9.",
            "We do not claim that this evidence establishes cleavage "
            "causality.",
            "We never claim that the deletion was caused by Cas9.",
        ],
    )
    def test_governing_negative_claim_framing_still_allowed(
        self, tmp_path: Path, sentence: str
    ) -> None:
        """Round 4 tightens the frame without banning genuine negative
        framing: 'no claim that' / 'no evidence that' / 'do not claim
        that' / 'never claim that' reaching the assertion with no
        contrastive barrier still exempt the match."""
        root = tmp_path / "repo"
        _write(root, "docs/report.md", f"{sentence}\n")
        assert MODULE.check_scientific_language(root) == []

    def test_immediate_negation_earlier_clause_does_not_excuse(
        self, tmp_path: Path
    ) -> None:
        """Round 4: the immediate verbal/adjectival negation is scoped
        to the immediate comma/colon-delimited clause — a negated
        earlier clause cannot whitewash a later affirmative clause."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "It is not established: the deletion was caused by Cas9.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "was caused by Cas9" in findings[0]
        assert "direct causality claim" in findings[0]

    def test_infinitive_establish_causality_flagged(
        self, tmp_path: Path
    ) -> None:
        """Bare infinitive 'establish cleavage causality' without negation is flagged."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "This study will establish cleavage causality across all target sites.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "establish cleavage causality" in findings[0]
        assert "establishment claim" in findings[0]

    def test_infinitive_establish_causality_negated_allowed(
        self, tmp_path: Path
    ) -> None:
        """Bare infinitive 'establish cleavage causality' with immediate negation passes."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "This study will not establish cleavage causality across these sites.\n",
        )
        assert MODULE.check_scientific_language(root) == []

    # --- Codex Sol Round 5 ---------------------------------------------
    # 1. Punctuation barriers (colon, Chinese comma) block clause crossing.
    # 2. Affirmative endorsements of 'any claim' are flagged.

    @pytest.mark.parametrize(
        ("sentence", "expected_fragment"),
        [
            (
                "There is no evidence of contamination: this evidence establishes cleavage causality.",
                "establishment claim",
            ),
            (
                "There is no evidence of contamination，the deletion was caused by Cas9.",
                "direct causality claim",
            ),
            (
                "We support any claim that the deletion was caused by Cas9.",
                "direct causality claim",
            ),
            (
                "We endorse any claim that this evidence establishes cleavage causality.",
                "establishment claim",
            ),
        ],
    )
    def test_round5_causality_bypass_probes_flagged(
        self, tmp_path: Path, sentence: str, expected_fragment: str
    ) -> None:
        """Codex Sol Round 5 counterexamples: colon/CJK punctuation blocks
        framing negation across clauses, and affirmative endorsements of 'any claim'
        must be flagged."""
        root = tmp_path / "repo"
        _write(root, "docs/report.md", f"{sentence}\n")
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "docs/report.md" in findings[0]
        assert expected_fragment in findings[0]

    # --- Codex Sol Round 6 ---------------------------------------------
    # Unrelated scope prohibitions, cross-sentence leaks, semicolon clause boundaries,
    # and assertion verbs must not govern 'any claim that'.

    @pytest.mark.parametrize(
        ("sentence", "expected_fragment"),
        [
            (
                "Contamination is out of scope. We confirm any claim that this evidence establishes cleavage causality.",
                "establishment claim",
            ),
            (
                "Out of scope: contamination; we assert any claim that the deletion was caused by Cas9.",
                "direct causality claim",
            ),
            (
                "Out of scope: contamination, we assert any claim that the deletion was caused by Cas9.",
                "direct causality claim",
            ),
            (
                "Out of scope: contamination, the team accepts any claim that the deletion was caused by Cas9.",
                "direct causality claim",
            ),
            (
                "Out of scope: contamination, the team accepts any claim that this evidence establishes cleavage causality.",
                "establishment claim",
            ),
            (
                "Out of scope: contamination, our position accepts any claim that the deletion was caused by Cas9.",
                "direct causality claim",
            ),
            (
                "超出范围：污染。我们确认任何声称 the deletion was caused by Cas9 的主张。",
                "direct causality claim",
            ),
            (
                "超出范围：污染；我们坚持任何声称 the deletion was caused by Cas9 的主张。",
                "direct causality claim",
            ),
            (
                "Out of scope: any claim that contamination exists, the deletion was caused by Cas9.",
                "direct causality claim",
            ),
            (
                "Out of scope: any claim that contamination exists, this evidence establishes cleavage causality.",
                "establishment claim",
            ),
            (
                "超出范围：任何声称存在污染，the deletion was caused by Cas9。",
                "direct causality claim",
            ),
            (
                "超出范围：任何声称存在污染，this evidence establishes cleavage causality。",
                "establishment claim",
            ),
        ],
    )
    def test_round6_scope_prohibition_bypass_probes_flagged(
        self, tmp_path: Path, sentence: str, expected_fragment: str
    ) -> None:
        """Codex Sol Round 6 & 7 counterexamples: unrelated scope prohibitions across
        sentences or semicolons, comma-separated clauses, and assertion verbs in lists must be flagged."""
        root = tmp_path / "repo"
        _write(root, "docs/report.md", f"{sentence}\n")
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "docs/report.md" in findings[0]
        assert expected_fragment in findings[0]

    def test_round7_wrapped_comma_clause_bypass_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 7: scope prohibitions wrapped after comma onto another line
        must not excuse subsequent affirmative assertions."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "Out of scope: any claim that contamination exists,\n"
            "the deletion was caused by Cas9.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "direct causality claim" in findings[0]

        _write(
            root,
            "docs/report2.md",
            "Out of scope: any claim that contamination exists,\n"
            "this evidence establishes cleavage causality.\n",
        )
        findings2 = MODULE.check_scientific_language(root)
        assert len(findings2) == 2

    def test_round6_legitimate_scope_prohibitions_pass(
        self, tmp_path: Path
    ) -> None:
        """Legitimate colon-governed scope prohibitions in English and Chinese pass."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "Out of scope for this closure: CAST/IS110, long reads, and any claim that the deletion was caused by Cas9.\n"
            "超出范围：PRJNA1088182、长读长测序，以及任何声称 the deletion was caused by Cas9 的结论。\n"
            "Out of scope: any claim that the deletion was caused by Cas9.\n"
            "Out of scope: any claim that contamination exists, or any claim that the deletion was caused by Cas9.\n"
            "超出范围：任何声称存在污染，或任何声称 the deletion was caused by Cas9 的主张。\n",
        )
        assert MODULE.check_scientific_language(root) == []

    def test_hedged_candidate_language_passes(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "docs/report.md",
            "Finding: candidate guide-dependent off-target event near SITE_000001.\n"
            f"{REQUIRED_DISCLAIMER}\n",
        )
        assert MODULE.check_scientific_language(root) == []

    def test_guide_dependent_discussion_requires_disclaimer(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol P1: documents discussing candidate guide-dependent
        associations must carry the required disclaimer."""
        root = tmp_path / "repo"
        _write(
            root,
            "docs/review.md",
            "# Review\n\nWe observed several candidate guide-dependent "
            "off-target events near the target locus.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "candidate guide-dependent" in findings[0]
        assert "lacks the required disclaimer" in findings[0]

    def test_generated_product_reports_exempt_from_disclaimer_requirement(
        self, tmp_path: Path
    ) -> None:
        """Engine-generated reports under benchmark/**/results/ are worded
        by the frozen engine template and are exempt from the Tier-C
        disclaimer requirement (direct-causality policing still applies)."""
        root = tmp_path / "repo"
        _write(
            root,
            "benchmark/results/m9_fixture_0001/prokadiff/report.md",
            "## 6. Candidate Guide-dependent Off-target Events\n\n"
            "No candidate guide-dependent off-target events were identified.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert findings == []

    def test_rule6_all_ncbi_differences_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "docs/legacy.md",
            "We report all differences vs NCBI as unintended edits.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "Rule 6" in findings[0]

    def test_rule6_unintended_vs_ncbi_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "benchmark/real_world/legacy_report.md",
            "The table lists unintended mutations relative to NCBI.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "Rule 6" in findings[0]

    def test_rule6_negated_prohibition_allowed(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "docs/policy.md",
            "Never report all differences vs NCBI as unintended edits.\n",
        )
        assert MODULE.check_scientific_language(root) == []

    def test_rule8_cast_claim(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "docs/roadmap.md",
            "CAST editors show no genome-wide mutagenesis, so screening is "
            "unnecessary.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "Rule 8" in findings[0]

    def test_rule8_is110_claim(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "benchmark/real_world/editors.md",
            "`is110` does not cause any genome-wide mutagenesis.\n",
        )
        findings = MODULE.check_scientific_language(root)
        assert len(findings) == 1
        assert "Rule 8" in findings[0]

    def test_rule8_roadmap_mention_allowed(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(
            root,
            "docs/schema_notes.md",
            "CAST/IS110 support is deferred to the roadmap; `cast` / `is110` "
            "error and point to STAGE_ROADMAP.md.\n",
        )
        assert MODULE.check_scientific_language(root) == []


# ---------------------------------------------------------------------------
# Gate 1: fail-closed scan guarantees
# ---------------------------------------------------------------------------


class TestScientificLanguageFailClosed:
    def test_zero_documents_fail_closed(self, tmp_path: Path) -> None:
        """Codex Sol counterexample: an empty scan scope is itself a
        failure, never a clean pass."""
        root = tmp_path / "repo"
        root.mkdir()
        (root / "docs").mkdir()
        (root / "benchmark").mkdir()
        findings = MODULE.check_scientific_language(root)
        assert any("zero scannable documents" in f for f in findings)

    def test_missing_trees_fail_closed(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        root.mkdir()
        findings = MODULE.check_scientific_language(root)
        assert any("zero scannable documents" in f for f in findings)

    def test_oversize_markdown_flagged_not_skipped(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _write(root, "docs/small.md", "Benign note.\n")
        huge = root / "docs" / "huge.md"
        huge.write_text("# huge\n" + ("x" * (MODULE.SCAN_SIZE_LIMIT_BYTES + 1)))
        findings = MODULE.check_scientific_language(root)
        assert any(
            "huge.md" in f and "SCAN_SIZE_LIMIT_BYTES" in f for f in findings
        )
        assert not any("forbidden" in f for f in findings)

    def test_oversize_metadata_table_flagged_not_skipped(
        self, tmp_path: Path
    ) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        huge = root / "benchmark/real_world/manifests/huge_manifest.tsv"
        huge.parent.mkdir(parents=True, exist_ok=True)
        huge.write_text("sample\tdv1_event_id\n" + ("x" * (MODULE.SCAN_SIZE_LIMIT_BYTES + 1)))
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "huge_manifest.tsv" in f and "SCAN_SIZE_LIMIT_BYTES" in f
            for f in findings
        )

    def test_oversize_bulk_artifact_table_remains_out_of_scope(
        self, tmp_path: Path
    ) -> None:
        """Raw engine bulk tables (e.g. repeat_ambiguous_junctions.tsv)
        above the size bound stay out of scope — documented exemption."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        bulk = root / "benchmark/results/m6/repeat_ambiguous_junctions.tsv"
        bulk.parent.mkdir(parents=True, exist_ok=True)
        bulk.write_text("junction\n" + ("x" * (MODULE.SCAN_SIZE_LIMIT_BYTES + 1)))
        assert MODULE.check_evidence_chain_traceability(root) == []


# ---------------------------------------------------------------------------
# Gate 2: evidence-chain traceability
# ---------------------------------------------------------------------------


class TestEvidenceChain:
    def test_real_repo_passes(self) -> None:
        assert MODULE.check_evidence_chain_traceability(REPO_ROOT) == []

    def test_clean_v2_product_set_passes(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        assert MODULE.check_evidence_chain_traceability(root) == []

    def test_v1_product_set_locus_reference_passes(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        base = "benchmark/real_world/results/bl21_user/B21_3_1_diff"
        _write(root, f"{base}/post_edit_variants.tsv", V1_POST_EDIT_VARIANTS)
        assert MODULE.check_evidence_chain_traceability(root) == []

    @pytest.mark.parametrize(
        "bad_event_id",
        [
            "EVENT_1",
            "DV1_0123",  # too short
            "DV1_ABCDEF0123456789ABCDEF0123456789",  # uppercase hex
            "VAR_0001",  # display label reused as identity
            "mut_3",  # legacy identifier reused as identity
        ],
    )
    def test_malformed_event_id_flagged(self, tmp_path: Path, bad_event_id: str) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "post_edit_variants.tsv"
        )
        rows = target.read_text().splitlines()
        rows[1] = rows[1].rsplit("\t", 1)[0] + f"\t{bad_event_id}"
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("post_edit_variants row" in f and bad_event_id in f for f in findings)

    def test_empty_event_id_cell_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "post_edit_variants.tsv"
        )
        rows = target.read_text().splitlines()
        rows[2] = rows[2].rsplit("\t", 1)[0] + "\t"
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("empty event_id cell" in f for f in findings)

    def test_na_event_id_with_valid_locus_passes(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "post_edit_variants.tsv"
        )
        rows = target.read_text().splitlines()
        rows[2] = rows[2].rsplit("\t", 1)[0] + "\tNA"
        target.write_text("\n".join(rows) + "\n")
        assert MODULE.check_evidence_chain_traceability(root) == []

    def test_site_unknown_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "mutation_offtarget_links.tsv"
        )
        rows = target.read_text().splitlines()
        rows[1] = rows[1].replace("SITE_000001", "SITE_UNKNOWN")
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("SITE_UNKNOWN" in f for f in findings)

    def test_unresolved_site_reference_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "mutation_offtarget_links.tsv"
        )
        rows = target.read_text().splitlines()
        rows[1] = rows[1].replace("SITE_000001", "SITE_999999")
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("does not resolve" in f for f in findings)

    def test_orphan_site_with_missing_sibling_table_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol counterexample: a site reference with no sibling
        ``offtarget_sites.tsv`` in the same directory is an orphan."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        (root / "benchmark/real_world/results/bl21_user/m5_p3/validation_1"
         / "BL21_G1" / "offtarget_sites.tsv").unlink()
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "SITE_000001" in f and "orphan reference" in f and
            "offtarget_sites.tsv" in f
            for f in findings
        )

    def test_orphan_event_with_missing_sibling_table_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol counterexample: a DV1 event reference with no sibling
        ``post_edit_variants.tsv`` in the same directory is an orphan."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        (root / "benchmark/real_world/results/bl21_user/m5_p3/validation_1"
         / "BL21_G1" / "post_edit_variants.tsv").unlink()
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            VALID_DV1 in f and "orphan reference" in f and
            "post_edit_variants.tsv" in f
            for f in findings
        )

    def test_orphan_event_id_not_in_sibling_flagged(self, tmp_path: Path) -> None:
        """A well-formed DV1 id that the sibling variants table does not
        declare is an orphan event reference."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "unintended.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("event_id")] = GHOST_DV1
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "does not resolve" in f and GHOST_DV1 in f for f in findings
        )

    def test_mechanistic_associations_orphan_site_flagged(
        self, tmp_path: Path
    ) -> None:
        """The generic product-table rule (e.g.
        ``mechanistic_associations.tsv``) polices site/event references."""
        root = tmp_path / "repo"
        base = _make_v2_product_set(root)
        _write(
            root,
            f"{base.relative_to(root).as_posix()}/mechanistic_associations.tsv",
            f"association_id\tsite_id\tevent_id\trationale\n"
            f"assoc_1\tSITE_777777\t{VALID_DV1}\tspatial\n",
        )
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "mechanistic_associations.tsv" in f and "SITE_777777" in f
            for f in findings
        )

    def test_mechanistic_associations_clean_refs_pass(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        base = _make_v2_product_set(root)
        _write(
            root,
            f"{base.relative_to(root).as_posix()}/mechanistic_associations.tsv",
            "association_id\tsite_id\tevent_id\trationale\n"
            f"assoc_1\tSITE_000001\t{VALID_DV1}\tspatial\n",
        )
        assert MODULE.check_evidence_chain_traceability(root) == []

    @pytest.mark.parametrize(
        ("column", "bad_value"),
        [
            ("position", "999999999"),
            ("end", "2000000000"),
        ],
    )
    def test_absurd_coordinate_flagged(
        self, tmp_path: Path, column: str, bad_value: str
    ) -> None:
        """Codex Sol P1: coordinates beyond bacterial-contig sanity bounds
        (pos <= 10,000,000) are flagged as absurd."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "post_edit_variants.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index(column)] = bad_value
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("sanity bound" in f or "absurd coordinate" in f for f in findings)

    def test_ecoli_scale_coordinate_passes(self, tmp_path: Path) -> None:
        """A 4.6 Mb E. coli-scale coordinate is perfectly sane."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        base = root / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        target = base / "post_edit_variants.tsv"
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("position")] = "4600000"
        row[header.index("end")] = "4600000"
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")

        # Also update sibling tables referencing DEV1 so coordinates remain reconciled
        unintended = base / "unintended.tsv"
        u_rows = unintended.read_text().splitlines()
        u_header = u_rows[0].split("\t")
        u_row = u_rows[1].split("\t")
        u_row[u_header.index("position")] = "4600000"
        if "end" in u_header:
            u_row[u_header.index("end")] = "4600000"
        u_rows[1] = "\t".join(u_row)
        unintended.write_text("\n".join(u_rows) + "\n")

        links = base / "mutation_offtarget_links.tsv"
        l_rows = links.read_text().splitlines()
        l_header = l_rows[0].split("\t")
        l_row = l_rows[1].split("\t")
        l_row[l_header.index("mutation_position")] = "4600000"
        l_rows[1] = "\t".join(l_row)
        links.write_text("\n".join(l_rows) + "\n")

        assert MODULE.check_evidence_chain_traceability(root) == []

    @pytest.mark.parametrize(
        ("column", "bad_value"),
        [
            ("position", "-1"),  # negative integer
            ("position", "abc"),  # non-numeric
            ("end", "10000001"),  # one past the sanity bound
            ("position", ""),  # blank cell
        ],
    )
    def test_non_numeric_or_invalid_coordinate_flagged(
        self, tmp_path: Path, column: str, bad_value: str
    ) -> None:
        """Codex Sol Round 2: a coordinate cell must be clean 'NA' or a
        positive integer within 1..10,000,000; '-1', 'abc', '10000001'
        and blank cells are errors, never silent skips."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "post_edit_variants.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index(column)] = bad_value
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            f"column '{column}'" in f and f"'{bad_value}'" in f
            for f in findings
        )

    def test_na_coordinate_cell_still_passes(self, tmp_path: Path) -> None:
        """The clean missing literal 'NA' remains accepted in a coordinate
        column (e.g. end for a SNP-scale variant)."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "offtarget_sites.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("end")] = "NA"
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        assert MODULE.check_evidence_chain_traceability(root) == []

    def test_blank_site_id_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 2: a blank site_id cell is an error, never a
        silent pass."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "mutation_offtarget_links.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("site_id")] = ""
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "site_id cell is empty" in f and "mutation_offtarget_links.tsv" in f
            for f in findings
        )

    def test_empty_zero_byte_product_table_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 2: a 0-byte post_edit_variants.tsv (no header
        row) is an EVIDCHAIN failure instead of a silent skip."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "post_edit_variants.tsv"
        )
        target.write_text("")  # 0 bytes
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "post_edit_variants.tsv" in f and "product table is empty" in f
            for f in findings
        )

    def test_mechanistic_associations_coordinates_validated(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 2: the generic product-table branch (e.g.
        mechanistic_associations.tsv) validates coordinate columns too."""
        root = tmp_path / "repo"
        base = _make_v2_product_set(root)
        _write(
            root,
            f"{base.relative_to(root).as_posix()}/mechanistic_associations.tsv",
            "association_id\tsite_id\tevent_id\tsite_start\tsite_end\n"
            f"assoc_1\tSITE_000001\t{VALID_DV1}\t-1\tabc\n",
        )
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "column 'site_start'" in f and "'-1'" in f for f in findings
        )
        assert any(
            "column 'site_end'" in f and "'abc'" in f for f in findings
        )

    @pytest.mark.parametrize(
        "bad_value", ["-1", "abc", "10000001"]
    )
    def test_mutation_position_invalid_flagged(
        self, tmp_path: Path, bad_value: str
    ) -> None:
        """Codex Sol Round 3: mutation_position is a policed coordinate
        column; '-1', 'abc', and '10000001' must trigger findings."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "mutation_offtarget_links.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("mutation_position")] = bad_value
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "column 'mutation_position'" in f and f"'{bad_value}'" in f
            for f in findings
        )

    @pytest.mark.parametrize(
        "bad_value", ["-1", "abc", "10000001"]
    )
    def test_side2_position_invalid_on_jc_flagged(
        self, tmp_path: Path, bad_value: str
    ) -> None:
        """Codex Sol Round 3: a populated side2_position on a JC row
        must be a valid coordinate; '-1', 'abc', and '10000001' trigger
        findings."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        _write(
            root,
            "benchmark/real_world/results/bl21_user/m5_p3/validation_1"
            "/BL21_G1/unintended.tsv",
            _unintended_with_side2("", bad_value),
        )
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "column 'side2_position'" in f and f"'{bad_value}'" in f
            for f in findings
        )

    @pytest.mark.parametrize(
        "bad_value", ["-1", "abc", "10000001"]
    )
    def test_side2_position_invalid_on_non_jc_flagged(
        self, tmp_path: Path, bad_value: str
    ) -> None:
        """Codex Sol Round 3: side2_position populated on a non-JC row
        is validated like any other coordinate; '-1', 'abc', and
        '10000001' trigger findings."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        _write(
            root,
            "benchmark/real_world/results/bl21_user/m5_p3/validation_1"
            "/BL21_G1/unintended.tsv",
            _unintended_with_side2(bad_value, "4000000"),
        )
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "column 'side2_position'" in f and f"'{bad_value}'" in f
            for f in findings
        )

    @pytest.mark.parametrize("blank", ["", ".", "NA"])
    def test_side2_position_blank_on_non_jc_permitted(
        self, tmp_path: Path, blank: str
    ) -> None:
        """Codex Sol Round 3: the schema permits blank second-side
        fields on non-JC rows; ''/'.'/'NA' side2_position cells pass
        while the sibling JC row carries a valid coordinate."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        _write(
            root,
            "benchmark/real_world/results/bl21_user/m5_p3/validation_1"
            "/BL21_G1/unintended.tsv",
            _unintended_with_side2(blank, "500000"),
        )
        assert MODULE.check_evidence_chain_traceability(root) == []

    @pytest.mark.parametrize("blank", ["", ".", "NA"])
    def test_side2_position_blank_on_jc_flagged(
        self, tmp_path: Path, blank: str
    ) -> None:
        """Codex Sol Round 3: a JC row must carry a real second-side
        coordinate; a blank side2_position on a JC row is a finding."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        _write(
            root,
            "benchmark/real_world/results/bl21_user/m5_p3/validation_1"
            "/BL21_G1/unintended.tsv",
            _unintended_with_side2("1234", blank),
        )
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "JC row" in f and "side2_position" in f for f in findings
        )

    def test_placeholder_pam_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "offtarget_sites.tsv"
        )
        rows = target.read_text().splitlines()
        rows[1] = rows[1].replace("\tCGG\t", "\tNGG\t")
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("query-pattern placeholder" in f for f in findings)

    @pytest.mark.parametrize("score_col", ["cfd_score", "hsu_score"])
    def test_unvalidated_score_flagged(self, tmp_path: Path, score_col: str) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "offtarget_sites.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index(score_col)] = "0.45"
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(score_col in f for f in findings)

    def test_association_event_id_cannot_be_na(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "mutation_offtarget_links.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("event_id")] = "NA"
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("association row" in f for f in findings)

    def test_edit_outcomes_dv1_list_validated(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "edit_outcomes.tsv"
        )
        rows = target.read_text().splitlines()
        rows[1] = rows[1].replace(VALID_DV1, "MAYBE_AN_EVENT")
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("matched_event_ids_dv1" in f for f in findings)

    def test_edit_outcomes_orphan_event_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "edit_outcomes.tsv"
        )
        rows = target.read_text().splitlines()
        rows[1] = rows[1].replace(VALID_DV1, GHOST_DV1)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(GHOST_DV1 in f and "does not resolve" in f for f in findings)

    def test_missing_locus_reference_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "post_edit_variants.tsv"
        )
        rows = target.read_text().splitlines()
        row = rows[1].split("\t")
        header = rows[0].split("\t")
        row[header.index("seq_id")] = ""
        row[header.index("position")] = "0"
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("empty seq_id" in f for f in findings)

    def test_ragged_product_row_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        base = "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        _make_v2_product_set(root)
        target = root / base / "post_edit_variants.tsv"
        rows = target.read_text().splitlines()
        rows[1] = "\t".join(rows[1].split("\t")[:-1])  # drop last cell
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("ragged row" in f for f in findings)

    @pytest.mark.parametrize(
        "bad_value", ["", "TBD", "DV1_placeholder", "UNKNOWN_EVENT"]
    )
    def test_manifest_dirty_event_values_flagged(
        self, tmp_path: Path, bad_value: str
    ) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        _write(
            root,
            "benchmark/real_world/manifests/event_recovery.tsv",
            "sample\tdv1_event_id\n"
            f"Cas9_1\t{bad_value}\n",
        )
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("dv1_event_id" in f for f in findings)

    def test_manifest_clean_sentinels_pass(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        _write(
            root,
            "benchmark/real_world/manifests/event_recovery.tsv",
            "sample\tdv1_event_id\tobserved_event_id\n"
            f"Cas9_1\tNOT_EVALUABLE\tNA\n"
            f"Cas9_2\t.\t{VALID_DV1}\n",
        )
        assert MODULE.check_evidence_chain_traceability(root) == []

    def test_non_product_table_site_ids_not_policed(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        _write(
            root,
            "benchmark/offtarget/cas_offinder/prokadiff.tsv",
            "site_id\tchr\tpos\nsite_1\tchr1\t100\n",
        )
        assert MODULE.check_evidence_chain_traceability(root) == []

    def test_unintended_mismatched_position_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 9: unintended row at pos 200 referencing event declared at pos 100 is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        unintended = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "unintended.tsv"
        )
        rows = unintended.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("position")] = "200"
        rows[1] = "\t".join(row)
        unintended.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "is inconsistently referenced with position '200'" in f
            for f in findings
        )

    def test_unintended_mismatched_seq_id_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 9: unintended row at seq_id chr2 referencing event declared at chr1 is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        unintended = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "unintended.tsv"
        )
        rows = unintended.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("seq_id")] = "chr2"
        rows[1] = "\t".join(row)
        unintended.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "is inconsistently referenced with seq_id 'chr2'" in f
            for f in findings
        )

    def test_post_edit_variants_no_supporting_gd_evidence_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 9: post_edit_variants row declaring an event_id without supporting GD evidence is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "post_edit_variants.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("evidence")] = "RA=NA;MC=NA;JC=NA"
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "no supporting GD evidence" in f and "RA=NA;MC=NA;JC=NA" in f
            for f in findings
        )

    @pytest.mark.parametrize(
        "bad_evidence",
        [
            "RA=garbage",
            "FAKE=yes",
            "999",
            "RA=1;;MC=NA;JC=NA",
            "RA=0;RA=1;MC=NA;JC=NA",
        ],
    )
    def test_post_edit_variants_malformed_evidence_grammar_flagged(
        self, tmp_path: Path, bad_evidence: str
    ) -> None:
        """Codex Sol Round 10/11: malformed evidence tokens (RA=garbage, FAKE=yes, bare 999, consecutive semicolons, duplicate keys) must be flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        target = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "post_edit_variants.tsv"
        )
        rows = target.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("evidence")] = bad_evidence
        rows[1] = "\t".join(row)
        target.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "invalid evidence grammar" in f and bad_evidence in f
            for f in findings
        )

    def test_post_edit_variants_unresolved_gd_mutation_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 10: post_edit_variants row not present in sibling .gd file is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        base = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        )
        # Sibling .gd has a mutation at position 42, but NOT chr1:100
        _write(base, "output.gd", "SNP\t1\t.\tchr1\t42\tA\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "does not resolve to any record in sibling GD file(s)" in f
            for f in findings
        )

    def test_post_edit_variants_type_mismatch_gd_mutation_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 11: post_edit_variants row declaring SNP chr1:100 is flagged if sibling GD has DEL at chr1:100."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        base = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        )
        _write(base, "output.gd", "DEL\t1\t.\tchr1\t100\t1\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "does not resolve to any record in sibling GD file(s)" in f
            for f in findings
        )

    def test_post_edit_variants_evidence_only_ra_gd_record_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 11: post_edit_variants row declaring SNP chr1:100 is flagged if sibling GD only has RA evidence record."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        base = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        )
        _write(base, "output.gd", "RA\t1\t.\tchr1\t100\t0\t10\t10\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "does not resolve to any record in sibling GD file(s)" in f
            for f in findings
        )

    def test_unintended_mismatched_end_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 10: unintended row with end diverging from declared end is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        unintended = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "unintended.tsv"
        )
        rows = unintended.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("end")] = "200"
        rows[1] = "\t".join(row)
        unintended.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "inconsistently referenced with end '200'" in f
            for f in findings
        )

    def test_unintended_mismatched_alt_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 10: unintended row with alt diverging from declared alt is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        unintended = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "unintended.tsv"
        )
        rows = unintended.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("alt")] = "G"
        rows[1] = "\t".join(row)
        unintended.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "inconsistently referenced with alt 'G'" in f
            for f in findings
        )

    def test_unintended_mismatched_ref_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 10: unintended row with ref diverging from declared ref is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        unintended = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
            / "unintended.tsv"
        )
        rows = unintended.read_text().splitlines()
        header = rows[0].split("\t")
        row = rows[1].split("\t")
        row[header.index("ref")] = "T"
        rows[1] = "\t".join(row)
        unintended.write_text("\n".join(rows) + "\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "inconsistently referenced with ref 'T'" in f
            for f in findings
        )

    def test_report_md_orphan_event_finding_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 9: report.md finding referencing an undeclared event ID is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        base = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        )
        _write(base, "report.md", f"Finding at chr1:100 for undeclared {GHOST_DV1}\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "orphan event finding" in f and GHOST_DV1 in f
            for f in findings
        )

    def test_report_md_mismatched_position_finding_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 9: report.md finding with coordinates diverging from declared event is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        base = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        )
        _write(base, "report.md", f"Finding at chr1:999 for event {VALID_DV1}\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "inconsistent with position '100'" in f and VALID_DV1 in f
            for f in findings
        )

    def test_report_md_mismatched_contig_finding_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 10: report.md finding referencing event at wrong contig (chr2:100 vs declared chr1:100) is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        base = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        )
        _write(base, "report.md", f"| {VALID_DV1} | chr2:100 |\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "inconsistent with contig 'chr1'" in f and VALID_DV1 in f
            for f in findings
        )

    def test_report_md_heading_finding_flagged(self, tmp_path: Path) -> None:
        """Codex Sol Round 10: report.md heading (# Finding ...) with coordinates diverging from declared event is flagged."""
        root = tmp_path / "repo"
        _make_v2_product_set(root)
        base = (
            root
            / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        )
        _write(base, "report.md", f"# Finding {VALID_DV1} at chr1:999\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "inconsistent with position '100'" in f and VALID_DV1 in f
            for f in findings
        )

    def test_report_md_orphan_with_empty_variants_table_flagged(
        self, tmp_path: Path
    ) -> None:
        """Codex Sol Round 10: report.md with DV1 event must be flagged even when sibling post_edit_variants.tsv is empty."""
        root = tmp_path / "repo"
        base = root / "benchmark/real_world/results/bl21_user/m5_p3/validation_1/BL21_G1"
        base.mkdir(parents=True, exist_ok=True)
        # Create minimal required product structure with 0 events
        _write(base, "offtarget_sites.tsv", "site_id\tseq_id\tstart\tend\tstrand\ttarget_seq\tpam\tmismatches\tbulge_type\tbulge_size\tsearch_backend\tcfd_score\thsu_score\n")
        _write(base, "post_edit_variants.tsv", "variant_id\tseq_id\tposition\tend\tgd_type\tref\talt\torigin_status\tintended_relation\tsize_class\treview_priority\tguide_relation\tguide_mismatches\tpam\tdist_to_site\tmobile_element\tevidence\tgene\tlocus_tag\tfeature_type\tlegacy_class\tevent_id\n")
        _write(base, "report.md", f"| {VALID_DV1} | chr1:100 |\n")
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any(
            "orphan event finding" in f and VALID_DV1 in f
            for f in findings
        )


# ---------------------------------------------------------------------------
# Gate 2: fail-closed scan guarantee
# ---------------------------------------------------------------------------


class TestEvidenceChainFailClosed:
    def test_zero_tables_fail_closed(self, tmp_path: Path) -> None:
        """Codex Sol counterexample: an evidence-chain scan with zero
        tables is a failure, never a clean pass."""
        root = tmp_path / "repo"
        root.mkdir()
        (root / "benchmark").mkdir()
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("zero evidence tables" in f for f in findings)

    def test_missing_benchmark_tree_fail_closed(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        root.mkdir()
        findings = MODULE.check_evidence_chain_traceability(root)
        assert any("zero evidence tables" in f for f in findings)


# ---------------------------------------------------------------------------
# Gate 3: quarantine and G1 hold integrity
# ---------------------------------------------------------------------------


class TestQuarantineHold:
    def test_real_repo_passes(self) -> None:
        assert MODULE.check_quarantine_and_hold_integrity(REPO_ROOT) == []

    def test_clean_registry_and_hold_pass(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_registry(root)
        _make_hold(root)
        assert MODULE.check_quarantine_and_hold_integrity(root) == []

    def test_missing_registry_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_hold(root)
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any("VALIDATION_REGISTRY.tsv" in f and "missing" in f for f in findings)

    def test_missing_hold_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_registry(root)
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any("WAVE11E_HOLD.md" in f for f in findings)

    @pytest.mark.parametrize("bad_status", ["ACTIVE", "COMPLETE", "VALIDATED"])
    def test_registry_status_tampered_flagged(
        self, tmp_path: Path, bad_status: str
    ) -> None:
        root = tmp_path / "repo"
        _make_registry(root, prjna884016_status=bad_status)
        _make_hold(root)
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any(
            "prjna884016" in f and bad_status in f and "PENDING" in f for f in findings
        )

    def test_registry_flag_true_while_pending_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_registry(
            root, prjna1088182_flags=("true", "false", "false", "false")
        )
        _make_hold(root)
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any("variant_validated=true" in f for f in findings)

    @pytest.mark.parametrize(
        "truthy", ["yes", "1", "t", "complete", "TRUE", "Y", "partial"]
    )
    def test_registry_flag_nonfalse_literal_rejected(
        self, tmp_path: Path, truthy: str
    ) -> None:
        """Codex Sol counterexample: ANY non-false quarantine flag value —
        including 'yes', '1', 't', 'complete' — must be rejected."""
        root = tmp_path / "repo"
        _make_registry(
            root, prjna1088182_flags=(truthy, "false", "false", "false")
        )
        _make_hold(root)
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any(f"variant_validated={truthy}" in f for f in findings)
        assert any("rejected" in f for f in findings)

    @pytest.mark.parametrize("falsy", ["false", "no", "0", "f", "n", "NA", ""])
    def test_registry_flag_accepted_false_literals_pass(
        self, tmp_path: Path, falsy: str
    ) -> None:
        root = tmp_path / "repo"
        _make_registry(
            root, prjna1088182_flags=(falsy, "false", "false", "false")
        )
        _make_hold(root)
        assert MODULE.check_quarantine_and_hold_integrity(root) == []

    @pytest.mark.parametrize(
        "drop", ["prjna1088182", "prjna884016"]
    )
    def test_registry_row_missing_flagged(self, tmp_path: Path, drop: str) -> None:
        root = tmp_path / "repo"
        _make_registry(
            root,
            include_prjna1088182=(drop != "prjna1088182"),
            include_prjna884016=(drop != "prjna884016"),
        )
        _make_hold(root)
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any(f"'{drop}'" in f and "missing" in f for f in findings)

    def test_hold_baseline_tampered_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_registry(root)
        _make_hold(root, HOLD_BODY.replace("OVER 24", "OVER 12"))
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any("baseline" in f for f in findings)

    def test_hold_pass_claim_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_registry(root)
        _make_hold(root, HOLD_BODY + "\nG1 is PASSED and the gate is closed.\n")
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any("G1 pass claim" in f for f in findings)

    def test_hold_gate_passed_true_claim_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_registry(root)
        _make_hold(root, HOLD_BODY + "\nFinal state: G1 gate_passed=true.\n")
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any("G1 pass claim" in f for f in findings)

    def test_hold_marker_removed_flagged(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_registry(root)
        _make_hold(root, HOLD_BODY.replace("G1 remains `PROJECT_ACCEPTANCE_PENDING`.", ""))
        findings = MODULE.check_quarantine_and_hold_integrity(root)
        assert any("non-pass marker" in f for f in findings)

    def test_hold_copy_inside_lock_archive_ignored(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        _make_registry(root)
        _make_hold(root)
        _write(
            root,
            "benchmark/results/parity_closure/state/lock_archive_20261103"
            "/WAVE11E_HOLD.md",
            "# tampered copy\nEverything passed: G1 QUALIFIED, baseline "
            "matched 10, OVER 0, UNDER 0.\n",
        )
        assert MODULE.check_quarantine_and_hold_integrity(root) == []


# ---------------------------------------------------------------------------
# Aggregate runner and CLI
# ---------------------------------------------------------------------------


class TestRunAllAndCli:
    def test_real_repo_run_all_passes(self) -> None:
        passed, results = MODULE.run_all_release_gates(REPO_ROOT)
        assert passed is True
        assert set(results.keys()) == {
            "scientific_language",
            "evidence_chain_traceability",
            "quarantine_and_hold_integrity",
        }
        assert all(v == [] for v in results.values())

    def test_run_all_detects_tampered_registry(self, tmp_path: Path) -> None:
        root = _make_clean_repo(tmp_path)
        _make_registry(root, prjna1088182_status="COMPLETE")  # overwrite
        passed, results = MODULE.run_all_release_gates(root)
        assert passed is False
        assert results["quarantine_and_hold_integrity"]
        assert results["scientific_language"] == []
        assert results["evidence_chain_traceability"] == []

    def test_run_all_detects_forbidden_phrase(self, tmp_path: Path) -> None:
        root = _make_clean_repo(tmp_path)
        _write(root, "docs/sneaky.md", "One confirmed off-target event.\n")
        passed, results = MODULE.run_all_release_gates(root)
        assert passed is False
        assert results["scientific_language"]

    def test_run_all_empty_repo_fails_closed(self, tmp_path: Path) -> None:
        """Codex Sol counterexample: zero documents / zero tables must
        fail closed across every gate."""
        root = tmp_path / "repo"
        root.mkdir()
        (root / "docs").mkdir()
        (root / "benchmark").mkdir()
        passed, results = MODULE.run_all_release_gates(root)
        assert passed is False
        assert any("zero scannable documents" in f for f in results["scientific_language"])
        assert any("zero evidence tables" in f for f in results["evidence_chain_traceability"])

    def _run_cli(self, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, str(SCRIPT_PATH), *args],
            capture_output=True,
            text=True,
            cwd=str(cwd) if cwd else None,
        )

    def test_cli_clean_repo_exits_zero(self, tmp_path: Path) -> None:
        root = _make_clean_repo(tmp_path)
        proc = self._run_cli("--repo-root", str(root))
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert "| Gate | Status | Findings |" in proc.stdout
        assert "**Overall: PASS**" in proc.stdout

    def test_cli_tampered_repo_exits_one(self, tmp_path: Path) -> None:
        root = _make_clean_repo(tmp_path)
        _make_hold(root, HOLD_BODY.replace("UNDER 8", "UNDER 3"))
        proc = self._run_cli("--repo-root", str(root))
        assert proc.returncode == 1
        assert "**Overall: FAIL**" in proc.stdout
        assert "Findings — quarantine_and_hold_integrity" in proc.stdout

    def test_cli_empty_repo_exits_one_fail_closed(self, tmp_path: Path) -> None:
        root = tmp_path / "repo"
        root.mkdir()
        (root / "docs").mkdir()
        (root / "benchmark").mkdir()
        proc = self._run_cli("--repo-root", str(root))
        assert proc.returncode == 1
        assert "zero scannable documents" in proc.stdout
        assert "zero evidence tables" in proc.stdout

    def test_cli_default_repo_root_is_cwd(self, tmp_path: Path) -> None:
        root = _make_clean_repo(tmp_path)
        proc = self._run_cli(cwd=root)
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert "**Overall: PASS**" in proc.stdout

    def test_cli_writes_optional_report(self, tmp_path: Path) -> None:
        root = _make_clean_repo(tmp_path)
        report_path = tmp_path / "out" / "gate_report.md"
        proc = self._run_cli("--repo-root", str(root), "--output", str(report_path))
        assert proc.returncode == 0
        assert report_path.is_file()
        assert "ProkaDiff Release Gate Audit" in report_path.read_text()

    def test_cli_missing_root_exits_two(self, tmp_path: Path) -> None:
        proc = self._run_cli("--repo-root", str(tmp_path / "does_not_exist"))
        assert proc.returncode == 2

    def test_cli_repo_root_passes_on_real_repo(self) -> None:
        proc = self._run_cli("--repo-root", str(REPO_ROOT))
        assert proc.returncode == 0, proc.stdout
        assert "**Overall: PASS**" in proc.stdout
