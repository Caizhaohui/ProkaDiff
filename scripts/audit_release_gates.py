#!/usr/bin/env python3
"""ProkaDiff M9 automated release-gate auditor.

Deterministic, pure-standard-library gate checks over a ProkaDiff repository
checkout.  The tool implements the defensive release gates defined by
``M9_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`` (Gate C scientific language,
Gate B evidence-chain traceability, and the quarantine/G1-hold invariants),
hardened after the independent Codex Sol review (2026-11-04) and its
Round 2 and Round 3 findings (2026-11-04), and its Round 4 findings
(2026-11-04: affirmative "claim that" framing is not negation, and a
claim-framing negation separated from the assertion by a contrastive
barrier governed the *prior* clause, not the assertion).

Gates
-----
1. ``check_scientific_language``
   Scans text documents (``.md`` / ``.txt`` / ``.tsv`` / ``.csv``) under
   ``docs/`` and ``benchmark/`` (ignoring ``docs/development/archive/`` and
   the parity-closure lock archives) for:
   * forbidden phrases ``confirmed off-target`` / ``confirmed offtarget`` /
     ``confirmed_off_target`` (exempt: explanatory negative occurrences in
     ``scientific_model.md``);
   * **direct cleavage-causality assertions** ("was caused by Cas9",
     "Cas9 caused ...", "induced by Cas9", "Cas9-induced", ...).  These
     are policed per occurrence, and only the **prefix that precedes
     the match inside its own sentence** (delimited by ``.`` ``;``
     ``!`` ``?`` newline, or their CJK equivalents) can excuse them
     (Codex Sol Round 3): a disclaimer elsewhere in the document does
     NOT excuse them, a negation in a *different* sentence of the same
     paragraph does not (Codex Sol Round 2), and a negation in the
     *suffix after the claim* ("... without ambiguity", "..., not by
     contamination") does not either (Codex Sol Round 3).  Only a
     verbal/adjectival negation directly before the claim ("was not
     caused", "is not Cas9-induced", "does not establish", "cannot
     establish", "fails to establish", "never establishes",
     "non-Cas9-induced") or a governing claim-framing negation in the
     clause prefix ("no claim that", "cannot prove that", "no evidence
     that", "any claim that") exempts a match.  Codex Sol Round 4
     tightened both prefix forms: an affirmative claim frame ("We
     claim that ...") is NOT a negation and never exempts a match;
     the immediate verbal/adjectival negation is evaluated on the
     *immediate* comma/colon-delimited clause of the prefix only; and
     a framing negation exempts the match only when no contrastive
     barrier ("but", "however", "yet", "nevertheless", "nonetheless",
     "although", "whereas", "though", or a comma) intervenes between
     the framing and the assertion — with a barrier, the framing
     governed the prior clause and does not excuse the assertion;
   * affirmative causality-*establishment* claims ("establishes cleavage
     causality", "cleavage causality is established").  A document-wide
     disclaimer does NOT excuse these (Codex Sol Round 2), and neither
     does a negation in the *suffix after the claim* (Codex Sol Round
     3); only a negation directly preceding the claim inside its own
     sentence ("does not establish cleavage causality") does.  The
     required disclaimer sentence
     "Spatial association is consistent with a possible guide-dependent
     event but does not establish cleavage causality." itself passes
     because its own clause prefix carries the negation;
   * curated documents that discuss candidate guide-dependent associations
     (or "mechanistic associations") without carrying the required
     disclaimer.  Machine-generated product reports under
     ``benchmark/results/`` and ``benchmark/real_world/results/`` are
     exempt from this requirement (their hedging boilerplate is fixed by
     the engine template); direct-causality policing still applies there;
   * Rule 6 violations (framing all NCBI-reference differences as
     unintended edits);
   * Rule 8 violations (claiming CAST / IS110 cause no genome-wide
     mutagenesis).

   Fail closed: a scan that finds zero scannable documents (missing or
   empty ``docs/`` + ``benchmark/`` trees) is itself a failure, and
   unreadable or oversize in-scope documents are reported instead of
   silently skipped.

2. ``check_evidence_chain_traceability``
   Audits TSV reports and manifests under ``benchmark/`` so that every
   variant row carries a valid ``DV1_<32 hex>`` event identity or a valid
   locus reference, missing data is the clean literal ``NA`` (or the
   documented manifest sentinels), site references resolve, no fabricated
   identifiers (``SITE_UNKNOWN``, invented scores or placeholder event
   ids) appear, and — per the Codex Sol review — **sibling integrity**
   holds: any product table that references a ``site_id`` requires a
   sibling ``offtarget_sites.tsv`` in the same directory containing that
   id, and any product table that references a ``DV1_`` event id requires
   a sibling ``post_edit_variants.tsv`` containing that id (orphan
   references and missing sibling tables are TRACEABILITY failures).
   Genomic coordinates in product tables must respect bacterial-contig
   sanity bounds (1 .. ``MAX_BACTERIAL_COORD``); every coordinate cell
   must be the clean literal ``NA`` or a positive integer inside those
   bounds — negative ("-1"), non-numeric ("abc"), blank, and
   above-bound ("10000001") cells are errors, never silent skips, and
   blank ``site_id`` reference cells are errors (Codex Sol Round 2).
   Coordinate policing covers ``position`` / ``start`` / ``end`` /
   ``site_start`` / ``site_end`` / ``expected_start`` / ``expected_end``
   and — per Codex Sol Round 3 — ``mutation_position`` and
   ``side2_position``.  ``side2_position`` additionally honours the
   schema's second-side contract: blank sentinels (""/"."/"NA") are
   permitted on non-JC rows (``gd_type`` present and != "JC"), a
   populated value must be a valid coordinate on any row, and a JC row
   must carry a real second-side coordinate (blank ``side2_position``
   on a JC row is itself a failure).  A
   file whose filename is in ``PRODUCT_FILENAMES`` but that carries no
   header row at all (0 bytes / empty) is itself an EVIDCHAIN failure
   instead of a silent skip.

   Fail closed: zero scannable evidence tables is itself a failure, and
   unreadable / oversize product or metadata tables are reported instead
   of silently skipped.  Oversize *non-metadata* bulk tables (raw engine
   artefacts such as ``repeat_ambiguous_junctions.tsv``) remain out of
   scope by the documented size bound.

3. ``check_quarantine_and_hold_integrity``
   Fails closed unless ``benchmark/real_world/VALIDATION_REGISTRY.tsv``
   keeps ``prjna1088182`` and ``prjna884016`` rows ``PENDING`` with all
   validation flags false, and ``WAVE11E_HOLD.md`` keeps the G1 strict
   gate un-passed with the frozen baseline ``matched 10, OVER 24,
   UNDER 8``.  Quarantine flags are accepted as false only for the
   literals ``false / no / 0 / f / n / <empty> / na`` (case-insensitive);
   ANY other value ("true", "yes", "1", "t", "complete", ...) is
   rejected.

Exit status of the CLI: ``0`` when every gate passes, ``1`` on any finding.

Limitations (documented, deterministic): files above
``SCAN_SIZE_LIMIT_BYTES`` are skipped by the tree scans **only** when they
are non-metadata bulk tables; in-scope markdown/text documents and
metadata/product tables above the limit are flagged (fail closed) rather
than silently skipped.  Raw engine artefacts (``.gd``, ``.html``, ``.log``,
``.err``, indexes, FASTQ/BAM) are out of scope by extension.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Iterable, Iterator

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Text extensions scanned by the scientific-language gate.
TEXT_EXTENSIONS = {".md", ".txt", ".tsv", ".csv"}

#: Files above this size are non-metadata bulk tables and are skipped by
#: the tree scans (the walk stays deterministic and bounded).  In-scope
#: documents (.md/.txt) and metadata/product tables above the limit are
#: *flagged* instead of skipped (fail closed, no silent whitewashing).
SCAN_SIZE_LIMIT_BYTES = 2 * 1024 * 1024

#: Directories never descended into.
SKIPPED_DIR_NAMES = {"__pycache__", ".git"}

#: Subtrees never scanned (relative posix path fragments).
SKIPPED_PATH_FRAGMENTS = (
    "docs/development/archive",
    "benchmark/results/parity_closure/state/lock_archive",
)

#: The canonical scientific disclaimer required alongside causality language.
REQUIRED_DISCLAIMER = (
    "Spatial association is consistent with a possible guide-dependent event "
    "but does not establish cleavage causality."
)

DISCLAIMER_CLAUSE_RE = re.compile(
    r"does\s+not\s+establish\s+cleavage\s+causality", re.IGNORECASE
)

#: Forbidden off-target finality phrases (whitespace/punctuation tolerant).
CONFIRMED_OFFTARGET_RE = re.compile(r"confirmed[\s_-]+off[\s_-]?target", re.IGNORECASE)

#: Negation / prohibition markers that make an occurrence "explanatory
#: negative" usage (allowed only inside ``scientific_model.md``).
NEGATION_MARKERS_RE = re.compile(
    r"严禁|禁止|不得|不应|不能|避免|禁止使用"
    r"|\bnever\b|\bnot\b|\bforbidden\b|\bprohibited\b|\bavoid\b|\banned\b",
    re.IGNORECASE,
)

#: Codex Sol Round 3: a causality match is excused only by its own
#: clause **prefix** — the text between the start of the match's
#: sentence (delimited by ``SENTENCE_DELIMITER_RE``) and the match
#: itself.  A negation in the *suffix* after the match ("without
#: ambiguity", ", not by contamination"), in another sentence of the
#: paragraph, or elsewhere in the document never exempts it.  Two
#: prefix forms count.
#:
#: 1. Immediate verbal/adjectival negation anchored directly before the
#:    match (optional auxiliary + negator + bounded run of
#:    particles/adverbs): "was not caused", "is not Cas9-induced",
#:    "does not establish", "cannot establish", "fails to establish",
#:    "without establishing", "never establishes", "non-Cas9-induced",
#:    or the Chinese 并非由 Cas9 引起 / 不建立切割因果 forms.
CAUSALITY_IMMEDIATE_NEGATION_RE = re.compile(
    r"\b(?:(?:am|is|are|was|were|be|been|being|do|does|did|has|have|had"
    r"|will|would|can|could|may|might|must|shall|should)\s+)?"
    r"(?:not|never|no|cannot|neither|nor|without|unproven|unsupported"
    r"|undetermined|unclear|unconfirmed|lacks?|lacking|absent|fails?"
    r"|failed|failing|refutes?|refuted|refuting|contradicts?"
    r"|contradicted|contradicting|disputes?|disputed|disputing"
    r"|excludes?|excluded|excluding|precludes?|precluded|precluding"
    r"|denies?|denied|denying|rejects?|rejected|rejecting|avoids?"
    r"|avoided|avoiding|prohibits?|prohibited|prohibiting|forbids?"
    r"|forbidden|forbidding|bans?|banned|banning)"
    r"(?:\s+(?:to|be|been|being|directly|necessarily|specifically"
    r"|strictly|formally|definitively|convincingly|rigorously|clearly"
    r"|explicitly|currently|yet|thereby))*"
    r"\s*$"
    r"|\bnon[- ]?\s*$"
    r"|(?:并非|并不|不是|不会|不能|不得|没有|无法|未曾|未能|未|非|不)"
    r"(?:由|是|因|被|建立|证明|证实|能|会|再|可)?\s*$",
    re.IGNORECASE,
)

#: 2. Governing claim-framing negation anywhere in the clause prefix:
#:    "no claim that", "no claim is made that", "cannot prove that",
#:    "no evidence that", "any claim that" (scope-limiting enumeration,
#:    as in the frozen M5_CLOSEOUT wording), "do not / never claim
#:    that", prohibition framing, and the corresponding Chinese framing
#:    forms.  Codex Sol Round 4: a **bare affirmative claim frame**
#:    ("We claim that ...", "claim of", "claim to") is NOT negation and
#:    must never match here — only genuinely negative framing forms
#:    are accepted.
CAUSALITY_FRAMING_NEGATION_RE = re.compile(
    r"\bno\s+(?:claim|evidence|proof|demonstration|assertion|statement"
    r"|suggestion|indication|support|basis|data|reason)\b"
    r"|\bnot\s+(?:a\s+|the\s+|any\s+)?(?:claim|proof|evidence"
    r"|demonstration|assertion)\b"
    r"|\b(?:do\s+not|does\s+not|did\s+not|don'?t|doesn'?t|didn'?t"
    r"|never)\s+claim(?:s|ed|ing)?\b"
    r"|\blittle\s+(?:evidence|support|reason)\b"
    r"|\b(?:cannot|can'?t|could\s+not|couldn'?t|will\s+not|won'?t"
    r"|shall\s+not|should\s+not|must\s+not|may\s+not|might\s+not"
    r"|fail(?:s|ed|ing)?\s+to|unable\s+to)\s+"
    r"(?:prove|establish|show|demonstrate|confirm|verify|determine"
    r"|conclude|rule\s+out|support|certify|guarantee)\b"
    r"|\b(?:unproven|unsupported|unsubstantiated|false|disputed)\s+(?:claim|assertion|hypothesis)\s+that\b"
    r"|\bprohibit\w*|\bforbid\w*|\binterdict\w*"
    r"|\b(?:reject\w*|dismiss\w*|den(?:y|ies|ied))\s+(?:any\s+claim\s+that|the\s+claim\s+that)\b"
    r"|\bwithout\s+(?:establishing|proving|claiming|asserting)\b"
    r"|无法证明|没有证据|并无证据|并非声称|不声称|未证实|未建立|非因果"
    r"|无法建立|禁止|严禁|不得|不应主张",
    re.IGNORECASE,
)

#: Codex Sol Round 4 & 5: an adversative contrast between a claim-framing
#: negation and the causality assertion means the framing governed the
#: *prior* clause and does NOT govern the assertion ("There is no
#: evidence of contamination, but the deletion was caused by Cas9.",
#: "...: this evidence establishes...", "..., the deletion was caused...").
#: Any contrastive conjunction/adverb, or punctuation (comma, semicolon,
#: colon, or CJK equivalents) counts as the barrier.
CONTRASTIVE_BARRIER_RE = re.compile(
    r"\b(?:but|however|yet|nevertheless|nonetheless|although|whereas"
    r"|though)\b"
    r"|[,;:，：；、]",
    re.IGNORECASE,
)

#: Sentence delimiters bounding the "exact sentence" used for guard
#: evaluation (ASCII ``. ; ! ?`` and newline, plus CJK ``。！？；``).
SENTENCE_DELIMITER_RE = re.compile(r"[.;!?\n。！？；]")

#: Codex Sol Round 4 & 5: clause delimiters bounding the *immediate* clause
#: on which a verbal/adjectival negation is evaluated (ASCII comma,
#: semicolon, colon, plus CJK ``，`` ``：`` ``；`` ``、``).  A negation in
#: an earlier clause of the same sentence never whitewashes a later affirmative
#: clause.
CLAUSE_DELIMITER_RE = re.compile(r"[,:;，：；、]")

#: Codex Sol Round 6 & 7: explicit scope-prohibition lead-in pattern anchored with colon:
#: "Out of scope for this closure:", "Out of scope:", "Scope limitation:", "超出范围："
SCOPE_PROHIBITION_LEADIN_RE = re.compile(
    r"(?:\b(?:out\s+of\s+scope|outside\s+(?:the\s+)?scope|excluded\s+from\s+(?:the\s+)?scope"
    r"|not\s+(?:with)?in\s+(?:the\s+)?scope|scope\s+limitation|scope\s+exclusion)\b"
    r"|超出范围|不在范围|范围之外|不在此次.*范围)"
    r"(?:\s+for\s+[^:：\n]+)?"
    r"\s*[:：]",
    re.IGNORECASE,
)

#: Codex Sol Round 6 & 7: conservative accepted grammar for scope-prohibition enumerations.
#: In an out-of-scope list, 'any claim that' (or Chinese 任何声称/主张) must be an enumeration
#: item directly introduced by the lead-in colon, or by a list delimiter (comma / enumeration comma)
#: optionally followed by a coordination conjunction (and / or / as well as / 以及 / 及 / 或).
#: Arbitrary text (e.g. 'the team accepts', 'we assert', 'our position confirms') between delimiter
#: and claim is structurally rejected without relying on an open-ended keyword blacklist.
ITEM_CLAIM_LEADIN_RE = re.compile(
    r"(?:[:：]|[,，、]\s*(?:and|or|as\s+well\s+as|以及|及|或)?)\s*"
    r"(?:any\s+claim\s+that\b|任何(?:声称|主张))",
    re.IGNORECASE,
)

#: Cheap literal prefilter: a line containing none of these tokens cannot
#: match any DIRECT_CAUSALITY_RES pattern (perf guard for large TSV scans).
DIRECT_CAUSALITY_PREFILTER_RE = re.compile(
    r"caused|mediated|introduced|induced|cas\d|grna|guide", re.IGNORECASE
)

#: Direct causality assertion patterns (English).  Policed per occurrence
#: at the exact-sentence level: a disclaimer elsewhere in the document —
#: or a negation in another sentence of the paragraph — does not excuse
#: them; only a guard inside the match's own sentence does.  (Codex Sol
#: P1: a document-wide disclaimer must not whitewash a direct "was caused
#: by Cas9" claim; Codex Sol Round 2: neither may "No contamination was
#: detected." excuse "The deletion was caused by Cas9." in the next
#: sentence, and the adjective form "Cas9-induced" must be caught.)
DIRECT_CAUSALITY_RES: tuple[re.Pattern[str], ...] = (
    re.compile(
        r"\b(?:was|were|is|are|been|be)?\s*(?:caused|mediated|introduced"
        r"|induced)\s+by\s+"
        r"(?:Cas9|Cas12a|CAS\d+|the\s+guide(?:\s+RNA)?|gRNA)\b",
        re.IGNORECASE,
    ),
    re.compile(
        r"\bCas(?:9|12a)?\s+(?:nuclease\s+)?(?:cleavage\s+)?(?:directly\s+)?"
        r"(?:caused|mediated)\s+(?:this|the|a|an|these)\s+"
        r"(?:mutation|variant|deletion|insertion|junction|event)s?\b",
        re.IGNORECASE,
    ),
    re.compile(r"\bCas(?:9|12a)\s+(?:directly\s+)?caused\b", re.IGNORECASE),
    re.compile(r"\binduced\s+by\s+Cas(?:9|12a)?\b", re.IGNORECASE),
    # Codex Sol Round 2: the adjective form ("The deletion is
    # Cas9-induced.").
    re.compile(r"\bCas(?:9|12a)?[- ]induced\b", re.IGNORECASE),
)

#: Causality-*establishment* claims (English).  An affirmative
#: establishment claim is itself a violation: a document-wide disclaimer
#: does NOT excuse it; only a negation inside the claim's own sentence
#: ("does not establish cleavage causality") does (Codex Sol Round 2).
CAUSALITY_ESTABLISHMENT_RES: tuple[re.Pattern[str], ...] = (
    re.compile(
        r"\b(?:establish(?:es)?|established|establishing|confirms?|confirmed|confirming|proves?|proven|proving"
        r"|demonstrates?|demonstrated|demonstrating)\s+cleavage\s+causality\b",
        re.IGNORECASE,
    ),
    re.compile(
        r"\bcleavage\s+causality\s+(?:is|was)\s+"
        r"(?:established|confirmed|proven|demonstrated)\b",
        re.IGNORECASE,
    ),
)

#: Discussion of candidate guide-dependent / mechanistic associations.
#: Curated documents matching this must carry the required disclaimer.
GUIDE_ASSOCIATION_DISCUSSION_RE = re.compile(
    r"candidate\s+guide-?dependent"
    r"|guide-?dependent\s+(?:off[\s_-]?target|event|association)s?"
    r"|mechanistic\s+associations?",
    re.IGNORECASE,
)

#: Machine-generated product-report subtrees: exempt from the Tier-C
#: "discussion requires disclaimer" rule (engine-template boilerplate
#: already hedges them), but NOT from direct-causality policing.
GENERATED_PRODUCT_REPORT_ROOTS = (
    "benchmark/results/",
    "benchmark/real_world/results/",
)

#: Rule 6: reporting all differences relative to NCBI as unintended edits.
RULE6_RES: tuple[re.Pattern[str], ...] = (
    re.compile(
        r"\ball\s+(?:genome\s+)?differences?\s+"
        r"(?:vs\.?|versus|against|relative\s+to|compared\s+(?:to|with))\s+"
        r"(?:the\s+)?NCBI\b",
        re.IGNORECASE,
    ),
    re.compile(
        r"\b(?:unintended|unexpected)\s+(?:edits?|mutations?|variants?"
        r"|differences?|calls?)\s+"
        r"(?:vs\.?|versus|against|relative\s+to|compared\s+(?:to|with))\s+"
        r"(?:the\s+)?NCBI\b",
        re.IGNORECASE,
    ),
)

#: Rule 8: editor tokens whose mutagenesis claims are policed.
RULE8_EDITOR_RE = re.compile(
    r"\bCAST\b|\bIS110\b|`cast`|`is110`|--editor[= ](?:cast|is110)\b"
)

#: Rule 8: assertions that the editors cause no genome-wide mutagenesis.
RULE8_NO_MUTAGENESIS_RE = re.compile(
    r"\bno\s+genome-?wide\s+mutagen\w*"
    r"|\bmutagen\w*-free\b"
    r"|(?:does\s+not|do\s+not|doesn'?t|never|cannot)\s+"
    r"(?:cause|produce|generate|perform|undergo|introduce)[^\n.]{0,60}"
    r"\bmutagen(?:esis|ts)\b",
    re.IGNORECASE,
)

#: Rule 8 guard: meta-statements that *prohibit* the claim, not make it.
RULE8_PROHIBITION_GUARD_RE = re.compile(
    r"(?:never|not|严禁|禁止|不得)\s*(?:claim|say|state|assert|report|hard-?code"
    r"|write|label)",
    re.IGNORECASE,
)

# --- Evidence-chain identifiers ---------------------------------------------

DV1_RE = re.compile(r"^DV1_[0-9a-f]{32}$")
DV1_SEARCH_RE = re.compile(r"\bDV1_[0-9a-f]{32}\b")
SITE_ID_RE = re.compile(r"^SITE_\d+$")
VARIANT_ID_RE = re.compile(r"^VAR_\d+$")
LEGACY_MUTATION_ID_RE = re.compile(r"^mut_\d+$")

#: Clean "missing" literals for product TSV event-identity columns.
PRODUCT_MISSING = {"NA"}

#: Clean "missing" literals for manifest event-identity columns.
MANIFEST_MISSING = {"NA", "NOT_EVALUABLE", "."}

#: Header names policed by the generic manifest event-column rule.
EVENT_COLUMN_EXACT = {"event_id", "event_id_if_applicable", "dv1_event_id"}
EVENT_COLUMN_SUFFIXES = ("_event_id", "_event_ids", "_dv1", "event_ids_dv1")

#: Query-pattern PAM placeholders that must never fill an observed PAM cell.
PAM_PLACEHOLDERS = {"NGG", "NAG", "NRG", "NNN", "TTV", "TTTV", "N"}

#: Product TSV filenames with dedicated audit rules.  Includes
#: ``mechanistic_associations.tsv`` (Codex Sol example): it receives the
#: generic product-table treatment (site/event reference integrity).
PRODUCT_FILENAMES = {
    "post_edit_variants.tsv",
    "unintended.tsv",
    "mutation_offtarget_links.tsv",
    "mechanistic_associations.tsv",
    "offtarget_sites.tsv",
    "edit_outcomes.tsv",
}

#: Metadata tables that must never be skipped by the size bound.
METADATA_TABLE_NAMES = PRODUCT_FILENAMES | {"validation_registry.tsv"}

#: Sibling tables required to resolve references inside product tables.
SITES_SIBLING_NAME = "offtarget_sites.tsv"
VARIANTS_SIBLING_NAME = "post_edit_variants.tsv"

#: Sanity bound for genomic coordinates on standard bacterial contigs
#: (largest known bacterial chromosomes are ~13 Mb; 10 Mb comfortably
#: covers E. coli / BL21 / REL606 while flagging absurd values).
MAX_BACTERIAL_COORD = 10_000_000

#: Coordinate-bearing column names policed in product tables: every cell
#: must be clean ``NA`` or a positive integer within
#: ``1 .. MAX_BACTERIAL_COORD`` (anything else is an error).  Codex Sol
#: Round 3 added ``mutation_position`` and ``side2_position``.
COORDINATE_COLUMN_NAMES = (
    "position",
    "start",
    "end",
    "site_start",
    "site_end",
    "expected_start",
    "expected_end",
    "mutation_position",
    "side2_position",
)

#: Blank sentinels permitted for ``side2_position`` on non-JC rows
#: (``gd_type`` present and != "JC"): the schema leaves the second-side
#: fields empty for single-side calls.  A JC row must instead carry a
#: real second-side coordinate (a blank is itself a failure).
SIDE2_BLANK_LITERALS = ("", ".", "NA")

# --- Quarantine / G1 hold ----------------------------------------------------

REGISTRY_RELPATH = Path("benchmark/real_world/VALIDATION_REGISTRY.tsv")
HOLD_FILENAME = "WAVE11E_HOLD.md"

QUARANTINED_DATASETS = ("prjna1088182", "prjna884016")
REGISTRY_FLAG_COLUMNS = (
    "variant_validated",
    "intended_validated",
    "mob_validated",
    "jc_validated",
)

#: The only literals accepted as "false" for quarantine validation flags.
#: ANY other value ("true", "yes", "1", "t", "complete", ...) is rejected.
FALSE_FLAG_LITERALS = {"false", "no", "0", "f", "n", "", "na"}

#: Frozen G1 strict-gate baseline recorded in ``WAVE11E_HOLD.md``.
G1_BASELINE_RE = re.compile(
    r"matched\s*10,\s*OVER\s*24,\s*UNDER\s*8", re.IGNORECASE
)
G1_NONPASS_MARKER_RE = re.compile(
    r"PROJECT_ACCEPTANCE_PENDING|pass\s*=\s*false|\bpass\s+false\b",
    re.IGNORECASE,
)
G1_PASS_CLAIM_RE = re.compile(
    r"\bG1\b[^.\n]{0,200}?"
    r"(?:PASS\s*=\s*true|gate_passed\s*[:=]\s*true|\bPASSED\b|\bQUALIFIED\b)",
    re.IGNORECASE,
)

GATE_SCIENTIFIC_LANGUAGE = "scientific_language"
GATE_EVIDENCE_CHAIN = "evidence_chain_traceability"
GATE_QUARANTINE_HOLD = "quarantine_and_hold_integrity"

GATE_ORDER = (
    GATE_SCIENTIFIC_LANGUAGE,
    GATE_EVIDENCE_CHAIN,
    GATE_QUARANTINE_HOLD,
)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _under_skipped_fragment(rel_posix: str) -> bool:
    """True when a repo-relative posix path equals or is under a skipped subtree."""
    return any(
        rel_posix == frag or rel_posix.startswith(frag + "/")
        for frag in SKIPPED_PATH_FRAGMENTS
    )


def _is_skipped_dir(rel_parent_posix: str, dir_name: str) -> bool:
    """True when a child directory must not be descended into."""
    if dir_name in SKIPPED_DIR_NAMES:
        return True
    if dir_name.startswith("lock_archive"):
        return True
    return _under_skipped_fragment(f"{rel_parent_posix}/{dir_name}")


def _is_metadata_table(rel_posix: str) -> bool:
    """True for curated tables that the size bound must never skip."""
    parts = Path(rel_posix).parts
    if Path(rel_posix).name.lower() in METADATA_TABLE_NAMES:
        return True
    return any(part == "manifests" for part in parts[:-1])


def _is_generated_product_report(rel_posix: str) -> bool:
    """True inside machine-generated product-report subtrees."""
    return any(
        rel_posix == root.rstrip("/") or rel_posix.startswith(root)
        for root in GENERATED_PRODUCT_REPORT_ROOTS
    )


def _walk_files(
    repo_root: Path,
    base: str,
    extensions: set[str],
    findings: list[str] | None = None,
) -> Iterator[tuple[Path, str]]:
    """Yield ``(absolute_path, repo-relative posix path)`` for scan targets.

    The walk is deterministic (sorted directories and files), skips the
    archive and lock subtrees, and ignores dot files.  Files above
    ``SCAN_SIZE_LIMIT_BYTES`` are skipped **only** when they are
    non-metadata bulk tables; in-scope markdown/text documents and
    metadata/product tables above the limit (or with unreadable stat) are
    reported as fail-closed findings when a ``findings`` list is given.
    """
    base_dir = repo_root / base
    if not base_dir.is_dir():
        return
    for root, dirs, files in os.walk(base_dir):
        dirs.sort()
        files.sort()
        rel_root_posix = Path(root).relative_to(repo_root).as_posix()
        if _under_skipped_fragment(rel_root_posix) and rel_root_posix != base:
            continue
        dirs[:] = [d for d in dirs if not _is_skipped_dir(rel_root_posix, d)]
        for name in files:
            if Path(name).suffix.lower() not in extensions:
                continue
            if name.startswith("."):
                continue
            abs_path = Path(root) / name
            rel_posix = f"{rel_root_posix}/{name}"
            try:
                size = abs_path.stat().st_size
            except OSError:
                if findings is not None:
                    findings.append(
                        _finding(
                            "FAILCLOSED",
                            rel_posix,
                            "in-scope file cannot be stat-ed; the audit "
                            "refuses to silently skip it",
                        )
                    )
                continue
            if size > SCAN_SIZE_LIMIT_BYTES:
                document_like = Path(name).suffix.lower() in {".md", ".txt"}
                if document_like or _is_metadata_table(rel_posix):
                    if findings is not None:
                        findings.append(
                            _finding(
                                "FAILCLOSED",
                                rel_posix,
                                f"file size {size} exceeds "
                                f"SCAN_SIZE_LIMIT_BYTES "
                                f"({SCAN_SIZE_LIMIT_BYTES}); in-scope "
                                f"documents and metadata tables must not "
                                f"be silently skipped",
                            )
                        )
                    continue
                continue  # non-metadata bulk table: out-of-scope artefact
            yield abs_path, rel_posix


def _read_text(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None


def _read_tsv_rows(path: Path) -> list[list[str]] | None:
    """Read a TSV as rows (``None`` = unreadable, ``[]`` = empty)."""
    text = _read_text(path)
    if text is None:
        return None
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        return []
    return [row for row in csv.reader(lines, delimiter="\t")]


def _finding(gate_tag: str, rel_path: str, detail: str) -> str:
    return f"[{gate_tag}] {rel_path}: {detail}"


def _sentence_prefix_at(lines: list[str], lineno_idx: int, col: int) -> str:
    """Return the reconstructed text of the current sentence preceding ``col`` on line ``lineno_idx``.

    Walks backwards within the current paragraph/block to reconstruct multi-line sentences
    (e.g. Markdown line-wrapped lists), stopping at any sentence delimiter (. ; ! ?
    。 ； ！ ？) or paragraph/header boundary.
    """
    current_line = lines[lineno_idx]
    delims = list(SENTENCE_DELIMITER_RE.finditer(current_line, 0, col))
    if delims:
        return current_line[delims[-1].end() : col]

    accum = [current_line[:col]]
    for j in range(lineno_idx - 1, -1, -1):
        prev_line = lines[j]
        prev_stripped = prev_line.strip()
        if not prev_stripped or prev_stripped.startswith(("#", "---", "===")):
            break
        prev_delims = list(SENTENCE_DELIMITER_RE.finditer(prev_line))
        if prev_delims:
            accum.insert(0, prev_line[prev_delims[-1].end() :])
            break
        accum.insert(0, prev_line)

    return " ".join(part.strip() for part in accum if part.strip())


def _is_scope_prohibition_governing(prefix: str) -> bool:
    """True when an explicit scope prohibition lead-in governs an 'any claim that'
    enumeration item in prefix under a conservative accepted grammar (Codex Sol Round 6 & 7).

    Governing requires:
    1. A scope prohibition lead-in with colon ("Out of scope for this closure:", "Scope limitation:", "超出范围：").
    2. The claim phrase ("any claim that", "任何声称/主张") must be a direct enumeration item
       introduced by the colon or by a list delimiter (comma optionally followed by and/or/as well as).
       Arbitrary subject/verb clauses between the delimiter and the claim phrase are structurally rejected.
    3. Between lead-in colon and item delimiter: no sentence delimiters or contrastive barriers.
    4. Remainder between claim phrase and causal assertion: only the that-clause noun phrase subject
       (no clause delimiters, list commas, semicolons, colons, or contrastive barriers).
    """
    lead = SCOPE_PROHIBITION_LEADIN_RE.search(prefix)
    if not lead:
        return False

    claim_items = list(ITEM_CLAIM_LEADIN_RE.finditer(prefix, lead.start()))
    if not claim_items:
        return False
    claim_item = claim_items[-1]

    lead_colon_pos = lead.end()
    item_delim_pos = claim_item.start()
    if item_delim_pos < lead_colon_pos - 1:
        return False

    between = prefix[lead_colon_pos:item_delim_pos]
    if re.search(r"[.;!?。！？；]", between):
        return False
    if re.search(
        r"\b(?:but|however|yet|nevertheless|nonetheless|although|whereas|though)\b",
        between,
        re.IGNORECASE,
    ):
        return False

    remainder = prefix[claim_item.end() :]
    if re.search(
        r"[,，、.;:!?。！？；：]|\b(?:but|however|yet|nevertheless|nonetheless|although|whereas|though)\b",
        remainder,
        re.IGNORECASE,
    ):
        return False

    return True


def _is_causality_assertion_negated(
    line: str,
    start: int,
    end: int,
    sentence_prefix: str = "",
) -> bool:
    """True when the causality assertion at ``line[start:end]`` is negated
    or claim-framed by its own sentence prefix.

    Codex Sol Round 3: only the text that **precedes** the match inside
    its sentence (delimited by ``.`` ``;`` ``!`` ``?`` a newline, or
    their CJK equivalents) can excuse it.  Two prefix forms count:

    * an immediate verbal/adjectival negation anchored directly before
      the match ("was not caused", "is not Cas9-induced", "does not
      establish", "cannot establish", "fails to establish", "without
      establishing", "never establishes", "non-Cas9-induced", and the
      Chinese 并非由 Cas9 引起 / 不建立切割因果 forms); or
    * a governing claim-framing negation anywhere in the prefix ("no
      claim that", "cannot prove that", "no evidence that").

    A negation in the *suffix* after the match ("... without ambiguity",
    "..., not by contamination") never excuses the affirmative claim,
    and neither does a guard in a different sentence of the paragraph
    (Codex Sol Round 2).  ``end`` is accepted for call-site symmetry;
    the suffix is deliberately never consulted.

    Codex Sol Round 4 tightened both prefix forms:

    * the immediate verbal/adjectival negation is evaluated on the
      **immediate clause** of the prefix only (``CLAUSE_DELIMITER_RE``
      splits at ASCII/CJK commas and colons), so a negation sitting in
      an earlier clause of the same sentence ("It is not established:
      the deletion was caused by Cas9.") cannot whitewash a later
      affirmative clause; and
    * a claim-framing negation only excuses the match when it governs
      the assertion, i.e. when **no contrastive barrier**
      (``CONTRASTIVE_BARRIER_RE``: "but", "however", "yet",
      "nevertheless", "nonetheless", "although", "whereas", "though",
      or a comma/semicolon) intervenes between the framing match and
      the causal assertion.  With a barrier the framing governed the
      *prior* clause ("There is no evidence of contamination, but the
      deletion was caused by Cas9.") and does not excuse the
      assertion.  An affirmative claim frame ("We claim that ...") is
      not in ``CAUSALITY_FRAMING_NEGATION_RE`` at all and never
      excuses anything.

    Codex Sol Round 5 & 6: scope-limiting enumerations ("Out of scope: ... and
    any claim that ...") are governed by ``_is_scope_prohibition_governing``,
    which requires an explicit colon-anchored lead-in and strictly rejects
    unrelated scope prohibitions, cross-sentence leakage, clause boundaries,
    and affirmative endorsement/assertion clauses.
    """
    del end  # the suffix after the match must never excuse the claim
    if not sentence_prefix:
        prefix_start = 0
        for delimiter in SENTENCE_DELIMITER_RE.finditer(line, 0, start):
            prefix_start = delimiter.end()
        prefix = line[prefix_start:start]
    else:
        prefix = sentence_prefix

    # (a) Immediate verbal/adjectival negation, evaluated on the
    #     immediate comma/colon-delimited clause only (Codex Sol
    #     Round 4).
    immediate_clause = CLAUSE_DELIMITER_RE.split(prefix)[-1].strip()
    if CAUSALITY_IMMEDIATE_NEGATION_RE.search(immediate_clause):
        return True

    # (b) Governing claim-framing negation: valid only when it reaches
    #     the assertion with no contrastive barrier in between (Codex
    #     Sol Round 4).
    for framing in CAUSALITY_FRAMING_NEGATION_RE.finditer(prefix):
        if not CONTRASTIVE_BARRIER_RE.search(prefix, framing.end()):
            return True

    # (c) Governing scope-prohibition enumeration (Codex Sol Round 5 & 6).
    if _is_scope_prohibition_governing(prefix):
        return True

    return False


# ---------------------------------------------------------------------------
# Gate 1: scientific language
# ---------------------------------------------------------------------------


def check_scientific_language(repo_root: Path) -> list[str]:
    """Scan ``docs/`` and ``benchmark/`` for forbidden scientific wording.

    Returns a sorted list of human-readable findings; an empty list means the
    gate passes.  Fails closed when zero documents can be scanned.
    """
    repo_root = Path(repo_root)
    findings: list[str] = []
    scanned = 0
    for base in ("docs", "benchmark"):
        for abs_path, rel_path in _walk_files(
            repo_root, base, TEXT_EXTENSIONS, findings
        ):
            text = _read_text(abs_path)
            if text is None:
                findings.append(
                    _finding(
                        "FAILCLOSED",
                        rel_path,
                        "in-scope document is unreadable; the audit refuses "
                        "to silently skip it",
                    )
                )
                continue
            scanned += 1
            findings.extend(_scan_document_language(rel_path, text))
    if scanned == 0:
        findings.append(
            _finding(
                "SCILANG",
                "docs/, benchmark/",
                "zero scannable documents found (fail closed): an empty or "
                "missing scan scope cannot be certified as clean",
            )
        )
    return sorted(set(findings))


def _scan_document_language(rel_path: str, text: str) -> list[str]:
    findings: list[str] = []
    lines = text.splitlines()
    is_science_model = Path(rel_path).name == "scientific_model.md"

    # (a) Forbidden off-target finality phrases.
    for lineno, line in enumerate(lines, start=1):
        if "confirmed" not in line.lower():
            continue
        for match in CONFIRMED_OFFTARGET_RE.finditer(line):
            if is_science_model and NEGATION_MARKERS_RE.search(line):
                continue
            findings.append(
                _finding(
                    "SCILANG",
                    rel_path,
                    f"line {lineno}: forbidden phrase "
                    f"'{match.group(0)}' (use 'candidate guide-dependent "
                    f"off-target'); only negative explanatory occurrences in "
                    f"scientific_model.md are allowed",
                )
            )

    # (b1) Direct causality assertions: clause-prefix guard only.  A
    #      document-wide disclaimer, a negation in another sentence of
    #      the same paragraph (Codex Sol Round 2), or a negation in the
    #      suffix after the claim (Codex Sol Round 3) never excuses
    #      these; only a negation/claim-framing prefix directly
    #      preceding the match inside its own sentence does.
    for lineno, line in enumerate(lines, start=1):
        low = line.lower()
        if not (
            "caused" in low
            or "mediated" in low
            or "introduced" in low
            or "induced" in low
            or "guide" in low
            or "grna" in low
            or "cas" in low
        ):
            continue
        if not DIRECT_CAUSALITY_PREFILTER_RE.search(line):
            continue
        accepted_spans: list[tuple[int, int]] = []
        for pattern in DIRECT_CAUSALITY_RES:
            for match in pattern.finditer(line):
                if any(
                    match.start() < end and start < match.end()
                    for start, end in accepted_spans
                ):
                    continue  # already covered by another pattern
                accepted_spans.append((match.start(), match.end()))
                sentence_prefix = _sentence_prefix_at(
                    lines, lineno - 1, match.start()
                )
                if _is_causality_assertion_negated(
                    line, match.start(), match.end(), sentence_prefix
                ):
                    continue
                findings.append(
                    _finding(
                        "SCILANG",
                        rel_path,
                        f"line {lineno}: cleavage-causality assertion "
                        f"'{match.group(0)}' is a direct causality claim; a "
                        f"disclaimer elsewhere in the document, a negation "
                        f"in another sentence, or a negation in the suffix "
                        f"after the claim does not excuse it (only a "
                        f"negation directly preceding the claim in its own "
                        f"clause — such as 'was not caused', 'does not "
                        f"establish' — or framing such as 'no claim that' / "
                        f"'cannot prove' does); required disclaimer for "
                        f"association language: '{REQUIRED_DISCLAIMER}'",
                    )
                )

    # (b2) Causality-establishment claims: an affirmative establishment
    #      claim is itself a violation.  A document-wide disclaimer does
    #      NOT excuse it (Codex Sol Round 2), and neither does a negation
    #      in the suffix after the claim (Codex Sol Round 3); only a
    #      negation directly preceding the claim in its own sentence
    #      ("does not establish cleavage causality") does.
    for lineno, line in enumerate(lines, start=1):
        if "causality" not in line.lower():
            continue
        accepted_spans = []
        for pattern in CAUSALITY_ESTABLISHMENT_RES:
            for match in pattern.finditer(line):
                if any(
                    match.start() < end and start < match.end()
                    for start, end in accepted_spans
                ):
                    continue  # already covered by another pattern
                accepted_spans.append((match.start(), match.end()))
                sentence_prefix = _sentence_prefix_at(
                    lines, lineno - 1, match.start()
                )
                if _is_causality_assertion_negated(
                    line, match.start(), match.end(), sentence_prefix
                ):
                    continue
                findings.append(
                    _finding(
                        "SCILANG",
                        rel_path,
                        f"line {lineno}: affirmative cleavage-causality "
                        f"establishment claim '{match.group(0)}' is a "
                        f"violation; a document-wide disclaimer does not "
                        f"excuse it, and neither does a negation in the "
                        f"suffix after the claim (only a negation directly "
                        f"preceding the claim in its own sentence, e.g. "
                        f"'does not establish cleavage causality', does); "
                        f"required hedged wording: "
                        f"'{REQUIRED_DISCLAIMER}'",
                    )
                )

    # (b3) Curated documents discussing candidate guide-dependent /
    #      mechanistic associations must carry the required disclaimer.
    #      Generated product reports are exempt (engine boilerplate).
    has_disclaimer = bool(DISCLAIMER_CLAUSE_RE.search(text))
    if not has_disclaimer and not _is_generated_product_report(rel_path):
        for lineno, line in enumerate(lines, start=1):
            low = line.lower()
            if not (
                "candidate" in low or "guide" in low or "mechanistic" in low
            ):
                continue
            match = GUIDE_ASSOCIATION_DISCUSSION_RE.search(line)
            if match:
                findings.append(
                    _finding(
                        "SCILANG",
                        rel_path,
                        f"line {lineno}: document discusses candidate "
                        f"guide-dependent associations ('{match.group(0)}') "
                        f"but lacks the required disclaimer "
                        f"'{REQUIRED_DISCLAIMER}'",
                    )
                )
                break  # one finding per document is sufficient signal

    # (c) Rule 6 / Rule 8 wording violations, line by line.
    for lineno, line in enumerate(lines, start=1):
        low = line.lower()
        if "ncbi" in low:
            for pattern in RULE6_RES:
                match = pattern.search(line)
                if match and not NEGATION_MARKERS_RE.search(line):
                    findings.append(
                        _finding(
                            "SCILANG",
                            rel_path,
                            f"line {lineno}: Rule 6 violation "
                            f"'{match.group(0)}' (starter-strain WGS is "
                            f"mandatory; NCBI-reference differences must "
                            f"not be reported as unintended edits)",
                        )
                    )
        if "cast" in low or "is110" in low:
            if RULE8_EDITOR_RE.search(line) and RULE8_NO_MUTAGENESIS_RE.search(line):
                if RULE8_PROHIBITION_GUARD_RE.search(line):
                    continue
                findings.append(
                    _finding(
                        "SCILANG",
                        rel_path,
                        f"line {lineno}: Rule 8 violation: CAST/IS110 "
                        f"genome-wide mutagenesis claim "
                        f"(hypothesis columns are optional annotations, "
                        f"not facts)",
                    )
                )
    return findings


# ---------------------------------------------------------------------------
# Gate 2: evidence-chain traceability
# ---------------------------------------------------------------------------


def check_evidence_chain_traceability(repo_root: Path) -> list[str]:
    """Audit product TSVs and manifests under ``benchmark/``.

    Every variant row must carry a valid ``DV1_<32 hex>`` event identity or a
    valid locus reference, missing data must be the clean literal ``NA`` (or
    documented manifest sentinels), site references must resolve inside the
    sibling ``offtarget_sites.tsv``, event references must resolve inside
    the sibling ``post_edit_variants.tsv``, and no fabricated identifiers
    or unvalidated scores may appear.  Fails closed when zero evidence
    tables can be scanned.
    """
    repo_root = Path(repo_root)
    findings: list[str] = []
    sibling_cache: dict[tuple[Path, str], set[str] | None] = {}
    sibling_event_cache: dict[Path, dict[str, dict[str, str]] | None] = {}
    gd_cache: dict[Path, set[str] | None] = {}
    audited_reports: set[Path] = set()
    scanned = 0

    for abs_path, rel_path in _walk_files(
        repo_root, "benchmark", {".tsv", ".csv"}, findings
    ):
        rows = _read_tsv_rows(abs_path)
        if rows is None:
            findings.append(
                _finding(
                    "FAILCLOSED",
                    rel_path,
                    "in-scope evidence table is unreadable; the audit "
                    "refuses to silently skip it",
                )
            )
            continue
        name = Path(rel_path).name
        if not rows:
            # Codex Sol Round 2: a product table with no header row at
            # all (0 bytes / empty) is a failure, never a silent skip.
            if name in PRODUCT_FILENAMES:
                findings.append(
                    _finding(
                        "EVIDCHAIN",
                        rel_path,
                        "product table is empty (0 bytes / no header row); "
                        "product tables must at least carry their header "
                        "schema, and an empty table cannot be certified as "
                        "clean",
                    )
                )
            continue
        scanned += 1
        parent = abs_path.parent
        if name in PRODUCT_FILENAMES:
            sibling_sites = _sibling_ids(
                parent, SITES_SIBLING_NAME, "site_id", sibling_cache
            )
            sibling_events = _sibling_event_records(
                parent, sibling_event_cache
            )
            gd_mutations = _load_parent_gd_mutations(parent, gd_cache)
            findings.extend(
                _audit_product_table(
                    rel_path, name, rows, parent, sibling_sites, sibling_events, gd_mutations
                )
            )
            if parent not in audited_reports:
                audited_reports.add(parent)
                report_file = parent / "report.md"
                if report_file.is_file():
                    findings.extend(
                        _audit_report_findings(repo_root, report_file, sibling_events)
                    )
        else:
            findings.extend(_audit_generic_event_columns(rel_path, rows))
    # Ensure any remaining report.md files under benchmark/ are audited even if directory has no TSVs
    benchmark_dir = repo_root / "benchmark"
    if benchmark_dir.is_dir():
        for report_file in benchmark_dir.glob("**/report.md"):
            if report_file.is_file() and report_file.parent not in audited_reports:
                audited_reports.add(report_file.parent)
                sibling_events = _sibling_event_records(
                    report_file.parent, sibling_event_cache
                )
                findings.extend(
                    _audit_report_findings(repo_root, report_file, sibling_events)
                )
    if scanned == 0:
        findings.append(
            _finding(
                "EVIDCHAIN",
                "benchmark/",
                "zero evidence tables found (fail closed): an empty or "
                "missing scan scope cannot be certified as clean",
            )
        )
    return sorted(set(findings))


def _sibling_ids(
    parent: Path,
    sibling_name: str,
    column: str,
    cache: dict[tuple[Path, str], set[str] | None],
    value_re: re.Pattern[str] | None = None,
) -> set[str] | None:
    """Return ids declared by ``sibling_name`` in ``parent``.

    Returns ``None`` when the sibling file is absent (references against a
    missing sibling are orphan failures); an empty set means the sibling
    exists but declares no ids of that column.
    """
    key = (parent, sibling_name)
    if key in cache:
        return cache[key]
    sibling_path = parent / sibling_name
    ids: set[str] | None = None
    if sibling_path.is_file():
        ids = set()
        rows = _read_tsv_rows(sibling_path)
        if rows:
            header = rows[0]
            if column in header:
                idx = header.index(column)
                for row in rows[1:]:
                    if idx < len(row):
                        value = row[idx]
                        if value_re is None or value_re.match(value):
                            ids.add(value)
    cache[key] = ids
    return ids


EVIDENCE_TOKEN_RE = re.compile(r"^(RA|MC|JC)=(1|0|NA)$")
GD_MUTATION_TYPES = {"SNP", "SUB", "DEL", "INS", "MOB", "AMP", "CON", "INV", "JC"}


def _validate_evidence_grammar(evidence: str) -> tuple[bool, bool, str | None]:
    """Validate evidence column against documented schema.

    Returns (is_valid_grammar, has_positive_support, error_detail).
    """
    if not evidence or evidence.strip() in ("", "NA", "NONE", "."):
        return False, False, "evidence cell is empty or missing"
    raw_tokens = evidence.split(";")
    seen_keys: set[str] = set()
    has_positive = False
    for tok in raw_tokens:
        tok = tok.strip()
        if not tok:
            return False, False, "empty token in evidence string (e.g. consecutive or trailing semicolon)"
        m = EVIDENCE_TOKEN_RE.match(tok)
        if not m:
            return (
                False,
                False,
                f"malformed evidence token '{tok}'; expected (RA|MC|JC)=(1|0|NA)",
            )
        key, val = m.group(1), m.group(2)
        if key in seen_keys:
            return False, False, f"duplicate evidence key '{key}'"
        seen_keys.add(key)
        if val == "1":
            has_positive = True
    if not has_positive:
        return True, False, "no positive supporting GD tokens (all are 0 or NA)"
    return True, True, None


def _load_parent_gd_mutations(
    parent: Path, gd_cache: dict[Path, set[tuple[str, str, str]] | None]
) -> set[tuple[str, str, str]] | None:
    """Load mutation locus identities from sibling *.gd files in parent.

    Returns None if no .gd file exists in parent. Only actual mutation records
    (not evidence records like RA/MC) are eligible.
    """
    if parent in gd_cache:
        return gd_cache[parent]
    gd_files = [f for f in parent.glob("*.gd") if f.is_file()]
    if not gd_files:
        gd_cache[parent] = None
        return None
    mutations: set[tuple[str, str, str]] = set()
    for gd_path in gd_files:
        try:
            content = gd_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for line in content.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 5:
                # parts[0]: type (SNP, INS, DEL, SUB, JC, RA, MC)
                gd_type = parts[0]
                if gd_type in GD_MUTATION_TYPES:
                    mutations.add((gd_type, parts[3], parts[4]))
                    if gd_type == "JC" and len(parts) >= 8:
                        mutations.add(("JC", parts[6], parts[7]))
    gd_cache[parent] = mutations
    return mutations


def _sibling_event_records(
    parent: Path,
    cache: dict[Path, dict[str, dict[str, str]] | None],
) -> dict[str, dict[str, str]] | None:
    """Return a mapping of event_id -> dict of record fields for post_edit_variants.tsv in parent.

    Returns None if post_edit_variants.tsv is absent.
    """
    if parent in cache:
        return cache[parent]
    sibling_path = parent / VARIANTS_SIBLING_NAME
    if not sibling_path.is_file():
        cache[parent] = None
        return None
    records: dict[str, dict[str, str]] = {}
    rows = _read_tsv_rows(sibling_path)
    if rows:
        header = rows[0]
        col_map = {c: i for i, c in enumerate(header)}
        if "event_id" in col_map:
            for row in rows[1:]:
                ev_id = row[col_map["event_id"]] if col_map["event_id"] < len(row) else ""
                if ev_id and DV1_RE.match(ev_id):
                    rec: dict[str, str] = {}
                    for fld in ("seq_id", "position", "end", "gd_type", "ref", "alt", "evidence"):
                        if fld in col_map and col_map[fld] < len(row):
                            rec[fld] = row[col_map[fld]]
                        else:
                            rec[fld] = ""
                    records[ev_id] = rec
    cache[parent] = records
    return records


def _audit_report_findings(
    repo_root: Path,
    report_path: Path,
    sibling_events: dict[str, dict[str, str]] | None,
) -> list[str]:
    """Audit report.md findings: each referenced DV1 event must resolve and match coordinates."""
    findings: list[str] = []
    rel_path = report_path.relative_to(repo_root).as_posix()
    try:
        text = report_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return findings
    events = sibling_events if sibling_events is not None else {}
    for lineno, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        for match in DV1_SEARCH_RE.finditer(line):
            ev_id = match.group(0)
            if ev_id not in events:
                findings.append(
                    _finding(
                        "EVIDCHAIN", rel_path,
                        f"line {lineno}: finding references event id '{ev_id}' which "
                        f"does not resolve in sibling {VARIANTS_SIBLING_NAME} (orphan event finding)",
                    )
                )
            else:
                decl = events[ev_id]
                locus_match = re.search(r"\b([A-Za-z0-9_.-]+):(\d+)\b", line)
                if locus_match:
                    locus_seq, locus_pos = locus_match.group(1), locus_match.group(2)
                    if decl.get("seq_id") and locus_seq != decl["seq_id"]:
                        findings.append(
                            _finding(
                                "EVIDCHAIN", rel_path,
                                f"line {lineno}: finding for event id '{ev_id}' at "
                                f"contig '{locus_seq}' is inconsistent with contig "
                                f"'{decl['seq_id']}' declared in {VARIANTS_SIBLING_NAME}",
                            )
                        )
                    if decl.get("position") and locus_pos != decl["position"]:
                        findings.append(
                            _finding(
                                "EVIDCHAIN", rel_path,
                                f"line {lineno}: finding for event id '{ev_id}' at "
                                f"position '{locus_pos}' is inconsistent with position "
                                f"'{decl['position']}' declared in {VARIANTS_SIBLING_NAME}",
                            )
                        )
    return findings


def _audit_product_table(
    rel_path: str,
    name: str,
    rows: list[list[str]],
    parent: Path,
    sibling_sites: set[str] | None,
    sibling_events: dict[str, dict[str, str]] | set[str] | None,
    gd_mutations: set[tuple[str, ...]] | set[str] | None = None,
) -> list[str]:
    findings: list[str] = []
    header = rows[0]
    data = rows[1:]
    col = {c: i for i, c in enumerate(header)}

    def cell(row: list[str], column: str) -> str:
        i = col.get(column)
        if i is None or i >= len(row):
            return ""
        return row[i]

    def check_site_reference(row_kind: str, site_id: str) -> None:
        findings.extend(
            _check_site_reference(rel_path, row_kind, site_id, sibling_sites)
        )

    def check_event_reference(
        row_kind: str, value: str, row: list[str] | None = None, *, allow_na: bool
    ) -> None:
        findings.extend(
            _check_event_reference(
                rel_path, row_kind, value, sibling_events, row, col, allow_na=allow_na
            )
        )

    def check_coordinates(row: list[str]) -> None:
        findings.extend(_check_coordinate_cells(rel_path, row, col))

    if name == VARIANTS_SIBLING_NAME:
        for required in ("variant_id", "seq_id", "position"):
            if required not in col:
                findings.append(
                    _finding("EVIDCHAIN", rel_path, f"missing required column '{required}'")
                )
        for row in data:
            if len(row) != len(header):
                findings.append(
                    _finding("EVIDCHAIN", rel_path, f"ragged row {row[:3]!r}")
                )
            variant_id = cell(row, "variant_id")
            if variant_id and not VARIANT_ID_RE.match(variant_id):
                findings.append(
                    _finding(
                        "EVIDCHAIN", rel_path,
                        f"variant_id '{variant_id}' is not a VAR_#### display label",
                    )
                )
            findings.extend(_check_locus_reference(rel_path, row, col))
            check_coordinates(row)
            if "event_id" in col:  # schema V2 only; self-declaring table
                ev_id = cell(row, "event_id")
                findings.extend(
                    _check_event_cell(
                        rel_path, "post_edit_variants row", ev_id,
                        allow_na=True,
                    )
                )
                if ev_id and DV1_RE.match(ev_id):
                    evidence_val = cell(row, "evidence")
                    is_valid, has_pos, err = _validate_evidence_grammar(evidence_val)
                    if not is_valid:
                        findings.append(
                            _finding(
                                "EVIDCHAIN", rel_path,
                                f"post_edit_variants row: event id '{ev_id}' has "
                                f"invalid evidence grammar: {err} (evidence='{evidence_val}')",
                            )
                        )
                    elif not has_pos:
                        findings.append(
                            _finding(
                                "EVIDCHAIN", rel_path,
                                f"post_edit_variants row: event id '{ev_id}' has "
                                f"no supporting GD evidence: {err} (evidence='{evidence_val}')",
                            )
                        )
                    elif gd_mutations is not None:
                        decl_type = cell(row, "gd_type")
                        decl_seq = cell(row, "seq_id")
                        decl_pos = cell(row, "position")
                        matched = (decl_type, decl_seq, decl_pos) in gd_mutations
                        if not matched:
                            findings.append(
                                _finding(
                                    "EVIDCHAIN", rel_path,
                                    f"post_edit_variants row: event id '{ev_id}' "
                                    f"({decl_type} at {decl_seq}:{decl_pos}) does not resolve "
                                    f"to any record in sibling GD file(s)",
                                )
                            )

    elif name == "unintended.tsv":
        for required in ("seq_id", "position"):
            if required not in col:
                findings.append(
                    _finding("EVIDCHAIN", rel_path, f"missing required column '{required}'")
                )
        for row in data:
            if len(row) != len(header):
                findings.append(
                    _finding("EVIDCHAIN", rel_path, f"ragged row {row[:3]!r}")
                )
            findings.extend(_check_locus_reference(rel_path, row, col))
            check_coordinates(row)
            if "event_id" in col:  # schema V2 only
                value = cell(row, "event_id")
                findings.extend(
                    _check_event_cell(
                        rel_path, "unintended row", value, allow_na=True
                    )
                )
                check_event_reference("unintended row", value, row, allow_na=True)

    elif name == "mutation_offtarget_links.tsv":
        for required in ("mutation_id", "site_id"):
            if required not in col:
                findings.append(
                    _finding("EVIDCHAIN", rel_path, f"missing required column '{required}'")
                )
        for row in data:
            if len(row) != len(header):
                findings.append(
                    _finding("EVIDCHAIN", rel_path, f"ragged row {row[:3]!r}")
                )
            check_site_reference("association row", cell(row, "site_id"))
            check_coordinates(row)
            mutation_id = cell(row, "mutation_id")
            if mutation_id and not LEGACY_MUTATION_ID_RE.match(mutation_id):
                findings.append(
                    _finding(
                        "EVIDCHAIN", rel_path,
                        f"mutation_id '{mutation_id}' is not a mut_<id> identifier",
                    )
                )
            if "event_id" in col:
                value = cell(row, "event_id")
                findings.extend(
                    _check_event_cell(
                        rel_path, "association row", value, allow_na=False
                    )
                )
                check_event_reference("association row", value, row, allow_na=False)

    elif name == SITES_SIBLING_NAME:
        if "site_id" not in col:
            findings.append(
                _finding("EVIDCHAIN", rel_path, "missing required column 'site_id'")
            )
        for row in data:
            if len(row) != len(header):
                findings.append(
                    _finding("EVIDCHAIN", rel_path, f"ragged row {row[:3]!r}")
                )
            site_id = cell(row, "site_id")
            if site_id and not SITE_ID_RE.match(site_id):
                findings.append(
                    _finding(
                        "EVIDCHAIN", rel_path,
                        f"site_id '{site_id}' is not a SITE_###### identifier",
                    )
                )
            check_coordinates(row)
            pam = cell(row, "pam")
            if pam in PAM_PLACEHOLDERS:
                findings.append(
                    _finding(
                        "EVIDCHAIN", rel_path,
                        f"PAM cell '{pam}' is a query-pattern placeholder; the "
                        f"observed PAM contract forbids it",
                    )
                )
            for score_col in ("cfd_score", "hsu_score"):
                value = cell(row, score_col)
                if value not in ("", "NA"):
                    findings.append(
                        _finding(
                            "EVIDCHAIN", rel_path,
                            f"{score_col} '{value}' is a score without a validated "
                            f"external oracle; it must remain NA "
                            f"(DISABLED_UNVALIDATED_ORACLE)",
                        )
                    )

    elif name == "edit_outcomes.tsv":
        for row in data:
            if len(row) != len(header):
                findings.append(
                    _finding("EVIDCHAIN", rel_path, f"ragged row {row[:3]!r}")
                )
            check_coordinates(row)
            for column in ("matched_event_ids_dv1", "unexpected_event_ids_dv1"):
                if column not in col:
                    continue
                value = cell(row, column)
                findings.extend(_check_dv1_list(rel_path, column, value))
                for token in re.split(r"[;,]", value):
                    token = token.strip()
                    if DV1_RE.match(token):
                        check_event_reference(
                            f"{column} token", token, row, allow_na=False
                        )
            if "matched_event_ids" in col:
                value = cell(row, "matched_event_ids")
                if value and value != "NONE":
                    for token in re.split(r"[;,]", value):
                        if not token.strip().isdigit():
                            findings.append(
                                _finding(
                                    "EVIDCHAIN", rel_path,
                                    f"matched_event_ids '{value}' is neither NONE "
                                    f"nor numeric GD ids",
                                )
                            )
                            break

    else:
        # Generic product table (e.g. mechanistic_associations.tsv):
        # police site/event reference integrity and coordinate sanity
        # without schema knowledge (Codex Sol Round 2: coordinate columns
        # are validated here too).
        for row in data:
            if len(row) != len(header):
                findings.append(
                    _finding("EVIDCHAIN", rel_path, f"ragged row {row[:3]!r}")
                )
            check_coordinates(row)
            if "site_id" in col:
                check_site_reference(f"{name} row", cell(row, "site_id"))
            if "event_id" in col:
                value = cell(row, "event_id")
                findings.extend(
                    _check_event_cell(rel_path, f"{name} row", value, allow_na=True)
                )
                check_event_reference(f"{name} row", value, row, allow_na=True)
            for column in ("matched_event_ids_dv1", "unexpected_event_ids_dv1"):
                if column in col:
                    value = cell(row, column)
                    findings.extend(_check_dv1_list(rel_path, column, value))
                    for token in re.split(r"[;,]", value):
                        token = token.strip()
                        if DV1_RE.match(token):
                            check_event_reference(
                                f"{column} token", token, row, allow_na=False
                            )
    return findings


def _check_site_reference(
    rel_path: str,
    row_kind: str,
    site_id: str,
    sibling_sites: set[str] | None,
) -> list[str]:
    """Site references must be well-formed and resolve in the sibling table."""
    if site_id == "SITE_UNKNOWN":
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"{row_kind}: fabricated 'SITE_UNKNOWN' site reference is "
                f"forbidden",
            )
        ]
    if not site_id:
        # Codex Sol Round 2: a blank site_id cell is an error, never a
        # silent pass.
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"{row_kind}: site_id cell is empty; blank site references "
                f"are forbidden — reference a resolving SITE_###### "
                f"identifier",
            )
        ]
    if not SITE_ID_RE.match(site_id):
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"{row_kind}: site_id '{site_id}' is not a SITE_###### "
                f"identifier",
            )
        ]
    if sibling_sites is None:
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"{row_kind}: site_id '{site_id}' references sibling "
                f"'{SITES_SIBLING_NAME}' but that table is missing in the "
                f"same directory (orphan reference)",
            )
        ]
    if site_id not in sibling_sites:
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"{row_kind}: site_id '{site_id}' does not resolve in sibling "
                f"{SITES_SIBLING_NAME}",
            )
        ]
    return []


def _check_event_reference(
    rel_path: str,
    row_kind: str,
    value: str,
    sibling_events: dict[str, dict[str, str]] | set[str] | None,
    row: list[str] | None = None,
    col: dict[str, int] | None = None,
    *,
    allow_na: bool,
) -> list[str]:
    """DV1 event references must resolve in the sibling variants table and be content-consistent."""
    if not DV1_RE.match(value):
        return []  # format handled by _check_event_cell
    if sibling_events is None:
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"{row_kind}: event id '{value}' references sibling "
                f"'{VARIANTS_SIBLING_NAME}' but that table is missing in "
                f"the same directory (orphan reference)",
            )
        ]
    if value not in sibling_events:
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"{row_kind}: event id '{value}' does not resolve in sibling "
                f"{VARIANTS_SIBLING_NAME} (orphan event reference)",
            )
        ]
    findings: list[str] = []
    if isinstance(sibling_events, dict) and row is not None and col is not None:
        decl = sibling_events[value]
        # Verify mutation coordinates consistency (seq_id, position)
        if "seq_id" in col:
            row_seq = row[col["seq_id"]] if col["seq_id"] < len(row) else ""
            if row_seq and decl.get("seq_id") and row_seq != decl["seq_id"]:
                findings.append(
                    _finding(
                        "EVIDCHAIN", rel_path,
                        f"{row_kind}: event id '{value}' declared on seq_id "
                        f"'{decl['seq_id']}' in {VARIANTS_SIBLING_NAME} is "
                        f"inconsistently referenced with seq_id '{row_seq}'",
                    )
                )
        pos_idx = col.get("position") if "position" in col else col.get("mutation_position")
        if pos_idx is not None:
            row_pos = row[pos_idx] if pos_idx < len(row) else ""
            if row_pos and decl.get("position") and row_pos != decl["position"]:
                col_name = "position" if "position" in col else "mutation_position"
                findings.append(
                    _finding(
                        "EVIDCHAIN", rel_path,
                        f"{row_kind}: event id '{value}' declared at position "
                        f"'{decl['position']}' in {VARIANTS_SIBLING_NAME} is "
                        f"inconsistently referenced with {col_name} '{row_pos}'",
                    )
                )
        if "end" in col:
            row_end = row[col["end"]] if col["end"] < len(row) else ""
            if row_end and decl.get("end") and row_end != decl["end"]:
                findings.append(
                    _finding(
                        "EVIDCHAIN", rel_path,
                        f"{row_kind}: event id '{value}' declared with end "
                        f"'{decl['end']}' in {VARIANTS_SIBLING_NAME} is "
                        f"inconsistently referenced with end '{row_end}'",
                    )
                )
        dna_re = re.compile(r"^[ACGTNacgtn]+$")
        decl_type = decl.get("gd_type", "")
        if "ref" in col:
            row_ref = row[col["ref"]] if col["ref"] < len(row) else ""
            decl_ref = decl.get("ref", "")
            if row_ref and decl_ref:
                should_match = (decl_type == "SNP") or bool(dna_re.match(row_ref) and dna_re.match(decl_ref))
                if should_match and row_ref != decl_ref:
                    findings.append(
                        _finding(
                            "EVIDCHAIN", rel_path,
                            f"{row_kind}: event id '{value}' declared with ref "
                            f"'{decl_ref}' in {VARIANTS_SIBLING_NAME} is "
                            f"inconsistently referenced with ref '{row_ref}'",
                        )
                    )
        if "alt" in col:
            row_alt = row[col["alt"]] if col["alt"] < len(row) else ""
            decl_alt = decl.get("alt", "")
            if row_alt and decl_alt:
                should_match = (decl_type == "SNP") or bool(dna_re.match(row_alt) and dna_re.match(decl_alt))
                if should_match and row_alt != decl_alt:
                    findings.append(
                        _finding(
                            "EVIDCHAIN", rel_path,
                            f"{row_kind}: event id '{value}' declared with alt "
                            f"'{decl_alt}' in {VARIANTS_SIBLING_NAME} is "
                            f"inconsistently referenced with alt '{row_alt}'",
                        )
                    )
        for type_col in ("gd_type", "mutation_type"):
            if type_col in col:
                row_type = row[col[type_col]] if col[type_col] < len(row) else ""
                if row_type and decl.get("gd_type") and row_type != decl["gd_type"]:
                    findings.append(
                        _finding(
                            "EVIDCHAIN", rel_path,
                            f"{row_kind}: event id '{value}' declared as {decl['gd_type']} "
                            f"in {VARIANTS_SIBLING_NAME} is inconsistently referenced "
                            f"with {type_col} '{row_type}'",
                        )
                    )
    return findings


def _check_locus_reference(
    rel_path: str, row: list[str], col: dict[str, int]
) -> list[str]:
    seq_id = row[col["seq_id"]] if col.get("seq_id") is not None and col["seq_id"] < len(row) else ""
    position = row[col["position"]] if col.get("position") is not None and col["position"] < len(row) else ""
    if not seq_id:
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"variant row {row[:3]!r} has an empty seq_id locus reference",
            )
        ]
    if not position.isdigit() or int(position) < 1:
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"variant row {row[:3]!r} has invalid 1-based position '{position}'",
            )
        ]
    if int(position) > MAX_BACTERIAL_COORD:
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"variant row {row[:3]!r} position {position} exceeds the "
                f"bacterial-contig sanity bound {MAX_BACTERIAL_COORD} "
                f"(absurd coordinate)",
            )
        ]
    return []


def _check_coordinate_cells(
    rel_path: str, row: list[str], col: dict[str, int]
) -> list[str]:
    """Validate product-table coordinate cells.

    Codex Sol Round 2: a coordinate cell must be the clean literal ``NA``
    or a positive integer within ``1 .. MAX_BACTERIAL_COORD``.  Any other
    value — negative ("-1"), non-numeric ("abc"), blank, or above-bound
    ("10000001") — is an error, never a silent skip.

    Codex Sol Round 3: ``mutation_position`` and ``side2_position`` are
    policed as coordinate columns too.  ``side2_position`` follows the
    schema's second-side contract: blank sentinels are permitted on
    non-JC rows (``gd_type`` present and != "JC"), a populated value
    must be a valid coordinate on any row (non-JC included), and a JC
    row must carry a real second-side coordinate — blank
    ``side2_position`` on a JC row is itself a failure.
    """
    findings: list[str] = []
    gd_type_index = col.get("gd_type")
    gd_type = (
        row[gd_type_index]
        if gd_type_index is not None and gd_type_index < len(row)
        else ""
    )
    for column in COORDINATE_COLUMN_NAMES:
        i = col.get(column)
        if i is None or i >= len(row):
            continue
        value = row[i]
        if column == "side2_position" and gd_type_index is not None:
            if gd_type != "JC" and value in SIDE2_BLANK_LITERALS:
                # Schema-permitted blank second-side field on a non-JC
                # row (single-side call).
                continue
            if gd_type == "JC" and value in SIDE2_BLANK_LITERALS:
                findings.append(
                    _finding(
                        "EVIDCHAIN", rel_path,
                        f"row {row[:3]!r} is a JC row but column "
                        f"'side2_position' is blank ('{value}'); JC rows "
                        f"must carry a real second-side coordinate within "
                        f"1..{MAX_BACTERIAL_COORD}",
                    )
                )
                continue
        if value == "NA":
            continue
        if not value.isdigit():
            findings.append(
                _finding(
                    "EVIDCHAIN", rel_path,
                    f"row {row[:3]!r} column '{column}' coordinate cell "
                    f"'{value}' is not the clean literal 'NA' or an "
                    f"integer within the bacterial-contig sanity bounds "
                    f"1..{MAX_BACTERIAL_COORD}",
                )
            )
            continue
        coordinate = int(value)
        if coordinate < 1 or coordinate > MAX_BACTERIAL_COORD:
            findings.append(
                _finding(
                    "EVIDCHAIN", rel_path,
                    f"row {row[:3]!r} column '{column}' coordinate "
                    f"'{value}' is outside the bacterial-contig sanity "
                    f"bounds 1..{MAX_BACTERIAL_COORD}",
                )
            )
    return findings


def _check_event_cell(
    rel_path: str, row_kind: str, value: str, *, allow_na: bool
) -> list[str]:
    if value == "":
        return [
            _finding(
                "EVIDCHAIN", rel_path,
                f"{row_kind} has an empty event_id cell; missing data must be "
                f"the literal 'NA', never blank",
            )
        ]
    if DV1_RE.match(value):
        return []
    if allow_na and value in PRODUCT_MISSING:
        return []
    return [
        _finding(
            "EVIDCHAIN", rel_path,
            f"{row_kind} event identity '{value}' is not a valid "
            f"DV1_<32 lowercase hex> id or clean 'NA'",
        )
    ]


def _check_dv1_list(rel_path: str, column: str, value: str) -> list[str]:
    if value == "" or value == "NONE":
        if value == "":
            return [
                _finding(
                    "EVIDCHAIN", rel_path,
                    f"{column} is empty; use the literal NONE for an empty list",
                )
            ]
        return []
    for token in re.split(r"[;,]", value):
        token = token.strip()
        if not DV1_RE.match(token):
            return [
                _finding(
                    "EVIDCHAIN", rel_path,
                    f"{column} token '{token}' is not a valid DV1 id or NONE",
                )
            ]
    return []


def _audit_generic_event_columns(rel_path: str, rows: list[list[str]]) -> list[str]:
    """Police event-identity columns in every non-product TSV/CSV table."""
    findings: list[str] = []
    header = rows[0]
    event_indexes: list[tuple[int, str]] = []
    for i, column in enumerate(header):
        if column in EVENT_COLUMN_EXACT or column.endswith(EVENT_COLUMN_SUFFIXES):
            event_indexes.append((i, column))
    if not event_indexes:
        return findings
    for row in rows[1:]:
        for i, column in event_indexes:
            value = row[i] if i < len(row) else ""
            if DV1_RE.match(value) or value in MANIFEST_MISSING:
                continue
            findings.append(
                _finding(
                    "EVIDCHAIN", rel_path,
                    f"column '{column}' value '{value}' is neither a valid "
                    f"DV1_<32 hex> id nor a clean missing sentinel "
                    f"{sorted(MANIFEST_MISSING)}",
                )
            )
    return findings


# ---------------------------------------------------------------------------
# Gate 3: quarantine and G1 hold integrity
# ---------------------------------------------------------------------------


def check_quarantine_and_hold_integrity(repo_root: Path) -> list[str]:
    """Validate the frozen quarantine registry and the G1 strict-gate hold."""
    repo_root = Path(repo_root)
    findings: list[str] = []
    findings.extend(_check_validation_registry(repo_root))
    findings.extend(_check_wave11e_hold(repo_root))
    return sorted(set(findings))


def _check_validation_registry(repo_root: Path) -> list[str]:
    findings: list[str] = []
    registry_path = repo_root / REGISTRY_RELPATH
    if not registry_path.is_file():
        return [
            _finding(
                "QUARANTINE",
                REGISTRY_RELPATH.as_posix(),
                "validation registry is missing; release evidence cannot be "
                "verified (fail closed)",
            )
        ]
    rows = _read_tsv_rows(registry_path)
    if rows is None:
        return [
            _finding(
                "QUARANTINE",
                REGISTRY_RELPATH.as_posix(),
                "validation registry is unreadable (fail closed)",
            )
        ]
    if not rows:
        return [
            _finding(
                "QUARANTINE", REGISTRY_RELPATH.as_posix(), "registry is empty"
            )
        ]
    header = rows[0]
    col = {c: i for i, c in enumerate(header)}
    for required in ("dataset", "status", *REGISTRY_FLAG_COLUMNS):
        if required not in col:
            findings.append(
                _finding(
                    "QUARANTINE",
                    REGISTRY_RELPATH.as_posix(),
                    f"missing required column '{required}'",
                )
            )
    if {"dataset", "status"} - set(col):
        return findings  # cannot evaluate further

    seen: dict[str, int] = {}
    for row in rows[1:]:
        dataset = row[col["dataset"]] if col["dataset"] < len(row) else ""
        status = row[col["status"]] if col["status"] < len(row) else ""
        if dataset not in QUARANTINED_DATASETS:
            continue
        seen[dataset] = seen.get(dataset, 0) + 1
        if status != "PENDING":
            findings.append(
                _finding(
                    "QUARANTINE",
                    REGISTRY_RELPATH.as_posix(),
                    f"dataset '{dataset}' status is '{status}'; it must remain "
                    f"'PENDING' until versioned calls and comparison results "
                    f"exist",
                )
            )
        for flag in REGISTRY_FLAG_COLUMNS:
            i = col.get(flag)
            value = row[i] if i is not None and i < len(row) else ""
            if value.strip().lower() not in FALSE_FLAG_LITERALS:
                findings.append(
                    _finding(
                        "QUARANTINE",
                        REGISTRY_RELPATH.as_posix(),
                        f"dataset '{dataset}' has {flag}={value}; quarantined "
                        f"datasets must keep all validation flags false "
                        f"(only {sorted(nonempty for nonempty in FALSE_FLAG_LITERALS if nonempty)} "
                        f"are accepted as false; any other value — 'true', "
                        f"'yes', '1', 't', 'complete', ... — is rejected)",
                    )
                )
    for dataset in QUARANTINED_DATASETS:
        if dataset not in seen:
            findings.append(
                _finding(
                    "QUARANTINE",
                    REGISTRY_RELPATH.as_posix(),
                    f"required dataset row '{dataset}' is missing",
                )
            )
    return findings


def _check_wave11e_hold(repo_root: Path) -> list[str]:
    findings: list[str] = []
    benchmark_dir = repo_root / "benchmark"
    hold_paths: list[Path] = []
    if benchmark_dir.is_dir():
        for root, dirs, files in os.walk(benchmark_dir):
            dirs.sort()
            files.sort()
            rel_root_posix = Path(root).relative_to(repo_root).as_posix()
            if _under_skipped_fragment(rel_root_posix) and rel_root_posix != "benchmark":
                continue
            dirs[:] = [
                d for d in dirs if not _is_skipped_dir(rel_root_posix, d)
            ]
            if HOLD_FILENAME in files:
                hold_paths.append(Path(root) / HOLD_FILENAME)
    if not hold_paths:
        return [
            _finding(
                "QUARANTINE",
                f"benchmark/**/{HOLD_FILENAME}",
                "G1 hold document is missing; the frozen strict-gate state "
                "cannot be verified (fail closed)",
            )
        ]
    for hold_path in hold_paths:
        rel_path = hold_path.relative_to(repo_root).as_posix()
        text = _read_text(hold_path)
        if text is None:
            findings.append(
                _finding("QUARANTINE", rel_path, "hold document is unreadable")
            )
            continue
        if not G1_BASELINE_RE.search(text):
            findings.append(
                _finding(
                    "QUARANTINE",
                    rel_path,
                    "frozen G1 strict-gate baseline 'matched 10, OVER 24, "
                    "UNDER 8' is absent or altered",
                )
            )
        if not G1_NONPASS_MARKER_RE.search(text):
            findings.append(
                _finding(
                    "QUARANTINE",
                    rel_path,
                    "G1 non-pass marker (PROJECT_ACCEPTANCE_PENDING / "
                    "PASS=false) is missing",
                )
            )
        for match in G1_PASS_CLAIM_RE.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            findings.append(
                _finding(
                    "QUARANTINE",
                    rel_path,
                    f"line {line}: G1 pass claim '{match.group(0)}' contradicts "
                    f"the frozen PASS=false hold",
                )
            )
    return findings


# ---------------------------------------------------------------------------
# Aggregate runner and CLI
# ---------------------------------------------------------------------------


def run_all_release_gates(repo_root: Path) -> tuple[bool, dict[str, list[str]]]:
    """Run every release gate; return ``(all_passed, findings_by_gate)``."""
    repo_root = Path(repo_root)
    results: dict[str, list[str]] = {
        GATE_SCIENTIFIC_LANGUAGE: check_scientific_language(repo_root),
        GATE_EVIDENCE_CHAIN: check_evidence_chain_traceability(repo_root),
        GATE_QUARANTINE_HOLD: check_quarantine_and_hold_integrity(repo_root),
    }
    return all(not v for v in results.values()), results


def _render_markdown(repo_root: Path, results: dict[str, list[str]]) -> str:
    all_pass = all(not v for v in results.values())
    lines: list[str] = []
    lines.append("# ProkaDiff Release Gate Audit")
    lines.append("")
    lines.append(f"- Repository: `{repo_root}`")
    lines.append(f"- Gates run: {len(GATE_ORDER)}")
    lines.append("")
    lines.append("| Gate | Status | Findings |")
    lines.append("| --- | --- | --- |")
    for gate in GATE_ORDER:
        findings = results.get(gate, [])
        status = "PASS" if not findings else "FAIL"
        lines.append(f"| {gate} | {status} | {len(findings)} |")
    lines.append("")
    lines.append(f"**Overall: {'PASS' if all_pass else 'FAIL'}**")
    lines.append("")
    for gate in GATE_ORDER:
        findings = results.get(gate, [])
        if findings:
            lines.append(f"## Findings — {gate}")
            lines.append("")
            for item in findings:
                lines.append(f"- {item}")
            lines.append("")
    return "\n".join(lines)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="audit_release_gates",
        description=(
            "Deterministic ProkaDiff release-gate audit (scientific language, "
            "evidence-chain traceability, quarantine/G1 hold integrity)."
        ),
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="ProkaDiff repository root to audit (default: current directory)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="optional path to also write the markdown report",
    )
    args = parser.parse_args(argv)

    repo_root = Path(args.repo_root).resolve()
    if not repo_root.is_dir():
        print(f"error: repository root is not a directory: {repo_root}", file=sys.stderr)
        return 2

    passed, results = run_all_release_gates(repo_root)
    report = _render_markdown(repo_root, results)
    print(report)

    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(report + "\n", encoding="utf-8")

    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
