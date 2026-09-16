#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("layer2_gd_compare.py")
SPEC = importlib.util.spec_from_file_location("layer2_gd_compare", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {MODULE_PATH}")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def row(kind: str, *fields: str) -> dict[str, object]:
    return {"kind": kind, "fields": list(fields), "line": ""}


class ComparatorTests(unittest.TestCase):
    def test_red_annotations_do_not_change_mutation_identity(self) -> None:
        prokadiff = [row("INS", "REL606", "3875632", "T")]
        breseq = [
            row(
                "INS",
                "REL606",
                "3875632",
                "T",
                "repeat_length=1",
                "repeat_new_copies=8",
            )
        ]

        result = MODULE.compare(prokadiff, breseq, tol=5)

        self.assertEqual(result["over_red_n"], 0)
        self.assertEqual(result["under_red_n"], 0)

    def test_red_allele_difference_remains_a_failure(self) -> None:
        prokadiff = [row("DEL", "REL606", "4126706", "1")]
        breseq = [row("DEL", "REL606", "4126706", "2")]

        result = MODULE.compare(prokadiff, breseq, tol=5)

        self.assertEqual(result["over_red_n"], 1)
        self.assertEqual(result["under_red_n"], 1)

    def test_jc_still_uses_tolerant_coordinate_matching(self) -> None:
        prokadiff = [
            row("JC", "REL606", "599", "1", "REL606", "2558", "-1")
        ]
        breseq = [
            row(
                "JC",
                "REL606",
                "602",
                "1",
                "REL606",
                "2555",
                "-1",
                "reject=COVERAGE_EVENNESS_SKEW",
            )
        ]

        result = MODULE.compare(prokadiff, breseq, tol=5)

        self.assertEqual(result["jc_matched"], 1)
        self.assertEqual(result["over_jc_unmatched"], 0)
        self.assertEqual(result["under_jc_unmatched"], 0)


if __name__ == "__main__":
    unittest.main()
