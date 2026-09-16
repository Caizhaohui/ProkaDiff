#!/usr/bin/env python3
"""
benchmark/real_world/scripts/validate_reports.py
Audits human-readable reports (report.md) for section completeness,
scientific neutrality (no unverified causal claims), and P0 errors.
"""

import sys
import os
import re
import argparse


def validate_report(report_path, check_neutrality=True):
    """
    Validates report.md against publication-readiness guidelines.
    Returns error count by level (P0, P1, P2, P3) and a list of issues.
    """
    issues = []
    p0_count = 0
    p1_count = 0
    p2_count = 0
    p3_count = 0

    if not os.path.exists(report_path):
        issues.append(("P0", f"Report file not found: {report_path}"))
        return 1, 0, 0, 0, issues

    with open(report_path, "r", encoding="utf-8") as f:
        content = f.read()

    # 1. Check required headers / sections
    required_sections = [
        ("ProkaDiff Genome Audit Report", "P1", "Title missing"),
        ("Executive Summary", "P1", "Executive Summary section missing"),
        ("Provenance", "P1", "Provenance metadata section missing"),
    ]

    for sec, lvl, msg in required_sections:
        if sec.lower() not in content.lower():
            issues.append((lvl, msg))
            if lvl == "P0":
                p0_count += 1
            elif lvl == "P1":
                p1_count += 1

    # 2. Check for scientific neutrality / unverified causal assertions
    if check_neutrality:
        # Should not assert SOS / DSB-stress as established causal fact in default mode
        causal_red_flags = [
            (r"mutations? caused by SOS response", "Unverified causal attribution to SOS response"),
            (r"DSB-induced stress mutagenesis confirmed", "Unverified causal claim of stress mutagenesis"),
            (r"quality score: \d+", "Unauthorized genome quality score present"),
        ]
        for pattern, desc in causal_red_flags:
            if re.search(pattern, content, re.IGNORECASE):
                issues.append(("P0", f"Violation of scientific neutrality: {desc}"))
                p0_count += 1

    # 3. Check for stable IDs
    # Variants and candidate sites should use structured IDs (VAR_*, SITE_*)
    has_var_id = bool(re.search(r"VAR_\d+", content))
    has_site_id = bool(re.search(r"SITE_\d+", content))

    # 4. Check for provenance completeness
    if "Commit" in content or "Git" in content:
        if "UNKNOWN" in content and "dev" in content:
            issues.append(("P3", "Provenance commit hash is dev/UNKNOWN"))
            p3_count += 1

    return p0_count, p1_count, p2_count, p3_count, issues


def main():
    parser = argparse.ArgumentParser(description="Validate human-readable report.md for correctness and neutrality")
    parser.add_argument("report", help="Path to report.md")
    parser.add_argument("--strict", action="store_true", help="Fail on any P0 or P1 errors")
    args = parser.parse_args()

    p0, p1, p2, p3, issues = validate_report(args.report)

    print(f"=== Report Validation: {args.report} ===")
    print(f"P0 Errors (Critical Release Blockers): {p0}")
    print(f"P1 Errors (Missing Sections/Structural): {p1}")
    print(f"P2 Errors (Annotation Issues):          {p2}")
    print(f"P3 Errors (Formatting / Minor):         {p3}")

    if issues:
        print("\nIssues Found:")
        for lvl, msg in issues:
            print(f"  [{lvl}] {msg}")
    else:
        print("\nAll validation checks PASSED cleanly.")

    if args.strict and (p0 > 0 or p1 > 0):
        sys.exit(1)


if __name__ == "__main__":
    main()
