#!/usr/bin/env python3
"""
benchmark/real_world/scripts/compare_intended_edits.py
Evaluates intended edit outcome classifications against ground truth.
Monitors the critical release gate: false_complete == 0.
"""

import sys
import json
import argparse
import os


def parse_audit_json(audit_path):
    """Parses audit.json or intended edit TSV produced by ProkaDiff."""
    if not os.path.exists(audit_path):
        raise FileNotFoundError(f"Audit file not found: {audit_path}")

    with open(audit_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    # In ProkaDiff audit.json:
    # "intended_edits": [
    #    {
    #       "edit": { "seq_id": ..., "start": ..., "end": ... },
    #       "status": "Complete" | "Partial" | "Missing" | "UnexpectedStructure",
    #       "observed_pos": ...
    #    }
    # ]
    results = []
    for item in data.get("intended_edits", []):
        status = item.get("status")
        # Normalize status string
        if isinstance(status, dict):
            status_name = list(status.keys())[0]
        else:
            status_name = str(status)
        results.append({
            "seq_id": item.get("edit", {}).get("seq_id", ""),
            "start": item.get("edit", {}).get("start", 0),
            "end": item.get("edit", {}).get("end", 0),
            "status": status_name,
            "evidence": item.get("evidence", ""),
        })
    return results


def parse_curated_truth(truth_path):
    """
    Parses truth TSV file.
    Expected columns: truth_id, seq_id, start, end, expected_status
    or if bl21_curated_truth.tsv: sample_id or truth_id with notes indicating expected outcome.
    """
    records = []
    with open(truth_path, "r", encoding="utf-8") as f:
        headers = None
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if headers is None:
                headers = [h.strip() for h in parts]
                continue
            row = dict(zip(headers, parts))
            records.append(row)
    return records


def evaluate_outcomes(observed_edits, expected_status):
    """
    Compares observed edit status vs expected status.
    Returns status counts and checks for false_complete.
    """
    metrics = {
        "complete_true_positive": 0,
        "partial_true_positive": 0,
        "missing_true_positive": 0,
        "unexpected_structure_true_positive": 0,
        "false_complete": 0,
        "false_partial": 0,
        "false_missing": 0,
        "total_evaluated": len(observed_edits),
    }

    for edit in observed_edits:
        obs = edit["status"].upper()
        exp = expected_status.upper()

        if exp in ("COMPLETE", "PASS"):
            if obs in ("COMPLETE", "PASS"):
                metrics["complete_true_positive"] += 1
            elif obs in ("PARTIAL",):
                metrics["false_partial"] += 1
            else:
                metrics["false_missing"] += 1
        elif exp in ("UNEXPECTED_STRUCTURE", "UNEXPECTEDSTRUCTURE", "ABERRANT", "DELETION"):
            if obs in ("UNEXPECTEDSTRUCTURE", "UNEXPECTED_STRUCTURE", "ABERRANT"):
                metrics["unexpected_structure_true_positive"] += 1
            elif obs in ("COMPLETE", "PASS"):
                metrics["false_complete"] += 1
            else:
                # Other non-complete status
                pass
        elif exp in ("PARTIAL",):
            if obs in ("PARTIAL",):
                metrics["partial_true_positive"] += 1
            elif obs in ("COMPLETE", "PASS"):
                metrics["false_complete"] += 1
            else:
                metrics["false_missing"] += 1
        elif exp in ("MISSING",):
            if obs in ("MISSING",):
                metrics["missing_true_positive"] += 1
            elif obs in ("COMPLETE", "PASS"):
                metrics["false_complete"] += 1
            else:
                pass

    return metrics


def main():
    parser = argparse.ArgumentParser(description="Evaluate intended edit classification against truth")
    parser.add_argument("--audit-json", required=True, help="ProkaDiff audit.json file")
    parser.add_argument("--expected-status", required=True, help="Expected status (COMPLETE, PARTIAL, MISSING, UNEXPECTED_STRUCTURE)")
    parser.add_argument("--tsv", action="store_true", help="Output tab-separated format")
    args = parser.parse_args()

    edits = parse_audit_json(args.audit_json)
    metrics = evaluate_outcomes(edits, args.expected_status)

    is_gate_passed = (metrics["false_complete"] == 0)

    if args.tsv:
        print("metric\tvalue")
        for k, v in metrics.items():
            print(f"{k}\t{v}")
        print(f"gate_passed\t{str(is_gate_passed).lower()}")
    else:
        print("=== Intended Edit Outcome Evaluation ===")
        print(f"Total Edits Evaluated:               {metrics['total_evaluated']}")
        print(f"Complete True Positive:              {metrics['complete_true_positive']}")
        print(f"Unexpected Structure True Positive:  {metrics['unexpected_structure_true_positive']}")
        print(f"Partial True Positive:               {metrics['partial_true_positive']}")
        print(f"Missing True Positive:               {metrics['missing_true_positive']}")
        print(f"False Complete (CRITICAL GATE):       {metrics['false_complete']}")
        print(f"False Partial:                       {metrics['false_partial']}")
        print(f"False Missing:                       {metrics['false_missing']}")
        print("---------------------------------------")
        print(f"Release Gate (false_complete == 0):   {'PASSED' if is_gate_passed else 'FAILED'}")

    if not is_gate_passed:
        sys.exit(1)


if __name__ == "__main__":
    main()
