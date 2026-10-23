#!/usr/bin/env python3
"""
benchmark/real_world/scripts/aggregate_results.py
Aggregates benchmark results across datasets into a standardized summary.tsv.

M7-P0b: comparator TSV keys are aliased to canonical keys when each metric
file is folded in, and values that conflict across files merge fail-closed to
NOT_EVALUABLE — a later file never replaces an earlier NOT_EVALUABLE sentinel
(or any conflicting concrete value) with its own number.
"""

import sys
import os
import glob
import argparse

NOT_EVALUABLE = "NOT_EVALUABLE"

# Comparator TSV keys are not the aggregator's canonical keys (M7-P0b):
# compare_variant_calls.py --tsv emits type-prefixed uppercase keys, and
# compare_mob_truth.py --tsv emits unprefixed precision/recall. Aliasing
# happens per file at fold time; the comparators' own TSV keys are unchanged.
_PREFIXED_ALIASES = {
    "MOB_precision": "mob_precision",
    "MOB_recall": "mob_recall",
    "JC_precision": "jc_precision",
    "JC_recall": "jc_recall",
    "SNP_precision": "snp_precision",
    "SNP_recall": "snp_recall",
}

# Only compare_mob_truth.py --tsv carries exact_fraction / median_error_bp;
# their presence is what licenses bare precision/recall to mean MOB metrics.
_MOB_MARKER_KEYS = ("exact_fraction", "median_error_bp")


def metric_value(ds_metrics, key, default="N/A"):
    """NOT_EVALUABLE propagates verbatim; it is never coerced to 0/1/0.0/1.0 (M7-P0b)."""
    value = ds_metrics.get(key, default)
    if value == NOT_EVALUABLE:
        return NOT_EVALUABLE
    return value


def parse_metric_tsv(filepath):
    metrics = {}
    if not os.path.exists(filepath):
        return metrics
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                metrics[parts[0].strip()] = parts[1].strip()
    return metrics


def merge_metric_value(existing, new):
    """Fail-closed conflict policy for one key across metric files (M7-P0b).

    Identical values pass through. If either side is NOT_EVALUABLE, or two
    concrete values disagree, the merged value is NOT_EVALUABLE; a later file
    never overwrites an earlier value with its own number.
    """
    if existing is None:
        return new
    if existing == new:
        return existing
    return NOT_EVALUABLE


def alias_metric_keys(metrics):
    """Adds canonical lowercase keys for comparator-shaped TSV keys (M7-P0b).

    Aliasing only adds (or fail-closes) canonical keys; the file's original
    keys are preserved. Bare precision/recall map onto mob_precision/mob_recall
    only when the same file also contains exact_fraction or median_error_bp
    (the compare_mob_truth.py --tsv output shape). A bare recall is never
    mapped onto variant_recall.
    """
    aliased = dict(metrics)
    for src, dst in _PREFIXED_ALIASES.items():
        if src in metrics:
            aliased[dst] = merge_metric_value(aliased.get(dst), metrics[src])
    if any(marker in metrics for marker in _MOB_MARKER_KEYS):
        for bare, mob_key in (("precision", "mob_precision"), ("recall", "mob_recall")):
            if bare in metrics:
                aliased[mob_key] = merge_metric_value(aliased.get(mob_key), metrics[bare])
    return aliased


def merge_metric_dicts(accumulated, incoming):
    """Merges one file's aliased metrics into the dataset accumulator,
    applying the fail-closed conflict policy per canonical key."""
    merged = dict(accumulated)
    for key, value in incoming.items():
        merged[key] = merge_metric_value(accumulated.get(key), value)
    return merged


def aggregate_dir(results_dir):
    """Aggregates all subdirectories under results_dir."""
    headers = [
        "dataset",
        "samples",
        "variant_precision",
        "variant_recall",
        "snp_precision",
        "snp_recall",
        "indel_precision",
        "indel_recall",
        "mob_precision",
        "mob_recall",
        "jc_precision",
        "jc_recall",
        "false_complete",
        "report_p0_errors",
        "median_runtime_sec",
        "median_rss_mb",
    ]

    rows = []
    datasets = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]

    for ds in sorted(datasets):
        ds_path = os.path.join(results_dir, ds)
        metric_files = glob.glob(os.path.join(ds_path, "**/*metrics.tsv"), recursive=True)
        summary_file = os.path.join(ds_path, "summary.tsv")
        if os.path.exists(summary_file):
            metric_files.append(summary_file)

        # Parse metrics, alias comparator-shaped keys to canonical keys, and
        # merge fail-closed so a later file can never overwrite an earlier
        # NOT_EVALUABLE sentinel with a number (M7-P0b).
        ds_metrics = {}
        sample_count = 0
        for mf in metric_files:
            m = alias_metric_keys(parse_metric_tsv(mf))
            ds_metrics = merge_metric_dicts(ds_metrics, m)
            sample_count += 1

        row = {
            "dataset": ds,
            "samples": str(sample_count if sample_count > 0 else 1),
            "variant_precision": metric_value(ds_metrics, "variant_precision"),
            "variant_recall": metric_value(ds_metrics, "variant_recall"),
            "snp_precision": metric_value(ds_metrics, "snp_precision"),
            "snp_recall": metric_value(ds_metrics, "snp_recall"),
            "indel_precision": metric_value(ds_metrics, "indel_precision"),
            "indel_recall": metric_value(ds_metrics, "indel_recall"),
            "mob_precision": metric_value(ds_metrics, "mob_precision"),
            "mob_recall": metric_value(ds_metrics, "mob_recall"),
            "jc_precision": metric_value(ds_metrics, "jc_precision"),
            "jc_recall": metric_value(ds_metrics, "jc_recall"),
            "false_complete": metric_value(ds_metrics, "false_complete", "0"),
            "report_p0_errors": metric_value(ds_metrics, "report_p0_errors", "0"),
            "median_runtime_sec": metric_value(ds_metrics, "median_runtime_sec"),
            "median_rss_mb": metric_value(ds_metrics, "median_rss_mb"),
        }
        rows.append(row)

    print("\t".join(headers))
    for r in rows:
        print("\t".join(r.get(h, "N/A") for h in headers))


def main():
    parser = argparse.ArgumentParser(description="Aggregate benchmark results into summary TSV")
    parser.add_argument("results_dir", help="Path to benchmark/real_world/results")
    args = parser.parse_args()

    aggregate_dir(args.results_dir)


if __name__ == "__main__":
    main()
