#!/usr/bin/env python3
"""
benchmark/real_world/scripts/aggregate_results.py
Aggregates benchmark results across datasets into a standardized summary.tsv.
"""

import sys
import os
import glob
import argparse


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

        # Parse metrics
        ds_metrics = {}
        sample_count = 0
        for mf in metric_files:
            m = parse_metric_tsv(mf)
            ds_metrics.update(m)
            sample_count += 1

        row = {
            "dataset": ds,
            "samples": str(sample_count if sample_count > 0 else 1),
            "variant_precision": ds_metrics.get("variant_precision", "N/A"),
            "variant_recall": ds_metrics.get("variant_recall", "N/A"),
            "snp_precision": ds_metrics.get("snp_precision", "N/A"),
            "snp_recall": ds_metrics.get("snp_recall", "N/A"),
            "indel_precision": ds_metrics.get("indel_precision", "N/A"),
            "indel_recall": ds_metrics.get("indel_recall", "N/A"),
            "mob_precision": ds_metrics.get("mob_precision", "N/A"),
            "mob_recall": ds_metrics.get("mob_recall", "N/A"),
            "jc_precision": ds_metrics.get("jc_precision", "N/A"),
            "jc_recall": ds_metrics.get("jc_recall", "N/A"),
            "false_complete": ds_metrics.get("false_complete", "0"),
            "report_p0_errors": ds_metrics.get("report_p0_errors", "0"),
            "median_runtime_sec": ds_metrics.get("median_runtime_sec", "N/A"),
            "median_rss_mb": ds_metrics.get("median_rss_mb", "N/A"),
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
