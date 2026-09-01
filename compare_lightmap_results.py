#!/usr/bin/env python3
"""Collect standard-reconstruction and GlobalFit metrics across map branches."""

import argparse
import csv
import os
import statistics
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/lightmap-matplotlib")
import matplotlib.pyplot as plt
import ROOT


def mean(values):
    return statistics.fmean(values) if values else float("nan")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("comparison_dir", type=Path)
    args = parser.parse_args()
    rows = []
    for branch in sorted(path for path in args.comparison_dir.iterdir()
                         if path.is_dir()):
        reconstruction = branch / "plots" / "reconstruction_efficiency.csv"
        analysis = branch / "global_fit_analysis.root"
        fit_file = branch / "global_fit.root"
        if not reconstruction.is_file() or not analysis.is_file() or not fit_file.is_file():
            continue
        with reconstruction.open(newline="") as stream:
            standard = list(csv.DictReader(stream))
        root_analysis = ROOT.TFile.Open(str(analysis), "READ")
        directions = root_analysis.Get("direction_comparison")
        angles = [float(entry.mc_fit_angle) for entry in directions]
        root_analysis.Close()
        root_fit = ROOT.TFile.Open(str(fit_file), "READ")
        fits = root_fit.Get("global_fit")
        statuses = [int(entry.status) for entry in fits]
        root_fit.Close()
        rows.append({
            "tag": branch.name,
            "events": len(standard),
            "mean_hit_voxel_efficiency": mean([
                float(row["hit_voxel_efficiency"]) for row in standard]),
            "mean_node_voxel_efficiency": mean([
                float(row["node_voxel_efficiency"]) for row in standard]),
            "mean_reco_offtrack_voxels": mean([
                float(row["reco_offtrack_voxels"]) for row in standard]),
            "mean_reco_unique_voxels": mean([
                float(row["reco_unique_voxels"]) for row in standard]),
            "global_fit_entries": len(statuses),
            "global_fit_converged_fraction": (
                sum(status == 0 for status in statuses) / len(statuses)
                if statuses else float("nan")),
            "global_fit_mean_angle_deg": mean(angles),
            "global_fit_median_angle_deg": (
                statistics.median(angles) if angles else float("nan")),
        })
    if not rows:
        raise RuntimeError("No complete comparison branches were found")
    output = args.comparison_dir / "comparison_summary.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    labels = [row["tag"] for row in rows]
    figure, axes = plt.subplots(1, 3, figsize=(max(15, 2.3 * len(rows)), 5),
                               constrained_layout=True)
    axes[0].bar(labels, [row["mean_hit_voxel_efficiency"] for row in rows])
    axes[0].set(ylabel="Mean efficiency", title="Standard reco: matched MC voxels")
    axes[1].bar(labels, [row["mean_reco_offtrack_voxels"] for row in rows])
    axes[1].set(ylabel="Mean voxels/event", title="Standard reco: off-track voxels")
    axes[2].bar(labels, [row["global_fit_median_angle_deg"] for row in rows])
    axes[2].set(ylabel="Median angle [deg]", title="GlobalFit direction error")
    for axis in axes:
        axis.tick_params(axis="x", rotation=55)
        axis.grid(axis="y", alpha=.2)
    figure.savefig(args.comparison_dir / "comparison_summary.png", dpi=170)
    plt.close(figure)
    print(f"Wrote {output}")
    print(f"Wrote {args.comparison_dir / 'comparison_summary.png'}")


if __name__ == "__main__":
    main()
