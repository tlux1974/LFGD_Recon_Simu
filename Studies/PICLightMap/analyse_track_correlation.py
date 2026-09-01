#!/usr/bin/env python3
"""Integrate per-fibre light-map response along 10 mm axis-parallel tracks."""

import argparse
import csv
import math
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/lightmap-matplotlib")
import matplotlib.pyplot as plt
import ROOT

from analyse_fibre_correlation import AXES, TRANSVERSE, fibre_direction, load_response


def analyse(filename, output_dir, label):
    root_file = ROOT.TFile.Open(str(filename), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")
    centre = root_file.Get("cubeCentre")
    global_centre = (float(centre.X()), float(centre.Y()), float(centre.Z()))
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    containment_rows = []
    figure, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    containment_figure, containment_axes = plt.subplots(
        1, 3, figsize=(17, 5), constrained_layout=True)
    colours = ("#4c78a8", "#f58518")

    for plot_axis, track_name in zip(axes, "xyz"):
        track_axis = AXES[track_name]
        source_positions = []
        integrated = defaultdict(float)
        for coordinate in [value - 4.5 for value in range(10)]:
            relative = [0.5, 0.5, 0.5]
            relative[track_axis] = coordinate
            source = tuple(global_centre[index] + relative[index]
                           for index in range(3))
            _, bin_number, response = load_response(root_file, source)
            source_positions.append((source, bin_number))
            for fibre_position, fraction in response:
                direction = fibre_direction(fibre_position)
                if direction is None:
                    continue
                axes_for_fibre = TRANSVERSE["xyz"[direction]]
                key = (direction,) + tuple(fibre_position[axis]
                                           for axis in axes_for_fibre)
                integrated[key] += fraction

        # Shortest distance between the (effectively infinite) central track
        # axis and each fibre axis.  For perpendicular lines this is the
        # separation along their one common transverse axis; for the parallel
        # fibre family it is the usual two-dimensional transverse distance.
        midpoint = (0.5, 0.5, 0.5)
        response_by_radius = defaultdict(float)
        for key, response_sum in integrated.items():
            fibre_axis = key[0]
            transverse_axes = TRANSVERSE["xyz"[fibre_axis]]
            fixed = dict(zip(transverse_axes, key[1:]))
            if fibre_axis == track_axis:
                radius = math.sqrt(sum(
                    (fixed[axis] - midpoint[axis]) ** 2
                    for axis in transverse_axes))
            else:
                common_axis = next(axis for axis in range(3)
                                   if axis not in (fibre_axis, track_axis))
                radius = abs(fixed[common_axis] - midpoint[common_axis])
            response_by_radius[radius] += response_sum

        total_response = sum(response_by_radius.values())
        cumulative = 0.0
        r99 = None
        radial_points = []
        for radius in sorted(response_by_radius):
            shell_response = response_by_radius[radius]
            cumulative += shell_response
            fraction = cumulative / total_response
            radial_points.append((radius, shell_response, cumulative, fraction))
            if r99 is None and fraction >= 0.99:
                r99 = radius
        containment_rows.extend({
            "track_axis": track_name,
            "radius_mm": radius,
            "shell_integrated_fraction": shell,
            "cumulative_integrated_fraction": cumulative_value,
            "cumulative_response_fraction": fraction,
            "r99_mm": r99,
        } for radius, shell, cumulative_value, fraction in radial_points)
        containment_axis = containment_axes["xyz".index(track_name)]
        containment_axis.step(
            [point[0] for point in radial_points],
            [point[3] for point in radial_points], where="post",
            color="#4c78a8", linewidth=2)
        containment_axis.axhline(.99, color="#e45756", linestyle="--",
                                 label="99% containment")
        containment_axis.axvline(r99, color="black", linestyle=":",
                                 label=f"R99 = {r99:.2f} mm")
        containment_axis.set(
            title=f"Track along {track_name}",
            xlabel="Shortest distance from track to fibre axis [mm]",
            ylabel="Cumulative recorded fibre response", ylim=(0, 1.015))
        containment_axis.grid(alpha=.2)
        containment_axis.legend(fontsize=8)

        fibre_axes = [axis for axis in range(3) if axis != track_axis]
        family_curves = []
        for family_index, fibre_axis in enumerate(fibre_axes):
            transverse_axes = TRANSVERSE["xyz"[fibre_axis]]
            fibres = {key[1:]: value for key, value in integrated.items()
                      if key[0] == fibre_axis}
            # Reference the fibre whose axis is closest to the track midpoint.
            midpoint = (0.5, 0.5, 0.5)
            nearest_key = min(
                fibres,
                key=lambda key: math.sqrt(sum(
                    (key[index] - midpoint[axis]) ** 2
                    for index, axis in enumerate(transverse_axes))))
            moving_index = transverse_axes.index(track_axis)
            reference = fibres[nearest_key]
            curve = {}
            for key, response_sum in fibres.items():
                if any(abs(key[index] - nearest_key[index]) > 1e-6
                       for index in range(2) if index != moving_index):
                    continue
                offset = (key[moving_index] - nearest_key[moving_index]) / 10.0
                if abs(offset - round(offset)) <= 1e-6:
                    curve[int(round(offset))] = response_sum
            family_curves.append((fibre_axis, nearest_key, reference, curve))
            offsets = sorted(curve)
            plot_axis.plot(offsets, [curve[n] / reference for n in offsets],
                           "o-", color=colours[family_index],
                           label=f"{'xyz'[fibre_axis].upper()}-directed fibres")

        offsets = sorted(set(family_curves[0][3]) | set(family_curves[1][3]))
        pair_reference = family_curves[0][2] + family_curves[1][2]
        combined = {}
        for offset in offsets:
            combined[offset] = sum(family[3].get(offset, 0.0)
                                   for family in family_curves) / pair_reference
        plot_axis.plot(offsets, [combined[n] for n in offsets], "s--",
                       color="black", linewidth=2, label="Combined pair")
        plot_axis.axhline(1.0, color="0.5", linestyle=":")
        plot_axis.set(title=f"Track along {track_name}",
                      xlabel="Fibre offset n [10 mm pitch]",
                      ylabel="Integrated response / reference", yscale="log")
        plot_axis.grid(alpha=.2, which="both")
        plot_axis.legend(fontsize=8)

        for offset in offsets:
            record = {"track_axis": track_name, "offset_n": offset,
                      "offset_mm": 10.0 * offset, "subvoxels": 10}
            for index, (fibre_axis, key, reference, curve) in enumerate(
                    family_curves, 1):
                response_sum = curve.get(offset, 0.0)
                record[f"fibre{index}_direction"] = "xyz"[fibre_axis]
                record[f"fibre{index}_reference_index"] = str(key)
                record[f"fibre{index}_reference_integrated_fraction"] = reference
                record[f"fibre{index}_integrated_fraction"] = response_sum
                record[f"fibre{index}_ratio"] = response_sum / reference
            record["combined_pair_ratio"] = combined[offset]
            rows.append(record)

    root_file.Close()
    figure.suptitle(
        f"{label}\n10-subvoxel tracks; fixed transverse coordinates = +0.5 mm")
    figure.savefig(output_dir / "integrated_track_fibre_correlation.png", dpi=170)
    plt.close(figure)
    containment_figure.suptitle(
        f"{label}\ntransverse containment of recorded fibre response")
    containment_figure.savefig(
        output_dir / "integrated_track_transverse_containment.png", dpi=170)
    plt.close(containment_figure)
    with (output_dir / "integrated_track_fibre_correlation.csv").open(
            "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (output_dir / "integrated_track_transverse_containment.csv").open(
            "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(containment_rows[0]))
        writer.writeheader()
        writer.writerows(containment_rows)
    print("Integrated ten subvoxels at -4.5, ..., +4.5 mm for each track axis")
    print(f"Wrote {output_dir / 'integrated_track_fibre_correlation.csv'}")
    print(f"Wrote {output_dir / 'integrated_track_fibre_correlation.png'}")
    for track_name in "xyz":
        r99 = next(row["r99_mm"] for row in containment_rows
                   if row["track_axis"] == track_name)
        print(f"  track {track_name}: R99 = {r99:.6g} mm")
    print(f"Wrote {output_dir / 'integrated_track_transverse_containment.csv'}")
    print(f"Wrote {output_dir / 'integrated_track_transverse_containment.png'}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="Merged ROOT light map")
    parser.add_argument("--label")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    analyse(args.input,
            args.output_dir or Path("lightmap_analysis") / args.input.stem,
            args.label or args.input.stem)


if __name__ == "__main__":
    main()
