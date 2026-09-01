#!/usr/bin/env python3
"""Measure the central-subvoxel response falloff along perpendicular fibres."""

import argparse
import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/lightmap-matplotlib")
import matplotlib.pyplot as plt
import ROOT


AXES = {"x": 0, "y": 1, "z": 2}
# Fibre direction -> the two coordinates indexing that fibre.
TRANSVERSE = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}


def load_response(root_file, source):
    histogram = root_file.Get("Starting Vertex Distribution")
    centre = root_file.Get("cubeCentre")
    positions_dir = root_file.Get("hitPositions")
    fractions_dir = root_file.Get("lightFractions")
    if not all((histogram, centre, positions_dir, fractions_dir)):
        raise RuntimeError("Incomplete light map")
    relative = tuple(source[i] - (centre.X(), centre.Y(), centre.Z())[i]
                     for i in range(3))
    bin_number = histogram.FindBin(*relative)
    if histogram.GetBinContent(bin_number) <= 0:
        raise RuntimeError(f"Selected source {source} is not a populated map bin")
    positions = positions_dir.Get(f"Bin_{bin_number}_hit_positions")
    fractions = fractions_dir.Get(f"Bin_{bin_number}_light_fractions")
    if not positions or fractions is None:
        raise RuntimeError(f"Missing per-fibre response for bin {bin_number}")
    response = []
    for index in range(min(positions.GetEntries(), len(fractions))):
        position = positions.At(index)
        response.append(((float(position.X()), float(position.Y()),
                          float(position.Z())), float(fractions[index])))
    return relative, bin_number, response


def load_central_response(filename, source=(0.5, 30.5, 910.5)):
    root_file = ROOT.TFile.Open(str(filename), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")
    result = load_response(root_file, source)
    root_file.Close()
    return result


def fibre_direction(position, tolerance=1e-6):
    zeros = [index for index, value in enumerate(position)
             if abs(value) < tolerance]
    return zeros[0] if len(zeros) == 1 else None


def analyse(filename, output_dir, label):
    source, bin_number, response = load_central_response(filename)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    distance_plot_data = {}
    colours = ("#4c78a8", "#f58518")
    figure, plot_axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)

    for plot_axis, track_name in zip(plot_axes, "xyz"):
        track_axis = AXES[track_name]
        fibre_axes = [axis for axis in range(3) if axis != track_axis]
        family_curves = []
        for family_index, fibre_axis in enumerate(fibre_axes):
            transverse_axes = TRANSVERSE["xyz"[fibre_axis]]
            fibres = {}
            for position, fraction in response:
                if fibre_direction(position) != fibre_axis:
                    continue
                key = tuple(position[axis] for axis in transverse_axes)
                fibres[key] = fraction
            nearest_key = min(
                fibres,
                key=lambda key: math.sqrt(sum(
                    (key[index] - source[axis]) ** 2
                    for index, axis in enumerate(transverse_axes))))
            moving_index = transverse_axes.index(track_axis)
            reference = fibres[nearest_key]
            curve = {}
            for key, fraction in fibres.items():
                if any(abs(key[index] - nearest_key[index]) > 1e-6
                       for index in range(2) if index != moving_index):
                    continue
                offset = (key[moving_index] - nearest_key[moving_index]) / 10.0
                if abs(offset - round(offset)) > 1e-6:
                    continue
                curve[int(round(offset))] = fraction
            family_curves.append((fibre_axis, nearest_key, reference, curve))
            offsets = sorted(curve)
            plot_axis.plot(
                offsets, [curve[offset] / reference for offset in offsets],
                "o-", color=colours[family_index],
                label=f"{'xyz'[fibre_axis].upper()}-directed fibres")

        all_offsets = sorted(set(family_curves[0][3]) |
                             set(family_curves[1][3]))
        pair_reference = family_curves[0][2] + family_curves[1][2]
        pair_ratios = {}
        for offset in all_offsets:
            pair = (family_curves[0][3].get(offset, 0.0)
                    + family_curves[1][3].get(offset, 0.0))
            pair_ratios[offset] = pair / pair_reference
        plot_axis.plot(all_offsets, [pair_ratios[n] for n in all_offsets],
                       "s--", color="black", linewidth=2,
                       label="Combined pair")
        plot_axis.axhline(1.0, color="0.5", linestyle=":")
        plot_axis.set(title=f"Track along {track_name}", xlabel="Fibre offset n [10 mm pitch]",
                      ylabel="Response / nearest-fibre response", yscale="lin")
        plot_axis.grid(alpha=.2, which="both")
        plot_axis.legend(fontsize=8)

        for offset in all_offsets:
            record = {"track_axis": track_name, "offset_n": offset,
                      "offset_mm": 10.0 * offset}
            fibre_distances = []
            for index, (fibre_axis, key, reference, curve) in enumerate(family_curves, 1):
                transverse_axes = TRANSVERSE["xyz"[fibre_axis]]
                moving_index = transverse_axes.index(track_axis)
                target_key = list(key)
                target_key[moving_index] += 10.0 * offset
                distance = math.sqrt(sum(
                    (target_key[key_index] - source[axis]) ** 2
                    for key_index, axis in enumerate(transverse_axes)))
                fibre_distances.append(distance)
                record[f"fibre{index}_direction"] = "xyz"[fibre_axis]
                record[f"fibre{index}_nearest_index"] = str(key)
                record[f"fibre{index}_distance_mm"] = distance
                record[f"fibre{index}_reference_fraction"] = reference
                fraction = curve.get(offset, 0.0)
                record[f"fibre{index}_fraction"] = fraction
                record[f"fibre{index}_ratio"] = fraction / reference
            record["mean_pair_distance_mm"] = sum(fibre_distances) / 2.0
            record["combined_pair_ratio"] = pair_ratios.get(offset, "")
            rows.append(record)
            distance_plot_data.setdefault(track_name, []).append(
                (record["mean_pair_distance_mm"], record["combined_pair_ratio"], offset))

    figure.suptitle(f"{label}\ncentral source = (0.5, 30.5, 910.5) mm; bin {bin_number}")
    figure.savefig(output_dir / "central_fibre_correlation.png", dpi=170)
    plt.close(figure)

    distance_figure, distance_axes = plt.subplots(
        1, 3, figsize=(17, 5), constrained_layout=True)
    for axis, track_name in zip(distance_axes, "xyz"):
        points = sorted(distance_plot_data[track_name])
        axis.plot([point[0] for point in points], [point[1] for point in points],
                  "o-", color="#4c78a8")
        for distance, ratio, offset in points:
            axis.annotate(f"n={offset:+d}", (distance, ratio), xytext=(3, 4),
                          textcoords="offset points", fontsize=8)
        axis.set(title=f"Track along {track_name}",
                 xlabel="Mean distance to fibre pair [mm]",
                 ylabel="Combined pair response / reference", yscale="lin")
        axis.grid(alpha=.2, which="both")
    distance_figure.suptitle(
        f"{label}\ncentral source = (0.5, 30.5, 910.5) mm")
    distance_figure.savefig(
        output_dir / "central_fibre_correlation_vs_distance.png", dpi=170)
    plt.close(distance_figure)
    fields = list(rows[0])
    with (output_dir / "central_fibre_correlation.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Central source: global (0.5, 30.5, 910.5) mm, relative {source}, bin {bin_number}")
    print(f"Wrote {output_dir / 'central_fibre_correlation.csv'}")
    print(f"Wrote {output_dir / 'central_fibre_correlation.png'}")
    print(f"Wrote {output_dir / 'central_fibre_correlation_vs_distance.png'}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="Merged ROOT light map")
    parser.add_argument("--label", help="Setting shown in the plot")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    analyse(args.input,
            args.output_dir or Path("lightmap_analysis") / args.input.stem,
            args.label or args.input.stem)


if __name__ == "__main__":
    main()
