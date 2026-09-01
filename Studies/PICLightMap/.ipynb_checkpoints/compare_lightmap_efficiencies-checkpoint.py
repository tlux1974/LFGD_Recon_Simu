#!/usr/bin/env python3
"""Compare total, per-view, and nearest-fibre light-map efficiencies."""

import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/lightmap-matplotlib")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

from analyse_lightmap import read_efficiencies
from analyse_fibre_correlation import TRANSVERSE, fibre_direction, load_response


AXIS_NAMES = ("x", "y", "z")
CONTAINMENT_FRACTIONS = (0.5, 0.7, 0.8, 0.9)


def _population_rms(values):
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean((values - values.mean()) ** 2)))


def read_detailed_efficiencies(filename):
    """Read total/view sums and the geometrically nearest fibre in each view."""
    filename = Path(filename)
    totals = read_efficiencies(filename)
    root_file = ROOT.TFile.Open(str(filename), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")
    starts = root_file.Get("Starting Vertex Distribution")
    centre = root_file.Get("cubeCentre")
    positions_dir = root_file.Get("hitPositions")
    fractions_dir = root_file.Get("lightFractions")
    if not all((starts, centre, positions_dir, fractions_dir)):
        root_file.Close()
        raise RuntimeError(f"Incomplete light map: {filename}")
    centre_xyz = np.asarray((centre.X(), centre.Y(), centre.Z()), dtype=float)
    records = []
    for x, y, z, total in totals:
        source = np.asarray((x, y, z), dtype=float)
        relative = source - centre_xyz
        bin_number = starts.FindBin(*relative)
        positions = positions_dir.Get(f"Bin_{bin_number}_hit_positions")
        fractions = fractions_dir.Get(f"Bin_{bin_number}_light_fractions")
        if not positions or fractions is None:
            root_file.Close()
            raise RuntimeError(f"Missing response for bin {bin_number}")
        view_sum = np.zeros(3)
        nearest_distance = np.full(3, np.inf)
        nearest_fraction = np.zeros(3)
        for index in range(min(positions.GetEntries(), len(fractions))):
            position = positions.At(index)
            fibre = np.asarray((position.X(), position.Y(), position.Z()),
                               dtype=float)
            zero_axes = np.flatnonzero(np.abs(fibre) < 1e-6)
            if len(zero_axes) != 1:
                root_file.Close()
                raise RuntimeError(f"Cannot identify fibre direction: {fibre}")
            axis = int(zero_axes[0])
            fraction = float(fractions[index])
            view_sum[axis] += fraction
            transverse = [coordinate for coordinate in range(3)
                          if coordinate != axis]
            # Stored fibre coordinates and the light-map histogram are both
            # relative to cubeCentre; x/y/z in position_efficiency are global.
            distance = float(np.linalg.norm(
                fibre[transverse] - relative[transverse]))
            if distance < nearest_distance[axis]:
                nearest_distance[axis] = distance
                nearest_fraction[axis] = fraction
        view_mean = float(view_sum.mean())
        nearest_mean = float(nearest_fraction.mean())
        records.append({
            "x_mm": x, "y_mm": y, "z_mm": z, "total": total,
            **{f"view_{name}": float(view_sum[axis])
               for axis, name in enumerate(AXIS_NAMES)},
            "view_rms": float(np.sqrt(np.mean((view_sum-view_mean)**2))),
            **{f"nearest_{name}": float(nearest_fraction[axis])
               for axis, name in enumerate(AXIS_NAMES)},
            **{f"nearest_distance_{name}_mm": float(nearest_distance[axis])
               for axis, name in enumerate(AXIS_NAMES)},
            "nearest_view_rms": float(np.sqrt(np.mean(
                (nearest_fraction-nearest_mean)**2))),
        })
    root_file.Close()
    return records


def _statistics(values):
    values = np.asarray(values, dtype=float)
    return {"mean": float(values.mean()), "min": float(values.min()),
            "max": float(values.max()), "rms": _population_rms(values)}


def track_containment_radii(filename, fractions=CONTAINMENT_FRACTIONS):
    """Containment radii for central ten-subvoxel tracks along x, y, and z."""
    root_file = ROOT.TFile.Open(str(filename), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")
    centre = root_file.Get("cubeCentre")
    global_centre = np.asarray((centre.X(), centre.Y(), centre.Z()), dtype=float)
    result = {}
    for track_axis, track_name in enumerate(AXIS_NAMES):
        integrated = {}
        for coordinate in (value - 4.5 for value in range(10)):
            relative = np.asarray((0.5, 0.5, 0.5), dtype=float)
            relative[track_axis] = coordinate
            source = tuple(global_centre + relative)
            _, _, response = load_response(root_file, source)
            for fibre_position, response_fraction in response:
                direction = fibre_direction(fibre_position)
                if direction is None:
                    continue
                transverse_axes = TRANSVERSE[AXIS_NAMES[direction]]
                key = (direction,) + tuple(fibre_position[axis]
                                           for axis in transverse_axes)
                integrated[key] = integrated.get(key, 0.0) + response_fraction

        response_by_radius = {}
        midpoint = (0.5, 0.5, 0.5)
        for key, response_sum in integrated.items():
            fibre_axis = key[0]
            transverse_axes = TRANSVERSE[AXIS_NAMES[fibre_axis]]
            fixed = dict(zip(transverse_axes, key[1:]))
            if fibre_axis == track_axis:
                radius = math.sqrt(sum((fixed[axis]-midpoint[axis])**2
                                       for axis in transverse_axes))
            else:
                common_axis = next(axis for axis in range(3)
                                   if axis not in (fibre_axis, track_axis))
                radius = abs(fixed[common_axis] - midpoint[common_axis])
            response_by_radius[radius] = (
                response_by_radius.get(radius, 0.0) + response_sum)

        total_response = sum(response_by_radius.values())
        cumulative = 0.0
        radial = []
        for radius in sorted(response_by_radius):
            cumulative += response_by_radius[radius]
            radial.append((radius, cumulative))
        axis_result = {"track_total_efficiency": total_response / 10.0}
        for fraction in fractions:
            suffix = str(int(round(100*fraction)))
            relative_target = fraction * total_response
            absolute_target = fraction * 10.0
            axis_result[f"relative_r{suffix}_mm"] = next(
                radius for radius, value in radial if value >= relative_target)
            axis_result[f"absolute_r{suffix}_mm"] = next(
                (radius for radius, value in radial if value >= absolute_target),
                math.nan)
        result[track_name] = axis_result
    root_file.Close()
    return result


def compare_lightmaps(settings, output_dir):
    """Analyse ``[(scatter_length_mm, label, ROOT_path), ...]``."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    detailed = {}
    quantities = (["total"] + [f"view_{name}" for name in AXIS_NAMES]
                  + ["view_rms"]
                  + [f"nearest_{name}" for name in AXIS_NAMES]
                  + ["nearest_view_rms"])
    for scatter_length, label, filename in sorted(settings):
        records = read_detailed_efficiencies(filename)
        containment = track_containment_radii(filename)
        detailed[label] = records
        row = {"scatter_length_mm": scatter_length, "label": label,
               "input": str(filename), "positions": len(records)}
        for quantity in quantities:
            stats = _statistics([record[quantity] for record in records])
            for statistic, value in stats.items():
                row[f"{quantity}_{statistic}"] = value
        for normalization in ("relative", "absolute"):
            for fraction in CONTAINMENT_FRACTIONS:
                suffix = str(int(round(100*fraction)))
                axis_values = [containment[axis][
                    f"{normalization}_r{suffix}_mm"] for axis in AXIS_NAMES]
                for axis, value in zip(AXIS_NAMES, axis_values):
                    row[f"{normalization}_r{suffix}_{axis}_mm"] = value
                row[f"{normalization}_r{suffix}_mean_mm"] = (
                    float(np.mean(axis_values))
                    if np.all(np.isfinite(axis_values)) else math.nan)
        summaries.append(row)

    with (output_dir / "lightmap_efficiency_comparison.csv").open(
            "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    lengths = np.asarray([row["scatter_length_mm"] for row in summaries])
    colours = {"x": "#e45756", "y": "#4c78a8", "z": "#54a24b"}

    figure, axis = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    axis.plot(lengths, [row["total_mean"] for row in summaries], "o-",
              color="black", label="Mean")
    axis.plot(lengths, [row["total_min"] for row in summaries], "v--",
              color="#e45756", label="Minimum")
    axis.plot(lengths, [row["total_max"] for row in summaries], "^--",
              color="#54a24b", label="Maximum")
    axis.set(xlabel="Scattering length [mm]", ylabel="Total light fraction",
             title="Total light collection versus scattering length")
    axis.grid(alpha=.25); axis.legend()
    figure.savefig(output_dir / "total_efficiency_vs_scattering_length.png",
                   dpi=170)
    plt.close(figure)

    def plot_view_statistics(prefix, ylabel, filename, title):
        figure, axes = plt.subplots(2, 2, figsize=(13, 9),
                                    constrained_layout=True)
        for axis, statistic in zip(axes.flat, ("mean", "min", "max", "rms")):
            for name in AXIS_NAMES:
                axis.plot(lengths,
                          [row[f"{prefix}_{name}_{statistic}"]
                           for row in summaries], "o-", color=colours[name],
                          label=f"{name.upper()}-directed fibres")
            axis.set(title=statistic.upper(), xlabel="Scattering length [mm]",
                     ylabel=ylabel)
            axis.grid(alpha=.25); axis.legend(fontsize=8)
        figure.suptitle(title)
        figure.savefig(output_dir / filename, dpi=170)
        plt.close(figure)

    plot_view_statistics(
        "view", "Summed light fraction in view",
        "view_efficiency_statistics_vs_scattering_length.png",
        "Per-view summed efficiency across source subvoxels")
    plot_view_statistics(
        "nearest", "Light fraction reaching nearest fibre",
        "nearest_fibre_statistics_vs_scattering_length.png",
        "Nearest-fibre response across source subvoxels")

    figure, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    for axis, prefix, title in (
            (axes[0], "view_rms", "RMS among summed view efficiencies"),
            (axes[1], "nearest_view_rms", "RMS among nearest-fibre fractions")):
        for statistic, marker, curve_label in (
                ("mean", "o", "MEAN"), ("min", "v", "MIN"),
                ("max", "^", "MAX"),
                ("rms", "s", "STD OF VIEW RMS")):
            axis.plot(lengths, [row[f"{prefix}_{statistic}"]
                                for row in summaries], marker+"-",
                      label=curve_label)
        axis.set(title=title, xlabel="Scattering length [mm]",
                 ylabel="Light-fraction RMS")
        axis.grid(alpha=.25); axis.legend(fontsize=8)
    figure.savefig(output_dir / "view_imbalance_vs_scattering_length.png",
                   dpi=170)
    plt.close(figure)

    containment_colours = ("#4c78a8", "#f58518", "#54a24b", "#e45756")
    for normalization, qualifier in (
            ("relative", "fraction of collected light"),
            ("absolute", "fraction of all launched photons")):
        figure, axis = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
        for fraction, colour in zip(CONTAINMENT_FRACTIONS,
                                    containment_colours):
            suffix = str(int(round(100*fraction)))
            axis.plot(lengths,
                      [row[f"{normalization}_r{suffix}_mean_mm"]
                       for row in summaries], "o-", color=colour,
                      label=f"{fraction:.1f}")
        axis.set(xlabel="Scattering length [mm]",
                 ylabel="Mean containment radius over X/Y/Z tracks [mm]",
                 title=f"Track-light containment: {qualifier}")
        axis.grid(alpha=.25); axis.legend(title="Target fraction")
        figure.savefig(output_dir /
                       f"track_containment_{normalization}_vs_scattering_length.png",
                       dpi=170)
        plt.close(figure)
    return summaries, detailed
