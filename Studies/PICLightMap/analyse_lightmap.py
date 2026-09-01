#!/usr/bin/env python3
"""Summarise light-map efficiency and draw 1D and 2D diagnostics."""

import argparse
import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/lightmap-matplotlib")
import matplotlib.pyplot as plt
import numpy as np
import ROOT


def read_efficiencies(filename: Path):
    if filename.suffix.lower() == ".csv":
        with filename.open(newline="") as stream:
            return [
                (float(r["x_mm"]), float(r["y_mm"]), float(r["z_mm"]),
                 float(r["fibre_hit_efficiency"]))
                for r in csv.DictReader(stream)
            ]

    root_file = ROOT.TFile.Open(str(filename), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")
    tree = root_file.Get("position_efficiency")
    if not tree:
        raise RuntimeError(
            f"{filename} has no position_efficiency tree; export it first with "
            "ExportLightMapEfficiency.C or pass its efficiency CSV")
    rows = [(float(e.x_mm), float(e.y_mm), float(e.z_mm),
             float(e.fibre_hit_efficiency)) for e in tree]
    root_file.Close()
    return rows


def read_view_efficiencies(filename: Path, rows=None):
    """Return total and x/y/z-directed-fibre efficiencies per source voxel."""
    if filename.suffix.lower() == ".csv":
        raise RuntimeError(
            "Per-view efficiencies require the merged ROOT light map; the "
            "efficiency CSV contains only the total")
    rows = rows if rows is not None else read_efficiencies(filename)
    root_file = ROOT.TFile.Open(str(filename), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")
    histogram = root_file.Get("Starting Vertex Distribution")
    centre = root_file.Get("cubeCentre")
    positions_dir = root_file.Get("hitPositions")
    fractions_dir = root_file.Get("lightFractions")
    if not all((histogram, centre, positions_dir, fractions_dir)):
        root_file.Close()
        raise RuntimeError(f"{filename} has incomplete per-fibre response data")
    centre_xyz = (float(centre.X()), float(centre.Y()), float(centre.Z()))
    result = []
    for x, y, z, total_efficiency in rows:
        relative = (x - centre_xyz[0], y - centre_xyz[1], z - centre_xyz[2])
        bin_number = histogram.FindBin(*relative)
        positions = positions_dir.Get(f"Bin_{bin_number}_hit_positions")
        fractions = fractions_dir.Get(f"Bin_{bin_number}_light_fractions")
        if not positions or fractions is None:
            root_file.Close()
            raise RuntimeError(f"Missing per-fibre response for bin {bin_number}")
        view = [0.0, 0.0, 0.0]
        for index in range(min(positions.GetEntries(), len(fractions))):
            position = positions.At(index)
            coordinates = (float(position.X()), float(position.Y()),
                           float(position.Z()))
            zero_axes = [axis for axis, value in enumerate(coordinates)
                         if abs(value) < 1e-6]
            if len(zero_axes) != 1:
                root_file.Close()
                raise RuntimeError(
                    f"Cannot identify fibre direction for {coordinates}")
            view[zero_axes[0]] += float(fractions[index])
        mean_view = sum(view) / 3.0
        view_rms = math.sqrt(sum((value - mean_view) ** 2 for value in view) / 3.0)
        result.append((x, y, z, total_efficiency, *view, view_rms))
    root_file.Close()
    return result


def analyse(rows, label: str, output_dir: Path, view_rows=None):
    if not rows:
        raise RuntimeError("The light map contains no valid scintillator positions")
    output_dir.mkdir(parents=True, exist_ok=True)
    values = np.asarray([r[3] for r in rows])
    mean = float(values.mean())
    std = float(values.std())
    standardized = (values - mean) / std
    skewness = float(np.mean(standardized ** 3))
    excess_kurtosis = float(np.mean(standardized ** 4) - 3.0)
    result = {
        "setting": label,
        "positions": len(rows),
        "mean": mean,
        "std": std,
        "coefficient_of_variation": std / mean,
        "skewness": skewness,
        "excess_kurtosis": excess_kurtosis,
        "min": float(values.min()),
        "max": float(values.max()),
        "relative_spread": (float(values.max()) - float(values.min())) / mean,
    }
    minimum = min(rows, key=lambda row: row[3])
    result.update(min_x_mm=minimum[0], min_y_mm=minimum[1],
                  min_z_mm=minimum[2])

    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=result.keys())
        writer.writeheader()
        writer.writerow(result)

    def draw_distribution(data, centre, width, xlabel, filename):
        figure, axis = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
        counts, edges, _ = axis.hist(
            data, bins="auto", color="#4c78a8", edgecolor="white",
            alpha=.85, label=f"{len(data)} positions")
        x_values = np.linspace(edges[0], edges[-1], 600)
        bin_width = float(np.mean(np.diff(edges)))
        gaussian = (len(data) * bin_width / (width * math.sqrt(2.0 * math.pi))
                    * np.exp(-0.5 * ((x_values - centre) / width) ** 2))
        axis.plot(x_values, gaussian, color="#e45756", linewidth=2.2,
                  label="Gaussian from mean and population std")
        axis.axvline(centre, color="black", linestyle="--", linewidth=1,
                     label="Mean")
        axis.set(xlabel=xlabel, ylabel="Number of positions", title=label)
        axis.grid(alpha=.2)
        axis.legend()
        axis.text(
            .98, .96,
            f"skewness = {skewness:.3f}\nexcess kurtosis = {excess_kurtosis:.3f}",
            transform=axis.transAxes, ha="right", va="top",
            bbox={"facecolor": "white", "alpha": .8, "edgecolor": "0.8"})
        figure.savefig(output_dir / filename, dpi=160)
        plt.close(figure)

    draw_distribution(values, mean, std, "Fibre-hit efficiency",
                      "efficiency_distribution.png")
    relative_deviation = (values - mean) / mean * 100.0
    draw_distribution(relative_deviation, 0.0, std / mean * 100.0,
                      "Relative deviation from mean [%]",
                      "relative_efficiency_distribution.png")

    xs = sorted({r[0] for r in rows})
    ys = sorted({r[1] for r in rows})
    zs = sorted({r[2] for r in rows})
    lookup = {(x, y, z): efficiency for x, y, z, efficiency in rows}

    # Infer the three fibre axes from the standard grid positions missing from
    # the 970-position map.  A fibre line is a group of ten missing positions
    # with two fixed coordinates and one coordinate spanning the grid.
    all_positions = {(x, y, z) for x in xs for y in ys for z in zs}
    missing_positions = all_positions - set(lookup)
    fibre_lines = []
    for axis in range(3):
        fixed_axes = [index for index in range(3) if index != axis]
        groups = {}
        for position in missing_positions:
            key = tuple(position[index] for index in fixed_axes)
            groups.setdefault(key, []).append(position)
        for fixed, positions in groups.items():
            if len(positions) == len((xs, ys, zs)[axis]):
                fibre_lines.append((axis, fixed))
    fibre_lines = list(dict.fromkeys(fibre_lines))
    if len(fibre_lines) != 3:
        raise RuntimeError(
            f"Expected three fibre lines from the 30 excluded positions, "
            f"inferred {len(fibre_lines)}")

    def distance_to_line(position, line):
        axis, fixed = line
        fixed_axes = [index for index in range(3) if index != axis]
        return math.sqrt(sum(
            (position[index] - value) ** 2
            for index, value in zip(fixed_axes, fixed)))

    distance_rows = []
    for x, y, z, efficiency in rows:
        distance = min(distance_to_line((x, y, z), line)
                       for line in fibre_lines)
        distance_rows.append((distance, efficiency))
    with (output_dir / "efficiency_vs_fibre_distance.csv").open(
            "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("distance_to_nearest_fibre_axis_mm", "efficiency"))
        writer.writerows(distance_rows)

    distances = np.asarray([row[0] for row in distance_rows])
    distance_efficiencies = np.asarray([row[1] for row in distance_rows])
    unique_distances = sorted(set(distances))
    distance_means = np.asarray([
        distance_efficiencies[distances == distance].mean()
        for distance in unique_distances])
    distance_stds = np.asarray([
        distance_efficiencies[distances == distance].std()
        for distance in unique_distances])
    distance_counts = np.asarray([
        np.count_nonzero(distances == distance) for distance in unique_distances])
    figure, axis = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    axis.scatter(distances, distance_efficiencies, s=13, alpha=.18,
                 color="#4c78a8", label="Individual positions")
    axis.errorbar(unique_distances, distance_means, yerr=distance_stds,
                  fmt="o-", color="#e45756", capsize=3,
                  label="Mean ± population std at each distance")
    for distance, average, count in zip(
            unique_distances, distance_means, distance_counts):
        axis.annotate(f"n={count}", (distance, average), xytext=(3, 6),
                      textcoords="offset points", fontsize=7)
    axis.set(xlabel="Distance to nearest fibre axis [mm]",
             ylabel="Fibre-hit efficiency",
             title=f"{label}: fibre-distance dependence")
    axis.grid(alpha=.2)
    axis.legend()
    figure.savefig(output_dir / "efficiency_vs_fibre_distance.png", dpi=160)
    plt.close(figure)

    vmin, vmax = float(values.min()), float(values.max())
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("lightgray")
    figure, axes = plt.subplots(2, 5, figsize=(18, 7), constrained_layout=True)
    overview_image = None
    for index, z in enumerate(zs):
        grid = np.full((len(ys), len(xs)), np.nan)
        for iy, y in enumerate(ys):
            for ix, x in enumerate(xs):
                grid[iy, ix] = lookup.get((x, y, z), np.nan)
        axis = axes.flat[index]
        overview_image = axis.imshow(
            grid, origin="lower", interpolation="nearest", cmap=cmap,
            vmin=vmin, vmax=vmax,
            extent=(xs[0] - .5, xs[-1] + .5, ys[0] - .5, ys[-1] + .5))
        axis.set(title=f"z = {z:g} mm", xlabel="x [mm]", ylabel="y [mm]")

        single, single_axis = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
        image = single_axis.imshow(
            grid, origin="lower", interpolation="nearest", cmap=cmap,
            vmin=vmin, vmax=vmax,
            extent=(xs[0] - .5, xs[-1] + .5, ys[0] - .5, ys[-1] + .5))
        single_axis.set(title=f"{label}: z = {z:g} mm",
                        xlabel="x [mm]", ylabel="y [mm]")
        single.colorbar(image, ax=single_axis, label="Fibre-hit efficiency")
        single.savefig(output_dir / f"efficiency_z_{z:g}mm.png", dpi=160)
        plt.close(single)
    figure.colorbar(overview_image, ax=axes, label="Fibre-hit efficiency", shrink=.9)
    figure.suptitle(f"{label}: ten efficiency slices")
    figure.savefig(output_dir / "efficiency_10_slices.png", dpi=160)
    plt.close(figure)

    if view_rows is not None:
        if len(view_rows) != len(rows):
            raise RuntimeError("Total and per-view efficiency rows do not match")
        fields = ("x_mm", "y_mm", "z_mm", "total_efficiency",
                  "x_fibre_efficiency", "y_fibre_efficiency",
                  "z_fibre_efficiency", "view_efficiency_rms")
        with (output_dir / "efficiency_per_view.csv").open(
                "w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(view_rows)

        per_view_lookup = {
            (row[0], row[1], row[2]): row[3:] for row in view_rows}
        summed_difference = np.asarray([
            sum(row[4:7]) - row[3] for row in view_rows])
        tolerance = max(1e-10, 1e-6 * float(np.max(np.abs(values))))
        if float(np.max(np.abs(summed_difference))) > tolerance:
            raise RuntimeError(
                "Per-view efficiencies do not reproduce total efficiency: "
                f"maximum difference {np.max(np.abs(summed_difference)):.6g}")

        quantities = (
            (0, "total_efficiency", "Total fibre-hit efficiency"),
            (1, "x_fibre_efficiency", "Efficiency: X-directed fibres"),
            (2, "y_fibre_efficiency", "Efficiency: Y-directed fibres"),
            (3, "z_fibre_efficiency", "Efficiency: Z-directed fibres"),
            (4, "view_efficiency_rms", "RMS among the three view efficiencies"),
        )
        for value_index, stem, colour_label in quantities:
            quantity_values = np.asarray([row[3 + value_index]
                                          for row in view_rows])
            quantity_min = float(quantity_values.min())
            quantity_max = float(quantity_values.max())
            slice_figure, slice_axes = plt.subplots(
                2, 5, figsize=(18, 7), constrained_layout=True)
            slice_image = None
            for slice_index, z in enumerate(zs):
                grid = np.full((len(ys), len(xs)), np.nan)
                for iy, y in enumerate(ys):
                    for ix, x in enumerate(xs):
                        record = per_view_lookup.get((x, y, z))
                        if record is not None:
                            grid[iy, ix] = record[value_index]
                slice_axis = slice_axes.flat[slice_index]
                slice_image = slice_axis.imshow(
                    grid, origin="lower", interpolation="nearest", cmap=cmap,
                    vmin=quantity_min, vmax=quantity_max,
                    extent=(xs[0] - .5, xs[-1] + .5,
                            ys[0] - .5, ys[-1] + .5))
                slice_axis.set(title=f"z = {z:g} mm", xlabel="x [mm]",
                               ylabel="y [mm]")
            slice_figure.colorbar(slice_image, ax=slice_axes,
                                  label=colour_label, shrink=.9)
            slice_figure.suptitle(f"{label}: {colour_label}, ten slices")
            slice_figure.savefig(output_dir / f"{stem}_10_slices.png", dpi=160)
            plt.close(slice_figure)

    print(f"{label}: {len(rows)} positions")
    print(f"  mean={mean:.10g}, std={std:.10g}, CV={std / mean:.6%}")
    print(f"  skewness={skewness:.6g}, excess kurtosis={excess_kurtosis:.6g}")
    print(f"  min={values.min():.10g} at ({minimum[0]:g}, {minimum[1]:g}, {minimum[2]:g}) mm")
    print(f"  max={values.max():.10g}, relative spread (max-min)/mean={result['relative_spread']:.6%}")
    print(f"  plots and summary: {output_dir}")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="Merged ROOT map or efficiency CSV")
    parser.add_argument("--label", help="Setting shown in plots")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    label = args.label or args.input.stem
    output_dir = args.output_dir or Path("lightmap_analysis") / args.input.stem
    rows = read_efficiencies(args.input)
    view_rows = (read_view_efficiencies(args.input, rows)
                 if args.input.suffix.lower() != ".csv" else None)
    analyse(rows, label, output_dir, view_rows)


if __name__ == "__main__":
    main()
