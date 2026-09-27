#!/usr/bin/env python3
"""Charge-aware side-by-side event displays for matched detector samples."""

import argparse
from collections import defaultdict
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
import ROOT


VIEWS = (("z", "x"), ("z", "y"), ("x", "y"))
CUBE_PITCH_MM = 10.0
# Global active-volume lower edges for the common HFGD/LFGD placement.
CUBE_GRID_ORIGIN = {"x": -1100.0, "y": -1120.0, "z": 310.0}


def draw_cube_grid(panel, axes, limits):
    """Draw projected cube boundaries aligned to the detector matrix."""
    for coordinate, limit, draw in (
            (axes[0], limits[0], panel.axvline),
            (axes[1], limits[1], panel.axhline)):
        origin = CUBE_GRID_ORIGIN[coordinate]
        first = origin + np.ceil((limit[0]-origin)/CUBE_PITCH_MM)*CUBE_PITCH_MM
        for value in np.arange(first, limit[1]+.5*CUBE_PITCH_MM,
                               CUBE_PITCH_MM):
            draw(value, color="0.78", lw=.45, alpha=.55, zorder=0)


def rows(tree, fields):
    result = []
    for row in tree:
        result.append({field: getattr(row, field) for field in fields})
    return result


def load(path):
    source = ROOT.TFile.Open(str(path))
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    data = {
        "fibres": rows(source["fiber_hits"],
                       ("event", "projection", "x", "y", "z", "charge")),
        "hits2d": rows(source["hits2d"],
                       ("event", "view", "x", "y", "z", "charge", "used")),
        "cubes": rows(source["hits3d"],
                      ("event", "x", "y", "z", "charge")),
        "tracks": rows(source["track_nodes"],
                       ("event", "track", "node", "x", "y", "z")),
        "truth": rows(source["mc_track_points"],
                      ("event", "track_id", "parent_id", "point",
                       "x", "y", "z")),
        "truth_segments": rows(source["mc_virtual_segments"],
                      ("event", "start_x", "start_y", "start_z",
                       "stop_x", "stop_y", "stop_z")),
    }
    source.Close()
    coordinate_ranges = {view: {axis: [np.inf, -np.inf]
                                for axis in "xyz"} for view in range(3)}
    for row in data["fibres"]:
        view = int(row["projection"])
        for axis in "xyz":
            coordinate_ranges[view][axis][0] = min(
                coordinate_ranges[view][axis][0], float(row[axis]))
            coordinate_ranges[view][axis][1] = max(
                coordinate_ranges[view][axis][1], float(row[axis]))
    data["fibre_axes"] = {
        view: tuple(sorted("xyz", key=lambda axis:
                    coordinate_ranges[view][axis][1]
                    - coordinate_ranges[view][axis][0], reverse=True)[:2])
        for view in range(3)
    }
    data["view_axes"] = dict(data["fibre_axes"])
    return data


def selected(data, category, event):
    return [row for row in data[category] if int(row["event"]) == event]


def truth_tracks(data, event, axes):
    grouped = defaultdict(list)
    for row in selected(data, "truth", event):
        if int(row["parent_id"]) != 0:
            continue
        grouped[int(row["track_id"])].append(
            (int(row["point"]), float(row[axes[0]]), float(row[axes[1]])))
    return [[(x, y) for _, x, y in sorted(points)]
            for points in grouped.values()]


def reco_tracks(data, event, axes):
    grouped = defaultdict(list)
    for row in selected(data, "tracks", event):
        grouped[int(row["track"])].append(
            (int(row["node"]), float(row[axes[0]]), float(row[axes[1]])))
    return [[(x, y) for _, x, y in sorted(points)]
            for points in grouped.values()]


def truth_segments(data, event, axes):
    return [[(float(row[f"start_{axes[0]}"]),
              float(row[f"start_{axes[1]}"])),
             (float(row[f"stop_{axes[0]}"]),
              float(row[f"stop_{axes[1]}"]))]
            for row in selected(data, "truth_segments", event)]


def truth_segments_3d(data, event):
    return [[(float(row["start_x"]), float(row["start_y"]),
              float(row["start_z"])),
             (float(row["stop_x"]), float(row["stop_y"]),
              float(row["stop_z"]))]
            for row in selected(data, "truth_segments", event)]


def charge_norm(values):
    positive = np.asarray([value for value in values if value > 0], dtype=float)
    if not len(positive):
        return LogNorm(1.0, 2.0)
    low = max(float(np.min(positive)), float(np.percentile(positive, 2)))
    high = max(low*1.01, float(np.percentile(positive, 99)))
    return LogNorm(low, high)


def weighted_quantile(values, weights, quantile):
    order = np.argsort(values)
    values = np.asarray(values, dtype=float)[order]
    weights = np.asarray(weights, dtype=float)[order]
    cumulative = np.cumsum(weights)
    if not len(values) or cumulative[-1] <= 0:
        return np.nan
    return float(np.interp(quantile*cumulative[-1], cumulative, values))


def lfgd_column_hits(data, event):
    """Make one charge-CoG 2D hit per dominant-axis slice and view."""
    fibres = selected(data, "fibres", event)
    spans = {}
    for axis in "xyz":
        values, weights = [], []
        for row in fibres:
            view = int(row["projection"])
            if axis not in data["fibre_axes"][view]:
                continue
            values.append(float(row[axis]))
            weights.append(max(float(row["charge"]), 0.0))
        spans[axis] = (weighted_quantile(values, weights, .95)
                       - weighted_quantile(values, weights, .05))
    dominant = max(spans, key=lambda axis: spans[axis])
    accumulators = defaultdict(lambda: [0.0, 0.0])
    for row in fibres:
        view = int(row["projection"])
        measured = data["fibre_axes"][view]
        if dominant not in measured:
            continue
        transverse = measured[0] if measured[1] == dominant else measured[1]
        charge = max(float(row["charge"]), 0.0)
        if charge <= 0:
            continue
        column = int(np.floor((float(row[dominant])
                               - CUBE_GRID_ORIGIN[dominant])/CUBE_PITCH_MM))
        accumulator = accumulators[(view, column, transverse)]
        accumulator[0] += charge
        accumulator[1] += charge*float(row[transverse])
    result = []
    for (view, column, transverse), (charge, weighted_position) in sorted(
            accumulators.items()):
        position = {"x": np.nan, "y": np.nan, "z": np.nan}
        position[dominant] = (CUBE_GRID_ORIGIN[dominant]
                              + (column+.5)*CUBE_PITCH_MM)
        position[transverse] = weighted_position/charge
        result.append({"event": event, "dominant_axis": dominant,
                       "view": view, "column_index": column,
                       "transverse_axis": transverse, "charge": charge,
                       **position})
    return result


def combine_column_hits(hits):
    """Combine the two perpendicular-view measurements into one 3D point."""
    grouped = defaultdict(list)
    for hit in hits:
        grouped[(hit["dominant_axis"], hit["column_index"])].append(hit)
    points = []
    for (dominant, column), measurements in sorted(grouped.items()):
        position = {"x": np.nan, "y": np.nan, "z": np.nan}
        position[dominant] = (CUBE_GRID_ORIGIN[dominant]
                              + (column+.5)*CUBE_PITCH_MM)
        charge = 0.0
        for measurement in measurements:
            transverse = measurement["transverse_axis"]
            position[transverse] = float(measurement[transverse])
            charge += float(measurement["charge"])
        if all(np.isfinite(position[axis]) for axis in "xyz"):
            points.append({"event": int(measurements[0]["event"]),
                           "dominant_axis": dominant,
                           "column_index": column, "charge": charge/2.0,
                           **position})
    return points


def fit_column_track(points):
    """Charge-weighted orthogonal straight-line fit using 3D PCA."""
    if len(points) < 2:
        return None
    coordinates = np.asarray([[row[axis] for axis in "xyz"]
                              for row in points], dtype=float)
    weights = np.asarray([max(float(row["charge"]), 1e-9)
                          for row in points], dtype=float)
    centre = np.average(coordinates, axis=0, weights=weights)
    centered = coordinates-centre
    covariance = (centered*weights[:, None]).T @ centered/weights.sum()
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    direction = eigenvectors[:, np.argmax(eigenvalues)]
    dominant_index = "xyz".index(points[0]["dominant_axis"])
    if direction[dominant_index] < 0:
        direction *= -1
    projection = centered @ direction
    residual = centered-np.outer(projection, direction)
    rms = np.sqrt(np.average(np.sum(residual*residual, axis=1),
                             weights=weights))
    return {"point": centre, "direction": direction, "rms": rms,
            "minimum": float(projection.min()),
            "maximum": float(projection.max())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hfgd")
    parser.add_argument("lfgd")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--event-range", nargs=2, type=int, default=(0, 9))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    samples = (("HFGD standard reconstruction", load(args.hfgd)),
               ("LFGD standard reconstruction", load(args.lfgd)))
    lfgd_columns = {event: lfgd_column_hits(samples[1][1], event)
                    for event in range(args.event_range[0],
                                       args.event_range[1]+1)}
    lfgd_points = {event: combine_column_hits(hits)
                   for event, hits in lfgd_columns.items()}
    lfgd_fits = {event: fit_column_track(points)
                 for event, points in lfgd_points.items()}
    with (args.output_dir / "lfgd_column_2d_hits.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        fields = ("event", "dominant_axis", "view", "column_index",
                  "transverse_axis", "charge", "x", "y", "z")
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for event in sorted(lfgd_columns):
            writer.writerows(lfgd_columns[event])
    with (args.output_dir / "lfgd_column_3d_points.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        fields = ("event", "dominant_axis", "column_index", "charge",
                  "x", "y", "z")
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for event in sorted(lfgd_points):
            writer.writerows(lfgd_points[event])
    with (args.output_dir / "lfgd_column_track_fits.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("event", "n_points", "point_x", "point_y", "point_z",
                         "direction_x", "direction_y", "direction_z",
                         "transverse_rms_mm"))
        for event, fit in sorted(lfgd_fits.items()):
            if fit is not None:
                writer.writerow((event, len(lfgd_points[event]), *fit["point"],
                                 *fit["direction"], fit["rms"]))

    for event in range(args.event_range[0], args.event_range[1] + 1):
        for axes in VIEWS:
            truth = [point for _, data in samples
                     for segment in truth_segments(data, event, axes)
                     for point in segment]
            if not truth:
                truth = [point for _, data in samples
                         for track in truth_tracks(data, event, axes)
                         for point in track]
            if not truth:
                continue
            xs, ys = zip(*truth)
            xmargin = max(20.0, 0.08*(max(xs)-min(xs)))
            ymargin = max(20.0, 0.08*(max(ys)-min(ys)))
            limits = ((min(xs)-xmargin, max(xs)+xmargin),
                      (min(ys)-ymargin, max(ys)+ymargin))

            hit2d_charge = []
            for _, data in samples:
                hit2d_charge.extend(float(row["charge"])
                    for row in selected(data, "hits2d", event)
                    if set(data["view_axes"][int(row["view"])])
                       == set(axes))
            hit2d_norm = charge_norm(hit2d_charge)

            fig, panels = plt.subplots(1, 2, figsize=(15, 6.3), sharex=True,
                                       sharey=True, constrained_layout=True)
            hit2d_artist = None
            for panel, (label, data) in zip(panels, samples):
                hits2d = [row for row in selected(data, "hits2d", event)
                          if set(data["view_axes"][int(row["view"])])
                             == set(axes)]
                used = [row for row in hits2d if bool(row["used"])]
                unused = [row for row in hits2d if not bool(row["used"])]
                if used:
                    hit2d_artist = panel.scatter(
                        [float(row[axes[0]]) for row in used],
                        [float(row[axes[1]]) for row in used],
                        c=[float(row["charge"]) for row in used],
                        cmap="viridis", norm=hit2d_norm, marker="o", s=30,
                        edgecolors="black", linewidths=.25,
                        label="used 2D hit")
                if unused:
                    panel.scatter(
                        [float(row[axes[0]]) for row in unused],
                        [float(row[axes[1]]) for row in unused],
                        c=[float(row["charge"]) for row in unused],
                        cmap="viridis", norm=hit2d_norm, marker="x", s=42,
                        linewidths=1.0, label="unused 2D hit")
                if label.startswith("LFGD"):
                    column_points = lfgd_points[event]
                    if column_points:
                        panel.scatter(
                            [float(row[axes[0]]) for row in column_points],
                            [float(row[axes[1]]) for row in column_points],
                            facecolors="none", edgecolors="deepskyblue",
                            marker="s", s=34, linewidths=1.1,
                            label="one CoG point / virtual-cube column")
                    fit = lfgd_fits[event]
                    if fit is not None:
                        parameters = np.linspace(fit["minimum"], fit["maximum"], 2)
                        line = fit["point"][None, :] + parameters[:, None]*fit["direction"][None, :]
                        indices = ("xyz".index(axes[0]), "xyz".index(axes[1]))
                        panel.plot(line[:, indices[0]], line[:, indices[1]],
                                   "--", color="deepskyblue", lw=1.8,
                                   label="column-point line fit")
                segments = truth_segments(data, event, axes)
                if segments:
                    collection = LineCollection(
                        segments, colors="limegreen", linewidths=1.8,
                        label="MC truth hit segments", zorder=4)
                    panel.add_collection(collection)
                else:
                    for index, track in enumerate(truth_tracks(data, event, axes)):
                        tx, ty = zip(*track)
                        panel.plot(tx, ty, "-", color="limegreen", lw=2.0,
                                   marker=".", ms=3,
                                   label="primary MC truth" if index == 0 else None)
                panel.set(title=label, xlabel=f"{axes[0]} [mm]",
                          xlim=limits[0], ylim=limits[1])
                draw_cube_grid(panel, axes, limits)
                panel.grid(alpha=.18, linewidth=.8)
                panel.legend(loc="best", fontsize=8)
            panels[0].set_ylabel(f"{axes[1]} [mm]")
            fig.suptitle(f"Matched event {event}: reconstructed 2D hits")
            if hit2d_artist is not None:
                fig.colorbar(hit2d_artist, ax=panels, shrink=.82,
                             label="2D-hit charge")
            fig.savefig(args.output_dir / f"event{event}_{axes[0]}{axes[1]}.png",
                        dpi=150)
            plt.close(fig)

        segments3d = [segment for _, data in samples
                      for segment in truth_segments_3d(data, event)]
        if not segments3d:
            continue
        truth_points3d = [point for segment in segments3d for point in segment]
        limits3d = []
        for coordinate in range(3):
            values = [point[coordinate] for point in truth_points3d]
            margin = max(20.0, .08*(max(values)-min(values)))
            limits3d.append((min(values)-margin, max(values)+margin))
        cube_charge = [float(row["charge"]) for _, data in samples
                       for row in selected(data, "cubes", event)]
        cube_norm = charge_norm(cube_charge)
        fig = plt.figure(figsize=(15, 7), constrained_layout=True)
        panels = [fig.add_subplot(1, 2, index+1, projection="3d")
                  for index in range(2)]
        cube_artist = None
        for panel, (label, data) in zip(panels, samples):
            cubes = selected(data, "cubes", event)
            if cubes:
                cube_artist = panel.scatter(
                    [float(row["x"]) for row in cubes],
                    [float(row["y"]) for row in cubes],
                    [float(row["z"]) for row in cubes],
                    c=[float(row["charge"]) for row in cubes],
                    cmap="autumn_r", norm=cube_norm, s=20, alpha=.72,
                    edgecolors="black", linewidths=.2,
                    label="reconstructed 3D hit")
            mc_segments = truth_segments_3d(data, event)
            if mc_segments:
                panel.add_collection3d(Line3DCollection(
                    mc_segments, colors="limegreen", linewidths=1.7,
                    label="MC truth hit segments"))
            grouped = defaultdict(list)
            for row in selected(data, "tracks", event):
                grouped[int(row["track"])].append(
                    (int(row["node"]), float(row["x"]), float(row["y"]),
                     float(row["z"])))
            for index, nodes in enumerate(grouped.values()):
                ordered = sorted(nodes)
                panel.plot([row[1] for row in ordered],
                           [row[2] for row in ordered],
                           [row[3] for row in ordered], color="red", lw=.8,
                           label="reconstructed tracks" if index == 0 else None)
            if label.startswith("LFGD") and lfgd_points[event]:
                points = lfgd_points[event]
                panel.scatter([row["x"] for row in points],
                              [row["y"] for row in points],
                              [row["z"] for row in points],
                              facecolors="none", edgecolors="deepskyblue",
                              marker="s", s=28, linewidths=1.0,
                              label="column CoG points")
                fit = lfgd_fits[event]
                if fit is not None:
                    parameters = np.linspace(fit["minimum"], fit["maximum"], 2)
                    line = fit["point"][None, :] + parameters[:, None]*fit["direction"][None, :]
                    panel.plot(line[:, 0], line[:, 1], line[:, 2], "--",
                               color="deepskyblue", lw=2.0,
                               label="column-point line fit")
            panel.set(title=label, xlabel="x [mm]", ylabel="y [mm]",
                      zlabel="z [mm]", xlim=limits3d[0], ylim=limits3d[1],
                      zlim=limits3d[2])
            panel.legend(loc="best", fontsize=8)
        fig.suptitle(f"Matched event {event}: reconstructed 3D hits")
        if cube_artist is not None:
            fig.colorbar(cube_artist, ax=panels, shrink=.72,
                         label="3D-hit charge")
        fig.savefig(args.output_dir / f"event{event}_3d.png", dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()
