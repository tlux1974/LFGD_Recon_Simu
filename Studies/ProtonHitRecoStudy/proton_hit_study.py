"""Analysis helpers for the isolated 1 GeV proton hit-position study."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot


VIEW_NAMES = {0: "XZ (fibres along y)", 1: "YZ (fibres along x)",
              2: "XY (fibres along z)"}
PROJECTIONS = (("x", "y"), ("x", "z"), ("y", "z"))


def _tree(path, name, columns=None):
    path = Path(path)
    with uproot.open(path) as source:
        if name not in source:
            return pd.DataFrame(columns=columns or [])
        tree = source[name]
        available = set(tree.keys())
        selected = list(available if columns is None else
                        [column for column in columns if column in available])
        arrays = tree.arrays(selected, library="np")
    return pd.DataFrame({key: np.asarray(value) for key, value in arrays.items()})


def load_run(path):
    path = Path(path)
    common = ["event", "x", "y", "z", "charge", "geom_id"]
    return {
        "path": path,
        "hits": _tree(path, "hits3d", common),
        "hits2d": _tree(path, "hits2d",
                        ["event", "view", "x", "y", "z", "charge",
                         "geom_id", "used", "rejection_code",
                         "rejection_reason"]),
        "views": _tree(path, "hit3d_views",
                       ["event", "hit", "view", "view_charge", "view_x",
                        "view_y", "view_z", "fiber_count",
                        "fiber_charge_sum", "geom_id", "view_geom_id"]),
        "tracks": _tree(path, "track_nodes",
                        ["event", "track", "node", "x", "y", "z"]),
        "truth_tracks": _tree(path, "mc_track_points",
                              ["event", "track_id", "parent_id", "pdg",
                               "point", "x", "y", "z"]),
        "truth_segments": _tree(path, "mc_virtual_segments",
                                ["event", "geom_id", "primary_pdg", "start_x",
                                 "start_y", "start_z", "stop_x", "stop_y",
                                 "stop_z", "energy_deposit", "track_length"]),
    }


def load_comparison(default_path, alternative_path):
    return {"Standard": load_run(default_path),
            "Dominant-axis": load_run(alternative_path)}


def summary_table(runs):
    rows = []
    for label, run in runs.items():
        events = set(run["hits"].event.astype(int)) if not run["hits"].empty else set()
        hits2d = run["hits2d"]
        rows.append({"method": label, "events_with_3d_hits": len(events),
                     "3d_hits": len(run["hits"]),
                     "2d_hits": len(hits2d) if not hits2d.empty else len(run["views"]),
                     "used_2d_hits": (int(hits2d.used.sum())
                                      if not hits2d.empty else len(run["views"])),
                     "unused_2d_hits": (int((~hits2d.used.astype(bool)).sum())
                                        if not hits2d.empty else np.nan),
                     "tracks": (len(run["tracks"][["event", "track"]].drop_duplicates())
                                if not run["tracks"].empty else 0),
                     "track_nodes": len(run["tracks"])})
    return pd.DataFrame(rows)


def reconstruction_counts(run):
    """Return event-level counts at the saved reconstruction stages.

    The 2D count is the number of distinct reconstructed view hits attached to
    at least one 3D hit. A view hit reused by several candidate 3D hits is
    counted once. Rejected/unmatched 2D candidates are not in the flat tree.
    """
    event_sources = [frame.event.astype(int) for frame in
                     (run["hits"], run["hits2d"], run["views"], run["tracks"])
                     if not frame.empty and "event" in frame]
    columns = ["event", "used_2d_hits", "unused_2d_hits", "view_0_hits",
               "view_1_hits", "view_2_hits", "3d_hits", "tracks"]
    if not event_sources:
        return pd.DataFrame(columns=columns)
    result = pd.DataFrame({"event": np.unique(np.concatenate(event_sources))})

    if not run["hits2d"].empty:
        unique_views = run["hits2d"]
        used_views = unique_views[unique_views.used.astype(bool)]
    else:
        views = run["views"].copy()
        identity = [column for column in
                    ("event", "view", "view_geom_id", "view_x", "view_y",
                     "view_z", "view_charge") if column in views]
        unique_views = views.drop_duplicates(identity) if identity else views
        used_views = unique_views
    total_2d = used_views.groupby("event").size()
    result["used_2d_hits"] = result.event.map(total_2d).fillna(0).astype(int)
    unused_2d = (unique_views[~unique_views.used.astype(bool)].groupby("event").size()
                 if "used" in unique_views else pd.Series(dtype=int))
    result["unused_2d_hits"] = result.event.map(unused_2d).fillna(0).astype(int)
    for view in range(3):
        counts = used_views[used_views.view == view].groupby("event").size()
        result[f"view_{view}_hits"] = result.event.map(counts).fillna(0).astype(int)

    hits3d = (run["hits"].groupby("event").size()
              if not run["hits"].empty else pd.Series(dtype=int))
    result["3d_hits"] = result.event.map(hits3d).fillna(0).astype(int)
    tracks = run["tracks"]
    track_counts = (tracks[["event", "track"]].drop_duplicates()
                    .groupby("event").size()) if not tracks.empty else pd.Series(dtype=int)
    result["tracks"] = result.event.map(track_counts).fillna(0).astype(int)
    return result[columns]


def event_count_table(runs, event):
    rows = []
    for label, run in runs.items():
        selected = reconstruction_counts(run)
        selected = selected[selected.event == event]
        row = selected.iloc[0].to_dict() if len(selected) else {"event": event}
        row["method"] = label
        rows.append(row)
    columns = ["method", "event", "view_0_hits", "view_1_hits",
               "view_2_hits", "used_2d_hits", "unused_2d_hits", "3d_hits",
               "tracks"]
    return pd.DataFrame(rows).reindex(columns=columns).fillna(0)


def count_summary_table(runs):
    rows = []
    for label, run in runs.items():
        counts = reconstruction_counts(run)
        for quantity in ("view_0_hits", "view_1_hits", "view_2_hits",
                         "used_2d_hits", "unused_2d_hits", "3d_hits", "tracks"):
            values = counts[quantity]
            rows.append({"method": label, "quantity": quantity,
                         "events": len(values), "total": int(values.sum()),
                         "mean_per_event": values.mean(),
                         "median_per_event": values.median(),
                         "minimum": int(values.min()) if len(values) else 0,
                         "maximum": int(values.max()) if len(values) else 0})
    return pd.DataFrame(rows)


def plot_count_distributions(runs, bins=50):
    quantities = (("used_2d_hits", "Used 2D hits"),
                  ("unused_2d_hits", "Unused 2D hits"),
                  ("3d_hits", "3D hits"), ("tracks", "Tracks"))
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.5), constrained_layout=True)
    colours = {"Standard": "tab:blue", "Dominant-axis": "tab:orange"}
    for ax, (quantity, title) in zip(axes, quantities):
        for label, run in runs.items():
            values = reconstruction_counts(run)[quantity]
            ax.hist(values, bins=bins, histtype="step", lw=2,
                    color=colours.get(label), label=label)
        ax.set(xlabel=f"{title} per event", ylabel="Events", title=title)
        ax.legend()
    return fig


def rejection_summary_table(runs, event=None):
    rows = []
    for label, run in runs.items():
        frame = run["hits2d"]
        if frame.empty:
            continue
        if event is not None:
            frame = frame[frame.event == event]
        unused = frame[~frame.used.astype(bool)]
        for (view, reason), group in unused.groupby(["view", "rejection_reason"]):
            rows.append({"method": label, "event": event if event is not None else "all",
                         "view": int(view), "view_name": VIEW_NAMES.get(int(view), str(view)),
                         "rejection_reason": reason, "hits": len(group)})
    return pd.DataFrame(rows)


def _pca_line(points):
    points = np.asarray(points, dtype=float)
    if len(points) < 2:
        return None
    centre = points.mean(axis=0)
    _, _, vt = np.linalg.svd(points-centre, full_matrices=False)
    return centre, vt[0]


def _line_endpoints(points, line):
    centre, direction = line
    scale = (points-centre) @ direction
    return centre+scale.min()*direction, centre+scale.max()*direction


def plot_event(runs, event=0, charge_scale=(1.0, 150.0)):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    truth = next(iter(runs.values()))["truth_tracks"]
    truth = truth[(truth.event == event) & (truth.pdg == 2212) &
                  (truth.parent_id == 0)].sort_values(["track_id", "point"])
    colours = {"Standard": "tab:blue", "Dominant-axis": "tab:orange"}
    markers = {"Standard": "o", "Dominant-axis": "x"}
    for ax, (horizontal, vertical) in zip(axes, PROJECTIONS):
        if not truth.empty:
            for _, trajectory in truth.groupby("track_id"):
                ax.plot(trajectory[horizontal], trajectory[vertical],
                        color="green", lw=2, label="MC proton")
        for label, run in runs.items():
            hits = run["hits"][run["hits"].event == event]
            if hits.empty:
                continue
            sizes = 18+50*np.sqrt(np.maximum(hits.charge.to_numpy(), 0)/
                                  max(float(hits.charge.max()), 1.0))
            ax.scatter(hits[horizontal], hits[vertical], s=sizes,
                       marker=markers[label], facecolors="none" if label == "Standard" else None,
                       color=colours[label], alpha=.8, label=f"{label} 3D hits")
            line = _pca_line(hits[["x", "y", "z"]].to_numpy())
            if line is not None:
                first, last = _line_endpoints(hits[["x", "y", "z"]].to_numpy(), line)
                index = {"x": 0, "y": 1, "z": 2}
                ax.plot([first[index[horizontal]], last[index[horizontal]]],
                        [first[index[vertical]], last[index[vertical]]],
                        ls="--", color=colours[label], lw=1.5,
                        label=f"{label} PCA")
            nodes = run["tracks"][run["tracks"].event == event]
            for track_index, nodes_one in nodes.groupby("track"):
                nodes_one = nodes_one.sort_values("node")
                ax.plot(nodes_one[horizontal], nodes_one[vertical], color=colours[label],
                        lw=2.5, alpha=.55,
                        label=f"{label} standard fitter" if track_index == nodes.track.min() else None)
        ax.set(xlabel=f"{horizontal} [mm]", ylabel=f"{vertical} [mm]",
               title=f"Event {event}: {horizontal.upper()}{vertical.upper()}")
        ax.set_aspect("equal", adjustable="datalim")
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), fontsize=8)
    return fig


def plot_2d_views(runs, event=0):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    colours = {"Standard": "tab:blue", "Dominant-axis": "tab:orange"}
    for view in range(3):
        for label, run in runs.items():
            frame = run["views"]
            needed = {"view_x", "view_y", "view_z"}
            if frame.empty or not needed.issubset(frame.columns):
                continue
            frame = frame[(frame.event == event) & (frame.view == view)]
            if frame.empty:
                continue
            if view == 0: horizontal, vertical = "view_x", "view_z"
            elif view == 1: horizontal, vertical = "view_y", "view_z"
            else: horizontal, vertical = "view_x", "view_y"
            axes[0, view].scatter(frame[horizontal], frame[vertical], s=28,
                                  alpha=.65, color=colours[label], label=label)
            axes[1, view].hist(frame.view_charge, bins=40, histtype="step",
                               lw=2, color=colours[label], label=label)
            diagnostics = run["hits2d"]
            if not diagnostics.empty:
                diagnostics = diagnostics[(diagnostics.event == event)
                                          & (diagnostics.view == view)
                                          & (~diagnostics.used.astype(bool))]
                if not diagnostics.empty:
                    if view == 0: dx, dy = "x", "z"
                    elif view == 1: dx, dy = "y", "z"
                    else: dx, dy = "x", "y"
                    axes[0, view].scatter(
                        diagnostics[dx], diagnostics[dy], s=42, marker="x",
                        color=colours[label], alpha=.85,
                        label=f"{label} unused 2D")
        axes[0, view].set_title(f"Event {event}: {VIEW_NAMES[view]}")
        axes[0, view].set_xlabel(horizontal.replace("view_", "")+" [mm]")
        axes[0, view].set_ylabel(vertical.replace("view_", "")+" [mm]")
        axes[0, view].legend()
        axes[1, view].set(xlabel="2D-cluster charge", ylabel="Clusters")
        axes[1, view].legend()
    return fig


def plot_hit_level_overlay(run, event=0):
    """Overlay saved 2D inputs, reconstructed 3D hits, and primary truth."""
    projections = ((0, "x", "z"), (1, "y", "z"), (2, "x", "y"))
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    hits2d = run["hits2d"]
    hits2d = hits2d[hits2d.event == event] if not hits2d.empty else hits2d
    hits3d = run["hits"]
    hits3d = hits3d[hits3d.event == event] if not hits3d.empty else hits3d
    truth = run["truth_tracks"]
    truth = truth[(truth.event == event) & (truth.pdg == 2212)
                  & (truth.parent_id == 0)].sort_values(["track_id", "point"])
    for ax, (view, horizontal, vertical) in zip(axes, projections):
        for _, trajectory in truth.groupby("track_id"):
            ax.plot(trajectory[horizontal], trajectory[vertical], color="green",
                    lw=2, label="MC proton")
        selected = hits2d[hits2d.view == view] if not hits2d.empty else hits2d
        used = selected[selected.used.astype(bool)] if not selected.empty else selected
        unused = selected[~selected.used.astype(bool)] if not selected.empty else selected
        if not used.empty:
            ax.scatter(used[horizontal], used[vertical], s=30,
                       facecolors="none", edgecolors="tab:blue",
                       label="used 2D hits")
        if not unused.empty:
            ax.scatter(unused[horizontal], unused[vertical], s=42, marker="x",
                       color="tab:red", label="unused 2D hits")
        if not hits3d.empty:
            ax.scatter(hits3d[horizontal], hits3d[vertical], s=24, marker="s",
                       color="tab:orange", alpha=.7, label="3D hits")
        ax.set(xlabel=f"{horizontal} [mm]", ylabel=f"{vertical} [mm]",
               title=f"Event {event}: {VIEW_NAMES[view]}")
        ax.set_aspect("equal", adjustable="datalim")
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), fontsize=8)
    return fig


def plot_position_distributions(runs, bins=80):
    fig, axes = plt.subplots(2, 3, figsize=(16, 8), constrained_layout=True)
    colours = {"Standard": "tab:blue", "Dominant-axis": "tab:orange"}
    for column, axis in enumerate(("x", "y", "z")):
        for label, run in runs.items():
            axes[0, column].hist(run["hits"][axis], bins=bins, histtype="step",
                                 lw=2, color=colours[label], label=label)
        axes[0, column].set(xlabel=f"3D-hit {axis} [mm]", ylabel="Hits")
        axes[0, column].legend()
        for view, name in VIEW_NAMES.items():
            coordinate = f"view_{axis}"
            frame = next(iter(runs.values()))["views"]
            if coordinate in frame:
                selected = frame[frame.view == view]
                axes[1, column].hist(selected[coordinate], bins=bins,
                                     histtype="step", label=name)
        axes[1, column].set(xlabel=f"Standard 2D-view {axis} [mm]",
                            ylabel="View clusters")
        axes[1, column].legend(fontsize=8)
    return fig


def matched_truth(run):
    segments = run["truth_segments"]
    hits = run["hits"]
    if segments.empty or hits.empty:
        return pd.DataFrame()
    segments = segments[segments.primary_pdg == 2212].copy()
    for axis in "xyz":
        segments[f"truth_{axis}"] = .5*(segments[f"start_{axis}"]+
                                              segments[f"stop_{axis}"])
    grouped = []
    for (event, geom_id), group in segments.groupby(["event", "geom_id"]):
        energy = group.energy_deposit.sum()
        weights = group.energy_deposit.to_numpy()
        if energy > 0:
            position = [np.average(group[f"truth_{axis}"], weights=weights)
                        for axis in "xyz"]
        else:
            position = [group[f"truth_{axis}"].mean() for axis in "xyz"]
        grouped.append((event, geom_id, energy, *position,
                        group.track_length.sum()))
    truth = pd.DataFrame(grouped, columns=["event", "geom_id", "truth_energy",
                                           "truth_x", "truth_y", "truth_z",
                                           "truth_length"])
    return hits.merge(truth, on=["event", "geom_id"], how="inner")


def plot_resolution(runs, bins=80, position_range=(-10, 10),
                    fractional_charge_range=(-2, 2)):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    colours = {"Standard": "tab:blue", "Dominant-axis": "tab:orange"}
    summaries = []
    for label, run in runs.items():
        matched = matched_truth(run)
        for column, axis in enumerate("xyz"):
            residual = matched[axis]-matched[f"truth_{axis}"]
            axes[0, column].hist(residual, bins=bins, range=position_range,
                                 histtype="step", lw=2, color=colours[label],
                                 label=label)
            summaries.append({"method": label, "quantity": f"{axis}_residual_mm",
                              "entries": len(residual), "mean": residual.mean(),
                              "rms": residual.std()})
        valid = matched[matched.truth_energy > 0].copy()
        scale = np.median(valid.charge/valid.truth_energy) if len(valid) else np.nan
        valid["fractional_charge_residual"] = (
            valid.charge/(scale*valid.truth_energy)-1.0)
        axes[1, 0].scatter(valid.truth_energy, valid.charge, s=4, alpha=.2,
                           color=colours[label], label=f"{label}; q/E={scale:.3g}")
        axes[1, 1].hist(valid.fractional_charge_residual, bins=bins,
                        range=fractional_charge_range, histtype="step", lw=2,
                        color=colours[label], label=label)
        axes[1, 2].hist(valid.charge, bins=bins, histtype="step", lw=2,
                        color=colours[label], label=label)
        summaries.append({"method": label, "quantity": "fractional_charge_residual",
                          "entries": len(valid),
                          "mean": valid.fractional_charge_residual.mean(),
                          "rms": valid.fractional_charge_residual.std()})
    for column, axis in enumerate("xyz"):
        axes[0, column].set(xlabel=f"reco {axis} − MC proton {axis} [mm]",
                            ylabel="Matched hits")
        axes[0, column].legend()
    axes[1, 0].set(xlabel="MC deposited energy in virtual cube",
                    ylabel="Reconstructed 3D-hit charge")
    axes[1, 1].set(xlabel="q/(median(q/E) E) − 1", ylabel="Matched hits")
    axes[1, 2].set(xlabel="Reconstructed 3D-hit charge", ylabel="Hits")
    for ax in axes[1]: ax.legend(fontsize=8)
    return fig, pd.DataFrame(summaries)
