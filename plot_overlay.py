#!/usr/bin/env python3
"""Overlay MC truth, fibre hits, reconstructed 3D hits, and track nodes."""

import argparse
from collections import Counter, defaultdict
import csv
from pathlib import Path
import ROOT


TREE_NAMES = ("fiber_hits", "hits3d", "track_nodes", "mc_track_points")
AXES = (("z", "x"), ("z", "y"), ("x", "y"))

# Global ND280 coordinates for the baseline-2024-plusplus detector envelopes.
# The configured local z ranges are shifted by -890 mm when the
# PlusPlusTracker is placed in the global geometry.
COMMON_ENVELOPE = {"x": (-1160.0, 1160.0), "y": (-1177.0, 1237.0)}
DETECTOR_BOUNDS = {
    "SWD": {**COMMON_ENVELOPE, "z": (-890.0, 300.0)},
    "LFGD": {"x": (-1000.0, 1000.0), "y": (-970.0, 1030.0),
             "z": (310.0, 1510.0)},
    "HFGD": {**COMMON_ENVELOPE, "z": (310.0, 1510.0)},
    "SciFi": {**COMMON_ENVELOPE, "z": (1520.0, 1820.0)},
}
BOUNDARY_COLORS = {
    "SWD": ROOT.kCyan + 2,
    "LFGD": ROOT.kOrange + 7,
    "HFGD": ROOT.kOrange + 7,
    "SciFi": ROOT.kMagenta + 1,
}


def event_counts(tree):
    return Counter(int(row.event) for row in tree)


def reconstructed_track_counts(tree):
    """Count tracks rather than nodes in each event."""
    tracks = defaultdict(set)
    for row in tree:
        tracks[int(row.event)].add(int(row.track))
    return Counter({event: len(track_ids)
                    for event, track_ids in tracks.items()})


def points_by_event(tree, selected_events, axes):
    result = defaultdict(list)
    for row in tree:
        event = int(row.event)
        if event in selected_events:
            result[event].append(
                (float(getattr(row, axes[0])), float(getattr(row, axes[1])))
            )
    return result


def transverse_fiber_points_by_event(tree, selected_events, axes,
                                     coordinate_axes):
    """Return fibres only when both displayed coordinates are measured.

    A fibre hit fixes the two coordinates transverse to the fibre.  Its third
    coordinate is merely a geometry reference point along the fibre and must
    not be plotted as though it were a measured hit position.
    """
    result = defaultdict(list)
    displayed_axes = set(axes)
    for row in tree:
        event = int(row.event)
        projection = int(row.projection)
        if event not in selected_events or projection not in coordinate_axes:
            continue
        if displayed_axes != set(coordinate_axes[projection]):
            continue
        result[event].append(
            (float(getattr(row, axes[0])), float(getattr(row, axes[1])))
        )
    return result


def mc_tracks_by_event(tree, selected_events, axes):
    tracks = defaultdict(lambda: defaultdict(list))
    for row in tree:
        event = int(row.event)
        if event not in selected_events:
            continue
        tracks[event][int(row.track_id)].append(
            (int(row.point), float(getattr(row, axes[0])),
             float(getattr(row, axes[1])))
        )
    return {
        event: [[(x, y) for _, x, y in sorted(values)]
                for values in event_tracks.values()]
        for event, event_tracks in tracks.items()
    }


def reco_tracks_by_event(tree, selected_events, axes):
    tracks = defaultdict(lambda: defaultdict(list))
    for row in tree:
        event = int(row.event)
        if event not in selected_events:
            continue
        tracks[event][int(row.track)].append(
            (int(row.node), float(getattr(row, axes[0])),
             float(getattr(row, axes[1])))
        )
    return {
        event: [[(x, y) for _, x, y in sorted(values)]
                for values in event_tracks.values()]
        for event, event_tracks in tracks.items()
    }


def charges_by_event(tree, selected_events):
    result = defaultdict(list)
    for row in tree:
        event = int(row.event)
        if event in selected_events:
            result[event].append(float(row.charge))
    return result


def fiber_coordinate_axes(tree):
    """Find the two transverse coordinates represented by each projection."""
    ranges = {projection: {axis: [float("inf"), float("-inf")]
                           for axis in ("x", "y", "z")}
              for projection in range(3)}
    for row in tree:
        projection = int(row.projection)
        if projection not in ranges:
            continue
        for axis in ("x", "y", "z"):
            value = float(getattr(row, axis))
            ranges[projection][axis][0] = min(ranges[projection][axis][0], value)
            ranges[projection][axis][1] = max(ranges[projection][axis][1], value)
    result = {}
    for projection in range(3):
        # The coordinate along an extended fibre is only a reference value;
        # the two transverse coordinates vary across the fibre grid.
        result[projection] = tuple(sorted(
            ("x", "y", "z"),
            key=lambda axis: ranges[projection][axis][1]
                             - ranges[projection][axis][0],
            reverse=True,
        )[:2])
    return result


def fiber_positions_by_event(tree, selected_events, coordinate_axes):
    result = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in tree:
        event = int(row.event)
        projection = int(row.projection)
        if event not in selected_events or projection not in coordinate_axes:
            continue
        for axis in coordinate_axes[projection]:
            result[event][projection][axis].append(float(getattr(row, axis)))
    return result


def graph(values, color, marker, size):
    # TGraph(0) is ambiguous in PyROOT and can select the TH1-pointer
    # constructor with a null pointer. Construct an empty graph explicitly.
    result = ROOT.TGraph()
    for index, (x, y) in enumerate(values):
        result.SetPoint(index, x, y)
    result.SetMarkerColor(color)
    result.SetLineColor(color)
    result.SetMarkerStyle(marker)
    result.SetMarkerSize(size)
    return result


def draw_detector_boundaries(axes, detector_name):
    """Draw projected detector envelopes and retain their ROOT objects."""
    objects = []
    names = ("SWD", detector_name, "SciFi")
    for name in names:
        bounds = DETECTOR_BOUNDS[name]
        box = ROOT.TBox(
            bounds[axes[0]][0], bounds[axes[1]][0],
            bounds[axes[0]][1], bounds[axes[1]][1])
        box.SetFillStyle(0)
        box.SetLineColor(BOUNDARY_COLORS[name])
        box.SetLineStyle(2)
        box.SetLineWidth(2)
        box.Draw("L SAME")
        objects.append((name, box))
    return objects


def write_index(output_dir, events, counts, track_counts, charge_sums):
    rows = []
    for event in events:
        images = "".join(
            f'<a href="event{event}_{a}{b}.png"><img '
            f'src="event{event}_{a}{b}.png" width="23%"></a>'
            for a, b in AXES
        )
        images += (
            f'<a href="event{event}_fiber_charge.png"><img '
            f'src="event{event}_fiber_charge.png" width="23%"></a>'
        )
        images += (
            f'<a href="event{event}_fiber_positions.png"><img '
            f'src="event{event}_fiber_positions.png" width="47%"></a>'
        )
        images += (
            f'<a href="event{event}_hit3d_charge.png"><img '
            f'src="event{event}_hit3d_charge.png" width="23%"></a>'
            f'<a href="event{event}_hit3d_positions.png"><img '
            f'src="event{event}_hit3d_positions.png" width="47%"></a>'
        )
        rows.append(
            f"<h2 id='event{event}'>Event {event}</h2>"
            f"<p>fibres={counts['fiber_hits'][event]}, "
            f"3D hits={counts['hits3d'][event]}, "
            f"tracks={track_counts[event]}, "
            f"track nodes={counts['track_nodes'][event]}, "
            f"MC points={counts['mc_track_points'][event]}; "
            f"sum fibre charge={charge_sums['fiber_hits'][event]:.3g}, "
            f"sum 3D-hit charge={charge_sums['hits3d'][event]:.3g}, "
            f"3D/fibre={charge_sums['ratio'][event]:.3g}</p>{images}"
        )
    document = (
        "<!doctype html><meta charset='utf-8'><title>LFGD event scan</title>"
        "<style>body{font-family:sans-serif}img{margin:0.5%;height:auto}</style>"
        "<h1>LFGD event scan</h1>"
        '<p><a href="fiber_charge_selected.png">Selected-event aggregate '
        "fibre-charge distribution</a> | "
        '<a href="fiber_positions_selected.png">Selected-event aggregate '
        "fibre-position distributions</a> | "
        '<a href="hit3d_charge_selected.png">Selected-event aggregate '
        "3D-hit charge distribution</a> | "
        '<a href="hit3d_positions_selected.png">Selected-event aggregate '
        "3D-hit position distributions</a> | "
        '<a href="charge_sums_selected.png">Event charge-sum comparison</a> | '
        '<a href="charge_sums.csv">Charge sums as CSV</a></p>'
        + "\n".join(rows)
    )
    (output_dir / "index.html").write_text(document)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--output-dir", type=Path, default=Path("plots"))
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--event", type=int, help="plot one event (default: 0)")
    selection.add_argument("--event-range", type=int, nargs=2,
                           metavar=("FIRST", "LAST"),
                           help="plot an inclusive event range")
    parser.add_argument(
        "--only-with-data",
        choices=("any", "fibers", "hits3d", "tracks"),
        help="within the selection, keep only events with this reco data",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root_file = ROOT.TFile.Open(args.input)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {args.input}")

    trees = {}
    counts = {}
    for name in TREE_NAMES:
        trees[name] = root_file[name]
        if not trees[name]:
            raise RuntimeError(f"Missing TTree '{name}' in {args.input}")
        counts[name] = event_counts(trees[name])
        print(f"{name}: {trees[name].GetEntries()} total entries, "
              f"{len(counts[name])} events")
    track_counts = reconstructed_track_counts(trees["track_nodes"])
    print(f"Reconstructed tracks: {sum(track_counts.values())} total in "
          f"{len(track_counts)} events")

    if args.event_range:
        first, last = args.event_range
        if last < first:
            parser.error("LAST must be greater than or equal to FIRST")
        events = list(range(first, last + 1))
    else:
        events = [args.event if args.event is not None else 0]

    filter_tree = {
        "any": None,
        "fibers": "fiber_hits",
        "hits3d": "hits3d",
        "tracks": "track_nodes",
    }.get(args.only_with_data)
    if args.only_with_data == "any":
        events = [event for event in events if any(
            counts[name][event] for name in
            ("fiber_hits", "hits3d", "track_nodes")
        )]
    elif filter_tree:
        events = [event for event in events if counts[filter_tree][event]]

    if not events:
        raise RuntimeError("No events match the requested selection")
    selected_events = set(events)
    if len(events) <= 20:
        print(f"Plotting {len(events)} event(s): "
              + ", ".join(str(event) for event in events))
    else:
        print(f"Plotting {len(events)} event(s): "
              f"{events[0]} through {events[-1]}")

    ROOT.gROOT.SetBatch(True)
    fiber_charges = charges_by_event(trees["fiber_hits"], selected_events)
    hit3d_charges = charges_by_event(trees["hits3d"], selected_events)
    charge_sums = {
        "fiber_hits": defaultdict(float),
        "hits3d": defaultdict(float),
        "ratio": defaultdict(float),
    }
    for event in events:
        charge_sums["fiber_hits"][event] = sum(fiber_charges[event])
        charge_sums["hits3d"][event] = sum(hit3d_charges[event])
        fiber_sum = charge_sums["fiber_hits"][event]
        charge_sums["ratio"][event] = (
            charge_sums["hits3d"][event] / fiber_sum if fiber_sum else 0.0)

    with (args.output_dir / "charge_sums.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("event", "fiber_charge_sum", "hit3d_charge_sum",
                         "hit3d_over_fiber", "hit3d_exceeds_fiber"))
        for event in events:
            fiber_sum = charge_sums["fiber_hits"][event]
            hit3d_sum = charge_sums["hits3d"][event]
            writer.writerow((event, fiber_sum, hit3d_sum,
                             charge_sums["ratio"][event],
                             int(hit3d_sum > fiber_sum + 1e-9)))

    charge_sum_graph = ROOT.TGraph()
    maximum_charge_sum = 1.0
    violations = []
    for index, event in enumerate(events):
        fiber_sum = charge_sums["fiber_hits"][event]
        hit3d_sum = charge_sums["hits3d"][event]
        charge_sum_graph.SetPoint(index, fiber_sum, hit3d_sum)
        maximum_charge_sum = max(maximum_charge_sum, fiber_sum, hit3d_sum)
        if hit3d_sum > fiber_sum + 1e-9:
            violations.append(event)
    maximum_charge_sum *= 1.05
    charge_sum_canvas = ROOT.TCanvas(
        "c_charge_sums_selected", "c_charge_sums_selected", 850, 750)
    charge_sum_canvas.DrawFrame(
        0.0, 0.0, maximum_charge_sum, maximum_charge_sum,
        "Selected events: charge conservation;"
        "sum of fibre-hit charges;sum of reconstructed 3D-hit charges")
    equality = ROOT.TLine(0.0, 0.0, maximum_charge_sum, maximum_charge_sum)
    equality.SetLineColor(ROOT.kRed + 1)
    equality.SetLineStyle(2)
    equality.SetLineWidth(2)
    equality.Draw("SAME")
    charge_sum_graph.SetMarkerStyle(20)
    charge_sum_graph.SetMarkerSize(0.7)
    charge_sum_graph.SetMarkerColor(ROOT.kBlue + 1)
    charge_sum_graph.Draw("P SAME")
    charge_sum_canvas.SaveAs(str(args.output_dir / "charge_sums_selected.png"))
    charge_sum_canvas.Close()
    print(f"Charge-sum check: {len(violations)} event(s) have a 3D-hit "
          "charge sum larger than the fibre-hit charge sum")
    if violations and len(violations) <= 20:
        print("  Events: " + ", ".join(str(event) for event in violations))
    selected_charges = [charge for event in events
                        for charge in fiber_charges[event]]
    def save_charge_histogram(charges, name, title, output,
                              bins=25, low=0.0, high=100.0):
        histogram = ROOT.TH1D(name, title, bins, low, high)
        histogram.SetDirectory(0)
        histogram.SetLineColor(ROOT.kGray + 2)
        histogram.SetFillColor(ROOT.kAzure - 9)
        histogram.SetLineWidth(2)
        for charge in charges:
            histogram.Fill(charge)
        canvas = ROOT.TCanvas(f"c_{name}", f"c_{name}", 900, 700)
        histogram.Draw("HIST")
        canvas.SaveAs(str(output))
        canvas.Close()

    save_charge_histogram(
        selected_charges,
        "fiber_charge_selected",
        "Selected events: fibre-hit charge;charge;fibre hits",
        args.output_dir / "fiber_charge_selected.png",
    )
    for event in events:
        save_charge_histogram(
            fiber_charges[event],
            f"fiber_charge_event{event}",
            f"Event {event}: fibre-hit charge;charge;fibre hits",
            args.output_dir / f"event{event}_fiber_charge.png",
        )

    selected_hit3d_charges = [charge for event in events
                              for charge in hit3d_charges[event]]
    hit3d_charge_high = max(selected_hit3d_charges, default=1.0)
    hit3d_charge_high = max(1.0, 1.05 * hit3d_charge_high)
    save_charge_histogram(
        selected_hit3d_charges,
        "hit3d_charge_selected",
        "Selected events: reconstructed 3D-hit charge;charge;3D hits",
        args.output_dir / "hit3d_charge_selected.png",
        bins=60, high=hit3d_charge_high,
    )
    for event in events:
        save_charge_histogram(
            hit3d_charges[event],
            f"hit3d_charge_event{event}",
            f"Event {event}: reconstructed 3D-hit charge;charge;3D hits",
            args.output_dir / f"event{event}_hit3d_charge.png",
            bins=60, high=hit3d_charge_high,
        )

    hit3d_positions = {
        axis: defaultdict(list) for axis in ("x", "y", "z")
    }
    for row in trees["hits3d"]:
        event = int(row.event)
        if event not in selected_events:
            continue
        for axis in ("x", "y", "z"):
            hit3d_positions[axis][event].append(float(getattr(row, axis)))

    def save_hit3d_positions(values_by_axis, name, title, output):
        canvas = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1500, 500)
        canvas.Divide(3, 1)
        histograms = []
        for pad, axis in enumerate(("x", "y", "z"), 1):
            values = values_by_axis[axis]
            low = min(values, default=-1.0)
            high = max(values, default=1.0)
            if high <= low:
                low -= 0.5
                high += 0.5
            margin = 0.03 * (high - low)
            histogram = ROOT.TH1D(
                f"{name}_{axis}",
                f"{title}: {axis};{axis} [mm];3D hits",
                60, low - margin, high + margin)
            histogram.SetDirectory(0)
            histogram.SetLineColor(ROOT.kBlue + 2)
            histogram.SetFillColor(ROOT.kAzure - 9)
            for value in values:
                histogram.Fill(value)
            canvas.cd(pad)
            histogram.Draw("HIST")
            histograms.append(histogram)
        canvas.SaveAs(str(output))
        canvas.Close()

    selected_hit3d_positions = {
        axis: [value for event in events
               for value in hit3d_positions[axis][event]]
        for axis in ("x", "y", "z")
    }
    save_hit3d_positions(
        selected_hit3d_positions,
        "hit3d_positions_selected",
        "Selected events: reconstructed 3D-hit position",
        args.output_dir / "hit3d_positions_selected.png",
    )
    for event in events:
        save_hit3d_positions(
            {axis: hit3d_positions[axis][event] for axis in ("x", "y", "z")},
            f"hit3d_positions_event{event}",
            f"Event {event}: reconstructed 3D-hit position",
            args.output_dir / f"event{event}_hit3d_positions.png",
        )

    coordinate_axes = fiber_coordinate_axes(trees["fiber_hits"])
    # HFGD has the wider transverse fibre grid.  This avoids relying on the
    # directory name when a flat tree has been copied or renamed.
    transverse_extent = max(
        abs(float(getattr(row, axis)))
        for row in trees["fiber_hits"] for axis in ("x", "y"))
    detector_name = "HFGD" if transverse_extent > 1050.0 else "LFGD"
    print(f"Detector boundary overlay: {detector_name}")
    fiber_positions = fiber_positions_by_event(
        trees["fiber_hits"], selected_events, coordinate_axes)

    def save_fiber_positions(position_values, name, title, output):
        canvas = ROOT.TCanvas(
            f"c_{name}", f"c_{name}", 1200, 1200)
        canvas.Divide(2, 3)
        histograms = []
        pad = 1
        for projection in range(3):
            for axis in coordinate_axes[projection]:
                values = position_values[projection][axis]
                low = min(values, default=-1.0)
                high = max(values, default=1.0)
                if high <= low:
                    low -= 0.5
                    high += 0.5
                margin = 0.03 * (high - low)
                histogram = ROOT.TH1D(
                    f"{name}_p{projection}_{axis}",
                    f"{title}, fibre projection {projection}: {axis};"
                    f"{axis} [mm];fibre hits",
                    60, low - margin, high + margin)
                histogram.SetDirectory(0)
                histogram.SetLineColor(ROOT.kGray + 2)
                histogram.SetFillColor(ROOT.kOrange - 3)
                for value in values:
                    histogram.Fill(value)
                canvas.cd(pad)
                histogram.Draw("HIST")
                histograms.append(histogram)
                pad += 1
        canvas.SaveAs(str(output))
        canvas.Close()

    selected_positions = defaultdict(lambda: defaultdict(list))
    for event in events:
        for projection in range(3):
            for axis in coordinate_axes[projection]:
                selected_positions[projection][axis].extend(
                    fiber_positions[event][projection][axis])
    save_fiber_positions(
        selected_positions,
        "fiber_positions_selected",
        "Selected events",
        args.output_dir / "fiber_positions_selected.png",
    )
    for event in events:
        save_fiber_positions(
            fiber_positions[event],
            f"fiber_positions_event{event}",
            f"Event {event}",
            args.output_dir / f"event{event}_fiber_positions.png",
        )

    for axes in AXES:
        fibers_by_event = transverse_fiber_points_by_event(
            trees["fiber_hits"], selected_events, axes, coordinate_axes)
        hits_by_event = points_by_event(trees["hits3d"], selected_events, axes)
        tracks_by_event = reco_tracks_by_event(
            trees["track_nodes"], selected_events, axes)
        truth_by_event = mc_tracks_by_event(
            trees["mc_track_points"], selected_events, axes)

        for event in events:
            fibers = graph(fibers_by_event[event], ROOT.kGray + 1, 7, 0.5)
            hits = graph(hits_by_event[event], ROOT.kBlue + 1, 20, 0.8)
            tracks = [graph(values, ROOT.kRed + 1, 24, 1.0)
                      for values in tracks_by_event.get(event, [])]
            truth = [graph(values, ROOT.kGreen + 2, 1, 0.0)
                     for values in truth_by_event.get(event, [])]

            canvas_name = f"c_{event}_{axes[0]}{axes[1]}"
            canvas = ROOT.TCanvas(canvas_name, canvas_name, 1000, 750)
            title = (
                f"Event {event}: fibres shown={len(fibers_by_event[event])}/"
                f"{counts['fiber_hits'][event]}, "
                f"3D hits={counts['hits3d'][event]}, "
                f"tracks={track_counts[event]}, "
                f"track nodes={counts['track_nodes'][event]};"
                f"{axes[0]} [mm];{axes[1]} [mm]"
            )
            canvas.DrawFrame(-3000, -1500, 3000, 1500, title)
            boundaries = draw_detector_boundaries(axes, detector_name)
            for trajectory in truth:
                trajectory.SetLineWidth(2)
                trajectory.Draw(
                    ("L" if trajectory.GetN() > 1 else "P") + " SAME")
            if fibers.GetN():
                fibers.Draw("P SAME")
            if hits.GetN():
                hits.Draw("P SAME")
            for track in tracks:
                if track.GetN():
                    track.Draw(("LP" if track.GetN() > 1 else "P") + " SAME")

            legend = ROOT.TLegend(0.68, 0.72, 0.91, 0.90)
            if fibers.GetN():
                legend.AddEntry(fibers, "fibre hits", "p")
            if hits.GetN():
                legend.AddEntry(hits, "3D reconstructed hits", "p")
            if tracks:
                legend.AddEntry(tracks[0], "reconstructed tracks", "lp")
            if truth:
                legend.AddEntry(truth[0], "MC truth trajectories", "l")
            for name, boundary in boundaries:
                legend.AddEntry(boundary, f"{name} envelope", "l")
            legend.Draw()
            output = args.output_dir / f"event{event}_{axes[0]}{axes[1]}.png"
            canvas.SaveAs(str(output))
            canvas.Close()

    write_index(args.output_dir, events, counts, track_counts, charge_sums)
    print(f"Browse {args.output_dir / 'index.html'}")


if __name__ == "__main__":
    main()
