#!/usr/bin/env python3
"""Overlay MC truth, fibre hits, reconstructed 3D hits, and track nodes."""

import argparse
from collections import Counter, defaultdict
import csv
import math
from pathlib import Path
import statistics
import ROOT


TREE_NAMES = ("fiber_hits", "hits3d", "hit3d_views", "track_nodes",
              "track_node_hits", "mc_virtual_segments")
AXES = (("z", "x"), ("z", "y"), ("x", "y"))

# Global ND280 coordinates for the baseline-2024-plusplus detector envelopes.
# The configured local z ranges are shifted by -890 mm when the
# PlusPlusTracker is placed in the global geometry.
COMMON_ENVELOPE = {"x": (-1160.0, 1160.0), "y": (-1177.0, 1237.0)}
DETECTOR_BOUNDS = {
    "SWD": {**COMMON_ENVELOPE, "z": (-890.0, 300.0)},
    # LFGD and HFGD are alternative 2200 x 2300 x 1200 mm active volumes
    # at the same placement.  The tracker contributes the +30 mm global-Y
    # offset and -890 mm global-Z offset.
    "LFGD": {"x": (-1100.0, 1100.0), "y": (-1120.0, 1180.0),
             "z": (310.0, 1510.0)},
    "HFGD": {"x": (-1100.0, 1100.0), "y": (-1120.0, 1180.0),
             "z": (310.0, 1510.0)},
    "SciFi": {**COMMON_ENVELOPE, "z": (1520.0, 1820.0)},
}
BOUNDARY_COLORS = {
    "SWD": ROOT.kCyan + 2,
    "LFGD": ROOT.kOrange + 7,
    "HFGD": ROOT.kOrange + 7,
    "SciFi": ROOT.kMagenta + 1,
}
VOXEL_PITCH_MM = 10.0


def clip_segment_to_box(start, stop, bounds):
    """Clip a segment to an axis-aligned detector box (Liang-Barsky)."""
    delta = tuple(stop[i] - start[i] for i in range(3))
    first, last = 0.0, 1.0
    for index, axis in enumerate(("x", "y", "z")):
        low, high = bounds[axis]
        if abs(delta[index]) < 1e-12:
            if start[index] < low or start[index] > high:
                return None
            continue
        enter = (low - start[index]) / delta[index]
        leave = (high - start[index]) / delta[index]
        if enter > leave:
            enter, leave = leave, enter
        first = max(first, enter)
        last = min(last, leave)
        if first > last:
            return None
    return (
        tuple(start[i] + first * delta[i] for i in range(3)),
        tuple(start[i] + last * delta[i] for i in range(3)),
    )


def position_to_voxel(position, bounds, pitch=VOXEL_PITCH_MM):
    indices = []
    for value, axis in zip(position, ("x", "y", "z")):
        low, high = bounds[axis]
        if value < low - 1e-7 or value > high + 1e-7:
            return None
        count = int(round((high - low) / pitch))
        index = int(math.floor((value - low) / pitch))
        # A point on the positive envelope belongs to the final voxel.
        index = min(max(index, 0), count - 1)
        indices.append(index)
    return tuple(indices)


def segment_voxels(start, stop, bounds, pitch=VOXEL_PITCH_MM):
    """Return every voxel crossed by a clipped straight trajectory segment."""
    clipped = clip_segment_to_box(start, stop, bounds)
    if not clipped:
        return set()
    start, stop = clipped
    delta = tuple(stop[i] - start[i] for i in range(3))
    crossings = {0.0, 1.0}
    for index, axis in enumerate(("x", "y", "z")):
        if abs(delta[index]) < 1e-12:
            continue
        low, high = bounds[axis]
        plane = low + pitch
        while plane < high - 1e-9:
            fraction = (plane - start[index]) / delta[index]
            if 1e-10 < fraction < 1.0 - 1e-10:
                crossings.add(fraction)
            plane += pitch
    crossings = sorted(crossings)
    voxels = set()
    for first, last in zip(crossings, crossings[1:]):
        fraction = 0.5 * (first + last)
        position = tuple(start[i] + fraction * delta[i] for i in range(3))
        voxel = position_to_voxel(position, bounds, pitch)
        if voxel is not None:
            voxels.add(voxel)
    return voxels


def mc_voxels_by_event(tree, selected_events, bounds):
    tracks = defaultdict(lambda: defaultdict(list))
    for row in tree:
        event = int(row.event)
        if event not in selected_events or int(row.parent_id) != 0:
            continue
        tracks[event][int(row.track_id)].append(
            (int(row.point), (float(row.x), float(row.y), float(row.z))))
    result = defaultdict(set)
    for event, event_tracks in tracks.items():
        for points in event_tracks.values():
            ordered = [position for _, position in sorted(points)]
            for start, stop in zip(ordered, ordered[1:]):
                result[event].update(segment_voxels(start, stop, bounds))
    return result


def reconstructed_voxels_by_event(tree, selected_events, bounds):
    result = defaultdict(set)
    for row in tree:
        event = int(row.event)
        if event not in selected_events:
            continue
        voxel = position_to_voxel(
            (float(row.x), float(row.y), float(row.z)), bounds)
        if voxel is not None:
            result[event].add(voxel)
    return result


def mc_virtual_voxels_by_event(tree, selected_events):
    """Use the cube indices stored with the exact virtual truth segments."""
    result = defaultdict(set)
    for row in tree:
        event = int(row.event)
        if event in selected_events:
            result[event].add(
                (int(row.cube_x), int(row.cube_y), int(row.cube_z)))
    return result


def mc_virtual_segments_by_event(tree, selected_events, axes):
    """Return projected entry/exit pairs for each virtual-cube segment."""
    result = defaultdict(list)
    for row in tree:
        event = int(row.event)
        if event not in selected_events:
            continue
        result[event].append((
            (float(getattr(row, f"start_{axes[0]}")),
             float(getattr(row, f"start_{axes[1]}"))),
            (float(getattr(row, f"stop_{axes[0]}")),
             float(getattr(row, f"stop_{axes[1]}"))),
        ))
    return result


def mc_cube_truth_by_event(tree, selected_events):
    """Combine exact segments into one energy-weighted truth point per cube."""
    accumulators = defaultdict(lambda: defaultdict(lambda: {
        "energy": 0.0, "length": 0.0,
        "energy_position": [0.0, 0.0, 0.0],
        "length_position": [0.0, 0.0, 0.0],
        "segments": 0,
    }))
    for row in tree:
        event = int(row.event)
        if event not in selected_events:
            continue
        voxel = (int(row.cube_x), int(row.cube_y), int(row.cube_z))
        energy = max(0.0, float(row.energy_deposit))
        length = max(0.0, float(row.track_length))
        midpoint = tuple(
            0.5 * (float(getattr(row, f"start_{axis}"))
                   + float(getattr(row, f"stop_{axis}")))
            for axis in ("x", "y", "z")
        )
        value = accumulators[event][voxel]
        value["energy"] += energy
        value["length"] += length
        value["segments"] += 1
        for index in range(3):
            value["energy_position"][index] += energy * midpoint[index]
            value["length_position"][index] += length * midpoint[index]

    result = defaultdict(dict)
    for event, voxels in accumulators.items():
        for voxel, value in voxels.items():
            if value["energy"] > 0.0:
                position = tuple(
                    coordinate / value["energy"]
                    for coordinate in value["energy_position"])
                weighting = "energy"
            elif value["length"] > 0.0:
                position = tuple(
                    coordinate / value["length"]
                    for coordinate in value["length_position"])
                weighting = "length"
            else:
                position = (0.0, 0.0, 0.0)
                weighting = "none"
            result[event][voxel] = {
                "position": position,
                "energy": value["energy"],
                "length": value["length"],
                "segments": value["segments"],
                "weighting": weighting,
            }
    return result


def reco_hits_by_voxel(tree, selected_events, bounds):
    """Return the closest continuous reconstructed hit for every voxel."""
    result = defaultdict(dict)
    for row in tree:
        event = int(row.event)
        if event not in selected_events:
            continue
        position = (float(row.x), float(row.y), float(row.z))
        voxel = position_to_voxel(position, bounds)
        if voxel is None:
            continue
        centre = tuple(
            bounds[axis][0] + (voxel[index] + 0.5) * VOXEL_PITCH_MM
            for index, axis in enumerate(("x", "y", "z")))
        distance2 = sum((position[index]-centre[index])**2
                        for index in range(3))
        previous = result[event].get(voxel)
        if previous is None or distance2 < previous["centre_distance2"]:
            result[event][voxel] = {
                "position": position,
                "charge": float(row.charge),
                "centre_distance2": distance2,
            }
    return result


def observed_fibres_by_event(tree, selected_events):
    result = defaultdict(lambda: defaultdict(set))
    for row in tree:
        event = int(row.event)
        if event in selected_events:
            result[event][int(row.projection)].add((int(row.u), int(row.v)))
    return result


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


def voxel_projection_points(voxels, bounds, axes):
    """Project virtual-cube indices to their global centre coordinates."""
    points = set()
    for voxel in voxels:
        coordinates = {
            axis: bounds[axis][0]
                  + (voxel[index] + 0.5) * VOXEL_PITCH_MM
            for index, axis in enumerate(("x", "y", "z"))
        }
        points.add((coordinates[axes[0]], coordinates[axes[1]]))
    return sorted(points)


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


def write_index(output_dir, events, counts, track_counts, charge_sums,
                diagnostics):
    rows = []
    for event in events:
        images = "".join(
            f'<a target="_blank" href="event{event}_{a}{b}.png"><img '
            f'src="event{event}_{a}{b}.png" width="23%"></a>'
            for a, b in AXES
        )
        images += "".join(
            f'<a target="_blank" href="event{event}_{a}{b}_voxel_match.png"><img '
            f'src="event{event}_{a}{b}_voxel_match.png" width="23%"></a>'
            for a, b in AXES
        )
        images += "".join(
            f'<a target="_blank" href="event{event}_{a}{b}_track_match.png"><img '
            f'src="event{event}_{a}{b}_track_match.png" width="23%"></a>'
            for a, b in AXES
        )
        images += (
            f'<a target="_blank" href="event{event}_fiber_charge.png"><img '
            f'src="event{event}_fiber_charge.png" width="23%"></a>'
        )
        images += (
            f'<a target="_blank" href="event{event}_fiber_positions.png"><img '
            f'src="event{event}_fiber_positions.png" width="47%"></a>'
        )
        images += (
            f'<a target="_blank" href="event{event}_hit3d_charge.png"><img '
            f'src="event{event}_hit3d_charge.png" width="23%"></a>'
            f'<a target="_blank" href="event{event}_hit3d_track_charge.png"><img '
            f'src="event{event}_hit3d_track_charge.png" width="23%"></a>'
            f'<a target="_blank" href="event{event}_hit3d_view_charge.png"><img '
            f'src="event{event}_hit3d_view_charge.png" width="47%"></a>'
            f'<a target="_blank" href="event{event}_hit3d_fiber_count.png"><img '
            f'src="event{event}_hit3d_fiber_count.png" width="47%"></a>'
            f'<a target="_blank" href="event{event}_hit3d_view_reuse.png"><img '
            f'src="event{event}_hit3d_view_reuse.png" width="47%"></a>'
            f'<a target="_blank" href="event{event}_hit3d_positions.png"><img '
            f'src="event{event}_hit3d_positions.png" width="47%"></a>'
        )
        rows.append(
            f"<h2 id='event{event}'>Event {event}</h2>"
            f"<p>fibres={counts['fiber_hits'][event]}, "
            f"3D hits={counts['hits3d'][event]}, "
            f"tracks={track_counts[event]}, "
            f"track nodes={counts['track_nodes'][event]}, "
            f"MC virtual-cube segments="
            f"{counts['mc_virtual_segments'][event]}; "
            f"measured fibre-charge sum="
            f"{charge_sums['fiber_hits'][event]:.3g}, "
            f"fitted 3D-deposit sum={charge_sums['hits3d'][event]:.3g}, "
            f"deposit/measured={charge_sums['ratio'][event]:.3g} "
            f"(response-model dependent); "
            f"MC voxels={diagnostics[event]['expected']}, "
            f"matched 3D voxels={diagnostics[event]['matched']}, "
            f"missing MC voxels={diagnostics[event]['missing']}, "
            f"off-track 3D voxels={diagnostics[event]['offtrack']}, "
            f"node-associated 3D-hit voxels={diagnostics[event]['nodes']}, "
            f"3D voxels represented by nodes="
            f"{len(diagnostics[event]['node_matched_reco_voxels'])}, "
            f"3D voxels without nodes="
            f"{len(diagnostics[event]['reco_without_node_voxels'])}</p>{images}"
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
        '<a href="hit3d_track_charge_selected.png">Used versus unused 3D-hit '
        "charge</a> | "
        '<a href="hit3d_view_charge_selected.png">3D-hit charge composition '
        "by view</a> | "
        '<a href="hit3d_fiber_count_selected.png">Contributing fibres by '
        "view</a> | "
        '<a href="hit3d_view_reuse_selected.png">2D-view reuse by 3D '
        "candidates</a> | "
        '<a href="hit3d_view_composition.csv">3D-hit view composition CSV</a> | '
        '<a href="hit3d_positions_selected.png">Selected-event aggregate '
        "3D-hit position distributions</a> | "
        '<a href="charge_sums_selected.png">Measured-charge versus fitted-'
        'deposit comparison</a> | '
        '<a href="charge_sums.csv">Charge sums as CSV</a> | '
        '<a href="reconstruction_efficiency.png">MC voxel/reconstruction '
        'comparison</a> | '
        '<a href="reconstruction_efficiency.csv">Reconstruction diagnostics '
        'as CSV</a> | '
        '<a href="voxel_position_residuals.png">MC/reco position residuals</a> | '
        '<a href="voxel_charge_response.png">MC-energy/3D-deposit response</a> | '
        '<a href="voxel_residuals.csv">Per-voxel residuals as CSV</a></p>'
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
        writer.writerow(("event", "measured_fiber_charge_sum",
                         "fitted_3d_deposit_sum", "deposit_over_measured"))
        for event in events:
            fiber_sum = charge_sums["fiber_hits"][event]
            hit3d_sum = charge_sums["hits3d"][event]
            writer.writerow((event, fiber_sum, hit3d_sum,
                             charge_sums["ratio"][event]))

    charge_sum_graph = ROOT.TGraph()
    maximum_charge_sum = 1.0
    for index, event in enumerate(events):
        fiber_sum = charge_sums["fiber_hits"][event]
        hit3d_sum = charge_sums["hits3d"][event]
        charge_sum_graph.SetPoint(index, fiber_sum, hit3d_sum)
        maximum_charge_sum = max(maximum_charge_sum, fiber_sum, hit3d_sum)
    maximum_charge_sum *= 1.05
    charge_sum_canvas = ROOT.TCanvas(
        "c_charge_sums_selected", "c_charge_sums_selected", 850, 750)
    charge_sum_canvas.DrawFrame(
        0.0, 0.0, maximum_charge_sum, maximum_charge_sum,
        "Measured fibre charge versus fitted 3D deposit;"
        "sum of measured fibre-hit charge;"
        "sum of fitted reconstructed deposit")
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
    ratios = [charge_sums["ratio"][event] for event in events]
    print("Charge comparison: fitted 3D values include the configured fibre "
          "response model; deposit/measured is directly comparable only "
          "when that response is unity; "
          f"range={min(ratios):.3g}--{max(ratios):.3g}")
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

    # Classify reconstructed hits using the exact hit identifiers attached to
    # node measurements.  Fitted node coordinates are intentionally not used.
    node_hit_ids = defaultdict(set)
    for row in trees["track_node_hits"]:
        event = int(row.event)
        if event in selected_events:
            node_hit_ids[event].add(int(row.geom_id))
    used_hit3d_charges = defaultdict(list)
    unused_hit3d_charges = defaultdict(list)
    for row in trees["hits3d"]:
        event = int(row.event)
        if event not in selected_events:
            continue
        target = (used_hit3d_charges if int(row.geom_id) in node_hit_ids[event]
                  else unused_hit3d_charges)
        target[event].append(float(row.charge))

    def save_hit3d_track_charge(found, not_found, name, title, output):
        high = 1.05 * max([1.0] + found + not_found)
        found_hist = ROOT.TH1D(
            f"{name}_found", f"{title};3D-hit charge;3D hits", 60, 0.0, high)
        missing_hist = ROOT.TH1D(
            f"{name}_not_found", f"{title};3D-hit charge;3D hits", 60, 0.0,
            high)
        for histogram in (found_hist, missing_hist):
            histogram.SetDirectory(0)
            histogram.SetLineWidth(2)
            histogram.SetStats(0)
        found_hist.SetLineColor(ROOT.kGreen + 2)
        missing_hist.SetLineColor(ROOT.kRed + 1)
        for charge in found:
            found_hist.Fill(charge)
        for charge in not_found:
            missing_hist.Fill(charge)
        found_hist.SetMaximum(1.18 * max(
            found_hist.GetMaximum(), missing_hist.GetMaximum(), 1.0))
        canvas = ROOT.TCanvas(f"c_{name}", f"c_{name}", 900, 700)
        found_hist.Draw("HIST")
        missing_hist.Draw("HIST SAME")
        legend = ROOT.TLegend(0.55, 0.73, 0.90, 0.89)
        legend.AddEntry(
            found_hist, f"used by track nodes ({len(found)})", "l")
        legend.AddEntry(
            missing_hist, f"not used by track nodes ({len(not_found)})", "l")
        legend.Draw()
        canvas.SaveAs(str(output))
        canvas.Close()

    selected_used_charges = [charge for event in events
                             for charge in used_hit3d_charges[event]]
    selected_unused_charges = [charge for event in events
                               for charge in unused_hit3d_charges[event]]
    save_hit3d_track_charge(
        selected_used_charges, selected_unused_charges,
        "hit3d_track_charge_selected",
        "Selected events: 3D-hit use by reconstructed tracks",
        args.output_dir / "hit3d_track_charge_selected.png")
    for event in events:
        save_hit3d_track_charge(
            used_hit3d_charges[event], unused_hit3d_charges[event],
            f"hit3d_track_charge_event{event}",
            f"Event {event}: 3D-hit use by reconstructed tracks",
            args.output_dir / f"event{event}_hit3d_track_charge.png")

    view_composition = defaultdict(lambda: defaultdict(lambda: {
        "used_charge": [], "unused_charge": [],
        "used_fibres": [], "unused_fibres": [],
        "used_reuse": [], "unused_reuse": []}))
    view_count = defaultdict(lambda: defaultdict(set))
    view_reuse = defaultdict(int)
    for row in trees["hit3d_views"]:
        event = int(row.event)
        if event in selected_events:
            view_reuse[(event, int(row.view), int(row.view_geom_id))] += 1
    with (args.output_dir / "hit3d_view_composition.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("event", "hit", "geom_id", "used_by_track", "view",
                         "view_geom_id", "hit3d_charge", "view_charge",
                         "fiber_count", "fiber_charge_sum",
                         "view_candidate_multiplicity"))
        for row in trees["hit3d_views"]:
            event = int(row.event)
            if event not in selected_events:
                continue
            geom_id = int(row.geom_id)
            view = int(row.view)
            used = geom_id in node_hit_ids[event]
            prefix = "used" if used else "unused"
            reuse = view_reuse[(event, view, int(row.view_geom_id))]
            view_composition[event][view][f"{prefix}_charge"].append(
                float(row.view_charge))
            view_composition[event][view][f"{prefix}_fibres"].append(
                int(row.fiber_count))
            view_composition[event][view][f"{prefix}_reuse"].append(reuse)
            view_count[event][geom_id].add(view)
            writer.writerow((
                event, int(row.hit), geom_id, int(used), view,
                int(row.view_geom_id),
                float(row.hit_charge), float(row.view_charge),
                int(row.fiber_count), float(row.fiber_charge_sum), reuse))

    def save_composition(values, field, name, title, output):
        canvas = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1500, 500)
        canvas.Divide(3, 1)
        histograms = []
        for pad, view in enumerate(range(3), 1):
            found = values[view][f"used_{field}"]
            not_found = values[view][f"unused_{field}"]
            if field == "charge":
                low = 0.0
                high = 1.05 * max([1.0] + found + not_found)
                bins = 60
                x_title = "2D-view cluster charge"
            elif field == "fibres":
                maximum = max([1] + found + not_found)
                low, high, bins = -0.5, maximum + 0.5, maximum + 1
                x_title = "contributing fibres in 2D cluster"
            else:
                maximum = max([1] + found + not_found)
                low, high, bins = 0.5, maximum + 0.5, maximum
                x_title = "3D candidates sharing this 2D cluster"
            used_hist = ROOT.TH1D(
                f"{name}_used_{view}",
                f"{title}: view {view};{x_title};3D-hit views",
                bins, low, high)
            unused_hist = ROOT.TH1D(
                f"{name}_unused_{view}",
                f"{title}: view {view};{x_title};3D-hit views",
                bins, low, high)
            for histogram in (used_hist, unused_hist):
                histogram.SetDirectory(0)
                histogram.SetStats(0)
                histogram.SetLineWidth(2)
            used_hist.SetLineColor(ROOT.kGreen + 2)
            unused_hist.SetLineColor(ROOT.kRed + 1)
            for value in found:
                used_hist.Fill(value)
            for value in not_found:
                unused_hist.Fill(value)
            used_hist.SetMaximum(1.2 * max(
                used_hist.GetMaximum(), unused_hist.GetMaximum(), 1.0))
            canvas.cd(pad)
            used_hist.Draw("HIST")
            unused_hist.Draw("HIST SAME")
            legend = ROOT.TLegend(0.48, 0.74, 0.89, 0.89)
            legend.AddEntry(used_hist, f"used ({len(found)})", "l")
            legend.AddEntry(unused_hist, f"unused ({len(not_found)})", "l")
            legend.Draw()
            histograms.extend((used_hist, unused_hist, legend))
        canvas.SaveAs(str(output))
        canvas.Close()

    selected_composition = defaultdict(lambda: {
        "used_charge": [], "unused_charge": [],
        "used_fibres": [], "unused_fibres": [],
        "used_reuse": [], "unused_reuse": []})
    for event in events:
        for view in range(3):
            for field in ("used_charge", "unused_charge", "used_fibres",
                          "unused_fibres", "used_reuse", "unused_reuse"):
                selected_composition[view][field].extend(
                    view_composition[event][view][field])
    for field, stem, label in (
            ("charge", "view_charge", "2D-view charge composition"),
            ("fibres", "fiber_count", "2D-view fibre multiplicity"),
            ("reuse", "view_reuse", "2D-view reuse by 3D candidates")):
        save_composition(
            selected_composition, field, f"hit3d_{stem}_selected",
            f"Selected events: {label}",
            args.output_dir / f"hit3d_{stem}_selected.png")
        for event in events:
            save_composition(
                view_composition[event], field,
                f"hit3d_{stem}_event{event}", f"Event {event}: {label}",
                args.output_dir / f"event{event}_hit3d_{stem}.png")

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

    # Compare reconstruction with the exact virtual-cube truth segments.
    # Unlike interpolation between sparse trajectory points, these indices
    # come directly from the Geant4 step clipping at every 10 mm boundary.
    detector_bounds = DETECTOR_BOUNDS[detector_name]
    mc_voxels = mc_virtual_voxels_by_event(
        trees["mc_virtual_segments"], selected_events)
    reco_voxels = reconstructed_voxels_by_event(
        trees["hits3d"], selected_events, detector_bounds)
    # Use the measurements attached to nodes for hit/node association.  The
    # positions in track_nodes are fitted states and may move into adjacent
    # voxels; they remain the correct coordinates for drawing the red track.
    node_voxels = reconstructed_voxels_by_event(
        trees["track_node_hits"], selected_events, detector_bounds)
    observed_fibres = observed_fibres_by_event(
        trees["fiber_hits"], selected_events)
    axis_index = {"x": 0, "y": 1, "z": 2}

    diagnostics = {}
    with (args.output_dir / "reconstruction_efficiency.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow((
            "event", "mc_voxels", "reco_3d_entries", "reco_unique_voxels",
            "reco_matched_mc_voxels", "reco_missing_mc_voxels",
            "reco_offtrack_voxels", "track_nodes", "track_node_voxels",
            "track_node_matched_mc_voxels", "hit_voxel_efficiency",
            "node_voxel_efficiency", "reco_voxels_with_track_node",
            "reco_voxels_without_track_node", "track_node_only_voxels",
            "expected_view0_rows", "observed_view0_fibres",
            "expected_view1_rows", "observed_view1_fibres",
            "expected_view2_rows", "observed_view2_fibres"))
        for event in events:
            expected = mc_voxels[event]
            reconstructed = reco_voxels[event]
            nodes = node_voxels[event]
            matched = reconstructed & expected
            node_matched = nodes & expected
            missing = expected - reconstructed
            offtrack = reconstructed - expected
            expected_views = {}
            for projection in range(3):
                axes = coordinate_axes[projection]
                indices = tuple(axis_index[axis] for axis in axes)
                expected_views[projection] = {
                    tuple(voxel[index] for index in indices)
                    for voxel in expected
                }
            hit_efficiency = len(matched) / len(expected) if expected else 0.0
            node_efficiency = (
                len(node_matched) / len(expected) if expected else 0.0)
            diagnostics[event] = {
                "expected": len(expected),
                "reco": len(reconstructed),
                "matched": len(matched),
                "missing": len(missing),
                "offtrack": len(offtrack),
                "nodes": len(nodes),
                "hit_efficiency": hit_efficiency,
                "node_efficiency": node_efficiency,
                "expected_voxels": expected,
                "matched_voxels": matched,
                "missing_voxels": missing,
                "offtrack_voxels": offtrack,
                "reco_voxels": reconstructed,
                "node_voxels": nodes,
                "node_matched_reco_voxels": nodes & reconstructed,
                "reco_without_node_voxels": reconstructed - nodes,
                "nodes_without_reco_voxels": nodes - reconstructed,
            }
            writer.writerow((
                event, len(expected), counts["hits3d"][event],
                len(reconstructed), len(matched), len(missing), len(offtrack),
                counts["track_nodes"][event], len(nodes), len(node_matched),
                hit_efficiency, node_efficiency,
                len(reconstructed & nodes), len(reconstructed - nodes),
                len(nodes - reconstructed),
                len(expected_views[0]), len(observed_fibres[event][0]),
                len(expected_views[1]), len(observed_fibres[event][1]),
                len(expected_views[2]), len(observed_fibres[event][2])))

    # Build an energy-weighted true crossing point in every virtual cube and
    # compare it with the continuous reconstructed position assigned to that
    # cube.  MC stores deposited energy (MeV), not PE charge.  A sample-wide
    # median PE/MeV response is therefore reported explicitly and used only
    # for the diagnostic charge residual.
    mc_cube_truth = mc_cube_truth_by_event(
        trees["mc_virtual_segments"], selected_events)
    reco_cube_hits = reco_hits_by_voxel(
        trees["hits3d"], selected_events, detector_bounds)
    charge_responses = []
    for event in events:
        for voxel in mc_cube_truth[event].keys() & reco_cube_hits[event].keys():
            energy = mc_cube_truth[event][voxel]["energy"]
            if energy > 0.0:
                charge_responses.append(
                    reco_cube_hits[event][voxel]["charge"] / energy)
    charge_scale = statistics.median(charge_responses) \
        if charge_responses else 0.0

    position_residuals = {axis: [] for axis in ("dx", "dy", "dz", "dr")}
    charge_points = []
    charge_fractional_residuals = []
    with (args.output_dir / "voxel_residuals.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow((
            "event", "cube_x", "cube_y", "cube_z", "mc_segments",
            "mc_position_weighting", "mc_x", "mc_y", "mc_z",
            "reco_x", "reco_y", "reco_z", "dx", "dy", "dz", "dr",
            "mc_energy_deposit_MeV", "reco_fitted_deposit",
            "sample_median_deposit_per_MeV", "expected_deposit_from_mc",
            "deposit_residual", "fractional_deposit_residual"))
        for event in events:
            shared = sorted(
                mc_cube_truth[event].keys() & reco_cube_hits[event].keys())
            for voxel in shared:
                truth = mc_cube_truth[event][voxel]
                reco = reco_cube_hits[event][voxel]
                delta = tuple(reco["position"][index]
                              - truth["position"][index]
                              for index in range(3))
                distance = math.sqrt(sum(value*value for value in delta))
                energy = truth["energy"]
                charge = reco["charge"]
                expected_charge = charge_scale * energy
                charge_residual = charge - expected_charge
                fractional = (charge_residual / expected_charge
                              if expected_charge > 0.0 else float("nan"))
                position_residuals["dx"].append(delta[0])
                position_residuals["dy"].append(delta[1])
                position_residuals["dz"].append(delta[2])
                position_residuals["dr"].append(distance)
                if energy > 0.0:
                    charge_points.append((energy, charge))
                    charge_fractional_residuals.append(fractional)
                writer.writerow((
                    event, *voxel, truth["segments"], truth["weighting"],
                    *truth["position"], *reco["position"], *delta, distance,
                    energy, charge, charge_scale, expected_charge,
                    charge_residual, fractional))

    position_canvas = ROOT.TCanvas(
        "c_voxel_position_residuals", "c_voxel_position_residuals", 1200, 900)
    position_canvas.Divide(2, 2)
    position_histograms = []
    for pad, (key, label) in enumerate((
            ("dx", "reco x - MC x [mm]"),
            ("dy", "reco y - MC y [mm]"),
            ("dz", "reco z - MC z [mm]"),
            ("dr", "|reco - MC| [mm]")), 1):
        values = position_residuals[key]
        low = min(values, default=-1.0)
        high = max(values, default=1.0)
        if high <= low:
            low -= 0.5
            high += 0.5
        histogram = ROOT.TH1D(
            f"voxel_residual_{key}",
            f"Matched virtual cubes;{label};cubes", 80, low, high)
        histogram.SetDirectory(0)
        for value in values:
            histogram.Fill(value)
        position_canvas.cd(pad)
        histogram.Draw("HIST")
        position_histograms.append(histogram)
    position_canvas.SaveAs(
        str(args.output_dir / "voxel_position_residuals.png"))
    position_canvas.Close()

    charge_canvas = ROOT.TCanvas(
        "c_voxel_charge_response", "c_voxel_charge_response", 1200, 550)
    charge_canvas.Divide(2, 1)
    charge_graph = graph(charge_points, ROOT.kBlue + 1, 20, 0.55)
    max_energy = max([1.0] + [point[0] for point in charge_points])
    max_charge = max([1.0] + [point[1] for point in charge_points])
    charge_canvas.cd(1)
    ROOT.gPad.DrawFrame(
        0.0, 0.0, 1.05*max_energy, 1.05*max_charge,
        "Matched virtual cubes;MC deposited energy [MeV];"
        "reconstructed fitted deposit")
    if charge_graph.GetN():
        charge_graph.Draw("P SAME")
    response_line = ROOT.TLine(
        0.0, 0.0, max_energy, charge_scale*max_energy)
    response_line.SetLineColor(ROOT.kRed + 1)
    response_line.SetLineWidth(2)
    response_line.Draw("SAME")
    finite_fractional = [value for value in charge_fractional_residuals
                         if math.isfinite(value)]
    # Suppress only the display range of extreme ratios caused by tiny MC
    # deposits.  Every untruncated residual remains available in the CSV.
    sorted_fractional = sorted(finite_fractional)
    frac_low = -1.0
    frac_high = (sorted_fractional[
        min(len(sorted_fractional)-1,
            int(0.98*len(sorted_fractional)))]
        if sorted_fractional else 1.0)
    frac_high = max(1.0, frac_high)
    if frac_high <= frac_low:
        frac_low -= 0.5
        frac_high += 0.5
    charge_residual_hist = ROOT.TH1D(
        "voxel_charge_fractional_residual",
        "Matched virtual cubes;"
        "(fitted deposit - scaled MC energy)/(scaled MC energy);cubes",
        80, frac_low, frac_high)
    charge_residual_hist.SetDirectory(0)
    for value in finite_fractional:
        charge_residual_hist.Fill(value)
    charge_canvas.cd(2)
    charge_residual_hist.Draw("HIST")
    charge_canvas.SaveAs(str(args.output_dir / "voxel_charge_response.png"))
    charge_canvas.Close()
    print(f"Voxel residuals: {len(position_residuals['dr'])} matched cubes; "
          f"median fitted-deposit/MC-energy scale={charge_scale:.6g}")

    efficiency_canvas = ROOT.TCanvas(
        "c_reconstruction_efficiency", "c_reconstruction_efficiency",
        1100, 750)
    frame_high = max(
        [1] + [max(diagnostics[event]["expected"],
                       counts["hits3d"][event],
                       counts["track_nodes"][event])
               for event in events]) * 1.10
    efficiency_canvas.DrawFrame(
        -0.5, 0.0, len(events)-0.5, frame_high,
        "Selected events: MC voxel and reconstruction counts;"
        "selected-event index;count")
    count_graphs = []
    for key, label, color, marker in (
            ("expected", "MC crossed voxels", ROOT.kGreen + 2, 20),
            ("matched", "3D voxels matched to MC", ROOT.kBlue + 1, 21),
            ("offtrack", "3D voxels off MC", ROOT.kOrange + 7, 22),
            ("nodes", "node-associated hit voxels", ROOT.kRed + 1, 24)):
        count_graph = ROOT.TGraph()
        for index, event in enumerate(events):
            count_graph.SetPoint(index, index, diagnostics[event][key])
        count_graph.SetMarkerColor(color)
        count_graph.SetLineColor(color)
        count_graph.SetMarkerStyle(marker)
        count_graph.Draw("PL SAME")
        count_graphs.append((label, count_graph))
    efficiency_legend = ROOT.TLegend(0.58, 0.68, 0.90, 0.89)
    for label, count_graph in count_graphs:
        efficiency_legend.AddEntry(count_graph, label, "lp")
    efficiency_legend.Draw()
    efficiency_canvas.SaveAs(
        str(args.output_dir / "reconstruction_efficiency.png"))
    efficiency_canvas.Close()
    print("MC voxel comparison written to reconstruction_efficiency.csv")
    if len(events) <= 20:
        for event in events:
            values = diagnostics[event]
            print(f"  Event {event}: MC voxels={values['expected']}, "
                  f"matched 3D={values['matched']}, "
                  f"missing={values['missing']}, "
                  f"off-track={values['offtrack']}, "
                  f"node-associated hit voxels={values['nodes']}")

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
        virtual_segments_by_event = mc_virtual_segments_by_event(
            trees["mc_virtual_segments"], selected_events, axes)

        for event in events:
            fibers = graph(fibers_by_event[event], ROOT.kGray + 1, 7, 0.5)
            hits = graph(hits_by_event[event], ROOT.kBlue + 1, 20, 0.8)
            tracks = [graph(values, ROOT.kRed + 1, 24, 1.0)
                      for values in tracks_by_event.get(event, [])]
            virtual_segments = []
            for start, stop in virtual_segments_by_event[event]:
                segment = ROOT.TLine(start[0], start[1], stop[0], stop[1])
                segment.SetLineColor(ROOT.kCyan + 2)
                segment.SetLineWidth(3)
                virtual_segments.append(segment)

            canvas_name = f"c_{event}_{axes[0]}{axes[1]}"
            canvas = ROOT.TCanvas(canvas_name, canvas_name, 1000, 750)
            title = (
                f"Event {event}: fibres shown={len(fibers_by_event[event])}/"
                f"{counts['fiber_hits'][event]}, "
                f"3D hits={counts['hits3d'][event]}, "
                f"MC cube segments={counts['mc_virtual_segments'][event]}, "
                f"tracks={track_counts[event]}, "
                f"track nodes={counts['track_nodes'][event]};"
                f"{axes[0]} [mm];{axes[1]} [mm]"
            )
            canvas.DrawFrame(-3000, -1500, 3000, 1500, title)
            boundaries = draw_detector_boundaries(axes, detector_name)
            for segment in virtual_segments:
                segment.Draw("SAME")
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
            if virtual_segments:
                legend.AddEntry(
                    virtual_segments[0], "MC virtual-cube segments", "l")
            for name, boundary in boundaries:
                legend.AddEntry(boundary, f"{name} envelope", "l")
            legend.Draw()
            output = args.output_dir / f"event{event}_{axes[0]}{axes[1]}.png"
            canvas.SaveAs(str(output))
            canvas.Close()

            # Detector-focused voxel classification.  This deliberately uses
            # only exact virtual-cube truth, not sparse MC trajectory points.
            voxel_truth = graph(voxel_projection_points(
                diagnostics[event]["expected_voxels"], detector_bounds,
                axes), ROOT.kGray + 2, 25, 0.55)
            voxel_matched = graph(voxel_projection_points(
                diagnostics[event]["matched_voxels"], detector_bounds,
                axes), ROOT.kGreen + 2, 20, 0.75)
            voxel_missing = graph(voxel_projection_points(
                diagnostics[event]["missing_voxels"], detector_bounds,
                axes), ROOT.kOrange + 7, 24, 0.75)
            voxel_offtrack = graph(voxel_projection_points(
                diagnostics[event]["offtrack_voxels"], detector_bounds,
                axes), ROOT.kRed + 1, 5, 0.8)

            focus_name = f"c_{event}_{axes[0]}{axes[1]}_voxel_match"
            focus = ROOT.TCanvas(focus_name, focus_name, 1000, 800)
            margin = 20.0
            focus.DrawFrame(
                detector_bounds[axes[0]][0] - margin,
                detector_bounds[axes[1]][0] - margin,
                detector_bounds[axes[0]][1] + margin,
                detector_bounds[axes[1]][1] + margin,
                f"Event {event}: LFGD virtual-voxel matching;"
                f"{axes[0]} [mm];{axes[1]} [mm]")
            focus.SetGrid()
            lfgd_box = ROOT.TBox(
                detector_bounds[axes[0]][0],
                detector_bounds[axes[1]][0],
                detector_bounds[axes[0]][1],
                detector_bounds[axes[1]][1])
            lfgd_box.SetFillStyle(0)
            lfgd_box.SetLineColor(ROOT.kBlack)
            lfgd_box.SetLineWidth(2)
            lfgd_box.Draw("L SAME")
            for voxel_graph in (
                    voxel_truth, voxel_missing, voxel_matched,
                    voxel_offtrack):
                if voxel_graph.GetN():
                    voxel_graph.Draw("P SAME")
            focus_legend = ROOT.TLegend(0.60, 0.72, 0.91, 0.90)
            if voxel_truth.GetN():
                focus_legend.AddEntry(
                    voxel_truth, "all MC truth voxels", "p")
            if voxel_matched.GetN():
                focus_legend.AddEntry(
                    voxel_matched, "matched reconstructed voxels", "p")
            if voxel_missing.GetN():
                focus_legend.AddEntry(
                    voxel_missing, "MC voxels without 3D hit", "p")
            if voxel_offtrack.GetN():
                focus_legend.AddEntry(
                    voxel_offtrack, "unmatched reconstructed voxels", "p")
            focus_legend.Draw()
            focus.SaveAs(str(
                args.output_dir
                / f"event{event}_{axes[0]}{axes[1]}_voxel_match.png"))
            focus.Close()

            all_reco = graph(voxel_projection_points(
                diagnostics[event]["reco_voxels"], detector_bounds, axes),
                ROOT.kGray + 2, 25, 0.55)
            node_reco = graph(voxel_projection_points(
                diagnostics[event]["node_matched_reco_voxels"],
                detector_bounds, axes), ROOT.kGreen + 2, 20, 0.75)
            reco_without_node = graph(voxel_projection_points(
                diagnostics[event]["reco_without_node_voxels"],
                detector_bounds, axes), ROOT.kRed + 1, 5, 0.8)
            node_without_reco = graph(voxel_projection_points(
                diagnostics[event]["nodes_without_reco_voxels"],
                detector_bounds, axes), ROOT.kMagenta + 1, 34, 0.75)
            track_focus_name = f"c_{event}_{axes[0]}{axes[1]}_track_match"
            track_focus = ROOT.TCanvas(
                track_focus_name, track_focus_name, 1000, 800)
            track_focus.DrawFrame(
                detector_bounds[axes[0]][0] - margin,
                detector_bounds[axes[1]][0] - margin,
                detector_bounds[axes[0]][1] + margin,
                detector_bounds[axes[1]][1] + margin,
                f"Event {event}: LFGD 3D-hit/node-measurement association;"
                f"{axes[0]} [mm];{axes[1]} [mm]")
            track_focus.SetGrid()
            track_box = ROOT.TBox(
                detector_bounds[axes[0]][0],
                detector_bounds[axes[1]][0],
                detector_bounds[axes[0]][1],
                detector_bounds[axes[1]][1])
            track_box.SetFillStyle(0)
            track_box.SetLineColor(ROOT.kBlack)
            track_box.SetLineWidth(2)
            track_box.Draw("L SAME")
            for voxel_graph in (
                    all_reco, node_reco, reco_without_node,
                    node_without_reco):
                if voxel_graph.GetN():
                    voxel_graph.Draw("P SAME")
            track_legend = ROOT.TLegend(0.58, 0.72, 0.91, 0.90)
            if all_reco.GetN():
                track_legend.AddEntry(all_reco, "all 3D-hit voxels", "p")
            if node_reco.GetN():
                track_legend.AddEntry(
                    node_reco, "3D voxels represented by track nodes", "p")
            if reco_without_node.GetN():
                track_legend.AddEntry(
                    reco_without_node, "3D voxels not used by tracks", "p")
            if node_without_reco.GetN():
                track_legend.AddEntry(
                    node_without_reco, "node-hit voxels absent from hits3d", "p")
            track_legend.Draw()
            track_focus.SaveAs(str(
                args.output_dir
                / f"event{event}_{axes[0]}{axes[1]}_track_match.png"))
            track_focus.Close()

    write_index(args.output_dir, events, counts, track_counts, charge_sums,
                diagnostics)
    print(f"Browse {args.output_dir / 'index.html'}")


if __name__ == "__main__":
    main()
