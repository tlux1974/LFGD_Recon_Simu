#!/usr/bin/env python3
"""Overlay curves from an arbitrary collection of CSV analysis outputs."""

import argparse
import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/lightmap-matplotlib")
import matplotlib.pyplot as plt


def parse_input(specification):
    if not specification.strip():
        raise ValueError(
            "Received a blank CSV argument. In a multiline Bash command, "
            "the backslash must be the final character on its line (no spaces after it).")
    if "=" in specification:
        label, filename = specification.split("=", 1)
    else:
        filename = specification
        label = Path(filename).parent.name or Path(filename).stem
    label = label.strip()
    filename = filename.strip()
    if not label or not filename:
        raise ValueError(f"Invalid CSV specification: {specification!r}")
    path = Path(filename).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"CSV input does not exist: {path}")
    return label, path


def read_rows(label, path, x_column, y_column, group_column):
    output = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {x_column, y_column}
        if group_column:
            required.add(group_column)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for line, row in enumerate(reader, 2):
            if not row[x_column].strip() or not row[y_column].strip():
                continue
            try:
                x_value = float(row[x_column])
                y_value = float(row[y_column])
            except ValueError as error:
                raise ValueError(f"{path}:{line}: non-numeric curve value") from error
            if not math.isfinite(x_value) or not math.isfinite(y_value):
                continue
            output.append((row[group_column] if group_column else "all",
                           x_value, y_value, label))
    return output


def main():
    parser = argparse.ArgumentParser(
        description="Overlay curves from CSV files. Inputs may be FILE.csv or LABEL=FILE.csv.")
    parser.add_argument("inputs", nargs="+", help="CSV inputs, optionally prefixed by LABEL=")
    parser.add_argument("--x-column", default="mean_pair_distance_mm")
    parser.add_argument("--y-column", default="combined_pair_ratio")
    parser.add_argument("--group-column", default="track_axis",
                        help="Column defining panels; pass an empty value for one panel")
    parser.add_argument("--xlabel")
    parser.add_argument("--ylabel")
    parser.add_argument("--title", default="Light-map curve comparison")
    parser.add_argument("--output", type=Path, default=Path("curve_comparison.png"))
    parser.add_argument("--xscale", choices=("linear", "log"), default="linear")
    parser.add_argument("--yscale", choices=("linear", "log"), default="linear")
    parser.add_argument("--connect", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--dpi", type=int, default=170)
    args = parser.parse_args()
    group_column = args.group_column or None

    all_rows, input_labels = [], []
    for specification in args.inputs:
        label, path = parse_input(specification)
        if label in input_labels:
            raise ValueError(f"Duplicate input label: {label}")
        input_labels.append(label)
        all_rows.extend(read_rows(label, path, args.x_column, args.y_column,
                                  group_column))
    if not all_rows:
        raise RuntimeError("No finite curve points were found")

    groups = sorted({row[0] for row in all_rows})
    figure, axes = plt.subplots(
        1, len(groups), figsize=(max(7, 5.5 * len(groups)), 5.3),
        squeeze=False, constrained_layout=True)
    markers = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "*")
    for group_index, group in enumerate(groups):
        axis = axes[0, group_index]
        for label_index, label in enumerate(input_labels):
            points = sorted((x, y) for row_group, x, y, row_label in all_rows
                            if row_group == group and row_label == label)
            if not points:
                continue
            x_values, y_values = zip(*points)
            axis.plot(x_values, y_values,
                      marker=markers[label_index % len(markers)],
                      linestyle="-" if args.connect else "none",
                      linewidth=1.7, markersize=5, label=label)
        axis.set_xscale(args.xscale)
        axis.set_yscale(args.yscale)
        axis.set_xlabel(args.xlabel or args.x_column.replace("_", " "))
        axis.set_ylabel(args.ylabel or args.y_column.replace("_", " "))
        axis.grid(alpha=.25, which="both")
        if group_column:
            axis.set_title(f"{group_column.replace('_', ' ')} = {group}")
        axis.legend()
    figure.suptitle(args.title)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=args.dpi)
    plt.close(figure)
    print(f"Wrote {args.output} from {len(input_labels)} CSV files and {len(groups)} group(s)")


if __name__ == "__main__":
    main()
