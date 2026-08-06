#!/usr/bin/env python3
"""Create a detector-independent, reproducible list of primary particles."""

import argparse
import csv
import math
import random
from pathlib import Path


PDG = {"mu-": 13, "mu+": -13, "pi+": 211, "pi-": -211,
       "proton": 2212, "e-": 11, "e+": -11, "gamma": 22}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--particle", required=True)
    parser.add_argument("--pdg", type=int)
    parser.add_argument("--energy-mev", type=float, required=True)
    parser.add_argument("--position-mm", nargs=3, type=float, required=True)
    parser.add_argument("--position-frame", choices=("global", "plusplus"),
                        default="plusplus")
    parser.add_argument("--direction-mode", choices=("fixed", "isotropic"),
                        required=True)
    parser.add_argument("--direction", nargs=3, type=float,
                        default=(0.0, 0.0, 1.0))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.events < 1:
        parser.error("--events must be positive")
    pdg = args.pdg if args.pdg is not None else PDG.get(args.particle)
    if pdg is None:
        parser.error("unknown PDG: pass --pdg for this particle name")

    x, y, z = args.position_mm
    if args.position_frame == "plusplus":
        # PlusPlusTracker placement in baseline-2024-plusplus.
        x, y, z = x, y + 30.0, z - 890.0

    rng = random.Random(args.seed)
    rows = []
    for event in range(args.events):
        if args.direction_mode == "isotropic":
            cos_theta = rng.uniform(-1.0, 1.0)
            phi = rng.uniform(0.0, 2.0*math.pi)
            sin_theta = math.sqrt(max(0.0, 1.0-cos_theta*cos_theta))
            dx = sin_theta*math.cos(phi)
            dy = sin_theta*math.sin(phi)
            dz = cos_theta
        else:
            dx, dy, dz = args.direction
            norm = math.sqrt(dx*dx+dy*dy+dz*dz)
            if norm == 0.0:
                parser.error("fixed direction must be nonzero")
            dx, dy, dz = dx/norm, dy/norm, dz/norm
        rows.append((event, args.particle, pdg, args.energy_mev,
                     x, y, z, dx, dy, dz))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("event", "particle", "pdg", "kinetic_energy_mev",
                         "x_mm", "y_mm", "z_mm", "dx", "dy", "dz"))
        writer.writerows(rows)
    print(args.output)


if __name__ == "__main__":
    main()
