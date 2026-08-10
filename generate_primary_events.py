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
    parser.add_argument("--direction-mode", choices=("fixed", "isotropic", "cone"),
                        required=True)
    parser.add_argument("--direction", nargs=3, type=float,
                        default=(0.0, 0.0, 1.0))
    parser.add_argument("--cone-half-angle-deg", type=float, default=5.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.events < 1:
        parser.error("--events must be positive")
    pdg = args.pdg if args.pdg is not None else PDG.get(args.particle)
    if pdg is None:
        parser.error("unknown PDG: pass --pdg for this particle name")
    if not 0.0 <= args.cone_half_angle_deg <= 180.0:
        parser.error("--cone-half-angle-deg must be between 0 and 180")

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
        elif args.direction_mode == "fixed":
            dx, dy, dz = args.direction
            norm = math.sqrt(dx*dx+dy*dy+dz*dz)
            if norm == 0.0:
                parser.error("fixed direction must be nonzero")
            dx, dy, dz = dx/norm, dy/norm, dz/norm
        else:
            # Uniform solid angle inside a cone around --direction.  Construct
            # a stable orthonormal basis transverse to the cone axis.
            ax, ay, az = args.direction
            norm = math.sqrt(ax*ax+ay*ay+az*az)
            if norm == 0.0:
                parser.error("cone direction must be nonzero")
            ax, ay, az = ax/norm, ay/norm, az/norm
            if abs(az) < 0.9:
                tx, ty, tz = -ay, ax, 0.0
            else:
                tx, ty, tz = 0.0, -az, ay
            norm_t = math.sqrt(tx*tx+ty*ty+tz*tz)
            tx, ty, tz = tx/norm_t, ty/norm_t, tz/norm_t
            ux, uy, uz = (ay*tz-az*ty, az*tx-ax*tz, ax*ty-ay*tx)
            cos_alpha = rng.uniform(
                math.cos(math.radians(args.cone_half_angle_deg)), 1.0)
            sin_alpha = math.sqrt(max(0.0, 1.0-cos_alpha*cos_alpha))
            phi = rng.uniform(0.0, 2.0*math.pi)
            dx = cos_alpha*ax + sin_alpha*(math.cos(phi)*tx + math.sin(phi)*ux)
            dy = cos_alpha*ay + sin_alpha*(math.cos(phi)*ty + math.sin(phi)*uy)
            dz = cos_alpha*az + sin_alpha*(math.cos(phi)*tz + math.sin(phi)*uz)
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
