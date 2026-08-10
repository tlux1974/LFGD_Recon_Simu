#!/usr/bin/env python3
"""Generate a fixed-point ND280 Geant4 GPS macro."""

import argparse
import csv
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--detector", choices=("homo", "hfg"), required=True)
    parser.add_argument("--baseline", default="baseline-2024-plusplus")
    parser.add_argument("--events", type=int, required=True)
    parser.add_argument("--particle", default="mu-")
    parser.add_argument("--energy-mev", type=float, default=700.0)
    parser.add_argument("--position-mm", nargs=3, type=float, required=True)
    parser.add_argument("--position-frame", choices=("global", "plusplus"),
                        default="global")
    parser.add_argument("--direction-mode", choices=("fixed", "isotropic", "cone"),
                        default="fixed")
    parser.add_argument("--direction", nargs=3, type=float,
                        default=(0.0, 0.0, 1.0))
    parser.add_argument("--cone-half-angle-deg", type=float, default=5.0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--primary-input", type=Path,
                        help="CSV with explicit per-event primaries")
    args = parser.parse_args()

    # baseline-2024-plusplus has the two alternative detectors at the same
    # location. Select exactly one before /t2k/update constructs the geometry.
    if args.detector == "homo":
        detector_commands = """\
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoHFGCMD true
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoHomoCMD false"""
    else:
        detector_commands = """\
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoHFGCMD false
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoHomoCMD true"""

    # This is a focused HOMO/HFG comparison.  Avoid constructing unrelated
    # upgrade detectors, particularly when the 6.07-million-node HOMO virtual
    # lattice is present.  Keeping them all enabled can exhaust memory before
    # the first event.  TPC3 deliberately remains enabled because the current
    # detector-response initialization expects its drift volumes to exist even
    # when TPC response is disabled.
    focused_geometry_commands = """\
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoInactiveWaterCMD true
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoTPC3CMD false
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoSWDCMD true
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoSciFiCMD true
/t2k/OA/Magnet/Basket/SFG/enable false
/t2k/OA/Magnet/Basket/TopHAT/enable false
/t2k/OA/Magnet/Basket/BottomHAT/enable false
/t2k/OA/Magnet/Basket/UTOF/enable false
/t2k/OA/Magnet/Basket/DTOF/enable false
/t2k/OA/Magnet/Basket/BTOF/enable false
/t2k/OA/Magnet/Basket/TTOF/enable false
/t2k/OA/Magnet/Basket/NTOF/enable false
/t2k/OA/Magnet/Basket/STOF/enable false"""

    input_px, input_py, input_pz = args.position_mm
    if args.position_frame == "plusplus":
        # Placement of PlusPlusTracker in baseline-2024-plusplus, verified
        # against the TGeoManager saved in g4.root.
        px, py, pz = input_px, input_py + 30.0, input_pz - 890.0
    else:
        px, py, pz = input_px, input_py, input_pz
    dx, dy, dz = args.direction
    primary_rows = []
    if args.primary_input:
        with args.primary_input.open(newline="", encoding="utf-8") as stream:
            primary_rows = list(csv.DictReader(stream))
        required = {"event", "particle", "pdg", "kinetic_energy_mev",
                    "x_mm", "y_mm", "z_mm", "dx", "dy", "dz"}
        if not primary_rows or not required.issubset(primary_rows[0]):
            parser.error("primary input has missing columns or no events")
        if len(primary_rows) < args.events:
            parser.error("primary input contains fewer rows than --events")
        primary_rows = primary_rows[:args.events]
        if [int(row["event"]) for row in primary_rows] != list(range(args.events)):
            parser.error("primary input event numbers must be consecutive from zero")
        angular_commands = "# Directions are set explicitly before each event."
    elif args.direction_mode == "isotropic":
        angular_commands = """\
/gps/ang/type iso
/gps/ang/mintheta 0 deg
/gps/ang/maxtheta 180 deg"""
    elif args.direction_mode == "fixed":
        angular_commands = f"/gps/direction {dx:g} {dy:g} {dz:g}"
    else:
        norm = math.sqrt(dx*dx + dy*dy + dz*dz)
        if norm == 0.0:
            parser.error("cone direction must be nonzero")
        ax, ay, az = dx/norm, dy/norm, dz/norm
        if abs(az) < 0.9:
            r1x, r1y, r1z = -ay, ax, 0.0
        else:
            r1x, r1y, r1z = 0.0, -az, ay
        r1norm = math.sqrt(r1x*r1x+r1y*r1y+r1z*r1z)
        r1x, r1y, r1z = r1x/r1norm, r1y/r1norm, r1z/r1norm
        r2x, r2y, r2z = (ay*r1z-az*r1y, az*r1x-ax*r1z,
                          ax*r1y-ay*r1x)
        angular_commands = f"""\
/gps/ang/type iso
/gps/ang/rot1 {r1x:g} {r1y:g} {r1z:g}
/gps/ang/rot2 {r2x:g} {r2y:g} {r2z:g}
/gps/ang/mintheta 0 deg
/gps/ang/maxtheta {args.cone_half_angle_deg:g} deg"""
    text = f"""\
/t2k/control {args.baseline} 1.0
/t2k/field 0.0 tesla   
{detector_commands}
{focused_geometry_commands}
/t2k/update
/gps/source/clear
/gps/source/multiplevertex true
/gps/source/add 1
/gps/particle {args.particle}
{angular_commands}
/gps/ene/type Mono
/gps/ene/mono {args.energy_mev:g} MeV
# Input vertex ({args.position_frame} frame): {input_px:g} {input_py:g} {input_pz:g} mm
# GPS vertex (global ND280 frame): {px:g} {py:g} {pz:g} mm
/gps/position {px:g} {py:g} {pz:g} mm
/gps/pos/type Point
/generator/add
"""
    if primary_rows:
        text += f"# Explicit primaries read from: {args.primary_input}\n"
        for row in primary_rows:
            text += f"""# Primary event {row['event']}, PDG {row['pdg']}
/gps/particle {row['particle']}
/gps/ene/mono {row['kinetic_energy_mev']} MeV
/gps/position {row['x_mm']} {row['y_mm']} {row['z_mm']} mm
/gps/direction {row['dx']} {row['dy']} {row['dz']}
/run/beamOn 1
"""
    else:
        text += f"/run/beamOn {args.events}\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
