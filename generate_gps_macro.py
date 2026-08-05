#!/usr/bin/env python3
"""Generate a fixed-point ND280 Geant4 GPS macro."""

import argparse
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
    parser.add_argument("--direction-mode", choices=("fixed", "isotropic"),
                        default="fixed")
    parser.add_argument("--direction", nargs=3, type=float,
                        default=(0.0, 0.0, 1.0))
    parser.add_argument("--output", type=Path, required=True)
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
    if args.direction_mode == "isotropic":
        angular_commands = """\
/gps/ang/type iso
/gps/ang/mintheta 0 deg
/gps/ang/maxtheta 180 deg"""
    else:
        angular_commands = f"/gps/direction {dx:g} {dy:g} {dz:g}"
    text = f"""\
/t2k/control {args.baseline} 1.0
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
/run/beamOn {args.events}
"""
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
