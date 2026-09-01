#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/config.sh"

usage() {
    cat <<EOF
Usage: $(basename "$0") ENERGY_MEV EVENTS OUTPUT_DIR

Generate one reusable HOMO Geant4 sample and stop after g4.root. Direction,
position, particle, and seed can be overridden through DIRECTION_MODE,
DIRECTION, CONE_HALF_ANGLE_DEG, POSITION_MM, POSITION_FRAME, PARTICLE, and SEED.

Example:
  DIRECTION_MODE=cone DIRECTION="1 0 0" CONE_HALF_ANGLE_DEG=5 $(basename "$0") 5000 1000 Studies/LightmapsStudie/5GeV/common
EOF
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then usage; exit 0; fi
(( $# == 3 )) || { usage >&2; exit 2; }
energy_mev=$1
events=$2
output_dir=$(realpath -m "$3")
[[ $energy_mev =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "ENERGY_MEV must be positive numeric text" >&2; exit 2; }
[[ $events =~ ^[1-9][0-9]*$ ]] || { echo "EVENTS must be a positive integer" >&2; exit 2; }
[[ ! -e $output_dir ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 2; }
for program in python3 ND280GEANT4SIM.exe; do command -v "$program" >/dev/null || { echo "Missing program: $program" >&2; exit 1; }; done

mkdir -p "$output_dir"
primary_csv="$output_dir/primary_events.csv"
macro="$output_dir/gps.mac"
g4="$output_dir/g4.root"
read -r px py pz <<< "$POSITION_MM"
read -r dx dy dz <<< "$DIRECTION"

python3 "$script_dir/generate_primary_events.py" --events "$events" --seed "$SEED" --particle "$PARTICLE" --energy-mev "$energy_mev" --position-mm "$px" "$py" "$pz" --position-frame "$POSITION_FRAME" --direction-mode "$DIRECTION_MODE" --direction "$dx" "$dy" "$dz" --cone-half-angle-deg "$CONE_HALF_ANGLE_DEG" --output "$primary_csv"
python3 "$script_dir/generate_gps_macro.py" --detector homo --baseline "$BASELINE" --events "$events" --particle "$PARTICLE" --energy-mev "$energy_mev" --position-mm "$px" "$py" "$pz" --position-frame "$POSITION_FRAME" --direction-mode "$DIRECTION_MODE" --direction "$dx" "$dy" "$dz" --cone-half-angle-deg "$CONE_HALF_ANGLE_DEG" --primary-input "$primary_csv" --output "$macro"

echo "Generating $events identical-replay primaries: $PARTICLE at $energy_mev MeV"
echo "Direction: $DIRECTION_MODE around $DIRECTION; cone half-angle $CONE_HALF_ANGLE_DEG deg"
echo "Seed: $SEED"
ND280GEANT4SIM.exe -s "$SEED" -o "${g4%.root}" "$macro" 2>&1 | tee "$output_dir/01_geant4.log"
[[ -s $g4 ]] || { echo "Geant4 did not produce $g4" >&2; exit 4; }
sha256sum "$g4" > "$output_dir/g4.sha256"
echo "Common Geant4 sample: $g4"
