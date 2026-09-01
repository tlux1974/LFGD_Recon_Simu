#!/usr/bin/env bash

set -euo pipefail

study_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
output=${1:-$study_dir/homo_response_octant_mirrored_10bin_5m_0p7mm_100kPhotons.root}
csv=${2:-${output%.root}_efficiency.csv}
mask=${3:-$study_dir/../../../SoftProj/detResponseSim/input/homo_response_250514_10bin_5m_1mm_100kPhotons.root}
map_exe="$study_dir/../../../SoftProj/nd280Geant4Sim/Linux-AlmaLinux_9.5-gcc_11-x86_64/bin/MAKEHOMOPHOTONMAP.exe"

mapfile -t inputs < <(
    find "$study_dir" -mindepth 3 -maxdepth 3 \
        -path "$study_dir/octant_chunk_*/results/*.root" -type f -print | sort -V
)
[[ ${#inputs[@]} == 5 ]] || {
    echo "Expected five octant chunk ROOT files, found ${#inputs[@]}" >&2
    exit 2
}
[[ -s $mask ]] || { echo "Missing valid-position mask: $mask" >&2; exit 2; }
[[ -x $map_exe ]] || { echo "Updated merger is not built: $map_exe" >&2; exit 3; }
command -v root >/dev/null || { echo "ROOT is not configured" >&2; exit 3; }

"$map_exe" -R -G "${inputs[0]}" -o "$output" \
    -O "cubeCentre=0 30 910 mm" \
    -O "cubeSide=10 mm" \
    -O "nBins=10" \
    -O "nPhotons=100000" \
    -O "useMirroring=1" \
    -O "validBinMask=$mask" \
    "${inputs[@]}"

root -l -b -q "$study_dir/ExportLightMapEfficiency.C(\"$output\",\"$csv\")"
entries=$(($(wc -l < "$csv") - 1))
[[ $entries == 970 ]] || {
    echo "Expected 970 valid mirrored positions, found $entries" >&2
    exit 4
}
echo "Mirrored full map: $output"
echo "Efficiency CSV:   $csv (970 positions)"
