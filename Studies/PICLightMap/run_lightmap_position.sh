#!/usr/bin/env bash

# Run one light-map grid position.  HTCondor starts this script inside the
# unmodified ND280 SIF via +SingularityImage.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_lightmap_position.sh POSITION_INDEX NPHOTONS SCAT_VALUE SCAT_UNIT ABS_VALUE ABS_UNIT SEED_BASE FIBRE_DIAMETER_MM [POSITIONS_FILE]

Example (first grid point, 100k photons, 0.5 mm scattering, 5 m absorption):
  ./run_lightmap_position.sh 0 100000 0.5 mm 5 m 730000 1.5 positions.txt

POSITION_INDEX is zero based: 0..999 for the supplied positions.txt.
EOF
}

if [[ ${1:-} == "--help" || ${1:-} == "-h" || $# -lt 8 || $# -gt 9 ]]; then
    usage
    [[ ${1:-} == "--help" || ${1:-} == "-h" ]] && exit 0
    exit 2
fi

[[ ${LIGHTMAP_DETECTOR:-} == "HOMO" ]] || {
    echo "Refusing to run: LIGHTMAP_DETECTOR=HOMO was not supplied by the submission file" >&2
    exit 2
}

position_index=$1
nphotons=$2
scat_value=$3
scat_unit=$4
abs_value=$5
abs_unit=$6
seed_base=$7
fibre_diameter_mm=$8
positions_file=${9:-positions.txt}

[[ $position_index =~ ^[0-9]+$ ]] || { echo "POSITION_INDEX must be a non-negative integer" >&2; exit 2; }
[[ $nphotons =~ ^[1-9][0-9]*$ ]] || { echo "NPHOTONS must be a positive integer" >&2; exit 2; }
[[ $seed_base =~ ^[1-9][0-9]*$ ]] || { echo "SEED_BASE must be a positive integer" >&2; exit 2; }
[[ $scat_value =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "SCAT_VALUE must be non-negative numeric text" >&2; exit 2; }
[[ $abs_value =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "ABS_VALUE must be non-negative numeric text" >&2; exit 2; }
[[ $fibre_diameter_mm =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "FIBRE_DIAMETER_MM must be positive numeric text" >&2; exit 2; }
awk -v d="$fibre_diameter_mm" 'BEGIN { exit !(d > 0 && d < 10) }' || { echo "FIBRE_DIAMETER_MM must be greater than 0 and less than the 10 mm pitch" >&2; exit 2; }
[[ $scat_unit =~ ^(nm|um|mm|cm|m)$ ]] || { echo "Unsupported SCAT_UNIT: $scat_unit" >&2; exit 2; }
[[ $abs_unit =~ ^(nm|um|mm|cm|m)$ ]] || { echo "Unsupported ABS_UNIT: $abs_unit" >&2; exit 2; }
[[ -r $positions_file ]] || { echo "Cannot read positions file: $positions_file" >&2; exit 2; }

derived_production_tag="fibre_${fibre_diameter_mm}mm/scat_${scat_value}mm"
production_tag=${LIGHTMAP_PRODUCTION_TAG:-$derived_production_tag}
[[ $production_tag == "$derived_production_tag" ]] || {
    echo "LIGHTMAP_PRODUCTION_TAG=$production_tag does not match SCAT_VALUE=$scat_value (expected $derived_production_tag)" >&2
    exit 2
}
job_directory="${production_tag}/position_${position_index}"

# Always create the declared Condor output, including before the simulation.
mkdir -p "$job_directory/results" "$job_directory/macros"

position_count=$(wc -l < "$positions_file")
(( position_index < position_count )) || {
    echo "POSITION_INDEX=$position_index is outside 0..$((position_count - 1))" >&2
    exit 2
}

position=$(sed -n "$((position_index + 1))p" "$positions_file")
[[ -n $position ]] || { echo "Position $position_index is empty" >&2; exit 2; }

# The SIF may or may not arrange the ND280 environment automatically.
if ! command -v ND280GEANT4SIM.exe >/dev/null 2>&1; then
    pilot_profile=/usr/local/t2k/current/nd280SoftwarePilot/nd280SoftwarePilot.profile
    setup_script=$(find /usr/local/t2k/current -maxdepth 4 -type f \
        -path '*/nd280SoftwareMaster*/bin/setup.sh' -print -quit 2>/dev/null || true)
    [[ -r $pilot_profile ]] || {
        echo "Cannot read the ND280 pilot profile: $pilot_profile" >&2
        exit 3
    }
    [[ -n $setup_script ]] || {
        echo "ND280GEANT4SIM.exe is unavailable and no ND280 master setup.sh was found" >&2
        exit 3
    }
    # Some package setup fragments return nonzero for advisory conditions
    # (notably a missing cached oaEvent ROOT geometry).  Do not let errexit
    # abort the environment setup; validate the executable explicitly below.
    set +e
    set +u
    # The pilot profile installs nd280-system and the platform helpers that
    # are prerequisites of the generated master setup.sh.
    # shellcheck disable=SC1090
    source "$pilot_profile"
    # shellcheck disable=SC1090
    source "$setup_script"
    set -u
    set -e
fi

command -v ND280GEANT4SIM.exe >/dev/null 2>&1 || {
    echo "ND280GEANT4SIM.exe is still unavailable after setup" >&2
    exit 3
}

index_tag=$(printf '%04d' "$position_index")
scat_tag=${scat_value//./p}${scat_unit}
abs_tag=${abs_value//./p}${abs_unit}
fibre_tag=${fibre_diameter_mm//./p}mm
stem="lightmap_fibre-${fibre_tag}_scat-${scat_tag}_abs-${abs_tag}_pos-${index_tag}"
macro="${job_directory}/macros/${stem}.mac"
output="${job_directory}/results/${stem}"
seed=$((seed_base + position_index))
fibre_radius_mm=$(awk -v d="$fibre_diameter_mm" 'BEGIN { printf "%.12g", d/2 }')
fibre_gap_mm=$(awk -v d="$fibre_diameter_mm" 'BEGIN { printf "%.12g", 10-d }')

cat > "$macro" <<EOF
/t2k/control baseline-2024-plusplus 1.0
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoHFGCMD true
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoHomoCMD false
/t2k/OA/Magnet/Basket/PlusPlusTracker/HOMO/Box/gapX ${fibre_gap_mm} mm
/t2k/OA/Magnet/Basket/PlusPlusTracker/HOMO/Box/gapY ${fibre_gap_mm} mm
/t2k/OA/Magnet/Basket/PlusPlusTracker/HOMO/Box/gapZ ${fibre_gap_mm} mm
/t2k/OA/Magnet/Basket/PlusPlusTracker/HOMO/Box/HomoFiber/radius ${fibre_radius_mm} mm
/t2k/phys/fullPhotonScintSim
/t2k/detector/liquidOscatLen ${scat_value} ${scat_unit}
/t2k/detector/liquidOabsLen ${abs_value} ${abs_unit}
/t2k/random/randomSeed ${seed}
/t2k/random/showRandomSeed
/t2k/update
/gps/source/clear
/gps/source/multiplevertex true
/gps/source/add 1
/gps/particle opticalphoton
/gps/ang/type iso
/gps/ang/maxtheta 180 deg
/gps/ang/rot1 -1 0 0
/gps/ang/rot2 0 1 0
/gps/number ${nphotons}
/gps/ene/type Mono
/gps/ene/mono 1 eV
/gps/position ${position}
/run/beamOn 1
/generator/add
EOF

echo "Host:              $(hostname)"
echo "ND280 executable:  $(command -v ND280GEANT4SIM.exe)"
echo "Position index:    $position_index of $position_count"
echo "Production dir:    $production_tag"
echo "Position:          $position"
echo "Photons:           $nphotons"
echo "Scattering length: $scat_value $scat_unit"
echo "Absorption length: $abs_value $abs_unit"
echo "Fibre diameter:    $fibre_diameter_mm mm"
echo "Fibre radius/gap:  $fibre_radius_mm/$fibre_gap_mm mm (10 mm pitch)"
echo "Random seed:       $seed"
echo "Output:            ${output}.root"

ND280GEANT4SIM.exe -o "$output" "$macro"

[[ -s ${output}.root ]] || {
    echo "Simulation returned without a non-empty ${output}.root" >&2
    exit 4
}

[[ -r ValidateHomoGeometry.C ]] || {
    echo "Missing required output validator ValidateHomoGeometry.C" >&2
    exit 5
}
root -l -b -q "ValidateHomoGeometry.C(\"${output}.root\")"

echo "Completed successfully:"
ls -lh "${output}.root"
