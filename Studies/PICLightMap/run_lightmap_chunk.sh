#!/usr/bin/env bash

# Run a contiguous chunk of light-map positions. HTCondor starts this inside the
# unmodified ND280 SIF via +SingularityImage.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_lightmap_chunk.sh START_INDEX POSITION_COUNT NPHOTONS SCAT_VALUE SCAT_UNIT ABS_VALUE ABS_UNIT SEED_BASE FIBRE_DIAMETER_MM SCINTILLATOR_RI SCATTERING_MODEL MIE_G MIE_FORWARD_FRACTION [POSITIONS_FILE]

Example (positions 0..99, 100k photons each, 0.30 mm Mie scattering, 5 m absorption):
  ./run_lightmap_chunk.sh 0 100 100000 0.30 mm 5 m 730000 1.0 1.48 double_hg 0.5 0.75 positions.txt

START_INDEX is zero based. POSITION_COUNT must be between 50 and 100, and
the requested range must fit inside positions.txt.
EOF
}

if [[ ${1:-} == "--help" || ${1:-} == "-h" || $# -lt 13 || $# -gt 14 ]]; then
    usage
    [[ ${1:-} == "--help" || ${1:-} == "-h" ]] && exit 0
    exit 2
fi

[[ ${LIGHTMAP_DETECTOR:-} == "HOMO" ]] || {
    echo "Refusing to run: LIGHTMAP_DETECTOR=HOMO was not supplied by the submission file" >&2
    exit 2
}

start_index=$1
chunk_count=$2
nphotons=$3
scat_value=$4
scat_unit=$5
abs_value=$6
abs_unit=$7
seed_base=$8
fibre_diameter_mm=$9
scintillator_ri=${10}
scattering_model=${11}
mie_g=${12}
mie_forward_fraction=${13}
positions_file=${14:-positions.txt}

[[ $start_index =~ ^[0-9]+$ ]] || { echo "START_INDEX must be a non-negative integer" >&2; exit 2; }
[[ $chunk_count =~ ^[1-9][0-9]*$ ]] || { echo "POSITION_COUNT must be a positive integer" >&2; exit 2; }
(( chunk_count >= 50 && chunk_count <= 100 )) || { echo "POSITION_COUNT must be between 50 and 100" >&2; exit 2; }
[[ $nphotons =~ ^[1-9][0-9]*$ ]] || { echo "NPHOTONS must be a positive integer" >&2; exit 2; }
[[ $seed_base =~ ^[1-9][0-9]*$ ]] || { echo "SEED_BASE must be a positive integer" >&2; exit 2; }
[[ $scat_value =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "SCAT_VALUE must be non-negative numeric text" >&2; exit 2; }
[[ $abs_value =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "ABS_VALUE must be non-negative numeric text" >&2; exit 2; }
[[ $fibre_diameter_mm =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "FIBRE_DIAMETER_MM must be positive numeric text" >&2; exit 2; }
[[ $scintillator_ri =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "SCINTILLATOR_RI must be positive numeric text" >&2; exit 2; }
[[ $scattering_model =~ ^(legacy_rayleigh|double_hg)$ ]] || { echo "SCATTERING_MODEL must be legacy_rayleigh or double_hg" >&2; exit 2; }
[[ $mie_g =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "MIE_G must be numeric text" >&2; exit 2; }
[[ $mie_forward_fraction =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "MIE_FORWARD_FRACTION must be numeric text" >&2; exit 2; }
awk -v d="$fibre_diameter_mm" 'BEGIN { exit !(d > 0 && d < 10) }' || { echo "FIBRE_DIAMETER_MM must be greater than 0 and less than the 10 mm pitch" >&2; exit 2; }
awk -v n="$scintillator_ri" 'BEGIN { exit !(n > 1) }' || { echo "SCINTILLATOR_RI must be greater than 1" >&2; exit 2; }
awk -v g="$mie_g" 'BEGIN { exit !(g >= 0 && g < 1) }' || { echo "MIE_G must satisfy 0 <= g < 1" >&2; exit 2; }
awk -v r="$mie_forward_fraction" 'BEGIN { exit !(r >= 0 && r <= 1) }' || { echo "MIE_FORWARD_FRACTION must satisfy 0 <= r_fb <= 1" >&2; exit 2; }
[[ $scat_unit =~ ^(nm|um|mm|cm|m)$ ]] || { echo "Unsupported SCAT_UNIT: $scat_unit" >&2; exit 2; }
[[ $abs_unit =~ ^(nm|um|mm|cm|m)$ ]] || { echo "Unsupported ABS_UNIT: $abs_unit" >&2; exit 2; }
[[ -r $positions_file ]] || { echo "Cannot read positions file: $positions_file" >&2; exit 2; }

# Separate simultaneous optical productions before separating their chunks.
# Condor supplies LIGHTMAP_PRODUCTION_TAG because it must know the transfer
# path before execution.  Direct/manual runs derive it from SCAT_VALUE.
scat_directory_tag="scat_${scat_value}mm"
fibre_directory_tag="fibre_${fibre_diameter_mm}mm"
ri_directory_tag="scin_ri_${scintillator_ri}"
model_directory_tag="model_${scattering_model}"
mie_directory_tag="g_${mie_g}_rfb_${mie_forward_fraction}"
derived_production_tag="${fibre_directory_tag}/${model_directory_tag}/${scat_directory_tag}/${mie_directory_tag}/${ri_directory_tag}"
production_tag=${LIGHTMAP_PRODUCTION_TAG:-$derived_production_tag}
[[ $production_tag == "$derived_production_tag" ]] || {
    echo "LIGHTMAP_PRODUCTION_TAG=$production_tag does not match the requested fibre/scattering/scintillator-index settings (expected $derived_production_tag)" >&2
    exit 2
}
job_directory="${production_tag}/chunk_${start_index}"
mkdir -p "$job_directory/results" "$job_directory/macros"

position_count=$(wc -l < "$positions_file")
end_index=$((start_index + chunk_count - 1))
(( end_index < position_count )) || {
    echo "Requested range $start_index..$end_index is outside 0..$((position_count - 1))" >&2
    exit 2
}

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

start_tag=$(printf '%04d' "$start_index")
end_tag=$(printf '%04d' "$end_index")
scat_tag=${scat_value//./p}${scat_unit}
abs_tag=${abs_value//./p}${abs_unit}
fibre_tag=${fibre_diameter_mm//./p}mm
ri_tag=${scintillator_ri//./p}
g_tag=${mie_g//./p}
rfb_tag=${mie_forward_fraction//./p}
stem="lightmap_fibre-${fibre_tag}_${scattering_model}_scat-${scat_tag}_g-${g_tag}_rfb-${rfb_tag}_abs-${abs_tag}_scin-ri-${ri_tag}_pos-${start_tag}-${end_tag}"
macro="${job_directory}/macros/${stem}.mac"
output="${job_directory}/results/${stem}"
seed=$((seed_base + start_index))
fibre_radius_mm=$(awk -v d="$fibre_diameter_mm" 'BEGIN { printf "%.12g", d/2 }')
fibre_gap_mm=$(awk -v d="$fibre_diameter_mm" 'BEGIN { printf "%.12g", 10-d }')

cat > "$macro" <<EOF
/t2k/control baseline-2024-plusplus 1.0
# The baseline disables HOMO and enables HFG by default.  A light-map job must
# select HOMO explicitly before /t2k/update constructs the geometry.
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoHFGCMD true
/t2k/OA/Magnet/Basket/PlusPlusTracker/NoHomoCMD false
/t2k/OA/Magnet/Basket/PlusPlusTracker/HOMO/Box/gapX ${fibre_gap_mm} mm
/t2k/OA/Magnet/Basket/PlusPlusTracker/HOMO/Box/gapY ${fibre_gap_mm} mm
/t2k/OA/Magnet/Basket/PlusPlusTracker/HOMO/Box/gapZ ${fibre_gap_mm} mm
/t2k/OA/Magnet/Basket/PlusPlusTracker/HOMO/Box/HomoFiber/radius ${fibre_radius_mm} mm
/t2k/phys/fullPhotonScintSim
EOF

# Only the new SIF understands the double-HG commands. The legacy branch is
# intentionally compatible with the historical ND280Reco_03.sif.
if [[ $scattering_model == double_hg ]]; then
    cat >> "$macro" <<EOF
/t2k/detector/liquidOScatteringModel double_hg
/t2k/detector/liquidOMieG ${mie_g}
/t2k/detector/liquidOMieForwardFraction ${mie_forward_fraction}
/t2k/detector/liquidOMieScatLen ${scat_value} ${scat_unit}
EOF
else
    cat >> "$macro" <<EOF
/t2k/detector/liquidOscatLen ${scat_value} ${scat_unit}
EOF
fi

cat >> "$macro" <<EOF
/t2k/detector/liquidOabsLen ${abs_value} ${abs_unit}
/t2k/detector/liquidORI ${scintillator_ri}
/t2k/random/randomSeed ${seed}
/t2k/random/showRandomSeed
/t2k/update
/gps/source/clear
/gps/source/multiplevertex true
EOF

for ((position_index=start_index; position_index<=end_index; ++position_index)); do
    position=$(sed -n "$((position_index + 1))p" "$positions_file")
    [[ -n $position ]] || { echo "Position $position_index is empty" >&2; exit 2; }
    cat >> "$macro" <<EOF
/gps/source/clear
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
EOF
done

cat >> "$macro" <<EOF
/generator/add
EOF

echo "Host:              $(hostname)"
echo "ND280 executable:  $(command -v ND280GEANT4SIM.exe)"
echo "Position range:    $start_index..$end_index of $position_count"
echo "Production dir:    $production_tag"
echo "Positions in job:  $chunk_count"
echo "Photons/position:  $nphotons"
echo "Scattering model:  $scattering_model"
echo "Scattering MFP:    $scat_value $scat_unit"
echo "Mie g / r_fb:      $mie_g / $mie_forward_fraction (r_fb is F/(F+B))"
echo "Absorption length: $abs_value $abs_unit"
echo "Fibre diameter:    $fibre_diameter_mm mm"
echo "Scintillator RI:   $scintillator_ri"
echo "Fibre radius/gap:  $fibre_radius_mm/$fibre_gap_mm mm (10 mm pitch)"
echo "Random seed:       $seed"
echo "Output:            ${output}.root"

simulation_log="${job_directory}/simulation.log"
timing_csv="${job_directory}/position_timing.csv"
set +e
ND280GEANT4SIM.exe -o "$output" "$macro" 2>&1 | tee "$simulation_log"
simulation_status=${PIPESTATUS[0]}
set -e

# Geant4 prints one "% User=... Real=... Sys=..." record after each
# /run/beamOn. Match it to the preceding GPS position and convert the local
# run number into the global zero-based light-map position index.
awk -v start="$start_index" '
BEGIN {
    print "position_index,local_run,x_mm,y_mm,z_mm,user_s,real_s,sys_s"
    local_run = -1
}
/^\/gps\/position[[:space:]]/ {
    ++local_run
    x = $2; y = $3; z = $4
    next
}
/User=[0-9.]+s[[:space:]]+Real=[0-9.]+s[[:space:]]+Sys=[0-9.]+s/ && local_run >= 0 {
    user = real = sys = ""
    for (field = 1; field <= NF; ++field) {
        split($field, value, "=")
        gsub(/s$/, "", value[2])
        if (value[1] == "User") user = value[2]
        else if (value[1] == "Real") real = value[2]
        else if (value[1] == "Sys") sys = value[2]
    }
    print start + local_run "," local_run "," x "," y "," z "," user "," real "," sys
}
' "$simulation_log" > "$timing_csv"

timed_positions=$(( $(wc -l < "$timing_csv") - 1 ))
echo "Per-position timings: $timing_csv ($timed_positions records)"

if (( simulation_status != 0 )); then
    echo "ND280GEANT4SIM.exe failed with status $simulation_status; partial timing data retained" >&2
    exit "$simulation_status"
fi

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
