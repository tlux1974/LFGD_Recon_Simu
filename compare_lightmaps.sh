#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
workspace_dir=$(cd -- "$script_dir/.." && pwd)
map_dir=${LIGHTMAP_DIR:-$script_dir/Studies/PICLightMap/lightmaps}
global_fit_dir="$script_dir/Studies/GlobalFitTest"

usage() {
    cat <<EOF
Usage: $(basename "$0") COMMON_G4.root OUTPUT_DIR [MAP.root ...]

Runs detector response, standard HFG reconstruction, flat-tree diagnostics,
and GlobalFit for every selected light map. COMMON_G4.root is read-only and is
reused verbatim, ensuring identical Geant4 events for every map.

If MAP.root arguments are omitted, all *.root files directly under
  $map_dir
are used except files whose names contain "octant"; those mirrored maps are
invalid for detector response. Environment overrides:
  LIGHTMAP_DIR, SEED (default 12345), GLOBAL_FIT_OPTIONS,
  DETRESPONSE_PARAMETER_FILE, HFGRECON_PARAMETER_FILE,
  FAST_HOMO_RESPONSE, DISABLE_HOMO_ATTENUATION, RESUME (default 0).

Set RESUME=1 to continue an interrupted study. Existing ROOT stages are reused
only when ROOT can open them without marking them zombie or recovered.
EOF
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then usage; exit 0; fi
(( $# >= 2 )) || { usage >&2; exit 2; }
common_g4=$(readlink -f "$1")
output_root=$(realpath -m "$2")
shift 2
[[ -s $common_g4 ]] || { echo "Missing common Geant4 file: $common_g4" >&2; exit 2; }

if (( $# )); then
    maps=("$@")
else
    mapfile -t maps < <(find "$map_dir" -maxdepth 1 -type f -name '*.root' ! -name '*octant*' -print | sort)
fi
(( ${#maps[@]} )) || { echo "No light maps selected" >&2; exit 2; }
resume=${RESUME:-0}
[[ $resume == 0 || $resume == 1 ]] || { echo "RESUME must be 0 or 1" >&2; exit 2; }

for program in DETRESPONSESIM.exe HFGRECON.exe LFGDFLATTREE.exe python3 root sha256sum; do
    command -v "$program" >/dev/null || { echo "Missing program: $program" >&2; exit 1; }
done
for program in global_light_fit; do
    [[ -x $global_fit_dir/$program ]] || { echo "Missing executable: $global_fit_dir/$program" >&2; exit 1; }
done
if [[ $resume == 0 && -d $output_root ]] && [[ -n $(find "$output_root" -mindepth 1 -print -quit) ]]; then
    echo "Refusing to overwrite non-empty output directory: $output_root" >&2
    exit 2
fi
[[ ! -e $output_root || -d $output_root ]] || {
    echo "Output path exists and is not a directory: $output_root" >&2
    exit 2
}
mkdir -p "$output_root"

manifest="$output_root/manifest.csv"
if [[ $resume == 0 || ! -f $manifest ]]; then
    printf 'tag,lightmap,lightmap_sha256,g4,g4_sha256,status\n' > "$manifest"
fi
g4_sha=$(sha256sum "$common_g4" | awk '{print $1}')
declare -A used_tags=()

valid_root() {
    local path=$1
    [[ -s $path ]] || return 1
    root -l -b -q -e "TFile f(\"$path\"); if(f.IsZombie() || f.TestBit(TFile::kRecovered)) gSystem->Exit(9);" >/dev/null 2>&1
}

slugify() {
    local value=${1%.root}
    value=${value#homo_response_}
    value=${value//[^A-Za-z0-9._-]/_}
    printf '%s' "$value"
}

for requested_map in "${maps[@]}"; do
    lightmap=$(readlink -f "$requested_map")
    [[ -s $lightmap ]] || { echo "Missing light map: $requested_map" >&2; exit 2; }
    tag=$(slugify "$(basename "$lightmap")")
    [[ -z ${used_tags[$tag]:-} ]] || { echo "Duplicate output tag: $tag" >&2; exit 2; }
    used_tags[$tag]=1
    branch="$output_root/$tag"
    mkdir -p "$branch"

    map_check=$(root -l -b -q -e "TFile f(\"$lightmap\"); auto h=(TH3*)f.Get(\"Starting Vertex Distribution\"); auto p=(TDirectory*)f.Get(\"hitPositions\"); auto q=(TDirectory*)f.Get(\"lightFractions\"); if(f.IsZombie()||!h||!p||!q||p->GetListOfKeys()->GetSize()!=970||q->GetListOfKeys()->GetSize()!=970) gSystem->Exit(9);" 2>&1) || {
        printf '%s\n' "$map_check" >&2
        echo "Invalid 970-bin light map: $lightmap" >&2
        exit 2
    }
    lightmap_sha=$(sha256sum "$lightmap" | awk '{print $1}')
    if [[ $resume == 1 ]] && valid_root "$branch/global_fit.root" \
            && [[ $(<"$branch/lightmap_path.txt") == "$lightmap" ]] \
            && [[ $(<"$branch/common_g4_path.txt") == "$common_g4" ]]; then
        sed -i "\|^$tag,|d" "$manifest"
        printf '%s,%s,%s,%s,%s,complete\n' "$tag" "$lightmap" "$lightmap_sha" "$common_g4" "$g4_sha" >> "$manifest"
        echo "[$tag] Complete; skipping"
        continue
    fi
    sed -i "\|^$tag,|d" "$manifest"
    printf '%s,%s,%s,%s,%s,running\n' "$tag" "$lightmap" "$lightmap_sha" "$common_g4" "$g4_sha" >> "$manifest"
    printf '%s\n' "$lightmap" > "$branch/lightmap_path.txt"
    printf '%s\n' "$common_g4" > "$branch/common_g4_path.txt"

    parameter_file="$branch/detresponse.parameters.dat"
    printf '< detResponseSim.LiquidO.Response.File = %s >\n' "$lightmap" > "$parameter_file"
    if [[ -n ${DETRESPONSE_PARAMETER_FILE:-} ]]; then cat "$DETRESPONSE_PARAMETER_FILE" >> "$parameter_file"; fi
    if [[ ${DISABLE_HOMO_ATTENUATION:-0} == 1 ]]; then cat "$script_dir/detresponse_no_homo_attenuation.parameters.dat" >> "$parameter_file"; fi
    if [[ ${FAST_HOMO_RESPONSE:-0} == 1 ]]; then printf '< detResponseSim.Homo.FastFiberResponse = 1 >\n' >> "$parameter_file"; fi

    if [[ $resume == 1 ]] && valid_root "$branch/detresponse.root"; then
        echo "[$tag] Reusing detector response"
    else
        echo "[$tag] Detector response"
        DETRESPONSESIM.exe -R -G "$common_g4" -O "randomseed=${SEED:-12345}" -O upgrade-nd280plus-homo -O disable-sfg -O disable-swd -O disable-scifi -O disable-ecal -O disable-smrd -O disable-tof -O disable-tpc -O disable-hat -O "par_override=$parameter_file" -o "$branch/detresponse.root" "$common_g4" 2>&1 | tee "$branch/02_detresponse.log"
    fi

    if [[ $resume == 1 ]] && valid_root "$branch/reco.root"; then
        echo "[$tag] Reusing standard reconstruction"
    else
        echo "[$tag] Standard reconstruction"
        reco_args=(-R -o "$branch/reco.root")
        if [[ -n ${HFGRECON_PARAMETER_FILE:-} ]]; then reco_args+=(-O "par_override=$HFGRECON_PARAMETER_FILE"); elif [[ ${DISABLE_HOMO_ATTENUATION:-0} == 1 ]]; then reco_args+=(-O "par_override=$script_dir/hfgrecon_no_homo_attenuation.parameters.dat"); fi
        HFGRECON.exe "${reco_args[@]}" "$branch/detresponse.root" 2>&1 | tee "$branch/03_hfgrecon.log"
    fi
    if [[ $resume == 1 ]] && valid_root "$branch/flat.root"; then
        echo "[$tag] Reusing flat tree"
    else
        LFGDFLATTREE.exe -R -O outfile="$branch/flat.root" "$branch/reco.root" 2>&1 | tee "$branch/04_flat_tree.log"
    fi

    echo "[$tag] GlobalFit"
    echo "GlobalFit light map: $lightmap" | tee "$branch/05_global_fit.log"
    echo "GlobalFit light-map SHA-256: $lightmap_sha" | tee -a "$branch/05_global_fit.log"
    read -r -a fit_options <<< "${GLOBAL_FIT_OPTIONS:-EVENT=all TREE=fiber_hits MIN_CHARGE=10 DBSCAN=1 DBSCAN_EPSILON_MM=14.2 DBSCAN_MIN_POINTS=2}"
    "$global_fit_dir/global_light_fit" "$branch/flat.root" "$branch/global_fit.root" "$lightmap" "${fit_options[@]}" 2>&1 | tee -a "$branch/05_global_fit.log"
    sed -i "s|^$tag,\(.*\),running$|$tag,\1,complete|" "$manifest"
done

echo "ROOT production complete: $output_root"
echo "Manifest: $manifest"
echo "Run analyse_lightmap_study.sh outside the production container to create diagnostics and plots."
