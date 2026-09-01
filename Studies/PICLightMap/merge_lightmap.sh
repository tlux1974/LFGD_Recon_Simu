#!/usr/bin/env bash

set -euo pipefail

study_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
map_exe="$study_dir/../../../SoftProj/nd280Geant4Sim/Linux-AlmaLinux_9.5-gcc_11-x86_64/bin/MAKEHOMOPHOTONMAP.exe"

usage() {
    cat <<EOF
Usage: $(basename "$0") --scat-length VALUE --fibre-diameter VALUE --scintillator-ri VALUE --scattering-model MODEL --mie-g VALUE --mie-forward-fraction VALUE [--output MAP.root] [--csv EFFICIENCY.csv]

VALUE is the required scattering length in mm, for example 0.7. The merger
selects only chunk ROOT files whose generated filename contains the matching
scattering-length, fibre-diameter, and scintillator-index tags. None of these
physical settings has a default.
EOF
}

scat_length=
fibre_diameter=
scintillator_ri=
scattering_model=
mie_g=
mie_forward_fraction=
output=
csv=
while (( $# )); do
    case $1 in
        --scat-length)
            (( $# >= 2 )) || { echo "--scat-length requires a value" >&2; exit 2; }
            scat_length=$2
            shift 2
            ;;
        --fibre-diameter|--fiber-diameter)
            (( $# >= 2 )) || { echo "$1 requires a value" >&2; exit 2; }
            fibre_diameter=$2
            shift 2
            ;;
        --scintillator-ri)
            (( $# >= 2 )) || { echo "--scintillator-ri requires a value" >&2; exit 2; }
            scintillator_ri=$2
            shift 2
            ;;
        --scattering-model)
            (( $# >= 2 )) || { echo "--scattering-model requires a value" >&2; exit 2; }
            scattering_model=$2
            shift 2
            ;;
        --mie-g)
            (( $# >= 2 )) || { echo "--mie-g requires a value" >&2; exit 2; }
            mie_g=$2
            shift 2
            ;;
        --mie-forward-fraction)
            (( $# >= 2 )) || { echo "--mie-forward-fraction requires a value" >&2; exit 2; }
            mie_forward_fraction=$2
            shift 2
            ;;
        --output)
            (( $# >= 2 )) || { echo "--output requires a path" >&2; exit 2; }
            output=$2
            shift 2
            ;;
        --csv)
            (( $# >= 2 )) || { echo "--csv requires a path" >&2; exit 2; }
            csv=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -n $scat_length ]] || { echo "Missing mandatory --scat-length VALUE" >&2; usage >&2; exit 2; }
[[ -n $fibre_diameter ]] || { echo "Missing mandatory --fibre-diameter VALUE" >&2; usage >&2; exit 2; }
[[ -n $scintillator_ri ]] || { echo "Missing mandatory --scintillator-ri VALUE" >&2; usage >&2; exit 2; }
[[ -n $scattering_model ]] || { echo "Missing mandatory --scattering-model MODEL" >&2; usage >&2; exit 2; }
[[ -n $mie_g ]] || { echo "Missing mandatory --mie-g VALUE" >&2; usage >&2; exit 2; }
[[ -n $mie_forward_fraction ]] || { echo "Missing mandatory --mie-forward-fraction VALUE" >&2; usage >&2; exit 2; }
[[ $scat_length =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "Invalid scattering length: $scat_length" >&2; exit 2; }
[[ $fibre_diameter =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "Invalid fibre diameter: $fibre_diameter" >&2; exit 2; }
[[ $scintillator_ri =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "Invalid scintillator refractive index: $scintillator_ri" >&2; exit 2; }
[[ $scattering_model =~ ^(legacy_rayleigh|double_hg)$ ]] || { echo "Invalid scattering model: $scattering_model" >&2; exit 2; }
[[ $mie_g =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "Invalid Mie g: $mie_g" >&2; exit 2; }
[[ $mie_forward_fraction =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "Invalid Mie forward fraction: $mie_forward_fraction" >&2; exit 2; }
awk -v g="$mie_g" 'BEGIN { exit !(g >= 0 && g < 1) }' || { echo "Mie g must satisfy 0 <= g < 1" >&2; exit 2; }
awk -v r="$mie_forward_fraction" 'BEGIN { exit !(r >= 0 && r <= 1) }' || { echo "Mie forward fraction must satisfy 0 <= r_fb <= 1" >&2; exit 2; }
scat_tag=${scat_length//./p}mm
fibre_tag=${fibre_diameter//./p}mm
ri_tag=${scintillator_ri//./p}
g_tag=${mie_g//./p}
rfb_tag=${mie_forward_fraction//./p}
production_dir="$study_dir/fibre_${fibre_diameter}mm/model_${scattering_model}/scat_${scat_length}mm/g_${mie_g}_rfb_${mie_forward_fraction}/scin_ri_${scintillator_ri}"
output=${output:-$study_dir/homo_response_10bin_5m_${scattering_model}_${scat_tag}_g-${g_tag}_rfb-${rfb_tag}_fibre-${fibre_tag}_scin-ri-${ri_tag}_100kPhotons.root}
csv=${csv:-${output%.root}_efficiency.csv}

# Existing 1 mm-fibre productions predate the diameter-aware directory layout.
# Keep them mergeable, while requiring explicit --fibre-diameter 1.
input_pattern="lightmap_fibre-${fibre_tag}_${scattering_model}_scat-${scat_tag}_g-${g_tag}_rfb-${rfb_tag}_*.root"
if [[ ! -d $production_dir && $scattering_model == legacy_rayleigh ]] && awk -v d="$fibre_diameter" -v n="$scintillator_ri" 'BEGIN { exit !(d == 1 && n == 1.48) }'; then
    production_dir="$study_dir/${scat_length//./p}"
    input_pattern="lightmap_scat-${scat_tag}_*.root"
    echo "Using legacy 1 mm-fibre chunk layout: $production_dir"
fi

mapfile -t inputs < <(
    find "$production_dir" -mindepth 3 -maxdepth 3 \
        -path "$production_dir/chunk_*/results/*.root" -type f \
        -name "$input_pattern" -print | sort -V
)
if (( ${#inputs[@]} < 10 || ${#inputs[@]} > 20 )); then
    echo "Expected 10 to 20 chunk ROOT files (50--100 positions each), found ${#inputs[@]}" >&2
    exit 2
fi

# Refuse partial, overlapping, duplicated, or mixed productions. The generated
# filename is authoritative for the global position interval in each chunk.
expected_start=0
for input in "${inputs[@]}"; do
    filename=$(basename -- "$input")
    if [[ $filename =~ _pos-([0-9]{4})-([0-9]{4})[.]root$ ]]; then
        chunk_start=$((10#${BASH_REMATCH[1]}))
        chunk_end=$((10#${BASH_REMATCH[2]}))
    else
        echo "Cannot read position range from chunk filename: $input" >&2
        exit 2
    fi
    chunk_size=$((chunk_end - chunk_start + 1))
    (( chunk_size >= 50 && chunk_size <= 100 )) || {
        echo "Chunk $input contains $chunk_size positions; expected 50--100" >&2
        exit 2
    }
    (( chunk_start == expected_start )) || {
        echo "Position coverage is not contiguous: expected start $expected_start, found $chunk_start in $input" >&2
        exit 2
    }
    expected_start=$((chunk_end + 1))
done
(( expected_start == 1000 )) || {
    echo "Position coverage ends at $((expected_start - 1)); expected complete range 0--999" >&2
    exit 2
}
[[ -x $map_exe ]] || {
    echo "Updated merger is not built: $map_exe" >&2
    exit 3
}
command -v root >/dev/null 2>&1 || {
    echo "ROOT is not configured; source the ND280 environment first" >&2
    exit 3
}

echo "Scattering length: $scat_length mm"
echo "Scattering model:  $scattering_model"
echo "Mie g / r_fb:      $mie_g / $mie_forward_fraction"
echo "Fibre diameter:    $fibre_diameter mm"
echo "Scintillator RI:   $scintillator_ri"
echo "Merging ${#inputs[@]} matching chunks into $output"
"$map_exe" -R -G "${inputs[0]}" \
    -o "$output" \
    -O "cubeCentre=0 30 910 mm" \
    -O "cubeSide=10 mm" \
    -O "nBins=10" \
    -O "nPhotons=100000" \
    -O "fibreDiameter=$fibre_diameter" \
    "${inputs[@]}"

root -l -b -q "$study_dir/ExportLightMapEfficiency.C(\"$output\",\"$csv\")"

# The valid source-position count depends on fibre diameter.  Count grid
# centres outside all three HOMO fibre cylinders using the same geometry as
# MAKEHOMOPHOTONMAP.exe, and reject a silently wrong denominator.
expected_positions=$(awk -v diameter="$fibre_diameter" 'BEGIN {
    radius2=(diameter/2)^2; count=0
    for (ix=0; ix<10; ++ix) for (iy=0; iy<10; ++iy) for (iz=0; iz<10; ++iz) {
        x=-4.5+ix; y=-4.5+iy; z=-4.5+iz
        dx2=(y+2.5)^2+(z-2.5)^2
        dy2=(x-2.5)^2+(z+2.5)^2
        dz2=(x+2.5)^2+(y-2.5)^2
        if (dx2>radius2+1e-9 && dy2>radius2+1e-9 && dz2>radius2+1e-9) ++count
    }
    print count
}')
actual_positions=$(root -l -b "$output" -e 'auto* tree=(TTree*)_file0->Get("position_efficiency"); std::cout << "LIGHTMAP_POSITION_COUNT " << (tree ? tree->GetEntries() : -1) << std::endl;' -q 2>/dev/null | awk '/LIGHTMAP_POSITION_COUNT/{print $2}')
[[ $actual_positions == "$expected_positions" ]] || {
    echo "Merged map has $actual_positions valid positions; expected $expected_positions for fibre diameter $fibre_diameter mm" >&2
    exit 4
}

echo "Merged light map: $output"
echo "Per-position efficiency: $csv ($actual_positions valid scintillator positions)"
