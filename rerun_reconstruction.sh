#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<EOF
Usage: $(basename "$0") DETRESPONSE_ROOT OUTPUT_DIR [HFGRECON_OPTIONS...]

Rerun the fibre-to-3D-hit reconstruction from an existing detector-response
file.  Geant4 and detector response are not rerun.  OUTPUT_DIR receives:

  reco.root
  flat.root
  03_hfgrecon.log
  04_flat_tree.log
  plots/

Examples:
  $(basename "$0") \
    output/homo_mu700_isotropic/detresponse.root \
    output/homo_mu700_isotropic/reco_test

Additional arguments are passed directly to HFGRECON.exe.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
if (( $# < 2 )); then
    usage
    exit 2
fi

input="$1"
output_dir="$2"
shift 2

if [[ ! -f "$input" ]]; then
    echo "Detector-response input does not exist: $input" >&2
    exit 2
fi

for program in HFGRECON.exe LFGDFLATTREE.exe python3; do
    command -v "$program" >/dev/null || {
        echo "Missing $program. Source ND280 setup and select/build local hfgRecon first." >&2
        exit 1
    }
done

mkdir -p "$output_dir"
reco="$output_dir/reco.root"
flat="$output_dir/flat.root"

echo "[1/3] Fibre-to-3D-hit reconstruction"
echo "Input:  $input"
echo "Output: $reco"
HFGRECON.exe -R "$@" -o "$reco" "$input" \
    2>&1 | tee "$output_dir/03_hfgrecon.log"

echo "[2/3] Flat diagnostic tree"
LFGDFLATTREE.exe -R -O outfile="$flat" "$reco" \
    2>&1 | tee "$output_dir/04_flat_tree.log"

echo "[3/3] Event display (event 0)"
python3 "$SCRIPT_DIR/plot_overlay.py" "$flat" \
    --output-dir "$output_dir/plots"

echo "Done: $output_dir"
echo "Browse: $output_dir/plots/index.html"
