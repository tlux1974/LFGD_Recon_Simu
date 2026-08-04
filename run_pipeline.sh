#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

usage() {
    cat <<EOF
Usage: $(basename "$0") [homo|hfg] [number-of-events]

Runs 700 MeV muons from a fixed point through:
  GPS macro -> Geant4 -> detector response -> HFGRECON -> flat tree -> plots

Defaults come from config.sh and can be overridden, for example:
  POSITION_FRAME=plusplus POSITION_MM="0 0 1800" DIRECTION="1 0 1" $0 homo 20
  POSITION_FRAME=global POSITION_MM="0 30 910" $0 homo 20
  DIRECTION_MODE=isotropic $0 homo 1000

Run this after sourcing the normal ND280++ environment and then selecting the
desired hfgRecon build with switch-hfgrecon.sh.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
DETECTOR="${1:-$DETECTOR}"
EVENTS="${2:-$EVENTS}"
if [[ "$DETECTOR" != "homo" && "$DETECTOR" != "hfg" ]]; then
    echo "Detector must be 'homo' or 'hfg'." >&2
    exit 2
fi
if [[ "$DIRECTION_MODE" != "fixed" && "$DIRECTION_MODE" != "isotropic" ]]; then
    echo "DIRECTION_MODE must be 'fixed' or 'isotropic'." >&2
    exit 2
fi
if [[ "$POSITION_FRAME" != "plusplus" && "$POSITION_FRAME" != "global" ]]; then
    echo "POSITION_FRAME must be 'plusplus' or 'global'." >&2
    exit 2
fi

for program in ND280GEANT4SIM.exe DETRESPONSESIM.exe HFGRECON.exe LFGDFLATTREE.exe python3; do
    command -v "$program" >/dev/null || {
        echo "Missing $program. Source ND280 setup and run build_flat_treemaker.sh." >&2
        exit 1
    }
done

OUTPUT_DIR="${SCRIPT_DIR}/output/${DETECTOR}_${SIM_NAME}"
mkdir -p "$OUTPUT_DIR"
macro="${OUTPUT_DIR}/gps.mac"
g4="${OUTPUT_DIR}/g4.root"
detresp="${OUTPUT_DIR}/detresponse.root"
reco="${OUTPUT_DIR}/reco.root"
flat="${OUTPUT_DIR}/flat.root"

read -r px py pz <<<"$POSITION_MM"
read -r dx dy dz <<<"$DIRECTION"
python3 "${SCRIPT_DIR}/generate_gps_macro.py" \
    --detector "$DETECTOR" --baseline "$BASELINE" --events "$EVENTS" \
    --particle "$PARTICLE" --energy-mev "$ENERGY_MEV" \
    --position-mm "$px" "$py" "$pz" --position-frame "$POSITION_FRAME" \
    --direction-mode "$DIRECTION_MODE" \
    --direction "$dx" "$dy" "$dz" \
    --output "$macro"

echo "[1/5] Geant4: $DETECTOR, $EVENTS x $PARTICLE at ${ENERGY_MEV} MeV, ${DIRECTION_MODE} directions, vertex $POSITION_MM mm ($POSITION_FRAME frame), seed $SEED"
ND280GEANT4SIM.exe -s "$SEED" -o "${g4%.root}" "$macro" \
    2>&1 | tee "${OUTPUT_DIR}/01_geant4.log"

echo "[2/5] Detector response"
detector_option="upgrade-nd280plus-homo"
[[ "$DETECTOR" == "hfg" ]] && detector_option="upgrade-nd280plus"
# The broad upgrade options also enable response simulation for SciFi, SWD,
# SFG, ECAL, SMRD, TOF, TPC, and HAT.  They are unrelated to this focused
# Homo/HFG comparison and can produce errors from hits in those subsystems
# (notably TPlusPlusScint's "Hit segment is not in a fiber").
detresponse_args=(
    -R
    -O "randomseed=$SEED"
    -O "$detector_option"
    -O disable-sfg
    -O disable-swd
    -O disable-scifi
    -O disable-ecal
    -O disable-smrd
    -O disable-tof
    -O disable-tpc
    -O disable-hat
    -o "$detresp"
)
DETRESPONSESIM.exe "${detresponse_args[@]}" "$g4" \
    2>&1 | tee "${OUTPUT_DIR}/02_detresponse.log"

echo "[3/5] HomoFGD/HFGD 3D-hit and track reconstruction"
HFGRECON.exe -R -o "$reco" "$detresp" \
    2>&1 | tee "${OUTPUT_DIR}/03_hfgrecon.log"

echo "[4/5] Flat diagnostic tree"
LFGDFLATTREE.exe -R -O outfile="$flat" "$reco" \
    2>&1 | tee "${OUTPUT_DIR}/04_flat_tree.log"

echo "[5/5] Overlay plots"
python3 "${SCRIPT_DIR}/plot_overlay.py" "$flat" --output-dir "${OUTPUT_DIR}/plots"

echo "Done: $OUTPUT_DIR"
