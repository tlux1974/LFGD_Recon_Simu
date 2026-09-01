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
  PRIMARY_INPUT_FILE=input/primary_mu700_center_isotropic_seed12345_1000.csv $0 homo 10

Run this after sourcing the normal ND280++ environment and then selecting the
desired simulation and reconstruction builds, for example:
  source "${SCRIPT_DIR}/switch-nd280geant4sim.sh" local
  source "${SCRIPT_DIR}/switch-hfgrecon.sh" local
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
if [[ "$DIRECTION_MODE" != "fixed" && "$DIRECTION_MODE" != "isotropic" && "$DIRECTION_MODE" != "cone" ]]; then
    echo "DIRECTION_MODE must be 'fixed', 'isotropic', or 'cone'." >&2
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
gps_args=(
    --detector "$DETECTOR" --baseline "$BASELINE" --events "$EVENTS"
    --particle "$PARTICLE" --energy-mev "$ENERGY_MEV"
    --position-mm "$px" "$py" "$pz" --position-frame "$POSITION_FRAME"
    --direction-mode "$DIRECTION_MODE"
    --direction "$dx" "$dy" "$dz"
    --cone-half-angle-deg "${CONE_HALF_ANGLE_DEG:-5}"
    --output "$macro"
)
if [[ -n "${PRIMARY_INPUT_FILE:-}" ]]; then
    [[ -f "$PRIMARY_INPUT_FILE" ]] || {
        echo "Missing primary input file: $PRIMARY_INPUT_FILE" >&2
        exit 2
    }
    gps_args+=(--primary-input "$PRIMARY_INPUT_FILE")
fi
python3 "${SCRIPT_DIR}/generate_gps_macro.py" "${gps_args[@]}"

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
    # Explicit primaries are generated with one /run/beamOn command per
    # event.  Pin the geometry from this input file so oaEvent does not reload
    # the identical multi-million-node HOMO geometry at every run boundary.
    -G "$g4"
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
detresponse_parameter_inputs=()
if [[ -n "${HOMO_LIGHTMAP_FILE:-}" ]]; then
    [[ "$DETECTOR" == "homo" ]] || {
        echo "HOMO_LIGHTMAP_FILE is only valid for the homo detector" >&2
        exit 2
    }
    HOMO_LIGHTMAP_FILE=$(readlink -f "$HOMO_LIGHTMAP_FILE")
    [[ -f "$HOMO_LIGHTMAP_FILE" ]] || {
        echo "Missing HOMO light map: $HOMO_LIGHTMAP_FILE" >&2
        exit 2
    }
    lightmap_parameter_file="${OUTPUT_DIR}/detresponse_lightmap.parameters.dat"
    printf '< detResponseSim.LiquidO.Response.File = %s >\n' \
        "$HOMO_LIGHTMAP_FILE" > "$lightmap_parameter_file"
    detresponse_parameter_inputs+=("$lightmap_parameter_file")
    echo "HOMO light map: $HOMO_LIGHTMAP_FILE"
fi
if [[ -n "${DETRESPONSE_PARAMETER_FILE:-}" ]]; then
    [[ -f "$DETRESPONSE_PARAMETER_FILE" ]] || {
        echo "Missing detector-response parameter file: $DETRESPONSE_PARAMETER_FILE" >&2
        exit 2
    }
    detresponse_parameter_inputs+=("$DETRESPONSE_PARAMETER_FILE")
fi
if [[ "${DISABLE_HOMO_ATTENUATION:-0}" == "1" ]]; then
    echo "Detector-response HOMO fibre attenuation disabled for validation"
    detresponse_parameter_inputs+=(
        "${SCRIPT_DIR}/detresponse_no_homo_attenuation.parameters.dat")
fi
if [[ "${FAST_HOMO_RESPONSE:-0}" == "1" ]]; then
    echo "Detector-response HOMO fast aggregate response enabled"
fi
if (( ${#detresponse_parameter_inputs[@]} > 0 )) \
   || [[ "${FAST_HOMO_RESPONSE:-0}" == "1" ]]; then
    detresponse_parameter_file="${OUTPUT_DIR}/detresponse.parameters.dat"
    : > "$detresponse_parameter_file"
    for input in "${detresponse_parameter_inputs[@]}"; do
        cat "$input" >> "$detresponse_parameter_file"
        printf '\n' >> "$detresponse_parameter_file"
    done
    if [[ "${FAST_HOMO_RESPONSE:-0}" == "1" ]]; then
        printf '< detResponseSim.Homo.FastFiberResponse = 1 >\n' \
            >> "$detresponse_parameter_file"
    fi
    echo "Detector-response parameters: $detresponse_parameter_file"
    detresponse_args+=(-O "par_override=$detresponse_parameter_file")
fi
DETRESPONSESIM.exe "${detresponse_args[@]}" "$g4" \
    2>&1 | tee "${OUTPUT_DIR}/02_detresponse.log"

echo "[3/5] HomoFGD/HFGD 3D-hit and track reconstruction"
hfgrecon_args=(-R -o "$reco")
hfgrecon_parameter_file="${HFGRECON_PARAMETER_FILE:-}"
if [[ "${DISABLE_HOMO_ATTENUATION:-0}" == "1"
      && -z "$hfgrecon_parameter_file" ]]; then
    echo "hfgRecon HOMO attenuation correction disabled for validation"
    hfgrecon_parameter_file="${SCRIPT_DIR}/hfgrecon_no_homo_attenuation.parameters.dat"
fi
if [[ -n "$hfgrecon_parameter_file" ]]; then
    [[ -f "$hfgrecon_parameter_file" ]] || {
        echo "Missing hfgRecon parameter file: $hfgrecon_parameter_file" >&2
        exit 2
    }
    echo "hfgRecon parameters: $hfgrecon_parameter_file"
    hfgrecon_args+=(-O "par_override=$hfgrecon_parameter_file")
fi
HFGRECON.exe "${hfgrecon_args[@]}" "$detresp" \
    2>&1 | tee "${OUTPUT_DIR}/03_hfgrecon.log"

echo "[4/5] Flat diagnostic tree"
LFGDFLATTREE.exe -R -O outfile="$flat" "$reco" \
    2>&1 | tee "${OUTPUT_DIR}/04_flat_tree.log"

echo "[5/5] Overlay plots"
plot_args=("$flat" --output-dir "${OUTPUT_DIR}/plots")
if [[ -n "${PLOT_EVENT_RANGE:-}" ]]; then
    read -r plot_first plot_last <<<"$PLOT_EVENT_RANGE"
    plot_args+=(--event-range "$plot_first" "$plot_last")
fi
python3 "${SCRIPT_DIR}/plot_overlay.py" "${plot_args[@]}"

echo "Done: $OUTPUT_DIR"
