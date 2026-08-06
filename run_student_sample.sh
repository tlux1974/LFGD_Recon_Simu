#!/usr/bin/env bash
set -euo pipefail

# One-command entry point for the three detector/reconstruction configurations
# used in the LFGD validation study. Run this inside ND280ppCont.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ND280_SYSTEM="${ND280_SYSTEM:-Linux-AlmaLinux_9.5-gcc_11-x86_64}"

usage() {
    cat <<EOF
Usage: $(basename "$0") CONFIG EVENTS [RUN_NAME]

CONFIG is one of:
  hfg-standard      HFGD with its standard reconstruction
  lfg-original-low Historical LFGD peak method with 10 PE thresholds
  lfg-best          Local LFGD 2D clustering plus view-average charge (mode 2)

Examples:
  ./run_student_sample.sh hfg-standard 100
  ./run_student_sample.sh lfg-original-low 100 comparison_a
  ./run_student_sample.sh lfg-best 100 comparison_a

The default particle sample is reproducible isotropic 700 MeV muons from the
detector centre, with random seed 12345. Optional environment overrides include:
  SEED=6789 ENERGY_MEV=1000 PARTICLE=pi+ DIRECTION_MODE=isotropic
  DIRECTION_MODE=fixed DIRECTION="0 0 1"

Every run gets a new directory containing g4.root, detresponse.root,
reco.root, flat.root, logs, and plots/index.html. The first ten events are
plotted; all events are written to every ROOT file.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
if (( $# < 2 )); then
    usage >&2
    exit 2
fi

configuration="$1"
events="$2"
run_name="${3:-$(date +%Y%m%d_%H%M%S)}"
[[ "$events" =~ ^[1-9][0-9]*$ ]] || {
    echo "EVENTS must be a positive integer." >&2
    exit 2
}
[[ "$run_name" =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "RUN_NAME may contain only letters, digits, dot, underscore, and dash." >&2
    exit 2
}

case "$configuration" in
    hfg-standard)
        detector=hfg
        parameter_file=""
        ;;
    lfg-original-low)
        detector=homo
        parameter_file="${SCRIPT_DIR}/student_lfg_original_low.parameters.dat"
        ;;
    lfg-best)
        detector=homo
        parameter_file="${SCRIPT_DIR}/student_lfg_best.parameters.dat"
        ;;
    *)
        echo "Unknown CONFIG: $configuration" >&2
        usage >&2
        exit 2
        ;;
esac

if [[ ! -d /usr/local/t2k/current ]]; then
    echo "Run this script inside the ND280++ Singularity container." >&2
    exit 1
fi

export ND280_SYSTEM
export PATH="/usr/local/t2k/current/nd280SoftwarePilot/scripts:/usr/local/t2k/current/nd280SoftwarePilot/${ND280_SYSTEM}/bin:${PATH}"
# The installed detResponseSim setup checks this even when HAT response is not
# used.  Match the value configured by the container's software pilot.
export ND280_DOWNLOADS="${ND280_DOWNLOADS:-http://nd280.lancs.ac.uk/downloads}"
# Point the setup fragment at the writable local build.  Otherwise it tries to
# refresh templates.index inside the read-only installed container package.
export DETRESPONSESIMROOT="${WORKSPACE_DIR}/SoftProj/detResponseSim"
export DETRESPONSESIMCONFIG="${ND280_SYSTEM}"
set +u
source /usr/local/t2k/current/nd280SoftwareMaster_14.36-plusplus.0.3/bin/setup.sh
set -u

# Use the coordinated local geometry, response, and reconstruction packages.
# Set these paths directly so a fresh checkout does not depend on private
# switch scripts outside this repository.
select_local_package() {
    local root_name="$1"
    local package_dir="$2"
    local package_root="${WORKSPACE_DIR}/SoftProj/${package_dir}"
    [[ -d "$package_root" ]] || {
        echo "Missing local package: $package_root" >&2
        echo "Follow REPOSITORY_SETUP.md before running samples." >&2
        exit 1
    }
    export "${root_name}ROOT=${package_root}"
    export "${root_name}CONFIG=${ND280_SYSTEM}"
    PATH="${package_root}/${ND280_SYSTEM}/bin:${PATH}"
    LD_LIBRARY_PATH="${package_root}/${ND280_SYSTEM}/lib:${LD_LIBRARY_PATH:-}"
}
select_local_package OAEVENT oaEvent
select_local_package OAGEOMINFO oaGeomInfo
select_local_package ND280GEANT4SIM nd280Geant4Sim
select_local_package DETRESPONSESIM detResponseSim
select_local_package HFGRECON hfgrecon
export PATH LD_LIBRARY_PATH
unset -f select_local_package
export PATH="${SCRIPT_DIR}/${ND280_SYSTEM}/bin:${PATH}"
export LD_LIBRARY_PATH="${SCRIPT_DIR}/${ND280_SYSTEM}/lib:${LD_LIBRARY_PATH:-}"

for program in ND280GEANT4SIM.exe DETRESPONSESIM.exe HFGRECON.exe LFGDFLATTREE.exe; do
    command -v "$program" >/dev/null || {
        echo "Missing $program. Follow STUDENT_GUIDE.md section 'Build once'." >&2
        exit 1
    }
done

plot_last=$((events - 1))
(( plot_last > 9 )) && plot_last=9

export DETECTOR="$detector"
export EVENTS="$events"
export SIM_NAME="student_${configuration}_${run_name}"
export DIRECTION_MODE="${DIRECTION_MODE:-isotropic}"
export POSITION_FRAME="${POSITION_FRAME:-plusplus}"
export POSITION_MM="${POSITION_MM:-0 0 1800}"
export SEED="${SEED:-12345}"
export HFGRECON_PARAMETER_FILE="$parameter_file"
export PLOT_EVENT_RANGE="0 $plot_last"

# Freeze the generator state event by event. A shared Geant4 seed alone is
# insufficient: different detector transport consumes different numbers of
# random values and would change later isotropic GPS directions.
actual_output_dir="${SCRIPT_DIR}/output/${detector}_${SIM_NAME}"
if [[ -e "$actual_output_dir" ]]; then
    echo "Refusing to overwrite existing output: $actual_output_dir" >&2
    echo "Choose a different RUN_NAME." >&2
    exit 2
fi
mkdir -p "$actual_output_dir"
if [[ -n "${PRIMARY_INPUT_FILE:-}" ]]; then
    PRIMARY_INPUT_FILE="$(readlink -f "$PRIMARY_INPUT_FILE")"
    [[ -f "$PRIMARY_INPUT_FILE" ]] || {
        echo "Missing primary input file: $PRIMARY_INPUT_FILE" >&2
        exit 2
    }
else
    PRIMARY_INPUT_FILE="${actual_output_dir}/primary_events.csv"
    read -r primary_x primary_y primary_z <<<"$POSITION_MM"
    read -r primary_dx primary_dy primary_dz <<<"${DIRECTION:-0 0 1}"
    python3 "${SCRIPT_DIR}/generate_primary_events.py" \
        --events "$events" --seed "$SEED" --particle "${PARTICLE:-mu-}" \
        --energy-mev "${ENERGY_MEV:-700}" \
        --position-mm "$primary_x" "$primary_y" "$primary_z" \
        --position-frame "$POSITION_FRAME" --direction-mode "$DIRECTION_MODE" \
        --direction "$primary_dx" "$primary_dy" "$primary_dz" \
        --output "$PRIMARY_INPUT_FILE"
fi
export PRIMARY_INPUT_FILE

echo "Configuration: $configuration"
echo "Detector:      $detector"
echo "Events:        $events"
echo "Parameters:    ${parameter_file:-standard package parameters}"
echo "Output:        $actual_output_dir"

"${SCRIPT_DIR}/run_pipeline.sh" "$detector" "$events"

echo
echo "Finished successfully."
echo "ROOT analysis file: $actual_output_dir/flat.root"
echo "Event display:      $actual_output_dir/plots/index.html"
