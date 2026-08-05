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
set +u
source /usr/local/t2k/current/nd280SoftwareMaster_14.36-plusplus.0.3/bin/setup.sh
set -u

# Use the coordinated local geometry, response, and reconstruction packages.
source "${WORKSPACE_DIR}/switch-oaevent.sh" local
source "${WORKSPACE_DIR}/switch-oageominfo.sh" local
source "${WORKSPACE_DIR}/switch-nd280geant4sim.sh" local
source "${WORKSPACE_DIR}/switch-detresponsesim.sh" local
source "${WORKSPACE_DIR}/switch-hfgrecon.sh" local
export PATH="${SCRIPT_DIR}/${ND280_SYSTEM}/bin:${PATH}"
export LD_LIBRARY_PATH="${SCRIPT_DIR}/${ND280_SYSTEM}/lib:${LD_LIBRARY_PATH:-}"

for program in ND280GEANT4SIM.exe DETRESPONSESIM.exe HFGRECON.exe LFGDFLATTREE.exe; do
    command -v "$program" >/dev/null || {
        echo "Missing $program. Follow STUDENT_GUIDE.md section 'Build once'." >&2
        exit 1
    }
done

output_dir="${SCRIPT_DIR}/output/student_${configuration}_${run_name}"
if [[ -e "$output_dir" ]]; then
    echo "Refusing to overwrite existing output: $output_dir" >&2
    echo "Choose a different RUN_NAME." >&2
    exit 2
fi

plot_last=$((events - 1))
(( plot_last > 9 )) && plot_last=9

export DETECTOR="$detector"
export EVENTS="$events"
export OUTPUT_DIR="$output_dir"
export SIM_NAME="student_${configuration}_${run_name}"
export DIRECTION_MODE="${DIRECTION_MODE:-isotropic}"
export POSITION_FRAME="${POSITION_FRAME:-plusplus}"
export POSITION_MM="${POSITION_MM:-0 0 1800}"
export SEED="${SEED:-12345}"
export HFGRECON_PARAMETER_FILE="$parameter_file"
export PLOT_EVENT_RANGE="0 $plot_last"

echo "Configuration: $configuration"
echo "Detector:      $detector"
echo "Events:        $events"
echo "Parameters:    ${parameter_file:-standard package parameters}"
echo "Output:        $output_dir"

"${SCRIPT_DIR}/run_pipeline.sh" "$detector" "$events"

echo
echo "Finished successfully."
echo "ROOT analysis file: $output_dir/flat.root"
echo "Event display:      $output_dir/plots/index.html"
