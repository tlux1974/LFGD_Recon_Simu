#!/usr/bin/env bash
set -euo pipefail

# Run a fresh 100-event HOMO/LFGD sample with the complete local virtual-cube
# geometry chain. This script is intended to be run inside ND280ppCont.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ND280_SYSTEM="${ND280_SYSTEM:-Linux-AlmaLinux_9.5-gcc_11-x86_64}"

if [[ ! -d /usr/local/t2k/current ]]; then
    echo "Run this script inside the ND280++ Singularity container." >&2
    echo "For example: singularity shell ND280ppCont" >&2
    exit 1
fi

# Make nd280-system available before sourcing the standard master setup.
export ND280_SYSTEM
export PATH="/usr/local/t2k/current/nd280SoftwarePilot/scripts:/usr/local/t2k/current/nd280SoftwarePilot/${ND280_SYSTEM}/bin:${PATH}"
# Some installed ND280 setup fragments inspect optional variables without
# defining them first (notably TEMPLATESFILES). Do not let our nounset policy
# turn those optional checks into fatal errors.
set +u
source /usr/local/t2k/current/nd280SoftwareMaster_14.36-plusplus.0.3/bin/setup.sh
set -u

# Select the coordinated local versions. hfgRecon was configured and rebuilt
# against the local oaEvent and oaGeomInfo packages.
source "${WORKSPACE_DIR}/switch-oaevent.sh" local
source "${WORKSPACE_DIR}/switch-oageominfo.sh" local
source "${WORKSPACE_DIR}/switch-nd280geant4sim.sh" local
source "${WORKSPACE_DIR}/switch-detresponsesim.sh" local
source "${WORKSPACE_DIR}/switch-hfgrecon.sh" local

# The diagnostic flat-tree executable is built directly in LFGD_Recon_Simu.
export PATH="${SCRIPT_DIR}/${ND280_SYSTEM}/bin:${PATH}"
export LD_LIBRARY_PATH="${SCRIPT_DIR}/${ND280_SYSTEM}/lib:${LD_LIBRARY_PATH:-}"

run_name="${1:-mu700_virtual_geometry_100ev_$(date +%Y%m%d_%H%M%S)}"
export OUTPUT_DIR="${SCRIPT_DIR}/output/${run_name}"

if [[ -e "${OUTPUT_DIR}" ]]; then
    echo "Refusing to overwrite existing output: ${OUTPUT_DIR}" >&2
    echo "Pass a different folder name as the first argument." >&2
    exit 2
fi

export SIM_NAME="${run_name}"
export DIRECTION_MODE=isotropic
export POSITION_FRAME=plusplus
export POSITION_MM="0 0 1800"
export SEED=12345

echo "Running 100 isotropic 700 MeV muons with the local virtual geometry."
echo "Output: ${OUTPUT_DIR}"

"${SCRIPT_DIR}/run_pipeline.sh" homo 100

echo "Completed. Results are in: ${OUTPUT_DIR}"
