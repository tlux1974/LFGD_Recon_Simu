#!/usr/bin/env bash
# Source this file: source switch-detresponsesim.sh {local|original|status}
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "This script must be sourced." >&2
    exit 2
fi
_nd280_switch_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
_nd280_workspace="$(cd -- "${_nd280_switch_dir}/.." && pwd)"
source "${_nd280_switch_dir}/switch-nd280-package-common.sh"
_nd280_switch_package detResponseSim DETRESPONSESIM \
    "${_nd280_workspace}/SoftProj/detResponseSim" \
    "${_nd280_workspace}/ND280ppCont/usr/local/t2k/current/detResponseSim_7.20-plusplus.0.1" \
    "${1:-help}"
_switch_result=$?
unset -f _nd280_switch_package
unset _switch_result _nd280_switch_dir _nd280_workspace
