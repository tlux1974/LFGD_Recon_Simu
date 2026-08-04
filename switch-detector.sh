#!/usr/bin/env bash
# Source this script so DETECTOR remains set in the current shell:
#   source switch-detector.sh lfgd
#   source switch-detector.sh hfgd
#   source switch-detector.sh status

_lfgd_detector_help() {
    cat <<EOF
Select the detector design used by the reconstruction simulation.

Usage:
  source ${BASH_SOURCE[0]} lfgd     Select LFGD (internal model: homo)
  source ${BASH_SOURCE[0]} hfgd     Select HFGD (internal model: hfg)
  source ${BASH_SOURCE[0]} status   Show the current selection
  source ${BASH_SOURCE[0]} help     Show this help

Both designs use the same configured position. The generated Geant4 macro
enables the selected detector and disables the alternative design.
EOF
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "This script must be sourced so it can modify the current shell." >&2
    echo "Use: source ${BASH_SOURCE[0]} {lfgd|hfgd|status}" >&2
    exit 2
fi

case "${1:-help}" in
    lfgd|homo)
        export DETECTOR=homo
        echo "Selected detector: LFGD (ND280 model: homo)"
        ;;
    hfgd|hfg)
        export DETECTOR=hfg
        echo "Selected detector: HFGD (ND280 model: hfg)"
        ;;
    status)
        case "${DETECTOR:-homo}" in
            homo) echo "Selected detector: LFGD (ND280 model: homo)" ;;
            hfg)  echo "Selected detector: HFGD (ND280 model: hfg)" ;;
            *)    echo "Unknown DETECTOR value: ${DETECTOR}" >&2; return 1 ;;
        esac
        ;;
    help|-h|--help)
        _lfgd_detector_help
        ;;
    *)
        echo "Unknown detector: $1" >&2
        _lfgd_detector_help >&2
        unset -f _lfgd_detector_help
        return 2
        ;;
esac

unset -f _lfgd_detector_help
