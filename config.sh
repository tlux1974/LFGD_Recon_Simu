#!/usr/bin/env bash
# Shared defaults. Every value can be overridden in the environment.

DETECTOR="${DETECTOR:-homo}"       # homo (LFGD) or hfg
EVENTS="${EVENTS:-10}"
PARTICLE="${PARTICLE:-mu-}"
ENERGY_MEV="${ENERGY_MEV:-700}"
# Use the familiar PlusPlusTracker layout frame by default.  The macro
# generator converts its detector centre (0,0,1800) to global (0,30,910) mm,
# since Geant4 GPS itself always consumes global ND280 coordinates.
POSITION_FRAME="${POSITION_FRAME:-plusplus}" # plusplus or global
POSITION_MM="${POSITION_MM:-0 0 1800}"
DIRECTION_MODE="${DIRECTION_MODE:-fixed}" # fixed or isotropic
DIRECTION="${DIRECTION:-0 0 1}"
SIM_NAME="${SIM_NAME:-mu700_${DIRECTION_MODE}}"
BASELINE="${BASELINE:-baseline-2024-plusplus}"
SEED="${SEED:-12345}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR:-${SCRIPT_DIR}/output/${DETECTOR}_${SIM_NAME}}"
