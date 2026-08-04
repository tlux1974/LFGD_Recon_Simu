#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SYSTEM="${ND280_SYSTEM:-Linux-AlmaLinux_9.5-gcc_11-x86_64}"
BUILD_DIR="${SCRIPT_DIR}/build/${SYSTEM}"
BIN_DIR="${SCRIPT_DIR}/${SYSTEM}/bin"

if [[ -z "${HFGRECONROOT:-}" ]]; then
    echo "HFGRECONROOT is unset. Source the ND280 setup and switch-hfgrecon.sh first." >&2
    exit 1
fi
if ! command -v nd280-system >/dev/null 2>&1; then
    echo "nd280-system is unavailable. Enter the normal ND280++ container/login environment first." >&2
    echo "Do not rely on a package bin/setup.sh alone; it expects the base environment." >&2
    exit 1
fi
if [[ -z "${CMAKE_PREFIX_PATH:-}" ]]; then
    echo "CMAKE_PREFIX_PATH is unset. Source the normal ND280++ environment first." >&2
    exit 1
fi

cmake -S "${SCRIPT_DIR}/cmake" -B "$BUILD_DIR"
cmake --build "$BUILD_DIR" -j"${JOBS:-4}"
echo "Built under $BUILD_DIR"
echo "Add its bin directory to PATH:"
echo "  export PATH=${BIN_DIR}:\$PATH"
