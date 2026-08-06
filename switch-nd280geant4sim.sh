#!/usr/bin/env bash
# Source this file to select the nd280Geant4Sim package used by this shell:
#   source switch-nd280geant4sim.sh local
#   source switch-nd280geant4sim.sh original
#   source switch-nd280geant4sim.sh status
#   source switch-nd280geant4sim.sh help

_g4_switch_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
_g4_workspace="$(cd -- "${_g4_switch_dir}/.." && pwd)"
_g4_local_root="${_g4_workspace}/SoftProj/nd280Geant4Sim"
_g4_original_root="${_g4_workspace}/ND280ppCont/usr/local/t2k/current/nd280Geant4Sim_7.17-plusplus.1"

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "This script must be sourced so it can modify the current shell." >&2
    echo "Usage: source ${BASH_SOURCE[0]} {local|original|status}" >&2
    exit 2
fi

_g4_system="${ND280_SYSTEM:-}"
if [[ -z "$_g4_system" ]]; then
    for _g4_candidate in "$_g4_original_root"/Linux-*; do
        if [[ -d "$_g4_candidate" ]]; then
            _g4_system="${_g4_candidate##*/}"
            break
        fi
    done
fi

if [[ -z "$_g4_system" ]]; then
    echo "Cannot determine ND280_SYSTEM. Set it before sourcing this script." >&2
    unset _g4_local_root _g4_original_root _g4_system _g4_candidate
    return 1
fi

_g4_remove_path() {
    local _g4_value=":${1:-}:"
    local _g4_entry
    shift
    for _g4_entry in "$@"; do
        _g4_value="${_g4_value//:${_g4_entry}:/:}"
    done
    _g4_value="${_g4_value#:}"
    _g4_value="${_g4_value%:}"
    printf '%s' "$_g4_value"
}

_g4_show_status() {
    echo "ND280GEANT4SIMROOT=${ND280GEANT4SIMROOT:-<unset>}"
    echo "ND280_SYSTEM=$_g4_system"
    if [[ -n "${ND280GEANT4SIMROOT:-}" ]]; then
        echo "bin=${ND280GEANT4SIMROOT}/${_g4_system}/bin"
        echo "lib=${ND280GEANT4SIMROOT}/${_g4_system}/lib"
    fi
    if command -v ND280GEANT4SIM.exe >/dev/null 2>&1; then
        echo "executable=$(command -v ND280GEANT4SIM.exe)"
    else
        echo "executable=<not found>"
    fi
}

_g4_show_help() {
    cat <<EOF
Switch between the local and installed nd280Geant4Sim packages.

This script must be sourced so it can update the current shell:

  source ${BASH_SOURCE[0]} local
      Use the development clone at:
      $_g4_local_root

  source ${BASH_SOURCE[0]} original
      Use the installed container package at:
      $_g4_original_root

  source ${BASH_SOURCE[0]} status
      Show the selected package, build paths, and resolved executable.

The script updates ND280GEANT4SIMROOT, ND280GEANT4SIMCONFIG,
ND280_SYSTEM, PATH, and LD_LIBRARY_PATH. The selected package must be
built under <package>/\${ND280_SYSTEM}/lib.
EOF
}

case "${1:-help}" in
    local)
        _g4_selected_root="$_g4_local_root"
        ;;
    original|installed)
        _g4_selected_root="$_g4_original_root"
        ;;
    status)
        _g4_show_status
        unset -f _g4_remove_path _g4_show_status _g4_show_help
        unset _g4_local_root _g4_original_root _g4_system _g4_candidate
        return 0
        ;;
    help|-h|--help)
        _g4_show_help
        unset -f _g4_remove_path _g4_show_status _g4_show_help
        unset _g4_local_root _g4_original_root _g4_system _g4_candidate
        return 0
        ;;
    *)
        echo "Unknown selection: $1" >&2
        _g4_show_help >&2
        unset -f _g4_remove_path _g4_show_status _g4_show_help
        unset _g4_local_root _g4_original_root _g4_system _g4_candidate
        return 2
        ;;
esac

if [[ ! -d "$_g4_selected_root" ]]; then
    echo "nd280Geant4Sim directory not found: $_g4_selected_root" >&2
    unset -f _g4_remove_path _g4_show_status _g4_show_help
    unset _g4_local_root _g4_original_root _g4_selected_root
    unset _g4_system _g4_candidate
    return 1
fi

_g4_local_bin="${_g4_local_root}/${_g4_system}/bin"
_g4_local_lib="${_g4_local_root}/${_g4_system}/lib"
_g4_original_bin="${_g4_original_root}/${_g4_system}/bin"
_g4_original_lib="${_g4_original_root}/${_g4_system}/lib"
_g4_selected_bin="${_g4_selected_root}/${_g4_system}/bin"
_g4_selected_lib="${_g4_selected_root}/${_g4_system}/lib"

if [[ ! -d "$_g4_selected_lib" ]]; then
    echo "Warning: selected package has not been built for $_g4_system:" >&2
    echo "  missing $_g4_selected_lib" >&2
fi

PATH="$(_g4_remove_path "${PATH:-}" \
    "$_g4_local_bin" "$_g4_original_bin")"
LD_LIBRARY_PATH="$(_g4_remove_path "${LD_LIBRARY_PATH:-}" \
    "$_g4_local_lib" "$_g4_original_lib")"

PATH="${_g4_selected_bin}${PATH:+:${PATH}}"
LD_LIBRARY_PATH="${_g4_selected_lib}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
ND280GEANT4SIMROOT="$_g4_selected_root"
ND280GEANT4SIMCONFIG="$_g4_system"
ND280_SYSTEM="$_g4_system"

export PATH LD_LIBRARY_PATH ND280GEANT4SIMROOT ND280GEANT4SIMCONFIG ND280_SYSTEM

echo "Selected nd280Geant4Sim: $ND280GEANT4SIMROOT"
echo "Library path:            $_g4_selected_lib"

unset -f _g4_remove_path _g4_show_status _g4_show_help
unset _g4_local_root _g4_original_root _g4_selected_root
unset _g4_local_bin _g4_local_lib _g4_original_bin _g4_original_lib
unset _g4_selected_bin _g4_selected_lib _g4_system _g4_candidate
