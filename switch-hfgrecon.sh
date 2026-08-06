#!/usr/bin/env bash
# Source this file to select the hfgRecon package used by the current shell:
#   source switch-hfgrecon.sh local
#   source switch-hfgrecon.sh original
#   source switch-hfgrecon.sh status
#   source switch-hfgrecon.sh help

_hfg_switch_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
_hfg_workspace="$(cd -- "${_hfg_switch_dir}/.." && pwd)"
_hfg_local_root="${_hfg_workspace}/SoftProj/hfgrecon"
_hfg_original_root="${_hfg_workspace}/ND280ppCont/usr/local/t2k/current/hfgRecon_0.0-plusplus"

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "This script must be sourced so it can modify the current shell." >&2
    echo "Usage: source ${BASH_SOURCE[0]} {local|original|status}" >&2
    exit 2
fi

_hfg_system="${ND280_SYSTEM:-}"
if [[ -z "$_hfg_system" ]]; then
    for _hfg_candidate in "$_hfg_original_root"/Linux-*; do
        if [[ -d "$_hfg_candidate" ]]; then
            _hfg_system="${_hfg_candidate##*/}"
            break
        fi
    done
fi

if [[ -z "$_hfg_system" ]]; then
    echo "Cannot determine ND280_SYSTEM. Set it before sourcing this script." >&2
    unset _hfg_local_root _hfg_original_root _hfg_system _hfg_candidate
    return 1
fi

_hfg_remove_path() {
    local _hfg_value=":${1:-}:"
    local _hfg_entry
    shift
    for _hfg_entry in "$@"; do
        _hfg_value="${_hfg_value//:${_hfg_entry}:/:}"
    done
    _hfg_value="${_hfg_value#:}"
    _hfg_value="${_hfg_value%:}"
    printf '%s' "$_hfg_value"
}

_hfg_show_status() {
    echo "HFGRECONROOT=${HFGRECONROOT:-<unset>}"
    echo "ND280_SYSTEM=$_hfg_system"
    if [[ -n "${HFGRECONROOT:-}" ]]; then
        echo "bin=${HFGRECONROOT}/${_hfg_system}/bin"
        echo "lib=${HFGRECONROOT}/${_hfg_system}/lib"
    fi
}

_hfg_show_help() {
    cat <<EOF
Switch between the local and installed hfgRecon packages.

This script must be sourced so it can update the current shell:

  source ${BASH_SOURCE[0]} local
      Use the development clone at:
      $_hfg_local_root

  source ${BASH_SOURCE[0]} original
      Use the installed package at:
      $_hfg_original_root

  source ${BASH_SOURCE[0]} status
      Show the currently selected package and build paths.

  source ${BASH_SOURCE[0]} help
      Show this help text.

The script updates HFGRECONROOT, ND280_SYSTEM, PATH, and
LD_LIBRARY_PATH. The selected package must be built under:

  <package>/\${ND280_SYSTEM}/lib

The script warns when that library directory is missing.
EOF
}

case "${1:-help}" in
    local)
        _hfg_selected_root="$_hfg_local_root"
        ;;
    original|installed)
        _hfg_selected_root="$_hfg_original_root"
        ;;
    status)
        _hfg_show_status
        unset -f _hfg_remove_path _hfg_show_status _hfg_show_help
        unset _hfg_local_root _hfg_original_root _hfg_system _hfg_candidate
        return 0
        ;;
    help|-h|--help)
        _hfg_show_help
        unset -f _hfg_remove_path _hfg_show_status _hfg_show_help
        unset _hfg_local_root _hfg_original_root _hfg_system _hfg_candidate
        return 0
        ;;
    *)
        echo "Unknown selection: $1" >&2
        _hfg_show_help >&2
        unset -f _hfg_remove_path _hfg_show_status _hfg_show_help
        unset _hfg_local_root _hfg_original_root _hfg_system _hfg_candidate
        return 2
        ;;
esac

if [[ ! -d "$_hfg_selected_root" ]]; then
    echo "hfgRecon directory not found: $_hfg_selected_root" >&2
    unset -f _hfg_remove_path _hfg_show_status _hfg_show_help
    unset _hfg_local_root _hfg_original_root _hfg_selected_root
    unset _hfg_system _hfg_candidate
    return 1
fi

_hfg_local_bin="${_hfg_local_root}/${_hfg_system}/bin"
_hfg_local_lib="${_hfg_local_root}/${_hfg_system}/lib"
_hfg_original_bin="${_hfg_original_root}/${_hfg_system}/bin"
_hfg_original_lib="${_hfg_original_root}/${_hfg_system}/lib"
_hfg_selected_bin="${_hfg_selected_root}/${_hfg_system}/bin"
_hfg_selected_lib="${_hfg_selected_root}/${_hfg_system}/lib"

if [[ ! -d "$_hfg_selected_lib" ]]; then
    echo "Warning: selected package has not been built for $_hfg_system:" >&2
    echo "  missing $_hfg_selected_lib" >&2
fi

PATH="$(_hfg_remove_path "${PATH:-}" \
    "$_hfg_local_bin" "$_hfg_original_bin")"
LD_LIBRARY_PATH="$(_hfg_remove_path "${LD_LIBRARY_PATH:-}" \
    "$_hfg_local_lib" "$_hfg_original_lib")"

PATH="${_hfg_selected_bin}${PATH:+:${PATH}}"
LD_LIBRARY_PATH="${_hfg_selected_lib}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
HFGRECONROOT="$_hfg_selected_root"
ND280_SYSTEM="$_hfg_system"

export PATH LD_LIBRARY_PATH HFGRECONROOT ND280_SYSTEM

echo "Selected hfgRecon: $HFGRECONROOT"
echo "Library path:      $_hfg_selected_lib"

unset -f _hfg_remove_path _hfg_show_status _hfg_show_help
unset _hfg_local_root _hfg_original_root _hfg_selected_root
unset _hfg_local_bin _hfg_local_lib _hfg_original_bin _hfg_original_lib
unset _hfg_selected_bin _hfg_selected_lib _hfg_system _hfg_candidate
