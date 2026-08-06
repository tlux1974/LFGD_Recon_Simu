#!/usr/bin/env bash
# Shared implementation for the package switch scripts. Source a public
# switch-*.sh wrapper instead of using this file directly.

_nd280_switch_package() {
    local package_name="$1"
    local env_prefix="$2"
    local local_root="$3"
    local original_root="$4"
    local selection="${5:-help}"
    local system="${ND280_SYSTEM:-}"
    local candidate selected_root local_bin local_lib original_bin original_lib
    local selected_bin selected_lib value entry executable

    if [[ -z "$system" ]]; then
        for candidate in "$original_root"/Linux-*; do
            if [[ -d "$candidate" ]]; then
                system="${candidate##*/}"
                break
            fi
        done
    fi
    if [[ -z "$system" ]]; then
        echo "Cannot determine ND280_SYSTEM. Set it before sourcing the switch script." >&2
        return 1
    fi

    case "$selection" in
        local) selected_root="$local_root" ;;
        original|installed) selected_root="$original_root" ;;
        status)
            eval "echo ${env_prefix}ROOT=\${${env_prefix}ROOT:-<unset>}"
            echo "ND280_SYSTEM=$system"
            eval "selected_root=\${${env_prefix}ROOT:-}"
            if [[ -n "$selected_root" ]]; then
                echo "bin=$selected_root/$system/bin"
                echo "lib=$selected_root/$system/lib"
            fi
            return 0
            ;;
        help|-h|--help)
            echo "Usage: source switch-${package_name,,}.sh {local|original|status}"
            echo "  local:    $local_root"
            echo "  original: $original_root"
            return 0
            ;;
        *)
            echo "Unknown selection: $selection" >&2
            return 2
            ;;
    esac

    if [[ ! -d "$selected_root" ]]; then
        echo "$package_name directory not found: $selected_root" >&2
        return 1
    fi

    local_bin="$local_root/$system/bin"
    local_lib="$local_root/$system/lib"
    original_bin="$original_root/$system/bin"
    original_lib="$original_root/$system/lib"
    selected_bin="$selected_root/$system/bin"
    selected_lib="$selected_root/$system/lib"

    value=":${PATH:-}:"
    for entry in "$local_bin" "$original_bin"; do value="${value//:${entry}:/:}"; done
    value="${value#:}"; PATH="${value%:}"
    value=":${LD_LIBRARY_PATH:-}:"
    for entry in "$local_lib" "$original_lib"; do value="${value//:${entry}:/:}"; done
    value="${value#:}"; LD_LIBRARY_PATH="${value%:}"

    PATH="$selected_bin${PATH:+:$PATH}"
    LD_LIBRARY_PATH="$selected_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    printf -v "${env_prefix}ROOT" '%s' "$selected_root"
    printf -v "${env_prefix}CONFIG" '%s' "$system"
    ND280_SYSTEM="$system"
    export PATH LD_LIBRARY_PATH ND280_SYSTEM
    export "${env_prefix}ROOT" "${env_prefix}CONFIG"

    if [[ ! -d "$selected_lib" ]]; then
        echo "Warning: $package_name has not been built for $system:" >&2
        echo "  missing $selected_lib" >&2
    fi
    executable="$(find "$selected_bin" -maxdepth 1 -type f -perm -u+x 2>/dev/null | head -1)"
    echo "Selected $package_name: $selected_root"
    echo "Library path:          $selected_lib"
    [[ -n "$executable" ]] && echo "Example executable:    $executable"
}
