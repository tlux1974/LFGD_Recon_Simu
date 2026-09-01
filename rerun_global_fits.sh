#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fit_exe="$script_dir/Studies/GlobalFitTest/global_light_fit"

usage() {
    cat <<EOF
Usage: $(basename "$0") STUDY_DIR OUTPUT_NAME KEY=value [KEY=value ...]

Reruns only GlobalFit for every completed light-map branch. OUTPUT_NAME is a
new filename such as global_fit_median_seed.root; existing outputs are never
overwritten. The KEY=value arguments are passed directly to global_light_fit.
EOF
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then usage; exit 0; fi
(( $# >= 3 )) || { usage >&2; exit 2; }
study_dir=$(realpath "$1")
output_name=$2
shift 2
[[ $output_name == *.root && $output_name != */* ]] || { echo "OUTPUT_NAME must be a ROOT filename without a directory" >&2; exit 2; }
[[ -x $fit_exe ]] || { echo "Missing executable: $fit_exe" >&2; exit 1; }
for option in "$@"; do [[ $option == *=* ]] || { echo "Fit option must use KEY=value: $option" >&2; exit 2; }; done

found=0
for branch in "$study_dir"/*; do
    [[ -d $branch && -s $branch/flat.root && -s $branch/lightmap_path.txt ]] || continue
    found=1
    output="$branch/$output_name"
    [[ ! -e $output ]] || { echo "Refusing to overwrite: $output" >&2; exit 2; }
    lightmap=$(<"$branch/lightmap_path.txt")
    [[ -s $lightmap ]] || { echo "Missing light map for $branch: $lightmap" >&2; exit 2; }
    log="$branch/${output_name%.root}.log"
    echo "[$(basename "$branch")] $output_name"
    echo "GlobalFit light map: $lightmap" | tee "$log"
    "$fit_exe" "$branch/flat.root" "$output" "$lightmap" "$@" 2>&1 | tee -a "$log"
done
(( found )) || { echo "No completed flat-tree branches found in: $study_dir" >&2; exit 2; }
