#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
study_dir=${1:-$script_dir/Studies/LightmapsStudie}
global_fit_analysis="$script_dir/Studies/GlobalFitTest/global_fit_analysis"

[[ -d $study_dir ]] || { echo "Missing study directory: $study_dir" >&2; exit 2; }
for program in python3 root; do
    command -v "$program" >/dev/null || { echo "Missing program: $program" >&2; exit 1; }
done
[[ -x $global_fit_analysis ]] || { echo "Missing executable: $global_fit_analysis" >&2; exit 1; }

found=0
for branch in "$study_dir"/*; do
    [[ -d $branch && -s $branch/flat.root && -s $branch/global_fit.root ]] || continue
    found=1
    echo "[$(basename "$branch")] Standard-reconstruction diagnostics"
    mkdir -p "$branch/plots"
    python3 "$script_dir/plot_overlay.py" "$branch/flat.root" \
        --output-dir "$branch/plots" --all-events
    echo "[$(basename "$branch")] GlobalFit diagnostics"
    (
        cd "$branch"
        "$global_fit_analysis" global_fit.root OUTPUT=global_fit_analysis.root \
            2>&1 | tee 06_global_fit_analysis.log
    )
done
(( found )) || { echo "No complete ROOT-production branches found in: $study_dir" >&2; exit 2; }

python3 "$script_dir/compare_lightmap_results.py" "$study_dir"
