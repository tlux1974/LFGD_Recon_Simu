#!/usr/bin/env bash

set -euo pipefail

[[ $# == 8 ]] || {
    echo "Usage: $0 OCTANT_START POSITION_COUNT NPHOTONS SCAT_VALUE SCAT_UNIT ABS_VALUE ABS_UNIT SEED_BASE" >&2
    exit 2
}

# HTCondor transfers outputs only from its per-job scratch directory.  The
# executable itself can resolve to the shared initialdir on PIC, so do not use
# the executable's directory as the working directory.
job_work_dir=${_CONDOR_SCRATCH_DIR:-$PWD}
cd "$job_work_dir"

# The +++ octant contains the five bin centres above the cube centre on each
# axis: 5 x 5 x 5 = 125 photon-bomb positions.  Keep their original ordering.
awk '$1 > 0 && $2 > 30 && $3 > 910' positions.txt > positions_octant_ppp.txt
[[ $(wc -l < positions_octant_ppp.txt) == 125 ]] || {
    echo "Expected 125 +++ octant positions" >&2
    exit 3
}

octant_start=$1
position_count=$2
[[ $octant_start =~ ^(0|25|50|75|100)$ ]] || {
    echo "OCTANT_START must be one of 0, 25, 50, 75, 100" >&2
    exit 3
}
[[ $position_count == 25 ]] || { echo "POSITION_COUNT must be 25" >&2; exit 3; }

export LIGHTMAP_JOB_TAG="octant_chunk_${octant_start}"
exec ./run_lightmap_chunk.sh "$octant_start" "$position_count" \
    "$3" "$4" "$5" "$6" "$7" "$8" \
    positions_octant_ppp.txt
