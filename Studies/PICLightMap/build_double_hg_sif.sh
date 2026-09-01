#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: build_double_hg_sif.sh OUTPUT.sif

Install the already-built local nd280Geant4Sim executable/library into the
writable ND280ppCont sandbox, preserve the previous files as *.pre-double-hg,
then convert that sandbox into a new immutable SIF. OUTPUT.sif must not exist.

The original ND280Reco_03.sif is never modified.
EOF
}

if [[ ${1:-} == --help || ${1:-} == -h || $# -ne 1 ]]; then
    usage
    [[ ${1:-} == --help || ${1:-} == -h ]] && exit 0
    exit 2
fi

study_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
workspace=$(cd -- "$study_dir/../../.." && pwd)
sandbox="$workspace/ND280ppCont"
build="$workspace/SoftProj/nd280Geant4Sim/Linux-AlmaLinux_9.5-gcc_11-x86_64"
install_dir="$sandbox/usr/local/t2k/current/nd280Geant4Sim_7.17-plusplus.1/Linux-AlmaLinux_9.5-gcc_11-x86_64"
output=$1
[[ $output = /* ]] || output="$PWD/$output"

[[ -d $sandbox/.singularity.d ]] || { echo "Not a Singularity sandbox: $sandbox" >&2; exit 3; }
[[ -x $build/bin/ND280GEANT4SIM.exe ]] || { echo "Missing local build: $build/bin/ND280GEANT4SIM.exe" >&2; exit 3; }
[[ -s $build/lib/libnd280Geant4Sim.so ]] || { echo "Missing local build: $build/lib/libnd280Geant4Sim.so" >&2; exit 3; }
[[ -d $install_dir/bin && -d $install_dir/lib ]] || { echo "Missing installed package in sandbox: $install_dir" >&2; exit 3; }
[[ ! -e $output ]] || { echo "Refusing to overwrite existing output: $output" >&2; exit 4; }

for relative in bin/ND280GEANT4SIM.exe lib/libnd280Geant4Sim.so; do
    installed="$install_dir/$relative"
    backup="${installed}.pre-double-hg"
    [[ -e $backup ]] || cp -p "$installed" "$backup"
    install -m 775 "$build/$relative" "$installed"
done

echo "Installed local double-HG build in sandbox:"
sha256sum "$install_dir/bin/ND280GEANT4SIM.exe" "$install_dir/lib/libnd280Geant4Sim.so"
grep -a -q liquidOMieForwardFraction "$install_dir/lib/libnd280Geant4Sim.so" || {
    echo "Installed library does not contain the new Mie command" >&2
    exit 5
}

echo "Building immutable SIF: $output"
singularity build --fakeroot "$output" "$sandbox"

echo "Verifying files inside the SIF"
singularity exec "$output" sha256sum \
    /usr/local/t2k/current/nd280Geant4Sim_7.17-plusplus.1/Linux-AlmaLinux_9.5-gcc_11-x86_64/bin/ND280GEANT4SIM.exe \
    /usr/local/t2k/current/nd280Geant4Sim_7.17-plusplus.1/Linux-AlmaLinux_9.5-gcc_11-x86_64/lib/libnd280Geant4Sim.so
echo "Created: $output"
