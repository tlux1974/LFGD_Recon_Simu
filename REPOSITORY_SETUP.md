# Obtaining the LFGD/HFGD development repositories

This manual reproduces the coordinated local software used by the LFGD
validation workflow. It assumes access to the T2K GitLab repositories and the
ND280++ Singularity container.

## 1. Access and directory layout

Check SSH access before cloning:

```bash
ssh -T git@git.t2k.org
ssh -T git@github.com
```

The T2K command requires a GitLab account with ND280 repository access. If it
reports `Permission denied (publickey)`, register the correct SSH public key.

Keep this relative layout; the absolute workspace path may differ:

```text
ND280++/
  LFGD_Recon_Simu/
  SoftProj/
    oaEvent/
    oaGeomInfo/
    nd280Geant4Sim/
    detResponseSim/
    hfgrecon/
  ND280ppCont
```

`ND280ppCont` is the existing ND280++ Singularity image or link and is not a
Git repository cloned below.

## 2. Clone the correct branches

```bash
mkdir -p "$HOME/HK/ND280++/SoftProj"
cd "$HOME/HK/ND280++"

git clone --branch main --single-branch \
  git@github.com:tlux1974/ND280pp_Tools.git LFGD_Recon_Simu

git clone --branch LFGD_tlux --single-branch \
  git@git.t2k.org:nd280/base/oaEvent.git SoftProj/oaEvent

git clone --branch LFGD_tlux --single-branch \
  git@git.t2k.org:nd280/base/oaGeomInfo.git SoftProj/oaGeomInfo

git clone --branch LFGD_tlux --single-branch \
  git@git.t2k.org:nd280/sim/nd280Geant4Sim.git SoftProj/nd280Geant4Sim

git clone --branch LFGD_tlux --single-branch \
  git@git.t2k.org:nd280/sim/detResponseSim.git SoftProj/detResponseSim

git clone --branch LFGDReco_tlux --single-branch \
  git@git.t2k.org:nd280/recon/hfgrecon.git SoftProj/hfgrecon
```

Do not substitute package `master` branches. The changes span geometry IDs,
geometry lookup, simulation, detector response, and reconstruction. Mixing
feature and standard package versions is unsupported.

## 3. Verify branches and revisions

The coordinated snapshot verified on 5 August 2026 is:

| Repository | Branch | Verified commit |
|---|---|---|
| `LFGD_Recon_Simu` | `main` | `3415a9f73add7aa474ab1f023c7888f15b68707c` |
| `oaEvent` | `LFGD_tlux` | `fe07b585160b4217e00a96e25110938b987de8dc` |
| `oaGeomInfo` | `LFGD_tlux` | `6a162e025fe16041c7b6f864b774b4c5900c6452` |
| `nd280Geant4Sim` | `LFGD_tlux` | `4aa1f881a308c31daadf1c268a16cb71b201ed40` |
| `detResponseSim` | `LFGD_tlux` | `2677c318b88146814b2cb95e14492bc66428` |
| `hfgrecon` | `LFGDReco_tlux` | `ea64bc5cb8c40436a28debce4e3c8d7388895829` |

Normally use the latest tips of these feature branches. Use the hashes only
to freeze an exact historical comparison. Print the active snapshot with:

```bash
cd "$HOME/HK/ND280++"
for repo in LFGD_Recon_Simu SoftProj/oaEvent SoftProj/oaGeomInfo \
            SoftProj/nd280Geant4Sim SoftProj/detResponseSim \
            SoftProj/hfgrecon; do
  printf '\n%-35s branch=%-20s commit=' "$repo" \
    "$(git -C "$repo" branch --show-current)"
  git -C "$repo" rev-parse --short HEAD
done
```

Update later without creating accidental merge commits:

```bash
for repo in LFGD_Recon_Simu SoftProj/oaEvent SoftProj/oaGeomInfo \
            SoftProj/nd280Geant4Sim SoftProj/detResponseSim \
            SoftProj/hfgrecon; do
  git -C "$repo" pull --ff-only
done
```

## 4. Enter and initialize the container

```bash
cd "$HOME/HK/ND280++"
singularity shell ND280ppCont
```

Inside the container:

```bash
export ND280_SYSTEM=Linux-AlmaLinux_9.5-gcc_11-x86_64
export PATH="/usr/local/t2k/current/nd280SoftwarePilot/scripts:/usr/local/t2k/current/nd280SoftwarePilot/${ND280_SYSTEM}/bin:${PATH}"
set +u
source /usr/local/t2k/current/nd280SoftwareMaster_14.36-plusplus.0.3/bin/setup.sh
set -u
```

The container/master version is part of the validated snapshot. Record any
version change before comparing samples made in different environments.

## 5. Build in dependency order

Run the following inside the container. It builds a package and then selects
that local build while configuring its downstream dependants.

```bash
cd "$HOME/HK/ND280++"
export LFGD_WORKSPACE="$PWD"
export JOBS="${JOBS:-4}"

select_package() {
  local root_name="$1" package="$2"
  local root="$LFGD_WORKSPACE/SoftProj/$package"
  export "${root_name}ROOT=$root"
  export "${root_name}CONFIG=$ND280_SYSTEM"
  export PATH="$root/$ND280_SYSTEM/bin:$PATH"
  export LD_LIBRARY_PATH="$root/$ND280_SYSTEM/lib:${LD_LIBRARY_PATH:-}"
}

build_package() {
  local package="$1"
  cmake -S "$LFGD_WORKSPACE/SoftProj/$package/cmake" \
        -B "$LFGD_WORKSPACE/SoftProj/$package/$ND280_SYSTEM"
  cmake --build "$LFGD_WORKSPACE/SoftProj/$package/$ND280_SYSTEM" \
        -j"$JOBS"
}

build_package oaEvent
select_package OAEVENT oaEvent

build_package oaGeomInfo
select_package OAGEOMINFO oaGeomInfo

build_package nd280Geant4Sim
select_package ND280GEANT4SIM nd280Geant4Sim

build_package detResponseSim
select_package DETRESPONSESIM detResponseSim

build_package hfgrecon
select_package HFGRECON hfgrecon

cd "$LFGD_WORKSPACE/LFGD_Recon_Simu"
./build_flat_treemaker.sh
```

Warnings about package directory names, Git version strings, modified
installed dependencies, or missing optional PYTORCH support are expected in
the present environment. A linker error or missing executable is not.

## 6. Verify and run

Check that local executables take precedence:

```bash
command -v ND280GEANT4SIM.exe
command -v DETRESPONSESIM.exe
command -v HFGRECON.exe
command -v LFGDFLATTREE.exe
```

The first three paths must be below `SoftProj`; the last must be below
`LFGD_Recon_Simu`. Then run the three comparable configurations:

```bash
cd "$HOME/HK/ND280++/LFGD_Recon_Simu"
./run_student_sample.sh hfg-standard 100 test01
./run_student_sample.sh lfg-original-low 100 test01
./run_student_sample.sh lfg-best 100 test01
```

The runner selects all five local packages directly; it does not require
private switch scripts outside the Git checkout. See `STUDENT_GUIDE.md` for
the configuration meanings, output files, plots, and systematic tests.

## 7. Common problems

- `could not load cache`: configure first with `cmake -S ... -B ...` and pass
  the same build directory to `cmake --build`.
- An executable resolves below `/usr/local/t2k/current`: the local build was
  not selected or did not finish successfully.
- `Remote branch ... not found`: run `git fetch origin` and check access and
  spelling. Only hfgrecon uses `LFGDReco_tlux`; the other ND280 packages use
  `LFGD_tlux`.
- Detached HEAD after testing a commit hash: return with `git switch` followed
  by the branch name in the table.
- Results differ between users: record all six commit hashes, container
  version, parameter file, random seed, particle settings, and event count.
