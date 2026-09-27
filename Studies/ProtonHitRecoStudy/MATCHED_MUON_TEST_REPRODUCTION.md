# Reproduce the matched HFGD/LFGD 1 GeV muon test

This document focuses only on the ten-event comparison performed on
2026-09-26.  It is intended both as a runnable recipe and as source material
for a future manual.

## Purpose

Generate one detector-independent CSV containing ten isotropic 1 GeV muon
primaries, use it for separate HFGD and LFGD Geant4/response productions, and
then compare these reconstruction variants on unchanged detector-response
files:

1. standard HFGD, charge sharing enabled;
2. standard HFGD, charge sharing disabled;
3. standard LFGD, charge sharing enabled;
4. standard LFGD, charge sharing disabled;
5. column LFGD plus the official fitter, charge sharing enabled;
6. column LFGD plus the official fitter, charge sharing disabled.

The standard 1 mm-fibre HOMO light map supplied by the container is used.  Do
not set `HOMO_LIGHTMAP_FILE` for this test.  The magnetic field is disabled by
the generated GPS macro used by this study.

## 1. Enter and prepare the local ND280 environment

Run from `/home/tlux/HK/ND280++` inside `ND280ppCont`.  The normal student
runner selects the local coordinated packages automatically.  If hfgRecon was
modified, rebuild it first:

```bash
source /usr/local/t2k/current/nd280SoftwarePilot/nd280SoftwarePilot.profile
source /usr/local/t2k/current/nd280SoftwareMaster_14.36-plusplus.0.3/bin/setup.sh
export CMAKE_PREFIX_PATH="$ND280SOFTWAREPOLICYROOT:$CMAKE_PREFIX_PATH"
cmake --build /home/tlux/HK/ND280++/SoftProj/hfgrecon/Linux-AlmaLinux_9.5-gcc_11-x86_64 -j2
cd /home/tlux/HK/ND280++/LFGD_Recon_Simu
./build_flat_treemaker.sh
cd /home/tlux/HK/ND280++
```

The explicit `CMAKE_PREFIX_PATH` addition is needed because the current master
setup exposes `ND280SOFTWAREPOLICYROOT` but may leave the policy directory out
of CMake's prefix path.

## 2. Create the common primary CSV once

Choose a permanent input name rather than placing the only copy inside one
detector's output directory:

```bash
python3 LFGD_Recon_Simu/generate_primary_events.py --events 10 --seed 12345 --particle mu- --energy-mev 1000 --position-mm 0 0 1800 --position-frame plusplus --direction-mode isotropic --output LFGD_Recon_Simu/input/matched_muon1GeV_isotropic_seed12345_10.csv
```

Inspect it before starting the expensive simulations:

```bash
head LFGD_Recon_Simu/input/matched_muon1GeV_isotropic_seed12345_10.csv
```

The PlusPlus-frame source `(0,0,1800)` is written as global
`(0,30,910)` mm.  Every row stores the exact direction.  Reusing merely the
same random seed is not equivalent.

## 3. Produce the standard HFGD sample

```bash
PRIMARY_INPUT_FILE="$PWD/LFGD_Recon_Simu/input/matched_muon1GeV_isotropic_seed12345_10.csv" PARTICLE=mu- ENERGY_MEV=1000 DIRECTION_MODE=isotropic LFGD_Recon_Simu/run_student_sample.sh hfg-standard 10 muon1GeV_matched_reproduction
```

Expected directory:

```text
LFGD_Recon_Simu/output/hfg_student_hfg-standard_muon1GeV_matched_reproduction
```

It contains `gps.mac`, `g4.root`, `detresponse.root`, `reco.root`, `flat.root`,
the four numbered logs and `plots/index.html`.

## 4. Produce the standard LFGD sample from the same CSV

```bash
PRIMARY_INPUT_FILE="$PWD/LFGD_Recon_Simu/input/matched_muon1GeV_isotropic_seed12345_10.csv" PARTICLE=mu- ENERGY_MEV=1000 DIRECTION_MODE=isotropic LFGD_Recon_Simu/run_student_sample.sh lfg-default-position 10 muon1GeV_matched_reproduction
```

Expected directory:

```text
LFGD_Recon_Simu/output/homo_student_lfg-default-position_muon1GeV_matched_reproduction
```

The two Geant4 files are necessarily different because their detector
geometries differ, but event `N` starts with the same particle, energy, vertex
and direction in both.

## 5. Select the local reconstruction executable for reruns

The following commands assume the normal ND280 setup has already been
sourced:

```bash
cd /home/tlux/HK/ND280++/LFGD_Recon_Simu
source ./switch-hfgrecon.sh local
export PATH="$PWD/Linux-AlmaLinux_9.5-gcc_11-x86_64/bin:$PATH"
export LD_LIBRARY_PATH="$PWD/Linux-AlmaLinux_9.5-gcc_11-x86_64/lib:${LD_LIBRARY_PATH:-}"
```

## 6. Rerun HFGD with charge sharing disabled

The standard HFGD run from section 3 is the sharing-enabled result.  Produce
the no-sharing variant without rerunning Geant4 or detector response:

```bash
./rerun_reconstruction.sh output/hfg_student_hfg-standard_muon1GeV_matched_reproduction/detresponse.root output/hfg_muon1GeV_matched_reproduction_no_charge_sharing -O "par_override=$PWD/Studies/HFGDHitRecoStudy/no_charge_sharing.parameters.dat"
```

That parameter file sets `hfgRecon.Hits3D.SharingAlgo = 0` and changes no
other HFGD reconstruction option.

## 7. Rerun standard LFGD with charge sharing disabled

The standard LFGD run from section 4 is the sharing-enabled result.  Produce
the corresponding no-sharing control:

```bash
./rerun_reconstruction.sh output/homo_student_lfg-default-position_muon1GeV_matched_reproduction/detresponse.root output/homo_muon1GeV_matched_reproduction_standard_no_charge_sharing -O "par_override=$PWD/Studies/ProtonHitRecoStudy/default_position_no_charge_sharing.parameters.dat"
```

This retains local-direction 2D clustering and position mode 2, while setting
`hfgRecon.Hits3D.SharingAlgo.homo = 0`.

## 8. Produce column LFGD with the official fitter

With normal HOMO charge sharing enabled:

```bash
./rerun_reconstruction.sh output/homo_student_lfg-default-position_muon1GeV_matched_reproduction/detresponse.root output/homo_muon1GeV_matched_reproduction_column_official -O "par_override=$PWD/Studies/ProtonHitRecoStudy/column_official.parameters.dat"
```

With HOMO charge sharing disabled:

```bash
./rerun_reconstruction.sh output/homo_student_lfg-default-position_muon1GeV_matched_reproduction/detresponse.root output/homo_muon1GeV_matched_reproduction_column_official_no_charge_sharing -O "par_override=$PWD/Studies/ProtonHitRecoStudy/column_official_no_charge_sharing.parameters.dat"
```

Both use position mode 4 to collapse preliminary candidates to one point per
dominant-axis column.  The complete downstream official chain—selection,
clustering, spanning tree, kink finding, track growth, merging and stochastic
fitting—then runs normally.

## 9. Verify that all variants use the intended input

Check the override recorded by each reconstruction log:

```bash
grep -H "runtime parameter override file" output/*muon1GeV_matched_reproduction*/03_hfgrecon.log
```

Check that both initial productions embedded the same directions:

```bash
diff -u <(grep '^/gps/direction' output/hfg_student_hfg-standard_muon1GeV_matched_reproduction/gps.mac) <(grep '^/gps/direction' output/homo_student_lfg-default-position_muon1GeV_matched_reproduction/gps.mac)
```

No output from `diff` is expected.  The macros can differ in detector setup
commands and geometry.

Check whether charge sharing ran.  Sharing-enabled logs contain messages such
as `Augmented Cubes` and optimization summaries.  No-sharing logs should not
run those optimization steps for the relevant detector model.

## 10. Count reconstructed tracks consistently

The following sums the official track count reported once per event by the
vertex stage:

```bash
for d in output/hfg_student_hfg-standard_muon1GeV_matched_reproduction output/hfg_muon1GeV_matched_reproduction_no_charge_sharing output/homo_student_lfg-default-position_muon1GeV_matched_reproduction output/homo_muon1GeV_matched_reproduction_standard_no_charge_sharing output/homo_muon1GeV_matched_reproduction_column_official output/homo_muon1GeV_matched_reproduction_column_official_no_charge_sharing; do printf '%s ' "$d"; grep -E 'THFGPairwiseVertices:: [0-9]+ tracks with required length' "$d/03_hfgrecon.log" | awk '{sum += $4} END {print sum+0}'; done
```

For the original ten-event run performed today, the totals were:

| detector/method | sharing | official tracks | tracks/event |
|---|---:|---:|---:|
| HFGD standard | enabled | 20 | 2.0 |
| HFGD standard | disabled | 30 | 3.0 |
| LFGD standard | enabled | 95 | 9.5 |
| LFGD standard | disabled | 123 | 12.3 |
| LFGD column | enabled | 40 | 4.0 |
| LFGD column | disabled | 40 | 4.0 |

These counts include genuine delta-electron and secondary branches.  They are
not the number of falsely reconstructed tracks alone.

## 11. Event displays and flat trees

Each output has an HTML page:

```text
OUTPUT_DIRECTORY/plots/index.html
```

For detailed comparisons use `flat.root`.  In particular:

- `hits2d` contains used and unused 2D hits and rejection information;
- `hits3d` contains the reconstructed points;
- `track_nodes` and `track_node_hits` contain the official fitter output;
- `mc_track_points` and `mc_virtual_segments` contain the MC reference.

The existing matched-display utility is:

```bash
python3 /home/tlux/HK/ND280++/LFGD_Recon_Simu/plot_matched_event_displays.py
```

Its configured input paths may need to be changed to the new
`matched_reproduction` directories before execution.

## 12. Interpretation boundary

The test isolates the consequences of hit construction and charge sharing on
the same detector response.  It does not yet:

- classify primary-muon pieces versus delta electrons;
- optimize standard fitter parameters for LFGD column hits;
- implement a physical LFGD-specific column charge likelihood;
- prove that mode 4 is suitable without an event-level direction estimate.

The meaningful immediate conclusion is that collapsing duplicated LFGD 3D
candidates substantially reduces track fragmentation, while disabling the
standard charge-sharing fit worsens the standard HFGD and standard LFGD track
multiplicities in this small sample.

