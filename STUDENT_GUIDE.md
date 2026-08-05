# Student guide: HFGD/LFGD simulation and reconstruction

This guide assumes no detailed knowledge of the ND280++ reconstruction.  Its
purpose is to produce comparable samples without editing commands or parameter
files by hand.

## What was changed, in plain language

The HFGD is made from physical scintillator cubes.  The LFGD is homogeneous,
but its signals still have to be addressed on a virtual 10 mm cube/fibre grid
so that the same reconstruction framework can be tested.  The coordinated
local packages now provide that grid from simulation through reconstruction:

1. `nd280Geant4Sim` builds the compact homogeneous detector and records truth
   segments split at virtual-cube boundaries.
2. `oaEvent` defines deterministic virtual cube and fibre identifiers.
3. `oaGeomInfo` calculates virtual cube centres, fibre positions, endpoints,
   and transit distances analytically.
4. `detResponseSim` finds the virtual LFGD fibres without requiring physical
   fibre volumes in the ROOT geometry.
5. `hfgRecon` reconstructs LFGD 2D fibre clusters, 3D hits, clusters, tracks,
   and track nodes. Several LFGD algorithms remain selectable for comparison.
6. `LFGDFLATTREE.exe` writes simple TTrees containing truth, fibre hits, 3D
   hits, their view composition, reconstructed tracks, and exact node-to-hit
   associations. `plot_overlay.py` makes browsable event displays and summary
   plots from this flat file.

None of these local changes is needed to change the standard HFGD algorithm;
the HFGD sample is the reference configuration.

## The three supported comparisons

### `hfg-standard`

The physical HFGD and its standard reconstruction parameters. Use this as the
reference sample.

### `lfg-original-low`

The historical LFGD peak-fibre method, with both relevant LFGD thresholds set
to 10 PE. The lower post-3D threshold was important: the older high threshold
discarded many LFGD hits before track building.

Important settings are:

```text
PeakClusteringMode.homo = 0
Min2DHitCharge.homo = 10
ChargeReconstructionMode.homo = 0
SelectHits3D.ChargeCut.homo = 10
VoxelPositionMode.homo = 2
```

### `lfg-best`

The best-looking experimental LFGD configuration so far:

- local 2D clustering with unique physical-fibre ownership;
- charge-weighted 2D positions;
- one reconstructed candidate per virtual cube ID while retaining the
  continuous reconstructed position;
- a 3D charge equal to the mean of the accepted 2D-view charges;
- a locally fibre-parallel view may be omitted, leaving a two-view mean;
- 10 PE 2D and selected-3D thresholds.

This mode produces continuous tracks and a nearly Gaussian 3D-hit charge peak
near 100 PE in the current 700 MeV muon sample. It is still experimental:
the same projected 2D cluster can be compatible with several 3D candidates,
so the event-level sum of 3D charge is not yet guaranteed to be less than the
fibre-charge sum. Do not use the current total 3D charge as a final calibrated
dE/dx measurement.

The rejected strict-ownership experiment is retained internally as charge mode
3 for debugging only. Do not use it for production comparisons; it removed
valid repeated projected cells and severely reduced efficiency.

## Build once

Enter the ND280++ Singularity container and run the normal environment setup.
The coordinated local packages should already be built. Build the flat-tree
maker once, or whenever `src/LFGDFlatTree.cxx` changes:

```bash
cd /home/tlux/HK/ND280++/LFGD_Recon_Simu
./build_flat_treemaker.sh
```

The student run script checks that every required executable is available and
stops with an explanation if a local package has not been built.

## Run a sample

From inside the container:

```bash
cd /home/tlux/HK/ND280++/LFGD_Recon_Simu

./run_student_sample.sh hfg-standard 100 test01
./run_student_sample.sh lfg-original-low 100 test01
./run_student_sample.sh lfg-best 100 test01
```

The three commands use the same default particle settings and random seed, but
write separate directories. `100` is the number of events and `test01` is an
arbitrary label. Existing output directories are never overwritten.

Defaults are:

```text
particle:       mu-
kinetic energy: 700 MeV
vertex:         (0,0,1800) mm in PlusPlusTracker coordinates
directions:     isotropic
random seed:    12345
```

Examples of controlled variations:

```bash
SEED=6789 ./run_student_sample.sh lfg-best 100 seed6789

ENERGY_MEV=1000 ./run_student_sample.sh lfg-best 100 mu1000

PARTICLE=pi+ ENERGY_MEV=700 \
  ./run_student_sample.sh lfg-best 100 pion700

DIRECTION_MODE=fixed DIRECTION="0 0 1" \
  ./run_student_sample.sh lfg-best 100 fixed_z
```

For a fair comparison, use the same event count, particle, energy, direction
mode, vertex, and seed for all three configurations.

## What the script runs

The five stages are:

```text
Geant4 -> detector response -> 3D reconstruction/tracking -> flat tree -> plots
```

Every output directory contains:

- `gps.mac`: generated particle-gun configuration;
- `g4.root`: Geant4 truth and detector energy deposits;
- `detresponse.root`: simulated fibre/electronics response;
- `reco.root`: full ND280 event containing 3D hits and reconstructed objects;
- `flat.root`: analysis-friendly diagnostic trees;
- `01_geant4.log` through `04_flat_tree.log`;
- `plots/index.html`: clickable displays for the first ten events and summary
  plots over those events.

All requested events are stored in the ROOT files. Only the first ten are
plotted automatically to avoid creating thousands of PNG files.

## Flat-tree contents

`flat.root` contains:

| Tree | Meaning |
|---|---|
| `fiber_hits` | Measured fibre hits, position, charge, projection and indices |
| `hits3d` | Reconstructed 3D-hit position, charge and virtual/physical cube ID |
| `hit3d_views` | View charge, retained fibre composition and view reuse for each 3D hit |
| `track_nodes` | Fitted track-node positions used to draw reconstructed tracks |
| `track_node_hits` | Exact 3D hits attached to each track node |
| `mc_virtual_segments` | Truth entry/exit segment, energy and length in each LFGD virtual cube or HFGD physical cube |
| `mc_track_points` | Saved Geant4 trajectory points for additional debugging |

Use `track_node_hits`, not the fitted `track_nodes` position, to decide whether
a reconstructed 3D hit was used by a track. A fit is allowed to move a node
position into a neighbouring voxel.

Despite its historical name, `mc_virtual_segments` is filled for both detector
configurations. For LFGD, `cube_x/y/z` identify the 10 mm virtual cube crossed
by the Geant4 segment. For HFGD they are the indices of the actual physical
cube. The `detector` branch is `0` for LFGD/HOMO and `1` for HFGD. The HFGD
mapping is taken from `hfg_truth`, which stores the exact cube-to-Geant4-
segment association; it is not obtained by rounding the segment position.

To calculate a position residual per cube, group rows by
`(event,cube_x,cube_y,cube_z)`, calculate the energy-weighted segment midpoint
from `start_*`, `stop_*`, and `energy_deposit`, and compare it with the
continuous `x/y/z` in `hits3d` having the same cube ID. More than one Geant4
segment can legitimately occur in a cube, for example after a hard
interaction. The plotting script already performs this aggregation and writes
the residual tables.

## Looking at more events without rerunning reconstruction

```bash
python3 plot_overlay.py output/student_lfg-best_test01/flat.root \
  --event-range 10 29 \
  --output-dir output/student_lfg-best_test01/plots_10_29
```

This reads only `flat.root`; it does not rerun simulation or reconstruction.

## Rerunning only reconstruction

To compare a reconstruction parameter change using an existing detector-
response file:

```bash
./rerun_reconstruction.sh \
  output/student_lfg-best_test01/detresponse.root \
  output/student_lfg-best_test01/reco_alternative \
  -O par_override=/absolute/path/to/parameters.dat
```

Never write the alternative into the original directory. Separate directories
make comparisons reproducible.

## Quantities to compare systematically

At minimum record, per configuration and particle direction:

- number of fibre hits, 3D hits, tracks and track nodes;
- MC virtual voxels crossed;
- matched, missing and off-track reconstructed voxels;
- 3D hits represented by track nodes and those left unused;
- 3D-hit charge distribution for used and unused hits;
- 2D-view charge and view-reuse multiplicity;
- sum of fibre charge, sum of 3D charge and their ratio;
- MC/reconstructed position residuals per virtual cube;
- visible track continuity in all three projections.

The CSV files in `plots/` allow these quantities to be combined across many
jobs without parsing PNG files.

## Known limitations

- Mode 2 improves the individual-hit charge distribution but does not yet
  solve repeated use of a projected 2D measurement by several 3D candidates.
- Strictly assigning every non-parallel 2D cluster to one voxel was tested and
  failed because finite 10 mm cells can legitimately repeat in projection.
- A future charge algorithm should use the other two projections to infer how
  many voxel crossings share an extended view and divide that view charge.
- The flat-tree and plots are validation tools, not an official analysis data
  format.

## Where detailed package notes live

- `SoftProj/nd280Geant4Sim/README_homo_virtual_cube_truth.md`
- `SoftProj/oaEvent/README_homo_virtual_geometry.md`
- `SoftProj/oaGeomInfo/README_homo_virtual_geometry.md`
- `SoftProj/detResponseSim/README_homo_virtual_geometry.md`
- `SoftProj/hfgrecon/README_local_direction_clustering.md`
- `SoftProj/hfgrecon/README_3d_hit_position_charge_tests.md`
