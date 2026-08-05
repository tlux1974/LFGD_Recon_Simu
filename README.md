# LFGD reconstruction and simulation

This directory provides a reproducible diagnostic chain for fixed 700 MeV
muons originating inside either detector design.

For the supported three-way student comparison and a plain-language summary
of all coordinated local changes, start with [STUDENT_GUIDE.md](STUDENT_GUIDE.md)
and `run_student_sample.sh`.

## Setup

```bash
cd ~/HK/ND280++/LFGD_Recon_Simu
source <ND280 setup>
source ~/HK/ND280++/switch-nd280geant4sim.sh local
source ~/HK/ND280++/switch-hfgrecon.sh local
./build_flat_treemaker.sh
export PATH="$PWD/$ND280_SYSTEM/bin:$PATH"
```

Use `original` instead of `local` with either switch to select the installed
container package for a direct comparison.

The package-level `bin/setup.sh` is not sufficient in a fresh shell: the base
container/login environment must already provide `nd280-system` and
`CMAKE_PREFIX_PATH`. The build script checks both and stops with an explanation
instead of producing a misleading partial configuration.

The tree maker normally only needs rebuilding after changing
`src/LFGDFlatTree.cxx` or its dependencies.

## Select the detector

The detector-switch script must be sourced because it exports `DETECTOR` in
the current shell:

```bash
source ./switch-detector.sh lfgd
source ./switch-detector.sh status
```

For the HFGD comparison, use:

```bash
source ./switch-detector.sh hfgd
source ./switch-detector.sh status
```

The internal ND280 model names are:

| Name used in this study | ND280 model name |
|---|---|
| LFGD | `homo` |
| HFGD | `hfg` |

Both designs are generated at the same detector location. The generated
Geant4 macro enables the selected design and disables the alternative.

## Run the complete chain

After selecting a detector:

```bash
./run_pipeline.sh
```

You can also specify the detector and number of events directly:

```bash
./run_pipeline.sh homo 10   # LFGD (the code calls this model "homo")
./run_pipeline.sh hfg 10    # HFGD comparison
```

The arguments are:

```text
./run_pipeline.sh <detector-model> <number-of-events>
```

Therefore, the `10` means **simulate ten events**. With the current particle
gun, each event contains one independent 700 MeV muon starting at the same
position and travelling in the same initial direction. For a quick check use
one event:

```bash
./run_pipeline.sh homo 1
```

For a larger sample of 1000 muons:

```bash
./run_pipeline.sh homo 1000
```

The command-line values override the detector selected with
`switch-detector.sh` and the event count in `config.sh` for that run only.

## Simulation settings

The default gun is `mu-`, 700 MeV kinetic energy, starting at
`(0,0,1800) mm` in fixed direction `(0,0,1)`. Edit `config.sh` or override
variables:

```bash
POSITION_FRAME=plusplus POSITION_MM="0 0 1800" \
    DIRECTION="1 0 1" ./run_pipeline.sh homo 20
```

For a full-sphere isotropic direction distribution from the same point:

```bash
DIRECTION_MODE=isotropic ./run_pipeline.sh homo 1000
```

For an explicit fixed direction:

```bash
DIRECTION_MODE=fixed DIRECTION="1 0 1" ./run_pipeline.sh homo 1000
```

`SEED=12345` is the default for both Geant4 and detector response, so repeated
runs with identical settings reproduce the same isotropic directions and
detector fluctuations. Override it explicitly with, for example, `SEED=6789`.
The default output names are `homo_mu700_isotropic` and `homo_mu700_fixed`, so
the two direction modes do not overwrite each other.

In this example, `20` means twenty events. Environment overrides apply only to
that command unless they are exported beforehand.

`baseline-2024-plusplus` places both alternative detector designs between
z=1200 and 2400 mm. The generated macro enables exactly the requested design.

## Products

Each run writes to `output/<detector>_mu700_fixed/`:

- `g4.root`: Geant4 truth and energy deposits;
- `detresponse.root`: calibrated fibre/MPPC hits;
- `reco.root`: `hfg_3d` hits and `THFGRecon` tracks;
- `flat.root`: `fiber_hits`, `hits3d`, `hit3d_views`, `track_nodes`,
  `track_node_hits`, `mc_virtual_segments`, and `mc_track_points` TTrees;
- `plots/`: ZX, ZY, and XY overlays;
- numbered logs for each processing stage.

The reconstruction executable is named HFGRECON for both models. Internally it
selects `THomoHits3D` for the `homo` collection and `THFGHits3D` for `hfg`, then
runs the shared HFG clustering and tracking chain.

`mc_virtual_segments` keeps its old name for file compatibility, but is now a
common per-cube truth table. It stores LFGD Geant4 segments in 10 mm virtual
cubes and HFGD segments in their physical cubes. `detector=0` means LFGD/HOMO
and `detector=1` means HFGD. The HFGD cube ID is read from the exact
`hfg_truth` cube-to-segment relationship rather than reconstructed by rounding
coordinates. Consequently the same flat-tree and plotting workflow can be
used for voxel-by-voxel MC/reconstructed-hit position and charge comparisons
for either detector.

The detector-response command disables response simulation for unrelated
SciFi, SWD, SFG, ECAL, SMRD, TOF, TPC, and HAT subsystems.  The broad ND280+
presets enable these by default, but they are not inputs to this focused
Homo/HFG comparison and can otherwise produce unrelated response warnings.

`POSITION_FRAME` makes the input coordinate system explicit.  Its default is
`plusplus`, for which `POSITION_MM="0 0 1800"` is the LFGD/HFGD centre.  The
macro generator converts this to `(0,30,910) mm` in the global ND280 frame
used by Geant4 GPS, hits, and trajectories.  Alternatively,
`POSITION_FRAME=global POSITION_MM="0 30 910"` specifies that same vertex
directly.  The generated macro records both input and converted coordinates.
The detector extends about 1200 mm in z, from global 310 to 1510 mm.

The fibre positions in the overlay are reference points for extended fibres;
they are useful diagnostically but are not independent 3D particle positions.
In each two-dimensional overlay, a fibre hit is therefore shown only when
both displayed coordinates are transverse to that fibre.  Views involving
the arbitrary along-fibre reference coordinate omit that fibre projection.
The overlays also draw dashed projected envelopes for SWD, the selected
LFGD/HFGD detector, and SciFi in their global positions.

To scan a range of events and create a browsable `plots/index.html`:

```bash
python3 plot_overlay.py output/homo_mu700_fixed/flat.root \
    --output-dir output/homo_mu700_fixed/plots --event-range 0 99
```

Add `--only-with-data fibers`, `--only-with-data hits3d`, or
`--only-with-data tracks` to retain only events containing that layer.
The plotter also creates a fibre-hit charge distribution for every selected
event and an aggregate `fiber_charge_selected.png` for the complete selection.
Charge plots cover 0--100 in bins of width 4. Each event also gets a six-panel
`eventN_fiber_positions.png`: the two transverse position coordinates are
shown separately for each of the three fibre projections.
An aggregate `fiber_positions_selected.png` shows the same six distributions
over all events retained by the event selection and data filter.
The corresponding reconstructed-hit diagnostics are
`eventN_hit3d_charge.png` and `eventN_hit3d_positions.png` (separate x, y,
and z panels), with aggregate `hit3d_charge_selected.png` and
`hit3d_positions_selected.png` files for all selected events.
For charge-conservation checks, `charge_sums_selected.png` compares each
event's summed reconstructed 3D-hit charge with its summed original fibre-hit
charge and draws the equality line. `charge_sums.csv` records both sums, their
ratio, and whether the reconstructed sum exceeds the fibre sum. The HTML event
summary shows the same values.

`mc_track_points` contains the Geant4 truth path for every saved particle,
including `track_id`, `parent_id`, `pdg`, particle name, point number, position,
time, and momentum.  The overlay draws each MC trajectory independently in
green, so primary and secondary truth tracks are not accidentally connected.

## Typical session

```bash
cd ~/HK/ND280++/LFGD_Recon_Simu

# Enter/source the normal ND280++ environment first.
source ~/HK/ND280++/switch-hfgrecon.sh local

# Build once, then expose the diagnostic executable.
./build_flat_treemaker.sh
export PATH="$PWD/$ND280_SYSTEM/bin:$PATH"

# Simulate and reconstruct ten LFGD events.
source ./switch-detector.sh lfgd
./run_pipeline.sh

To rerun reconstruction without repeating Geant4 or detector response, give
the existing detector-response file and a separate output directory:

```bash
./rerun_reconstruction.sh \
    output/homo_mu700_isotropic/detresponse.root \
    output/homo_mu700_isotropic/reco_test
```

The new directory receives `reco.root`, `flat.root`, reconstruction logs, and
event displays. Additional arguments are forwarded to `HFGRECON.exe`, which
allows reconstruction parameter variations to be kept in separate folders.

# Run the same ten-muon configuration with HFGD.
source ./switch-detector.sh hfgd
./run_pipeline.sh
```

The LFGD and HFGD results go into separate output directories, so running the
comparison does not overwrite the first detector's files.
