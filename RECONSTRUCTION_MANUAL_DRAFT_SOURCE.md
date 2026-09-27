# Source material for the complete LFGD/HFGD simulation and reconstruction manual

This is a project-wide document, not material specific to the proton-hit or
matched-muon studies.  It is deliberately an inventory rather than a polished
user manual.  It records the available programs, switches, reconstruction
modes, parameter interactions, outputs, and known limitations.  It can be
supplied to ChatGPT together with the repository when drafting the complete
manual.

## 1. Scope and reconstruction chain

The local study runs this chain:

1. `generate_primary_events.py` optionally creates a detector-independent CSV
   of primary particles.
2. `generate_gps_macro.py` converts either that CSV or a statistical source
   description into an ND280 Geant4 macro.
3. `ND280GEANT4SIM.exe` transports the primary particles through either the
   HFGD or HOMO/LFGD geometry.
4. `DETRESPONSESIM.exe` converts Geant4 deposits into fibre/MPPC response.  For
   HOMO this includes a selectable ROOT light map.
5. `HFGRECON.exe` constructs 2D fibre clusters, 3D hits, clusters, spanning
   trees and fitted tracks.
6. `LFGDFLATTREE.exe` writes diagnostic flat trees.
7. `plot_overlay.py` produces ROOT-based control plots and event displays.

The HFGD and LFGD currently share most of the downstream reconstruction.  The
important study question is which parts of an algorithm designed for physical
scintillator cubes remain meaningful when LFGD light is distributed over
multiple fibres and neighbouring virtual columns.

## 2. Main entry points

### `run_student_sample.sh`

Usage:

```bash
LFGD_Recon_Simu/run_student_sample.sh CONFIG EVENTS [RUN_NAME]
```

It creates a new output directory and refuses to overwrite an existing run.
It performs the complete chain from primary generation through plots.

Recognized configurations:

- `hfg-standard`: physical HFGD cubes and standard HFGD reconstruction.
- `lfg-original-low`: historical HOMO peak finder, with 10 PE thresholds.
- `lfg-best`: local-direction 2D clustering and view-average charge.
- `lfg-default-position`: current local 2D clustering, standard charge
  sharing and continuous reconstructed position; used as the control.
- `lfg-dominant-axis-position`: same chain, but the coordinate along the
  dominant axis is placed at the slice centre and only the two perpendicular
  views determine transverse position.

Relevant environment variables:

- `PARTICLE`, default `mu-`.
- `ENERGY_MEV`, default 700.
- `SEED`, default 12345.
- `DIRECTION_MODE`: `fixed`, `isotropic`, or `cone`.
- `DIRECTION`: three-vector used by fixed and cone modes.
- `CONE_HALF_ANGLE_DEG`, default 5 degrees.
- `POSITION_MM`, default `0 0 1800`.
- `POSITION_FRAME`: `plusplus` or `global`.  In the PlusPlus frame the default
  centre maps to global `(0,30,910)` mm.
- `PRIMARY_INPUT_FILE`: explicit event CSV.  This is required for exact
  event-by-event comparisons between detector configurations.
- `REPLAY_PRIMARY_EVENTS=1`, default: explicit primaries are replayed one by
  one.  Setting it to zero uses one Geant4 run but no longer guarantees equal
  later directions in different detector geometries.
- `HOMO_LIGHTMAP_FILE`: override the installed HOMO light map.  It is invalid
  for HFGD.
- `DETRESPONSE_PARAMETER_FILE`: additional detector-response parameters.
- `DISABLE_HOMO_ATTENUATION=1`: coordinated validation switch disabling HOMO
  fibre attenuation in detector response and reconstruction.
- `FAST_HOMO_RESPONSE=1`: aggregate photons per fibre to reduce memory use.

### `rerun_reconstruction.sh`

Usage:

```bash
LFGD_Recon_Simu/rerun_reconstruction.sh DETRESPONSE_ROOT OUTPUT_DIR [HFGRECON_OPTIONS...]
```

This reuses an existing detector-response file.  It reruns only HFGRECON, the
flat-tree maker and plotting.  A reconstruction override is supplied as:

```bash
-O "par_override=/absolute/path/to/parameters.dat"
```

This is the preferred way to compare hit-building and fitting options because
the Geant4 deposits and fibre response remain identical.

## 3. Exact primary-event replay

`generate_primary_events.py` writes these columns:

`event, particle, pdg, kinetic_energy_mev, x_mm, y_mm, z_mm, dx, dy, dz`.

The CSV contains global coordinates after conversion from the requested
frame.  It fixes the primary state, not the random numbers used during later
particle transport.  Therefore HFGD and LFGD have the same initial MC truth,
while their subsequent interactions can legitimately differ.

A common Geant4 seed without a CSV is not sufficient.  Different detector
geometries consume different numbers of random values, so subsequent
isotropic directions diverge after the first event.

## 4. Detector-response options

`run_pipeline.sh` selects `upgrade-nd280plus` for HFGD and
`upgrade-nd280plus-homo` for LFGD.  Unrelated detector responses are disabled
for this focused comparison.

For HOMO, `HOMO_LIGHTMAP_FILE` creates a runtime override for
`detResponseSim.LiquidO.Response.File`.  If omitted, the standard light map
installed with the container is used.  The present matched-muon test uses only
that standard 1 mm-fibre light map.

Detector-response output is `detresponse.root`.  Reconstruction comparisons
must reuse this file rather than rerun response simulation.

## 5. HFGRECON 2D and 3D hit options

Parameters without `.homo` apply to HFGD.  HOMO/LFGD-specific overrides use
the `.homo` suffix.

### Timing and initial fibre selection

- `hfgRecon.TimeSlice.GapCut`: split fibre signals at large time gaps;
  standard value 50 ns.
- `hfgRecon.TimeSlice.ChargeThreshold`: initial calibrated fibre threshold;
  standard value 2.5 PE.
- `hfgRecon.Hits3D.Min2DHitCharge.homo`: minimum HOMO 2D fibre/cluster charge;
  current study value 10 PE.
- `hfgRecon.Hits3D.Allow2dHits[.homo]`: allow 3D candidates with only two
  measured views.
- `hfgRecon.Hits3D.AllowedDeadFibers`: allowed absent/dead fibres.

### HOMO 2D peak clustering

`hfgRecon.Hits3D.FindPeaks.homo` enables the peak stage.

`hfgRecon.Hits3D.PeakClusteringMode.homo` selects:

- `0`: historical peak finder.  Its stored charge is the peak-fibre charge.
- `1`: experimental local-direction clustering.  Nearby peaks estimate a
  local tangent; fibres are assigned in longitudinal/transverse windows and
  each fibre is owned by at most one 2D cluster.  Cluster charge is the sum of
  its fibres.

Mode-0 option:

- `OriginalPeakNeighbourDistance.homo`: peak comparison distance in fibre
  cells.

Mode-1 options:

- `LocalPeakNeighbourDistance.homo`.
- `LocalTangentPeakCount.homo`.
- `LocalMaximumTangentDistance.homo`.
- `LocalLongitudinalWindow.homo`.
- `LocalTransverseWindow.homo`.
- `LocalLongitudinalWeight.homo`.

`UsePeakWeightedPosition.homo=1` uses charge-weighted 2D positions.
`PositionMatchTolerance.homo` controls agreement when clusters from the three
views are combined into 3D candidates.

### Charge sharing

`hfgRecon.Hits3D.SharingAlgo` controls HFGD and
`hfgRecon.Hits3D.SharingAlgo.homo` controls HOMO:

- `1`: `OptimizeCubes`, the normal/default HFGD/SFGD charge-sharing fit.
- `2`: constrained charge sharing.
- `3`: maximum-entropy method; slow and numerically problematic.
- Any other value, conventionally `0`: disable charge sharing.

The default algorithm assumes cube deposits whose signals are shared by
fibres.  That model is not automatically correct for LFGD virtual columns.
Turning it off is useful as a topology diagnostic, but it changes reconstructed
charges and therefore can also change selection, clustering and tracks.

`ConserveCharge` constrains total normalization.  HOMO attenuation correction
is controlled separately by `ApplyAttenuationCorrection.homo` and should only
be disabled when response simulation also omitted attenuation.

### HOMO charge modes

`hfgRecon.Hits3D.ChargeReconstructionMode.homo` selects:

- `0`: original SFGD-style sharing using every available view.
- `1`: direction-aware selection of views followed by the existing sharing
  algorithm; requires sharing algorithm 1 or 2.
- `2`: direction-aware arithmetic mean of accepted view-cluster charges,
  without the sharing fit.
- `3`: mode 2 plus strict ownership of localizing 2D hits.  This is diagnostic
  because a projected cell can legitimately be shared by successive virtual
  voxels.

Direction-aware charge options include `ChargeDirectionNeighbours.homo`,
`ChargeMinimumDirectionNeighbours.homo`,
`ChargeDirectionMaximumDistance.homo`,
`ChargeViewMaximumParallelCosine.homo`, and `ChargeMinimumViews.homo`.

### HOMO position modes

`hfgRecon.Hits3D.VoxelPositionMode.homo` selects:

- `0`: historical continuous charge-weighted position.
- `1`: snap to the nearest 10 mm virtual-cube centre, assign the cube ID and
  retain one candidate per voxel.
- `2`: assign the virtual-cube ID and remove duplicate voxel candidates, but
  retain the chosen continuous position.  This is the current control.
- `3`: estimate a local direction, select its dominant detector axis, fix that
  coordinate at the 10 mm slice centre, and use the two perpendicular views
  for the transverse coordinates.  It falls back to mode 2 if the direction
  is inadequate or farther than 45 degrees from every axis.
- `4`: experimental event-level column collapse used in the matched-muon
  test.  A charge-weighted PCA of preliminary candidates selects the dominant
  detector axis.  All candidates in each 10 mm slice are collapsed into one
  charge-weighted point, with the dominant coordinate at the slice centre.
  Downstream official selection and fitting remain unchanged.  This is a
  diagnostic prototype, not yet a direction-free local reconstruction.

Mode-3 controls are `PositionDirectionNeighbours.homo`,
`PositionMinimumDirectionNeighbours.homo`,
`PositionDirectionMaximumDistance.homo`, and
`PositionDominantAxisMinimumCosine.homo`.

## 6. Hit selection and official track reconstruction

`SelectHits3D.Enable` is 2 for standard HFGD.  Mode 2 applies charge and
neighbour requirements.  HOMO normally uses `Enable.homo=1`, a direct charge
cut.  Setting an unrecognized/disabled algorithm copies input hits.

Relevant selection controls:

- `SelectHits3D.ChargeCut[.homo]`.
- `SelectHits3D.ChargeCut2`.
- `SelectHits3D.NeighborCut`.
- `SelectHits3D.UseValidFibers`.

Downstream stages are:

1. `THFGClusterHits` with `ClusterHits.HomoNeighbourhood` for HOMO.
2. `THFGSpanningTree`; `DistanceType` selects edge weighting and
   `GhostHitCut` removes unphysical tiny-charge hits.
3. `THFGFindKinks`; controls include `ScanLength`, `KinkThreshold`, and
   `LengthFraction`.
4. Cluster growth; important controls include `MaxLineHits`,
   `Chi2Threshold`, and `MinChargePerHit`.
5. Track growth/merging; controls include `MergeDistance`,
   `MinimumDistance`, `AllowedKink`, `GoodnessCut`, `MatchedAngleCut`, and
   `MatchedPositionCut`.
6. `THFGMergeXTalk`, whose DBSCAN distance is controlled by
   `MergeXTalk.MaxClusterDist[.homo]`.
7. PID and vertex creation.

The stochastic official track fitter is used after these stages.  Track
multiplicity therefore depends not only on the 3D positions, but also on hit
charges, thresholds, graph connectivity, kink splitting and merging.

## 7. Parameter sets used in current studies

- `default_position.parameters.dat`: local 2D clustering, sharing enabled,
  position mode 2.
- `dominant_axis_position.parameters.dat`: position mode 3.
- `column_official.parameters.dat`: position mode 4 followed by the official
  chain.
- `default_position_no_charge_sharing.parameters.dat`: standard LFGD control
  with HOMO sharing disabled.
- `column_official_no_charge_sharing.parameters.dat`: column mode with HOMO
  sharing disabled.
- `Studies/HFGDHitRecoStudy/no_charge_sharing.parameters.dat`: standard HFGD
  reconstruction with sharing disabled.

## 8. Flat-tree output and diagnostics

Important trees in `flat.root`:

- `fiber_hits`: detector-response fibre signals.
- `homo_raw`: raw HOMO-specific response information.
- `homo_truth`: truth information associated with HOMO response.
- `hits2d`: both used and unused 2D candidates, including rejection codes and
  reasons where available.
- `hits3d`: reconstructed 3D hits.
- `hit3d_views`: views contributing to each 3D hit.
- `track_nodes`: nodes from the official fitted tracks.
- `track_node_hits`: hit association for track nodes.
- `mc_track_points`: MC trajectory points.
- `mc_virtual_segments`: truth segments through virtual voxels.

`plot_overlay.py` creates charge, fibre-position, 3D-hit, view-composition,
truth-matching and per-event projection plots.  Its HTML index is written to
`plots/index.html`.

## 9. Current matched-muon observations

For ten identical 1 GeV isotropic muon primaries:

| reconstruction | charge sharing | 2D hits | 3D hits | official tracks |
|---|---:|---:|---:|---:|
| HFGD standard | enabled | — | — | 20 |
| LFGD standard | enabled | 2638 | 2455 | 95 |
| LFGD column mode | enabled | 2638 | 870 | 40 |
| HFGD standard | disabled | — | 1571 | 30 |
| LFGD standard | disabled | 2638 | — | 123 |
| LFGD column mode | disabled | 2638 | 870 | 40 |

These are ten-event diagnostic numbers, not performance measurements.  Delta
electrons and true secondary branches have not been classified.  No tracking
parameters have been optimized for the LFGD column representation.  The
column mode nevertheless demonstrates that duplicated light-sharing 3D
candidates are a major source of LFGD track fragmentation.

## 10. Known limitations and next decisions

- Mode 4 uses an event-level PCA to choose an axis.  A future reconstruction
  must work locally and should not require a completed track seed.
- The current column collapse uses preliminary 3D candidates rather than a
  direct two-view likelihood within each slice.
- Its collapsed charge is a temporary definition, not an LFGD-specific
  physical charge-sharing result.
- Delta electrons, forks and separated collinear pieces must be classified
  before tuning merging aggressively.
- A future local method should identify a cluster, choose a local dominant
  axis, define pitch-wide perpendicular slices, reconstruct transverse CoG or
  a restricted light-map likelihood from the two perpendicular views, and
  then hand those points to the official fitter.
- Track-fitter parameters should be scanned only after hit charge and point
  construction are sufficiently stable.
