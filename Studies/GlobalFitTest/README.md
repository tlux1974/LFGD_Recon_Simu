# Global HOMO light-distribution fit test

This is a standalone ROOT-based test. It does not read `mc_virtual_segments`
or any other MC track/segment information.

## Build

Enter the ND280 container, source the standard ND280 setup, and run:

```sh
cd LFGD_Recon_Simu/Studies/GlobalFitTest
make
```

## Produce the initial fit seed

```sh
./global_fit_seed flat.root EVENT=0 TREE=homo_truth
```

Optional arguments use the `KEY=value` format. Run `./global_fit_seed --help`
or run it without arguments for the complete list. For example:

```sh
./global_fit_seed flat.root EVENT=12 TREE=homo_raw OUTPUT_PREFIX=raw_seed_event_12
./global_fit_seed flat.root EVENT=12 TREE=fiber_hits OUTPUT_PREFIX=fibre_seed_event_12
./global_fit_seed flat.root EVENT=12 TREE=homo_truth OUTPUT_PREFIX=truth_seed_event_12 MIN_CHARGE=10
```

The same reader and seed calculation are used for `homo_truth`, `homo_raw`
and `fiber_hits`. All three trees expose the required `event`, `x`, `y`, `z`,
`charge`, and `projection` branches.

The start-point seed is local `(0,0,0)` mm, stored in the flat-tree global
frame as the HOMO centre `(0,30,910)` mm, and remains free transversely in the
likelihood fit. For each Cartesian coordinate, the direction seed uses the fibre hit with
the largest absolute displacement from the start point; its coordinate sign
sets the direction sign. The resulting three-vector is normalized.

A fibre does not constrain position along its own axis. Consequently, the X
distribution uses Y- and Z-directed fibres, Y uses X- and Z-directed fibres,
and Z uses X- and Y-directed fibres. The program writes charge-weighted
coordinate distributions to `seed_event_N.png` and `seed_event_N.root`.

Only fibres with `charge >= 10` are used by default. The cut is applied before
both the coordinate plots and direction-seed calculation, removing weak,
distant light-map tails. In `homo_raw` and `fiber_hits` this threshold is in the
stored photoelectron-like charge units. In `homo_truth` it applies to the
pre-fibre expected photon charge, rather than literal detected photoelectrons.

## Run the light-pattern likelihood fit

Fit one event and write a new, self-contained ROOT file:

```sh
./global_light_fit flat.root flat_with_global_fit.root \
  homo_response_250514_10bin_5m_1mm_100kPhotons.root \
  EVENT=0 TREE=homo_truth MIN_CHARGE=10
```

Fit every event:

```sh
./global_light_fit flat.root flat_with_global_fit.root \
  homo_response_250514_10bin_5m_1mm_100kPhotons.root \
  EVENT=all TREE=homo_truth MIN_CHARGE=10
```

DBSCAN is enabled by default. All options use `KEY=value`; run
`./global_light_fit --help` or run it without arguments for the complete list.
The default epsilon of 14.2 mm is slightly
larger than `10*sqrt(2)` mm, so fibres touching horizontally, vertically, or
diagonally belong to the same neighbourhood. DBSCAN runs independently in XY,
XZ, and YZ; the cluster with the largest summed charge in each projection is
retained. For the intended detector-level test:

```sh
./global_light_fit flat.root fibre_fit.root \
  homo_response_250514_10bin_5m_1mm_100kPhotons.root \
  EVENT=all TREE=fiber_hits MIN_CHARGE=10 DBSCAN=1 \
  DBSCAN_EPSILON_MM=14.2 DBSCAN_MIN_POINTS=2
```

To compare without clustering, set the DBSCAN argument to zero:

```sh
./global_light_fit flat.root fibre_fit_no_dbscan.root \
  homo_response_250514_10bin_5m_1mm_100kPhotons.root \
  EVENT=all TREE=fiber_hits MIN_CHARGE=10 DBSCAN=0
```

The peak corridor is the next optional pair of arguments: enabled and
half-width in fibres. It is disabled by default. This enables DBSCAN and then
keeps the peak fibre plus one transverse neighbour on either side in every
10 mm longitudinal slice:

```sh
./global_light_fit flat.root fibre_fit_corridor.root \
  homo_response_250514_10bin_5m_1mm_100kPhotons.root \
  EVENT=all TREE=fiber_hits MIN_CHARGE=10 DBSCAN=1 \
  DBSCAN_EPSILON_MM=14.2 DBSCAN_MIN_POINTS=2 CORRIDOR=1 \
  CORRIDOR_HALF_WIDTH_FIBRES=1
```

The provisional post-DBSCAN direction defines a fixed straight corridor in
XY, XZ, and YZ. Fibres farther than `10*half_width+0.5 mm` from the projected
seed line are removed; unlike the earlier peak-following version, this cannot
turn and follow a delta-electron branch. The direction seed is recalculated
after applying the corridor. The output records
`observations_before_clustering`, `observations_after_dbscan`, final
`observations`, `corridor_enabled`, and `corridor_half_width`.

The output also contains `global_fit_fibres`, with one row per selected fibre
and the copied `event`, `x`, `y`, `z`, `time`, `charge`, `geom_id`,
`projection`, `u`, and `v` values. This makes the exact fit selection available
to the display and downstream checks without rerunning the clustering.

The input and output names must differ. The output
clones the newest cycle of every input tree and adds a `global_fit` tree; the
input is never modified.

The likelihood uses only the selected fibre charges and the light map. It does
not read `mc_virtual_segments`, `mc_track_points`, reconstructed tracks or 3D
hits. It is a shape likelihood: predicted light is normalized over all fibres,
so no light-yield normalization is fitted. Constant dE/dx therefore cancels.
The fitted point is the line's point of closest approach to the origin. This
removes the unobservable freedom to slide an infinite straight line along its
own direction while still allowing its transverse position to move.

The current integration samples the track every 2 mm. `status == 0` denotes a
converged ROOT minimization; nonzero statuses are retained in the output and
must not be treated as successful fits.

## Event display

With an X display available in the container, run:

```sh
./global_fit_display flat_with_global_fit.root
```

Charge colouring is logarithmic by default. It can be selected explicitly or
changed to a linear scale or uniform markers:

```sh
./global_fit_display flat_with_global_fit.root EVENT=0 MODE=log
./global_fit_display flat_with_global_fit.root EVENT=0 MODE=linear
./global_fit_display flat_with_global_fit.root EVENT=0 MODE=all
./global_fit_display flat_with_global_fit.root EVENT=0 MODE=all-log
./global_fit_display flat_with_global_fit.root EVENT=0 MODE=all-linear
./global_fit_display flat_with_global_fit.root EVENT=0 MODE=selected
```

`all` is the classical display: it reads the original fit-input fibre tree,
draws every fibre passing the charge cut uniformly, and completely ignores the
`global_fit_fibres` selection. `selected` draws only fit-selected fibres. In
the charge-coloured `log` and `linear` modes, black selection outlines are
automatically suppressed when every charge-selected fibre was retained.

The canvas contains XY, XZ and YZ fibre views plus a rotatable ROOT 3D view.
All fibres passing the charge cut are charge-coloured, fibres actually selected
for the fit have open black markers, the fitted line is red, and primary MC
segments are green. MC is display-only
and never enters the likelihood. In the terminal use `n`, `p`, an event number, or `q`
for next, previous, jump and quit. Only events present in `global_fit` are
offered. The 3D points are representative fibre positions; the three 2D views
are the geometrically meaningful fibre measurements.

## Direction analysis

```sh
./global_fit_analysis flat_with_global_fit.root OUTPUT=direction_analysis.root
```

This writes component correlations for MC versus seed and MC versus fitted
directions, signed MC–seed and MC–fit angular residuals, and a per-event
`direction_comparison` tree. It also creates `direction_comparison.png`. The MC
direction is obtained from the dominant HOMO primary's ordered segment
displacements and is used only by this validation program.

The fitter also stores data-variance Poisson `chi2 = sum((q-mu)^2/q)`, `ndof`,
and hence `chi2/ndof`. The comparison is conditioned on fibres passing the
charge cut, with predicted charge normalized to their observed total. This
form remains finite when a selected fibre currently has zero predicted light.
There are four fitted line
parameters and one normalization constraint, so `ndof = N_fibres - 5`.
`global_fit_analysis` writes `fit_chi2`, `fit_ndof`, and `fit_chi2_ndof`
histograms and the summary image `fit_quality.png`. The event display prints
all three values for the selected event.

The final optional fitter argument is a minimum light-map fraction. For
example, this keeps only map contributions of at least `1e-4`:

```sh
./global_light_fit flat.root fibre_fit_fraction.root \
  homo_response_250514_10bin_5m_1mm_100kPhotons.root \
  EVENT=all TREE=fiber_hits MIN_CHARGE=10 DBSCAN=1 \
  DBSCAN_EPSILON_MM=14.2 DBSCAN_MIN_POINTS=2 CORRIDOR=1 \
  CORRIDOR_HALF_WIDTH_FIBRES=1 MIN_MAP_FRACTION=1e-4
```

The value is stored as `minimum_map_fraction` in `global_fit`. Zero retains the
complete map. This selection changes the forward model, not the measured-fibre
charge cut, so the zero and `1e-4` results should be compared explicitly.

The analysis also associates `mc_virtual_segments.primary_id` with muon
trajectories (`abs(pdg)==13`) from `mc_track_points`. It deduplicates shared
segment start/end points at 1 micrometre precision and calculates the minimum
perpendicular distance of every unique point to the infinite fitted line. The
distances are stored in `muon_fit_distances`, histogrammed as
`muon_fit_point_distance`, and fitted around the histogram peak with
`muon_distance_gaussian`. The plot is written as `muon_fit_distance.png`.

Because a 3D line has a two-dimensional transverse plane rather than one
intrinsic positive side, the analysis additionally stores signed perpendicular
residuals in XY, XZ, and YZ. Positive is the left-hand normal
`(-d_b,d_a)` when looking along the fitted projected direction. The branches
are `residual_xy`, `residual_xz`, and `residual_yz`; corresponding histograms
and Gaussian fits are written to the analysis ROOT file and plotted in
`muon_fit_signed_residuals.png`. The unsigned 3D distance is retained.

Unsigned minimum distances are also stored separately as `distance_xy`,
`distance_xz`, and `distance_yz`. Their histograms and Gaussian peak fits are
written to the ROOT output, with the combined plot `muon_fit_view_distances.png`.
