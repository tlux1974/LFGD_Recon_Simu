# Global HOMO light-distribution fit test

This is a standalone ROOT-based test. The normal likelihood and
`SEED_DIRECTION=FIBRE` reconstruction use only fibre information. MC segments
are copied for display and diagnostic benchmarks. The explicitly
truth-assisted `SEED_DIRECTION=MC_SEGMENT_10` mode also reads them and must not
be treated as reconstruction performance.

## Build

Enter the ND280 container, source the standard ND280 setup, and run:

```sh
cd LFGD_Recon_Simu/Studies/GlobalFitTest
make
```

## Produce the initial fit seed

For interactive event-by-event threshold studies, open
`global_fit_seed_study.ipynb` in Jupyter. Set `INPUT_FILE` to a branch's
`flat.root`, run all cells, and use the event and charge-threshold controls.
The notebook shows XY, XZ, and YZ fibre distributions plus a combined 3D
panel, overlays MC segments in green, draws the current data-only seed in red,
and reports its unsigned angular difference from the MC direction. MC is not
used to construct the data seed.

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

Both fitters restrict the predicted line to a data-derived 3D range by default
rather than integrating it across the complete detector. For X, Y, and Z they
use only fibre views that measure the respective coordinate, trim the lower
and upper 1% of charge, add 10 mm padding, and clamp the result to the detector.
This is `FIT_RANGE=1`, the default. `FIT_RANGE_QUANTILE=0.01` trims 1% of the
charge independently from each end of each measured coordinate distribution;
it is a charge fraction, not a distance. `FIT_RANGE_PADDING_MM=10` then extends
both ends by 10 mm to retain endpoint light. Set `FIT_RANGE=0` to recover the
old detector-spanning line; quantile and padding are then ignored. Each event's
exact settings and bounds are stored as `fit_range_enabled`, `fit_low_*`,
`fit_high_*`, `fit_range_quantile`, and `fit_range_padding`.

The forward light-map integral is sampled at 1 mm in both the straight
geometry fit and the column-charge fit, matching the light map's subvoxel
pitch. The RL refinement uses the same 1 mm default through
`RESPONSE_STEP_MM=1`.

```sh
./global_light_fit flat.root fitted.root lightmap.root EVENT=all FIT_RANGE=1 FIT_RANGE_QUANTILE=0.01 FIT_RANGE_PADDING_MM=10
./global_light_fit flat.root fitted_detector_span.root lightmap.root EVENT=all FIT_RANGE=0
```

### Experimental two-stage column fit

`global_light_fit_columns` is a separate experimental executable; the baseline
`global_light_fit` behaviour is unchanged. The dominant axis is chosen from
the data seed. The geometry likelihood integrates the line in exact 10 mm
dominant-axis columns while retaining constant charge per path length. After
the geometry fit, its parameters are fixed and a non-negative Poisson/EM fit
releases one charge contribution per crossed column. The result is stored in
`global_fit_column_charges` with the axis, column bounds, path length, fitted
selected-fibre charge, charge per millimetre, and response scale.
The `global_fit` row also records `column_fit_converged`,
`column_fit_iterations`, `column_fit_max_relative_change`,
`column_fit_active_columns`, `column_nll`, `column_chi2`, and `column_ndof`.
The convergence flag requires the largest relative column-charge update to be
below `1e-6` before the 2000-update limit.

An optional one-sided MIP prior can suppress implausibly low internal-column
charges without penalizing the Landau tail above the MIP most-probable value.
It is applied after the ordinary Poisson/EM column solution and minimizes the
combined light likelihood and squared lower-side penalty with backtracking.
The floor and MPV scale with each column's path length relative to 10 mm. The
first and last columns are always excluded; other short columns are excluded
using `COLUMN_MIP_MIN_PATH_FRACTION` relative to the median internal path.
The required `COLUMN_CHARGE_PER_ENERGY` must be fixed from an independent
calibration sample, not the sample whose reconstruction performance is being
measured. The feature is disabled by default.

```sh
Studies/GlobalFitTest/global_light_fit_columns flat.root global_fit_columns_mip_1mm.root lightmap.root EVENT=all TREE=fiber_hits MIN_CHARGE=3 DBSCAN=1 DBSCAN_EPSILON_MM=14.2 DBSCAN_MIN_POINTS=2 SEED_DIRECTION=VIEW_MEDIAN_DIAMETER MAX_FUNCTION_CALLS=5000 TOLERANCE=1e-2 FIT_RANGE=1 COLUMN_MIP_PRIOR=1 COLUMN_MIP_FLOOR_ENERGY=0.5 COLUMN_MIP_MPV_ENERGY=1.5 COLUMN_MIP_STRENGTH=0.1 COLUMN_CHARGE_PER_ENERGY=65.98 COLUMN_MIP_MIN_PATH_FRACTION=0.8
```

The column tree retains both solutions: ordinary fields contain the
regularized result, while `unregularized_fitted_selected_fibre_charge` and
`unregularized_response_scale` retain the initial EM result. It also records
eligibility and the path-scaled floor and MPV per column. The `global_fit`
tree records prior convergence/stalling, iterations, eligible-column count,
penalty, combined objective, calibration and all prior settings.

### Offline RL segment refinement

`rl_segment_refinement` is an optional third stage and never reruns the seed,
geometry fit, or fixed-geometry column fit. The second-stage column-charge EM
update is already Richardson--Lucy for fixed segment positions. This program
extends its components to a transverse position grid around every fitted-line
column. The lightmap remains the forward response operator.

The starting position weights are an aggregate two-dimensional KDE of the
primary-muon reco--MC line residuals over all events in the input file. This is
a closure-test prior, not event-specific MC truth. For a final performance
measurement, derive and freeze it using an independent MC training sample.
Because an independent displacement PDF does not prevent column-to-column
zigzags, every RL iteration also applies a discrete-curvature prior based on
`offset[k+1] - 2*offset[k] + offset[k-1]`. It permits a changing local slope
while suppressing unsupported rapid oscillation. Both priors preserve each
column's total response-scale amplitude during their regularization step. The
accepted objective is the charge-normalized light-pattern NLL plus the
weighted displacement and curvature penalties. Backtracking accepts an update
only when that combined objective decreases.

Build only this independent executable:

```sh
make -C Studies/GlobalFitTest rl_segment_refinement
```

Initial one-event closure test for the 0.4 mm lightmap (one-line prompt form):

```sh
d=Studies/LightmapsStudie/5GeV/maps/10bin_5m_0p4mm_100kPhotons; Studies/GlobalFitTest/rl_segment_refinement "$d/global_fit_columns_particle_ids.root" "$d/rl_segments_event0.root" "$(cat "$d/lightmap_path.txt")" EVENT=0 GRID_RADIUS_MM=3 GRID_STEP_MM=1 ITERATIONS=50 CONVERGENCE=1e-4 PRIOR_BANDWIDTH_MM=0.5 PRIOR_STRENGTH=0.02 CURVATURE_SIGMA_MM=0.75 CURVATURE_STRENGTH=0.20 RESPONSE_STEP_MM=1
```

After inspecting that output, replace `EVENT=0` with `EVENT=all` and use a new
output filename. The result contains `rl_fit`, `rl_iterations`, `rl_columns`,
and (unless `SAVE_CANDIDATES=0`) `rl_candidates`. It records NLL evolution,
amplitude convergence, curvature RMS, refined transverse offsets and centres,
column response scales, predicted selected-fibre charge, curvature components,
the final candidate distribution, and the prior/regularization settings.
`rl_iterations` stores the NLL before an update, immediately after pure RL,
and after accepted regularization, together with both penalties, the combined
objective, accepted line-search fraction, amplitude change, centroid-position
change, and curvature RMS. A rejected line search marks the event as stalled.
The output preserves two solutions: the final/best accepted penalized-objective
solution in the ordinary `rl_columns` and `rl_candidates` fields, and the
lowest raw-light-NLL solution in fields prefixed with `best_nll_`. `rl_fit`
records both selected iteration numbers and both NLL values.

Overlay either or both reconstructed RL paths on the fibre display and MC
truth using the original column-fit file as the primary input:

```sh
d=Studies/LightmapsStudie/5GeV/maps/10bin_5m_0p4mm_100kPhotons; Studies/GlobalFitTest/global_fit_display "$d/global_fit_columns_particle_ids.root" EVENT=0 MODE=all-log RL="$d/rl_segments_all.root" RL_SOLUTION=both
```

The straight global fit is red, primary-muon MC cube segments are green, the
regularized RL path is magenta, and the lowest-raw-NLL path is orange dashed.
`RL_SOLUTION` also accepts `regularized` and `best-nll`.

The main controls are `GRID_RADIUS_MM`, `GRID_STEP_MM`, `ITERATIONS`,
`CONVERGENCE`, `POSITION_CONVERGENCE_MM`, `OBJECTIVE_CONVERGENCE`,
`PRIOR_BANDWIDTH_MM`, `PRIOR_MAX_RESIDUAL_MM`,
`PRIOR_STRENGTH`, `CURVATURE_SIGMA_MM`, `CURVATURE_STRENGTH`, and
`RESPONSE_STEP_MM`. Run the executable without arguments for its complete
help. Multiplicative RL cannot populate an exactly zero component; therefore
the aggregate residual PDF initializes every supported grid candidate with a
strictly positive amplitude. `CONVERGENCE` is applied to the largest
per-column L1 fractional change of the candidate-amplitude distribution; this
avoids tiny tail candidates falsely preventing convergence.

The column charge is normalized on the selected fibre set; it is not yet an
absolute deposited-energy calibration. Internal columns are defined by their
dominant-axis span, avoiding independent parameters for short pieces created
when a diagonal line clips a cube. The first and last detector columns can
still have short path lengths.

Build this variant explicitly when ready:

```sh
make -B -C Studies/GlobalFitTest global_light_fit_columns
```

It accepts the same options as `global_light_fit`:

```sh
Studies/GlobalFitTest/global_light_fit_columns INPUT.root OUTPUT.root LIGHTMAP.root EVENT=all TREE=fiber_hits MIN_CHARGE=3 DBSCAN=1 DBSCAN_EPSILON_MM=14.2 DBSCAN_MIN_POINTS=2 CORRIDOR=0 SEED_DIRECTION=VIEW_MEDIAN_DIAMETER SEED_MEDIAN_FACTOR=1 MAX_FUNCTION_CALLS=5000 TOLERANCE=1e-3
```

For interactive comparisons, open `global_fit_column_analysis.ipynb` with the
Jupyter Python environment. Its first configuration cell accepts any labelled
set of column-fit ROOT files. It overlays convergence, column-charge/MC-energy
profiles and residuals across light maps, and includes a pure-MC straight-line
PCA reference made from the primary muon's `mc_track_points`. The helper module
`column_fit_notebook_tools.py` uses PyROOT, pandas, matplotlib and ipywidgets.

For a truth-assisted seed diagnostic, the direction can be initialized from
the primary muon's first MC segment start toward MC segment 10:

```sh
./global_light_fit flat.root fitted_mc_seed.root lightmap.root EVENT=all TREE=fiber_hits SEED_DIRECTION=MC_SEGMENT_10
```

This option is intentionally a cheat for diagnosing seed sensitivity. It must
not be used when quoting reconstruction performance. The default remains the
data-driven fibre seed, `SEED_DIRECTION=FIBRE`.

The light-sharing seed developed for the scattering-length comparison uses a
moderate absolute charge threshold first, then DBSCAN independently in each
view. Within each retained main cluster it calculates the per-view median,
keeps fibres above `median * SEED_MEDIAN_FACTOR`, finds the maximum-distance
pair, and combines the projected XY/XZ/YZ directions into a 3D seed:

```sh
./global_light_fit flat.root fitted_median_seed.root lightmap.root EVENT=all TREE=fiber_hits MIN_CHARGE=10 DBSCAN=1 DBSCAN_EPSILON_MM=14.2 DBSCAN_MIN_POINTS=2 SEED_DIRECTION=VIEW_MEDIAN_DIAMETER SEED_MEDIAN_FACTOR=1
```

The median-filtered subset is used only to determine the seed. The likelihood
continues to use the broader DBSCAN-selected fibre set. The output records the
seed method, per-view medians, median factor, number of seed fibres, and NLL
and chi-squared evaluated at the seed before minimization.

The default minimizer settings are `MAX_FUNCTION_CALLS=5000` and
`TOLERANCE=1e-2`. Both can be overridden explicitly. The `global_fit` output
tree records `edm`, `function_calls`, `maximum_function_calls`, and `tolerance`
for every event, allowing status-3 failures to be compared with the requested
convergence threshold:

```sh
./global_light_fit flat.root fitted.root lightmap.root EVENT=all TREE=fiber_hits MIN_CHARGE=3 DBSCAN=1 DBSCAN_EPSILON_MM=14.2 DBSCAN_MIN_POINTS=2 CORRIDOR=0 SEED_DIRECTION=VIEW_MEDIAN_DIAMETER SEED_MEDIAN_FACTOR=1 MAX_FUNCTION_CALLS=5000 TOLERANCE=1e-2
```

To rerun only GlobalFit across an existing light-map study while preserving
the previous outputs, use `rerun_global_fits.sh` from `LFGD_Recon_Simu` and
give it a new output filename.

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

### Two-view fits

The `FIT_VIEWS` option controls which physical fibre projections enter hit
selection, direction seeding, DBSCAN and corridor processing, and the
likelihood normalization:

| Setting | Views retained | Fibre direction excluded | Intended track axis |
| --- | --- | --- | --- |
| `ALL` | XY, XZ, YZ | none | general 3-view fit |
| `2D_X` | XZ, XY | X-directed fibres (YZ view) | approximately X |
| `2D_Y` | YZ, XY | Y-directed fibres (XZ view) | approximately Y |
| `2D_Z` | XZ, YZ | Z-directed fibres (XY view) | approximately Z |

For example, for tracks approximately parallel to X, exclude the YZ view:

```sh
./global_light_fit flat.root fibre_fit_2d_x.root \
  homo_response_250514_10bin_5m_1mm_100kPhotons.root \
  EVENT=all TREE=fiber_hits MIN_CHARGE=10 FIT_VIEWS=2D_X
```

The fit still returns a 3D point and 3D direction; "2D" means that two fibre
projections constrain that line. Equivalent explicit lists are accepted, for
example `FIT_VIEWS=XZ,XY`. The chosen canonical view list is stored as
`fit_views` in the `global_fit` tree, while `global_fit_fibres` contains only
fibres from those views. The interactive display reads this stored selection;
there is no separate display-time `FIT_VIEWS` switch.

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

## MC-segment straight-line benchmark

Every `global_light_fit` run also performs an independent best-case benchmark
that never reads fibre charges or the light map. For each requested event it
associates `mc_virtual_segments.primary_id` with `abs(pdg)==13` trajectories in
`mc_track_points`, restricts the segments to detector 0, and, if more than one
muon ID exists, selects the ID with the largest number of segments. A 3D
orthogonal-regression/PCA line is fitted using only the segment **start**
positions. Its sign is chosen to follow the summed start-to-stop direction.

The benchmark is written into the same output ROOT file as the fibre
likelihood result:

- `mc_segment_line_fit`: one row per event with `primary_id`, number of start
  points, fitted point and direction, `sum_squared_residual_mm2`,
  `rms_residual_mm`, and `ndof=2*N-4`;
- `mc_segment_line_residuals`: one row per selected segment, with separate
  `start_distance` and `end_distance`, signed Cartesian residual-vector
  components, and signed residuals in XY, XZ, and YZ for both endpoints;
- `mc_segment_start_distance` and `mc_segment_end_distance`: separate ROOT
  histograms of the unsigned 3D perpendicular distances.

The Cartesian residual vector is the point-to-line displacement after
removing its component parallel to the fitted line. A single 3D distance has
no intrinsic sign, so `start_distance` and `end_distance` are non-negative;
the per-axis and projected-view residual branches retain signed information.
Comparing starts and ends tests extrapolation within each virtual-cube segment,
while the width versus track length gives a direct view of the limitation from
multiple Coulomb scattering.

Running `global_fit_analysis` on this file additionally creates
`fibre_vs_mc_line_start_residuals.png` and
`fibre_vs_mc_line_end_residuals.png`. Each canvas overlays normalized signed
residuals from the fibre/light-pattern line and the MC-start PCA line in XY,
XZ, and YZ over `-10..+10 mm`. Independent Gaussian fits are made around each
central peak, and their fitted mean and sigma are printed in the legends. The
underlying six fibre histograms, six MC-line histograms, Gaussian functions,
and canvases are also saved in the analysis ROOT output.

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
./global_fit_display flat_with_global_fit.root EVENT=0 MODE=dbscan
```

`all` is the classical display: it reads the original fit-input fibre tree,
draws every fibre passing the charge cut uniformly, and completely ignores the
`global_fit_fibres` selection. `selected` draws only fit-selected fibres. In
the charge-coloured `log` and `linear` modes, black selection outlines are
automatically suppressed when every charge-selected fibre was retained.
`dbscan` overlays every fibre passing the fit's charge threshold and outlines
the fibres retained immediately after DBSCAN in magenta. This is the selection
before any optional corridor cut. The fit output stores that intermediate set
in `global_fit_dbscan_fibres`; older fit files must be regenerated to use this
mode.

The canvas contains XY, XZ and YZ fibre views plus a rotatable ROOT 3D view.
All fibres passing the charge cut are charge-coloured, fibres actually selected
for the fit have open black markers, the fitted line is red, and primary MC
segments are green. MC is display-only
and never enters the likelihood. In the terminal use `n`, `p`, an event number, or `q`
for next, previous, jump and quit. Only events present in `global_fit` are
offered. The 3D points are representative fibre positions; the three 2D views
are the geometrically meaningful fibre measurements.

## Direction analysis

The current analysis explicitly compares initialization and minimization. Its
`direction_comparison` tree contains seed-to-MC, fit-to-MC, and seed-to-fit
angles; X/Y/Z direction-component residuals; Minuit status; NLL at the seed
and fit; and the corresponding NLL and chi-squared improvements.

For every primary-muon MC segment start and end, the
`mc_seed_fit_residuals` tree stores the X/Y/Z residual vector and perpendicular
distance to both the seed and fitted lines. The
`mc_seed_fit_xyz_residuals.png` summary plots these component distributions,
while `direction_comparison.png` includes the seed-to-fit angular movement.

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
