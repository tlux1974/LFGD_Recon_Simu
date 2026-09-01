# Light-map analyser manual

`analyse_lightmap.py` analyses a complete HOMO light map independently of how
the map was produced. It works with full-cube productions and maps reconstructed
from an octant.

## Requirements

The analyser uses Python 3, PyROOT, NumPy, and Matplotlib. Run it from an
environment where ROOT is configured. The input contains the valid
scintillator source positions in the 10 x 10 x 10 scan grid. Their number is
geometry-dependent: 970 for 1.0 mm fibres and 850 for 2.0 mm fibres with the
current grid and fibre offsets. Source centres inside fibres are absent.

## Accepted inputs

The positional input can be either:

- a merged ROOT light map containing the `position_efficiency` tree; or
- an efficiency CSV produced by `ExportLightMapEfficiency.C`.

Legacy ROOT maps without `position_efficiency` must first be exported to CSV:

```bash
root -l -b -q \
  'ExportLightMapEfficiency.C("legacy_map.root","legacy_map_efficiency.csv")'
```

For per-view analysis, retain the ROOT structure and add the missing tree to a
copy of the legacy map. The photon count must match the production:

```sh
cp lightmaps/homo_response_250514_10bin_5m_1mm_100kPhotons.root lightmaps/homo_response_250514_10bin_5m_1mm_100kPhotons_with_efficiency.root
root -l -b -q 'AddLightMapEfficiencyTree.C("lightmaps/homo_response_250514_10bin_5m_1mm_100kPhotons_with_efficiency.root",100000)'
```

The updater refuses to replace an existing `position_efficiency` tree.

## Create the merged ROOT input

Run the appropriate merger from `Studies/PICLightMap` after copying all chunk
directories into that directory.

For a full-cube production consisting of ten `chunk_*` directories:

```bash
./merge_lightmap.sh --scat-length 0.7 --fibre-diameter 1.0
```

For a one-octant production consisting of the five directories
`octant_chunk_0`, `octant_chunk_25`, `octant_chunk_50`, `octant_chunk_75`, and
`octant_chunk_100`:

```bash
./merge_lightmap_octant.sh \
  homo_response_octant_mirrored_10bin_5m_0p7mm_100kPhotons.root
```

The octant merger requires all five ROOT chunks and mirrors them into a full
map using the original map's valid-position mask. It refuses output unless the
result contains exactly 970 positions. Both merger scripts also write the
corresponding `_efficiency.csv` automatically.

The merger executable needs the configured ND280 environment. Perform this
step inside the writable local sandbox container used for development and
merging. The PIC SIF is only needed by submitted cluster jobs and is not used
for local aggregation. Inside the sandbox container, source its normal ND280
environment if necessary, then run:

```bash
source /usr/local/t2k/current/nd280SoftwarePilot/nd280SoftwarePilot.profile
source /usr/local/t2k/current/nd280SoftwareMaster_14.36-plusplus.0.3/bin/setup.sh
cd /home/tlux/HK/ND280++/LFGD_Recon_Simu/Studies/PICLightMap
./merge_lightmap_octant.sh \
  homo_response_octant_mirrored_10bin_5m_0p7mm_100kPhotons.root
```

Before merging an octant, verify the five inputs:

```bash
find octant_chunk_* -path '*/results/*.root' -type f | sort
```

Exactly five files should be listed.

## Command-line usage

From `Studies/PICLightMap`:

```bash
python3 analyse_lightmap.py MAP.root \
  --label 'scattering 0.7 mm, absorption 5 m' \
  --output-dir lightmap_analysis/0p7mm
```

The equivalent command using a CSV is:

```bash
python3 analyse_lightmap.py MAP_efficiency.csv \
  --label 'scattering 0.7 mm, absorption 5 m' \
  --output-dir lightmap_analysis/0p7mm
```

Arguments:

- `input` is the required ROOT or CSV input path.
- `--label` sets the descriptive name printed in the terminal and plots. If
  omitted, the input filename stem is used.
- `--output-dir` selects the result directory. If omitted, results go to
  `lightmap_analysis/<input filename stem>/`.

Display the built-in command summary with:

```bash
python3 analyse_lightmap.py --help
```

## Output files

The selected output directory contains:

- `summary.csv`: numerical statistics and the lowest-efficiency coordinates;
- `efficiency_distribution.png`: 1D distribution of the 970 efficiencies;
- `relative_efficiency_distribution.png`: the same distribution expressed as
  percentage deviation from the mean;
- `efficiency_vs_fibre_distance.png`: individual efficiencies and the mean
  plus population standard deviation versus distance to the nearest fibre
  axis;
- `efficiency_vs_fibre_distance.csv`: the underlying distance and efficiency
  for every valid position;
- `efficiency_z_905.5mm.png` through `efficiency_z_914.5mm.png`: ten separate
  x-y efficiency maps, one per z layer;
- `efficiency_10_slices.png`: all ten layers in one overview, using a common
  colour scale so that layers can be compared directly.

When the input is the merged ROOT map, the analyser additionally writes:

- `efficiency_per_view.csv`: total efficiency, efficiencies for X-, Y-, and
  Z-directed fibres, and their population RMS for every valid subvoxel;
- `total_efficiency_10_slices.png`: the ten total-efficiency slices;
- `x_fibre_efficiency_10_slices.png`, `y_fibre_efficiency_10_slices.png`, and
  `z_fibre_efficiency_10_slices.png`: ten slices for each fibre view;
- `view_efficiency_rms_10_slices.png`: the ten slices of the population RMS
  among the three view efficiencies.

The directional plots require the ROOT map because an exported efficiency CSV
contains only the total efficiency.

Grey cells in the plots are the excluded fibre-source positions, not measured
zero efficiencies.

Both 1D histograms overlay the Gaussian defined by the measured mean and
population standard deviation. This is a visual reference, not an assumption
or forced Gaussian fit. The plots also show skewness and excess kurtosis:
values near zero are consistent with a roughly symmetric, Gaussian-shaped
distribution, while sizeable values reveal asymmetry or non-Gaussian tails.

## Interpreting fibre-adjacent minima

The analyser infers the three fibre axes from the 30 missing grid positions
and plots efficiency against the shortest perpendicular distance to those
axes. This tests the observed tendency for efficiency to be lowest in the
roughly 1 mm neighbourhood of a fibre.

Such a minimum is physically plausible even though the fibre is nearby. The
reported efficiency is the total number of counted fibre hits divided by the
number of launched photons; reaching a nearby fibre boundary does not
guarantee capture as a recorded hit. Fibre and cladding volumes displace
scintillator, incidence angles can be unfavourable, and the nearest fibre can
shadow photon paths toward the other two views. With strong scattering, a
slightly more distant source can also sample more directions and gain more
opportunities to reach any fibre.

Positions inside fibres are absent rather than recorded with zero efficiency.
Consequently, a depression around a grey/missing fibre line appears in its
neighbouring cells. The distance plot helps quantify this feature. For ROOT
input, the analyser now also reads the per-fibre response objects and separates
the three fibre directions; the `position_efficiency` tree itself still stores
only their total.

## Reported statistics

For the 970 efficiencies, `summary.csv` reports:

- `mean`: arithmetic mean efficiency;
- `std`: population standard deviation;
- `coefficient_of_variation`: relative standard deviation, `std / mean`;
- `skewness`: asymmetry of the 1D distribution;
- `excess_kurtosis`: tail/peak measure relative to a Gaussian, whose value is
  zero;
- `min` and `max`: lowest and highest efficiency;
- `relative_spread`: full relative range, `(max - min) / mean`;
- `min_x_mm`, `min_y_mm`, and `min_z_mm`: location of the minimum.

The two relative-spread measures answer different questions. The coefficient
of variation measures the typical variation across the map, while the full
relative range is controlled by the two most extreme positions.

Values in `summary.csv` are stored as fractions. Multiply by 100 to express
them as percentages. For example, `0.00236` corresponds to `0.236%`.

## Jupyter notebook

`LightMapAnalysis.ipynb` provides an interactive front end to the same analysis
code. Start Jupyter in this directory, open the notebook, and edit:

```python
input_file = Path('lightmaps/homo_response_10bin_5m_0p7mm_100kPhotons.root')
label = 'scattering 0.7 mm, absorption 5 m'
```

Running the analysis cell calls the same `read_efficiencies`,
`read_view_efficiencies`, and `analyse` functions used by the command-line
program, so notebook and batch results are consistent. The following cell
displays the total, three directional, and RMS ten-slice overviews inline.

## Central-subvoxel fibre-correlation kernel

`analyse_fibre_correlation.py` analyses the individual fibre fractions rather
than their sum. It fixes the source at the recommended central valid subvoxel
`(0.5, 30.5, 910.5) mm`. For a hypothetical track along each coordinate axis,
it selects the two fibre families perpendicular to that axis and finds the
nearest fibre in each family.

For fibre family `k`, offset by `n` pitches along the track direction, it
calculates

```text
r_k(n) = F_k(n) / F_k(0).
```

It also reports the combined pair response

```text
C(n) = [F_1(n) + F_2(n)] / [F_1(0) + F_2(0)].
```

For every `n`, the CSV includes the perpendicular source-to-fibre-axis
distance for each fibre, plus their arithmetic mean. Because the fibre axes do
not pass through the subvoxel centre, `+n` and `-n` generally correspond to
different physical distances. The signed-offset plot therefore need not be
symmetric even for symmetric optical transport.

This measures how far the central subvoxel's light pattern extends along the
two fibre families relevant to a track. Positive and negative `n` are retained
separately, and one pitch is 10 mm. Run it on a merged ROOT map with:

```bash
python3 analyse_fibre_correlation.py MAP.root \
  --label 'scattering 0.4 mm, absorption 5 m' \
  --output-dir lightmap_analysis/0p4mm_full
```

It writes `central_fibre_correlation.csv`, the signed-pitch plot
`central_fibre_correlation.png`, and
`central_fibre_correlation_vs_distance.png`. The latter plots the combined
response against the mean physical distance to the fibre pair and is the more
appropriate view for judging whether the apparent signed asymmetry is real. A
ROOT map is required because the efficiency CSV contains only the total
response per source position, not the individual fibre fractions.

## Overlay curves from several settings

`plot_csv_curves.py` overlays any number of CSV results. Write inputs as
`LABEL=FILE.csv`; the label appears in the legend. Its defaults plot
`combined_pair_ratio` against `mean_pair_distance_mm`, with separate panels
for the x-, y-, and z-track directions:

```bash
python3 plot_csv_curves.py \
  '0.25 mm=lightmap_analysis/0p25mm_full/central_fibre_correlation.csv' \
  '0.4 mm=lightmap_analysis/0p4mm_full/central_fibre_correlation.csv' \
  '0.7 mm=lightmap_analysis/0p7mm_octant/central_fibre_correlation.csv' \
  --title 'Central-subvoxel fibre correlation' \
  --xlabel 'Mean distance to fibre pair [mm]' \
  --ylabel 'Combined response / nearest pair' \
  --output lightmap_analysis/fibre_correlation_comparison.png
```

The interface is generic for future CSV curves. Select other columns with
`--x-column` and `--y-column`, and another panel category with
`--group-column`. Use `--group-column ''` for a single panel, `--no-connect`
for unconnected points, and `--xscale`/`--yscale` to select linear or
logarithmic axes. Any number of inputs is accepted. See every option with:

```bash
python3 plot_csv_curves.py --help
```

## Integrated ten-subvoxel track correlation

`analyse_track_correlation.py` models a straight 10 mm track parallel to each
coordinate axis (and therefore parallel to one fibre family). For each track
axis it fixes the other two cube-relative coordinates at `+0.5 mm`, evaluates
the ten source positions from `-4.5 mm` through `+4.5 mm`, and sums every
individual fibre fraction over those ten subvoxels.

The two fibre families perpendicular to the track are then normalized to their
nearest-fibre responses. The combined result is

```text
C_track(n) = [sum_i F_1(i,n) + sum_i F_2(i,n)]
             / [sum_i F_1(i,0) + sum_i F_2(i,0)].
```

Run it with a merged ROOT map:

```bash
python3 analyse_track_correlation.py lightmaps/MAP.root --label '0.4 mm full production' --output-dir lightmap_analysis/0p4mm_full
```

It writes `integrated_track_fibre_correlation.csv` and
`integrated_track_fibre_correlation.png`. For the 0.4 mm full-production map,
the combined response at the adjacent fibres (`n=+/-1`) is about `0.17-0.21`,
larger than the central point-source value because the track endpoints are
closer to neighbouring fibre axes.

The same program also calculates the transverse containment of the complete
recorded fibre response. For every fibre it computes the shortest distance
between that fibre axis and the central track axis, sums responses in radial
shells, and constructs the cumulative fraction. It reports

```text
R99 = smallest radius containing at least 99% of recorded fibre response.
```

The outputs are `integrated_track_transverse_containment.csv` and
`integrated_track_transverse_containment.png`. Because fibres occupy discrete
axes, `R99` is a discrete containment radius rather than an interpolated fit.
Also note that this is 99% of photons recorded in fibres, not 99% of all
launched photons: the light map does not store the final absorption location
of photons that never reach a fibre. The even 10-bin grid has no source centre
at exactly zero, so the discretized central track uses transverse coordinates
`(+0.5,+0.5) mm`.

Overlay several integrated-track results with the generic plotter. Since this
CSV uses pitch offset rather than mean point-to-fibre distance, select
`offset_mm` explicitly:

```bash
python3 plot_csv_curves.py '0.25 mm=lightmap_analysis/0p25mm/integrated_track_fibre_correlation.csv' '0.4 mm=lightmap_analysis/0p4mm_full/integrated_track_fibre_correlation.csv' '0.7 mm=lightmap_analysis/0p7mm_octant/integrated_track_fibre_correlation.csv' '1.0 mm=lightmap_analysis/1p0mm/integrated_track_fibre_correlation.csv' --x-column offset_mm --y-column combined_pair_ratio --group-column track_axis --xlabel 'Fibre offset along track [mm]' --ylabel 'Integrated response / nearest pair' --title 'Ten-subvoxel track correlation' --output lightmap_analysis/integrated_track_correlation_comparison.png
```

## Example: existing 0.25 mm production

```bash
python3 analyse_lightmap.py \
  homo_response_10bin_5m_0p25mm_100kPhotons.root \
  --label 'scattering 0.25 mm, absorption 5 m' \
  --output-dir lightmap_analysis/0p25mm
```

This validated example contains 970 positions and produced a coefficient of
variation of `0.235996%`, a full relative range of `1.355226%`, and a minimum
efficiency of `0.80642` at `(3.5, 27.5, 908.5) mm`.

## Common errors

- An unexpected position count from `merge_lightmap.sh` means that the map
  maker's source mask does not agree with the requested fibre diameter.
- `has no position_efficiency tree`: export the legacy map with
  `ExportLightMapEfficiency.C` and analyse the resulting CSV.
- Python reports that `ROOT` or `matplotlib` is unavailable: source the ND280
  environment or run inside the configured project container.
