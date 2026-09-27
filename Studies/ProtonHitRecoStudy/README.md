# Isotropic 1 GeV proton hit-position study

Additional source documents:

- `../../RECONSTRUCTION_MANUAL_DRAFT_SOURCE.md` is the project-wide inventory
  of generation, response, hit-reconstruction and official-fitting options
  for preparing a complete manual.  It is not specific to this study.
- `MATCHED_MUON_TEST_REPRODUCTION.md` gives the complete matched HFGD/LFGD
  1 GeV muon workflow used in the current ten-event comparison.

This study is deliberately separate from `GlobalFitTest`.  It compares the
current HOMO reconstruction with an alternative position construction while
leaving the standard charge sharing, HFG clustering, and track fitter
unchanged.

The study uses only the standard 1 mm-fibre light map supplied by the ND280
container:

```text
${DETRESPONSESIMROOT}/input/homo_response_250514_10bin_5m_1mm_100kPhotons.root
```

Accordingly, the production command below deliberately does not set
`HOMO_LIGHTMAP_FILE`.  Setting that variable would override the standard map.

Both reconstructions use the standard configured local-direction 2D peak
clustering, its charge-centre-of-gravity positions, and charge-reconstruction
mode 0.  The default reconstruction combines
all available views.  The alternative estimates a local 3D direction, selects
the dominant detector axis, fixes that coordinate to the centre of its 10 mm
slice, and obtains the two transverse coordinates only from the two fibre
views perpendicular to that axis.  A direction farther than 45 degrees from
every detector axis retains the default continuous position.

The generated GPS macro explicitly contains `/t2k/field 0.0 tesla`; no extra
switch is required to disable the magnetic field.

## Build

Inside the normal ND280++ container, rebuild the locally modified packages:

```bash
singularity exec ND280ppCont bash -lc 'source /usr/local/t2k/current/nd280SoftwarePilot/nd280SoftwarePilot.profile >/dev/null 2>&1; source /usr/local/t2k/current/nd280SoftwareMaster_14.36-plusplus.0.3/bin/setup.sh >/dev/null 2>&1; cmake --build /home/tlux/HK/ND280++/SoftProj/hfgrecon/Linux-AlmaLinux_9.5-gcc_11-x86_64 -j2; cd /home/tlux/HK/ND280++/LFGD_Recon_Simu; ./build_flat_treemaker.sh'
```

## Generate one standard-lightmap proton sample

This command runs Geant4, detector response with the container's standard
light map, the default reconstruction, the standard track fitter, and the
flat tree. It does not overwrite an existing run directory.

```bash
PARTICLE=proton ENERGY_MEV=1000 DIRECTION_MODE=isotropic LFGD_Recon_Simu/run_student_sample.sh lfg-default-position 1000 proton1GeV_standard_default
```

The output directory is
`LFGD_Recon_Simu/output/homo_student_lfg-default-position_proton1GeV_standard_default`.

## Rerun only reconstruction with the alternative positions

This reuses exactly the same Geant4 and detector-response events:

```bash
base=LFGD_Recon_Simu/output/homo_student_lfg-default-position_proton1GeV_standard_default; LFGD_Recon_Simu/rerun_reconstruction.sh "$base/detresponse.root" LFGD_Recon_Simu/output/homo_proton1GeV_standard_dominant_axis -O "par_override=$PWD/LFGD_Recon_Simu/Studies/ProtonHitRecoStudy/dominant_axis_position.parameters.dat"
```

The standard fitter results are stored in `track_nodes` in both flat files.
The analysis notebook additionally performs an unweighted orthogonal/PCA line
fit directly to each set of reconstructed 3D hits.

To compare the one-point-per-column construction with the baseline using the
official track fitter (rather than the notebook's diagnostic PCA line), rerun
reconstruction with `column_official.parameters.dat`.  Only the hit-building
mode changes; hit selection and the complete downstream official tracking
chain remain enabled and unchanged.

## Analysis

Open `ProtonHitRecoComparison.ipynb`.  Its first code cell contains the two
flat-file paths and the event number.  It provides:

- event-level and aggregate counts of used 2D hits per view, reconstructed 3D
  hits, and reconstructed tracks;
- reconstructed 2D positions separately for XZ, YZ, and XY views;
- overlaid default and dominant-axis 3D hits with proton MC truth;
- standard-fitter track nodes and simple PCA lines;
- aggregate one-dimensional position distributions;
- hit multiplicity and charge distributions;
- per-coordinate reconstructed-minus-truth residuals, matched by virtual
  cube ID;
- reconstructed charge versus true deposited energy and a data-derived
  charge-per-energy residual.

The flat-tree extension adding `view_x`, `view_y`, and `view_z` requires the
flat tree to be regenerated.  Older flat files remain readable, but cannot
provide the new 2D-position plots.

The `hits2d` flat tree stores both used and unused reconstructed 2D peak/CoG
hits. Unused hits carry a rejection code and text reason: low charge, no
spatial three-view match, missing coordinate information, outside the voxel
grid, timing mismatch, compatible but not selected, or another combination
failure. The notebook reports these reasons per event and over the full
sample, and marks unused hits with crosses in the 2D event display.
