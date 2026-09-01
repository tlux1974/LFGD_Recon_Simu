# Reconstruct identical events with several light maps

The comparison starts from one existing Geant4 `g4.root`. This is stronger
than regenerating with the same seed: every light-map branch reads the exact
same Geant4 events, trajectories, energy deposits, and geometry. Only detector
response and its downstream reconstruction are rerun.

Inside the writable ND280 sandbox container, run this as one line:

```bash
./compare_lightmaps.sh /absolute/path/to/common/g4.root output/lightmap_comparison
```

With no additional arguments, all ROOT files directly inside
`Studies/PICLightMap/lightmaps/` are used, except the invalid octant-mirrored
maps. To select maps explicitly:

```bash
./compare_lightmaps.sh /absolute/path/to/common/g4.root output/lightmap_comparison Studies/PICLightMap/lightmaps/map_a.root Studies/PICLightMap/lightmaps/map_b.root
```

Each map gets an independent output directory containing detector response,
standard HFG reconstruction, flat tree, GlobalFit output, and logs.
`manifest.csv` records the
absolute input paths and SHA-256 checksums of the common Geant4 file and every
map. The runner validates that every selected map contains 970 response bins.

ROOT production is deliberately separate from analysis. After production,
leave the container, activate a Python environment containing Matplotlib and
PyROOT, and run this as one line:

```bash
./analyse_lightmap_study.sh /absolute/path/to/output/lightmap_comparison
```

This produces the reconstruction diagnostics, GlobalFit analysis, and the
cross-map `comparison_summary.csv` and `comparison_summary.png`. They compare
mean matched-voxel efficiency, mean off-track reconstructed voxels, GlobalFit
convergence, and GlobalFit direction error across all maps.

The default GlobalFit settings are `EVENT=all TREE=fiber_hits MIN_CHARGE=10
DBSCAN=1 DBSCAN_EPSILON_MM=14.2 DBSCAN_MIN_POINTS=2`. Override the complete
option list with `GLOBAL_FIT_OPTIONS`. The common detector-response seed is
controlled by `SEED` and defaults to 12345.

Do not compare a branch produced with a different `g4.root`; the manifest
checksum is the proof that truth events were held fixed.

## Current x-axis cone study

The selected common sample contains 1,000 events:

```text
output/homo_student_lfg-best_xaxis_cone/g4.root
```

Results go under `Studies/LightmapsStudie/`, with one subdirectory per map.
Use the `lfg-best` reconstruction parameters consistently for every branch:

```bash
HFGRECON_PARAMETER_FILE=/home/tlux/HK/ND280++/LFGD_Recon_Simu/student_lfg_best.parameters.dat /home/tlux/HK/ND280++/LFGD_Recon_Simu/compare_lightmaps.sh /home/tlux/HK/ND280++/LFGD_Recon_Simu/output/homo_student_lfg-best_xaxis_cone/g4.root /home/tlux/HK/ND280++/LFGD_Recon_Simu/Studies/LightmapsStudie
```

## Five-GeV low-multiple-scattering study

Generate one common 1,000-event sample inside the configured writable ND280
container. This command stops after producing `g4.root`:

```bash
DIRECTION_MODE=cone DIRECTION="1 0 0" CONE_HALF_ANGLE_DEG=5 ./generate_common_g4.sh 5000 1000 Studies/LightmapsStudie/5GeV/common
```

Then run every full light map against exactly that common sample:

```bash
HFGRECON_PARAMETER_FILE=/home/tlux/HK/ND280++/LFGD_Recon_Simu/student_lfg_best.parameters.dat ./compare_lightmaps.sh Studies/LightmapsStudie/5GeV/common/g4.root Studies/LightmapsStudie/5GeV/maps
```

Production and analysis remain separate. Outside the container, in the
Jupyter environment containing PyROOT and Matplotlib, analyse the result:

```bash
bash ./analyse_lightmap_study.sh Studies/LightmapsStudie/5GeV/maps
```

If production is interrupted, rebuild or correct the failed executable and
resume with the same inputs and output directory by prefixing the production
command with `RESUME=1`. The launcher validates existing ROOT files, reuses
completed detector-response, reconstruction, flat-tree, and GlobalFit stages,
and continues at the first missing stage.
