# HFGD 1 GeV isotropic-proton hit reconstruction

This isolated study uses the standard local HFGD reconstruction without
changing its hit builder, charge sharing, clustering, or track fitter. In the
standard HFGD chain, a calibrated fibre measurement is the 2D hit supplied to
the cube builder; there is no HOMO-style peak/CoG clustering stage.

The generated GPS macro contains `/t2k/field 0.0 tesla`, so the magnetic field
is off.

## Build

Inside the normal ND280++ container:

```bash
cd /home/tlux/HK/ND280++/LFGD_Recon_Simu && ./build_flat_treemaker.sh
```

The local `hfgrecon` package must also be rebuilt after changing its source.

## Produce the sample

```bash
PARTICLE=proton ENERGY_MEV=1000 DIRECTION_MODE=isotropic LFGD_Recon_Simu/run_student_sample.sh hfg-standard 1000 proton1GeV_standard
```

The output is
`LFGD_Recon_Simu/output/hfg_student_hfg-standard_proton1GeV_standard`.

## Analysis

Open `HFGDHitReco.ipynb`. For a selectable event it overlays:

- primary-proton MC truth in green;
- used HFGD 2D fibre hits as blue open circles;
- unused 2D fibre hits as red crosses;
- reconstructed 3D cube hits as orange squares;
- standard fitted tracks and a simple PCA line in the separate 3D display.

It also reports 2D hits per view, unused hits, 3D hits, reconstructed tracks,
position distributions, and 3D-position/charge residuals against MC truth.

For HFGD an unused 2D hit currently has the reason `not_used_in_3d`. More
detailed failure reasons would require instrumenting the standard exact fibre
intersection and timing decisions inside `THFGHits3D`; the reconstruction
itself is intentionally left unchanged for this baseline study.
