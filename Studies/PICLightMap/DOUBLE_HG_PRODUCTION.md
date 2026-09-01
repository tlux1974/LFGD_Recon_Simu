# Double-HG ND280 light-map production

## Definitions

For `double_hg`, `SCAT_VALUE` is the true Mie scattering mean free path
`l_s^Mie = 1/mu_s`: the exponential mean distance between physical scattering
events. It is not the transport mean free path. `MIE_FORWARD_FRACTION` is
`F/(F+B)`, not `F/B`.

The nominal starting point is:

```text
SCATTERING_MODEL = double_hg
SCAT_VALUE = 0.30
MIE_G = 0.5
MIE_FORWARD_FRACTION = 0.75
```

The first-moment anisotropy is
`g_eff = (2*MIE_FORWARD_FRACTION-1)*MIE_G = 0.25`. These are working guesses;
when the oLS producer provides `mu_s`, use `SCAT_VALUE = 1/mu_s` with consistent
units instead of converting a phenomenological isotropic length.

## Local source and commands

The local `nd280Geant4Sim` adds these pre-initialization commands:

```text
/t2k/detector/liquidOScatteringModel legacy_rayleigh
/t2k/detector/liquidOScatteringModel double_hg
/t2k/detector/liquidOMieScatLen 0.30 mm
/t2k/detector/liquidOMieG 0.5
/t2k/detector/liquidOMieForwardFraction 0.75
```

`legacy_rayleigh` is the default and preserves existing macros. In
`double_hg` mode only `MIEHG` is configured; Rayleigh is not configured at the
same time.

## Rebuild locally

Run from the workspace as one line:

```bash
singularity exec ND280ppCont bash -lc 'source /usr/local/t2k/current/nd280SoftwarePilot/nd280SoftwarePilot.profile >/dev/null 2>&1; source /usr/local/t2k/current/nd280SoftwareMaster_14.36-plusplus.0.3/bin/setup.sh >/dev/null 2>&1; cd /home/tlux/HK/ND280++/SoftProj/nd280Geant4Sim; make -C Linux-AlmaLinux_9.5-gcc_11-x86_64 -j2'
```

## Build a new SIF

The helper keeps `ND280Reco_03.sif` unchanged. It installs the local executable
and library into the existing writable `ND280ppCont` development sandbox,
keeping one `*.pre-double-hg` backup of each displaced file, and converts the
sandbox to a new SIF:

```bash
cd /home/tlux/HK/ND280++/LFGD_Recon_Simu/Studies/PICLightMap && ./build_double_hg_sif.sh /home/tlux/HK/ND280++/ND280Reco_doubleHG_fixed.sif
```

Do not run the helper until the local validation is satisfactory. Building the
SIF needs roughly the size of the existing image in free disk space.

Copy the new SIF and production files to PIC. `lightmap_chunks.sub` already
uses its cluster path. Required files are:

```text
ND280Reco_doubleHG_fixed.sif
run_lightmap_chunk.sh
lightmap_chunks.sub
chunks_50.txt
chunks_100.txt
positions.txt
ValidateHomoGeometry.C
```

Before submission, use the existing dry-run procedure and inspect the emitted
macro. It must contain `liquidOScatteringModel double_hg`,
`liquidOMieScatLen`, `liquidOMieG`, and `liquidOMieForwardFraction` before
`/t2k/update`.

## Merge

After all 10 or 20 chunks return, run locally in the ND280 environment:

```bash
./merge_lightmap.sh --scat-length 0.30 --fibre-diameter 1.0 --scintillator-ri 1.48 --scattering-model double_hg --mie-g 0.5 --mie-forward-fraction 0.75
```

The model, true Mie mean free path, `g`, and fractional forward weight are all
included in production directories, chunk filenames, logs, and merged-map
filenames so incompatible productions cannot be mixed silently.
