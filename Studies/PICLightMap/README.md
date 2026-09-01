# PIC light-map production

This directory generates a complete 10×10×10 HOMO light map as either ten
jobs of 100 positions or twenty jobs of 50 positions. The normal workflow uses
`lightmap_chunks.sub` and the separately built `ND280Reco_doubleHG_fixed.sif`.

## 1. Install or update the PIC working directory

The submit files use this fixed directory:

```text
/nfs/pic.es/user/t/tlux/tlux/PICLightMap
```

From the local `Studies/PICLightMap` directory, copy the small production
files using the same PIC login host for which your SSH key is configured:

```bash
scp README.md DOUBLE_HG_PRODUCTION.md lightmap_chunks.sub chunks_50.txt chunks_100.txt run_lightmap_chunk.sh merge_lightmap.sh ValidateHomoGeometry.C ExportLightMapEfficiency.C positions.txt tlux@ui.pic.es:/nfs/pic.es/user/t/tlux/tlux/PICLightMap/
```

Copy the approximately 4.7 GB SIF only if it is not already present:

```bash
scp /home/tlux/HK/ND280++/ND280Reco_doubleHG_fixed.sif tlux@ui.pic.es:/nfs/pic.es/user/t/tlux/tlux/PICLightMap/
```

On PIC:

```bash
cd /nfs/pic.es/user/t/tlux/tlux/PICLightMap
chmod +x run_lightmap_chunk.sh merge_lightmap.sh
mkdir -p logs
test -r ND280Reco_doubleHG_fixed.sif && echo 'SIF path is readable' || echo 'FIX THE SIF PATH BEFORE SUBMITTING'
```

Do not put the SIF in `transfer_input_files`; all jobs use the single copy on
the shared filesystem.

## 2. Submit a complete production

`lightmap_chunks.sub` defaults to ten 100-position jobs. Override
`CHUNK_SIZE = 50` to submit twenty 50-position jobs. Each job starts Geant4
once and produces one ROOT fragment. The worker refuses sizes outside 50--100.

The existing `lightmap_10x100.sub` remains the legacy-Rayleigh comparison
route. It uses `ND280Reco_03.sif`, defaults to `legacy_rayleigh`, and its
generated macro deliberately contains none of the new double-HG commands.
Keep it alongside, rather than replacing it with, `lightmap_chunks.sub` on PIC.

Set the scattering length, fibre diameter, and scintillator refractive index
with independent overrides:

- `SCAT_VALUE` is the numerical scattering length in mm;
- `FIBRE_DIAMETER_MM` is the physical fibre diameter in mm;
- `SCINTILLATOR_RI` is the `OpticalLiquidO` refractive index (standard value
  `1.48`). The existing fibre remains unchanged at `n=1.3`.

The worker derives all directory and filename tags from these values. It sets
the HOMO fibre radius to half the diameter and the edge-to-edge gap to
`10 mm - diameter`, preserving the existing 10 mm fibre pitch. Keep
`SCAT_UNIT=mm`, `NPHOTONS=100000`, and `ABS_VALUE=5` unless intentionally
starting a different study.

Dry-run the nominal double-HG model first with 100 positions per job:

```bash
condor_submit -dry-run double-hg-100.dryrun lightmap_chunks.sub
```

Submit it after checking the dry-run. If the first timing indicates that a
100-position job is too slow, use the second command instead:

```bash
condor_submit lightmap_chunks.sub
condor_submit -append 'CHUNK_SIZE = 50' lightmap_chunks.sub
```

Do not submit both chunk sizes for the same physical setting: that would
deliberately produce overlapping position ranges. Override physical settings
with `-append`, for example as a single command:

```bash
condor_submit -append 'CHUNK_SIZE = 50' -append 'SCAT_VALUE = 0.25' -append 'MIE_G = 0.5' -append 'MIE_FORWARD_FRACTION = 0.75' lightmap_chunks.sub
```

The scattering model, Mie mean free path, g, forward fraction, fibre diameter,
and scintillator index are all encoded in the output path.

## 3. Output layout

Productions with different scattering lengths cannot collide:

```text
fibre_1.0mm/
  model_double_hg/
    scat_0.30mm/
      g_0.5_rfb_0.75/
        scin_ri_1.48/
          chunk_0/{results,macros,simulation.log,position_timing.csv}
          chunk_100/{results,macros,simulation.log,position_timing.csv}
          ...
          chunk_900/{results,macros,simulation.log,position_timing.csv}
```

Log filenames also contain the tag, for example
`logs/fibre-1.0mm.double_hg.scat-0.30mm.g-0.5.rfb-0.75.chunk.100.CLUSTER.PROCESS.out`.

List the returned ROOT fragments for one production with:

```bash
find fibre_1.0mm/model_double_hg/scat_0.30mm/g_0.5_rfb_0.75/scin_ri_1.48/chunk_* -path '*/results/*.root' -type f -ls
```

Each timing CSV records the global position index, local run number, source
coordinates, and Geant4 user/real/system time. Show the twenty slowest 0.8 mm
positions with:

```bash
{ head -1 fibre_1.0mm/model_double_hg/scat_0.30mm/g_0.5_rfb_0.75/scin_ri_1.48/chunk_0/position_timing.csv; tail -q -n +2 fibre_1.0mm/model_double_hg/scat_0.30mm/g_0.5_rfb_0.75/scin_ri_1.48/chunk_*/position_timing.csv | sort -t, -k7,7nr | head -20; }
```

## 4. Monitor and diagnose

```bash
condor_q
condor_q -hold
```

Inspect a running job and cancel a cluster with:

```bash
condor_tail CLUSTER.PROCESS
condor_rm CLUSTER
```

Show detailed hold reasons with:

```bash
condor_q -const 'JobStatus == 5' -af ClusterId ProcId HoldReasonCode HoldReasonSubCode HoldReason
```

A successful worker log ends with `VALIDATION PASSED` and lists a non-empty
ROOT file. The validator rejects output without HOMO geometry, with HFG still
enabled, or without serialized optical-fibre hits.

## 5. Merge one complete production

Run the merger in the configured local/sandbox environment where
`MAKEHOMOPHOTONMAP.exe` and ROOT are available. The scattering length is
mandatory and selects only its tagged directory:

```bash
./merge_lightmap.sh --scat-length 0.30 --fibre-diameter 1.0 --scintillator-ri 1.48 --scattering-model double_hg --mie-g 0.5 --mie-forward-fraction 0.75
```

It accepts either 10--20 matching files, reads their position intervals from
the filenames, and refuses gaps, overlaps, chunks outside 50--100 positions,
or coverage other than exactly `0..999`. It excludes source centres inside the
requested fibre diameter and validates the geometry-dependent number of
scintillator positions (970 for 1.0 mm and 850 for 2.0 mm). It writes a normal HOMO response ROOT
file plus its efficiency CSV. The ROOT format is unchanged and can therefore
be used by the ND280 detector response and as the explicit light-map argument
to the global fit.

The default output basename includes `double_hg`, `0p30mm`, `g-0p5`,
`rfb-0p75`, the fibre diameter, and the scintillator refractive index.

The aggregation uses `cubeCentre="0 30 910 mm"`, `cubeSide="10 mm"`,
`nBins=10`, and `nPhotons=100000`.

## Worker safeguards

The submit file sets `LIGHTMAP_DETECTOR=HOMO`. The generated macro explicitly
disables HFG and enables HOMO before `/t2k/update`, enables full optical-photon
simulation, and uses 1 eV primary photons. The worker sources the ND280 pilot
and master setup layers, checks that `ND280GEANT4SIM.exe` is available, and
runs `ValidateHomoGeometry.C` before accepting the result.
