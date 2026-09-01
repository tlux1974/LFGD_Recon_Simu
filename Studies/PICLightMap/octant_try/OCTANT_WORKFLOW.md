# One-octant light-map production and analysis

> **Study only — do not use octant-mirrored maps for detector response or
> reconstruction.** This workflow was an investigation of whether one octant
> could be simulated and reflected into a full light map. The test showed that
> simple reflection about the cube centre is not valid because the physical
> fibre lattice is offset from that centre. It can place mirrored photon hits
> at coordinates where no fibre exists, causing detector response to discard
> most of the light. For production and quick scattering-length comparisons,
> simulate all valid positions and merge them without mirroring. To reduce the
> runtime of a quick cross-check, lower the number of photons per position
> instead of reducing the spatial coverage.

## Submit the octant

Copy these files to the PIC working directory together with the existing SIF:

```bash
scp positions.txt run_lightmap_chunk.sh run_lightmap_octant.sh \
    lightmap_octant.sub ValidateHomoGeometry.C \
    tlux@ui.pic.es:/nfs/pic.es/user/t/tlux/tlux/PICLightMap/
```

On PIC:

```bash
cd /nfs/pic.es/user/t/tlux/tlux/PICLightMap
chmod +x run_lightmap_chunk.sh run_lightmap_octant.sh
mkdir -p logs
condor_submit -dry-run octant.dryrun lightmap_octant.sub
condor_submit lightmap_octant.sub
```

This submits five jobs of 25 positions, covering all 125 positions in the
`+++` octant. Results are returned under `octant_chunk_0/`,
`octant_chunk_25/`, ..., `octant_chunk_100/`. The production parameters are
100,000 photons per position, 0.7 mm scattering length, and 5 m absorption
length.

The wrapper runs in HTCondor's `_CONDOR_SCRATCH_DIR`, ensuring that the five
declared output directories are available to Condor's output-transfer stage.

## Build the full map by mirroring

Copy all five `octant_chunk_*` directories back into this directory. In the
configured ND280++ container environment, run:

```bash
./merge_lightmap_octant.sh
```

The merger reflects each simulated response into all eight octants. It also
uses the original/default map as a geometry mask. This is necessary because
the 30 fibre-source positions are not octant-symmetric. The script refuses the
result unless the final efficiency CSV contains exactly 970 valid positions.

The default outputs are:

```text
homo_response_octant_mirrored_10bin_5m_0p7mm_100kPhotons.root
homo_response_octant_mirrored_10bin_5m_0p7mm_100kPhotons_efficiency.csv
```

Alternative output, CSV, or mask paths can be supplied as positional
arguments; see the variable assignments at the top of the script.

## Analyse a setting

The analyser is general and is not specific to octant production. Its inputs,
statistics, plots, command-line usage, notebook instructions, and troubleshooting
are documented separately in `LIGHTMAP_ANALYSER.md`.
