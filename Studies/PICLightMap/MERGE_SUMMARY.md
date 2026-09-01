# PIC light-map merge summary (2026-08-12)

## Run the full-map merger

The scattering length is mandatory and is specified in millimetres. From
`Studies/PICLightMap`, run this as one line, changing `0.7` to the production
setting:

```bash
./merge_lightmap.sh --scat-length 0.7 --fibre-diameter 1.0 --scintillator-ri 1.48 --scattering-model legacy_rayleigh --mie-g 0.5 --mie-forward-fraction 0.75
```

There is deliberately no default. The merger converts the value to the chunk
filename tag (for example, `0.7` becomes `0p7mm`), selects only those ten
matching chunk ROOT files, and writes
`homo_response_10bin_5m_0p7mm_100kPhotons.root` plus its `_efficiency.csv`.
This prevents an output filename for one scattering length from accidentally
being filled with chunks produced at another setting.

To choose explicit output paths:

```bash
./merge_lightmap.sh --scat-length 0.7 --fibre-diameter 1.0 --output lightmaps/homo_response_10bin_5m_0p7mm_100kPhotons.root --csv lightmaps/homo_response_10bin_5m_0p7mm_100kPhotons_efficiency.csv
```

## Outputs

- `homo_response_10bin_5m_0p25mm_100kPhotons.root`: merged response map from
  the ten production chunks (`0000-0099` through `0900-0999`).
- `homo_response_10bin_5m_0p25mm_100kPhotons_efficiency.csv`: per-position
  fibre-hit efficiency for the merged map.
- `original_1mm_efficiency.csv`: the same efficiency extraction from the
  installed original/default 1 mm-scattering map.

## Valid-position check

The merger read all 1,000 simulated events. Both the new and original maps
contain the same 970 populated positions. The same 30 grid points are absent
from both maps because those source positions fall on fibres rather than in
scintillator. Thus the populated-bin mask is exactly consistent with the
original map (970 common, zero present in only one map).

## Efficiency comparison with the original map

The files have different scattering lengths, so identical efficiencies are
not expected: the new map uses 0.25 mm and the original uses 1 mm.

| Quantity | New 0.25 mm | Original 1 mm |
|---|---:|---:|
| Mean efficiency | 0.81314857 | 0.93987470 |
| Population standard deviation | 0.00191900 | 0.00128232 |
| Minimum | 0.80642 | 0.93603 |
| Maximum | 0.81744 | 0.94280 |

Across matched positions, the mean new/original ratio is 0.86516705, its
standard deviation is 0.00169922, and its range is 0.85905425-0.87105763.
The per-position Pearson correlation is 0.55545.

## Octant-only consistency check

Octants are defined relative to the cube centre `(0, 30, 910) mm`. Counts are
120 or 125 because the 30 excluded fibre points lie on specific octants.

| x sign | y sign | z sign | Positions | New mean | Original mean | Ratio |
|---:|---:|---:|---:|---:|---:|---:|
| - | - | - | 125 | 0.81357792 | 0.94049936 | 0.86504888 |
| - | - | + | 120 | 0.81318442 | 0.93970617 | 0.86536031 |
| - | + | - | 120 | 0.81291208 | 0.93965692 | 0.86511584 |
| - | + | + | 120 | 0.81290933 | 0.93965717 | 0.86511268 |
| + | - | - | 120 | 0.81300808 | 0.93971767 | 0.86516207 |
| + | - | + | 120 | 0.81316467 | 0.93970950 | 0.86533622 |
| + | + | - | 120 | 0.81290292 | 0.93957617 | 0.86518044 |
| + | + | + | 125 | 0.81349672 | 0.94042568 | 0.86503031 |

The octant mean ratios span only 0.86503031-0.86536031 (0.038% relative
spread), so there is no meaningful octant-dependent inconsistency.

## Code and workflow changes

- Restricted `merge_lightmap.sh` input discovery to the ten production
  `chunk_*/results/*.root` files, excluding local two-position validation
  files.
- Updated `MakeHomoPhotonMap.cc` so hitless fibre-source positions are skipped,
  matching the original merger behavior, while valid positions remain in the
  new `position_efficiency` tree.
- Extended `ExportLightMapEfficiency.C` to export both new maps containing the
  efficiency tree and legacy/original maps containing only light-fraction
  objects.
- Rebuilt `MAKEHOMOPHOTONMAP.exe` and completed the full merge in the ND280++
  container environment.

The original/default ROOT file was read only and was not modified.
