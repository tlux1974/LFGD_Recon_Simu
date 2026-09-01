# Bug report: cumulative memory growth in HOMO light-map response

## Summary

`DETRESPONSESIM.exe` showed approximately linear resident-memory growth while
processing HOMO events with the compact light-map response. The process was
repeatedly killed by the operating system shortly after event 600 in a
1,000-event sample.

The problem was present in the original upstream HOMO light-map
implementation. It was not introduced by the local virtual-fibre geometry or
fast aggregate-response changes.

## Observed behaviour

The failure was reproducible with all of the following configurations:

- detailed simulation with individual `TWLSPhoton` objects;
- fast aggregate fibre response;
- HOMO fibre attenuation enabled or disabled;
- `homo_raw` and `homo_truth` enabled or disabled.

Events 600--604 completed normally when processed alone in a new process using
the event-loop `-s 600 -n 5` options. This ruled out a single pathological
event and demonstrated that the failure was cumulative.

With fast response and raw/truth collections disabled, the built-in oaEvent
memory logger measured:

| Event | Resident memory | Virtual memory |
|------:|----------------:|---------------:|
| 0 | 0.333 GB | 0.738 GB |
| 1 | 1.059 GB | 1.437 GB |
| 9 | 1.134 GB | 1.512 GB |
| 19 | 1.264 GB | 1.641 GB |
| 29 | 1.384 GB | 1.761 GB |
| 39 | 1.506 GB | 1.884 GB |
| 49 | 1.640 GB | 2.017 GB |

After geometry initialization, resident memory therefore increased by about
12.3 MB per event. Extrapolation gives roughly 8.4 GB near event 600,
consistent with the observed operating-system kill.

## Root cause

For every Geant4 hit segment, `THomoScint::GetLightMapValues()` called
`TDirectory::Get()` for both compact light-map objects:

```cpp
photonPositions = fphotonPositionDir->Get<TObjArray>(...);
lightFractions = fphotonFractionDir->Get<std::vector<double>>(...);
```

Repeated ROOT deserialization, particularly of the STL vector object, left
allocations resident after each lookup. Since the lookup occurs for every
Geant4 segment, memory increased with every event even when individual photon
simulation and diagnostic hit collections were disabled.

There was a second bug in the same code path:

```cpp
TVector3* hitPos = static_cast<TVector3*>(photonPositions->At(hitId));
*hitPos += offset;
```

This modified the light-map position object itself. Reusing a loaded map bin
could therefore accumulate detector offsets and corrupt subsequent fibre
positions.

## Git provenance

`git blame` attributes both the repeated compact-map reads and the in-place
position modification to upstream commit:

```text
22c614d  Merge branch 'feature/homo-hit-timing-info' into 'ND280-plusplus'
         Author: Ewan Miller
         Date: 2025-10-15
```

The local virtual-fibre commit `2677c31` changed the geometry lookup to
`GetFiberGeomIdFuzzy()` but retained the pre-existing light-map reads and
in-place position modification. It did not introduce either bug.

## Fix

The compact map is now read once during `THomoScint::InitModel()` and stored in
two caches indexed by the light-map global bin:

- `compactPositionCacheMap` for `TObjArray` positions;
- `lightFractionCacheMap` for light-fraction vectors.

`GetLightMapValues()` now returns pointers from these caches instead of calling
ROOT `Get()` for every segment.

Mapped positions are translated into a temporary value:

```cpp
const TVector3 hitPos = *mapHitPos + offset;
```

The cached light-map object is no longer modified.

Relevant implementation files are:

- `SoftProj/detResponseSim/inc/THomoScint.hxx`
- `SoftProj/detResponseSim/src/THomoScint.cxx`

## Verification

The fixed DetResponseSim package built successfully. Repeating the same
50-event fast/no-truth memory test gave:

```text
Maximum resident memory usage: 1.0523 GB
Maximum virtual memory usage:  1.4297 GB
Total Events Read:             50
Total Events Written:          50
```

The previous linear rise to 1.64 GB by event 50 was absent. The one-time cache
contains 970 populated compact light-map bins.

## Related but separate changes

The investigation also added:

- optional fast aggregate HOMO fibre response;
- explicit release of transient photon and avalanche data after each event;
- flat-tree output for `homo_raw` and `homo_truth`.

These changes reduce peak memory or expose diagnostic data, but they were not
the fix for the cumulative 12.3 MB/event leak described here.

## Separate general reconstruction bug found during validation

After the light-map correction, one LFGD event exposed an independent numeric
problem in `THFGStochTrackFit`. This fitter is shared by HFGD and LFGD, so the
bug was not detector-specific.

The stochastic prior estimated deposited charge per unit length using only
the separation of the first and last prior measurements:

```cpp
double length = (firstPosition - lastPosition).Mag();
avgEDep /= length;
```

Repeated endpoints or a folded track can have zero endpoint separation even
when several nodes are present. This generated an infinite energy-deposition
state (`S[0] = inf`) and terminated the entire reconstruction event loop with
`Numeric problem`.

The general fitter was made safe by:

- calculating sampled path length as the sum of consecutive node distances;
- requiring finite, nonzero path length and spatial extent;
- selecting spatially distinct measurements to initialize track direction,
  with the farthest pair as a deterministic fallback;
- returning an invalid fit handle for a degenerate candidate instead of
  throwing an exception that aborts the complete event loop.

The correction is in:

- `SoftProj/hfgrecon/src/THFGStochTrackFit.cxx`

The corrected `hfgrecon` package built successfully. The previously failing
input event (event index 231) was then reconstructed alone with the same
`lfg-best` parameter file. It completed with one event read and one event
written, without an infinite sample or numeric exception.

## Future key improvement: preserve fractional HOMO fibre gaps

The HOMO geometry stores the three fibre gaps as `double`, and its UI commands
correctly accept non-integral lengths. However, the corresponding getters in
`ND280HomoBoxConstructor.hh` return `int`:

```cpp
void SetGapX(double gap) { fGapX = gap; }
int GetGapX() { return fGapX; }
```

`GetGapY()` and `GetGapZ()` have the same defect. A fractional gap is therefore
silently truncated whenever the fibre pitch or placement is calculated.

This was exposed by the 1.5 mm-diameter fibre light-map study. Keeping the
standard 10 mm centre-to-centre pitch requires a radius of 0.75 mm and an
edge-to-edge gap of 8.5 mm. The getter truncated that gap to 8 mm, producing a
9.5 mm pitch. The number of fibres consequently increased (for example, 199
became 209 in one transverse direction), and ROOT geometry voxelisation
terminated with `std::bad_alloc` before the first simulated position. The
resulting 430--434 byte ROOT files are incomplete crash artifacts. Increasing
the job memory would not make this production valid because the geometry pitch
would still be wrong.

The required correction, now applied in the local source, is:

```cpp
double GetGapX() { return fGapX; }
double GetGapY() { return fGapY; }
double GetGapZ() { return fGapZ; }
```

in:

- `SoftProj/nd280Geant4Sim/inc/nd280-plusplus/homo/ND280HomoBoxConstructor.hh`

The local `ND280GEANT4SIM.exe` and `MAKEHOMOPHOTONMAP.exe` were rebuilt
successfully after this change. A PIC production still requires these local
changes to be installed in a new SIF; the historical SIF does not acquire them
from the host checkout automatically.

After rebuilding `nd280Geant4Sim`, geometry validation should explicitly check
both the requested fibre diameter and the resulting 10 mm fibre pitch before a
large light-map production is submitted. Integral-gap configurations do not
trigger the truncation: the available 2.0 mm-diameter test chunk used a 1.0 mm
radius plus an 8 mm gap, retained the expected fibre counts, and completed all
100 source positions successfully.

### Fibre-diameter-dependent source-position denominator

The historical map maker and analyser effectively retained the 970 valid-bin
convention of the 1.0 mm fibre geometry for every fibre diameter. All four
existing 2.0 mm ND280 maps were confirmed to contain 970
`position_efficiency` entries, whereas the same 10x10x10 source-centre grid
contains only 850 scintillator positions for 2.0 mm fibres. This biased
diameter comparisons by including source centres geometrically occupied by
the wider fibres.

`MAKEHOMOPHOTONMAP.exe` now accepts `-O fibreDiameter=...` and excludes source
centres inside any of the three HOMO fibre cylinders before filling the map.
`merge_lightmap.sh` passes the requested diameter and validates the resulting
entry count against the geometry (970 for 1.0 mm and 850 for 2.0 mm). The
analyser no longer hard-codes 970 and reports the actual valid-position count.
Existing Geant4 chunks can be re-merged with this correction; they do not need
to be simulated again for this particular bug.

## Light-map photon undercount from merged hit segments

The historical `MAKEHOMOPHOTONMAP.exe` counted one detected photon for every
stored `TG4HitSegment`:

```cpp
fiberHitsMap[geomId]++;
```

This is not one segment per photon. `ND280SegmentSD` deliberately merges
different tracks when their energy deposits occur in the same fibre, within
the configured 1 mm spatial tolerance and within 1 ns. `ND280HitSegment`
preserves all track IDs in `GetContributors()` and sums their deposited energy,
but the map builder discarded that multiplicity. A 100,000-photon bomb is
therefore undercounted, especially when a short scattering length concentrates
many photons into the same fibres and nearby positions.

Direct checks of existing raw chunks found:

| Sample | Stored segments | Contributor tracks | Deposited energy / 1 eV |
|---|---:|---:|---:|
| 0.4 mm scattering, scintillator RI 1.36 | 87,007 | 98,173 | 98,173 |
| 0.85 mm scattering | 93,220 | 98,847 | 98,847 |

Thus the old map reported efficiencies of 87.007% and 93.220%, whereas the
actual fractions absorbed in fibres were 98.173% and 98.847%. This accounting
bug also suppresses the response of the brightest fibres and biases apparent
scattering-length trends, nearest-fibre fractions, and view-to-view balance.

The local map builder now uses:

```cpp
fiberHitsMap[geomId] += seg->GetContributors().size();
```

It also checks contributor multiplicity against the deposited energy expected
from the monoenergetic 1 eV photon bombs. Existing Geant4 chunk files retain
the contributor lists, so corrected maps can be produced by rebuilding the map
maker and re-merging the existing chunks; the expensive Geant4 simulations do
not need to be repeated.

## Historical Mie-scattering typo

The original `OpticalLiquidO` material registered the intended Geant4 Mie
scattering length under:

```cpp
"MIEGH"
```

The Geant4 material-property key is `MIEHG`. Consequently the Mie process was
inactive and the historical light-map production used only its simultaneously
configured `RAYLEIGH` process. The local implementation preserves that exact
behavior under the explicit `legacy_rayleigh` model and uses the correctly
spelled `MIEHG` property only for the separately selectable `double_hg` model.
The standalone isotropic scan shows that the angular model does not explain
the large historical efficiency deficit; the merged-segment accounting above
does.
