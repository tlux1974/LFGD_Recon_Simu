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
