CONTAINER/PACKAGE CHANGE: nd280Geant4Sim
=======================================

Date: 2026-08-04
Submission scope: nd280Geant4Sim only

This change is separate from the local LFGD_Recon_Simu changes and should be
submitted later as an independent container/package change.


File changed
------------

Container path:
  /usr/local/t2k/current/nd280Geant4Sim_7.17-plusplus.1/
  src/ND280UserDetectorConstruction.cc

Workspace-backed path:
  ND280ppCont/usr/local/t2k/current/nd280Geant4Sim_7.17-plusplus.1/
  src/ND280UserDetectorConstruction.cc


Problem 1: missing ROOT geometry color
--------------------------------------

Observed error:

  Missing color for "OpticalFiberCore" in volume RotatedFiber/Fiber

FiberCore was registered with ND280RootGeometryManager, but OpticalFiberCore
was not. After constructing OpticalFiberCore, the following registration was
added:

  gMan->SetDrawAtt(opticalFiberCore,kCyan-4);

This gives OpticalFiberCore the same ROOT geometry color as FiberCore.


Problem 2: incorrect OpticalFiberCore density
---------------------------------------------

OpticalFiberCore is intended to be a copy of FiberCore. However, the mutable
local variable named density had already been changed to the fluorinated
cladding density (1.43 g/cm3) before OpticalFiberCore was constructed.

Construction was changed from using that stale density variable to:

  fiberCore->GetDensity()

This gives OpticalFiberCore the intended fibre-core density of 1.05 g/cm3 and
prevents an incorrect fibre mass in the detector geometry.


Configure and rebuild inside the container
------------------------------------------

The installed output directory may contain a generated Makefile but no
CMakeCache.txt. In that case, running only "cmake --build" reports:

  Error: could not load cache

Configure it first, then build the Geant4 simulation executable and its
library dependency:

  package=/usr/local/t2k/current/nd280Geant4Sim_7.17-plusplus.1

  cmake -S "$package/cmake" -B "$package/$ND280_SYSTEM"

  cmake --build "$package/$ND280_SYSTEM" \
    --target ND280GEANT4SIM.exe -j4

After rebuilding, rerun the Geant4 simulation. Do not reuse g4.root from the
failed run because it may be absent or incomplete.
