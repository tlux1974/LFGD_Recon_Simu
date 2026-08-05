LOCAL LFGD_Recon_Simu CHANGES
=============================

Date: 2026-08-04
Scope: /home/tlux/HK/ND280++/LFGD_Recon_Simu

Update 2026-08-05: the local hfgRecon default now uses
`VoxelPositionMode.homo = 2`.  It assigns each continuous reconstructed hit a
virtual-cube identity and resolves multiple candidates in that cube without
moving the selected hit to the cube centre.  The ten-event validation is in
`output/homo_mu700_virtual_geometry_100ev_20260805_145407/voxel_dedup_test10b`.
It reduced 3974 3D-hit entries in 2420 occupied voxels to 1997 entries in 1997
voxels, while retaining continuous reconstructed tracks.

The event overlays now use only exact `mc_virtual_segments` truth.  Sparse
Geant4 trajectory points are no longer drawn or required by the plotting
script.  The charge-sum plot is labelled as measured fibre charge versus
fitted, attenuation-corrected 3D deposit; their numerical ratio is not a
charge-conservation test.

Each standard wide event view is now accompanied by an LFGD-focused virtual-
voxel classification view (`eventN_zx_voxel_match.png`, and likewise for ZY
and XY).  It shows all exact MC voxels, matched reconstructed voxels, missing
MC voxels, and unmatched/off-track reconstructed voxels with separate marker
styles.  These images are linked directly from the generated HTML event scan.

Additional `eventN_*_track_match.png` views classify reconstructed 3D voxels
according to whether a fitted track node occupies the same voxel.  The CSV
also records the number of 3D voxels with and without track nodes and any
track-node-only voxels.

`voxel_residuals.csv` contains a row for every MC/reconstruction matched cube.
Multiple exact MC segments in one cube are combined into an energy-weighted
true crossing position (track-length weighting is the zero-energy fallback).
It stores the true and reconstructed positions, dx/dy/dz/dr, summed MC energy
deposit, fitted 3D deposit, and a sample-normalized charge-response residual.
Aggregate position and charge plots are written as
`voxel_position_residuals.png` and `voxel_charge_response.png`.

For controlled no-attenuation validation, `run_pipeline.sh` accepts
`DISABLE_HOMO_ATTENUATION=1`.  This applies coordinated parameter overrides to
both detResponseSim and hfgRecon.  Normal runs retain attenuation.

These changes belong to the local diagnostic study and can be submitted as
one set. They do not include the nd280Geant4Sim container/package change.


1. MC truth added to the flat tree
----------------------------------
File changed:
  LFGD_Recon_Simu/src/LFGDFlatTree.cxx

The tree maker already wrote:

  fiber_hits  - calibrated fibre hits
  hits3d      - reconstructed hfg_3d hits
  track_nodes - reconstructed THFGRecon/hfgRecon track nodes

A fourth TTree was added:

  mc_track_points - Geant4 truth trajectory points

Truth is read from:

  truth/G4Trajectories

Branches in mc_track_points:

  event      event number
  track_id   Geant4 trajectory identifier
  parent_id  parent trajectory identifier
  pdg        PDG particle code
  particle   particle name
  point      ordered point number within the trajectory
  x, y, z    true position
  t          true time
  px, py, pz true momentum at the trajectory point

All saved primary and secondary Geant4 trajectories are written. track_id and
parent_id preserve their relationships.


2. MC truth added to overlay plots
-----------------------------------
File changed:
  LFGD_Recon_Simu/plot_overlay.py

The plots now show:

  green line      MC truth trajectories
  grey points     calibrated fibre hits
  blue points     reconstructed 3D hits
  red line/points reconstructed track nodes

Each MC trajectory is drawn independently, preventing unrelated particles
from being connected. The modified Python file passed a syntax check.


3. CMake link failure fixed
----------------------------
File changed:
  LFGD_Recon_Simu/cmake/CMakeLists.txt

Observed failure:

  /usr/bin/ld: cannot find -lLFGDReconSimu

This package defines an executable but no LFGDReconSimu library. The following
declaration was added before ND280_EXECUTABLE:

  ND280_NO_LIBRARY(LFGDReconSimu)

Its position is important because ND280_EXECUTABLE decides at creation time
whether it should link the package's own library.


4. Executable PATH corrected
-----------------------------
Files changed:
  LFGD_Recon_Simu/build_flat_treemaker.sh
  LFGD_Recon_Simu/README.md

The old instructions incorrectly used:

  $PWD/build/$ND280_SYSTEM/bin

The executable is actually installed under:

  $PWD/$ND280_SYSTEM/bin

After building, use:

  export PATH="$PWD/$ND280_SYSTEM/bin:$PATH"

Expected location for the present system:

  LFGD_Recon_Simu/Linux-AlmaLinux_9.5-gcc_11-x86_64/bin/LFGDFLATTREE.exe


5. Documentation updated
-------------------------
File changed:
  LFGD_Recon_Simu/README.md

The documented flat.root content now includes mc_track_points and its truth
fields. The executable PATH examples were corrected as described above.


Local build and run
-------------------

Inside the configured Singularity environment:

  cd /home/tlux/HK/ND280++/LFGD_Recon_Simu
  . build_flat_treemaker.sh
  export PATH="$PWD/$ND280_SYSTEM/bin:$PATH"

The missing PYTORCHConfig.cmake message from hfgRecon configuration is a
warning and was unrelated to the LFGDFLATTREE link failure.

Run LFGD, for example:

  source ./switch-detector.sh lfgd
  ./run_pipeline.sh homo 10

Run HFGD, for example:

  source ./switch-detector.sh hfgd
  ./run_pipeline.sh hfg 10

Expected flat.root TTrees:

  fiber_hits
  hits3d
  track_nodes
  mc_track_points

Together these provide MC truth, fibre response, reconstructed 3D hits, and
reconstructed tracks for event-by-event diagnosis.


6. Detector response restricted to the detector under study
-------------------------------------------------------------
Files changed:
  LFGD_Recon_Simu/run_pipeline.sh
  LFGD_Recon_Simu/README.md

The broad upgrade-nd280plus and upgrade-nd280plus-homo options also enable
response simulation for several unrelated subsystems. This caused messages
from the SciFi/SWD TPlusPlusScint response path:

  Hit segment is not in a fiber

TPlusPlusScint is not the Homo response implementation; Homo uses THomoScint.
The pipeline now passes the supported disable options for SFG, SWD, SciFi,
ECAL, SMRD, TOF, TPC, and HAT. The selected Homo or HFG response remains
enabled. MC truth is carried through the event and is not removed by these
response-stage selections.


7. Empty plot layers made safe
-------------------------------
File changed:
  LFGD_Recon_Simu/plot_overlay.py

Constructing TGraph(0) is ambiguous in PyROOT 6.32 and selected a null
histogram-pointer overload. Drawing that invalid graph caused a segmentation
fault whenever an event had an empty fibre, 3D-hit, or track layer.

Graphs are now created with the unambiguous empty TGraph() constructor and
filled point by point. Empty graphs are not drawn or added to the legend.
Single-point MC trajectories use point drawing rather than an invalid line.
The deprecated TFile.tree attribute syntax was replaced with TFile["tree"].
The plotter also validates all required trees, prints their total entry
counts, and uses a unique canvas name for each projection.

The plotter additionally supports an inclusive --event-range and filters for
events containing fibre hits, 3D hits, or tracks. Multi-event output includes
an index.html page with counts and ZX, ZY, and XY thumbnails for each event.


8. Fixed and isotropic particle directions
-------------------------------------------
Files changed:
  LFGD_Recon_Simu/config.sh
  LFGD_Recon_Simu/generate_gps_macro.py
  LFGD_Recon_Simu/run_pipeline.sh
  LFGD_Recon_Simu/README.md

DIRECTION_MODE selects either "fixed" (the existing behavior using DIRECTION)
or "isotropic" (uniform over the full 4-pi solid angle using GPS angular type
iso with theta from 0 to 180 degrees). The source position remains controlled
by POSITION_MM and defaults to 0 0 1800 mm.

The existing SEED default of 12345 is passed to Geant4 and detector response,
making repeated isotropic runs reproducible. Default output names include the
direction mode (mu700_fixed or mu700_isotropic) to avoid overwriting results
from the other mode.


9. Fibre-hit charge distributions
----------------------------------
File changed:
  LFGD_Recon_Simu/plot_overlay.py

The plotter now creates eventN_fiber_charge.png for every selected event and
fiber_charge_selected.png containing the aggregate charge distribution for
all selected events. The histogram range is derived from the selected fibre
charges, and the HTML event index includes the per-event charge plot.

The charge-axis range is fixed to 0--100 with 25 bins (width 4). The plotter
also creates eventN_fiber_positions.png with six panels: two transverse
position-coordinate distributions for each of the three fibre projections.
The transverse coordinates are detected from the coordinates that vary over
the fibre grid, avoiding the arbitrary along-fibre reference coordinate.
An additional fiber_positions_selected.png contains the same six fibre
position distributions accumulated over all selected events and is linked
from the top of the HTML index.


10. Local hfgRecon fibre-boundary peak check
--------------------------------------------
File changed:
  SoftProj/hfgrecon/src/THomoPeakHit.cxx

THomoPeakHit::CheckIsPeak previously constructed both neighbouring fibre IDs
before checking whether they were inside the detector. Isotropic tracks made
boundary hits common, causing repeated Homo::Fiber errors with u=-1 or v=-1.

The local hfgRecon implementation now rejects negative neighbouring U/V
indices before constructing their geometry IDs. Neighbours are subsequently
looked up directly by geometry-ID value in the measured-hit map, avoiding an
ambiguous mapping between THomoGeom fibre-count names and GeomId U/V axes. A
boundary hit for which a two-sided peak test is undefined returns false for
that test. The installed hfgRecon copy in the container was not changed.

The temporary test that assigned the sum of the central and neighbouring
fibre charges to every Homo peak was reverted. It caused substantial charge
double counting: the summed reconstructed 3D-hit charge was typically about
2.4--2.7 times the sum of all original fibre charges. THomoPeakHit therefore
again carries only the central peak-fibre charge, matching the original
implementation. HFGD reconstruction was not affected by this test.

Overlay titles now also state the fibre-hit, 3D-hit, and track-node counts for
the displayed event, making empty reconstruction layers explicit.


11. Filter unmeasured along-fibre coordinates from overlays
------------------------------------------------------------
File changed:
  LFGD_Recon_Simu/plot_overlay.py

A fibre hit determines only the two coordinates transverse to the fibre. Its
GetPosition() nevertheless contains a geometry reference value for the third,
along-fibre coordinate. Plotting that value as a measured position produced
artificial straight grey bands (for example x=0 in one LFGD projection and
x=-1095 mm in one HFGD projection).

The ZX, ZY, and XY overlays now include a fibre only in the view whose two
axes are both transverse to that fibre. Overlay titles report the number of
fibre hits shown in that view as well as the total number in the event.


12. Correct global source centre and detector-boundary overlays
---------------------------------------------------------------
Files changed:
  LFGD_Recon_Simu/config.sh
  LFGD_Recon_Simu/run_pipeline.sh
  LFGD_Recon_Simu/plot_overlay.py
  LFGD_Recon_Simu/README.md

The previous default `(0,0,1800) mm` confused the detector centre in the
local PlusPlusTracker frame with the global coordinate system used by GPS.
The geometry places both LFGD and HFGD at global centre `(0,30,910) mm`, so
that is now the default source position.

Event overlays now include dashed projected boundaries for SWD, the enabled
LFGD or HFGD envelope, and SciFi. The target type is identified from the
transverse span of the fibre grid, independently of the input filename.

The source interface now distinguishes `POSITION_FRAME=plusplus` from
`POSITION_FRAME=global`. The default `(0,0,1800) mm` is interpreted in the
PlusPlusTracker layout frame and converted to `(0,30,910) mm` for the global
GPS command. A global position can still be supplied without conversion.

The detector-response invocation is assembled as one Bash argument array,
preventing any continued `-O` option from being parsed as a separate command.


13. Reconstruction-only rerun helper
------------------------------------
File added:
  LFGD_Recon_Simu/rerun_reconstruction.sh

This helper takes an existing detresponse.root and a distinct output folder,
then reruns HFGRECON, creates the flat diagnostic tree, and produces event
displays without repeating Geant4 or detector response. Extra command-line
arguments are forwarded to HFGRECON.exe for parameter studies.


14. Reconstructed 3D-hit position and charge plots
---------------------------------------------------
File changed:
  LFGD_Recon_Simu/plot_overlay.py

The event display now creates per-event reconstructed 3D-hit charge and
separate x/y/z position distributions. Aggregate versions contain all events
retained by the requested event selection. All new diagnostics are linked in
the generated HTML index.

An event-level charge-conservation diagnostic now sums all original fibre-hit
charges and all reconstructed 3D-hit charges. The aggregate scatter plot shows
the two sums against the equality line, while charge_sums.csv and the HTML
index report the values and 3D/fibre ratio for each event. Events above the
equality line are counted and reported by the script.


15. Open LFGD charge-model and parameter-scan studies
-----------------------------------------------------
No production reconstruction change is retained for cluster-charge summing.
The following items still need a controlled study:

* Determine how neighbouring LFGD fibre charge should contribute to the
  reconstructed energy and dE/dx without assigning the same fibre charge to
  overlapping peaks or multiple 3D candidates.
* Compare central-peak charge, unique-cluster charge, and a charge-sharing
  model that maintains an explicit mapping from original fibres to 3D hits.
* Scan hfgRecon.Hits3D.Min2DHitCharge.homo (including the original value 70
  and the current local value 20), since it currently acts on the central
  peak-fibre charge.
* Scan or compare FindPeaks.homo and UsePeakWeightedPosition.homo, together
  with matching efficiency, fake-hit rate, charge conservation, and dE/dx.
* Use charge_sums.csv and charge_sums_selected.png for every scan, but compare
  accepted/used fibre charge separately from the sum over all raw fibres when
  evaluating strict charge conservation.


16. Optional local-direction LFGD 2D clustering
------------------------------------------------
Files changed in SoftProj/hfgrecon:
  inc/THomoPeakHit.hxx
  inc/THomoHits3D.hxx
  src/THomoPeakHit.cxx
  src/THomoHits3D.cxx
  parameters/hfgRecon.parameters.dat

The original Homo peak finder remains available and is still the default.
The new hfgRecon.Hits3D.PeakClusteringMode.homo parameter selects the method:

  0: original peak finder and central-peak-fibre charge
  1: experimental local-direction clustering with unique fibre ownership

Mode 1 deliberately uses the original peak candidates so comparisons isolate
the clustering change. Around each peak, a straight local tangent is obtained
with a small PCA of that peak and up to four nearby peaks. Fibres in a narrow
longitudinal slice and a wider direction perpendicular to the tangent are
assigned to the best compatible peak. Each input fibre is considered once and
can therefore be owned by no more than one 2D cluster. The cluster position is
charge weighted and its stored charge is the sum of its uniquely owned fibres.

Each projection logs its input and assigned fibre counts and charge. This
provides a direct 2D conservation diagnostic before the still-independent
problem of selecting unique XZ/YZ/XY combinations in the 3D matcher is fixed.
The local hfgRecon library and executables compile successfully inside the
ND280++ Singularity environment.

Two small parameter override files were added to LFGD_Recon_Simu for direct
comparison. hfgrecon_original_2d.parameters.dat selects mode 0, while
hfgrecon_local_2d.parameters.dat selects mode 1 and enables use of its
charge-weighted 2D position in the 3D matching stage.
# Virtual truth overlay and reconstructed-coordinate validation (2026-08-05)

`LFGDFlatTree.cxx` now writes `mc_virtual_segments`, containing the exact
entry/exit points, energy deposit, track length, contributor information and
derived HOMO virtual-cube index for every `truth/g4Hits/homoVirtualCube`
segment.  `plot_overlay.py` overlays these segments in cyan and uses their cube
indices for the reconstruction-efficiency comparison instead of interpolating
the MC trajectory.

The overlay was used to isolate a HOMO projection-convention mismatch and an
invalid geometry-node transform in the hfgRecon 2D peak positions.  After the
oaEvent/hfgRecon fixes, events 16, 29 and 32 follow the MC trajectory in all
three views.  The exact cube-matching efficiency remains deliberately stricter
than visual alignment because reconstructed positions are charge weighted and
are not forced to voxel centres.
