Optional local-direction LFGD 2D clustering
============================================

Overview
--------

An experimental two-dimensional fibre-clustering method has been added to the
local hfgRecon package. The original method has not been overwritten and
remains the default, allowing both methods to be run on the same detector-
response input and compared directly.

The method is selected with:

  hfgRecon.Hits3D.PeakClusteringMode.homo

The available values are:

  0: Original peak finder. The 2D-hit charge is the central peak-fibre charge.

  1: Experimental local-direction clustering. Fibres are assigned uniquely
     to a peak and the 2D-hit charge is the sum of the assigned fibre charges.


Local clustering algorithm
--------------------------

The new method deliberately starts with the same peak candidates as the
original method. This makes it possible to isolate the effects of fibre
assignment and cluster charge when comparing the two modes.

For every peak, a local straight-line direction is estimated using a small
principal-component fit (PCA) containing the peak and up to four nearby peaks.
This is computationally inexpensive and allows the direction to change along
the track instead of requiring one global straight line or helix.

Input fibres are considered in a narrow slice along the local track direction
and a wider region perpendicular to it. If a fibre is compatible with more
than one peak, it is assigned to the best local match. Peak charge is used only
as a deterministic final tie-break.

Every input fibre can be owned by at most one 2D cluster. For each cluster:

  cluster position = charge-weighted position of the owned fibres
  cluster charge   = sum of the charges of the owned fibres

Fibres that are not compatible with any peak remain unassigned. They are
reported explicitly rather than silently added to more than one cluster.

This first implementation uses a local direction independently at each peak
instead of walking from a detector-boundary endpoint. This avoids dependence
on selecting the correct endpoint and avoids accumulating an incorrect
direction after one bad step. An explicit endpoint walker can still be added
later as another clustering mode.


Charge and ownership diagnostics
--------------------------------

For each XZ, YZ and XY projection, mode 1 reports:

  number of input fibres
  number of peak candidates
  number of clusters
  number of unassigned fibres
  sum of all input fibre charges
  sum of all assigned fibre charges

Because ownership is unique, the following relation holds separately in each
projection:

  assigned cluster charge + unassigned fibre charge = input fibre charge

This is the relevant charge-conservation test for the 2D clustering stage.
The subsequent 3D matching is a separate problem and is not changed by this
implementation.


Boundary correction
-------------------

The Homo peak neighbour check was also corrected. Negative neighbouring U/V
indices are rejected before a geometry ID is constructed, preventing the
repeated Homo::Fiber errors for u=-1 or v=-1.

Neighbouring measured fibres are looked up directly by their geometry-ID
value. This avoids relying on the ambiguous relationship between THomoGeom
fibre-count names and the GeomId U/V convention.


One-event validation
--------------------

The local hfgRecon library and executables compile successfully in the ND280++
Singularity environment. Event 0 of homo_mu700_isotropic completed with both
clustering modes.

The mode-1 two-dimensional diagnostics were:

  XZ: fibres=277, peaks=136, clusters=136, unassigned=13
      input charge=5855.57, assigned charge=5682.13

  YZ: fibres=186, peaks=91, clusters=91, unassigned=13
      input charge=4270.57, assigned charge=4072.59

  XY: fibres=316, peaks=94, clusters=94, unassigned=33
      input charge=6045.12, assigned charge=5371.21

With the current, unchanged 3D matcher, this event produced:

  original 2D mode: 391 reconstructed 3D candidates
  local 2D mode:    243 reconstructed 3D candidates

This difference must be evaluated with the event displays and efficiency
diagnostics. It is not yet evidence that either result is preferable.


Comparison parameter files
--------------------------

Two self-contained runtime-parameter override files are provided:

  hfgrecon_original_2d.parameters.dat
  hfgrecon_local_2d.parameters.dat

Both explicitly set FindPeaks and Min2DHitCharge=20 so that the comparison does
not depend on whether runtime parameters are loaded from the local package or
the installed container package. The local-mode file also enables
UsePeakWeightedPosition.


Running the comparison
----------------------

Enter the normal ND280++ environment and select the local hfgRecon build:

  source ~/HK/ND280++/switch-hfgrecon.sh local

From LFGD_Recon_Simu, run the original method into its own output directory:

  ./rerun_reconstruction.sh \
      output/homo_mu700_isotropic/detresponse.root \
      output/homo_mu700_isotropic/reco_original_2d \
      -O par_override="$PWD/hfgrecon_original_2d.parameters.dat"

Run the experimental local method into a different directory:

  ./rerun_reconstruction.sh \
      output/homo_mu700_isotropic/detresponse.root \
      output/homo_mu700_isotropic/reco_local_2d \
      -O par_override="$PWD/hfgrecon_local_2d.parameters.dat"

The two output directories contain independent reco.root and flat.root files,
logs, charge diagnostics and event displays.


Important environment check
---------------------------

Both the executable and shared library must come from the local build. Merely
checking command -v HFGRECON.exe is insufficient if LD_LIBRARY_PATH still puts
the installed container library first. The switch-hfgrecon.sh helper updates
both PATH and LD_LIBRARY_PATH.

The selected paths can be checked with:

  source ~/HK/ND280++/switch-hfgrecon.sh status
  command -v HFGRECON.exe
  ldd "$(command -v HFGRECON.exe)" | grep hfgRecon


Files changed in the local hfgRecon repository
----------------------------------------------

  inc/THomoHits3D.hxx
  inc/THomoPeakHit.hxx
  parameters/hfgRecon.parameters.dat
  src/THomoHits3D.cxx
  src/THomoPeakHit.cxx


Remaining 3D work
-----------------

The existing 3D matcher can still reuse the same XZ, YZ or XY 2D object in
several valid combinations. Consequently, unique 2D ownership does not by
itself guarantee a unique or charge-conserving 3D assignment.

The next stage should score competing XZ/YZ/XY combinations and either choose
one global best assignment or explicitly divide a 2D cluster charge among
competing 3D candidates.
