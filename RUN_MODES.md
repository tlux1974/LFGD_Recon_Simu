# Running LFGD response modes

Run these commands inside the ND280++ Singularity container from
`/home/tlux/HK/ND280++/LFGD_Recon_Simu`.

## Detailed response (default)

```bash
./run_student_sample.sh lfg-best 100 detailed_test
```

This creates individual WLS photons, simulates MPPC pixels and runs the
CITIROC electronics. Fibre attenuation is enabled. Transient photon and
avalanche objects are explicitly released at the end of every event.

## Fast aggregate response

```bash
FAST_HOMO_RESPONSE=1 \
  ./run_student_sample.sh lfg-best 100 fast_test
```

This accumulates trapped and surviving photons per fibre without creating one
object per photon. It applies the MPPC PDE and finite-pixel saturation
statistically. It does not simulate detailed photon timing, pixel recovery,
crosstalk, afterpulsing or CITIROC electronics. Fibre attenuation remains
enabled. Use this mode for large LFGD samples and light-pattern studies.

## Disable fibre attenuation

This is mainly a validation mode:

```bash
DISABLE_HOMO_ATTENUATION=1 \
  ./run_student_sample.sh lfg-best 100 no_attenuation_test
```

It can also be combined with aggregate response:

```bash
FAST_HOMO_RESPONSE=1 DISABLE_HOMO_ATTENUATION=1 \
  ./run_student_sample.sh lfg-best 100 fast_no_attenuation_test
```

By default, `run_student_sample.sh` first generates a CSV containing the
primary particle of every event and then replays those primaries explicitly.
To let Geant4 GPS generate all 100 events directly in a single run, use:

```bash
FAST_HOMO_RESPONSE=1 DISABLE_HOMO_ATTENUATION=1 REPLAY_PRIMARY_EVENTS=0 \
  ./run_student_sample.sh lfg-best 100 fast_no_attenuation_gps_test
```

The default direction mode is isotropic. For a fixed GPS direction, add for
example `DIRECTION_MODE=fixed DIRECTION="0 0 1"` before the command.

For muons starting at the detector centre and travelling within a 5 degree
cone around the ND280 +X axis:

```sh
DIRECTION_MODE=cone DIRECTION="1 0 0" CONE_HALF_ANGLE_DEG=5 \
  ./run_student_sample.sh lfg-best 1000 xaxis_cone5
```

Change `CONE_HALF_ANGLE_DEG=10` for a 10 degree half-opening angle. Directions
are uniform in solid angle inside the cone and are saved event by event in
`primary_events.csv`. The normal student default vertex `(0,0,1800) mm` in the
PlusPlus frame corresponds to the HOMO centre `(0,30,910) mm` globally.

Each run writes to a new directory under `output/`. The main analysis file is
`flat.root`, containing `fiber_hits`, `homo_raw` and `homo_truth`.

## Memory diagnostic without raw/truth collections

To test whether the diagnostic HOMO hit collections cause the remaining
memory growth, run:

```bash
FAST_HOMO_RESPONSE=1 \
DISABLE_HOMO_ATTENUATION=1 \
REPLAY_PRIMARY_EVENTS=0 \
DETRESPONSE_PARAMETER_FILE="$PWD/detresponse_no_raw_truth.parameters.dat" \
  ./run_student_sample.sh lfg-best 1000 fast_no_attenuation_gps_no_truth
```

In this mode, `homo_raw` and `homo_truth` in `flat.root` are empty.
