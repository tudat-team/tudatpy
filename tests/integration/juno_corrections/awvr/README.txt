Juno PJ3/PJ6: measured wet troposphere integration

This directory contains a snapshot of the two fixed-purpose Juno fitting
scripts. It uses the installed tudatpy package, including the native AMC
reader and measured tropospheric correction. There is no machine-specific
import bootstrap. The ordinary main function performs the original X-only
and combined X/Ka estimations; these are full fits, not unit tests.

Local setup and run

1. Activate an environment containing the matching tudatpy build, NumPy and
   Matplotlib. Alternatively, from the checkout root, select this build with:

   export PYTHONPATH="$PWD/build/src:$PWD/src${PYTHONPATH:+:$PYTHONPATH}"

2. Place the existing PJ3/PJ6 data cache beside these scripts as data/, or
   create a local symlink. Do not add the data or generated output to git:

   cd tests/integration/juno_corrections/awvr
   ln -s /path/to/existing/juno-data-cache data

   Keep the original cache layout, including perijove-times.txt, amc/, odf/,
   kernels/, tro/, ion/, sff/ and other supplied ancillary directories. The
   standard Tudat SPICE/resource cache must also be configured as for the
   baseline Juno example. No mission data is downloaded by these scripts.

3. Run locally:

   python juno_odf_doppler_estimation_spin.py

   The default AWVR mode writes two residual PDFs under output/awvr/.
   Run the data-independent screening/masking tests with the same package:

   python -m pytest -q test_awvr_preprocessing.py

Correction and comparison settings

The fixed settings near the top of the main script are use_awvr=True,
awvr_unit='01', and use_measured_dry_delay=False. Unit 01 and unit 02 are
redundant radiometers at DSS-25; they are not different tracking stations.
Only the selected unit is used, without smoothing or averaging the units.

The native reader returns UTC seconds since J2000 and zenith wet/dry path
lengths in metres. The example conservatively requires all six quality
flags to be zero, finite/nonnegative wet and dry paths, and no 999 cm
missing retrieval (9.99 m after conversion). It reports retained counts,
algorithm identifiers and overlapping rejection-reason counts. Algorithm
010 is preserved as metadata rather than interpreted as a quality flag;
its undocumented meaning is not inferred here. Unknown nonzero quality
flags are rejected by this conservative selection.

At DSS-25, measured wet zenith paths replace TSAC wet paths. TSAC dry paths
remain active by default. The native measured correction uses Niell
mapping for Earth and the duplicate EarthKa station body. Linear
interpolation is limited to 30.01 s between adjacent accepted samples.
Gaps and times outside coverage use the ordinary TSAC model. Other
stations retain their TSAC corrections. Optional measured dry replacement
is selected only by setting use_measured_dry_delay=True.

Before assigning the HGA reference point, the example removes complete
60 s Doppler counts that cross a measured/TSAC source boundary. Both Earth
transmission and reception events are checked using iterated reference
SPICE geometry, with a one-second guard around count endpoints. This
source-transition mask runs in both AWVR and TSAC comparison modes.
Intervals wholly inside a gap remain available with TSAC fallback.

Set use_awvr=False for the TSAC comparison with that same boundary mask;
its PDFs are written under output/. Retain the printed sample/count
summaries when comparing runs. Propagation, parameters, priors, weights,
iteration limits and the original residual-plot style are preserved.

Completed local validation, 2026-10-18

Using the supplied PJ03/PJ06 data and the original 45 parameters, fixed
priors and six iterations, exact link/band/epoch support agrees in both
TSAC/AWVR comparisons. The source-transition mask removes 11 X and 16 Ka
counts. X-only RMS is 1.807846 -> 1.721530 mHz (1328 counts). In the
combined fit (2210 counts), Ka RMS is 3.559845 -> 2.144437 mHz on PJ03
and 3.180732 -> 2.509926 mHz on PJ06. X also improves in both passes.
The largest parameter change is 0.244 baseline formal standard errors
for X-only and 0.421 for the combined fit. These are conditional fit
comparisons, not an external determination of gravity-field accuracy.
The selected unit-01 wet correction is retained on both passes.

validation_summary.json records the counts and per-pass/band comparisons.
Native model/reader tests, the Python reader tests and all 36 preprocessing
tests pass. Full local logs, exact support, vectors, residuals and PDFs
are retained in juno-corrections/validation/awvr outside this source tree.

References: Buccino et al. (2021), doi:10.1029/2021RS007387, describes
the Juno AWVR/Niell calibration; Zannoni's thesis, section 4.4.2, motivates
the same-pass TSAC comparison and checking rather than assuming a lag:
https://amsdottorato.unibo.it/id/eprint/6536/
AMC field/flag specification:
https://atmos.nmsu.edu/PDS/data/jnogrv_1001/DOCUMENT/DORS_002_020501.TXT
