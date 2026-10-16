Fixed-ephemeris lunar laser ranging diagnostic
==============================================

Current primary: paired March--April 2017 arc
--------------------------------------------

The current primary result uses fixed DE440/ICRF and the independently
preselected INPOP17a G+IR APOLLO/Grasse station pair. All **306 points** remain:
69 APOLLO points on nine UTC nights and 237 Grasse points on 24 nights, across
all five reflectors. Raw one-way RMS about zero is **58.488 mm** for APOLLO and
**359.583 mm** for Grasse. A retained **+5.43 m** A15 observation dominates the
Grasse RMS; its origin remains unexplained. These discrepancies do not measure
intrinsic observation noise. No offset, trend, uncertainty inflation or
residual-based cut is applied.

The `raw package <../analysis/laser_ranging/results/paired_long_arc_20170301_20170430/README.md>`_
retains every CSV row, source audit and provenance. The
`time plots <../analysis/laser_ranging/results/paired_long_arc_time_plots/README.rst>`_
retain the full range plus labeled central detail, wavelength markers and raw
nightly means. Two months do not establish a robust monthly period.

Prepare the hash-pinned inputs in
`llr_long_arc_dataset.json <../analysis/laser_ranging/llr_long_arc_dataset.json>`_
and the declared kernels, then run from the repository root::

    python analysis/laser_ranging/llr_long_arc.py \
      --data-root /path/to/longarc-20170304 \
      --kernel-dir /path/to/llr/kernels \
      --output /path/to/new-output

Native source ``a366f436c`` passed AtmosphereCorrection **15 cases / 9,880
assertions** and the Python settings check at actual 532.1/1064.2 nm wavelengths,
without rounding observed wavelengths or changing formulas/tolerances. Kernel
SHA-256 is
``cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef``.
The original driver runtime was **26.6935 s**. A subsequent controlled comparison
uses geocentric TT only for slow linear station drift: solver time falls from
**19.631 to 12.778 s**, and driver time from **27.066 to 18.994 s**. All 306 rows
(11,934 CSV fields), geometry and source metadata remain byte-identical. The
coordinate change from omitting the topocentric term in that drift calculation
is below 2e-15 m; the observation clock and simultaneity calculations are unchanged.
The `runtime control <../analysis/laser_ranging/results/drift_clock_runtime_control/README.rst>`_
records commands, source and binary provenance. Other optimization proposals
remain unimplemented.

Additional bounded checks find no missing or duplicated millisecond TDB-minus-TT
launch conversion: three absolute epochs agree with independent ERFA/SOFA within
2.5 ps. Native lunar displacement agrees with independently evaluated original
constants within 0.145 nm, including explicit radius and once-only application
checks. The daily Grasse archive repeats the same retained +5.43 m record.
Evidence is in the `absolute clock audit
<../analysis/laser_ranging/results/absolute_launch_clock_audit/README.rst>`_ and
the `deformation and source audit
<../analysis/laser_ranging/results/tide_parameter_plumbing/README.md>`_. These
controls narrow the possible causes; they do not explain the remaining residuals.

The `finite moving-body Shapiro control
<../analysis/laser_ranging/results/moving_monopole_control/README.md>`_
changes individual leg lengths by at most 2.282 micrometres at three saved
observations. Its nearly cancelling fixed-endpoint sums are not newly iterated
observables. The `INPOP17a time-kernel control
<../analysis/laser_ranging/results/inpop17a_time_and_prism_audit/README.rst>`_
agrees with SOFA within 1.393 ns at three launches. Neither check identifies
a centimetre error. The latter report also records the unresolved optical
reference point of INPOP's fitted reflector coordinates.

INPOP17a comparison
------------------

The optional INPOP17a control uses its own translation, principal-axis rotation,
GMs, published G+IR station/reflector coordinates, degree-2 surface tides and
lunar spatial scale. The kernel pool is cleared before loading the matched
products. Its original header constants support association with the published
G+IR fit, although the archive does not explicitly name that solution column.
The documented Lunokhod table misprints are reconstructed from the cited INPOP13c
baselines before any residual calculation.

The fitted reflector optical reference point remains unspecified in the retrieved
sources. Both optical conventions were declared before calculation; neither is
chosen by residual size. All original points remain, with no estimation or bias
correction. Raw one-way means and RMS values are:

.. csv-table:: Full March--April 2017 comparison, millimetres
   :header: "Separate glass term", "APOLLO mean / RMS (69)", "Grasse mean / RMS (237)"

   "JPL n*h per leg", "-47.439 / 49.408", "-37.637 / 358.749"
   "None", "-3.363 / 12.341", "+10.348 / 356.316"

The retained Grasse discrepancy remains about 5.47 m in the second case.
Using every point, its median absolute residual is 12.529 mm and its 95th
percentile is 37.777 mm; APOLLO's corresponding values are 8.212 and 22.850 mm.
These order statistics do not remove or recenter observations. The March 22
APOLLO spread between five reflector nightly means falls from 139.25 mm in
the DE440 case to 17.06 mm in the second INPOP case, but a shared trough and
recovery remain. Intrinsic measurement-noise residuals have not been achieved.

The published fit data end before this observation arc. Its station range biases,
including a Grasse term equivalent to 12.3 mm one-way, are not applied here.
Optical reference points, differing Earth loading/EOP reductions and small omitted
lunar terms prevent an exact reproduction of the fitted reduction. The `complete
comparison package <../analysis/laser_ranging/results/inpop17a_long_arc_20170301_20170430/README.md>`_
contains both raw time series, PNG/PDF figures, source references, an unchanged
eleven-point DE440 control and the independent actual-radius tide check.

Historical single-night setup
-----------------------------

``analysis/laser_ranging/llr_residuals.py`` compares all eleven APOLLO normal points
dated 2017-03-09 with a fixed DE440 model. It estimates no parameters and does not
remove a mean or trend, select points by residual, or inflate their uncertainties.
The input subset is checked by SHA-256 before use. Data and kernels remain outside
the repository. The adjacent ``analysis/laser_ranging/llr_dataset.json`` supplies
portable provider URLs, SHA-256 values, and the exact selection recipe: retain
all source lines beginning with the ASCII bytes ``512017 3 9`` in original order,
including their existing newlines, from the hash-verified full Group-e release.
The result must contain eleven rows (original lines 134–144) and match the pinned
subset hash. The model performs no downloads; supply ``--data-dir``,
``--kernel-dir`` and ``--output-dir`` explicitly after preparing these inputs.

This is a diagnostic, not a reproduction of JPL's fitted reduction. The APOLLO
station position comes from the 2005 WGS84 survey, whose stated accuracy is about
0.1 m; the fitted DE440 station position and velocity have not been recovered.
DE440's published centimetre residuals do not imply millimetre residuals with an
approximate station. Normal-point uncertainty is measurement precision, not an
allowance for modelling error.

Model conventions
-----------------

* The model uses ICRF axes throughout the environment.
  A hash-pinned analysis text frame, ``icrf_spice_alias.tf``, is an exact identity
  alias to SPICE's ``J2000``. Modern DE440 products labelled J2000 use ICRF axes,
  as documented in `NAIF's ICRF versus J2000 discussion
  <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html#ICRF%20vs%20J2000>`_.
  The native Earth rotation model interprets its J2000 option as the dynamical
  mean equator/equinox and applies a SOFA frame bias. The historical diagnostic
  combined that biased Earth rotation with ICRF translation and lunar rotation.
  Historical committed sources and the recorded source patch preserve that calculation.
  No native frame semantics or provider products are changed by the correction.

* The input UTC epoch is laser launch, constructed as a precise ``Time``. The
  selected Earth-orientation converter maps it to TDB. The native n-way model
  computes the atomic interval with ``time_scale_for_observable=tai_scale``;
  there is no second clock-transfer bias. No leap occurs during these intervals,
  so observed UTC interval and TAI interval agree. The observable is full
  ``c * RTT``; the reported one-way residual is
  ``(c * observed_RTT - computed_range) / 2``.
* Only the declared DE440 translation and ``MOON_PA_DE440`` rotation kernels are
  loaded. Reflector coordinates are Park et al. (2021), Table 1, in that PA
  frame. JPL's 2020 reduction evaluates lunar orientation directly in TDB.
* Earth has IERS2010 solid and pole displacement models. The Earth terminal
  spatial transformation is IERS Eq.11.19 with Sun potential. An analysis-local
  custom displacement also applies the simultaneity term
  ``-((v_E dot r)/c²) * rotational_station_velocity``. Deformation velocity
  over the approximately two-microsecond offset is omitted. The omitted
  acceleration terms are about 1.34 micrometres at Earth's surface.
* Lunar elastic displacement uses the existing degree-2/3 model, Earth and Sun,
  radius 1738 km, and JPL2020's ``h2=.0423, l2=.0107, h3=.0234, l3=.0030``.
  The full displaced PA vector is converted once with JPL2020 Eqs.16–17:
  ``delta_r = -(LM + U_sun/c²)r - .5(v_M dot r)v_M/c²``,
  ``LM=1.4825e-8``. Earth's TT scaling is not applied to the Moon. The newer TL
  convention is not substituted for this published DE440 reduction convention.
* Solar, Earth and lunar first-order Shapiro corrections are evaluated inside
  each leg's light-time solution. This differs from the DE440-fit JPL2020 model
  (Section 3.17), which omits the Moon's sub-millimetre term as well as Jupiter
  and other distant planets; retaining the lunar term is a physical improvement,
  not an exact replay of that fitted reduction. Mendes–Pavlis optical group delay with FCULa
  mapping uses the recorded wavelength and measured pressure, temperature and
  humidity. The existing Tetens utility supplies vapour pressure from humidity.
  Weather is interpolated between the measured UTC tags and held outside the
  first/last tag for the short transmit/receive interval. There is no standard
  atmosphere fallback.
* JPL2020's reflector coordinates refer to the array front face. The explicitly
  supplied per-leg internal path is ``n*h``, with ``n=1.4555`` and
  ``h=.0275 m`` (Apollo) or ``.0435 m`` (Lunokhod). The calibrated APOLLO point
  already represents a fitted array profile centre, so the script applies no
  second array-width or centroid correction. The native signature's zero-offset
  point supplies only this stated internal optical path.
* Light-time corrections are iterated. Failure to converge raises an exception.
  The native fractional tolerance is explicitly ``1e-13`` (the legacy Python
  binding calls its third argument ``absolute_tolerance``). No failed point is
  silently dropped or accepted at a looser tolerance.

Outputs and limitations
-----------------------

``--ephemeris de421-mer`` selects a separately pinned, coherent DE421 SPK,
DE421 binary lunar PCK and its ``moon_080317.tf`` frame definitions. Both frozen
dates lie inside their common 1900–2051 coverage. The lunar fixed frame is
explicitly ``MOON_ME_DE421`` and reflector coordinates are the DE421 MER column
of Park et al. (2021), Table 1; DE440 PA coordinates are not reused as MER
coordinates. Sun, Earth and Moon mass parameters come from the public DE421
header's GMS, GMB, EMRAT and AU, converted to SI, rather than the DE440 GM file.
``de421_dataset.json`` pins every product, the constants derivation and reference
point. Supply its cached files in the selected ``--kernel-dir``.

This is a fixed-station ephemeris/realization sensitivity. The published MER
coordinates agree with Eq.19-rotated DE440 PA coordinates within 0.758 mm from
printed precision; they are not an independently recovered original DE421 fit.
The actual DE421 versus DE440 orientation difference remains: at the frozen
events it changes the A15 inertial surface vector by about 0.319 m (APOLLO) and
0.337 m (Matera). Station coordinates, EOP, elastic tides, spatial conversions,
range corrections and iteration settings stay fixed across the two ephemerides.
The separate Hofmann station options retain DE440. No fitted transformation or
choice by residual improvement is introduced.

``normal_points.csv`` retains every row and its original uncertainty, together
with the one-way residual, normalized residual and sequential Shapiro,
atmospheric and glass increments. These increments come from re-solving the
same observation with successive corrections; they are not a fixed-event
additive decomposition. ``summary.json`` records raw mean, RMS about zero,
scatter about the mean, weighted mean, normalized RMS, unfitted chi-squared,
station coordinates and hashes of every loaded kernel.

``event_diagnostics.json`` records the native link-end epochs, Earth transmit
and receive elevations and registered per-leg correction values for each of the
four sequential solves. The event order is Earth transmit, Moon receive, Moon
retransmit, Earth receive. These native epoch dependent variables are doubles;
they cannot resolve all differences retained internally by ``Time``. The
separate launch-time split preserves the precise input epoch. Per-leg correction
values are in seconds, in registration order (Shapiro, atmosphere, glass as
enabled). They describe the corrections at each solution's events rather than
differences of separately solved ranges.

Unresolved differences include the station reference frame/tide convention and
velocity, ocean/atmospheric loading, lunar dissipative and thermal displacement,
and the exact calibration version used in the DE440 fit. Earth orientation uses
Tudat's selected public EOP configuration rather than a claimed copy of JPL's
historical solution. The full constant lunar tide is retained: JPL2020
Section 3.9 explicitly leaves Eq. 10's constant part unmodified, and Section 3.12
applies that same equation to the Moon, with constant components tabulated in
Appendix A2. These limitations must accompany the results rather than be absorbed
into a fitted range bias.

Published station-only sensitivity
----------------------------------

The independent Hofmann (2017) standard station solution is selected explicitly
with ``--station-profile hofmann2017-apollo-standard`` or, for the Matera case,
``--station-profile hofmann2017-matera-standard``. The separate
``hofmann2017_station_profiles.json`` pins the source tables, SRP/tide convention,
reference epoch and velocity interpretation. These options change station
coordinates and their published linear velocities only; they do not import
fitted biases, reflector coordinates or the source's other reduction models.

An independently justified axes correction was evaluated with the same frozen
observations, station parameters, tides, no loading, 30-second Time splines and
``1e-13`` fractional tolerance. All eleven APOLLO and both Matera points remain.
The explicit ICRF identity alias leaves the DE440 translations and lunar PA
rotation bitwise identical at every probed event; only the Earth frame bias is
removed. No residual mean or offset is fitted.

=========================  ==============  ==================
Station configuration      Mean (mm)       RMS about zero (mm)
=========================  ==============  ==================
APOLLO INPOP17a G+IR        -25.500         49.900
APOLLO survey2005           +114.175        121.085
Matera SLRF2020             -65.200         65.228
=========================  ==============  ==================

The computed half-range increases by 86.705–121.784 mm for APOLLO and by
241.035/189.799 mm for Matera. These full re-solves agree with the independently
projected station displacement at the original events within 21.2 micrometres
for APOLLO and 83.9 micrometres for Matera. The corrected INPOP17a case still has
42.893 mm scatter, with reflector means A11 -85.690, A14 +14.505, A15 -51.002
and L1 +25.689 mm. Correcting the frame mismatch does not establish a complete
centimetre reduction. The compact raw CSVs, per-point projected/full changes,
exact diagnostic source patch, alias and software/input hashes are retained in
``analysis/laser_ranging/results/icrf_frame_correction``. Historical results
remain separate and retain their original frame convention.

The same evidence directory contains an independent 80-digit coefficient audit.
Official DE440 ASCII coefficients identify fit dates 2020-06-25 and agree with
the pinned PCK librations within ``4.17e-13`` radians at all thirteen frozen
launches, below a micrometre at a reflector. The PCK's 2021 packaging date does
not indicate a different lunar fit. At these launches, direct SPICE's double
epoch conversion can move absolute Earth/Moon positions by up to 1.75 mm;
independent stored-polynomial evaluation at the already rounded epochs differs
by at most 27.1 micrometres. These are state controls, not measured range errors,
and do not explain the remaining centimetre reflector pattern.

The default ``--station-profile survey2005`` retains the approximate survey and
zero velocity. The explicit ``--station-profile inpop17a-gir`` uses the
independently preselected INPOP17a G+IR solution from Table 8, printed page 21:
``(-1463998.9155, -5166632.7690, 3435012.8659) m`` at J2000 and
``(-0.0130, +0.0001, 0.0) m/year``. The row expressions are algebraic offsets,
not unsigned coordinates. The caption specifies J2000 and ITRF2005 reference
values. The cited Manche reduction specifies TT for station drift; this
experiment explicitly interprets the epoch as JD 2451545.0 TT and the year as
365.25 days. The displacement is applied to both the intrinsic and transformed
Earth station. DE440 body states, lunar orientation and reflector coordinates
remain fixed. The report's APOLLO fit ends on 2016-11-25, before this selected
night; neither a station fit nor a choice between published solutions uses
these eleven residuals.

``llr_station_profiles.json`` pins the source URLs, PDF hashes, row algebra and
published component uncertainties. Those uncertainties are metadata, with no
change to normal-point weights. This is a station-only sensitivity, not a
matched INPOP or DE440 reduction. Fitted station coordinates are not guaranteed
to realize canonical ITRF2005 or the selected EOP frame, so no speculative frame
transformation is applied. A tide-free interpretation is supported by the
cited regularized station coordinates and IERS routines, but not explicitly
flagged in the 2017 report. The older reduction includes atmospheric and CSR4.0
ocean loading, whose precise 2017 choices are not restated; this experiment
adds neither. These limitations prevent treating the formal coordinate errors
as an accuracy claim for the mixed model.

With 30-second cubic tables, the published station changes the first launch's
kinematic position relative to the survey by
``(-0.241819751, +0.146447608, +0.299287507) m``. All eleven residuals remain
unfitted: mean ``+0.079877301 m``, RMS about zero ``0.089294128 m``, scatter
``0.039913133 m``, weighted mean ``+0.067742052 m``, normalized RMS ``36.989336``.
The residual changes range from ``-0.154586852`` to ``-0.122534931 m``.
A survey control with the same native kernel exactly reproduces the archived
30-second survey CSV, isolating the station change from intervening builds.

``analysis/laser_ranging/results/apollo_20170309`` retains exact raw CSV outputs
for direct SPICE and 30-second survey baselines, this sensitivity and its
survey control. Portable summaries retain every numeric field and provenance
hash, shorten runtime import paths and identify the original summary hashes.
The archive identifies Git commits/blobs for the source and input/profile
manifests, with minimal patches for the executed source before formatting.
Use ``git show REV:path`` (or ``git show BLOB``), apply the listed patch and
verify the original SHA-256; ``provenance/git_sources.json`` gives the records.
Large kernels and provider observations remain external.
``station_sensitivity.csv`` records every raw residual change and launch
kinematic station shift. The station experiment and its control used native
kernel SHA-256
``1864ed04530210fa056e78128bdbbfe1ad536cae6a9d3ade4f23295fbeda1744``.

Optional ocean loading and historical radial-basis controls
----------------------------------------------------------

``--got410c-harpos PATH`` enables the hash-pinned GOT4.10c HARPOS product from
``https://massloading.net/imsl/load_har_list/toc/got410c/toc_got410c_harmod.hps``.
The existing helper selects the exact case-sensitive ``Apollo`` site. Its
callback is included once in each Earth station's motion list, intrinsic and
final; lunar stations are unchanged. It evaluates the 17 supplied TT harmonics
and returns terrestrial XYZ displacement. This is a whole-Earth-centre-of-mass
STW105/REAR realization,
different from JPL FES2004 and the older INPOP CSR4.0 model. It includes neither
ocean-pole nor equilibrium terms. The current helper uses geodetic Up, following
the product narrative and the provider's source code; the analysis README records
the formal HARPOS specification's conflicting radial wording and provenance.
Loading remains disabled by default for both station profiles.

The following archived controls used the original radial-Up interpretation.
Their bytes and metadata are retained; they are not results from the corrected
geodetic helper. With all eleven points and 30-second cubic tables, survey plus
loading gives
mean ``+0.231725481 m``, RMS ``0.235149176 m`` and scatter ``0.039980448 m``.
INPOP17a station plus loading gives mean ``+0.092049003 m``, RMS
``0.100350376 m`` and scatter ``0.039964722 m``. The respective mean changes
from the recorded unloaded controls are ``+12.175132`` and ``+12.171702 mm``.
No station, reflector, uncertainty or correction parameter is fitted or chosen
using these residuals. Group means retain the night's sampling:

.. csv-table:: Mean one-way residual by reflector, millimetres
   :header: "Reflector (N)", "Survey", "Survey + loading", "INPOP station", "INPOP + loading"

   "A15 (4)", 196.049, 208.154, 55.557, 67.652
   "A11 (2)", 174.437, 186.501, 27.656, 39.698
   "A14 (3)", 251.387, 263.652, 115.450, 127.723
   "L1 (2)", 263.910, 276.203, 127.381, 139.684

The archive's ``loading_comparison.json`` records compact run summaries,
product metadata and provenance, including the different native build from
the unloaded controls. The input kernels and EOP resource hashes are unchanged.
The two additional CSVs preserve every raw result; ``loading_per_reflector.csv``
retains the unrounded group means. Launch loading changes from approximately
``(-1.648, -12.576, +7.533) mm`` to ``(-1.275, -11.366, +8.769) mm`` over the
night, identically for both profiles under the provider's fixed site basis.

SPICE epoch precision experiment
--------------------------------

The default ``--spice-knot-step 30`` uses the Sun, Earth and Moon translation
ephemerides from the existing Time-valued ``interpolated_spice`` route, using cubic splines and
integer TDB knots. Integer knot epochs are exact doubles, so sampling SPICE does
not lose their time tags. The table extends at least approximately one hour
outside the observations to suppress natural-spline boundary effects. Cubic
splines avoid the generic Lagrange interpolator's near-node handling. The
analysis-local spatial-correction callbacks retain their original direct-SPICE
queries. This experiment changes no station position, constants, measured
weather, uncertainty, correction formula or convergence tolerance.
The explicit option ``--spice-knot-step 0`` reproduces the original direct-SPICE
baseline; ``60`` retains the coarser convergence control. The default is selected
from the independent time-continuity and table-convergence controls below.

``ephemeris_controls.json`` compares each selected ephemeris against direct
SPICE at 22 exactly representable off-knot epochs, and compares displacement
over two precise epochs 20 ns apart that convert to the same double. At the
2017 observation date, one double epoch spacing is about 119 ns. Direct SPICE
therefore returns identical states for the two precise test epochs, missing
0.600 mm of Earth motion and 0.615 mm of Moon motion. The cubic tables recover
that motion within 24 micrometres at 60-second knots and 25 micrometres at
30-second knots. Maximum off-knot position differences from direct SPICE are
0.032/0.092 mm for Earth/Moon at 60 seconds and 0.063/0.062 mm at 30 seconds.
The latter comparison is a control at representable epochs, not evidence that
direct SPICE resolves arbitrary precise epochs.

All eleven points were retained in the direct, 60-second and 30-second runs.
Recording native dependent variables reproduced every original direct-SPICE CSV
field exactly. The following statistics are unfitted one-way residuals in metres:

.. list-table:: Unfitted one-way residual statistics
   :header-rows: 1

   * - Ephemeris
     - Mean
     - RMS about zero
     - Scatter
   * - Direct SPICE
     - 0.219974697
     - 0.223582017
     - 0.040000639
   * - Cubic, 60-second knots
     - 0.219550593
     - 0.223164011
     - 0.039996411
   * - Cubic, 30-second knots
     - 0.219550349
     - 0.223163915
     - 0.039997214

Changing from direct SPICE to 60-second cubic tables changes individual complete
half-ranges by up to 1.710 mm. Halving the table spacing changes complete
half-ranges by at most 0.0262 mm (RMS 0.0126 mm). The maximum 60-to-30-second
changes in geometric half-range, sequential Shapiro, atmospheric and glass
increments are respectively 0.0283, 0.0278, 0.0251 and 0.0214 mm. These are
observed numerical sensitivities for this fixed case, not adjusted tolerances
or a general error bound.

The final A14 point illustrates why a sequential range increment should not be
used alone to validate a correction formula. Its glass increment is
0.041121006 m with direct SPICE and 0.039928138/0.039907932 m with the
60/30-second tables, whereas the registered native contribution is exactly
``n*h/c`` on each leg for all three routes, equivalent to a 0.040026250 m
half-range contribution at fixed events. The direct-SPICE anomaly largely
disappears after restoring smooth translation, but the remaining approximately
0.1 mm sensitivity of the sequential solve is unresolved. Tightening only the
final point's fractional convergence tolerance from ``1e-13`` to ``1e-14``
leaves its atmospheric and complete ranges bit-for-bit unchanged; this control
does not implicate the stopping tolerance. The double-valued
event output cannot diagnose every precise event shift. This experiment does
not claim closure of all sub-millimetre numerical paths. The change of default
does not use residual reduction as its criterion.

These runs used the actual imported native kernel SHA-256
``c475a3e90a39eed4332442f58c631398cd608f2195d25bc08cd59e3cfd0bd866``.
Each output summary records its own script, kernel, input and Earth-orientation
resource identities, so later builds can be distinguished from this experiment.
The separate stricter-convergence control used kernel SHA-256
``029f3ea61f723093697a8b055e304189d3906ffb541d6fae929008890625deff``;
its same-process ``1e-13`` result also exactly matches the original 30-second run.

Frozen Matera prediction case
-----------------------------

``--case matera20220821`` retains both Matera 7941/Apollo15 CRD normal points
from 2022-08-21, outside the DE440 fitting span. The same solver, corrections,
30-second cubic ephemerides and convergence tolerance are reused. The native
CRD reader supplies the original transmit-tagged UTC seconds of day, RTT,
measured 532 nm weather and calibration metadata. Station delay and receive
amplitude corrections are already applied according to the pinned headers;
the recorded ``108680.1 ps`` calibration is not subtracted again. No second
array-profile or centre-of-mass adjustment is introduced.

``matera_dataset.json`` pins the observations and station products. The native
SINEX reader selects SLRF2020 release 250205, 7941/A/solution 1, DOMES 12734S008.
The documented optical-axis intersection has zero XYZ eccentricity over the
whole day. Native station-history motion is attached to both intrinsic and
final Earth stations, with the parsed 2015-01-01 reference vector and velocity.
Its declared nominal-civil-UTC convention uses native reference/event time
conversion and MJD differences divided by 365.25 for linear terms; no station
position is manually propagated. The frozen history conservatively uses the
open-ended ``SOLUTION/EPOCHS`` observation span rather than claiming it is a
discontinuity record. The pinned official ITRF2020 SLR PSD and discontinuity
products contain no Matera entries. No post-seismic term is assumed merely
from the presence of a station vector. Fitted seasonal terms and loading are
not added in this initial case.

The two one-way residuals are ``+0.177724421`` and ``+0.122627497 m``; mean
``+0.150175959 m``, RMS about zero ``0.152681814 m`` and scatter
``0.027548462 m``. CRD bin RMS values ``318`` and ``363 ps`` describe return
spread, not normal-point standard error. The CSV leaves uncertainty and
normalized residual blank; weighted statistics are null. No division by the
square root of return count or noise-level claim is made.

``analysis/laser_ranging/results/matera_20220821`` retains the exact raw CSV,
portable summary and direct references to its executed source Git blobs.
The native reference-position and one-Julian-year history
controls agree with the parsed solution to double-position precision. The
shared-script APOLLO control reproduces all eleven archived CSV rows exactly.

Optical input and observation-release controls
---------------------------------------------

An independent evaluation of IERS chapter 9 expressions at the archived native
events agrees with all 48 per-leg optical delays from the survey APOLLO,
INPOP17a APOLLO and Matera cases within ``9e-16 m``. This checks the actual
geodetic latitude, ellipsoidal height, local elevation, 532 nm wavelength,
UTC-interpolated pressure, temperature and humidity supplied to the model;
it does not establish the accuracy of the measured atmosphere. Nominal heights
are 2786.620 m, 2786.669 m and 536.982 m respectively. A 100 m height change
changes these delays by at most 0.097 mm. The
`Matera site log <https://ilrs.gsfc.nasa.gov/network/stations/active/MATM_sitelog.html>`_
records a barometer 2.4 m above the station reference point. Whether the supplied
CRD pressure was already reduced to that point remains undocumented; the
corresponding uncorrected height sensitivity is 0.77--0.91 mm per leg. No guessed
pressure reduction or change to the physical IERS mapping function is applied.

NASA's `archive release notice
<https://earth.gsfc.nasa.gov/geo/projects/sgp/apache-point/updates/complete-dataset-apache-point-lunar-laser-ranging-station>`_
documents reprocessing with ACS calibration in February 2024. A separate
data-version audit retrieved the current EDC CRDv2 files for the same four
reflectors on 2017-03-09, created 2025-06-26 and received 2025-06-27. These files
are served under ``quarantine/7045/2017``. Their documented
`WGL021 warning <https://edc.dgfi.tum.de/en/oc/crd/v2/#WGL021>`_
means that EDC lacks station XYZ coordinates for calculating an orbit fit.
Their quarantine status is retained; they are not substituted into the frozen
baseline or evaluated as a new accepted baseline.

The directory ``analysis/laser_ranging/results/apollo_20170309_release_comparison``
records source URLs and SHA-256 hashes, and every matched point's exact decimal
tags, RTTs, weather, photon counts and provider uncertainty. Matching uses only
reflector and acquisition-window overlap. Ten original points retain exactly
the same transmit tag; their RTT revisions correspond to half-range changes
between -3.687 and +12.052 mm. The third A15 observation moves from
01:56:10 to 01:56:15 UTC, so its RTTs are not subtracted. The current files also
contain an additional eight-photon A14 point at 02:05:05 UTC. It is recorded
only as unmatched release metadata. All eleven original points remain intact.

The current headers specify transmit epochs, station delay already applied,
and no troposphere, centre-of-mass, receive-amplitude or spacecraft-delay
correction. The legacy fixed-width rows do not expose equivalent CRD flags,
so identical correction status cannot be asserted. Current normal-point
uncertainties are taken from explicit provider comments, which retain a caveat
about a possible later scale factor or RSS term; CRD bin RMS is kept separately.
Weather and photon counts also differ. These are release differences, with no
residual calculation, data cuts or preferred-release selection.

Primary sources
---------------

* `Park et al., DE440/441 (2021) <https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf>`_
* `Williams and Boggs, JPL LLR Model 2020 <https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf>`_,
  Sections 3.11–3.18 and Appendix A.
* `APOLLO normal-point format and release <https://tmurphy.physics.ucsd.edu/apollo/norm_pts.html>`_
* `APOLLO survey note <https://tmurphy.physics.ucsd.edu/apollo/doc/APO_coords.pdf>`_
* `INPOP17a report <https://www.imcce.fr/content/medias/recherche/equipes/asd/inpop/inpop17a_v3_090817.pdf>`_,
  Tables 4 and 8; source hashes and interpretation in ``llr_station_profiles.json``.
* `Manche (2011) reduction model <https://www.imcce.fr/content/medias/publications/publications-recherche/theses-habilitations/docs/Manche_These.pdf>`_,
  Sections 9.1.1–9.1.5.
* `APOLLO data characterization (2023) <https://arxiv.org/pdf/2304.11128>`_ and
  `calibration (2023) <https://arxiv.org/pdf/2304.11174>`_
* `IERS chapter 9 and official FCUL routines <https://iers-conventions.obspm.fr/chapter9.php>`_
* `Turyshev et al., relativistic lunar transformations <https://arxiv.org/html/2406.16147v4>`_,
  Eq.36 for the size of Earth's omitted acceleration terms; its newer lunar
  coordinate convention is not silently substituted for JPL2020's model.
