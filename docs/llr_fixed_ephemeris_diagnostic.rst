Fixed-ephemeris lunar laser ranging diagnostic
==============================================

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

Primary sources
---------------

* `Park et al., DE440/441 (2021) <https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf>`_
* `Williams and Boggs, JPL LLR Model 2020 <https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf>`_,
  Sections 3.11–3.18 and Appendix A.
* `APOLLO normal-point format and release <https://tmurphy.physics.ucsd.edu/apollo/norm_pts.html>`_
* `APOLLO survey note <https://tmurphy.physics.ucsd.edu/apollo/doc/APO_coords.pdf>`_
* `APOLLO data characterization (2023) <https://arxiv.org/pdf/2304.11128>`_ and
  `calibration (2023) <https://arxiv.org/pdf/2304.11174>`_
* `IERS chapter 9 and official FCUL routines <https://iers-conventions.obspm.fr/chapter9.php>`_
* `Turyshev et al., relativistic lunar transformations <https://arxiv.org/html/2406.16147v4>`_,
  Eq.36 for the size of Earth's omitted acceleration terms; its newer lunar
  coordinate convention is not silently substituted for JPL2020's model.
