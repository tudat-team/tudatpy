GROOPS comparison for Earth rotation and station modelling
========================================================

Scope and provenance
--------------------

This source review was requested by Dominic on 14 October 2026. The inspected
`GROOPS revision
<https://github.com/groops-devs/groops/tree/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee>`_
is ``7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee``. It is a useful independent
implementation and a source of test scenarios. This review has not yet run
GROOPS or established numerical agreement. IERS conventions, identified
official updates and original reference routines remain the conventional
baseline. Existing checklist task numbers are retained.

The comparison concerns station positions and Earth orientation. Gravity
recovery and other GROOPS applications are outside the current programme.
GROOPS sources remain in a separate reference checkout; Tudat implementations
are developed from the scientific specifications. Any later executable
benchmark must pin its input data, build options and conventions as well as
the code revision.

Earth orientation and derivatives: tasks 07, 09--13, 25, 27
----------------------------------------------------------

The `IERS 2010 implementation
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/classes/earthRotation/earthRotationIers2010.cpp>`_
uses cubic interpolation, removes leap jumps from the interpolated UT1
ordinate, adds ocean/libration corrections, and uses ERFA precession/nutation
with observed celestial pole offsets. Its interpolation abscissae and several
high-frequency argument calls use UTC. Tudat's continuous TAI interpolation
and explicitly routed TT/UT1 arguments therefore require a documented
comparison policy, especially at leap seconds. Retain independent IERS/SOFA
tests for these differences; do not tune tolerances to hide them.

The `rotation-rate implementation
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/classes/earthRotation/earthRotation.cpp>`_
normally differentiates complete rotations with an eighth-order stencil at
30-minute spacing. If evaluation fails, it uses an LOD-based approximation.
This supplies an interior cross-check, with a convergence study. Tudat's rate
must instead follow the exact selected position interpolant and report
unavailable derivatives explicitly. Test changing polar motion, CPO, UT1
slope, interpolation knots, coverage endpoints and full station velocity
closure; a rotation-matrix match alone does not validate velocity.

The separate `high-frequency EOP configuration
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/classes/earthRotation/earthRotationIers2010b.cpp>`_
provides a useful later comparison for the optional Desai--Sibois model.
Keep that choice separate from the conventional default, and test that an
alternative ocean contribution replaces the corresponding default exactly
once. Task 13 remains deferred until the standards baseline is validated.

`GNSS Earth-orientation partials
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/gnss/gnssParametrization/gnssParametrizationEarthRotation.cpp>`_
and the corresponding SLR class separate spatial derivatives from temporal
basis functions. Include tests for angular versus UT1 units, independently
perturbed complete observations, temporal support boundaries and parameter
ordering. Approximate matrix expressions in another implementation are not
an oracle for exact Tudat derivatives.

Tides and loading: tasks 06--07, 14--17, 21--22, 26--28
------------------------------------------------------

`Pole tides
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/classes/tides/tidesPole.cpp>`_
and `ocean pole tides
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/classes/tides/tidesOceanPole.cpp>`_
use an explicit secular-pole model. This reinforces tests of zero wobble,
the sign of the second wobble component, angular units and the local-to-
Cartesian basis. GROOPS computes ocean pole loading from spherical harmonics;
Tudat currently uses the official displacement grid. Comparisons must align
the load Love numbers, truncation, constants and reference origin before
assigning a numerical tolerance.

`Harmonic loading
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/classes/tides/tidesDoodsonHarmonic.cpp>`_
supports selected constituents, explicit admittance and nodal-correction
choices. For Tudat, retain isolated-constituent sign/phase tests, major/minor
constituent identity checks, and tests preventing duplicate nodal or minor-tide
corrections. A different loading synthesis is a physical-model comparison,
not an exact HARDISP-equivalence test. GROOPS's spherical-harmonic path does
not establish the latitude convention of the Bernese RP03 displacement grid.

The `GNSS station loading assembly
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/gnss/gnssReceiverGenerator/gnssReceiverGeneratorStationNetwork.cpp>`_
distinguishes deformation Love numbers from potential Love numbers when a
field includes the induced potential. Record whether a provider product is
a displacement, loading potential or total potential, and its reference
origin and included signals. Tasks 21--22 must reject incompatible or duplicate
contributions. Test tide/loading sums component by component before combining
complete station trajectories.

Station histories and instrument geometry: tasks 18--20, 24, 27--28
------------------------------------------------------------------

`SINEX station trajectories
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/programs/conversion/gnss/sinex2StationPositions.cpp>`_
combine solution intervals, velocities, discontinuities, logarithmic and
exponential post-seismic terms, and annual/semiannual terms in XYZ or local
coordinates. This is a valuable complete-trajectory comparison target.

Tudat tests must retain station/point/solution identity, support reordered
parameter records, distinguish a 365.25-day velocity unit from calendar-year
seasonal phase, and exercise multiple earthquakes and solution changes.
Validate analytic rates, continuity where the model prescribes it, and
one-sided behaviour at actual discontinuities. Keep CM/CF seasonal products
and potential overlap with non-tidal loading explicit.

The inspected converter repairs some overlaps, fills gaps and can extrapolate
or extend seasonal terms. These are processing policies, not source-data
facts. Tudat's input layer must preserve the original intervals; consumers
must make any such extension explicit. Include rejection tests for ambiguous
solutions, gaps, unsupported extrapolation and mismatched station identities.

`SLR eccentricity histories
<https://github.com/groops-devs/groops/blob/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee/source/programs/conversion/slr/sinexEccentricties2SlrPlatform.cpp>`_
retain equipment validity and distinguish UNE components from Cartesian
coordinates. Add dated equipment changes, local-axis order, marker-to-
instrument reference-point offsets and overlap diagnostics to the existing
task 24 tests. Preserve these histories in neutral input records.

Planned independent execution
-----------------------------

Task 27 will run a pinned GROOPS configuration after the required Tudat
components are validated. Start with rotation-only and individually enabled
station terms, then complete trajectories at shared epochs and conventions.
Use matched inputs and explicitly list any algorithm differences. Save the
configuration, source/data hashes, generation command and results as offline
fixtures. Normal C++ regression tests must not require GROOPS, network access
or a Python-kernel build. Reserve build threads within the existing global
eight-thread limit when producing the reference results.
