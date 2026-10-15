Earth rotation values and physical time derivatives
==================================================

The native ``GcrsToItrsRotationModel`` computes the rotation rate with respect
to the input time scale declared at construction. The inverse
rate is the transpose of the forward rate. The constant GCRS-to-base matrix
premultiplies both. GCRS, SPICE J2000 and ICRF have identical axes; EME2000 uses
the selected SOFA frame bias and ECLIPJ2000 uses the fixed SPICE ecliptic axes.

The complete derivative is
``Qdot*A*W + Q*Adot*W + Q*A*Wdot``. Elementary rotations and the celestial pole
matrix are differentiated analytically. The celestial pole matrix uses
``a = 1/(1+sqrt(1-X*X-Y*Y))``, rather than a truncated expansion, so the matrix
and its quaternion describe the same orthogonal rotation. ERA has rate
``2*pi*1.00273781191135448/86400 * dUT1/dInput``. There is no empirical ERA
phase shift or independent LOD insertion. LOD remains observed metadata:
using a different LOD series in the rate would not differentiate the selected
UT1 position interpolant.

Both the TIO locator s' and its rate use TT. Nominal CIO locator s is computed
before adding observed dX and dY, as in the
`SOFA Earth attitude cookbook, sections 5.5/5.6 <https://www.iausofa.org/s/sofa_pn_c.pdf>`_.
The matrix product, TT argument and TIO locator follow
`IERS Conventions 2010 chapter 5, equations 5.1, 5.2 and 5.13 <https://iers-conventions.obspm.fr/content/chapter5/icc5.pdf>`_.

Clocks, knots and availability
-----------------------------

``EarthOrientationEvaluatedAnglesAndRates`` holds TAI, TT, UT1, the five small angles and, when
requested, their derivatives with respect to the input time scale. A double
input is first converted to ``Time``; time conversion then proceeds through TAI.
TT and UT1 retain extended precision throughout rotation and rate evaluation.
The extended-time virtual derivative APIs are overridden.
Intermediate-frame rotations use the same TT/UT1 evaluation.

Daily UT1-TAI, polar motion and observed celestial pole rates come from each
selected EOP interpolator. At a knot they use the right interval, except the
final sample uses the left interval, with the selected four-point Lagrange
window when requested. Node comparisons admit differences up to
``2*epsilon(long double)*TIME_NORMALIZATION_TERM`` seconds, accounting for
rounding in the fractional-hour part of ``Time`` when the same epoch is formed
through different arithmetic. This representation budget is separate from
the data-coverage margin. An unavailable slope, rejected gap or missing required
CPO raises an exception; it is never silently replaced by zero. The explicitly
selected missing-CPO zero policy still applies. Covered raw positions remain
usable when a neighboring gap makes the rate unavailable. Historical Delta-T
position conversion remains supported where selected, but its rate is not
provided by this implementation.

Clock factors are one for TAI and TT, the native geocentric dTT/dTDB for TDB,
and the inverse modeled dUT1/dTAI for UT1. Scalar UTC uses the SOFA DAT drift
within a civil-date segment, including pre-1972 drift. Its derivative is
undefined at an offset or drift boundary, including midnight immediately after
a leap second, and raises an exception there. This restriction does not affect
continuous TAI/TT/TDB/UT1 rotation through inserted seconds. The scalar UTC
representation itself cannot label an inserted second uniquely. The geocentric
TT/TDB rate also supports a selected TT/TDB interpolator: a five-point stencil
with fixed step 32 s differentiates its converted TT-minus-TDB values. The
stencil queries obey the selected interpolator's boundary policy.

Numerical components and error budget
------------------------------------

Short-period UT1 and polar motion differentiate their actual harmonic sums along
the local TT/UT1 tangent with a five-point stencil, default step 32 s. The daily
interpolation slope is never estimated by probing neighboring raw EOP. UT1's
short-period correction is evaluated at baseline UT1 in the position path, so
its rate uses that same baseline tangent. Polar motion uses final corrected
UT1 and its full rate.

Nominal SOFA X/Y/s use a five-point derivative with fixed step 32 s. If
nominal-CIO interpolation is configured, the returned positions use that
interpolator, but their rates are still computed from SOFA. The interpolating
polynomial is not differentiated. Nonfinite nominal rates are rejected.

For the nominal-coordinate step comparison, the test budgets eight double
machine epsilons of absolute error from unit-scale SOFA matrix arithmetic per
returned coordinate. The five-point weights amplify this by ``3/(2*h)``.
Comparing two stencils adds their budgets; the vector bound also includes
``sqrt(3)``. A conservative nominal amplitude of 1e-3 rad at the shortest nominal
CIP period of two days gives the truncation allowance
``amplitude*(2*pi/period)^5*h^4/30``. This is a numerical budget, not an accuracy
claim for the physical IAU model. The default 32 s roundoff allowance corresponds
to approximately 1e-9 m/s at Earth radius, below the complete station-velocity
closure target of 1e-7 m/s. Tests compare 8, 32 and 128 s steps and an independent
SOFA call using split dates at 64 s, across 1962, 2000, 2017 and 2026.

The independent complete SOFA matrix-sequence test uses explicit affine UT1,
polar motion and celestial pole inputs, SOFA ``iauC2ixys``, ``iauPom00`` and
``iauC2tcio``. It does not use the production quaternion or EOP evaluation path.
At one Earth station the 512, 128, 32 and 8 s five-point velocity errors are
approximately 1.97e-5, 7.67e-8, 7.42e-10 and 2.67e-9 m/s respectively. The coarse
steps show fourth-order convergence; the smallest steps show the expected
roundoff floor. A separate production-position velocity closure is also checked.

Construction
------------

The polar motion, precession-nutation and time scale objects are always
constructed from the same interpolated EOP data (``EarthOrientationDataEvaluator``),
so every angle and rate is computed in one way.

Native validation
-----------------

``test_earth_orientation_EarthRotationDerivatives`` checks independent SOFA
positions and rates for IAU 2000A, 2000B and 2006 and all supported base axes;
exact celestial-matrix orthogonality; controlled spin, PM and CPO slopes;
forward/inverse rates and the rotation tangent identity; double and ``Time``
routing; leap-second continuity; UTC drift/discontinuities; strict endpoints,
knot sides and missing rates; native short-period convergence; selected nominal
interpolation. Complete deformed-station state
closure is validated separately under task 12.
