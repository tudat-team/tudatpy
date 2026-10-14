Earth rotation values and physical time derivatives
==================================================

The native ``GcrsToItrsRotationModel`` differentiates its complete position
rotation, with respect to the input clock declared at construction. The inverse
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

``EarthOrientationAngles`` holds the common TAI, TT, UT1 and five small angles.
``EarthOrientationAnglesAndRates`` adds their input-clock derivatives. The
original double input is canonicalized through ``getTaiTime`` before conversion
to ``Time``; TT and UT1 then retain extended precision throughout rotation and
rate evaluation. The extended-time virtual derivative APIs are overridden.
Intermediate-frame rotations use the same TT/UT1 evaluation.

Daily UT1-TAI, polar motion and observed celestial pole rates come from each
selected EOP interpolator. At a knot they use the right interval, except the
final sample uses the left interval, with the selected four-point Lagrange
window when requested. An unavailable slope, rejected gap or missing required
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
representation itself cannot label an inserted second uniquely. A selected
custom TT/TDB clock interpolator has positions but no native clock-rate API.

Numerical components and error budget
------------------------------------

Short-period UT1 and polar motion differentiate their actual harmonic sums along
the local TT/UT1 tangent with a five-point stencil, default step 32 s. The daily
interpolation slope is never estimated by probing neighboring raw EOP. UT1's
short-period correction is evaluated at baseline UT1 in the position path, so
its rate uses that same baseline tangent. Polar motion uses final corrected
UT1 and its full rate.

Nominal SOFA X/Y/s use a five-point derivative with default step 32 s. The step
is available through ``getNominalDerivativeTimeStep`` and configurable through
``setNominalDerivativeTimeStep``; nonpositive, nonfinite, or numerically
unresolvable stencils are rejected. If nominal-CIO interpolation is configured,
its selected linear, cubic, Hermite or Lagrange polynomial branch is
differentiated consistently. Polynomial recovery samples only inside that
interval. Knot derivatives use the right branch, except the final endpoint
uses the left branch; configured cubic Lagrange boundary fallback is respected.
Rates outside nominal interpolation coverage or for unsupported interpolation
policies throw, independently of the position extrapolation setting.

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

Legacy constructor policy
-------------------------

Legacy/custom constructors remain usable. Their rate fallback differentiates
the actual supplied position/clock path with a 32 s five-point stencil. It
supports the known daily linear, jump-linear, cubic, Hermite, Lagrange and
piecewise-constant interpolator types only. It requires the entire stencil to
lie strictly within one daily interpolation interval (and one nominal-CIO
interval, when present). It rejects UTC jumps, interval knots and coverage
boundaries, and pre-1972 legacy clock conversions. These restrictions are
explicit; the fallback is not the strict evaluator's analytic daily-rate path.
Custom argument functions remain responsible for smoothness inside their
reported interval. The legacy orientation consumer updates all clocks together
to prevent a previous epoch's cached UT1 from surviving a partial TT update.

Native validation
-----------------

``test_earth_orientation_EarthRotationDerivatives`` checks independent SOFA
positions and rates for IAU 2000A, 2000B and 2006 and all supported base axes;
exact celestial-matrix orthogonality; controlled spin, PM and CPO slopes;
forward/inverse rates and the rotation tangent identity; double and ``Time``
routing; leap-second continuity; UTC drift/discontinuities; strict endpoints,
knot sides and missing rates; native short-period convergence; selected nominal
interpolation; and the guarded legacy fallback. Complete deformed-station state
closure is validated separately under task 12.
