Relativistic station spatial approximation
=========================================

The Earth convenience ``earthRelativisticStationMotionSettings`` implements
IERS Conventions (2010), Chapter 11, Eq. (11.19), from a TT-compatible local
station position to a TDB-compatible barycentric spatial offset. The generic
settings require explicit coordinate compatibility: unscaled TCG/TCB Eq. (11.18)
or TT/TDB Eq. (11.19). Omitting compatibility is rejected rather than guessing.
All coordinates and gravitational parameters remain in SI units.

**This is a spatial approximation evaluated at the supplied geocentre TDB epoch.
It is not a full IAU event transformation or a station worldline at a supplied
barycentric event time.** In particular, the station-dependent simultaneity term
relative to geocentre TT, approximately ``-v dot r / c²``, is omitted. For a
moving surface station its position effect can be about a millimetre. The full
IAU B1.3 spatial acceleration terms are also omitted. A velocity returned by this
model differentiates the declared approximation and does not close those omitted
parts of the relativistic event transformation.

With ``r`` the full intrinsic geocentric position in the rotation's inertial base
axes, ``v`` the host barycentric velocity in those same axes, and ``U`` the
external potential evaluated at the host centre, the spatial map is

::

    r_B = (1 - U/c² - L_C) r - 0.5 v (v dot r)/c²

``L_C`` is zero for unscaled compatibility and ``1.48082686741e-8`` for TT/TDB
compatibility. The correction is returned in body-fixed axes; the reference-point
factory rotates it and adds the host translation once. Native Earth rotation
models are viewed through TDB even when their standalone clock is TT, TAI, UTC
or UT1. Other rotation models retain their supplied clock contract.

Potential provenance and composition
------------------------------------

The external potential is the sum of ``GM_J / |x_host - x_J|`` over precisely the
configured body names, excluding the host. The Earth convenience defaults to
``["Sun"]``, an explicitly Sun-only potential approximation; pass the desired
external-body list to include other bodies. The model retains this list and its
callbacks through ``getExternalBodies``. Missing bodies, missing GM, duplicate
names, the host itself, zero separation and non-finite values are rejected.
Gravitational parameters are assumed constant. Use the environment-aware public
``createGroundStation`` overload to supply ``SystemOfBodies`` for those callbacks.
The external-body callback binds each selected external body, not the station host.

Intrinsic linear motion, displacements and body deformation are combined first.
The relativistic map then consumes that complete position/state exactly once,
regardless of its place in the settings list. Nested aggregators are flattened
and multiple relativistic maps are rejected. For a non-SSB target origin this
barycentric spatial correction is inactive.

Position-only access evaluates the central position, host velocity, rotation and
potential. It does not request intrinsic station velocities, orientation rates,
external-body velocities or host acceleration. Those unavailable rates therefore
do not remove an otherwise available position.

State derivative and fallback
-----------------------------

Write the correction in body-fixed axes as ``K q``. Then its local velocity is
``Kdot q + K qdot``. ``Kdot`` includes the rotation derivative, host acceleration
and external potential rate. The latter uses the supplied host and external-body
velocities analytically; native Earth orientation rates use the selected EOP
and clock policies. The complete intrinsic velocity is already a TDB derivative,
so no second clock factor or second intrinsic displacement is added.

Host acceleration is obtained from the supplied velocity callback with a
five-point numerical stencil. Legacy callback-only rotations and external
positions use the same fallback when analytic rates are absent. It compares
32-second and 16-second stencils, checks finite/distinct/representably uniform
samples and convergence, propagates callback/coverage errors, and restores
central-epoch caches while preserving the original failure type. Callbacks must
be smooth throughout the entire stencil; convergence alone cannot prove a black
box callback has no discontinuity. The fallback is a numerical derivative of
the supplied callback, not an invented analytic host dynamics model. Unsupported
or unavailable native orientation rates fail rather than using zero or replacing
the selected policy. Non-finite required velocities and rate outputs also fail.

Reference-point ephemeris metadata now records the rotation's actual inertial
base axes and the supplied global origin. ``createReferencePointEphemerisFromId``
passes the environment origin. Direct calls infer an origin only from an
unconverted ephemeris (or the host's explicit origin identity); after an origin
conversion, callers must supply the optional ``frameOrigin`` argument.

Validation
----------

``test_astro_ground_stations_RelativisticStationPosition`` checks original
70-decimal-digit Eq. (11.18)/(11.19) fixtures and potential rates, the exact
flat-spacetime Lorentz contraction of a stationary comoving point, terminal
composition/order and duplicate rejection, finite differences of a moving
station with changing host motion/potential/orientation, position-only access
with unavailable rates, coverage/representability failures, and native public
factory wiring for five Earth clocks with distinct Earth/Sun/Moon states.
Fixture inputs, source hashes and simultaneity conventions are recorded in
``tests/test_tudat/reference_files/relativistic_station/README.md``.

Primary sources:

* IERS Conventions Chapter 11, Eq. (11.18)/(11.19):
  https://iers-conventions.obspm.fr/content/chapter11/tn36_c11.pdf
* IERS Conventions Chapter 10, event/time-coordinate transformation:
  https://iers-conventions.obspm.fr/content/chapter10/tn36_c10.pdf
* Soffel et al. (2003), Appendix A, IAU 2000 B1.3 and Section 3.4:
  https://arxiv.org/abs/astro-ph/0303376
