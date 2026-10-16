LLR runtime audit: proposals only
================================

No implementation, interpolation, model or tolerance was changed. Evidence uses
the existing March 9 control and two complete normal points from its unchanged
source ``a453128ac8179fc9c3212648f78efc531ecc1114``, with native kernel
``80b924287707124a8efea66455b26233b5511b7856983d35c05d21f82734f215``
and one numerical thread. Native microprobes make 100 calls over a 1.3-second
window; these timings are indicative measurements, not precision benchmarks.

Measured costs
--------------

The thin-driver control reproduces all eleven raw residuals exactly: setup
0.0689 s, complete solves 0.6309 s, total 0.8970 s. The first long-arc attempt
completed all 69 APOLLO points over nine nights with setup about 0.594 s and
solves 3.925 s, then stopped on the Grasse 1064.2 nm input guard. This is partial
arc timing, not a completed 306-point run. A subsequent complete APOLLO-only
69-point run (source ``8d06418de``) took 0.6326 s setup, 4.2439 s solver and
5.3316 s total. These measurements indicate about 57--62 ms per point.
The completed 306-point two-station run (source ``6e7d4b807``, kernel SHA256
``cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef``)
subsequently took 26.6935 s total: 5.8127 s setup and 19.3589 s
solver across 79 groups. This is the unoptimized complete model, with the
separately reviewed actual-wavelength guard fix. Group count/setup contributes
materially on Grasse's smaller sessions; no additional timing probe was run.
The long-arc driver already evaluates the complete model once per point; the
four sequential component solves in the short diagnostic are not its cost.

Two-point cProfile wall time was 0.1162 s. Its 52 topocentric clock conversions
took 0.03236 s, about 28%; two simulator creations together took 0.000136 s,
under 0.2%. Python lunar-spatial callbacks accumulated 0.00125 s. Native
simulation self time was 0.0682 s; its internal clock/tide/orientation costs are
not individually resolved by Python profiling. Cumulative timings overlap and
must not be added to claim a total speedup.

Native advancing-epoch microprobes, microseconds per call:

========================================  ==========
Function                                  Time
========================================  ==========
Earth interpolated body state             0.50
Moon SPICE rotation                       2.85
Geocentric TDB to TT                      15.95
Earth rotation matrix                     46.74
Earth rotation derivative                 254.36
Topocentric TDB to TT (Time / double)     578 /584
Station linear-motion callback            580.97
Earth intrinsic position                  730.61
Earth final position                      1816.13
Moon intrinsic / final position           5.41/30.85
========================================  ==========

Earth matrix/rate timings at repeated identical epochs remained about 47/252
microseconds. The intrinsic Earth position includes solid/pole tides and the
linear station callback; it is not an isolated tide-only timing. The large
topocentric conversion cost is present for both Time and double inputs.

Ranked proposals, pending review
-------------------------------

1. Use geocentric TT only for the published linear station drift, as the simpler
   first candidate. The present callback performs the expensive topocentric
   inverse solely to multiply elapsed TT by a velocity of centimetres/year.
   At all 306 source epochs, replacing that TT by geocentric TT changes APOLLO's
   linear coordinates by at most 8.62e-16 m and Grasse's by 1.68e-15 m; maximum
   clock differences are 2.092 and 1.983 microseconds. The bound is
   ``norm(velocity) * abs(TT_geocentre - TT_station)``; it is a coordinate bound,
   not a measured range difference. The profile's J2000 TT/Julian-year
   interpretation remains explicit. Keep the generic converter, observable
   clock transfer and simultaneity exact. This targets the 28% clock-call cost
   visible here; replacement overhead and other native costs prevent claiming
   that speedup. No approximation has been implemented.
2. Reuse exact same-event results. The final Earth station evaluates its
   linear callback directly and again through the intrinsic station used for
   simultaneity: 52 kinematic callbacks versus 26 simultaneity callbacks in two
   points. Avoiding just the repeated half has a measured visible ceiling near
   14% of this profile, before overhead; eliminating all measured clock-call
   cost would be at most 28%. Use complete Time, station XYZ, converter/profile,
   origin and model identity as applicable cache keys. No rounded-epoch snapping
   or reuse across transmit/receive events. Reuse R/Rdot or intrinsic station
   results at identical events before introducing approximations. Simulator
   reuse by station/reflector/wavelength is possible, but its measured ceiling
   is small here and cannot address the dominant cost.
3. Investigate a safeguarded faster topocentric inverse. The current converter
   brackets raw EOP intervals and bisects in Time with a 1e-12 s residual stop.
   A bracketed Newton/secant method using justified clock derivatives, retaining
   bisection fallback, is a candidate. Preserve coverage, leap/EOP boundaries,
   monotonicity, full Time and the existing residual/convergence tolerances.
   This targets the measured clock path; no aggregate speedup is established.
4. Only after measurement, consider complete Earth rotation/state or intrinsic
   displacement tables. Existing CIO/TDB-TT generation settings are available,
   but CIO interpolation alone leaves four direct SOFA samples in each nominal
   CIP derivative. The existing geocentric clock table does not remove the
   measured station-position-dependent inverse. Native tabulated rotation uses
   quaternion/angular-velocity interpolation, but lacks an exposed precise-Time
   replacement path and does not preserve the GCRS model's converter lookup
   automatically. Python custom station motion also reduces Time to double and
   supplies zero additional velocity; it is not a full state/rate replacement.
   These constraints require design review rather than a driver substitution.
5. Independent station/night multiprocessing is a later option. Preserve one
   process-local SPICE/model/converter and deterministic point identities/output
   order, with numerical threads capped. Current APOLLO work is only seconds;
   launch/setup/I/O costs may offset gains. No parallel speedup was measured.

Any approximation needs a separately declared acceptance test: a 0.1 mm Earth
terminal budget corresponds to about 1.6e-11 rad (3.2 microarcseconds), and 1 mm
to about 32 microarcseconds. Validate full rotations and derivatives against the
direct model at withheld precise-Time events, including two epochs sharing one
double, quaternion sign/orthogonality and derivative consistency. Keep the ICRF
identity, exact clock transfer, separate retarded event epochs, padded boundaries
and all EOP/leap/station discontinuities. Validate intrinsic displacement and
full line-of-sight half-range errors separately at both endpoints; never apply
the original tide/spatial/simultaneity correction a second time. Knot spacing
and achievable accuracy remain unmeasured proposals, not accepted settings.

Source locations and retained evidence
-------------------------------------

* ``analysis/laser_ranging/llr_residuals.py``: station propagation, Earth
  simultaneity, lunar spatial callbacks, ``computed_range`` and model creation.
* ``analysis/laser_ranging/llr_long_arc.py``: group setup and per-point solver
  timers; source a453128ac. No timing labels were inserted by this audit.
* ``terrestrialTimeScaleConverter.cpp``: ``getTaiTime`` and
  ``invertTimeFromTai``; nonzero-position TDB uses raw-interval bisection.
* ``precessionNutationCalculator.cpp``: ``getNominalPositionAndRate``;
  ``itrsToGcrsRotationModel.h``: complete matrix derivative.
* ``createRotationModelEarthOrientation.cpp``: existing CIO/clock interpolation;
  ``tabulatedRotationalEphemeris.h`` and ``createGroundStations.h``: Time/rate
  limitations of possible existing table/callback routes.

Full profiling output, microprobe JSON and the measurement harness remain
external under ``research/llr-runtime-audit`` and ``verification``. This report
contains proposals only; no optimization was applied.
The measurement JSON SHA256 is
``6fa3974c713035535077261b369ffcd6d31261567e458ec608d54eb3cd5e1d49``;
the clock microprobe JSON SHA256 is
``2d625a7ac35930fdd1decfc3b27d33bc1a689d033e579fe99469ac883cef21aa``.
The coordinate-only 306-epoch probe is
``verification/audit_grasse_wavelength_and_drift_clock.py``, SHA256
``f16653da8d86aa8c5e31d8d1522de4bf74d0f3d6f77559bff6a2a49535672e98``;
its retained JSON SHA256 is
``2fd9df47e59c76aebb4207d53cb68283b3a16cfeac7a7ab768109b3ae3cf18bc``.
The active model source SHA256 is
``60a790cb9b25c6f45881afbbc62475a844a5e3e2eac2dfadaf0940ce9db49e2a``,
and the frozen long-arc manifest SHA256 is
``10906a8e5ebbf349b7d684844d968690cc500640aae0a634f6f82e722e72c86f``.
