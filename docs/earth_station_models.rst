Earth rotation and station modelling implementation
==================================================

Starting point and decisions
----------------------------

This work starts at PR #900 (``feature/data-refactor``), commit
``3070b20e83e18a749767fe0332a5280464f59b16``.  Changes are developed on
``fix/earth_deformation_rotation`` and committed as validated units. The session
was recovered on 14 October 2026 after an interrupted context compaction.
Work in the recovered session remains local; remote publication requires an
explicit instruction in the current conversation under Dominic's persistent
workspace policy.

* First reproduce accepted standards and correct existing models.  The
  published IERS Conventions (2010), explicitly identified official updates,
  and pinned reference implementations define the conventional baseline.
* SPICE ``J2000`` denotes ICRF axes.  It must not acquire an extra dynamical
  J2000 frame-bias rotation.  EME2000, when requested, is a separate explicit
  frame.  Frame orientation does not by itself specify origin or relativistic
  coordinate scaling.
* Data retrieval, parsing, validation and source metadata belong to
  ``data_input``.  This layer must have no dependency on ``dynamics`` or
  ``estimation``, including transitive imports, C++ includes and public types.
  Model construction and evaluation belong in consuming layers.  Shared
  time, units and mathematical utilities must be independent of both.
* Reference programs may be downloaded and executed to generate independent
  fixtures.  Record the source URL, version, checksums, generation command,
  input conventions and applicable licence.  Normal regression tests must
  use local fixtures and must not download data.
* All concurrent builds together are limited to eight worker threads.
  The coordinator assigns explicit allocations before concurrent builds and
  controls any reallocation across the integration build and isolated workers.
  Disable hidden nested build parallelism and multithreaded numerical runtimes.
* Every change includes relevant bindings, documentation and tests.  Tolerances
  come from reference accuracy and numerical error budgets, not from fitting
  the implementation's current output.
* Unit tests primarily use C++.  Python tests are reserved for Python-specific
  imports/bindings or cases that cannot straightforwardly be exercised in C++.
  Reference-data generation scripts may use Python to drive external programs.
  Compile only the native tests highly relevant to each change. Record unrelated
  unexpected failures for final integration instead of broadening every build.
  Reuse passing isolated-worktree results after integration unless a concrete
  merge, dependency or behaviour change creates a real risk that needs a rerun.
  Compile the Python kernel only to test Python-only functionality, such as
  binding registration, Python imports or Python-specific API behaviour. A full
  kernel build is not a general C++ integration requirement.
* Implementation agents work in separate worktrees and separate build directories.
  The coordinator integrates reviewed changes into the authoritative worktree.
* SINEX readers and ILRS interfaces do not need backward compatibility, as
  instructed during session recovery. Prioritize strict, complete history
  semantics and update affected consumers consistently.
* Do not add estimation of new parameters. Earth-orientation and station-motion/
  loading estimation extensions and synthetic parameter recovery are excluded.
  This does not exclude physical velocities or time derivatives of models.

Data architecture
-----------------

The Python entry points for Earth data live beneath
``tudatpy.data_input.environment_data``.  C++ parsers and source-data records
live in ``tudat::input_output`` beneath ``tudat/io``.  The records carry values,
units, epoch/time-scale identifiers, validity, uncertainties where supplied,
reference-frame/tide-system identifiers and source provenance.  They do not
contain bodies, rotation models, observation models or estimatable parameters.

Dynamics factories consume these records and own interpolation policies,
model selection and construction.  Estimation consumes the resulting models
and provides parameterisation and observation partials.  Existing convenience
interfaces may delegate loading to the input layer; they must not cause a
reverse dependency.  Binding registration and shared utility dependencies
must also respect this direction.

Task checklist
--------------

Unchecked items remain open even when a prerequisite or an initial subset has
been committed.  Larger items may have several focused commits.  The
validation ledger below records completed changes and remaining limitations.

Status: 17 of 25 tasks done; 3 active, 2 partial and 3 queued.
Reviewed 14 October 2026. This list is updated as relevant work is implemented,
validated and integrated.

.. code-block:: text

   [x] 01 Establish the PR #900 baseline, isolated build and task tracker.  [DONE]
   [x] 02 Define neutral data contracts and enforce the dependency boundary.  [DONE]
   [x] 03 Separate existing EOP and BLQ readers from model evaluation.  [DONE]
   [ ] 04 Establish pinned independent IERS/SOFA reference fixtures.  [PARTIAL — reference coverage grows with each model]
   [x] 05 Correct solid pole tides and validate spherical/ellipsoidal geometry.  [DONE]
   [x] 06 Validate solid Earth tides through the production time/argument path.  [DONE]
   [x] 07 Resolve inertial frames, time scales and displacement conventions.  [DONE]
   [x] 08 Implement robust legacy/current/rapid EOP ingestion with metadata.  [DONE]
   [x] 09 Implement EOP validity and leap-safe interpolation policies.  [DONE]
   [x] 10 Verify existing short-period EOP corrections independently.  [DONE]
   [x] 11 Complete Earth-rotation derivatives.  [DONE]
   [x] 12 Complete deformation velocities and full station-state closure.  [DONE]
   [x] 13 Add the optional Desai-Sibois high-frequency EOP configuration.  [DONE — opt-in only; replaces corresponding default ocean effects]
   [x] 14 Validate BLQ loading inputs and preserve their conventions.  [DONE]
   [x] 15 Implement HARDISP-equivalent ocean tidal loading.  [DONE]
   [x] 16 Add conventional ocean pole-tide loading.  [DONE]
   [x] 17 Add conventional atmospheric S1/S2 tidal loading.  [DONE]
   [x] 18 Preserve full SINEX solution/validity/eccentricity histories.  [DONE]
   [ ] 19 Implement post-seismic and seasonal station motion with rates.  [ACTIVE — explicit source clocks, PSD and seasonal evaluation]
   [ ] 20 Modernise ITRF transformations and station-catalogue inputs.  [PARTIAL — transformations validated and committed; catalogue integration pending]
   [ ] 21 Add non-tidal atmosphere/ocean/hydrology loading inputs.  [ACTIVE — strict provider inputs and metadata]
   [ ] 22 Integrate non-tidal loading with explicit overlap/validity policies.  [QUEUED]
   [ ] 23 Validate relativistic station position and velocity transformations.  [ACTIVE — spatial/event conventions and full-state composition]
   [ ] 24 Improve instrument reference-point and calibrated geometry support.  [QUEUED]
   [ ] 29 Consolidate public interfaces, diagnostics and reproducible examples.  [QUEUED]

Task numbers are stable identifiers rather than a rigid execution sequence.
Tasks 25--28 and 30 were removed at Dominic's request: new estimation
parameters/partials, cross-software complete-trajectory comparisons, observation
validation/controlled omissions and research-extension evaluation are outside
this programme. The retained scope comprises 25 tasks. Targeted independent
model-reference tests and physical state/time derivatives remain in scope.
Task 13 is explicitly authorized now as an opt-in; it must not become the
default. Selecting Desai--Sibois must replace the
corresponding default ocean-EOP contributions, with no duplicate application;
independent effects outside its scope remain explicitly selected.

Active work
-----------

Implementation workers use **GPT-6.1 Sol, high reasoning**, explicitly selected.
The coordinator handles architecture and integration. Current compiler allocations
are relativistic station transforms 4, station history 2, and non-tidal loading 2; the combined
limit remains eight. Earlier workers are stopped and their work is preserved.

* Task 04: independent reference packages are in use for solid Earth tides,
  short-period EOP and HARDISP; reference coverage grows with each added model.
* Task 11: complete rotation derivatives are integrated with all ten checkpoint
  files unchanged from the tested worktree. Thirteen focused native cases and
  334 assertions pass, including independent full SOFA matrix/state comparisons.
  The additive CMake overlap was resolved to the tested eleven-library list.
* Task 12: complete IERS solid-tide, HARDISP ocean, pole-tide and secular station
  state closure is integrated. Five focused cases/151 assertions pass, including
  the complete public station factory, nonzero body translation, four base axes,
  five declared rotation clocks, double/extended epochs, leap and coverage cases.
  The maximum native velocity closure error is 4.53e-9 m/s, below 1e-7 m/s.
* Task 13: Desai--Sibois is integrated as explicit opt-in. It replaces only ocean
  EOP terms and retains conventional libration once. The default is unchanged;
  ambiguous custom-file combinations reject. Seven native cases/276 assertions
  pass, with independent value/rate fixtures, full factory SOFA comparisons,
  provenance/dependency checks and syntax checks of both changed binding units.
* Task 17: complete atmospheric S1/S2 positions and velocities are integrated.
  All 19 relevant native atmospheric cases and the three UT1 prerequisite cases
  pass, together with input dependency checks. No Python kernel was compiled.
* Task 18: the neutral SINEX history parser, equipment records, full-history
  bindings and documentation are integrated. Sixteen native cases and 260
  assertions pass; all eight pinned complete IGN products parse. A malformed
  official IGS row is deliberately rejected and retained as a regression fixture.
  The merged CMake graph passes its raw-input dependency audit. Unchanged parser
  tests were not rerun. Compiled Python binding behaviour remains task 29.
* Task 19: station-motion evaluation is active in the history worktree. Physical
  validity uses TAI; the named IGN nominal-UTC profile preserves the published
  linear, post-seismic and seasonal source-time conventions and supplies rates
  only where their clock and validity semantics are defined.
* Task 20: neutral ITRF parameter ingestion and affine position/state transforms
  pass integrated comparisons against 75 independent high-precision fixtures.
  The transformation checkpoint includes bindings and API documentation.
  Station-catalogue inputs and their explicit physical-to-decimal-year mapping
  remain open with tasks 18 and 19; final compiled binding checks are task 29.
* Task 23: the rotation worker is now resolving the relativistic transformation
  contract and implementing correct spatial/state composition, external-body
  selection and explicit coordinate/time compatibility against primary standards.
* Task 21: strict EPHEDISP parsing and neutral preloaded loading-series records
  are being implemented from the pinned specification and a complete real IMLS
  provider file. Model evaluation, interpolation and overlap policy follow in 22.
* Task 29 will also simplify deformation-model entry points, following Dominic's
  AtmosphericTideDeformation example: consolidate redundant clock/geometry
  adapters and avoid convenience methods with subtly different output frames.

Other unchecked tasks are queued. Python-only binding/import validation is part of task 29;
kernel compilation is allowed only when needed for those Python-specific tests.

The pinned :doc:`GROOPS source comparison <earth_station_models_groops_review>`
records useful implementation and regression-test scenarios. Its former
estimation and complete-trajectory comparison proposals are outside the revised
scope; targeted reference checks remain useful.

Subsequent gravity-field programme
---------------------------------

Dominic requested this follow-on programme on 14 October 2026. It starts only
after **every retained item in the current 25-task checklist is complete**. Do not create its branch or begin its
implementation while any current item remains open.

At that point create ``fix/gravity_field_variations`` and review, implement
and independently validate detailed Earth gravity-field variation models used
in space geodesy, including solid Earth tides, ocean tides and pole tides.
Apply the same standards-first literature review, neutral input-layer
architecture, pinned reference data, primarily C++ validation and resource
limits. Estimation of gravity-field coefficients is explicitly outside this
follow-on programme's scope. Its detailed checklist will be established when
the current programme is finished.

Completion criteria
-------------------

For an implementation task, completion requires review of the final diff,
appropriate reference/regression tests, public-interface coverage where
affected, and a local commit.  Record pre-existing baseline failures separately.
Do not silently mark a model conformant when only equation-level or internal
consistency tests have run.  A numerical match to a specified model is not a
claim of equal physical accuracy for a real station.

Initial numerical targets, subject to the documented reference precision, are
1--10 micrometres for identical displacement algorithms, 1 microarcsecond for
equivalent rotations, and 0.05 microseconds for equivalent UT1 calculations.
Velocity closure initially targets 1e-7 m/s with a finite-difference convergence
study.  Actual discontinuities are tested as discontinuities.

Validation ledger
-----------------

2026-10-14, optional high-frequency EOP integration: task 13 commit ``8155c1145``
cleanly integrates ``1216f6bb1``. Three reader cases/32 assertions and four model
cases/244 assertions pass. Fourteen independent Decimal-polynomial phase/value
and analytic-rate fixtures cover 1962--2025; model tests include retained
libration, selected UT1/rotation rates and the full factory SOFA matrix sequence.
The embedded 159-term NASA mirror is byte-identical to its pinned 12586-byte
source. Three architecture/provenance checks and both binding translation-unit
syntax checks pass. The default conventional path remains unchanged; the explicit
Desai--Sibois selector replaces ocean contributions and rejects ambiguous custom
correction composition. Identical numerical tests were not rerun after integration.

2026-10-14, complete station-state integration: task 12 commit ``0eeccfcb2``
cleanly integrates ``5c1029c9a``. The station suite passes 5 cases/151 assertions;
the affected clock-routing suite passes 14/253 and the orientation suite 3/66.
The public reference-point factory now evaluates native Earth rotations using
TDB through the same angle calculator, preserving their EOP policies and axes
while differentiating with respect to the environment clock. Full intrinsic
motion, multiple deformations, rotation and nonzero body translation close
against independent position differences. The isolated build's stale borrowed
factory objects were rebuilt before validation. No post-merge numerical rerun
was needed because the implementation and dependencies were not altered by
integration. Evidence: ``rotation-derivative-build/task12-final.log``.

2026-10-14, rotation and station-history integration: task 11 commit
``e06188e6e`` preserves every file from the independently tested checkpoint
``431ed6160``. Its 13 cases and 334 assertions include a 7.42e-10 m/s full SOFA
station-velocity comparison at the selected 32-second stencil. Task 18 commit
``9719e2655`` preserves the strict parser, adapters, binding implementations and
fixtures from ``8bba93620``; its 16 cases/260 assertions, complete-provider-file
checks and fixture-regeneration evidence remain applicable. Conflicts only
combined additive library/module registrations. CMake regeneration checked the
combined target graph and neutral-input dependency closure. In accordance with
Dominic's latest instruction, identical numerical tests were not rerun after
integration; no Python kernel was compiled. Runtime Python binding checks are
reserved for task 29's Python-specific validation.

2026-10-14, session recovery: the original readable conversation, user requests,
last checklist and worker handoffs were recovered to
``/home/dominic/Tudat/earth-station-model-work/recovery-20261014``. Original session
records remain unchanged. Three workers survived interruption of their parent;
their turns were stopped through the local application server before new agents
took exclusive ownership of the existing worktrees. Their late edits were
preserved. The integration source matched all 38 frozen ITRF paths and passed
all 28 selected CTest checks plus both source-boundary checks. The ITRF
transformation checkpoint was subsequently reviewed and committed; its detailed
validation appears below. The independently tested UT1 value/rate prerequisite
was integrated separately. Tasks 11 and 17 remain open until their complete
models and their relevant tests pass. Compilation allocations are centrally
limited to eight jobs across all workers; Python kernel compilation is permitted
only for Python-specific validation.

2026-10-14: starting revision verified against GitHub PR #900.  Task 01 baseline
build completed with six build jobs: all 12 selected C++ test executables passed
(EarthOrientationCalculator, EopReader, PolarMotionCalculator,
TimeScaleConverter, ShortPeriodEopCorrections, ItrsToGcrsRotationModel,
Iers2010DeformationModel, OceanTideDisplacement, PoleTideDisplacement,
SofaEarthOrientationInterface, SofaFundamentalArguments, SofaTimeConversions).
Configuration used the local ``tudatpy-dev`` compiler/dependencies, Release with
``-O1``, serialization/PCH/explicit instantiations/extended precision disabled,
and no unity build or LTO.  These are targeted baseline results, not a complete
suite or a numerical-conformance claim.  Exact commands, compiler provenance,
resource hashes and logs are saved locally in
``/home/dominic/Tudat/earth-station-model-work/logs``.

Task 02 source-boundary checks pass for transitive local Python imports and C++
includes/associated implementations. The committed checks use only Python's
standard library and can run without a Tudat build. Independent input-library
linkage and compiled binding initialization also pass, as detailed below.

Task 05: corrected the pole-tide wobble units and east/north formulas, and use
the geocentric basis for the IERS displacement before projecting to a station's
topocentric frame.  The C++ PoleTideDisplacement executable passes eight frozen
reference cases evaluated independently from the original colatitude equations,
zero-wobble checks, and Cartesian/geodetic/spherical representation checks on
an oblate Earth.  The conventional secular pole remains the explicitly pinned
1 February 2018 IERS update.  Compilation also exposed and fixed a missing
direct constants include in the geodetic-coordinate header.

Task 07, inertial-frame portion: ``J2000`` and ``ICRF`` now share the GCRS axis
orientation, ``ECLIPJ2000`` uses the SPICE fixed ecliptic rotation, and the
previous dynamical-J2000 bias is available explicitly as ``EME2000``.  The C++
ItrsToGcrsRotationModel executable passes independent SOFA cookbook and fixed
frame-bias/ecliptic matrix cases, including forward/inverse transformations and
the double/extended-time interfaces.  Full time-scale and displacement-system
conventions remain open under task 07; derivative completeness is task 11.

Tasks 04/06: pinned and executed the unmodified corrected IERS DEHANTTIDEINEL
Fortran package, reproduced its two published examples, and generated 18
offline reference cases spanning station geometry, epochs and non-midnight
times.  Production solid-tide arguments now reproduce the reference's
precession-adjusted lunar longitude and UTC-hour phase; explicit custom
argument callbacks remain supported.  The C++ full-displacement comparisons
pass at 1 nm component tolerance through both the default production factory
and the independently implemented reference-argument callback.  The legacy
test's fractional-hour conversion was corrected; its rounded-century midnight
epoch was replaced by the exact published midnight without changing expected
values or tolerances.  The existing API's TDB-as-TT approximation remains
documented pending the broader time-scale work in task 07.  Task 04 remains
open for the other models' reference packages.

After these corrections, all 12 targeted C++ integration executables pass,
along with both source-dependency checks.  Python binding compilation follows
integration of the new EOP and ocean-loading input modules.

Task 10: executed pinned, unmodified PMSDNUT2, UTLIBR, ORTHO_EOP/CNMTX
and the Paris Observatory 71-wave ocean routine at 35 epochs spanning
1962--2025 and an 18.6-year nodal cycle. C++ tests distinguish exact reference
arguments from the production time/argument path, check default summed
corrections and single application of daily polar motion and observed CPO,
and cover positive amplitude cutoffs and empty results. Corrected a reversed
amplitude/multiplier argument order in coefficient filtering and updated the
legacy PMSDNUT2 check to the corrected 2011 reference values.

The published Table 5.1b T2 UT1 sine coefficient is 0.05 microseconds, whereas
UTLIBR uses 0.04. Tudat retains the published coefficient; the regression
explicitly accounts for that independently documented difference. The original
UTLIBR LOD convention is tested separately and is not represented as a new
production LOD API. Source hashes, licensing, generation instructions and
tolerance budgets are recorded in ``short_period_reference/README.md`` beside
the C++ suite. All nine short-period test cases pass; the complete set of
12 selected integration executables also passes after this change (1.63 s).
Broader interpolation, time-scale and rotation-rate work remains open.

Tasks 02/03/08/14: neutral EOP and BLQ records and readers now live in the input
layer, with separately linkable C++ targets and Python modules beneath
``data_input.environment_data``. EOP records retain SI values, optional missing
quantities, per-component Bulletin A status flags, source headers and supplied
frame metadata. Tests cover legacy 16-column and modern 21-column C04 plus
fixed-width finals2000A Bulletin A, using pinned real input extracts. BLQ records
retain the eleven standard constituents, up/west/south metre amplitudes,
degree phase lags, provider metadata and declared CE/CM origin. Malformed input
and duplicate station identifiers are rejected. Existing EOPReader/readBlqFile
interfaces delegate parsing to these readers; models can consume preloaded
records. Interpolation and full HARDISP synthesis remain separate tasks.

The raw I/O archive excludes the three adapters that depend on environment
models. GRAIL raw UTC-to-TDB conversions use the shared SOFA time utilities,
with C++ regression tests against the previous conversion results. A CMake
check traverses actual transitive links, including private static-library
dependencies, from all three input targets and rejects model dependencies.
All five reader/adapter test executables and all 12 selected scientific
integration executables pass, as do the two source-dependency tests and the
target-link audit. The EOP, ocean-loading and GRAIL suites contain seven, five
and three C++ cases, respectively.

The Python extension compiled successfully with input registration ordered
before dynamics and estimation types. Import/registration smoke checks and
both EOP binding-only Python tests pass against the newly built extension.
Four unused serialization-namespace aliases needed conditional guards for the
serialization-disabled build; that narrow build fix is a separate commit.
Subsequent full kernel builds are deferred to final integration at Dominic's
request. Logs are in ``earth-station-model-work/logs``; raw EOP validity and
the legacy adapter's artificial endpoint buffer remain task 09.

Task 15: the default ocean tidal-loading model now expands the eleven BLQ
constituents to the 342-constituent IERS HARDISP synthesis with complex admittance
interpolation. The eleven-constituent implementation remains an explicit legacy
option, and explicit custom argument callbacks retain that legacy behavior.
Filename and preloaded-data factories select the same model. Cartesian loading
uses the geocentric radial/north/east basis; public station ENU results are
projected into the station's actual topocentric frame.

The original corrected IERS Fortran package generates 36 frozen comparisons
covering six stations and non-midnight epochs over 25 years. Its regeneration
check passes, with source checksums and complete original notices preserved.
The largest C++ difference is 0.106 micrometres, within the independently
documented one-micrometre budget for the reference's single-precision synthesis.
The existing 24-hour Onsala and Reykjavik reference series now pass at one
micrometre in every component, replacing former 0.25--4 mm tolerances.

All six C++ ocean-loading cases pass, including factory selection, zero
amplitudes, legacy callbacks and ellipsoidal-coordinate geometry. The complete
15-executable scientific/input regression selection and raw-target dependency
audit pass (16 checks, 3.12 s). Logs are ``tests-hardisp-integration.log`` and
``hardisp-generator-check.log`` in the local work log directory. Binding source
and API documentation are included; final Python-kernel validation remains
task 29. Explicit environment TDB-to-model-TT routing and UTC-phase derivative
boundary handling remain tasks 07 and 12, respectively.

Task 09: production EOP evaluation now uses continuous TAI epochs and interpolates
UT1-TAI, with strict raw coverage, a default 36-hour gap limit and no artificial
endpoint padding. Linear interpolation remains the default; an explicit
four-point Lagrange option validates the entire selected stencil and its UT1
monotonicity. There is no silent extrapolation, order reduction or missing-CPO
fallback. Raw knot values remain available when a neighbouring gap prevents a
rate; interval rates use the right-hand interval, except at the final knot.

Twenty-seven frozen four-point comparisons were generated with the unmodified
Paris Observatory LAGINT and original SOFA date/UTC routines for modern epochs,
the 2016 inserted second and the 1962 drift era. Analytic cubic tests separately
validate slopes and reject an interior non-monotonic UT1 inverse. The original
LAGINT's undefined final-node branch is not executed; that endpoint is checked
against the raw record. Provenance and generation instructions are recorded in
``earth_orientation/interpolation_reference/README.md`` beside the C++ tests.

The rotation, polar-motion and time-conversion paths share the selected EOP
model. UT1 and observer-dependent TDB inversions stay within covered intervals.
Scalar UTC explicitly rejects the inserted second; continuous TAI/TT/UT1 paths
remain usable there. Exact rounded double images of coverage endpoints are
accepted without accepting later precise-Time epochs. Identity conversions
preserve their input without consulting unrelated tables. Geocentric atomic
clock conversion requires no EOP or UTC coverage. Historical Delta T is an
explicit clock-only policy, using the retained Julian-year table coordinate;
it does not invent historical polar motion, CPO or precise UTC.

The native geocentric TT/TDB rate differences the small SOFA clock offset rather
than large epochs. Independent 16/32/64-second convergence checks agree within
2.22e-16 in the rate. Selected clock interpolators retain position conversion
but explicitly reject unsupported rate evaluation. The final integrated EOP
suite passes 23 C++ cases, including exact endpoints, clock-interpolator coverage,
legacy custom paths and precise-Time representability regressions. Full Earth
rotation derivatives and deformation rate closure remain tasks 11 and 12.

Task 12, initial API: body deformation models now provide displacement velocities
and six-component states. The generic five-point difference uses a validated,
representable uniform stencil, preserves the original failure type and restores
central-epoch caches. Existing custom scalar callbacks retain their time basis.
Native solid/ocean-tide UTC phase models reject derivative intervals touching
an offset step or ambiguous inserted second; smooth historical UTC drift is
allowed. This is a phase-representation limitation, not a physical singularity.

Station-motion aggregation adds deformation velocities and propagates derivative
failures once the central displacement is available. Optional missing central
models can still be skipped. Position-only station access evaluates only the
central displacement, so an available position at a coverage boundary remains
usable when a derivative stencil is unavailable. Non-finite outputs are rejected.

All 12 integrated C++ velocity cases pass, including fourth-order convergence,
cache/exception behavior, large-epoch stencil representability, 2012/2016 leap
boundaries, historical UTC drift/steps and a controlled rotating-frame closure.
The IERS velocity changes by 1.22e-15 m/s between one- and half-second steps in
the smooth test case. This initial commit does not close task 12: native context
routing, selected-EOP-interval rates and closure with the complete physical
Earth-rotation derivative remain required under tasks 07, 11 and 12.

Task 16: the ocean pole-tide input reader preserves the dimensionless complex
radial/north/east grid, source headers and CM origin. Its independent input target
has no dynamics or estimation dependency. The model evaluates conventional IERS
Chapter 7 Eq. (29), the 2018 secular pole and IERS 2010 constants, with bilinear
geocentric-latitude interpolation and periodic longitude only for full global
grids. The production coefficient is K=5340.428562229677 m/rad. The official CM
grid already includes the geocentre contribution; no extra translation is added.
Geocentric displacements rotate to Cartesian before optional station-ENU output.

The raw grid stops at +/-89.75 degrees. Default evaluation requires grid coverage;
clamping to the nearest latitude is an explicitly selected approximation. Bad
headers, axis counts, duplicate/missing nodes and malformed rows produce source
and line diagnostics. Construction loads a file once through the input layer;
preloaded and file factories share the same evaluator and Earth polar-motion
model. Python input/model declarations and documentation accompany the C++ API.

All 731 published reference rows agree within 1 nm (maximum 0.7894 nm), using
only the published benchmark's explicitly older K=5339.4043696 m/rad and mean
pole for that comparison. These older constants do not change the production
model. The original full 720-by-360 grid was parsed independently; normal tests
use a checksum-pinned four-node extract and the unchanged published test series.
Source hashes, extraction instructions, constants and coordinate provenance are
recorded in ``reference_files/ocean_pole_tide/README.md``.

The final integrated ocean-pole model and input suites pass seven and four C++
cases, respectively, including file/preloaded parity, geometry, seam/coverage
policies and continuous-EOP evaluation through the 2016 inserted second. The
combined task 09/initial-12/16 milestone passes all 20 selected C++ executables
and the actual input-target dependency audit (21 CTest checks, 2.07 s), plus both
source-boundary checks. Exact-file formatting hooks pass. Evidence is recorded
in ``tests-combined-final-identity.log``, ``architecture-final-identity.log`` and
``combined-final-validation-summary.json`` in the local work log directory.
No Python kernel was rebuilt for this milestone; final binding validation is
still task 29.

Task 07: station motion now supplies a physical-epoch context that preserves
the original TDB input and lazily derives geocentric TT/TAI. Native solid tides
use TDB ephemerides, the rotational ephemeris's declared clock and TT/UTC tide
arguments. HARDISP and pole tides use their terrestrial clock explicitly.
UT1 and polar motion use the requesting model's own converter and coverage
policy. Generic/custom scalar callbacks retain their original forwarding.

Both native IERS and basic solid-tide factories honor a GCRS/ITRS rotation's
declared input time scale; SPICE and custom rotation inputs retain TDB. TT-only
displacement rates include dTT/dTDB, while mixed-clock solid tides differentiate
the complete routed displacement. Position evaluation remains independent of
rate availability. The public documentation states the conventional tide-free
nominal-coordinate requirement and the retained permanent degree-2 displacement.

All 14 integrated time-routing C++ cases pass, including direct-SOFA clock
fixtures, separate ephemeris/rotation/argument epochs, precise coverage limits,
leap-second behavior, custom callbacks, factory ownership and TT-versus-TDB
rotation equivalence. Existing solid/ocean/pole-tide and initial velocity tests
also pass in the combined 25-check integration selection. Full rotation rates
and station-state closure remain tasks 11 and 12. Logs are
``tests-time-routing-atmospheric-final.log`` and
``architecture-time-routing-atmospheric-final.log`` in the local work directory.

Task 17, position checkpoint: neutral atmospheric S1/S2 grid and station readers
preserve provider coefficients, source/origin metadata and explicit geometry
conventions. The separately linkable input target has no dynamics or estimation
dependency. Continuous UT1 drives the conventional midnight-referenced sine/cosine
synthesis. Grid lookup latitude and displacement basis are separate choices;
geodetic choices require the configured body shape. Cartesian displacement is
independent of the station's input coordinate representation.

The complete pinned RP03 provider grid was parsed independently. Its rounded
longitude-seam and polar coefficients do not satisfy exact consistency. Strict
evaluation rejects these inconsistencies; an explicit averaging approximation
reports the discrepancy and preserves the raw data. Source-origin preservation
is the default. Optional CE-to-CM translation uses IERS Table 7.6 only for an
explicitly identified RP03/no-IB/Farrell source and applies once in Cartesian
coordinates. An already-CM source receives no extra translation.

The six input, eight synthesis/interpolation and seven model/factory C++ cases
pass, including all lookup/basis combinations, oblate geometry, file/preloaded
parity, exact station identifiers and the 2016 inserted second. Independent
complex-Fourier fixtures check the published origin-translation equations;
this fixture generator is not presented as original IERS reference software.
A mismatching Bernese station example is recorded as an unsuitable oracle,
without widening a tolerance. Source provenance and these limitations are in
``reference_files/atmospheric_tides/README.md``.

The optimized integration build exposed a dangling Eigen expression in one
test's expected-value calculation. Materializing the vector fixed the test
without changing model equations, fixtures or tolerances. All 24 selected C++
executables and the actual input-target dependency audit now pass (25 CTest
checks, 1.98 s), together with both source-boundary checks. Exact-file formatting
also passes. The final source identities and evidence are recorded in
``time-routing-atmospheric-final-files.json`` and
``time-routing-atmospheric-validation-summary.json``. The full atmospheric
velocity/state interfaces explicitly reject evaluation until the selected UT1
rate is connected; task 17 remains open. No Python kernel was rebuilt.

Task 20, transformation checkpoint: the neutral IGN table reader preserves
published realization/release identities, decimal-year parameter epochs,
original records/units and source provenance. The standalone input target links
only Eigen. Explicit catalogue selection retains the original ITRF2014 and
superseding ITRF2020 tables without inventing datum changes for aligned releases.

Requested-epoch affine position and state evaluation use the published linearized
Helmert matrix, its exact inverse and complete derivative/composition terms.
The structured epoch supplies dDecimalYear/dInputSecond explicitly. Compatibility
functions retain their uniform-Julian-year seconds since J2000 and ITRF2014
catalogue, with corrected table transcription, same-frame identity and moving
states. Their previous fixed-reference-epoch propagation is removed.

All 25 published rows are checked at the reference epoch and +/-10.25 years:
75 independent 80-digit Decimal forward/inverse moving-state fixtures agree
within 2 nm and 2e-17 m/input-second. The generator verifies original and embedded
table hashes and reproduces the frozen fixtures. It implements the published
equations independently; it is not original IGN/IERS reference software.
Five model and three input C++ cases also check individual parameter/rate signs,
nonunit clock rates, inverse conditioning, identity, provenance and malformed
data. Existing reference-frame compatibility tests pass.

The recovered integration source matched all 38 files in its frozen manifest.
All 27 selected C++ executables and the actual input-target dependency audit pass
(28 CTest checks), together with both source-boundary checks. Exact-file
formatting and fixture regeneration pass. Binding declarations and discoverable
data/model API references accompany the implementation; no Python kernel was
rebuilt. Evidence is recorded in ``recovery-itrf-validation-summary.json`` and
the associated ``recovery-itrf-*``/``recovery-test-itrf-*`` logs in the local work
directory. Task 20 remains open for station-catalogue integration with tasks 18
and 19 and its chosen physical-time mapping; task 29 owns compiled binding checks.

Task 17, rate completion after session recovery: native atmospheric deformation
now differentiates the selected continuous UT1 model and applies dUT1/dTAI and
dTT/dTDB exactly once. Positions remain available where a selected slope is
unavailable. Model-owned time conversion preserves coverage endpoint policy;
origin translations are differentiated in Cartesian coordinates. The nine
synthesis and ten model cases pass in the actual integration CMake build,
including knots, endpoints, leap seconds, geometry and factory equivalence.
The three UT1 primitive cases and raw target-link audit also pass. Independent
real128 origin-rate fixtures regenerate identically. Evidence is in
``recovery-atmospheric-main-validation-summary.json``; Python-only runtime
binding validation is queued under task 29.
