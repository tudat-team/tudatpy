Earth rotation and station modelling implementation
==================================================

Starting point and decisions
----------------------------

This work starts at PR #900 (``feature/data-refactor``), commit
``3070b20e83e18a749767fe0332a5280464f59b16``.  Changes are developed on
``feature/earth-station-models`` and committed locally as validated units.

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
  Initially the Tudat build owner has six and the reference worker has two.
  No other worker may compile until the coordinator reallocates that budget.
  Disable hidden nested build parallelism and multithreaded numerical runtimes.
* Every change includes relevant bindings, documentation and tests.  Tolerances
  come from reference accuracy and numerical error budgets, not from fitting
  the implementation's current output.
* Unit tests primarily use C++.  Python tests are reserved for Python-specific
  imports/bindings or cases that cannot straightforwardly be exercised in C++.
  Reference-data generation scripts may use Python to drive external programs.

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

.. code-block:: text

   [x] 01 Establish the PR #900 baseline, isolated build and task tracker.
   [ ] 02 Define neutral data contracts and enforce the dependency boundary.
   [ ] 03 Separate existing EOP and BLQ readers from model evaluation.
   [ ] 04 Establish pinned independent IERS/SOFA reference fixtures.
   [x] 05 Correct solid pole tides and validate spherical/ellipsoidal geometry.
   [ ] 06 Validate solid Earth tides through the production time/argument path.
   [ ] 07 Resolve inertial frames, time scales and displacement conventions.
   [ ] 08 Implement robust legacy/current/rapid EOP ingestion with metadata.
   [ ] 09 Implement EOP validity and leap-safe interpolation policies.
   [ ] 10 Verify existing short-period EOP corrections independently.
   [ ] 11 Complete Earth-rotation derivatives.
   [ ] 12 Complete deformation velocities and full station-state closure.
   [ ] 13 Add the optional Desai-Sibois high-frequency EOP configuration.
   [ ] 14 Validate BLQ loading inputs and preserve their conventions.
   [ ] 15 Implement HARDISP-equivalent ocean tidal loading.
   [ ] 16 Add conventional ocean pole-tide loading.
   [ ] 17 Add conventional atmospheric S1/S2 tidal loading.
   [ ] 18 Preserve full SINEX solution/validity/eccentricity histories.
   [ ] 19 Implement post-seismic and seasonal station motion with rates.
   [ ] 20 Modernise ITRF transformations and station-catalogue inputs.
   [ ] 21 Add non-tidal atmosphere/ocean/hydrology loading inputs.
   [ ] 22 Integrate non-tidal loading with explicit overlap/validity policies.
   [ ] 23 Validate relativistic station position and velocity transformations.
   [ ] 24 Improve instrument reference-point and calibrated geometry support.
   [ ] 25 Add Earth-orientation estimation parameters and verified partials.
   [ ] 26 Add station-motion/loading parameters and verified partials.
   [ ] 27 Compare complete trajectories with independent geodetic software.
   [ ] 28 Validate observations, controlled omissions and synthetic recovery.
   [ ] 29 Consolidate public interfaces, diagnostics and reproducible examples.
   [ ] 30 Evaluate research extensions after the accepted-standard baseline.

Task numbers are stable identifiers rather than a rigid execution sequence.
Tasks 05 and 06 can proceed with the initial build and data architecture audit.
Tasks 13 and 30 follow verification and implementation of the accepted-standard
baseline; experimental alternatives do not become defaults implicitly.

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
includes/associated implementations.  The committed checks use only Python's
standard library and can run without a Tudat build.  Data contracts, independent
input-library linkage and compiled binding initialization remain to validate.

Task 05: corrected the pole-tide wobble units and east/north formulas, and use
the geocentric basis for the IERS displacement before projecting to a station's
topocentric frame.  The C++ PoleTideDisplacement executable passes eight frozen
reference cases evaluated independently from the original colatitude equations,
zero-wobble checks, and Cartesian/geodetic/spherical representation checks on
an oblate Earth.  The conventional secular pole remains the explicitly pinned
1 February 2018 IERS update.  Compilation also exposed and fixed a missing
direct constants include in the geodetic-coordinate header.

Solid-tide and inertial-frame corrections are under integration test.
