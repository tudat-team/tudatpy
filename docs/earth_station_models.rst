Earth rotation and station modelling implementation
==================================================

Starting point and decisions
----------------------------

This work starts at PR #900 (``feature/data-refactor``), commit
``3070b20e83e18a749767fe0332a5280464f59b16``.  Changes are developed on
``fix/earth_deformation_rotation`` and committed as validated units. The session
was recovered on 14 October 2026 after an interrupted context compaction.
Dominic explicitly authorized pushing the current milestone and subsequent major
milestones on 14 October 2026. The coordinator publishes reviewed integration
checkpoints to the existing ``origin/fix/earth_deformation_rotation`` branch.

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

Status: all 25 retained implementation tasks and the initial full-PR review are
complete. Dominic requested a further personal static review of complexity and
exposure on 14 October 2026; its findings and outstanding cleanup are recorded
below. Completion of the physical models does not close those review findings.

.. code-block:: text

   [x] 01 Establish the PR #900 baseline, isolated build and task tracker.  [DONE]
   [x] 02 Define neutral data contracts and enforce the dependency boundary.  [DONE]
   [x] 03 Separate existing EOP and BLQ readers from model evaluation.  [DONE]
   [x] 04 Establish pinned independent IERS/SOFA reference fixtures.  [DONE]
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
   [x] 19 Implement post-seismic and seasonal station motion with rates.  [DONE]
   [x] 20 Modernise ITRF transformations and station-catalogue inputs.  [DONE]
   [x] 21 Add non-tidal atmosphere/ocean/hydrology loading inputs.  [DONE]
   [x] 22 Integrate non-tidal loading with explicit overlap/validity policies.  [DONE]
   [x] 23 Validate relativistic station position and velocity transformations.  [DONE — explicit IERS spatial approximation and full intrinsic-state composition]
   [x] 24 Improve instrument reference-point and calibrated geometry support.  [DONE]
   [x] 29 Consolidate public interfaces, diagnostics and reproducible examples.  [DONE — simplified API and actual Python runtime validation]

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

All 25 retained tasks and the initial full-PR review are complete. Dominic's
comments on `draft PR 1092 <https://github.com/tudat-team/tudatpy/pull/1092>`_
are now being addressed.
Implementation workers use **GPT-6.1 Sol, high reasoning** in isolated worktrees.
The coordinator owns design, literature analysis, review and integration.
The three workers' documentation/formatting work is paused and preserved in
their separate worktrees. The coordinator personally reviewed the current PR
for unnecessary layers, duplicate code paths, numerical risks and exposure
placement. Gravity work remains paused at its checkpoints. Dominic's latest
instruction prohibits builds and test/application execution; no compilation,
tests or numerical experiments were run for this review. The external-review
handover is held until the additional cleanup is resolved.

.. code-block:: text

   [x] Remove the EOP adapter that derives from a mathematical interpolator.  [DONE — pushed e2be54597]
   [x] Remove the additional CMake dependency-check file and its registration.  [DONE]
   [x] Simplify duplicate polar-motion overloads and the short-period rate helper.  [DONE]
   [ ] Rewrite the reviewed comments and explain the remaining EOP/clock responsibilities.  [PAUSED — structural review first]
   [x] Remove the type aliases introduced by this PR.  [DONE]
   [ ] Add braces on separate lines around single-statement control flow introduced by the PR.  [PAUSED — edits preserved]
   [ ] Describe functions in headers and implementations, variable groups and computation steps.  [PAUSED — edits preserved]
   [ ] Add blank lines after variable definitions and function declarations.  [PAUSED — edits preserved]
   [x] Validate the affected native EOP and rotation paths; review and push the fixes.  [DONE]
   [x] Personally review complexity, reuse, numerical contracts and C++/Python exposure.  [DONE — static analysis of 5a69997a6]
   [ ] Resolve the structural findings below before the next review handover.  [OPEN — recommendations recorded, implementation not started]
   [ ] Push the comment/style checkpoint and create MODIFICATIONS-DONE.md for local review.  [HELD — no marker created]
   [ ] Check every ten minutes for REVIEW-DONE.md and assess its recommendations.  [QUEUED — follows handover]
   [ ] Implement accepted external-review changes and push the PR branch.  [QUEUED — execution remains prohibited until instructed otherwise]

The detailed personal review is saved locally at
``/home/dominic/Tudat/earth-station-model-work/pr1092-personal-static-review.rst``.
It identifies source locations, specific removals/consolidations, numerical
constraints, exposure placement, code to retain, and test reductions. These are
review recommendations, not claims of implemented or validated refactoring.

.. code-block:: text

   [ ] R01 Unify time-conversion routing while preserving custom clocks and precise Time.  [RECOMMENDED — numerical refactor]
   [ ] R02 Centralize EOP endpoint handling and correct the inconsistent scalar CPO route.  [RECOMMENDED — static correctness concern]
   [ ] R03 Remove duplicate daily-EOP record ownership and the factory's legacy-reader detour.  [RECOMMENDED]
   [ ] R04 Make orientation wrappers delegate to one precise calculation.  [RECOMMENDED]
   [ ] R05 Remove copied histories from transformation/reference-point wrappers.  [RECOMMENDED]
   [ ] R06 Remove the frame epoch wrapper and redundant clock-forwarding API.  [RECOMMENDED]
   [ ] R07 Move runtime bindings into existing frame_conversion/environment modules.  [RECOMMENDED]
   [ ] R08 Simplify Python construction and keep internal helper records out of normal workflows.  [RECOMMENDED]
   [ ] R09 Normalize file input into data-based settings at construction.  [RECOMMENDED — confirm file-loading timing]
   [ ] R10 Consolidate station adapters and repeated value/rate loading work.  [RECOMMENDED]
   [ ] R11 Normalize relativistic callbacks and share guarded sampling locally.  [RECOMMENDED — preserve fallback guarantees]
   [ ] R12 Distinguish unavailable deformation data from failed calculations.  [RECOMMENDED — inherited broad catch remains a concern]
   [ ] R13 Reduce loading-specific coupling and repeated compatibility work.  [RECOMMENDED — preserve active-model checks]
   [ ] R14 Reuse existing calendar/basis utilities and hide model-specific helpers.  [RECOMMENDED]
   [ ] R15 Simplify neutral-reader target registration.  [CONDITIONAL — preserve dependency/install contracts]

The review retains independent reference fixtures, precise continuous clocks,
selected-interpolant derivatives, HARDISP-specific spline behavior, explicit
frame/origin identities and separate position/state availability. The earlier
handover estimate no longer applies. Worker edits have not been discarded or
integrated during this personal review.

The EOP review preserves continuous TAI sample times, UT1-TAI interpolation,
rates of the selected interpolating function, and the existing coverage and
missing-data checks. The standard factories already use this data directly;
they do not need an additional interpolator adapter. Existing constructors
accepting user-supplied interpolators remain supported.

The first follow-up removes the CMake check and clarifies the existing custom
interpolator path, station position/state methods, derivative fallback, fitted
seasonal terms and relativistic coordinate conventions. These are comments and
build-check removal; they do not change numerical evaluation. The adapter removal
is a separate change and is not included in that first checkpoint.

The six added type aliases are removed at ``206ae7d45``. Their uses now spell
out the existing underlying types; affected native sources and tests compile.
Numerical tests and the Python kernel were not rerun for this spelling change.

``e2be54597`` removes ``EarthOrientationDataInterpolator`` and its factory
allocations. Models created from EOP records use the shared evaluator directly;
the old interpolator getters remain available for constructors supplied with
actual interpolators. The public clock-update method also uses the continuous
conversion functions, avoiding a null interpolator after this removal. Focused
checks pass: direct clocks 3 cases/110 assertions, reader/leap behavior 2/128,
factory/coverage 5/95 and polar motion 1/7. Original numerical tolerances remain.

``658bef7fa`` replaces the duplicate TT polar-motion methods with one template
and writes the short-period rate calculation as an explicit four-sample loop.
The paired-difference arithmetic and precise-time handling are retained. The
existing short-period rate and precise TT-boundary checks pass 2 cases/24 assertions.

The requested ``MODIFICATIONS-DONE.md`` review trigger will be written in the
monorepo parent after the documentation/formatting checkpoint is pushed. The
source ownership manifest is recorded locally in
``earth-station-model-work/logs/pr1092-style-ownership.json``.
Another local AI will write ``REVIEW-DONE.md`` beside the trigger. After writing
the trigger, check for the review every ten minutes. Assess its suggestions
seriously against Tudat's existing code and numerical requirements; do not copy
them without review. Sol agents implement accepted changes and run all relevant
tests. Fix failures before the final reviewed checkpoint is committed and pushed.

The station-history/parser/frame/reference-point, atmosphere/ocean/build,
rotation/core and cross-component reviews are complete. Their fixes are
integrated through source checkpoint ``e2335bee7``. The
independent tests listed below remain valid when their tested implementation and
dependencies survive integration unchanged; integration alone is not a reason
for another test run. The actual Python runtime validation is complete.

The completed-task milestone was pushed at ``474ef5bf3``. This reviewed checkpoint
and its final ledger are published to ``origin/fix/earth_deformation_rotation``
under Dominic's existing authorization.

Required full-PR review
-----------------------

.. code-block:: text

   [x] Review the complete station/rotation PR after all 25 tasks.  [DONE]
   [x] Review history/parser/frame/reference-point implementation and Python ownership.  [DONE]
   [x] Review rotation, time conversion, deformation core and relativistic composition.  [DONE]
   [x] Review atmospheric/ocean models, high-frequency EOP and build/data repetition.  [DONE]
   [x] Complete cross-component review and final evidence.  [DONE]

This review covers the entire diff from ``3070b20e83e18a749767fe0332a5280464f59b16``.
Beyond numerical correctness and ordinary regression review, it checks reuse of
existing functionality, proportional changes to generic code, unnecessary
abstractions and boilerplate, consistency with neighboring Python APIs, and
excessive tests. It seeks fewer added lines without compromising numerical quality
or significantly increasing runtime. Independent scientific fixtures are retained;
a compact loop with many numerical assertions is not itself excessive testing.

Integrated review changes:

* The SINEX state reader now delegates to the strict history parser. The duplicate
  headerless parser and its compatibility-only test were removed, as backwards
  compatibility was waived. Real SINEX fixtures and ILRS consumers remain covered.
  This removes 291 net lines; both affected native suites pass (13 cases/130
  assertions combined).
* Terrestrial frame transformations reuse the existing cross-product matrix.
  Atmospheric loading reuses the existing ENU quaternion after an explicit R/N/E
  to E/N/U component permutation. The respective existing suites pass 5/995 and
  10/351, with unchanged fixtures and tolerances. These remove eight net lines.
* Relativistic station rate fallback now rejects falsely agreeing tiny stencils
  when a callback internally narrows Time. It compares against a representable
  scalar-resolution stencil without narrowing the actual Time samples. The
  factory retains the selected converter and extended-time rotation derivative.
  Regressions fail before the fix and pass afterward; the affected suite is 7/103.
  The additional stencil also rejects nonfinite derivatives.
* A precise historical time-conversion cutoff defect is fixed with typed bounds.
  Five reproduced failures (including a roughly 29.858 ms wrong-branch jump) now
  pass in the affected 1-case/11-assertion transition regression.
* Six neutral-reader CMake blocks now use a local helper, removing 91 lines.
  Static/shared and imported/fallback Eigen contracts match, generated exports
  are byte-identical, and the dependency checks pass. No global build framework
  changed.
* Loading fixture provenance now uses exact source ranges, explicit count edits
  and whole-source/extract hashes, removing 874 net lines. Offline regeneration
  reproduces both fixture bytes and hashes. The provenance and dependency checks
  pass. Atmospheric latitude lookup directly uses the existing configured-shape
  coordinate converter. The unused new EOP-helper interpolation argument was
  removed; affected atmospheric and EOP tests pass 10/351 and 4/244.
* Rotation derivatives reuse the existing cross-product matrix, and nominal CIP
  interpolation nodes are copied once. Two duplicated routing fixture cases
  were removed; independent IERS/ocean reference coverage remains. The affected
  rotation and routing checks pass 4/169 and 9/65, removing 69 net lines.
* The common displacement derivative now checks sub-scalar-resolution Time
  stencils against a representable four-ULP stencil. The scalar callback
  counterexample fails before the fix; the precise-Time control and all 12/72
  affected checks pass afterward. Both intervals are validated and central
  caches restored. Normal-step evaluation cost is unchanged.

The integrated review removes more than 1,500 net lines relative to the completed
task checkpoint. The coordinator also checked the generic-code changes,
construction/registration dependencies, scientific conventions and Python API
consistency. No unresolved station/rotation correctness finding remains.

Completed behavior
------------------

Earth orientation and clocks
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Neutral readers preserve legacy/current C04 and finals2000A source conventions,
missing values, per-component status and provenance. Production evaluation uses
continuous TAI and UT1-TAI, strict raw coverage and a default 36-hour gap limit.
Linear interpolation is the default; explicit four-point Lagrange interpolation
validates its selected stencil and UT1 monotonicity. There is no hidden endpoint
padding, extrapolation, order reduction or missing-CPO fallback.

Rotation, polar motion and time conversion share the selected EOP model. TT/UT1
arguments and their rates are routed explicitly. J2000/ICRF share the GCRS axis
orientation; EME2000 applies its explicit bias and ECLIPJ2000 its fixed rotation.
This does not equate origins or relativistic coordinate scaling. The optional
Desai--Sibois 159-term model replaces the corresponding ocean EOP corrections;
conventional libration is retained once. The conventional default is unchanged.

Deformation and complete station states
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``EarthDeformationEpoch`` is removed. The common base exposes only displacement
and complete displacement state, taking TDB ``Time`` plus a station. Derived
models implement protected hooks; converters are owned/shared by the models.
Output is Cartesian SI, with state velocity differentiated with respect to TDB.
Station, motion and reference-point paths preserve native Time. Scalar endpoint
normalization is confined to real scalar API boundaries; precise Time is unsnapped.

Solid Earth tides retain their corrected production phase conventions. Ocean
loading defaults to the 342-constituent HARDISP synthesis from the eleven BLQ
constituents. The older eleven-constituent synthesis remains an explicit option.
Solid and ocean pole models share the 2018 secular pole. Atmospheric S1/S2 loading
uses continuous UT1, separate lookup-latitude and displacement-basis conventions,
and analytic rates. Loading origins and provider-specific translations are explicit.
Non-tidal loading has strict station/coverage selection and component-overlap checks.

Position-only access remains usable when a derivative is unavailable. Generic
finite differences validate representable uniform stencils, preserve failures
and restore central caches. Relativistic fallback also checks coarse/fine
convergence. The common displacement fallback checks scalar resolution when a
native Time step is tiny. Optional missing central data may be
skipped; invalid configuration and nonfinite numerical results are not suppressed.
Complete station states compose intrinsic motion, deformations, rotation and body
translation. Terminal relativistic composition consumes that complete state once.

Station histories, frames and instrument geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The strict SINEX reader retains solution validity, discontinuities, equipment,
eccentricities and station/point/solution identity. Motion supports cumulative
post-seismic and seasonal terms with physical TAI validity and an explicit IGN
nominal-UTC source clock. Julian-year velocity units and calendar seasonal phases
remain distinct. ITRF catalogue selection retains release/origin/tide conventions;
complete affine inverses and state derivatives reuse one transformation model.

Dated marker/ARP/calibrated reference points reuse retained equipment histories and
fixed local-axis conventions. The Python record getters return copies so callers
cannot mutate validated model records through a property. Public factories accept
preloaded neutral data and follow the existing settings/module conventions.

Validation evidence
-------------------

Normal tests are offline. Pinned fixture directories record sources, revisions,
hashes, licences, exact extraction/generation steps and the numerical budget.
Reference programs executed for fixtures are distinguished from independent
implementations of published equations. Logs and source manifests are retained in
``/home/dominic/Tudat/earth-station-model-work/logs`` and the isolated build trees.
The previous chronological ledger remains in Git history; the evidence below
states the final behavior instead of retaining superseded open-task claims.

* Solid Earth displacement: the unmodified corrected IERS DEHANTTIDEINEL package,
  its two published examples and 18 frozen cases; 1 nm component tolerance through
  the production factory. Final affected suite: 2 cases/153 assertions.
* Solid pole displacement: eight independent colatitude-equation vectors. The
  audited largest reference difference is 8.67e-19 m within a 1e-14 m budget.
  Final suite: 3/135.
* Short-period EOP: pinned PMSDNUT2, UTLIBR, ORTHO_EOP/CNMTX and the Paris 71-wave
  routine at 35 epochs. The published T2 sine coefficient (0.05 microseconds)
  differs from UTLIBR (0.04); that documented distinction is retained explicitly.
* EOP interpolation: 27 frozen LAGINT/SOFA cases, including the 2016 leap second
  and 1962 drift era, plus analytic interpolation/rate and inverse checks.
  The original LAGINT undefined final-node branch is not used as an oracle.
* Rotation/station closure: 13 original rotation cases/334 assertions; full
  SOFA station-velocity comparison error 7.42e-10 m/s. Complete station suite
  5/151 covers four base axes, five clocks, Time/scalar inputs and body translation;
  maximum closure error 4.53e-9 m/s against a 1e-7 m/s budget.
* Desai--Sibois: reader 3/32 and model 4/244, fourteen independent value/rate
  fixtures over 1962--2025, full factory/SOFA sequence, retained libration and
  ambiguous-overlap rejection. The embedded source matches its pinned NASA mirror.
* Ocean loading: 36 original IERS HARDISP cases, maximum difference 0.106
  micrometres within its documented 1 micrometre budget; existing 24-hour series
  also pass at 1 micrometre. Final suite: 6/164.
* Ocean pole loading: 731 unchanged published rows within 1 nm (maximum 0.7894 nm).
  Their historical constants/pole are used only for that comparison. Production
  uses IERS 2010 constants and the shared 2018 secular pole. Model suite 7/2287;
  reader suite four cases; full provider grid parsed independently.
* Atmospheric loading: independent Fourier/origin/rate fixtures and complete RP03
  input inspection; current model suite 10/351. Source imperfections and explicit
  approximation policies are tested without fitting tolerances to provider errors.
* SINEX/history: eight complete pinned IGN products and an intentionally rejected
  malformed official IGS record; motion 9/343 includes 126 independent original
  IGN Fortran rows. The strict-reader/ILRS cleanup passes 6/68 and 7/62.
* Terrestrial frames: 75 independent 80-digit Decimal forward/inverse state
  fixtures from all 25 published rows; agreement within 2 nm and 2e-17 m per input
  second. Affine suite 5/995; catalogue adapter 4/217 against 13 published rows.
* Non-tidal loading: reader 4/141 plus four dependency/provenance checks; model
  7/80. A provider file with eight duplicate station labels is deliberately
  rejected; positive and negative extracts retain the original evidence.
* Relativistic stations: independent 70-digit equation fixtures, an exact Lorentz
  limiting case, complete factory/origin conversion and analytic-rate controls.
  The current affected suite passes 7/103, including reproduced precision failures.
* Instrument geometry: 5/52 plus relevant history and frame-adapter cases.
* Core deformation/state: velocity 12/72, station aggregation 3/331 and clock
  routing 9/65 after duplicate reference cases were removed. Every changed native core
  translation unit compiled, including the actual loading-composition validator.
* Python-specific validation: a fresh coherent CMake kernel was built and its
  actual loaded path asserted. Seven cases passed initially; two exposed a mixed
  Time/scalar endpoint comparison and missing ocean-pole settings registration.
  Both affected cases passed after the fixes. The unchanged seven are retained;
  there is no claim of an unnecessary repeated full run.
* Architecture: transitive Python imports, C++ includes and actual CMake link
  closure enforce the neutral input boundary, including private static links.
  Binding registration orders neutral data before dynamics/estimation. The
  isolated stale EOP-factory ABI objects were rebuilt and the failing case passed.

Material limits and explicit policies
------------------------------------

* Model-equation agreement is not a claim of equivalent real-world physical
  accuracy. Declared reference constants, phases, clocks and origins must match.
* Scalar UTC cannot represent an inserted second. TAI/TT/UT1 remain continuous.
  UTC-phase derivative models reject intervals crossing an offset step or its
  ambiguous inserted second; smooth historical UTC drift is supported.
* EOP knots/coverage can provide positions without usable rates. Selected clock
  interpolators lacking a rate contract reject rate evaluation. Historical
  Delta T is an explicit clock-only policy, not invented historical EOP or UTC.
* The ocean-pole grid ends at +/-89.75 degrees. Latitude clamping is an explicit
  approximation. The CM grid already includes its geocentre contribution.
* RP03 longitude-seam/polar rounding defects reject by default. Explicit averaging
  reports the discrepancy. CE-to-CM translation is allowed only for the identified
  RP03/no-IB/Farrell profile; it is not applied to already-CM coefficients.
* Relativistic station conversion implements the documented IERS spatial
  approximation. It does not claim full event/simultaneity/worldline dynamics.
  Numerical differentiation requires smooth callbacks and sufficient coverage;
  convergence cannot prove arbitrary black-box smoothness.
* No new estimation parameters, complete cross-software trajectory campaign,
  observation omission study or deferred research-extension task is retained.

Parallel gravity-field programme
--------------------------------

Dominic explicitly authorized gravity work in parallel on a separate branch.
The initial station review is complete; the additional personal review and
cleanup status is recorded above. ``fix/gravity_field_variations`` and
``tudatpy-earth-gravity-variations`` contain the integration branch/worktree;
``docs/earth_gravity_field_variations.rst`` is its separate checklist, also open
in gedit. Primary-literature review and the existing-model inventory are complete.
Compatibility design continues; the static ICGEM/tide-system foundation is
integrated and validated. Gravity work is currently paused, with separate
worktrees and checkpoints preserved. Implementations
must share compatible rotation/EOP, clocks, frame axes, static tide conventions
and physical component inventories. Desai--Sibois ocean EOP corrections do not
replace ocean gravity tides or either pole tide. Gravity estimation is excluded.
