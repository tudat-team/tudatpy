Earth gravity-field variations: design and implementation tracker
================================================================

Scope and status
----------------

Branch: ``fix/gravity_field_variations``, created from the Earth rotation/station
checkpoint ``a3be1612a`` on 14 October 2026. Dominic authorized this programme to
start in parallel with the remaining station validation and review. The two
programmes have separate integration worktrees. The reviewed station fixes have
now been merged before extending gravity setup.

Status updated 15 October 2026: resumed after the rotation/deformation handover.
The static ICGEM adapter and solid-Earth, solid-pole and ocean-pole models are
integrated and have passing reference checks. Their conventional Earth factories
are also implemented locally; the earlier factory checkpoint passed 3 cases and
16 assertions. The tracker below supersedes the older "factory pending" status.
The harmonic ocean checkpoint is integrated. Its recovered original reader and
model tests pass 7 cases and 57 assertions; the largest of four COST-G M2 force
residuals is 1.57e-13 m/s^2 against the unchanged 1e-11 m/s^2 budget.

Authoritative checkout: ``tudatpy-earth-gravity-variations`` on
``fix/gravity_field_variations``. Last published gravity
checkpoint: ``98b2e67e1``. Merge ``2a79c232d`` includes reviewed rotation/deformation
revision ``d34b32718``; ``83f36358c`` integrates harmonic checkpoint ``353aa728b``.
Checkpoint ``0ad9cf79d`` adds harmonic consuming settings and removes the two
empty conventional preset classes. The later local fixes declare existing
Earth-orientation link dependencies and compare harmonic matrix dimensions
before coefficient equality.
Checkpoints ``39f58f9b6`` and ``d42704f96`` add AOD1B RL07 input, interpolation
and Earth-factory integration, including component overlap checks and serialization.
TiME22 source is integrated at ``2ae468a64`` and the first Python bindings at
``c2700aa53``. ``5e982c022`` moves final harmonic inventories into settings to avoid
an unnecessary large temporary copy.
``6a6fc5082`` documents the harmonic code and skips exact-zero admittance weights.
``b39c320fd`` integrates neutral temporal ICGEM records and the explicit static
consumer rejection. A separate comments-only pass covers the conventional and
ocean-pole models without changing numerical expressions.
``72ed3b5fc`` removes the raw ICGEM record include from the generic gravity setup
header; only code that consumes those records now includes their complete type.
The temporal evaluator and consuming setup are integrated at ``9202b8481``.
At Dominic's request, merge ``53284e544`` brings in the latest published
rotation/deformation revision ``d74af05bc``. Its shared text parser is reused by
the ICGEM and AOD1B readers in ``4eb777328``, removing forty net lines without
changing coefficient interpretation. The complete planned Python settings are
also integrated, including the temporal factory and explicit time conventions.
Gravity has not been pushed during this resumed work.
Local rotation/deformation commit ``3613e38d3`` adds its C++ fixtures, support
files, generators and seven resources to ``tests/test_tudat/data.zip``; the
existing 705 archive entries are unchanged. All 121 added source files passed
fresh-extraction byte comparisons. This separately requested change is on the
local rotation branch and merged into gravity at ``5b2be31cd``; it was not pushed.

The coordinator designs, reviews and handles major debugging. Three GPT-6.1 Sol
workers at high reasoning use separate worktrees: rotation_impl owns
dependency integration, the build and a separate Python exposure checkout;
desai_impl has completed the independent temporal numerical tests and evaluator
review and now owns Dominic's separate rotation/deformation ``data.zip`` task;
history_impl has completed the temporal evaluator, consuming factory and parser
reuse. rotation_impl has validated the initial targets, harmonic
setup, AOD1B and TiME22, and owns the remaining native and Python validation in
``/home/dominic/Tudat/gravity-resume-build``
from a consistent source tree and owns the entire eight-job allowance. Its
configuration follows the validated PR1092 build, including serialization.
Tests run serially. Only relevant native tests run; the Python kernel is built
solely for Python-specific validation. No new estimation parameters, gravity
recovery work or examples are included.

.. code-block:: text

   [x] G01 Review primary literature, provider conventions and existing implementations.  [DONE]
   [x] G02 Define rotation/gravity compatibility and reuse existing Tudat functionality.  [DONE — selected clocks, frames, tide systems, component overlap and temporal baseline conventions are explicit]
   [x] G03 Preserve neutral coefficient, constituent and source-convention records.  [DONE — static/temporal ICGEM, harmonic, ocean-pole and AOD1B readers pass]
   [x] G04 Complete conventional solid Earth gravity tides, including frequency dependence.  [DONE — model and factory checkpoint validated; dependency integration is G11]
   [x] G05 Add solid Earth pole-tide gravity corrections with the shared pole convention.  [DONE — model and factory checkpoint validated]
   [x] G06 Add ocean pole-tide gravity corrections with explicit truncation and origin.  [DONE — full model and factory checkpoint validated]
   [x] G07 Add harmonic ocean-tide synthesis and documented minor-constituent admittance.  [DONE — harmonic setup, TiME22 published admittance and exact-zero optimization pass]
   [x] G08 Support identified ocean-tide products without mixing their phase conventions.  [DONE — FES2004/FES2014b/EOT11a profiles and TiME22 1.1 pass focused native checks]
   [x] G09 Handle atmospheric/non-tidal coefficient products and component overlap explicitly.  [DONE — AOD1B reader/model/factory and TiME22 atmosphere pass]
   [x] G10 Connect secular/periodic/tabulated variations to their static-field baseline.  [DONE — temporal evaluator passes 3 cases/63 assertions against independent scalar references]
   [x] G11 Integrate Earth factories with the selected rotation and enforce compatibility.  [DONE — final Earth assembly/serialization passes 8 cases/57 assertions]
   [x] G12 Validate individual models and assembled forces against pinned references.  [DONE — 14 native executables, 68 cases and 3069 assertions pass, including the two narrowed links]
   [ ] G13 Complete Python settings and diagnostics.  [ACTIVE — all planned bindings integrated and reviewed; kernel/runtime/autosummary validation pending; examples excluded by Dominic]
   [ ] G14 Review the complete gravity diff for correctness, reuse, API clarity and size.  [ACTIVE — source/API review complete; verification-driven follow-up remains]

Resumed verification
--------------------

The first build of the isolated solid-tide test exposed missing link declarations
for the Earth-orientation library's existing SOFA and raw-reader dependencies.
The two-line owning-target fix links successfully; no numerical source change
was required. Logs retain the failed attempt and its successful replacement in
``/home/dominic/Tudat/gravity-resume-build/verification-20261015``.

The initial ten native executables pass: ICGEM, ocean-pole and harmonic-potential
readers; solid-Earth, ocean-pole and harmonic gravity models; ICGEM settings;
Earth variation setup; existing gravity variations; and environment model setup.
Together they pass 52 cases and 2452 assertions, with zero failures. This includes
the Earth-relative Moon and deferred-interpolation factory checks. Evidence is in
``verification-20261015/native-summary-initial.json``. Runtime tests are serial
with numerical-library threads set to one.

The later harmonic checkpoint also passes: harmonic reader 3 cases/30 assertions
and Earth setup 6 cases/39 assertions, including ordered admittance, selected
clocks, deferred sampling and exact serialization. AOD1B and TiME22 share the next
focused four-target build, which now passes: AOD1B reader 2 cases/49 assertions,
AOD1B model 2/12, updated Earth factory 7/49 and TiME22 setup 2/430. The build and
serial runs succeeded without numerical fixes or tolerance changes. Evidence is
in ``verification-20261015/aod-time22-build.log`` and
``verification-20261015/aod-time22-runtime.log`` with per-case reports.
The retained results through source ``2621e1b4e`` cover thirteen distinct native
executables, 62 cases and 2978 assertions, all passing. This total counts the
latest result for each executable, rather than adding repeated runs.

The subsequent exact-zero admittance optimization, temporal reader and include
cleanup pass their focused integration checks. Python source and documentation
now cover AOD1B, TiME22 and raw temporal records. Review fixed Python getters that
could otherwise expose mutable owned harmonic metadata/constituents: these return
copies, with assertions added to the existing provider-settings test. Kernel and
Python runtime validation remain pending.
The updated temporal ICGEM reader now also passes 6 cases and 50 assertions;
``verification-20261015/icgem-temporal-reader-runtime.log`` records that run.
That four-target build has now finished successfully at ``86c71fb1f``: ICGEM
reader 6 cases/50 assertions, harmonic model 4/29, static ICGEM adapter 2/28 and
TiME22 setup 2/430, all passing. The temporal evaluator has separate, independently generated scalar
references covering published multi-interval periodic terms, exact shared
boundaries, normalization and the two explicit definitions of a year. Those
model checks now pass: 3 cases and 63 assertions. Final Earth setup passes
8 cases/57 assertions, including the selected converter, adopted nominal
baseline, shared interval start and settings serialization. The static ICGEM
adapter also passes 2 cases/28 assertions. No numerical fixes or tolerance
changes were needed. The build log is ``verification-20261015/g10-final-build.log``.
The compiler checkout has integrated the reviewed source, including the updated
rotation dependency and shared parser reuse. One Python kernel build now precedes focused Python/architecture tests,
autosummary-name validation and the full Python suite against the established
offline-resource baseline. No broad native test build is planned.
At this final snapshot, the ICGEM reader (6 cases/50 assertions) and AOD1B reader
(2/49) already pass after shared parser reuse. The autosummary inventory covers
283 entries in 21 changed RST pages; runtime lookup still awaits the kernel.

Source review has checked the generic changes, baseline subtraction, selected
clocks, independent temporal origins, normalization and model overlap. The final
code reuses the existing gravity update, interpolation and derivative machinery;
it adds no epoch wrapper or universal configuration layer. Review removed two
empty settings classes, duplicated parsing, an unnecessary harmonic inventory
copy and mutable Python access to owned harmonic settings. Comments and new API
placement were checked against nearby Tudat code. The independent temporal
numerical checks use three cases, with frozen arithmetic budgets established
before running the implementation. ``eba17ab61`` also narrows the two new
standalone AOD1B/temporal model test links so they no longer request the full
propagation library list. This CMake-only cleanup is now applied in the compiler
checkout; the two links are checked alongside the kernel build. The tracked
source files match ``eba17ab61`` exactly; the compiler checkout retains the
small patch uncommitted until its build finishes, to avoid changing generated
Git metadata during verification. It does not change a numerical expression.
Both narrowed targets now link and pass: AOD1B 2 cases/12 assertions and temporal
ICGEM 3/63. The final retained native results cover fourteen distinct executables,
68 cases and 3069 assertions, with zero failures. The machine-readable record is
``verification-20261015/retained-latest-native-summary.json``. No further native
runs are pending. The kernel build is active; Python runtime validation remains.

Literature and model decisions
-----------------------------

The following is a design analysis, not a declaration that the listed models are
interchangeable or that all optional products will become defaults.

* The `official IERS Conventions 2010
  <https://iers-conventions.obspm.fr/conventions/content/tn36.pdf>`_ and separately
  identified `chapter 6 updates
  <https://iers-conventions.obspm.fr/chapter6.php>`_ supply the conventional
  baseline. Updated pages are explicitly described by IERS as under development;
  they must not be labelled a newly ratified edition. Record the exact selected
  update, especially the 2018 secular-pole change, rather than a vague "latest
  IERS" option. The FES2004 update history includes coefficient and phase fixes
  and an unresolved very-long-period phase discussion. Pin the data revision and
  expose that choice; do not silently repair the provider file.
* `IERS chapter 6, sections 6.2--6.5
  <https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf>`_ distinguishes
  direct solid-body response, degree-two forcing of degree-four coefficients,
  frequency corrections, permanent-tide treatment, ocean tides and two pole
  tides. The implementation must preserve these distinct contributions. A static
  tide-free field accepts the permanent solid tide; a zero-tide field requires
  its removal from the added correction. The full ocean pole response is a
  harmonic expansion, not merely its dominant C21/S21 approximation.
* The `IERS version history
  <https://iers-conventions.obspm.fr/conventions_versions.php>`_ identifies the
  transition from mean-pole to secular-pole conventions. Tudat's station models
  already use the 2018 secular-pole helper. Gravity pole tides should consume
  that same convention. An explicitly reproduced older reference uses its own
  declared pole model; its results must not become modern defaults by accident.
* `CNES/GRGS FES2014b spherical harmonics
  <https://grace.obs-mip.fr/dealiasing_and_tides/ocean-tides/>`_ distinguish water
  heights from potential coefficients and provide versions with and without S1.
  Their release notes document corrected MSf, land masks, Doodson identifiers and
  phase conventions. Therefore "FES2014b" alone is not sufficient provenance.
  Preserve the conversion provider and revision, reference constants, loading
  response and constituent inventory. This is the first modern published
  coefficient product to validate after the conventional reference.
* `FES2022b from AVISO
  <https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes/release-fes22.html>`_
  is a modern ocean atlas. Grid heights and the prediction library do not by
  themselves define the geopotential coefficient conversion. A later explicit
  FES2022 profile must identify its conversion, ellipsoidal treatment, load Love
  numbers and admittance. Reuse the harmonic evaluator; do not create a separate
  evaluator merely because the atlas name changes.
* The `TiME22 documentation, revision 1.1
  <https://datapub.gfz-potsdam.de/download/10.5880.GFZ.1.3.2022.006-VENonf/2022-006_Sulzbach-et-al_TiME22_README.pdf>`_
  separates atmospheric and oceanic components, includes loading in its potential
  coefficients, and absorbs the phase convention into its cosine/sine fields.
  Its revised admittance maps 57 input tides to 435 constituents. It must not
  receive another phase shift or load-response factor. Its atmospheric and oceanic
  S1 definitions differ. Treat it as an alternative product profile, retaining
  the supplied constituent identities and admittance rather than substituting
  the FES2004 table.
* `GFZ AOD1B
  <https://isdc.gfz.de/esmdata/aod1b/>`_ supplies non-tidal atmospheric/oceanic
  gravity variability. `Release notices
  <https://isdc.gfz.de/esmdata/aod1b/news>`_ document separate tidal products and
  corrected RL07 files. Record release and file revision. The
  `component definitions
  <https://www.gfz.de/en/section/earth-system-modelling/esmdata/esmdata/aod1b-product/low-degree-coefficients>`_
  make GLO the sum of ATM and OCN; these three must not be added together. OBA
  represents a different domain/physical combination, not a synonym for OCN.
* `ICGEM's 2023 format specification
  <https://icgem.gfz-potsdam.de/docs/ICGEM-Format-2023.pdf>`_ carries GM, radius,
  normalization, tide system and temporal coefficient records. Existing numeric
  matrix readers cannot supply this information after discarding the header.
  Preserve metadata in neutral input records and honour piecewise validity and
  reference epochs when adapting temporal terms to the existing variation API.

Existing-model comparison
-------------------------

Tudat already has basic complex Love-number response, real mode coupling,
tabulated, periodic and polynomial coefficient variations, a variation collection,
time-dependent spherical harmonics and force evaluation. The implementation uses
these existing interfaces. At the initial audit, the conventional solid, ocean
and pole enum names lacked their corresponding factory implementations; those
factories are now present. The local legacy ILR code supplied reusable numerical
structure for the solid, ocean-pole and harmonic models. Its unrelated framework
and obsolete build machinery were not imported.

The pinned `GROOPS source
<https://github.com/groops-devs/groops/tree/7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee>`_
separates solid, pole, ocean-pole and Doodson-harmonic models. Its harmonic path
binds an admittance matrix to an ordered major-constituent list and rescales
coefficients for the consumer's GM/radius. These are useful comparison cases.
Its documented/default data files must be pinned independently of the software.
GROOPS remains an external reference; its GPL source is not copied into Tudat.

`Orekit 13.1.8 SolidTidesField
<https://gitlab.orekit.org/orekit/orekit/-/blob/13.1.8/src/main/java/org/orekit/forces/gravity/SolidTidesField.java>`_
provides another implementation of the three-step solid-tide construction. It
explicitly receives the central field's tide system and body frame. Use that
separation to check our design, without importing its unrelated caching or API.
Reference execution has not yet been performed.

Setup compatibility rules
-------------------------

These enforcement rules follow the analysis above and the existing Tudat setup
paths. They are implemented in the Earth-specific consuming factory, including
the temporal baseline checks. Neutral readers do not depend on
bodies, dynamics or estimation.

.. list-table:: Compatibility decisions
   :header-rows: 1
   :widths: 25 75

   * - Boundary
     - Required behaviour
   * - Rotation and terrestrial axes
     - Coefficients and tide-generating body positions must use the selected
       gravity frame. Reject an Earth high-accuracy setup with an incompatible
       rotation target or inertial base. Reuse the selected GCRS/ITRS rotation;
       no independent hidden Earth rotation or extra J2000 frame bias.
   * - Rotation clock and EOP
     - Consume the rotation's converter and EOP/polar-motion calculation. Route
       TT and UT1 explicitly to fundamental arguments; do not use a convenience
       routine that substitutes UTC for UT1 or TDB for TT. The production Earth
       gravity factory requires the selected rotation to accept TDB, the existing
       default. Adapting only the gravity callback would leave the body's cached
       rotation, and hence the final force transformation, at the wrong epoch.
   * - Pole convention
     - Both gravity pole tides and station pole loading use the same secular-pole
       convention and EOP source. A conventional high-accuracy preset rejects a
       simple/SPICE rotation lacking the required EOP contract rather than
       silently creating a different polar-motion model.
   * - Static tide system
     - Require explicit known tide-free or zero-tide metadata for the conventional
       solid-tide preset. Unknown and mean-tide inputs need an explicit supported
       conversion; do not infer a tide system from the sign or magnitude of C20.
   * - GM, radius and normalization
     - Normalize once at ingestion and rescale once to the consuming field using
       the degree-dependent GM/radius relation. Keep model constants distinct
       from output scaling constants. Never apply a loading factor to coefficients
       that already include the induced potential.
   * - Coordinate scaling and origin
     - Identify the field/ephemeris coordinate convention and center of mass.
       Handle degree zero/one only when the product and dynamics explicitly
       support them. Do not turn a station displacement or water-height grid
       into a gravity correction without its prescribed physical conversion.
   * - Solid-tide overlap
     - A complete conventional Earth model replaces the basic Sun/Moon response
       for those same bodies. Reject duplicate complete/basic contributions;
       unrelated generic-body models remain available.
   * - Ocean EOP versus ocean gravity
     - Desai--Sibois replaces the corresponding ocean EOP terms only. It does
       not replace ocean gravity tides or either pole tide. Ocean atlas names
       alone are not a compatibility rule: compare actual conventions and
       inventories, not name equality.
   * - Ocean/atmosphere overlap
     - Admit one contribution per physical component and constituent. Distinguish
       ocean mass from direct air mass and radiational ocean response. Select
       S1 and de-aliasing components according to the provider release; reject
       GLO plus ATM/OCN and duplicated radiational ocean signals.
   * - Minor constituents and nodal effects
     - Bind an admittance matrix to its exact ordered major-wave inventory.
       Apply each minor/nodal contribution once. Do not combine explicit
       sidebands with an additional modulation that represents those sidebands.
   * - Secular/time-variable baseline
     - Temporal records state whether they replace a coefficient baseline or add
       an anomaly. Do not add a full time-variable field to its own static mean,
       duplicate an already included trend, or extrapolate beyond validity.

Keep the API close to existing ``gravity_field_variation`` settings and
``data_input.environment_data`` readers. Avoid a new universal configuration
framework or another epoch wrapper. The current generic variation/update API
uses double epochs; quantify the coefficient/force error before extending generic
interfaces for additional precision. Converters and internal continuous clocks
must retain ``Time`` where already supported. No precision claim may exceed the
actual public evaluation path.

Validation design and provenance
--------------------------------

The `COST-G benchmark publication
<https://adgeo.copernicus.org/articles/55/1/2020/>`_ provides independent force
values and shared auxiliary inputs, with a published 1e-11 m/s^2 comparison
criterion. Use selected model-specific samples, not a new trajectory propagation
campaign. Its `dataset README
<https://ftp.tugraz.at/pub/ITSG/COST-G/softwareComparison/00README_simulation.txt>`_
explicitly uses UTC instead of UT1 in its Doodson helper, includes permanent solid
tides, and specifies its secular pole and ocean-pole load Love numbers. A matched
input test can validate synthesis; it cannot validate a different clock policy.
Retain separate production-factory tests with correctly routed TT/UT1.

The `official ocean-pole test
<https://iers-conventions.obspm.fr/content/chapter6/additional_info/desaiscopole2004.txt>`_
contains explicit constants, a historical pole model, wobble and C21/S21 outputs.
Test the coefficient evaluation using its declared wobble. Do not compare its
historical pole prescription directly against the modern secular-pole default.

Tests should cover independent numerical values, phase/sign/normalization,
permanent-tide consistency, real factory wiring and material overlap rejection.
Retain a small set of meaningful cases rather than one test for every argument
check. Reuse existing spherical-harmonic force code and passing native evidence.
Python tests cover registration, settings ownership and imports only.

External reference cache:
``/home/dominic/Tudat/earth-station-model-work/references/gravity``.
``sources.json`` records retrieved URLs, SHA256 and byte sizes; ``costg/sources.json``
does the same for benchmark data. Downloaded reference sources are analysis inputs,
not vendored production dependencies. Any committed fixture extract must record
its exact source ranges, conventions, licence and reproduction command. Ordinary
tests must be entirely offline.

Remaining work at this checkpoint
---------------------------------

* Finish the temporal consuming adapter and validate its explicit date/year and
  nominal-baseline choices.
* Complete the focused native tests for the later AOD1B, TiME22 and temporal
  checkpoints, retaining the numerical tolerances set before execution.
* Finish the corresponding Python bindings, compile the kernel for their tests,
  and check documentation imports and the relevant architecture tests.
* Finish the review of added code, comments, API placement and offline test data.
  Examples remain excluded.

Implementation contracts
------------------------

The first implementation is an ICGEM static-field reader and explicit tide-system
metadata. This is independent of the Earth model formulae and prevents a later
factory from guessing the permanent-tide treatment from C20. Raw readers live in
``tudat::input_output`` and retain the published normalization, error columns,
header, source name, GM, radius and declared degree. Their data do not import
simulation or dynamics. The consuming settings adapter normalizes coefficients
once using the existing Legendre normalization utility. It must not manufacture
C00, silently discard temporal records, or confuse a missing coefficient with a
published zero. The initial static adapter rejects temporal records explicitly;
the later temporal adapter owns reference epochs and piecewise validity.

A small shared ``TideSystem`` enum (unknown, tide-free, zero-tide, mean-tide) is
carried by spherical-harmonic settings and the constructed field. Existing
constructors remain valid with an unknown default. The field's existing frame,
GM and radius remain authoritative; no parallel universal metadata/configuration
framework is introduced. Earth solid-tide construction reads the actual field's
tide system. Unknown/mean-tide values are allowed as input records but rejected
by that high-accuracy factory until an explicit conversion is supported.

The nearby gravity/rigid-body branch was inspected. Its broad derivative and
rigid-body changes are unnecessary for these models. The conventional degree-two
and degree-three forcing already exists here and is the first reuse choice.
Degree-four response to degree-two forcing must use the degree-two Legendre
functions in equation 6.7, not a new degree-four tide raiser. The unrelated
arbitrary-degree recurrence will only be imported if a concrete model needs it.

Factory helpers remain local to Earth gravity setup. They obtain the selected
``GcrsToItrsRotationModel``, its converter and polar-motion calculator, check the
field's target axes against the rotation, and require its input scale to be TDB.
This also protects the existing body rotation and force-evaluation path without
adding Earth-specific time machinery to the generic Body class. Body ephemeris base axes must agree with the
rotation's base axes; already supported explicit base-frame rotations remain
available. No extra frame-bias rotation is added. The full conventional preset
uses the selected corrected pole with the shared 2018 secular pole; matched
benchmark tests inject their published wobble instead of changing production
defaults. Setup rejects complete/basic solid-tide overlap for the same raiser.

The existing construction order is correct. Variation settings are attached with
correction initialization deferred. ``processBodyFrameDefinitions`` first
installs the global-frame ephemeris callbacks, then initializes the correction
functions and samples interpolated variations. An actual BodySettings test with
an Earth-relative lunar ephemeris passes 1 case/4 assertions. Preserve this
existing mechanism; moving frame initialization earlier would be incorrect.

Dominic identified the legacy implementation under
``/home/dominic/Software/ILR/trunk/ILR/ILR/Astrodynamics/Gravitation``.
It contains Earth solid, solid pole, ocean pole and ocean-tide implementations,
with IERS frequency tables and FES2004 data in the ILR DataFiles tree. The audit
is complete: all 71 frequency-table rows agree with the pinned primary tables.
The solid-tide implementation reuses the old degree-two-to-four construction
and frequency corrections alongside Tudat's current geometry and coefficient
helpers. The diurnal S21 sign is corrected, and the old implicit secular trends
are omitted. The audit also identified an old solid-pole coefficient, ocean-pole
degree-index/constructor defects and an incorrect six-character Doodson decoder.
Those defects are not carried into the new models. FES2004 data remain useful
for validating the harmonic reader's distinct product convention.

The solid-Earth and solid-pole backends are integrated at ``1ea8d4c6b`` and the
ocean-pole backend at ``9e7a2d936``. The solid/pole native suite passes 3 cases
and 34 assertions, including four COST-G force samples: maximum errors are
8.94e-17 and 3.45e-16 m/s^2, respectively, against a fixed 1e-11 m/s^2 budget.
The ocean-pole reader passes 2 cases/18 assertions and its model passes 3/830,
including all 366 published daily reference rows and four degree-180 COST-G
force samples (maximum error 2.57e-21 m/s^2). Independent table/fixture
reproduction and the isolated reader target checks pass. Integration retained
the tested numerical code unchanged; its only conflict combined adjacent test
registrations, so the numerical tests were not repeated.

The static ICGEM foundation is integrated at ``3d7c7d59c``. The reader preserves
sparse coefficient presence, declared normalization/error columns, source rows
and tide-system metadata. The consuming adapter requires an explicit frame and
complete selected coefficients unless zero-fill is requested; C00 is always
required. Normalization is applied once. Tests cover a real GFZ snapshot and its
comment preamble, the published COST-G duplicate-C00 failure, normalization,
static/time-dependent factory metadata and the existing numeric reader's
truncated-order defect: raw reader 4 cases/32 assertions, consuming factory
2/26, plus two architecture checks. The actual raw-only target links the new
standard-library-only reader. No Python kernel was built for these native tests.

The existing scalar coefficient-evaluation interface is retained. At epochs
within roughly 100 years of J2000 its double rounding is below 0.25 microseconds;
for a diurnal coefficient of amplitude 1e-8 this contributes less than 2e-19 to
the coefficient (using the local phase rate). This is below the 1e-11 m/s^2
benchmark force budget. This estimate does not justify extrapolation, event
snapping or a claim of native-Time public gravity evaluation. Internal physical
clock conversion and boundary selection retain Time where the existing API
provides it.

Conventional solid and pole models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conventional solid-tide implementation uses the anelastic parameters of
table 6.3: degree-two complex Love numbers 0.30190, 0.29830-0.00144i and
0.30102-0.00130i; degree-three values 0.093, 0.093, 0.093 and 0.094. The
degree-two-to-four responses are -0.00089, -0.00080 and -0.00057. Reuse the
existing per-body geometry and normalized single-coefficient forcing helper.
Do not import a new generic rate interface or disguise the induced response as
degree-four direct forcing. A small Earth-specific loop over these fixed terms
is preferable to changing the generic Love-number API.

Tables 6.5a/b/c supply the frequency corrections, once per Earth evaluation,
not once per tide raiser. For Doodson phase theta and table amplitudes ip/op,
the corrections are C20 += ip*cos(theta)-op*sin(theta),
C21 += ip*sin(theta)+op*cos(theta), S21 += ip*cos(theta)-op*sin(theta),
and C22 += ip*cos(theta), S22 -= ip*sin(theta). Units are 1e-12.
The tables contain 48 diurnal, 21 zonal and two semidiurnal terms. Pin the
numerical extraction against the primary IERS table; the cached Orekit 13.1.8
resource files provide a second transcription check with their Apache notices.

For a zero-tide nominal field, subtract the permanent C20 contribution
4.4228e-8 * (-0.31460) * 0.30190 exactly once. A tide-free field retains it.
Keep the Earth model's conventional constants explicit and apply the
GM_source/GM_field * (R_source/R_field)^n output conversion once where required.
Do not multiply by an additional normalization or loading factor.

Solid pole gravity uses the IERS section 6.4 C21/S21 equations with wobble in
arcseconds. Ocean pole gravity uses the full section 6.5 equation 6.23 with
the supplied complex A/B coefficient table, explicit load Love numbers and
wobble in radians. It must expose selected maximum degree/order and require
all needed load coefficients. Standard geocentric gravity uses degrees >=2;
degree-one geocentre physics is not silently added to this body-centered force.
The short degree-two equation 6.24 is labelled an approximation if offered.

Selected COST-G force samples provide an independent 1e-11 m/s^2 target using
their exact supplied rotation, pole, constants and astronomical arguments.
The solid-tide benchmark uses DE421: generate and freeze only the required
geometric Sun/Moon vectors, not a runtime ephemeris download. Separate factory
tests check proper production TT/UT1 routing and tide-system/overlap enforcement.

Harmonic ocean products
~~~~~~~~~~~~~~~~~~~~~~

Use one harmonic evaluator for declared potential coefficients, with canonical
cosine/sine amplitudes for each Cnm and Snm. For IERS potential-format prograde/
retrograde columns Cp, Sp, Cm, Sm, equation 6.15 gives
Ccos=Cp+Cm, Csin=Sp+Sm, Scos=Sp-Sm, Ssin=Cm-Cp. This transformation also handles
zonal products whose retrograde coefficients are zero and prograde coefficients
already doubled. Set physical Sn0 to zero; do not double the zonal terms again.
Potential-format files already absorb their published phase/loading conversion.
Water-height tables are a different input product and are not accepted as
potential coefficients merely because their columns look similar.

The FES2014b-v1 CNES/GRGS potential file is pinned separately from COST-G's own
FES2014b conversion. Its S1 is the radiational 164.555 constituent and Sa is
056.555; preserve those identifiers. The benchmark's EOT11a potential file has
a contradictory leftover 'in cm' header alongside its explicit 1e-11 Stokes
units; support it only with the identified provider profile, never by guessing
from the filename. Source GM/radius, normalization, component, revision and
included loading are explicit product metadata when absent from the header.

Admittance takes an explicit ordered major inventory, target inventory and
matrix. Validate the inventories and synthesize each output constituent once.
Do not mix such explicit sidebands with another nodal modulation. Published
FES2004 interpolation and TiME22-v1.1 matrices remain separate named profiles;
the latter already includes phase and loading in its cosine/sine fields.

The complete official TiME22 revision 1.1 archive was retrieved on 15 October
2026 into the local reference cache. Its SHA256 is
``79c9ef695b4c15d45f5a35dd706a7ed2e1df5acc4d7e509b5f6e404a71467ab3``;
the archive has 153 files, with individual checksums in
``references/gravity/time22-v1_1/files.json``. The provider permits reuse under
CC BY 4.0. Only necessary documented fixtures will enter the repository.
The ocean matrix has 57 input rows and 435 output columns and must be transposed
for the existing target-by-major helper. The atmospheric matrix is 16 by 16.
The supplied multiplier tables, rather than filename parsing alone, define
the output phases. Real S5/S6 identifiers include ``5A0.555`` and ``6BZ.555``;
section 4.3 defines A=10, B=11 and Z=-1 before the usual Doodson offsets.
The current setup work includes support for these verified identifier codes.

Non-tidal and temporal adaptation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AOD1B RL07 consumption keeps the provider GPS sampling grid and reuses
``TabulatedGravityFieldVariations`` with linear interpolation and no
extrapolation. Each query is converted from the selected rotation's TDB clock
to GPS before interpolation; converting only the sample knots to TDB would
change the interpolation weights. Missing three-hour intervals reject.
Published anomalies are added relative to the provider's 2007--2014 mean and
scaled once to the consuming field's GM and radius. Standard Earth setup starts
at degree two. GLO cannot accompany ATM or OCN, and OBA cannot be combined with
another AOD component without a supported explicit physical decomposition.

The ICGEM temporal audit confirms that ``gfct``, trend and periodic terms define
replacement coefficient histories. An adapter must subtract the consuming
nominal baseline once, without subtracting nearly equal complete coefficients
on every evaluation. The 2023 specification does not define the time scale or
duration of its year; these must be explicit consuming choices. Piecewise
coefficient intervals and their endpoints require a documented rule. This
neutral reader retains temporal records, while the static
consumer rejects them. The consuming adapter uses one
coefficient/interval model: existing polynomial and periodic models have one
global reference epoch and cannot directly represent the published independent
intervals. Tabulation is unsuitable across coefficient jumps or coverage gaps.

Consumption requires an explicit source time scale and year convention. The
Julian-year option uses elapsed TAI seconds divided by 31557600; the fractional
Gregorian-year option uses the calendar coordinate in the declared source scale.
Format 2.0 intervals are half-open, with the right-hand record owning a shared
boundary. Its date tokens require HHMM. Format 1.0 date-only records use midnight
and have no artificial validity limit.
An optional map of printed date tokens to source-scale epochs permits an explicit
caller correction of irregular provider dates. In particular, EIGEN's
``20041226.0060`` is preserved by the reader and is not silently rolled to 01:00.
The nominal field at construction defines the added anomaly. Later intentional
nominal-field changes act as offsets, as with existing additive variations.
No deferred interpolation or extra epoch wrapper is introduced for this model.

The primary EIGEN-6S4 version 2 file is now pinned in the reference cache
(``eigen6s4-v2-source.json``). Its independently dated periodic records confirm
that one cannot attach every term to an identical ``gfct`` interval: twelve
published periodic records have longer validity intervals than their associated
baseline records. This is part of the temporal adapter design, not a reason to
change or approximate the published coefficients.

Ocean-pole backend checkpoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``readOceanPoleGravityData`` reads the distinct IERS normalized complex A/B
integrals into neutral sparse records. It preserves headers/source lines,
published degree-one values and missing-entry identity; duplicate, malformed,
nonfinite or physically nonzero zonal sine integrals reject. Its independent
reader uses only the C++ standard library and is now part of the consolidated
``tudat_environment_data`` target.

``OceanPoleTideGravityFieldVariations`` implements equation 6.23 for explicitly
selected degrees >=2 and orders including zero/one. It requires all selected
integrals and degree-indexed load Love numbers, an explicit radians wobble
callback at TDB seconds, physical forcing constants, and output field GM/R.
The constructor freezes the scaled response; evaluations combine its four
matrices with the two wobble components. Above-triangle entries remain zero.
No factory/default pole model is selected at this backend checkpoint.

The local legacy ILR kernel supplied the same useful A/B response structure;
its uninitialized base dimensions, index-versus-size bounds and shifted load
response indices are corrected. All 366 official 2004 gravity-reference epochs
pass with budgets derived from published decimal precision. Four exact COST-G
samples use degree/order 180, the supplied historical secular pole/EOP,
Gegout97 load Love numbers and rotation/position records. Their maximum force
residual is 2.57e-21 m/s^2 against the fixed 1e-11 m/s^2 budget. Separate checks
cover selected blocks, GM/R scaling, immutable response and missing coverage.
Pinned source hashes, exact extraction ranges and an offline reproduction script
are in ``tests/test_tudat/reference_files/ocean_pole_gravity``.

Harmonic backend checkpoint
~~~~~~~~~~~~~~~~~~~~~~~~~~

``readHarmonicPotentialData`` preserves the ordered raw potential inventory,
source IDs, headers/units and sparse Cp/Sp/Cm/Sm records. Identified FES2004,
CNES/GRGS FES2014b-v1 and COST-G EOT11a profiles retain their published wave
extents; unknown units, duplicate records and incompatible profile inventories
reject. FES2004's 1e-12 scale remains distinct from the other profiles' 1e-11.
EOT11a's contradictory 'in cm' text requires its explicit identified profile.
No source GM/R is inferred. The neutral target uses only the standard library.

The single conversion helper applies equation 6.15 and physical Sn0=0. Selected
Doodson IDs are unique and must exist, and their explicit supplied order is
retained; the default selects every wave in raw order. Missing selected records
inside a declared extent reject. Only declared outer extents become zero for
identified profiles; generic sparse data needs an explicit zero-fill decision.
The common canonical evaluator receives declared source GM/R, normalization,
component, revision/loading and a Doodson callback, and applies output GM/R
scaling once. It adds no phase, loading or nodal modulation. Ordered admittance
synthesizes a replacement target inventory once.

Three reader cases and four evaluator cases pass. Closed-phase checks cover the
IERS signs, zonal convention, reference scaling and ordered admittance. Four
COST-G M2 degree/order-120 force samples agree within 1.57e-13 m/s^2 against the
fixed 1e-11 target with supplied benchmark arguments/positions/rotations and an
explicitly labelled benchmark GM/R convention. This does not supply missing
provider reference constants. Source ranges/hashes and offline extraction are
in ``tests/test_tudat/reference_files/harmonic_potential``. Consuming Earth setup
and Python registration remain separate checkpoints.
