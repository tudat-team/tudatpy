Earth gravity-field variations: design and implementation tracker
================================================================

Scope and status
----------------

Branch: ``fix/gravity_field_variations``, created from the Earth rotation/station
checkpoint ``a3be1612a`` on 14 October 2026. Dominic authorized this programme to
start in parallel with the remaining station validation and review. The two
programmes have separate integration worktrees. Reviewed station fixes will be
merged before gravity setup depends on them.

Status: primary-literature review and existing-code inventory are complete. The
first implementation contract below is ready; model-specific product checks
continue before their implementations. No new gravity model is claimed validated. The coordinator designs,
reviews and handles major debugging; GPT-6.1 Sol agents at high reasoning implement,
compile and handle minor debugging. The global limit is eight compiler jobs across
all worktrees. Only relevant native tests run; the Python kernel is built solely
for Python-specific tests. No new estimation parameters or gravity recovery work.

.. code-block:: text

   [x] G01 Review primary literature, provider conventions and existing implementations.  [DONE]
   [ ] G02 Define rotation/gravity compatibility and reuse existing Tudat functionality.  [ACTIVE]
   [ ] G03 Preserve neutral coefficient, constituent and source-convention records.  [ACTIVE — static ICGEM foundation]
   [ ] G04 Complete conventional solid Earth gravity tides, including frequency dependence.  [QUEUED]
   [ ] G05 Add solid Earth pole-tide gravity corrections with the shared pole convention.  [QUEUED]
   [ ] G06 Add ocean pole-tide gravity corrections with explicit truncation and origin.  [QUEUED]
   [ ] G07 Add harmonic ocean-tide synthesis and documented minor-constituent admittance.  [QUEUED]
   [ ] G08 Support identified ocean-tide products without mixing their phase conventions.  [QUEUED]
   [ ] G09 Handle atmospheric/non-tidal coefficient products and component overlap explicitly.  [QUEUED]
   [ ] G10 Connect secular/periodic/tabulated variations to their static-field baseline.  [QUEUED]
   [ ] G11 Integrate Earth factories with the selected rotation and enforce compatibility.  [QUEUED]
   [ ] G12 Validate individual models and assembled forces against pinned references.  [QUEUED]
   [ ] G13 Complete Python settings, diagnostics and a reproducible setup example.  [QUEUED]
   [ ] G14 Review the complete gravity diff for correctness, reuse, API clarity and size.  [QUEUED]

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
time-dependent spherical harmonics and force evaluation. Extend these seams.
The ``iers_2010_tidal``, ``ocean_tide`` and ``pole_tide`` enum names currently do
not have matching factory implementations. Do not present enum exposure as model
support. A separate local gravity/rigid-body branch contains candidate numerical
reuse; inspect individual routines rather than importing its unrelated machinery.

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

These are proposed enforcement rules derived from the analysis above and the
existing Tudat setup paths. They will be checked in the Earth-specific consuming
factory; neutral readers must not depend on bodies, dynamics or estimation.

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
       routine that substitutes UTC for UT1 or TDB for TT. A non-TDB rotation
       must be adapted at the factory, not called with an incorrectly labelled
       propagation epoch.
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

Open decisions before implementation
-----------------------------------

* Finish the local-branch reuse assessment and final Earth factory/API design.
* Verify coefficient-file licences and precise product headers before choosing
  fixture extracts and packaged optional profiles.
* Derive model-specific coefficient and acceleration tolerances from published
  reference precision and truncation, then freeze them before testing.
* Decide the smallest supported atmospheric/non-tidal input set after checking
  real product records. Do not manufacture missing scientific metadata.
* Preserve an explicit difference between a low-degree pole approximation and
  the full declared harmonic expansion.

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
field's target axes against the rotation, and adapt the propagation TDB epoch to
the rotation's declared input scale. Body ephemeris base axes must agree with the
rotation's base axes; already supported explicit base-frame rotations remain
available. No extra frame-bias rotation is added. The full conventional preset
uses the selected corrected pole with the shared 2018 secular pole; matched
benchmark tests inject their published wobble instead of changing production
defaults. Setup rejects complete/basic solid-tide overlap for the same raiser.

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
