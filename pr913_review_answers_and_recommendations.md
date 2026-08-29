# PR 913 review answers and recommendations

## Purpose and review principle

This report consolidates:

1. the answers to the review questions raised after the initial implementation;
2. the original list of proposed changes produced from those questions; and
3. the revised recommendations after independently assessing correctness, complexity, compatibility, and maintenance cost.

The review questions are treated as prompts for investigation, not automatic change requests. A change belongs in the eventual implementation plan only when it has a clear net benefit. In particular, moving `scaled_mean_moment_of_inertia` is a gated design item: it must not be implemented until there is a clear, concise solution with the required backwards compatibility.

## Answers to the review questions

### 1. Why does the gravity field own `scaled_mean_moment_of_inertia`?

It should not be part of the runtime gravity-field model. The current implementation stores the value in `SphericalHarmonicsGravityFieldSettings`, copies it into `SphericalHarmonicsGravityField`, and lets the gravity field decide whether it can return an inertia tensor. Although `RigidBodyProperties` caches the resulting tensor and derivative, the gravity model remains authoritative for an inertia property.

The intended ownership is:

- the gravity field owns gravitational parameter, reference radius, coefficients, reference frame, and gravity variations;
- gravity-linked rigid-body properties own the scaled mean moment, current inertia tensor, and inertia-tensor derivative; and
- the rigid-body properties retain a read-only link to the gravity field as the source of the coefficients needed for the conversion.

This migration is potentially broad. It also affects environment settings, automatic body creation, estimation of the mean moment, polyhedron handling, and C++ construction paths. It is therefore a deferred design item rather than an approved implementation change.

### 2. Is `docs/tudatpy/source/dynamics/propagation_setup/gravity_deformation.rst` a new Python module?

No. The new public entry points are functions in the existing modules:

- `tudatpy.dynamics.propagation_setup.create_gravity_deformation_models`; and
- `tudatpy.dynamics.propagation_setup.propagator.gravity_deformation` and `maxwell_deformation`.

The standalone page is a narrative feature guide, but it was added to a toctree captioned “Modules”. That incorrectly presents it as a Python module.

The API-specific material should be in the existing function docstrings. Broader lifecycle and conceptual material should be in the Tudat user guide. If a compact API overview is still useful, it should be integrated into the existing `propagator.rst` page rather than represented as a new module page.

### 3. Why does `GravityFieldModel` have `getInertiaTensor()`?

The method exists because the current implementation uses gravity fields as generic inertia providers:

- a spherical-harmonic gravity field derives inertia from its degree-two coefficients and scaled mean moment; and
- a polyhedron gravity field derives inertia from its geometry and density.

`FromGravityFieldRigidBodyProperties` then retrieves and stores the matrix. This is convenient, but it makes the gravity model responsible for calculating and exposing inertia.

Architecturally, gravity-linked rigid-body properties should instead retrieve the necessary read-only source data and perform the conversion themselves, possibly through free calculation functions. However, removing `hasInertiaTensor()` and `getInertiaTensor()` must be considered together with the scaled-mean-moment ownership redesign; it should not be performed independently.

### 4. Why did `coupledStateDerivativeSolver.h` exist?

There is a genuine same-epoch algebraic feedback loop:

```text
gravity-coefficient rates -> inertia-tensor derivative
                           -> angular acceleration
                           -> gravity-coefficient rates
```

More specifically:

- the Maxwell coefficient rates depend on angular acceleration through the derivative of the centrifugal equilibrium contribution;
- the inertia-tensor derivative is a linear function of the degree-two coefficient rates; and
- angular acceleration depends on the inertia-tensor derivative through the inertial torque.

Sequential evaluation would use a stale derivative and make the result depend on model evaluation order.

The removed solver represented this as a fixed point `d = F(d)`. It numerically reconstructed a scaled affine coupling matrix, attempted a dense direct solve, validated the result, and fell back to fixed-point iteration. It also returned an implicit multiplier intended to keep the variational equations consistent.

The physical need for a simultaneous solution is real. The generic global solver was nevertheless more general than required because the supported Maxwell/inertial-torque loop is affine and body-local. It has now been removed, and the earlier bounded derivative-update loop has been restored pending a focused design.

### 5. Why was `coupledStateDerivativeSolverSettings_` not allowed to be null?

There is no domain-level reason for the unconditional restriction. It was added to guarantee that the pointer can be dereferenced in the coupled branch.

Ordinary propagation normally has no coupled state-derivative updater and never uses these settings. In that case, a null pointer is valid and avoids allocating irrelevant configuration.

If such a generic solver were retained, a better contract would be:

- allow null settings when no coupled dependency exists;
- create internal defaults only after coupling is detected, or require settings only at that point; and
- do not expose a mandatory solver object on every propagator setting.

The member, mandatory non-null check, and Python API have now been removed with the generic solver.

### 6. Why was `stateDerivativeIndices` changed from an `IntegratedStateType` key to a `StateDerivativeDependency` key?

The new updater iterates dependency-specific environment update functions. Using the same dependency as the index-map key avoids converting back through `getStateTypeForDependency()` and would prevent two different dependency types sourced from the same integrated state type from sharing and potentially mixing one index vector.

For the currently supported dependencies, however, the mapping is one-to-one:

```text
inertia_tensor_derivative_dependency -> gravity_deformation_state
rotation_rate_derivative_dependency  -> rotational_state
```

The key change was therefore structurally defensive, not functionally necessary today. It has been reverted together with the generic solver, restoring the `IntegratedStateType` key used by the earlier updater.

### 7. What is the role of the large coupled block in `DynamicsStateDerivativeModel::computeStateDerivative()`?

The block is the integration layer around the generic coupled solver. It:

1. discovers the five gravity-rate and three angular-acceleration entries involved in coupling;
2. maps propagated and conventional state indices across every registered model and body;
3. builds an initial derivative and numerical component scales;
4. defines `F(d)` by injecting trial derivatives, updating the environment, refreshing dependent models, and recomputing the affected derivative blocks;
5. solves `d = F(d)`;
6. re-evaluates with the solution so the environment and full state derivative contain the same final values; and
7. retains `(I - dF/dd)^{-1}` for the subsequent variational-equation calculation.

The role was necessary under the generic solver design. Its implementation was intrusive: it repeatedly discovered relationships that are static during propagation, constructed a global dense problem, and duplicated some index bookkeeping already performed during setup.

The block and its implicit variational-equation multiplier have now been removed. The pre-solver bounded derivative-update loop is restored; any later replacement should be a focused physical design rather than another global generic solver.

### 8. Why are the propagated-gravity setters and `updateCurrentGravityField()` functions on `Body`?

The functions currently coordinate multiple objects:

- `setCurrentPropagatedGravityFieldVariation` finds the integrated variation, stores the live five-element correction, and refreshes the total field and gravity-linked rigid-body properties;
- `setCurrentPropagatedGravityFieldVariationDerivative` stores the live coefficient rate and forwards it to the gravity-linked rigid-body properties so the inertia derivative is updated; and
- `updateCurrentGravityField` updates the time-dependent gravity field, refreshes gravity-linked rigid-body properties, and outside propagation obtains a coefficient rate from installed history.

The synchronization is required somewhere. The functions are largely convenience/orchestration functions rather than fundamental new `Body` behavior, but that does not automatically mean they should be removed. `Body` can be a useful consistency boundary because it owns the participating environment objects.

The decision should be based on which alternative is simpler:

- retain one small internal/possibly private body-level synchronization operation; or
- create an internal per-body propagation interface holding direct links to the integrated variation, time-dependent gravity field, and gravity-linked rigid-body properties.

Replacing the methods is worthwhile only if the second design materially reduces coupling and duplication.

### 9. Was the old normalization and static-degree-two addition moved elsewhere?

The two parts were handled differently.

The unnormalised-to-geodesy-normalized conversion moved to `IntegratedGravityFieldVariations::calculateSphericalHarmonicsCorrections()`. The integrated state remains an unnormalised vector ordered as

```text
[Delta C20, Delta C21, Delta C22, Delta S21, Delta S22]
```

and the variation object converts it into normalized correction matrices for the time-dependent gravity field.

The explicit addition of `staticDegreeTwoCoefficients_` was not moved verbatim. Under the new variation semantics, the time-dependent field starts from its complete nominal coefficients and then adds the integrated variation and all other configured variations:

```text
total field = nominal environment field
            + integrated variation
            + all other gravity-field variations
```

A static contribution already present in the nominal field is retained automatically and must not be added again.

`static_coefficients` are now used inside the Maxwell constitutive model. They are converted to unnormalised form and subtracted from the current degree-two coefficients when defining the deformable baseline. This is more than a mechanical code move, so regression tests must establish whether historical Python workflows retain identical total-field behavior before any further change is made.

### 10. Requested comments on `IntegratedGravityFieldVariations`

Brief public-function documentation and focused implementation comments were added for:

- state ordering and normalization;
- live-state versus installed-history selection;
- default interpolation behavior;
- derivative evaluation from history; and
- propagation lifecycle behavior.

This comment-only change was committed as `066601762` (`Document integrated gravity field variation lifecycle`).

## Original proposed change list

The following was the original proposal produced directly after the review questions. It is retained here as historical context and is **not** the current approved plan.

1. Move inertia ownership fully into rigid-body properties.
   - Add `FromGravityFieldRigidBodyPropertiesSettings`.
   - Store the scaled mean moment in those settings and runtime properties.
   - Compute and cache inertia there using the linked gravity coefficients.
   - Move the estimatable mean-moment parameter away from the gravity model.

2. Remove inertia functionality from gravity-field models.
   - Remove `GravityFieldModel::hasInertiaTensor()` and `getInertiaTensor()`.
   - Remove the scaled mean moment from spherical-harmonic gravity constructors and runtime state.
   - Apply the same ownership rule to polyhedron-derived inertia.

3. Preserve Python compatibility through a deprecated forwarding property.
   - Add a canonical `rigid_body.from_gravity_field(...)` configuration path.
   - Retain the existing gravity-settings property as a forwarding compatibility alias.
   - Detect conflicting old and new configuration.

4. Remove propagation-specific gravity mutation methods from `Body`.
   - Replace the three methods with an internal per-body gravity-state interface.

5. Keep `IntegratedGravityFieldVariations` narrowly focused on live state, state rate, history, interpolation, and propagation status.

6. Replace the generic coupled solver with dependency-specific evaluation.
   - Use normal evaluation when there is no dependency.
   - Use ordered evaluation for one-way dependencies.
   - Use a direct per-body affine solve for a two-way Maxwell/inertial-torque cycle.

7. Remove generic coupled-solver configuration from the public API.
   - Remove the solver settings, failure policy, propagator member, and Python bindings if no nonlinear iterative solve remains.

8. Remove the large generic coupled block from `DynamicsStateDerivativeModel`.
   - Precompute body/model/index relationships during setup.

9. Simplify derivative-dependency bookkeeping.
   - Remove generic dependency components and parallel maps if per-body descriptors make them redundant.

10. Preserve and clearly define total-versus-variation and static-coefficient semantics.

11. Restructure the documentation.
    - Remove the false module page.
    - Move API material to function docstrings and conceptual material to the user guide.
    - Correct the scaled-mean-moment definition to `(Ixx + Iyy + Izz) / (3 M R^2)`.

12. Add focused tests for ownership, backwards compatibility, lifecycle, coupling, variational equations, multiple bodies, and historical Python workflows.

## Revised recommendations

### High-confidence changes

1. **Correct the documentation structure.**
   `gravity_deformation.rst` is not a Python module and should not appear in the module toctree. Its information should be redistributed to existing API docstrings and the user guide.

2. **Add regression tests before changing coefficient behavior.**
   Tests should explicitly verify nonzero nominal degree-two coefficients, nonzero `static_coefficients`, integrated variations, ordinary gravity variations, normalization, and post-propagation interpolation together.

3. **Do not require solver settings for uncoupled propagation.**
   If the current generic solver remains, null settings must be valid whenever no coupled dependency was created. This is a small correctness and API-cleanliness improvement, although it may become obsolete as part of the larger solver redesign.

4. **Keep the added `IntegratedGravityFieldVariations` comments.**
   They document behavior that is otherwise distributed across environment updating and result processing.

### Worthwhile redesigns, subject to confirming the full equations

1. **Represent the derivative dependency graph explicitly.**
   At setup, distinguish no dependency, a one-way dependency, and a true algebraic cycle for each body.

2. **Use ordered evaluation for one-way dependencies.**
   A solver is unnecessary when only coefficient rates depend on angular acceleration or only angular acceleration depends on the inertia derivative.

3. **Use a direct per-body affine solve for true cycles.**
   The currently supported Maxwell/inertial-torque equations appear affine in the five coefficient rates and three angular accelerations. This should be verified analytically. If confirmed, solve one compact block system per body and use the same block Jacobian for variational equations.

4. **Remove generic solver settings only after the direct formulation is complete.**
   The current solver and large dynamics block should not be deleted until the direct solution, final environment state, state-derivative modifiers, and variational equations are all demonstrably equivalent.

### Conditional changes—not justified automatically

1. **Do not automatically remove the body-level synchronization methods.**
   Keeping a small body-owned consistency operation may be cleaner than introducing another interface class. Change this only if the alternative substantially reduces dependencies and duplication.

2. **Do not automatically revert the dependency-keyed index map.**
   `StateDerivativeDependency` is a semantically valid key for dependency-specific consumers. It should remain if the generic updater remains and disappear only if the broader redesign makes it unnecessary.

3. **Do not change `static_coefficients` semantics based only on the code diff.**
   First reproduce historical workflows and determine whether the new variation state changes their total field or initial-state interpretation.

### Historical deferred, jointly designed inertia-ownership change

Moving `scaled_mean_moment_of_inertia` remains architecturally desirable but is not yet approved for implementation.

A promising two-stage design is:

1. Introduce gravity-derived rigid-body-property settings as the canonical owner.
2. Store the value in `FromGravityFieldRigidBodyProperties` at runtime.
3. Retain `SphericalHarmonicsGravityFieldSettings.scaled_mean_moment_of_inertia` in Python as a compatibility-only input carrier.
4. During body-settings reconciliation, transfer that legacy value into automatically created gravity-derived rigid-body-property settings.
5. Never copy the value into the instantiated gravity-field model.
6. Preserve existing explicit non-gravity rigid-body settings taking precedence over automatically generated gravity-linked settings.
7. Rebind mean-moment estimation to the rigid-body properties.
8. Treat removal of gravity-field inertia accessors and polyhedron changes as part of this same migration, not separate cleanup.

Before implementation, this design must satisfy all of the following:

- existing Python construction and property-assignment syntax produces the same environment;
- a gravity field remains valid without an inertia tensor;
- explicit rigid-body settings retain their existing precedence;
- estimation and torque partials remain correct;
- no duplicate source of truth exists at runtime;
- direct C++ construction and migration behavior have an agreed policy; and
- tests demonstrate compatibility rather than relying only on API shape.

This was the implementation gate at the time of the review. The design was subsequently approved with the additional reverse-link, update-callback, and Python-deprecation requirements recorded below.

## Current status

- The report records both the original analysis and subsequent review decisions.
- The generic coupled solver, its settings and Python API, dependency-keyed index map, large dynamics block, and solver-specific tests have been removed.
- The comments from commit `066601762` remain.
- The solver rollback was finalized separately in commit `420a3987a` (`Remove generic coupled state derivative solver`).
- The inertia-ownership migration was subsequently approved and implemented as a separate change.

## Implemented inertia-ownership design

The implemented design satisfies the approved ownership and compatibility requirements as follows:

1. `FromGravityFieldRigidBodyPropertiesSettings` is the canonical settings type, exposed in Python as `rigid_body.from_gravity_field(...)`.
2. `FromGravityFieldRigidBodyProperties` is the sole runtime owner of `scaled_mean_moment_of_inertia`, the derived inertia tensor, and its derivative.
3. `SphericalHarmonicsGravityFieldSettings.scaled_mean_moment_of_inertia` remains only as a deprecated Python-compatible input carrier. Body-settings reconciliation transfers it to automatically created gravity-derived rigid-body settings; no runtime gravity model receives a copy.
4. Explicit non-gravity rigid-body settings and the legacy explicit constant-mass input retain precedence. Gravity changes do not overwrite them.
5. Runtime gravity models have no inertia or center-of-mass accessor. Spherical-harmonic conversion is performed by the gravity-derived rigid-body properties; homogeneous polyhedron inertia moved there as part of the same migration. Ring and point-mass fields provide no inertia.
6. A gravity model stores a non-owning reverse pointer to the associated rigid-body properties. `Body` installs callbacks so changes to gravitational parameter update gravity-derived mass and inertia normalization, while spherical-harmonic coefficient changes update center of mass and inertia. A weak pointer is used deliberately to avoid a gravity/rigid-body ownership cycle.
7. Mean-moment estimation now reads and writes the gravity-derived rigid-body properties. The spherical-harmonic gravitational-torque partial also includes the gravitational-parameter contribution caused by the linked relation `M = mu/G`; this contribution is omitted when mass is explicitly configured and independent.
8. A spherical-harmonic gravity field remains valid without a finite scaled mean moment or inertia tensor. Ordinary instantaneous gravity variations update the current inertia tensor, but do not claim an inertia derivative because those variation models do not currently expose coefficient rates.
9. Direct C++ construction now has one explicit policy: construct the runtime gravity field without inertia configuration, then install `FromGravityFieldRigidBodyProperties(gravity_field, scaled_mean_moment)` on the same `Body` when gravity-derived inertia is required.
10. Python deprecation warnings are emitted for the legacy gravity-settings property and the old direct `GravityFieldModel(..., update_inertia_tensor=...)` callback argument. Both compatibility paths remain behaviorally functional.

## Validation outcome

- Full release build in `cmake-build-release`, using the `tudatpy-dev` environment and `-j8`: passed all 745 Ninja actions.
- Complete C++ CTest registry: 359 passed, 0 failed.
- Complete `tests/test_tudatpy` suite against the rebuilt extension: 510 passed, 12 skipped because optional external credentials/data were unavailable, and 1 existing expected-failure test unexpectedly passed.
- Formatting hooks and `git diff --check`: passed.
