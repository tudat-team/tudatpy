# PR 913 implementation plan

This document records the staged implementation plan for the coupled orbital,
rotational, and tidal gravity model introduced by PR 913.

## Branch update policy

The feature branch must be updated by merging `develop` into
`feature/coupled_model`. Do not rebase this feature branch.

## 1. Baseline and scope cleanup

Status: completed on the isolated feature branch after merging the 2026-08-29
`origin/develop` head. No build or tests were run, by explicit instruction.

- Merge the latest `origin/develop` into `feature/coupled_model` and resolve the
  conflicts without rebasing.
- Remove temporary debug output introduced by the feature work.
- Restore the complete zero-proper-mode rotational initial-state algorithm,
  including its backward-propagation phase and original behavioral contract.
- Identify unrelated API or implementation changes and revert them or separate
  them from the coupled-model work.

## 2. Rigid-body ownership

Status: completed. Rigid-body properties own the current inertia tensor and
derivative, with explicit availability and an optional gravity-field provider.

- Make `RigidBodyProperties` the sole owner of the current inertia tensor and
  its derivative.
- Add explicit availability semantics for inertia and inertia derivatives; an
  unavailable tensor must not be represented by a zero matrix.
- Introduce an optional, one-way provider from rigid-body properties to the
  degree-two gravity-field data needed to calculate inertia.
- Keep gravity-only bodies valid when no inertia tensor is configured.

## 3. Integrated gravity variation

Status: completed together with the propagation lifecycle in points 4 and 5.

- Add an integrated gravity-field variation object to the existing gravity
  variation system.
- Store the current integrated coefficient variation in this object during
  propagation.
- Store a post-propagation history/interpolator in the same object.
- Select the current value or interpolated value through an
  `isBodyInPropagation` flag analogous to synchronous rotation.
- Compose this variation additively with the nominal field and all other
  gravity-field variations.

## 4. Propagation and postprocessing

Status: completed. Propagation and postprocessing now route coefficient
variations through the integrated variation object.

- Replace direct writes to a body's spherical-harmonic coefficient matrices
  with updates to the integrated variation object.
- Use typed five-component degree-two state slices for each propagated body.
- Install each body's complete coefficient-variation history after propagation.
- Forward propagation-state lifecycle changes to the integrated variation
  object.

## 5. Physical derivative models

Status: completed for the supported degree-two Maxwell propagation model.

- Connect the Maxwell deformation model to the integrated variation state.
- Validate supported degree, order, coefficient count, and component ordering.
- Define one normalization convention and perform conversions at one explicit
  boundary.
- Keep nominal/static coefficients out of `Body` and out of the propagated
  variation state.

## 6. Coupled derivative solution

Status: open. The generic coupled-state-derivative solver that was initially
implemented has been removed pending a focused design for this physical loop.

- Formulate the dependency between gravity-coefficient rates, inertia
  derivatives, and angular acceleration explicitly.
- Prefer a direct per-body coupled algebraic solve when the dependency is
  linear.
- If iteration is required, provide scaled convergence criteria, configurable
  limits, and a defined non-convergence policy.
- Keep state-derivative modifiers and variational-equation evaluation
  consistent with the final coupled derivative.

## 7. Python API, documentation, and validation

Status: partially completed for the supported single-arc degree-two Maxwell
feature. The generic coupled-solver settings and their Python API have been
removed; validation of the final point-6 design remains open.

- Expose high-level propagation and environment settings rather than internal
  body mutation functions.
- Document total-versus-variation semantics, normalization, inertia
  availability, interpolation behavior, and supported propagation modes.
- Add unit, coupled-dynamics, regression, multi-body, lifecycle, and Python
  end-to-end tests before merging the complete feature.

## Open follow-up: derivatives of prescribed gravity variations

Existing periodic, polynomial, tabulated, tidal, and custom gravity-field
variation models provide instantaneous coefficient corrections, but not their
time derivatives.  Their contribution to the gravity-linked inertia tensor is
included, while their contribution to the inertia-tensor derivative cannot yet
be evaluated.  Adding explicit coefficient-rate interfaces to those models is
left open; no numerical differentiation is introduced by this PR.
