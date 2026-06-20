# Full Two-Body Variational Debug State

Date: 2026-06-20

This is a temporary debugging record for the Mars-Phobos full two-body variational-equations work.

## Current Scope

The active investigation is focused on the restricted equivalence case:

- legacy mutual spherical-harmonic acceleration + second-order gravitational torque
- restricted full two-body acceleration/torque using `getExtendedSinglePointMassInteractions( 2, 2, 2, 2 )`

These should be physically equivalent for the single-point interaction terms. The test file currently contains temporary diagnostic paths and has the other full two-body variational tests disabled with `#if 0` so the restricted comparison can be isolated.

## Confirmed Fix Already Present

`fullTwoBodySphericalHarmonicGravityPartial.cpp` contains the sign change in the Body1 relative-quaternion chain:

```cpp
-detail::getRightQuaternionMultiplicationMatrix( conjugatedQuaternionVectorOfBody2 ) *
 partialOfInertialToBodyQuaternionWrtBodyToInertialQuaternion;
```

This fixed `test_acceleration_partials_MutualExtendedBodySphericalHarmonicPartials` in the earlier debugging session. The analogous torque sign flip was tested earlier and rejected because it broke the full two-body torque partial tests.

## Current Test Diagnostics

The file `unitTestFullTwoBodyVariationalEquations.cpp` currently includes:

- no-propagation initial partial diagnostics using `SingleArcVariationalEquationsSolver(..., integrateEquationsOnCreation=false)`
- explicit `DynamicsStateDerivativeModel` evaluation at the initial state/time
- `simulation_setup::setAreBodiesInPropagation( setup.bodies, true )` around the derivative/partial extraction
- direct extraction of the full two-body acceleration partial wrt Phobos rotation and position
- split full-two-body interaction diagnostics:
  - `fullTwoBodyPhobosFigureSinglePointTerms`
  - `fullTwoBodyMarsFigureSinglePointTerms`
- a one-step normal propagation comparison that extracts final raw state, STM, sensitivities, and final RHS without creating state interpolators

The restricted full two-body torque now uses the same interaction list as the restricted acceleration:

```cpp
fullTwoBodySphericalHarmonicGravitationalTorque(
        getExtendedSinglePointMassInteractions( 2, 2, 2, 2 ) )
```

## Build/Test Command Used

Before rebuilds, ongoing C++ compilation was killed as requested. The target was built with:

```bash
/snap/clion/465/bin/ninja/linux/x64/ninja test_propagators_FullTwoBodyVariationalEquations -j6
```

from:

```bash
/home/dominic/Tudat/tudat-monorepo/tudatpy/cmake-build-release
```

The current diagnostic target builds successfully.

## Main Finding

At the initial state, with no propagation:

- acceleration values match between restricted full two-body and legacy mutual
- torque values match between restricted full two-body and legacy mutual
- full-two-body variational assembly equals the direct full-two-body acceleration partial
- full-two-body acceleration partial wrt Phobos position matches the mutual variational RHS position block
- the dominant mismatch is isolated to the full-two-body acceleration partial wrt Phobos quaternion

Representative initial-state no-propagation result for `d a / d q0`:

```text
restricted full two-body:
[ 0.50940124612360516, -0.0003277227830680296, 0.0012727106417694811 ]

legacy mutual:
[-3.0958541528919266e-07, -1.2745715784193753e-07, -3.3398995465139373e-07]

absolute error: 0.50940325177613677
relative error: 1.0771821782502872e6
```

The full-two-body direct partial is inserted correctly into the variational RHS. This is not a propagation problem and not an assembly problem.

## Split Interaction Diagnostic

The restricted full two-body interaction list was split into:

- Phobos-figure terms: keep `(l1,m1,0,0)`, remove `(0,0,l2,m2)` and `(0,0,0,0)`
- Mars-figure terms: keep `(0,0,l2,m2)`, remove `(l1,m1,0,0)` and `(0,0,0,0)`

For `d a / d q0`, the split sum does not reproduce the combined restricted full-two-body partial:

```text
Phobos-figure-only + Mars-figure-only:
[-0.00068986453081826271, -0.00032772278306861214, 0.0012727106417694551]

combined restricted full two-body:
[ 0.50940124612360516, -0.0003277227830680296, 0.0012727106417694811]
```

The split sum differs from the combined value almost entirely in the first component.

The Mars-figure-only partial wrt Phobos rotation is nonzero, even though this dependency is physically suspect:

```text
Mars-figure-only d a / d q0:
[-0.00068955494540297356, -0.00032759532591077018, 0.0012730446317241065]
```

This points to an error inside the full-two-body acceleration orientation partial calculation when both restricted single-point interaction groups are present. The most likely area is `FullTwoBodySphericalHarmonicsGravityPartial::updateCurrentOrientationPartials`, specifically the handling of effective coefficients / relative coefficient rotation / body orientation chain rules across mixed `(l1,m1,0,0)` and `(0,0,l2,m2)` terms.

## Acceleration vs Torque Isolation

Two temporary variants were tested:

1. legacy mutual acceleration + restricted full-two-body torque
   - large translational mismatches disappeared
   - only smaller rotational-row mismatches remained

2. restricted full-two-body acceleration + legacy second-order torque
   - large translational mismatches returned immediately

Conclusion: the dominant `r,v` mismatch is in the acceleration partial, not the torque partial.

Another temporary diagnostic removed full-two-body torque terms with body-undergoing-torque degree/order `(0,0)`. This did not change the large mismatch.

## Current Failure Pattern

The no-propagation diagnostic fails only in the acceleration partial wrt Phobos rotation comparison and split-consistency checks.

The one-step propagated comparison still shows large mismatches in:

- `d(r,v)/d(q)`
- `d(r,v)/d(omega)`
- selected translational sensitivity rows

But these are downstream of the initial acceleration partial mismatch.

## Important Caveats

The current commit is intentionally temporary and diagnostic-heavy:

- it adds temporary enum values for split interaction cases
- it adds a no-propagation diagnostic test path
- it keeps unrelated tests in this file disabled with `#if 0`
- it keeps the one-step comparison and RHS extraction diagnostics

Before final PR cleanup, these diagnostics should either be removed or converted into a concise permanent regression test.
