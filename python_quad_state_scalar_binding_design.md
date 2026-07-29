# Runtime-selectable state scalar in TudatPy

## Goal

Keep the existing Python interface based on Python `float` and NumPy
`float64`, while allowing a process-wide choice between C++ state
computations instantiated with:

- `double`; or
- `tudat::HighPrecisionStateScalar`, configured as
  `boost::multiprecision::cpp_bin_float_quad`.

The selected C++ specialization performs the propagation, interpolation,
ephemeris, observation, and estimation work. Conversion to or from `double`
occurs only at the Python boundary.

## Recommended user-facing semantics

The selection should be made before importing `tudatpy.kernel` or any TudatPy
module that imports it:

```python
import tudatpy

tudatpy.set_state_scalar_type("quad")

from tudatpy.dynamics import propagation_setup, simulator
```

The default remains `"double"` for Python compatibility. Useful top-level
functions are:

```python
tudatpy.set_state_scalar_type("double" | "quad")
tudatpy.get_state_scalar_type()
tudatpy.quad_precision_available()
```

The choice should be locked after `tudatpy.kernel` is imported. A later
attempt to change it should raise `RuntimeError`. Hot switching is unsafe
because already-created Python objects contain mutually incompatible C++
template specializations.

## Binding initialization

Compile both relevant specializations into the existing `tudatpy.kernel`
extension, but register only one when the module initializes:

```cpp
template< typename StateScalarType >
void exposeStateScalarBindingsImpl( py::module_& module );

void exposeStateScalarBindings( py::module_& module )
{
    if( getRequestedPythonStateScalarType( ) == PythonStateScalarType::quad_precision )
    {
        exposeStateScalarBindingsImpl< tudat::HighPrecisionStateScalar >( module );
    }
    else
    {
        exposeStateScalarBindingsImpl< double >( module );
    }
}
```

Registering only the selected specialization preserves the existing Python
class names and avoids duplicate pybind11 type registrations. A single
extension is preferable to separate double and quad kernels because it
preserves canonical module names and avoids duplicating all scalar-independent
bindings.

The Tudat core used for such a Python build must be configured with:

```text
TUDAT_HIGH_PRECISION_STATE_SCALAR=CPP_BIN_FLOAT_QUAD
```

If it is not, requesting `"quad"` must fail clearly instead of silently using
`long double`.

## Python boundary conversions

Boost multiprecision Eigen matrices must not be exposed directly through
pybind11's NumPy caster. Public arguments and return values stay `double`.

For a NumPy state vector:

```cpp
[]( const Eigen::VectorXd& pythonState )
{
    return templatedFunction< StateScalarType >(
            pythonState.template cast< StateScalarType >( ) );
}
```

For a returned state:

```cpp
return internalState.template cast< double >( ).eval( );
```

Histories require conversion of every mapped state while preserving the time
key. Python callbacks require conversions in both directions.

This means Python cannot inject or retrieve digits beyond float64, but
arithmetic and storage between the two boundaries can remain quad precision.
Sub-double temporal resolution should continue to use `tudat::Time`.

## `propagator.translational`

This is where the initial NumPy state enters the propagation object graph. The
factory should accept `Eigen::VectorXd`, explicitly cast it, and construct the
selected settings specialization:

```cpp
template< typename StateScalarType >
void exposePropagatorSettingsImpl( py::module_& module )
{
    module.def(
            "translational",
            []( const std::vector< std::string >& centralBodies,
                const basic_astrodynamics::AccelerationMap& accelerationModels,
                const std::vector< std::string >& bodiesToIntegrate,
                const Eigen::VectorXd& initialStates,
                const tudat::Time& initialTime,
                const std::shared_ptr< IntegratorSettings< tudat::Time > >& integratorSettings,
                const std::shared_ptr< PropagationTerminationSettings >& terminationSettings,
                const TranslationalPropagatorType propagator,
                const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >& outputVariables,
                const std::shared_ptr< SingleArcPropagatorProcessingSettings >& processingSettings )
            {
                return translationalStatePropagatorSettings<
                        StateScalarType,
                        tudat::Time >(
                        centralBodies,
                        accelerationModels,
                        bodiesToIntegrate,
                        initialStates.template cast< StateScalarType >( ),
                        initialTime,
                        integratorSettings,
                        terminationSettings,
                        propagator,
                        outputVariables,
                        processingSettings );
            },
            /* retain the existing py::arg declarations and documentation */ );
}
```

`PropagatorSettings<StateScalarType>` and all derived settings classes must be
registered with the same selected scalar. Their `initial_states` Python
property should return `Eigen::VectorXd` and accept `Eigen::VectorXd`, with an
explicit cast when resetting the internal state.

Although the immediate use case is translational propagation, the rotational,
mass, custom-state, multi-type, multi-arc, and hybrid-arc factories must use
the same specialization to keep the module's settings type graph consistent.

## `create_dynamics_simulator`

This factory receives no raw Python floating-point state. Its scalar-bearing
argument is already a typed propagator-settings object. Therefore no settings
conversion should be attempted here:

```cpp
template< typename StateScalarType >
void exposeDynamicsSimulatorImpl( py::module_& module )
{
    module.def(
            "create_dynamics_simulator",
            []( const SystemOfBodies& bodies,
                const std::shared_ptr< PropagatorSettings< StateScalarType > >& propagatorSettings,
                const bool simulateDynamicsOnCreation )
            {
                return createDynamicsSimulator< StateScalarType, tudat::Time >(
                        bodies,
                        propagatorSettings,
                        simulateDynamicsOnCreation );
            },
            py::arg( "bodies" ),
            py::arg( "propagator_settings" ),
            py::arg( "simulate_dynamics_on_creation" ) = true,
            /* retain the existing documentation */ );
}
```

The typed function signature prevents accidentally passing double settings to
a quad simulator. The related `DynamicsSimulator`, `SingleArcSimulator`,
`MultiArcSimulator`, and `HybridArcSimulator` Python classes must be
registered from the same specialization.

The following simulator interfaces need explicit float64 boundary adapters:

- reintegration initial states;
- state-derivative functions;
- state-derivative modifier callbacks;
- processed and raw state histories;
- any direct state getter that otherwise returns an Eigen matrix containing
  the selected scalar.

Dependent-variable histories that are intrinsically `Eigen::VectorXd` do not
need state-scalar conversion.

## Connected binding graph

Changing only these two factories is not sufficient for a usable quad mode.
The same selection must propagate through all mutually interacting exposed
types, including:

- propagation settings and result classes;
- dynamics and variational simulators;
- environment state creation and ephemeris factories;
- state interpolators and state-history utilities;
- estimatable parameter sets and orbit-determination managers;
- observation simulators, observation sets, and observation collections;
- tracking-data containers.

The current `STATE_SCALAR_TYPE` macro occurs in 24 binding source files. These
files should use a C++ template parameter in their scalar-dependent binding
sections rather than testing the runtime mode repeatedly.

## Callback limitation

Python callbacks remain float64 interfaces. A callback invoked at every
derivative evaluation therefore downcasts the current quad state and upcasts
the returned derivative every time. Quad-precision behavior can only be
claimed end-to-end when the relevant dynamics and observation models stay
inside C++.

## Validation

Double and quad Python tests must run in separate processes because the mode
is locked when `tudatpy.kernel` loads. Tests should verify:

- default double behavior is unchanged;
- module paths, class names, signatures, and NumPy dtypes are identical;
- the selected scalar reports 53 or 113 binary digits as appropriate;
- float64 initial states create quad-specialized propagator settings in quad
  mode;
- `create_dynamics_simulator` receives the matching settings specialization;
- propagation, state retrieval, reset ephemerides, range, n-way range, and
  DSN Doppler use the selected scalar internally;
- attempts to change mode after kernel import fail;
- Python callbacks still function, with their precision limitation
  documented.

## Prototype findings

A local, uncommitted prototype templated the propagator-settings and dynamics
simulator exposure functions and syntax-checked both their double and quad
instantiations. Directly exposing Boost multiprecision Eigen matrices produced
the expected pybind11 error:

```text
Attempt to use a non-POD or unimplemented POD type as a numpy dtype
```

Explicit `Eigen::VectorXd`/`Eigen::MatrixXd` adapters eliminated that error for
the prototype. The prototype also confirmed that settings and simulator
classes have to be selected together; converting an existing
`PropagatorSettings<double>` object inside `create_dynamics_simulator` is not a
reasonable design.
