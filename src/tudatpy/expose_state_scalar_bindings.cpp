#include <pybind11/pybind11.h>

#include "dynamics/environment/expose_environment_state_scalar.h"
#include "pythonStateScalar.h"

namespace py = pybind11;

#if TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE
#define TUDATPY_DISPATCH_STATE_SCALAR_BINDINGS( functionName, module )                               \
    do                                                                                               \
    {                                                                                                \
        if( tudatpy::getPythonStateScalarType( ) == tudatpy::PythonStateScalarType::quad_precision ) \
        {                                                                                            \
            quad_precision::functionName( module );                                                  \
        }                                                                                            \
        else                                                                                         \
        {                                                                                            \
            double_precision::functionName( module );                                                \
        }                                                                                            \
    } while( false )
#define TUDATPY_DECLARE_QUAD_STATE_SCALAR_BINDINGS( functionName ) \
    namespace quad_precision                                       \
    {                                                              \
    void functionName( py::module& module );                       \
    }
#else
#define TUDATPY_DISPATCH_STATE_SCALAR_BINDINGS( functionName, module ) double_precision::functionName( module )
#define TUDATPY_DECLARE_QUAD_STATE_SCALAR_BINDINGS( functionName )
#endif

#define TUDATPY_DECLARE_DOUBLE_STATE_SCALAR_BINDINGS( functionName ) \
    namespace double_precision                                       \
    {                                                                \
    void functionName( py::module& module );                         \
    }

#define TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( functionName ) \
    TUDATPY_DECLARE_DOUBLE_STATE_SCALAR_BINDINGS( functionName )        \
    TUDATPY_DECLARE_QUAD_STATE_SCALAR_BINDINGS( functionName )          \
    void functionName( py::module& module )                             \
    {                                                                   \
        TUDATPY_DISPATCH_STATE_SCALAR_BINDINGS( functionName, module ); \
    }

namespace tudatpy
{
namespace math
{
namespace interpolators
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_interpolators_state_scalar )

}  // namespace interpolators
}  // namespace math

namespace dynamics
{
namespace environment
{

namespace double_precision
{
void expose_environment_state_scalar( BodyPythonBinding& bodyBinding, SystemOfBodiesPythonBinding& systemOfBodiesBinding );
}
#if TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE
namespace quad_precision
{
void expose_environment_state_scalar( BodyPythonBinding& bodyBinding, SystemOfBodiesPythonBinding& systemOfBodiesBinding );
}
#endif

void expose_environment_state_scalar( BodyPythonBinding& bodyBinding, SystemOfBodiesPythonBinding& systemOfBodiesBinding )
{
#if TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE
    if( tudatpy::getPythonStateScalarType( ) == tudatpy::PythonStateScalarType::quad_precision )
    {
        quad_precision::expose_environment_state_scalar( bodyBinding, systemOfBodiesBinding );
        return;
    }
#endif
    double_precision::expose_environment_state_scalar( bodyBinding, systemOfBodiesBinding );
}

}  // namespace environment

namespace environment_setup
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_environment_setup_state_scalar )

namespace ephemeris
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_ephemeris_setup_state_scalar )

}  // namespace ephemeris
}  // namespace environment_setup

namespace parameters
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_parameters )

}  // namespace parameters

namespace parameters_setup
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_parameters_setup_state_scalar )

}  // namespace parameters_setup

namespace propagation
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_propagation_state_utility_bindings )
TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_propagation_results_bindings )

}  // namespace propagation

namespace propagation_setup
{
namespace propagator
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_propagator_setup )

}  // namespace propagator
}  // namespace propagation_setup

namespace simulator
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_simulator_dynamics_bindings )
TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_simulator_variational_bindings )

}  // namespace simulator
}  // namespace dynamics

namespace estimation
{
namespace estimation_analysis
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_estimation_analysis )
TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_estimation_analysis_orbit_determination_helpers )
TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_estimation_analysis_estimator )
TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_estimation_analysis_ephemeris_fit )

}  // namespace estimation_analysis

namespace observable_models
{
namespace observables_simulation
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_observables_simulation )

}  // namespace observables_simulation
}  // namespace observable_models

namespace observations
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_observations )

}  // namespace observations

namespace observations_setup
{
namespace observations_simulation_settings
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_observation_simulation_settings_factory_bindings )

}  // namespace observations_simulation_settings

namespace observations_wrapper
{

TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_observations_wrapper_io_bindings )
TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER( expose_observations_wrapper_simulation_bindings )

}  // namespace observations_wrapper
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy

#undef TUDATPY_DEFINE_STATE_SCALAR_BINDINGS_DISPATCHER
#undef TUDATPY_DECLARE_DOUBLE_STATE_SCALAR_BINDINGS
#undef TUDATPY_DECLARE_QUAD_STATE_SCALAR_BINDINGS
#undef TUDATPY_DISPATCH_STATE_SCALAR_BINDINGS
