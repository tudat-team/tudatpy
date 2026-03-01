#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_simulator_variational_helpers.h"

#include <pybind11/pybind11.h>

#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "kernelExternTemplates.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

void bind_create_variational_equations_solver_function( py::module& m )
{
    m.def( "create_variational_equations_solver",
           &tss::createVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "propagator_settings" ),
           py::arg( "parameters_to_estimate" ),
           py::arg( "simulate_dynamics_on_creation" ) = true,
           R"doc(Create a variational equations solver for the provided setup.)doc" );
}

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy
