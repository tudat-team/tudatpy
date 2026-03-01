#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_simulator_variational_helpers.h"

#include <pybind11/pybind11.h>

#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tep = tudat::estimatable_parameters;

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

void bind_single_arc_variational_integrate_methods( SingleArcVariationalSimulatorPyClass& singleArcClass )
{
    singleArcClass
            .def( py::init< const tss::SystemOfBodies&,
                            const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< TIME_TYPE > >,
                            const std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >,
                            const std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                            const bool,
                            const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >,
                            const bool,
                            const bool,
                            const bool,
                            const bool >(),
                  py::arg( "bodies" ),
                  py::arg( "integrator_settings" ),
                  py::arg( "propagator_settings" ),
                  py::arg( "estimated_parameters" ),
                  py::arg( "integrate_equations_concurrently" ) = true,
                  py::arg( "variational_only_integrator_settings" ) =
                          std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< TIME_TYPE > >( ),
                  py::arg( "clear_numerical_solutions" ) = false,
                  py::arg( "integrate_on_creation" ) = true,
                  py::arg( "set_integrated_result" ) = false,
                  py::arg( "print_dependent_variable_data" ) = true,
                  R"doc(Create a single-arc variational simulator.)doc" )
            .def( "integrate_equations_of_motion_only",
                  &SingleArcVariationalSolverType::integrateDynamicalEquationsOfMotionOnly,
                  py::arg( "initial_states" ),
                  R"doc(Integrate only the equations of motion.)doc" )
            .def( "integrate_full_equations",
                  &SingleArcVariationalSolverType::integrateVariationalAndDynamicalEquations,
                  py::arg( "initial_states" ),
                  py::arg( "integrate_equations_concurrently" ) = true,
                  R"doc(Integrate equations of motion and variational equations.)doc" )
            .def_property_readonly(
                    "parameter_vector",
                    &SingleArcVariationalSolverType::getParametersToEstimate,
                    R"doc(Estimatable parameter set used in the variational equations.)doc" );
}

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy
