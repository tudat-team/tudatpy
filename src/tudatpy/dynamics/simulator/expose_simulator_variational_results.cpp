#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_simulator_variational_helpers.h"

#include <pybind11/pybind11.h>

#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "kernelExternTemplates.h"

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

void bind_single_arc_variational_result_properties( SingleArcVariationalSimulatorPyClass& singleArcClass )
{
    singleArcClass
            .def_property_readonly(
                    "variational_equations_history",
                    &SingleArcVariationalSolverType::getNumericalVariationalEquationsSolution,
                    R"doc(Variational equations history [state transition, sensitivity].)doc" )
            .def_property_readonly(
                    "state_transition_matrix_history",
                    &SingleArcVariationalSolverType::getStateTransitionMatrixSolution,
                    R"doc(State transition matrix history.)doc" )
            .def_property_readonly(
                    "sensitivity_matrix_history",
                    &SingleArcVariationalSolverType::getSensitivityMatrixSolution,
                    R"doc(Sensitivity matrix history.)doc" )
            .def_property_readonly(
                    "state_history",
                    &SingleArcVariationalSolverType::getEquationsOfMotionSolutionDouble,
                    R"doc(Integrated state history.)doc" )
            .def_property_readonly(
                    "variational_propagation_results",
                    &SingleArcVariationalSolverType::getVariationalPropagationResults,
                    R"doc(Combined results object with dynamics and variational histories.)doc" )
            .def_property_readonly(
                    "dynamics_simulator",
                    &SingleArcVariationalSolverType::getDynamicsSimulator,
                    R"doc(Dynamics simulator used by this variational solver.)doc" );
}

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy
