#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_estimation_analysis_estimator_bindings.h"

#include <pybind11/pybind11.h>

#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tp = tudat::propagators;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

void bind_estimator_constructor_and_properties( EstimatorPyClass& estimatorClass )
{
    estimatorClass
            .def( py::init< const tss::SystemOfBodies&,
                            const std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                            const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&,
                            const std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >,
                            const bool >(),
                  py::arg( "bodies" ),
                  py::arg( "estimated_parameters" ),
                  py::arg( "observation_settings" ),
                  py::arg( "propagator_settings" ),
                  py::arg( "integrate_on_creation" ) = true,
                  R"doc(Create an estimator from environment, parameters, observation settings, and propagation settings.)doc" )
            .def_property_readonly(
                    "observation_simulators",
                    &EstimatorType::getObservationSimulators,
                    R"doc(Observation simulators owned by this estimator.)doc" )
            .def_property_readonly(
                    "observation_managers",
                    &EstimatorType::getObservationManagers,
                    R"doc(Observation managers owned by this estimator.)doc" )
            .def_property_readonly(
                    "state_transition_interface",
                    &EstimatorType::getStateTransitionAndSensitivityMatrixInterface,
                    R"doc(State-transition/sensitivity interface used during estimation.)doc" )
            .def_property_readonly(
                    "variational_solver",
                    &EstimatorType::getVariationalEquationsSolver,
                    R"doc(Variational solver used by this estimator.)doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
