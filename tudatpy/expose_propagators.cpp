//
// Created by ggarrett on 26-04-20.
//

#include <pybind11/pybind11.h>
#include <pybind11/chrono.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "expose_propagators.h"

//Did you forget to ``? Or <pybind11/complex.h>,
//<pybind11/functional.h>, <pybind11/chrono.h>, etc. Some automatic
//conversions are optional and require extra headers to be included
//        when compiling your pybind11 module.

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;

namespace tudatpy {
    void expose_propagators(py::module &m) {

        // Astrodynamics/Propagators/singeStateTypeDerivative.h
        m.def("get_single_integration_size",
              &tp::getSingleIntegrationSize,
              py::arg("state_type"));
        m.def("get_single_integration_differential_equation_order",
              &tp::getSingleIntegrationDifferentialEquationOrder,
              py::arg("state_type"));
        m.def("get_generalized_acceleration_size",
              &tp::getGeneralizedAccelerationSize,
              py::arg("state_type"));

        py::class_<
                tp::SingleStateTypeDerivative<double, double>,
                std::shared_ptr<tp::SingleStateTypeDerivative<double, double >>
        > SingleStateTypeDerivative_(m,
                                     "SingleStateTypeDerivative");
//                    .def(py::init<const tp::IntegratedStateType>(),
//                         py::arg("integrated_state_type"));

        py::class_<
                tp::NBodyStateDerivative<double, double>,
                std::shared_ptr<tp::NBodyStateDerivative<double, double >>
        > NBodyStateDerivative_(m, "NBodyStateDerivative");

        py::class_<
                tp::NBodyCowellStateDerivative<double, double>,
                std::shared_ptr<tp::NBodyCowellStateDerivative<double, double >>
        >(m, "NBodyCowellStateDerivative")
                .def(py::init<
                             const tudat::basic_astrodynamics::AccelerationMap &,
                             const std::shared_ptr<tp::CentralBodyData<double, double>>,
                             const std::vector<std::string> &
                     >(),
                     py::arg("acceleration_models_per_body"),
                     py::arg("central_body_data"),
                     py::arg("bodies_to_integrate"));

        py::class_<
                tp::SingleArcDynamicsSimulator<double, double>,
                std::shared_ptr<tp::SingleArcDynamicsSimulator<double, double >>>(m, "SingleArcDynamicsSimulator")
                .def(py::init<
                             const tudat::simulation_setup::NamedBodyMap &,
                             const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double>>,
                             const std::shared_ptr<tp::PropagatorSettings<double>>,
                             const bool,
                             const bool,
                             const bool,
                             const bool,
                             const std::chrono::steady_clock::time_point,
                             const std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>> &
                     >(),
                     py::arg("body_map"),
                     py::arg("integrator_settings"),
                     py::arg("propagator_settings"),
                     py::arg("are_equations_of_motion_to_be_integrated") = true,
                     py::arg("clear_numerical_solutions") = false,
                     py::arg("set_integrated_result") = false,
                     py::arg("print_number_of_function_evaluations") = false,
                     py::arg("initial_clock_time") = std::chrono::steady_clock::now(),
                     py::arg("state_derivative_models") =
                             std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double> > >()
                )
                .def("integrate_equations_of_motion",
                     &tp::SingleArcDynamicsSimulator<double, double>::integrateEquationsOfMotion,
                     py::arg("initial_states"))
                .def("get_equations_of_motion_numerical_solution",
                     &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolution)
                .def("get_equations_of_motion_numerical_solution_raw",
                     &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionRaw)
                .def("get_dependent_variable_history",
                     &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableHistory)
                .def("get_cumulative_computation_time_history",
                     &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeComputationTimeHistory)
                .def("get_cumulative_number_of_function_evaluations",
                     &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeNumberOfFunctionEvaluations)
                .def("get_equations_of_motion_numerical_solution_base",
                     &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionBase)
                .def("get_dependent_variable_numerical_solution_base",
                     &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableNumericalSolutionBase)
                .def("get_cumulative_computation_time_history_base",
                     &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeComputationTimeHistoryBase)
                .def("manually_set_and_process_raw_numerical_equations_of_motion_solution",
                     &tp::SingleArcDynamicsSimulator<double, double>::manuallySetAndProcessRawNumericalEquationsOfMotionSolution,
                     py::arg("equations_of_motion_numerical_solution"),
                     py::arg("dependent_variable_history"),
                     py::arg("process_solution"))
                .def("get_integrator_settings",
                     &tp::SingleArcDynamicsSimulator<double, double>::getIntegratorSettings)
                .def("get_state_derivative_function",
                     &tp::SingleArcDynamicsSimulator<double, double>::getStateDerivativeFunction)
                .def("get_double_state_derivative_function",
                     &tp::SingleArcDynamicsSimulator<double, double>::getDoubleStateDerivativeFunction)
                .def("get_environment_updater",
                     &tp::SingleArcDynamicsSimulator<double, double>::getEnvironmentUpdater)
                .def("get_dynamics_state_derivative",
                     &tp::SingleArcDynamicsSimulator<double, double>::getDynamicsStateDerivative)
                .def("get_propagation_termination_condition",
                     &tp::SingleArcDynamicsSimulator<double, double>::getPropagationTerminationCondition)
                .def("get_integrated_state_processors",
                     &tp::SingleArcDynamicsSimulator<double, double>::getIntegratedStateProcessors)
                .def("get_propagation_termination_reason",
                     &tp::SingleArcDynamicsSimulator<double, double>::getPropagationTerminationReason)
                .def("integration_completed_successfully",
                     &tp::SingleArcDynamicsSimulator<double, double>::integrationCompletedSuccessfully)
                .def("get_dependent_variable_ids",
                     &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableIds)
                .def("get_initial_propagation_time",
                     &tp::SingleArcDynamicsSimulator<double, double>::getInitialPropagationTime)
                .def("reset_initial_propagation_time",
                     &tp::SingleArcDynamicsSimulator<double, double>::resetInitialPropagationTime)
                .def("get_dependent_variables_functions",
                     &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariablesFunctions)
                .def("reset_propagation_termination_conditions",
                     &tp::SingleArcDynamicsSimulator<double, double>::resetPropagationTerminationConditions)
                .def("process_numerical_equations_of_motion_solution",
                     &tp::SingleArcDynamicsSimulator<double, double>::processNumericalEquationsOfMotionSolution);

        py::enum_<tp::TranslationalPropagatorType>(m, "TranslationalPropagatorType")
                .value("undefined_translational_propagator",
                       tp::TranslationalPropagatorType::undefined_translational_propagator)
                .value("cowell",
                       tp::TranslationalPropagatorType::cowell)
                .value("encke",
                       tp::TranslationalPropagatorType::encke)
                .value("gauss_keplerian",
                       tp::TranslationalPropagatorType::gauss_keplerian)
                .value("gauss_modified_equinoctial",
                       tp::TranslationalPropagatorType::gauss_modified_equinoctial)
                .value("unified_state_model_quaternions",
                       tp::TranslationalPropagatorType::unified_state_model_quaternions)
                .value("unified_state_model_modified_rodrigues_parameters",
                       tp::TranslationalPropagatorType::unified_state_model_modified_rodrigues_parameters)
                .value("unified_state_model_exponential_map",
                       tp::unified_state_model_exponential_map)
                .export_values();

        py::class_<tp::DependentVariableSaveSettings,
                std::shared_ptr<tp::DependentVariableSaveSettings >>(m, "DependentVariableSaveSettings")
                .def(py::init<
                             const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings>>,
                             const bool
                     >(),
                     py::arg("dependent_variables"),
                     py::arg("print_dependent_variable_types") = true);

        py::class_<
                tp::PropagatorSettings<double>,
                std::shared_ptr<tp::PropagatorSettings<double >>
        > PropagatorSettings_(m, "PropagatorSettings");

        py::class_<
                tp::SingleArcPropagatorSettings<double>,
                std::shared_ptr<tp::SingleArcPropagatorSettings<double >>,
                tp::PropagatorSettings<double>
        > SingleArcPropagatorSettings_(m, "SingleArcPropagatorSettings");

        py::class_<
                tp::TranslationalStatePropagatorSettings<double>,
                std::shared_ptr<tp::TranslationalStatePropagatorSettings<double >>,
                tp::SingleArcPropagatorSettings<double>
        >(m, "TranslationalStatePropagatorSettings")
                .def( // ctor 1
                        py::init<
                                const std::vector<std::string> &,
                                const tba::AccelerationMap &,
                                const std::vector<std::string> &,
                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                                const std::shared_ptr<tp::PropagationTerminationSettings>,
                                const tp::TranslationalPropagatorType,
                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
                                const double>(),
                        py::arg("central_bodies"),
                        py::arg("accelerations_map"),
                        py::arg("bodies_to_integrate"),
                        py::arg("initial_body_states"),
                        py::arg("termination_settings"),
                        py::arg("propagator") = tp::TranslationalPropagatorType::cowell,
                        py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                        py::arg("print_interval") = TUDAT_NAN)
                .def( // ctor 2
                        py::init<const std::vector<std::string> &,
                                const tss::SelectedAccelerationMap &,
                                const std::vector<std::string> &,
                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                                const std::shared_ptr<tp::PropagationTerminationSettings>,
                                const tp::TranslationalPropagatorType,
                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
                                const double>(),
                        py::arg("central_bodies"),
                        py::arg("acceleration_settings_map"),
                        py::arg("bodies_to_integrate"),
                        py::arg("initial_body_states"),
                        py::arg("termination_settings"),
                        py::arg("propagator") = tp::cowell,
                        py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                        py::arg("print_interval") = TUDAT_NAN)
                .def( // ctor 3
                        py::init<const std::vector<std::string> &,
                                const tba::AccelerationMap &,
                                const std::vector<std::string> &,
                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                                const double,
                                const tp::TranslationalPropagatorType,
                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
                                const double>(),
                        py::arg("central_bodies"),
                        py::arg("accelerations_map"),
                        py::arg("bodies_to_integrate"),
                        py::arg("initial_body_states"),
                        py::arg("end_time"),
                        py::arg("propagator") = tp::cowell,
                        py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                        py::arg("print_interval") = TUDAT_NAN)
                .def( // ctor 4
                        py::init<const std::vector<std::string> &,
                                const tss::SelectedAccelerationMap &,
                                const std::vector<std::string> &,
                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                                const double,
                                const tp::TranslationalPropagatorType,
                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
                                const double>(),
                        py::arg("central_bodies"),
                        py::arg("acceleration_settings_map"),
                        py::arg("bodies_to_integrate"),
                        py::arg("initial_body_states"),
                        py::arg("end_time"),
                        py::arg("propagator") = tp::cowell,
                        py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                        py::arg("print_interval") = TUDAT_NAN);
    }
}
