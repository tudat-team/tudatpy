/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_estimation_setup.h"
#include "expose_estimation_setup/expose_estimated_parameter_setup.h"
#include "expose_estimation_setup/expose_observation_setup.h"

#include "tudat/astro/propagators/propagateCovariance.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace ts = tudat::statistics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tom = tudat::observation_models;


namespace tudatpy {
namespace simulation {
namespace estimation_setup {


void expose_estimation_setup(py::module &m) {

    //TODO: Remove variationalOnlyIntegratorSettings
    py::class_<
            tp::SingleArcVariationalEquationsSolver<double, double>,
            std::shared_ptr<tp::SingleArcVariationalEquationsSolver<double, double>>>(m, "SingleArcVariationalEquationsSolver")
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies&,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings<double>>,
                 const std::shared_ptr< tp::PropagatorSettings<double>>,
                 const std::shared_ptr< tep::EstimatableParameterSet< double > >,
                 const bool,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >,
                 const bool,
                 const bool,
                 const bool >(),
                 py::arg("bodies"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("estimated_parameters"),
                 py::arg("integrate_equations_concurrently") = true,
                 py::arg("variational_only_integrator_settings") = std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >( ),
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("integrate_on_creation") = true,
                 py::arg("set_integrated_result") = false )
            .def("integrate_equations_of_motion_only",
                 &tp::SingleArcVariationalEquationsSolver<double, double>::integrateDynamicalEquationsOfMotionOnly,
                 py::arg("initial_states"))
            .def("integrate_full_equations",
                 &tp::SingleArcVariationalEquationsSolver<double, double>::integrateVariationalAndDynamicalEquations,
                 py::arg("initial_states"),
                 py::arg("integrate_equations_concurrently"))
            .def_property("parameter_vector",
                          &tp::SingleArcVariationalEquationsSolver<double, double>::getParametersToEstimate,
                          &tp::SingleArcVariationalEquationsSolver<double, double>::resetParameterEstimate)
            .def_property_readonly("variational_equations_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getNumericalVariationalEquationsSolution)
            .def_property_readonly("state_transition_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getStateTransitionMatrixSolution)
            .def_property_readonly("sensitivity_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getSensitivityMatrixSolution)
            .def_property_readonly("state_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getEquationsOfMotionSolution)
            .def_property_readonly("dynamics_simulator",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getDynamicsSimulator);

    py::class_<
            tp::CombinedStateTransitionAndSensitivityMatrixInterface,
            std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>>(
                m, "CombinedStateTransitionAndSensitivityMatrixInterface")
            .def("state_transition_sensitivity_at_epoch",
                 &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                 getCombinedStateTransitionAndSensitivityMatrix,
                 py::arg("time"))
            .def("full_state_transition_sensitivity_at_epoch",
                 &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                 getFullCombinedStateTransitionAndSensitivityMatrix,
                 py::arg("time"))
            .def_property_readonly(
                "state_transition_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getStateTransitionMatrixSize)
            .def_property_readonly(
                "sensitivity_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getSensitivityMatrixSize)
            .def_property_readonly(
                "full_parameter_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getFullParameterVectorSize);;

    m.def("propagate_covariance",
          py::overload_cast<
          const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >(
              &tp::propagateCovariance ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times") );

    m.def("propagate_formal_errors",
          py::overload_cast< const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >( &tp::propagateFormalErrors ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times") );

    m.def("propagate_covariance",
          py::overload_cast< const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >( &tp::propagateCovariance ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times") );

    py::class_<tep::EstimatableParameterSet<double>,
            std::shared_ptr<tep::EstimatableParameterSet<double>>>(m, "EstimatableParameterSet")
            .def_property_readonly( "parameter_set_size",
                                    &tep::EstimatableParameterSet<double>::getEstimatedParameterSetSize )
            .def_property_readonly( "initial_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalStateParameterSize )
            .def_property_readonly( "initial_multi_arc_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalSingleArcStateParameterSize )
            .def_property_readonly( "initial_multi_arc_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalMultiArcStateParameterSize )
            .def_property_readonly( "constraints_size",
                                    &tep::EstimatableParameterSet<double>::getConstraintSize )
            .def_property_readonly( "values",
                                    &tep::EstimatableParameterSet<double>::getFullParameterValues< double > )
            .def_property( "parameter_vector",
                           &tep::EstimatableParameterSet<double>::getFullParameterValues< double >,
                           &tep::EstimatableParameterSet<double>::resetParameterValues< double > )
            .def( "indices_for_parameter_type",
                  &tep::EstimatableParameterSet<double>::getIndicesForParameterType,
                  py::arg("parameter_type") );

    py::class_<
            tss::EstimationConvergenceChecker,
            std::shared_ptr<tss::EstimationConvergenceChecker>>(m, "EstimationConvergenceChecker")
            .def(py::init< const unsigned int,
                 const double,
                 const double,
                 const int >( ),
                 py::arg( "maximum_iterations" ) = 5,
                 py::arg( "minimum_residual_change" ) = 0.0,
                 py::arg( "minimum_residual" ) = 0.0,
                 py::arg( "number_of_iterations_without_improvement" ) = 2 );

    py::class_<
            tss::PodInput<double, double>,
            std::shared_ptr<tss::PodInput<double, double>>>(m, "PodInput")
            .def(py::init<
                 const std::shared_ptr< tom::ObservationCollection< > >&,
                 const int,
                 const Eigen::MatrixXd,
                 const Eigen::VectorXd >( ),
                 py::arg( "observations_and_times" ),
                 py::arg( "parameter_size" ),
                 py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 py::arg( "apriori_parameter_correction" ) = Eigen::VectorXd( 0 ) )
            .def( "set_constant_weight",
                  &tss::PodInput<double, double>::setConstantWeightsMatrix,
                  py::arg( "weight" ) )
            .def( "set_constant_weight_per_observable",
                  &tss::PodInput<double, double>::setConstantPerObservableWeightsMatrix,
                  py::arg( "weight_per_observable" ) )
            .def( "set_constant_weight_per_observable_and_link_end",
                  &tss::PodInput<double, double>::setConstantPerObservableAndLinkEndsWeights,
                  py::arg( "weight_per_observable_and_link" ) )
            .def( "define_estimation_settings",
                  &tss::PodInput<double, double>::defineEstimationSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "save_residuals_and_parameters_per_iteration" ) = true,
                  py::arg( "save_state_history_per_iteration" ) = false );

    py::class_<
            tss::PodOutput<double, double>,
            std::shared_ptr<tss::PodOutput<double, double>>>(m, "PodOutput")
            .def_property_readonly("inverse_covariance",
                                   &tss::PodOutput<double, double>::getUnnormalizedInverseCovarianceMatrix)
            .def_property_readonly("covariance",
                                   &tss::PodOutput<double, double>::getUnnormalizedCovarianceMatrix)
            .def_property_readonly("formal_errors",
                                   &tss::PodOutput<double, double>::getFormalErrorVector)
            .def_property_readonly("correlations",
                                   &tss::PodOutput<double, double>::getCorrelationMatrix)
            .def_property_readonly("residual_history",
                                   &tss::PodOutput<double, double>::getResidualHistoryMatrix)
            .def_property_readonly("parameter_history",
                                   &tss::PodOutput<double, double>::getParameterHistoryMatrix)
            .def_property_readonly("design_matrix",
                                   &tss::PodOutput<double, double>::getUnnormalizedInformationMatrix)
            .def_property_readonly("normalized_design_matrix",
                                   &tss::PodOutput<double, double>::getNormalizedInformationMatrix)
            .def_property_readonly("weighted_design_matrix",
                                   &tss::PodOutput<double, double>::getUnnormalizedWeightedInformationMatrix)
            .def_property_readonly("weighted_normalized_design_matrix",
                                   &tss::PodOutput<double, double>::getNormalizedWeightedInformationMatrix);


    py::class_<
            tss::OrbitDeterminationManager<double, double>,
            std::shared_ptr<tss::OrbitDeterminationManager<double, double>>>(m, "OrbitDeterminationManager")
            .def(py::init<const tss::SystemOfBodies&,
                 const std::shared_ptr< tep::EstimatableParameterSet< double > >,
                 const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&,
                 const std::shared_ptr< tni::IntegratorSettings< double > >,
                 const std::shared_ptr< tp::PropagatorSettings< double > >,
                 const bool >( ),
                 py::arg("bodies"),
                 py::arg("estimated_parameters"),
                 py::arg("observation_settings"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("integrate_on_creation") = true  )
            .def_property_readonly("observation_simulators",
                                   &tss::OrbitDeterminationManager<double, double>::getObservationSimulators)
            .def_property_readonly("observation_managers",
                                   &tss::OrbitDeterminationManager<double, double>::getObservationManagers)
            .def_property_readonly("state_transition_interface",
                                   &tss::OrbitDeterminationManager<double, double>::getStateTransitionAndSensitivityMatrixInterface)

            .def("perform_estimation",
                 &tss::OrbitDeterminationManager<double, double>::estimateParameters,
                 py::arg( "estimation_input" ),
                 py::arg( "convergence_checker" ) = std::make_shared< tss::EstimationConvergenceChecker >( ) );


    m.def("create_parameters_to_estimate",
          &tss::createParametersToEstimate< double >,
          py::arg("parameter_settings"),
          py::arg("bodies"),
          py::arg("propagator_settings") =
            std::shared_ptr< tp::PropagatorSettings< double > >( ) );

    auto parameter_setup = m.def_submodule("parameter");
    parameter::expose_estimated_parameter_setup(parameter_setup);

    auto observation_setup = m.def_submodule("observation");
    observation::expose_observation_setup(observation_setup);
}

}
}
}// namespace tudatpy
