/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_estimation_analysis.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"
#include "tudat/astro/propagators/propagateCovariance.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tp = tudat::propagators;
namespace trf = tudat::reference_frames;


namespace tudat
{

namespace propagators
{

std::map< double, Eigen::MatrixXd > propagateCovarianceRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE > >
                estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > >
                orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    tp::propagateCovariance(
            propagatedCovariance,
            estimationOutput->getUnnormalizedCovarianceMatrix( ),
            orbitDeterminationManager->getStateTransitionAndSensitivityMatrixInterface( ),
            evaluationTimes );

    tss::SystemOfBodies bodies = orbitDeterminationManager->getBodies( );

    std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > > parameterSet =
            orbitDeterminationManager->getParametersToEstimate( );

    std::map< int,
              std::shared_ptr< tep::EstimatableParameter<
                      Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > > > >
            initialStates = parameterSet->getInitialStateParameters( );
    std::map< std::pair< std::string, std::string >, std::vector< int > > transformationList;
    for( auto it: initialStates )
    {
        if( std::dynamic_pointer_cast<
                    tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > >( it.second ) )
        {
            std::shared_ptr< tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > >
                    currentInitialState = std::dynamic_pointer_cast<
                            tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > >(
                            it.second );
            transformationList[ std::make_pair(
                                        currentInitialState->getParameterName( ).second.first,
                                        currentInitialState->getCentralBody( ) ) ]
                    .push_back( it.first );
        }
        else if( std::dynamic_pointer_cast<
                         tep::ArcWiseInitialTranslationalStateParameter< STATE_SCALAR_TYPE > >(
                         it.second ) )
        {
            throw std::runtime_error(
                    "Error, multi-arc not yet supported in automatic "
                    "covariance conversion" );
        }
    }

    Eigen::Matrix3d currentInertialToRswPosition;
    Eigen::Matrix6d currentInertialToRswState;
    Eigen::MatrixXd currentFullInertialToRswState = Eigen::MatrixXd::Zero( 6, 6 );

    std::map< double, Eigen::MatrixXd > propagatedRswCovariance;
    for( auto it: propagatedCovariance )
    {
        double currentTime = static_cast< double >( it.first );
        Eigen::MatrixXd currentCovariance = it.second;
        currentFullInertialToRswState.setZero( );

        for( auto it_body: transformationList )
        {
            Eigen::Vector6d relativeState =
                    bodies.getBody( it_body.first.first )
                            ->getStateInBaseFrameFromEphemeris( currentTime ) -
                    bodies.getBody( it_body.first.second )
                            ->getStateInBaseFrameFromEphemeris( currentTime );
            currentInertialToRswPosition =
                    trf::getInertialToRswSatelliteCenteredFrameRotationMatrix( relativeState );
            currentInertialToRswState.block( 0, 0, 3, 3 ) = currentInertialToRswPosition;
            currentInertialToRswState.block( 3, 3, 3, 3 ) = currentInertialToRswPosition;
            for( unsigned int j = 0; j < it_body.second.size( ); j++ )
            {
                int currentStartIndex = it_body.second.at( j );
                currentFullInertialToRswState.block( currentStartIndex, currentStartIndex, 6, 6 ) =
                        currentInertialToRswState;
            }
        }
        propagatedRswCovariance[ currentTime ] = currentFullInertialToRswState * currentCovariance *
                currentFullInertialToRswState.transpose( );
    }
    return propagatedRswCovariance;
}

std::pair< std::vector< double >, std::vector< Eigen::MatrixXd > > propagateCovarianceVectorsRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE > >
                estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > >
                orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedRswCovariance =
            propagateCovarianceRsw( estimationOutput, orbitDeterminationManager, evaluationTimes );

    return std::make_pair( utilities::createVectorFromMapKeys( propagatedRswCovariance ),
                           utilities::createVectorFromMapValues( propagatedRswCovariance ) );
}

std::map< double, Eigen::VectorXd > propagateFormalErrorsRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE > >
                estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > >
                orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    std::map< double, Eigen::VectorXd > propagatedFormalErrors;

    propagatedCovariance =
            propagateCovarianceRsw( estimationOutput, orbitDeterminationManager, evaluationTimes );
    tp::convertCovarianceHistoryToFormalErrorHistory( propagatedFormalErrors,
                                                      propagatedCovariance );

    return propagatedFormalErrors;
}

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > propagateFormalErrorVectorsRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE > >
                estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > >
                orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::VectorXd > propagatedFormalErrors = propagateFormalErrorsRsw(
            estimationOutput, orbitDeterminationManager, evaluationTimes );
    tp::propagateFormalErrorsRsw( estimationOutput, orbitDeterminationManager, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys( propagatedFormalErrors ),
                           utilities::createVectorFromMapValues( propagatedFormalErrors ) );
}

std::pair< std::vector< double >, std::vector< Eigen::MatrixXd > > propagateCovarianceVectors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >
                stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    tp::propagateCovariance(
            propagatedCovariance, initialCovariance, stateTransitionInterface, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys( propagatedCovariance ),
                           utilities::createVectorFromMapValues( propagatedCovariance ) );
}

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > propagateFormalErrorVectors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >
                stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::VectorXd > propagatedFormalErrors;
    tp::propagateFormalErrors(
            propagatedFormalErrors, initialCovariance, stateTransitionInterface, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys( propagatedFormalErrors ),
                           utilities::createVectorFromMapValues( propagatedFormalErrors ) );
}

}  // namespace propagators

}  // namespace tudat


namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

void expose_estimation_analysis( py::module& m )
{

    // ************** Modules ***************

    // TO BE MODIFIED ACCORDINGLY - this was for the observable_models_setup module
    // auto biases = m.def_submodule( "biases" );
    // biases::expose_biases( biases );


    py::class_< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "Estimator",
            R"doc(

         Class for consolidating all estimation functionality.

         Class for consolidating all functionality required to perform an estimation.





      )doc" )
            .def( py::init< const tss::SystemOfBodies &,
                            const std::shared_ptr<
                                    tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                            const std::vector< std::shared_ptr< tom::ObservationModelSettings > > &,
                            const std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >,
                            const bool >( ),
                  py::arg( "bodies" ),
                  py::arg( "estimated_parameters" ),
                  py::arg( "observation_settings" ),
                  py::arg( "propagator_settings" ),
                  py::arg( "integrate_on_creation" ) = true,
                  R"doc(

         Class constructor.

         Constructor through which the user can create instances of this class.
         Defines environment, propagation and integrations models, as well as a number of settings related
         to the estimatable parameters and observation settings.

         .. note:: When using default settings, creating an object of
                   this type automatically triggers the propagation


         Parameters
         ----------
         bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
             Object defining the physical environment, with all
             properties of artificial and natural bodies.

         estimated_parameters : :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
             Object defining a consolidated set of estimatable parameters,
             linked to the environment and acceleration settings of the simulation.

         observation_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`
             List of settings objects, each object defining the observation model settings for one
             combination of observable and link geometry that is to be simulated.

         integrator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorSettings`
             Settings to create the numerical integrator that is to be
             used for the integration of the equations of motion

         propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
             Settings to create the propagator that is to be
             used for the propagation of dynamics

         integrate_on_creation : bool, default = True
             Boolean defining whether the propagation should be
             performed immediately (default), or at a later time
             (when calling the :func:`perform_estimation` member function.





     )doc" )
            .def_property_readonly(
                    "observation_simulators",
                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getObservationSimulators,
                    R"doc(

         **read-only**

         Observation simulators contained in the Estimator object. A single observation simulator hosts
         the functionality for simulating a given observable over the defined link geometry.


         :type: list[ :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` ]
      )doc" )
            .def_property_readonly(
                    "observation_managers",
                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getObservationManagers,
                    R"doc(

         **read-only**

         Observation managers contained in the Estimator object. A single observation manager can simulate observations and
         calculate observation partials for all link ends involved in the given observable type.


         :type: dict[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`, :class:`~tudatpy.numerical_simulation.estimation.ObservationManager` ]
      )doc" )
            .def_property_readonly(
                    "state_transition_interface",
                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getStateTransitionAndSensitivityMatrixInterface,
                    R"doc(

         **read-only**

         State transition and sensitivity matrix interface, setting the variational equations/dynamics in the
         Estimator object.


         :type: :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
      )doc" )
            .def( "perform_estimation",
                  &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                   TIME_TYPE >::estimateParameters,
                  py::arg( "estimation_input" ),
                  R"doc(

         Function to trigger the parameter estimation.


         Function to trigger the parameter estimation. Much of the process and requirements are similar to those described in the
         :func:`~tudatpy.numerical_simulation.Estimator.compute_covariance` function. This function uses an iterative least-squares
         estimate process to fit the data (inside ``estimation_input``) to the model defined by the inputs to the ``Estimator`` constructor.s


         Parameters
         ----------
         estimation_input : :class:`~tudatpy.numerical_simulation.estimation.EstimationInput`
             Object consolidating all relevant settings for the estimation
             This includes foremost the simulated observations, as well as a priori information about the estimatable parameters and convergence criteria for the least squares estimation.

         Returns
         -------
         :class:`~tudatpy.numerical_simulation.estimation.EstimationOutput`
             Object containing all outputs from the estimation process.





     )doc" )
            .def( "compute_covariance",
                  &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                   TIME_TYPE >::computeCovariance,
                  py::arg( "covariance_analysis_input" ),
                  R"doc(

         Function to perform a covariance analysis for the given observations and parameters


         Function to perform a covariance analysis for the given observations and parameters. The observations are provided through the
         ``covariance_analysis_input`` input, as are the weights :math:`\mathbf{W}` and inverse a priori covariance :math:`(\mathbf{P}_{0})^{-1}`.
         Calling this function uses the environment and propagator settings provided to the constructor of this `Estimator` class to simulate
         the dynamics of any relevant bodies for the observations (and associated variational equations). The observations are then
         computed using the observation models created by the settings provided to the constructor of this `Estimator` class, as is the
         associated design matrix :math:`\mathbf{H}`. This function then produces the covariance :math:`\mathbf{P}` (omitting the normalization used
         internally for numerical stability)

         .. math::
            \mathbf{P}=\left(\mathbf{H}^{T}\mathbf{W}\mathbf{H}+(\mathbf{P}_{0})^{-1}\right)^{-1}

         Note that, although the actual observations are formally not required for a covariance analysis, all additional data (e.g. observation time, type, link ends, etc.)
         are. And, as such, the ``covariance_analysis_input`` does require the full set of observations and associated information, for consistency purposes (e.g., same input as
         ``perform_estimation`` function) .


         Parameters
         ----------
         covariance_analysis_input : :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput`
             Object consolidating all relevant settings for the covariance analysis
             This includes foremost the simulated observations, as well as a priori information about the estimatable parameters

         Returns
         -------
         :class:`~tudatpy.numerical_simulation.estimation.vOutput`
             Object containing all outputs from the estimation process.





     )doc" )
            .def_property_readonly(
                    "variational_solver",
                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getVariationalEquationsSolver,
                    R"doc(

         **read-only**

         Variational equations solver, which is used to manage and execute the numerical integration of
         equations of motion and variational equations/dynamics in the Estimator object.

         :type: :class:`~tudatpy.numerical_simulation.SingleArcVariationalSimulator`
      )doc" );


    py::class_< tss::EstimationConvergenceChecker,
                std::shared_ptr< tss::EstimationConvergenceChecker > >(
            m,
            "EstimationConvergenceChecker",
            R"doc(

         Class defining the convergence criteria for an estimation.

         Class defining the convergence criteria for an estimation.
         The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation.estimation_convergence_checker` function.





      )doc" );

    m.def( "estimation_convergence_checker",
           &tss::estimationConvergenceChecker,
           py::arg( "maximum_iterations" ) = 5,
           py::arg( "minimum_residual_change" ) = 0.0,
           py::arg( "minimum_residual" ) = 0.0,
           py::arg( "number_of_iterations_without_improvement" ) = 2,
           R"doc(

 Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object.

 Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object, which is required for defining the convergence criteria of an estimation.


 Parameters
 ----------
 maximum_iterations : int, default = 5
     Maximum number of allowed iterations for estimation.
 minimum_residual_change : float, default = 0.0
     Minimum required change in residual between two iterations.
 minimum_residual : float, default = 0.0
     Minimum value of observation residual below which estimation is converged.
 number_of_iterations_without_improvement : int, default = 2
     Number of iterations without reduction of residual.
 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` class, defining the convergence criteria for an estimation.






     )doc" );

    py::class_< tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "CovarianceAnalysisInput",
            R"doc(

         Class for defining all specific inputs to a covariance analysis.





      )doc" )
            .def( py::init< const std::shared_ptr<
                                    tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >&,
                            const Eigen::MatrixXd,
                            const Eigen::MatrixXd >( ),
                  py::arg( "observations_and_times" ),
                  py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                  py::arg( "consider_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                  R"doc(

         Class constructor.

         Constructor through which the user can create instances of this class. Note that the weight are all initiated as 1.0, and the default settings of ``define_covariance_settings`` are used.


         Parameters
         ----------
         observations_and_times : ObservationCollection
             Total data structure of observations and associated times/link ends/type/etc.
         inverse_apriori_covariance : numpy.ndarray[numpy.float64[m, n]], default = [ ]
             A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
         Returns
         -------
         :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput`
             Instance of the :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput` class, defining the data and other settings to be used for the covariance analysis.





     )doc" )
            .def( "set_constant_weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE,
                                                 TIME_TYPE >::setConstantWeightsMatrix,
                  py::arg( "weight" ),
                  R"doc(

         Function to set a constant weight matrix for all observables.

         Function to set a constant weight matrix for all observables.
         The weights are applied to all observations managed by the given PodInput object.


         Parameters
         ----------
         constant_weight : float
             Constant weight factor that is to be applied to all observations.
         Returns
         -------
         None
             Function modifies the object in-place.





     )doc" )
            .def( "set_weights_from_observation_collection",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE,
                                                 TIME_TYPE >::setWeightsFromObservationCollection,
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_single_observable_weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE,
                                                 TIME_TYPE >::setConstantSingleObservableWeights,
                  py::arg( "observable_type" ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_single_observable_vector_weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::
                          setConstantSingleObservableVectorWeights,
                  py::arg( "observable_type" ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_single_observable_and_link_end_weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::
                          setConstantSingleObservableAndLinkEndsWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_single_observable_and_link_end_vector_"
                  "weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::
                          setConstantSingleObservableAndLinkEndsVectorWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_total_single_observable_and_link_end_vector_"
                  "weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::
                          setTabulatedSingleObservableAndLinkEndsWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight_vector" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_weight_per_observable",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE,
                                                 TIME_TYPE >::setConstantPerObservableWeightsMatrix,
                  py::arg( "weight_per_observable" ),
                  R"doc(

         Function to set a constant weight matrix for a given type of observable.

         Function to set a constant weight matrix for a given type of observable.
         The weights are applied to all observations of the observable type specified by the `weight_per_observable` parameter.


         Parameters
         ----------
         constant_weight : Dict[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`, float ]
             Constant weight factor that is to be applied to all observations.
         Returns
         -------
         None
             Function modifies the object in-place.





     )doc" )
            .def( "set_constant_vector_weight_per_observable",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::
                          setConstantPerObservableVectorWeightsMatrix,
                  py::arg( "weight_per_observable" ),
                  R"doc(No documentation found.)doc" )
            .def( "define_covariance_settings",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE,
                                                 TIME_TYPE >::defineCovarianceSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "limit_condition_number_for_warning" ) = 1.0E8,
                  R"doc(

         Function to define specific settings for covariance analysis process

         Function to define specific settings for covariance analysis process


         Parameters
         ----------
         reintegrate_equations : bool, default = True
             Boolean denoting whether the dynamics and variational equations are to be reintegrated
             or if existing values are to be used to perform first iteration.

         reintegrate_variational_equations : bool, default = True
             Boolean denoting whether the variational equations are to be reintegrated during estimation
             (if this is set to False, and ``reintegrate_equations`` to true, only the dynamics are re-integrated)

         save_design_matrix : bool, default = True
             Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\mathbf{H}` in the output. Setting this to false makes the
             :math:`\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.

         print_output_to_terminal : bool, default = True
             Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.

         Returns
         -------
         None
             Function modifies the object in-place.





     )doc" )
            .def_property( "weight_matrix_diagonal",
                           &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE,
                                                          TIME_TYPE >::getWeightsMatrixDiagonals,
                           &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE,
                                                          TIME_TYPE >::setWeightsMatrixDiagonals,
                           R"doc(

         **read-only**

         Complete diagonal of the weights matrix that is to be used

         :type: numpy.ndarray[numpy.float64[n, 1]]
      )doc" );

    py::class_< tss::EstimationInput< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::EstimationInput< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                "EstimationInput",
                                                                                R"doc(

         Class for defining all inputs to the estimation.





      )doc" )
            .def( py::init< const std::shared_ptr<
                                    tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >&,
                            const Eigen::MatrixXd,
                            std::shared_ptr< tss::EstimationConvergenceChecker >,
                            const Eigen::MatrixXd,
                            const Eigen::VectorXd,
                            const bool >( ),
                  py::arg( "observations_and_times" ),
                  py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                  py::arg( "convergence_checker" ) =
                          std::make_shared< tss::EstimationConvergenceChecker >( ),
                  py::arg( "consider_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                  py::arg( "consider_parameters_deviations" ) = Eigen::VectorXd::Zero( 0 ),
                  py::arg( "apply_final_parameter_correction" ) = true,
                  R"doc(

         Class constructor.

         Constructor through which the user can create instances of this class.


         Parameters
         ----------
         observations_and_times : ObservationCollection
             Total data structure of observations and associated times/link ends/type/etc.
         inverse_apriori_covariance : numpy.ndarray[numpy.float64[m, n]], default = [ ]
             A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
         convergence_checker : :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker`, default = :func:`~tudatpy.numerical_simulation.estimation.estimation_convergence_checker`
             Object defining when the estimation is converged.
         Returns
         -------
         :class:`~tudatpy.numerical_simulation.estimation.EstimationInput`
             Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationInput` class, defining the data and other settings to be used for the estimation.





     )doc" )
            .def( "define_estimation_settings",
                  &tss::EstimationInput< STATE_SCALAR_TYPE, TIME_TYPE >::defineEstimationSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "save_residuals_and_parameters_per_iteration" ) = true,
                  py::arg( "save_state_history_per_iteration" ) = false,
                  py::arg( "limit_condition_number_for_warning" ) = 1.0E8,
                  py::arg( "condition_number_warning_each_iteration" ) = true,
                  R"doc(

         Function to define specific settings for the estimation process

         Function to define specific settings for covariance analysis process


         Parameters
         ----------
         reintegrate_equations_on_first_iteration : bool, default = True
             Boolean denoting whether the dynamics and variational equations are to be reintegrated
             or if existing values are to be used to perform first iteration.

         reintegrate_variational_equations : bool, default = True
             Boolean denoting whether the variational equations are to be reintegrated during estimation
             (if this is set to False, and ``reintegrate_equations_on_first_iteration`` to true, only the dynamics are re-integrated)

         save_design_matrix : bool, default = True
             Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\mathbf{H}` in the output. Setting this to false makes the
             :math:`\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.

         print_output_to_terminal : bool, default = True
             Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.

         save_residuals_and_parameters_per_iteration : bool, default = True
             Boolean denoting whether the residuals and parameters from the each iteration are to be saved.

         save_state_history_per_iteration : bool, default = False
             Boolean denoting whether the state history and dependent variables are to be saved on each iteration.

         Returns
         -------
         None
             Function modifies the object in-place.





     )doc" );

     m.attr( "PodInput" ) = m.attr( "EstimationInput" );

    py::class_< tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "CovarianceAnalysisOutput",
            R"doc(

         Class collecting all outputs from the covariance analysis process.





      )doc" )
            .def_property_readonly( "inverse_covariance",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::
                                            getUnnormalizedInverseCovarianceMatrix,
                                    R"doc(

         **read-only**

         (Unnormalized) inverse estimation covariance matrix :math:`\mathbf{P}^{-1}`.

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )
            .def_property_readonly(
                    "covariance",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getUnnormalizedCovarianceMatrix,
                    R"doc(

         **read-only**

         (Unnormalized) estimation covariance matrix :math:`\mathbf{P}`.

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )
            .def_property_readonly( "inverse_normalized_covariance",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::
                                            getNormalizedInverseCovarianceMatrix,
                                    R"doc(

         **read-only**

         Normalized inverse estimation covariance matrix :math:`\mathbf{\tilde{P}}^{-1}`.

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )
            .def_property_readonly(
                    "normalized_covariance",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getNormalizedCovarianceMatrix,
                    R"doc(

         **read-only**

         Normalized estimation covariance matrix :math:`\mathbf{\tilde{P}}`.

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )
            .def_property_readonly(
                    "formal_errors",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getFormalErrorVector,
                    R"doc(

         **read-only**

         Formal error vector :math:`\boldsymbol{\sigma}` of the estimation result (e.g. square root of diagonal entries of covariance)s

         :type: numpy.ndarray[numpy.float64[m, 1]]s
      )doc" )
            .def_property_readonly(
                    "correlations",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getCorrelationMatrix,
                    R"doc(

         **read-only**

         Correlation matrix of the estimation result. Entry :math:`i,j` is equal to :math:`P_{i,j}/(\sigma_{i}\sigma_{j})`

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )
            .def_property_readonly(
                    "design_matrix",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getUnnormalizedDesignMatrix,
                    R"doc(

         **read-only**

         Matrix of unnormalized partial derivatives :math:`\mathbf{H}=\frac{\partial\mathbf{h}}{\partial\mathbf{p}}`.

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly(
                    "normalized_design_matrix",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getNormalizedDesignMatrix,
                    R"doc(

         **read-only**

         Matrix of normalized partial derivatives :math:`\tilde{\mathbf{H}}`.

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly( "weighted_design_matrix",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::
                                            getUnnormalizedWeightedDesignMatrix,
                                    R"doc(

         **read-only**

         Matrix of weighted partial derivatives, equal to :math:`\mathbf{W}^{1/2}{\mathbf{H}}`

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly(
                    "weighted_normalized_design_matrix",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getNormalizedWeightedDesignMatrix,
                    R"doc(

         **read-only**

         Matrix of weighted, normalized partial derivatives, equal to :math:`\mathbf{W}^{1/2}\tilde{\mathbf{H}}`

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly(
                    "consider_covariance_contribution",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getConsiderCovarianceContribution,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "normalized_covariance_with_consider_parameters",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::
                                            getNormalizedCovarianceWithConsiderParameters,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "unnormalized_covariance_with_consider_parameters",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::
                                            getUnnormalizedCovarianceWithConsiderParameters,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "normalized_design_matrix_consider_parameters",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::
                                            getNormalizedDesignMatrixConsiderParameters,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "consider_normalization_factors",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getConsiderNormalizationFactors,
                    R"doc(No documentation found.)doc" )
            .def_readonly( "normalization_terms",
                           &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::
                                   designMatrixTransformationDiagonal_,
                           R"doc(

         **read-only**

         Vector of normalization terms used for covariance and design matrix

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" );

    py::class_< tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                 "EstimationOutput",
                                                                                 R"doc(

         Class collecting all outputs from the iterative estimation process.





      )doc" )
            .def_property_readonly( "residual_history",
                                    &tss::EstimationOutput< STATE_SCALAR_TYPE,
                                                            TIME_TYPE >::getResidualHistoryMatrix,
                                    R"doc(

         **read-only**

         Residual vectors, concatenated per iteration into a matrix; the :math:`i^{th}` column has the residuals from the :math:`i^{th}` iteration.

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly( "parameter_history",
                                    &tss::EstimationOutput< STATE_SCALAR_TYPE,
                                                            TIME_TYPE >::getParameterHistoryMatrix,
                                    R"doc(

         **read-only**

         Parameter vectors, concatenated per iteration into a matrix. Column 0 contains pre-estimation values. The :math:`(i+1)^{th}` column has the residuals from the :math:`i^{th}` iteration.

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly(
                    "simulation_results_per_iteration",
                    &tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getSimulationResults,
                    R"doc(

         **read-only**

         List of complete numerical propagation results, with the :math:`i^{th}` entry of thee list thee results of the :math:`i^{th}` propagation

         :type: list[SimulationResults]
      )doc" )
            .def_readonly( "final_residuals",
                           &tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::residuals_,
                           R"doc(

         **read-only**

         Vector of post-fit observation residuals, for the iteration with the lowest rms residuals.

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_readonly(
                    "final_parameters",
                    &tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::parameterEstimate_,
                    R"doc(No documentation found.)doc" )
            .def_readonly( "best_iteration",
                           &tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::bestIteration_,
                           R"doc(No documentation found.)doc" );

    m.attr( "PodOutput" ) = m.attr( "EstimationOutput" );


    m.def( "create_best_fit_to_ephemeris",
           &tss::createBestFitToCurrentEphemeris< TIME_TYPE, STATE_SCALAR_TYPE >,
           py::arg( "bodies" ),
           py::arg( "acceleration_models" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "integrator_settings" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "data_point_interval" ),
           py::arg( "additional_parameter_names" ) = std::vector< std::shared_ptr< tep::EstimatableParameterSettings > >( ),
           py::arg( "number_of_iterations" ) = 3,
           py::arg( "reintegrate_variational_equations" ) = true,
           py::arg( "results_print_frequency" ) = 0.0,
           R"doc(No documentation found.)doc" );


    // PROPAGATE COVARIANCE

    m.def( "propagate_covariance_rsw_split_output",
           &tp::propagateCovarianceVectorsRsw,
           py::arg( "covariance_output" ),
           py::arg( "estimator" ),
           py::arg( "output_times" ),
           R"doc(No documentation found.)doc" );

    m.def( "propagate_formal_errors_rsw_split_output",
           &tp::propagateFormalErrorVectorsRsw,
           py::arg( "covariance_output" ),
           py::arg( "estimator" ),
           py::arg( "output_times" ),
           R"doc(No documentation found.)doc" );

    m.def( "propagate_covariance_split_output",
           py::overload_cast< const Eigen::MatrixXd,
                              const std::shared_ptr<
                                      tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
                              const std::vector< double > >( &tp::propagateCovarianceVectors ),
           py::arg( "initial_covariance" ),
           py::arg( "state_transition_interface" ),
           py::arg( "output_times" ),
           R"doc(No documentation found.)doc" );

    m.def( "propagate_covariance",
           py::overload_cast< const Eigen::MatrixXd,
                              const std::shared_ptr<
                                      tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
                              const std::vector< double > >( &tp::propagateCovariance ),
           py::arg( "initial_covariance" ),
           py::arg( "state_transition_interface" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system covariance through time.

 Function to propagate the covariance of a given system through time.
 The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.

 output_times : List[ float ]
     Times at which the propagated covariance matrix shall be reported.
     Note that this argument has no impact on the integration time-steps of the covariance propagation,
     which always adheres to the integrator settings that the `state_transition_interface` links to.
     Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 Dict[ float, numpy.ndarray[numpy.float64[m, n]] ]
     Dictionary reporting the propagated covariances at each output time.






     )doc" );

    m.def( "propagate_formal_errors_split_output",
           py::overload_cast< const Eigen::MatrixXd,
                              const std::shared_ptr<
                                      tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
                              const std::vector< double > >( &tp::propagateFormalErrorVectors ),
           py::arg( "initial_covariance" ),
           py::arg( "state_transition_interface" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system formal errors through time.

 Function to propagate the formal errors of a given system through time.
 Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
 The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.

 output_times : List[ float ]
     Times at which the propagated covariance matrix shall be reported.
     Note that this argument has no impact on the integration time-steps of the covariance propagation,
     which always adheres to the integrator settings that the `state_transition_interface` links to.
     Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 Dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
     Dictionary reporting the propagated formal errors at each output time.






     )doc" );

    m.def( "propagate_formal_errors",
           py::overload_cast< const Eigen::MatrixXd,
                              const std::shared_ptr<
                                      tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
                              const std::vector< double > >( &tp::propagateFormalErrors ),
           py::arg( "initial_covariance" ),
           py::arg( "state_transition_interface" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system formal errors through time.

 Function to propagate the formal errors of a given system through time.
 Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
 The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.

 output_times : List[ float ]
     Times at which the propagated covariance matrix shall be reported.
     Note that this argument has no impact on the integration time-steps of the covariance propagation,
     which always adheres to the integrator settings that the `state_transition_interface` links to.
     Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 Dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
     Dictionary reporting the propagated formal errors at each output time.






     )doc" );

}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
