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
#include "tudat/astro/propagators/propagateCovariance.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

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
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::cerr << "The propagation of covariances converted to RSW frame is deprecated as of v1.0 due to insufficient flexibility and lack "
                 "of use cases"
              << std::endl;

    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    tp::propagateCovariance( propagatedCovariance,
                             initialCovariance,
                             orbitDeterminationManager->getStateTransitionAndSensitivityMatrixInterface( ),
                             evaluationTimes );

    tss::SystemOfBodies bodies = orbitDeterminationManager->getBodies( );

    std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > > parameterSet =
            orbitDeterminationManager->getParametersToEstimate( );

    std::map< int, std::shared_ptr< tep::EstimatableParameter< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > > > > initialStates =
            parameterSet->getInitialStateParameters( );
    std::map< std::pair< std::string, std::string >, std::vector< int > > transformationList;
    for( auto it : initialStates )
    {
        if( std::dynamic_pointer_cast< tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > >( it.second ) )
        {
            std::shared_ptr< tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > > currentInitialState =
                    std::dynamic_pointer_cast< tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > >( it.second );
            transformationList[ std::make_pair( currentInitialState->getParameterName( ).second.first,
                                                currentInitialState->getCentralBody( ) ) ]
                    .push_back( it.first );
        }
        else if( std::dynamic_pointer_cast< tep::ArcWiseInitialTranslationalStateParameter< STATE_SCALAR_TYPE > >( it.second ) )
        {
            throw std::runtime_error(
                    "Error, multi-arc not yet supported in automatic "
                    "covariance conversion" );
        }
    }

    Eigen::Matrix3d currentInertialToRswPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix6d currentInertialToRswState = Eigen::Matrix6d::Zero( );
    Eigen::MatrixXd currentFullInertialToRswState = Eigen::MatrixXd::Zero( 6, 6 );

    std::map< double, Eigen::MatrixXd > propagatedRswCovariance;
    for( auto it : propagatedCovariance )
    {
        double currentTime = static_cast< double >( it.first );
        Eigen::MatrixXd currentCovariance = it.second;
        currentFullInertialToRswState.setZero( );

        for( auto it_body : transformationList )
        {
            Eigen::Vector6d relativeState = bodies.getBody( it_body.first.first )->getStateInBaseFrameFromEphemeris( currentTime ) -
                    bodies.getBody( it_body.first.second )->getStateInBaseFrameFromEphemeris( currentTime );
            currentInertialToRswPosition = trf::getInertialToRswSatelliteCenteredFrameRotationMatrix( relativeState );
            currentInertialToRswState.block( 0, 0, 3, 3 ) = currentInertialToRswPosition;
            currentInertialToRswState.block( 3, 3, 3, 3 ) = currentInertialToRswPosition;
            for( unsigned int j = 0; j < it_body.second.size( ); j++ )
            {
                int currentStartIndex = it_body.second.at( j );
                currentFullInertialToRswState.block( currentStartIndex, currentStartIndex, 6, 6 ) = currentInertialToRswState;
            }
        }
        propagatedRswCovariance[ currentTime ] =
                currentFullInertialToRswState * currentCovariance * currentFullInertialToRswState.transpose( );
    }
    return propagatedRswCovariance;
}

std::pair< std::vector< double >, std::vector< Eigen::MatrixXd > > propagateCovarianceVectorsRsw(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedRswCovariance =
            propagateCovarianceRsw( initialCovariance, orbitDeterminationManager, evaluationTimes );

    return std::make_pair( utilities::createVectorFromMapKeys( propagatedRswCovariance ),
                           utilities::createVectorFromMapValues( propagatedRswCovariance ) );
}

std::map< double, Eigen::VectorXd > propagateFormalErrorsRsw(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    std::map< double, Eigen::VectorXd > propagatedFormalErrors;

    propagatedCovariance = propagateCovarianceRsw( initialCovariance, orbitDeterminationManager, evaluationTimes );
    tp::convertCovarianceHistoryToFormalErrorHistory( propagatedFormalErrors, propagatedCovariance );

    return propagatedFormalErrors;
}

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > propagateFormalErrorVectorsRsw(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::cerr << "The propagate_covariance_rsw_split_output function is deprecated as of v1.0, use propagate_covariance_rsw instead"
              << std::endl;

    std::map< double, Eigen::VectorXd > propagatedFormalErrors =
            propagateFormalErrorsRsw( initialCovariance, orbitDeterminationManager, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys( propagatedFormalErrors ),
                           utilities::createVectorFromMapValues( propagatedFormalErrors ) );
}

std::pair< std::vector< double >, std::vector< Eigen::MatrixXd > > propagateCovarianceVectors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::cerr << "The propagate_covariance_split_output function is deprecated as of v1.0, use propagate_covariance instead" << std::endl;
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    tp::propagateCovariance( propagatedCovariance, initialCovariance, stateTransitionInterface, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys( propagatedCovariance ),
                           utilities::createVectorFromMapValues( propagatedCovariance ) );
}

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > propagateFormalErrorVectors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::cerr << "The propagate_formal_errors_split_output function is deprecated as of v1.0, use propagate_formal_errors instead"
              << std::endl;

    std::map< double, Eigen::VectorXd > propagatedFormalErrors;
    tp::propagateFormalErrors( propagatedFormalErrors, initialCovariance, stateTransitionInterface, evaluationTimes );
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

    py::class_< tss::EstimationConvergenceChecker, std::shared_ptr< tss::EstimationConvergenceChecker > >( m,
                                                                                                           "EstimationConvergenceChecker",
                                                                                                           R"doc(

         Class defining the convergence criteria for an estimation.

         Class defining the convergence criteria for an estimation.
         The user typically creates instances of this class via the :func:`~tudatpy.estimation.estimation_analysis.estimation_convergence_checker` function.





      )doc" );

    m.def( "estimation_convergence_checker",
           &tss::estimationConvergenceChecker,
           py::arg( "maximum_iterations" ) = 5,
           py::arg( "minimum_residual_change" ) = 0.0,
           py::arg( "minimum_residual" ) = 0.0,
           py::arg( "number_of_iterations_without_improvement" ) = 2,
           R"doc(

 Function for creating an :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker` object.

 Function for creating an :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker` object, which is required for defining the convergence criteria of an estimation.


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
 :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker`
     Instance of the :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker` class, defining the convergence criteria for an estimation.






     )doc" );

    py::class_< tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                   "CovarianceAnalysisInput",
                                                                                                   R"doc(

         Class for defining all inputs to a covariance analysis.


      )doc" )
            .def( py::init< const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >&,
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
         :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisInput`
             Instance of the :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisInput` class, defining the data and other settings to be used for the covariance analysis.





     )doc" )
            .def( "set_constant_weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightsMatrix,
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
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightsFromObservationCollection,
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_single_observable_weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservableWeights,
                  py::arg( "observable_type" ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_single_observable_vector_weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservableVectorWeights,
                  py::arg( "observable_type" ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_single_observable_and_link_end_weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservableAndLinkEndsWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_single_observable_and_link_end_vector_"
                  "weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservableAndLinkEndsVectorWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_total_single_observable_and_link_end_vector_"
                  "weight",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedSingleObservableAndLinkEndsWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight_vector" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_weight_per_observable",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantPerObservableWeightsMatrix,
                  py::arg( "weight_per_observable" ),
                  R"doc(

         Function to set a constant weight matrix for a given type of observable.

         Function to set a constant weight matrix for a given type of observable.
         The weights are applied to all observations of the observable type specified by the `weight_per_observable` parameter.


         Parameters
         ----------
         constant_weight : Dict[ :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`, float ]
             Constant weight factor that is to be applied to all observations.
         Returns
         -------
         None
             Function modifies the object in-place.





     )doc" )
            .def( "set_constant_vector_weight_per_observable",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantPerObservableVectorWeightsMatrix,
                  py::arg( "weight_per_observable" ),
                  R"doc(No documentation found.)doc" )
            .def( "define_covariance_settings",
                  &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::defineCovarianceSettings,
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
                           &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsMatrixDiagonals,
                           &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightsMatrixDiagonals,
                           R"doc(

         **read-only**

         Complete diagonal of the weights matrix that is to be used

         :type: numpy.ndarray[numpy.float64[n, 1]]
      )doc" )
            .def_property_readonly(
                    "inverse_apriori_covariance",
                    py::overload_cast<>( &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::getInverseOfAprioriCovariance ),
                    R"doc(

         **read-only**

         Inverse a-priori covariance matrix of the estimated parameters.

         :type: numpy.ndarray[numpy.float64[n, n]]
      )doc" )
            .def_property_readonly( "consider_covariance",
                                    &tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE >::getConsiderCovariance,
                                    R"doc(

         **read-only**

         A-priori covariance matrix of the considered parameters.

         :type: numpy.ndarray[numpy.float64[n, n]]
      )doc" );

    py::class_< tss::EstimationInput< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::EstimationInput< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tss::CovarianceAnalysisInput< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                "EstimationInput",
                                                                                R"doc(

         Class for defining all inputs to the estimation.





      )doc" )
            .def( py::init< const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >&,
                            const Eigen::MatrixXd,
                            std::shared_ptr< tss::EstimationConvergenceChecker >,
                            const Eigen::MatrixXd,
                            const Eigen::VectorXd,
                            const bool >( ),
                  py::arg( "observations_and_times" ),
                  py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                  py::arg( "convergence_checker" ) = std::make_shared< tss::EstimationConvergenceChecker >( ),
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
         convergence_checker : :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker`, default = :func:`~tudatpy.estimation.estimation_analysis.estimation_convergence_checker`
             Object defining when the estimation is converged.
         Returns
         -------
         :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
             Instance of the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput` class, defining the data and other settings to be used for the estimation.





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
                std::shared_ptr< tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                    "CovarianceAnalysisOutput",
                                                                                                    R"doc(

         Class collecting all outputs from the covariance analysis process (which takes a :class:`~CovarianceAnalysisInput` as input)

         This object is typically created by the :attr:`~tudatpy.estimation.estimation_analysis.Estimator.compute_covariance` or
         :attr:`~tudatpy.estimation.estimation_analysis.Estimator.perform_estimation` function of the :class:`~tudatpy.estimation.estimation_analysis.Estimator` class.
         The primary inputs used to create this object are (see `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/estimation-settings.html#covariance-analysis-settings>`_ for underlying models)

         * The partials matrix :math:`\mathbf{H}=\frac{\partial\mathbf{h}}{\partial\mathbf{p}}` of the observations w.r.t. the estimated parameters
         * The inverse covariance matrix :math:`\mathbf{P}^{-1}` of the estimated parameters (without influence of consider parameters). The inverse covariance is provided as input for situations where the inverse is unstable
         * The consider partials matrix :math:`\mathbf{H}_{c}=\frac{\partial\mathbf{h}}{\partial\mathbf{p}_{c}}` of the observations w.r.t. the consider parameters (if any)
         * The contribution :math:`\Delta \mathbf{P}_{c}` of the consider parameters to the estimated parameter covariance

         In the computation of the covariance  (see TODO), the columns of the :math:`H` matrices are normalized to reduce numerical instability
         that can result from the partials w.r.t. different parameters being of a very different order of magnitude. The normalization is achieved
         by computing a vector :math:`\mathbf{N}` (of the same size as the parameter vector :math:`\mathbf{p}`, such that for each column of the matrix :math:`\mathbf{H}`, we have:

         .. math::

             \max_{i}\left| \frac{H_{ij}}{N_{j}}\right|=1

         That is, the entries of :math:`\mathbf{N}` are chosen such that they normalize the corresponding column of :math:`\mathbf{H}` to be
         in the range :math:`[-1,1]`. We denote the normalized quantities with a tilde, so that:

         .. math::

             \tilde{H}_{ij}=\frac{H_{ij}}{N{j}}\\
             \tilde{P}_{ij}=P_{ij}N_{i}N_{j}

         When wanting to recreate the internal workings of the analysis, use the normalized quantities,
         when interested in the actual covariances, sensitivities, etc of the observations/parameters, use the unnormalized quantities.

      )doc" )
            .def_property_readonly( "inverse_covariance",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getUnnormalizedInverseCovarianceMatrix,
                                    R"doc(

         **read-only**

         (Unnormalized) inverse estimation covariance matrix :math:`\mathbf{P}^{-1}`. Note: if the any consider parameters are included in the analysis,
         this attribute does not include their contribution. The contribution of the consider parameters to the covariance can be retrieved
         from the :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.consider_covariance_contribution` attribute

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )
            .def_property_readonly( "covariance",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getUnnormalizedCovarianceMatrix,
                                    R"doc(

         **read-only**

         (Unnormalized) estimation covariance matrix :math:`\mathbf{P}`. Note: if the any consider parameters are included in the analysis,
         this attribute does not include their contribution. The contribution of the consider parameters to the covariance can be retrieved
         from the :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.consider_covariance_contribution` attribute

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )
            .def_property_readonly( "inverse_normalized_covariance",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getNormalizedInverseCovarianceMatrix,
                                    R"doc(

         **read-only**

         Normalized version of :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.inverse_covariance`

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )
            .def_property_readonly( "normalized_covariance",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getNormalizedCovarianceMatrix,
                                    R"doc(

         **read-only**

         Normalized version of :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.covariance`

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )
            .def_property_readonly( "formal_errors",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getFormalErrorVector,
                                    R"doc(

         **read-only**

         Formal error vector :math:`\boldsymbol{\sigma}` of the estimation result (e.g. square root of diagonal entries of covariance :math:`\mathbf{P}`)

         :type: numpy.ndarray[numpy.float64[m, 1]]s
      )doc" )
            .def_property_readonly( "correlations",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getCorrelationMatrix,
                                    R"doc(

         **read-only**

         Correlation matrix of the estimation result. Entry :math:`i,j` is equal to :math:`P_{i,j}/(\sigma_{i}\sigma_{j})`

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )

            .def_property_readonly( "consider_covariance",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getConsiderCovariance,
                                    R"doc(

         **read-only**

         Covariance :math:`\mathbf{C}` of consider parameters used in calculations

         :type: numpy.ndarray[numpy.float64[m, m]]
      )doc" )

            .def_property_readonly( "design_matrix",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getUnnormalizedDesignMatrix,
                                    R"doc(

         **read-only**

         Matrix of partial derivatives :math:`\mathbf{H}=\frac{\partial\mathbf{h}}{\partial\mathbf{p}}`.

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly( "normalized_design_matrix",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getNormalizedDesignMatrix,
                                    R"doc(

         **read-only**

         Normalized version of :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.design_matrix`

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly( "weighted_design_matrix",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getUnnormalizedWeightedDesignMatrix,
                                    R"doc(

         **read-only**

         Matrix of weighted partial derivatives, equal to :math:`\mathbf{W}^{1/2}{\mathbf{H}}`

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly( "weighted_normalized_design_matrix",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getNormalizedWeightedDesignMatrix,
                                    R"doc(

         **read-only**

         Normalized version of :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.weighted_design_matrix`

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly( "consider_covariance_contribution",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getConsiderCovarianceContribution,
                                    R"doc(

         **read-only**

         Contribution of the consider parameters to the estimated parameter covariance matrix, equal to :math:`(\mathbf{P} \mathbf{H}^{T} \mathbf{W}) (\mathbf{H}_c \mathbf{C} \mathbf{H}^{T}_c) (\mathbf{P} \mathbf{H}^{T} \mathbf{W})^{T}`

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly(
                    "covariance_with_consider_parameters",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getUnnormalizedCovarianceWithConsiderParameters,
                    R"doc(

         **read-only**

         Covariance matrix of the estimated parameters that includes the contribution of the consider parameters, 
         equal to the sum of the :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.covariance` matrix and the :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.consider_covariance_contribution` matrix.


         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly(
                    "normalized_covariance_with_consider_parameters",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getNormalizedCovarianceWithConsiderParameters,
                    R"doc(

         **read-only**

         Normalized version of :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.covariance_with_consider_parameters`

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly(
                    "design_matrix_consider_parameters",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getUnnormalizedDesignMatrixConsiderParameters,
                    R"doc(

         **read-only**

         Matrix of partial derivatives of observations w.r.t. consider parameters :math:`\mathbf{H}_c=\frac{\partial\mathbf{h}}{\partial\mathbf{c}}`.

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly(
                    "normalized_design_matrix_consider_parameters",
                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getNormalizedDesignMatrixConsiderParameters,
                    R"doc(

         **read-only**

         Normalized version of :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.design_matrix_consider_parameters`

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )

            .def_readonly( "normalization_terms",
                           &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::designMatrixTransformationDiagonal_,
                           R"doc(

         **read-only**

         Vector of estimated parameter normalization terms :math:`\math{N}`

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly( "consider_normalization_terms",
                                    &tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getConsiderNormalizationFactors,
                                    R"doc(

         **read-only**

         Vector of consider parameter normalization terms :math:`\mathbf{N}_{c}`

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
                                    &tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualHistoryMatrix,
                                    R"doc(

         **read-only**

         Residual vectors, concatenated per iteration into a matrix; the :math:`i^{th}` column has the residuals from the :math:`i^{th}` iteration.

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly( "parameter_history",
                                    &tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getParameterHistoryMatrix,
                                    R"doc(

         **read-only**

         Parameter vectors, concatenated per iteration into a matrix. Column 0 contains pre-estimation values. The :math:`(i+1)^{th}` column has the residuals from the :math:`i^{th}` iteration.

         :type: numpy.ndarray[numpy.float64[m, n]]
      )doc" )
            .def_property_readonly( "simulation_results_per_iteration",
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
            .def_readonly( "final_parameters",
                           &tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::parameterEstimate_,
                           R"doc(No documentation found.)doc" )
            .def_readonly( "best_iteration",
                           &tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::bestIteration_,
                           R"doc(No documentation found.)doc" );

    m.attr( "PodOutput" ) = m.attr( "EstimationOutput" );

    // PROPAGATE COVARIANCE

    m.def( "propagate_covariance",
           py::overload_cast< const Eigen::MatrixXd,
                              const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
                              const std::vector< double > >( &tp::propagateCovariance ),
           py::arg( "initial_covariance" ),
           py::arg( "state_transition_interface" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system covariance through time.

 Function to propagate the covariance of a given system through time. The covariance of the states and parameters at the reference epoch
 :math:`\mathbf{P}(t_{0})` is to be provided as input by the user. The definition of the parameters provided in the entries of this covariance must
 be consistent with those in the associated ``state_transition_interface`` input. Splitting the covariance into dynamical parameters
 (initial states :math:`\mathbf{x_{0}}=\mathbf{x}_{0}`) and static parameters :math:\mathbf{p}`, we can write:

 .. math::
    \mathbf{P}(t_{0}) &= \begin{pmatrix} \mathbf{P}_{x0,x0} & \mathbf{P}_{x0,p} \\ \mathbf{P}_{p,x0} & \mathbf{P}_{p,p} \end{pmatrix}\\
    \boldsymbol{\Psi(t,t_{0})} &= \begin{pmatrix} \boldsymbol{\Phi}(t,t_{0}) & \mathbf{S}(t) \\ \mathbf{0} & \mathbf{1} \end{pmatrix}

 with :math:`\boldsymbol{\Phi}(t,t_{0})` and :math:`\mathbf{S}(t)` the state transition and sensitivity matrices. The covariance is then propagated
 to any time :math:`t` using:

 .. math::
   \mathbf{P}(t) &= \boldsymbol{\Psi(t,t_{0})} \mathbf{P}(t_{0}) \boldsymbol{\Psi(t,t_{0})}^{T}
                 &= \begin{pmatrix} \mathbf{P}_{x,x} & \mathbf{P}_{x,p} \\ \mathbf{P}_{p,x} & \mathbf{P}_{p,p} \end{pmatrix}

 To automatically link the covariance computed by an :class:`~tudatpy.estimation.estimation_analysis.Estimator` object, use the
 :func:`~tudatpy.estimation.estimation_analysis.propagate_covariance_from_analysis_objects` function

 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[n, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by ``state_transition_interface``)

 state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time (typically retrieved from :attr:`~tudatpy.estimation.estimation_analysis.Estimator.state_transition_interface`).

 output_times : List[ astro.time_representation.Time ]
     Times at which the propagated covariance matrix shall be reported.
     Note that this argument has no impact on the integration time-steps of the variational equations,
     which happens before the call to this function, with results stored in the ``state_transition_interface`` input.
     Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 Dict[ astro.time_representation.Time, numpy.ndarray[numpy.float64[n, n]] ]
     Dictionary reporting the propagated covariances at each output time.


     )doc" );

    m.def( "propagate_covariance_from_analysis_objects",
           &tss::propagateCovarianceFromObjects< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "analysis_output" ),
           py::arg( "state_transition_interface" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system covariance through time.

 Function to propagate the covariance of a given system through time. Functionality is similar as :func:`~tudatpy.estimation.estimation_analysis.propagate_covariance`
 with the difference that the initial covariance is automatically extracted from the ``analysis_output`` input.

 This function also automatically deals with the presence of any consider parameters :math:`\mathbf{p}_{c}`. If these are present (with covariance :math:`\mathbf{C}`,  taken from :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.consider_covariance`)
 the system propagates the matrix :math:`\bar{\mathbf{P}}` constructed from the consider-parameter-free matrix :math:`\mathbf{P}`  (taken from :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.covariance`) as:

 .. math::
    \bar{\mathbf{P}}(t_{0}) &= \begin{pmatrix} \mathbf{P} + \Delta\mathbf{P}_{c} & \mathbf{0}  \\ \mathbf{0} & \mathbf{C} \end{pmatrix}

 where :math:`\Delta\mathbf{P}_{c}` (taken from :attr:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput.consider_covariance_contribution`) is the contribution of the consider parameters to the estimated parameter covariance at :math:`t_{0}`

 Parameters
 ----------
 analysis_output : CovarianceAnalysisOutput
     Object storing all results of a covariance analysis (using :meth:`~tudatpy.estimation.estimation_analysis.Estimator.compute_covariance`) or estimation (using :meth:`~tudatpy.estimation.estimation_analysis.Estimator.perform_estimation`)

 state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time (typically retrieved from :attr:`~tudatpy.estimation.estimation_analysis.Estimator.state_transition_interface`).

 output_times : List[ astro.time_representation.Time ]
     Times at which the propagated covariance matrix shall be reported.
     Note that this argument has no impact on the integration time-steps of the covariance propagation,
     which always adheres to the integrator settings that the `state_transition_interface` links to.
     Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 Dict[ astro.time_representation.Time, numpy.ndarray[numpy.float64[m, n]] ]
     Dictionary reporting the propagated covariances at each output time.


     )doc" );

    m.def( "propagate_formal_errors",
           py::overload_cast< const Eigen::MatrixXd,
                              const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
                              const std::vector< double > >( &tp::propagateFormalErrors ),
           py::arg( "initial_covariance" ),
           py::arg( "state_transition_interface" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system formal errors through time.

 Identical to :func:`~tudatpy.estimation.estimation_analysis.propagate_covariance`, but returning only the formal errors :math:\sigma_{i}` (from diagonal of covariance :math:`P_{ii}=\sigma_{i}^{2}`)

 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
    Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time (typically retrieved from :attr:`~tudatpy.estimation.estimation_analysis.Estimator.state_transition_interface`).

 output_times : List[ astro.time_representation.Time ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the variational equations,
    which happens before the call to this function, with results stored in the ``state_transition_interface`` input.
    Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 Dict[ astro.time_representation.Time, numpy.ndarray[numpy.float64[m, 1]] ]
     Dictionary reporting the propagated formal errors at each output time.






     )doc" );

    /************************** DEPRECATED ***************************/

    m.def( "propagate_covariance_split_output",
           py::overload_cast< const Eigen::MatrixXd,
                              const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
                              const std::vector< double > >( &tp::propagateCovarianceVectors ),
           py::arg( "initial_covariance" ),
           py::arg( "state_transition_interface" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system covariance through time.

 Identical to :func:`~tudatpy.estimation.estimation_analysis.propagate_covariance`, but with two lists rather than a dictionary as output

 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time (typically retrieved from :attr:`~tudatpy.estimation.estimation_analysis.Estimator.state_transition_interface`).

 output_times : List[ astro.time_representation.Time ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the variational equations,
    which happens before the call to this function, with results stored in the ``state_transition_interface`` input.
    Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 tuple[ list[astro.time_representation.Time], list[numpy.ndarray[numpy.float64[m, n]]] ]
     Tuple containing a list of output times, and a list of propagated covariances at each output time.



     )doc" );

    m.def( "propagate_covariance_rsw_split_output",
           &tp::propagateCovarianceVectorsRsw,
           py::arg( "initial_covariance" ),
           py::arg( "estimator" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system covariance through time and convert it to RSW frame.

 The covariance of a given system is propagated through time and afterwards converted to RSW frame.
 The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time (typically retrieved from :attr:`~tudatpy.estimation.estimation_analysis.Estimator.state_transition_interface`).

 output_times : List[ astro.time_representation.Time ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the variational equations,
    which happens before the call to this function, with results stored in the ``state_transition_interface`` input.
    Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 tuple[ list[float], list[numpy.ndarray[numpy.float64[m, n]]] ]
     Tuple containing a list of output times, and a list of propagated covariances in RSW frame at each output time.






     )doc" );

    m.def( "propagate_formal_errors_rsw_split_output",
           &tp::propagateFormalErrorVectorsRsw,
           py::arg( "initial_covariance" ),
           py::arg( "estimator" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system formal errors through time and convert to RSW frame.

 Function to propagate the formal errors of a given system through time.
 Note that in practice the entire covariance matrix is propagated and converted to RSW frame, but only the formal errors (variances) are reported at the output times.
 The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time (typically retrieved from :attr:`~tudatpy.estimation.estimation_analysis.Estimator.state_transition_interface`).


 output_times : List[ astro.time_representation.Time ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the variational equations,
    which happens before the call to this function, with results stored in the ``state_transition_interface`` input.
    Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 tuple[ list[astro.time_representation.Time], list[numpy.ndarray[numpy.float64[m, n]]] ]
     Tuple containing a list of output times, and a list of propagated formal errors in RSW frame at each output time.


     )doc" );

    m.def( "propagate_formal_errors_split_output",
           py::overload_cast< const Eigen::MatrixXd,
                              const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
                              const std::vector< double > >( &tp::propagateFormalErrorVectors ),
           py::arg( "initial_covariance" ),
           py::arg( "state_transition_interface" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system formal errors through time.

 Identical to :func:`~tudatpy.estimation.estimation_analysis.propagate_formal_errors, but with two lists rather than a dictionary as output

 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time (typically retrieved from :attr:`~tudatpy.estimation.estimation_analysis.Estimator.state_transition_interface`).


 output_times : List[ astro.time_representation.Time ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the variational equations,
    which happens before the call to this function, with results stored in the ``state_transition_interface`` input.
    Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 tuple[ list[astro.time_representation.Time], list[numpy.ndarray[numpy.float64[m, n]]] ]
     Tuple containing a list of output times, and a list of propagated formal errors at each output time.






     )doc" );

    m.def( "propagate_covariance_rsw_split_output",
           &tp::propagateCovarianceVectorsRsw,
           py::arg( "initial_covariance" ),
           py::arg( "estimator" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system covariance through time and convert it to RSW frame.

 The covariance of a given system is propagated through time and afterwards converted to RSW frame.
 The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time (typically retrieved from :attr:`~tudatpy.estimation.estimation_analysis.Estimator.state_transition_interface`).

 output_times : List[ astro.time_representation.Time ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the variational equations,
    which happens before the call to this function, with results stored in the ``state_transition_interface`` input.
    Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 tuple[ list[float], list[numpy.ndarray[numpy.float64[m, n]]] ]
     Tuple containing a list of output times, and a list of propagated covariances in RSW frame at each output time.






     )doc" );

    m.def( "propagate_formal_errors_rsw_split_output",
           &tp::propagateFormalErrorVectorsRsw,
           py::arg( "initial_covariance" ),
           py::arg( "estimator" ),
           py::arg( "output_times" ),
           R"doc(

 Function to propagate system formal errors through time and convert to RSW frame.

 Function to propagate the formal errors of a given system through time.
 Note that in practice the entire covariance matrix is propagated and converted to RSW frame, but only the formal errors (variances) are reported at the output times.
 The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


 Parameters
 ----------
 initial_covariance : numpy.ndarray[numpy.float64[m, n]]
     System covariance matrix (symmetric and positive semi-definite) at initial time.
     Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

 state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
     Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time (typically retrieved from :attr:`~tudatpy.estimation.estimation_analysis.Estimator.state_transition_interface`).


 output_times : List[ astro.time_representation.Time ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the variational equations,
    which happens before the call to this function, with results stored in the ``state_transition_interface`` input.
    Output times which do not coincide with integration time steps are calculated via interpolation.

 Returns
 -------
 tuple[ list[astro.time_representation.Time], list[numpy.ndarray[numpy.float64[m, n]]] ]
     Tuple containing a list of output times, and a list of propagated formal errors in RSW frame at each output time.


     )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
