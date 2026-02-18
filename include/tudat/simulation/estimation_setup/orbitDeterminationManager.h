/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ORBITDETERMINATIONMANAGER_H
#define TUDAT_ORBITDETERMINATIONMANAGER_H

#include <map>
#include <memory>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolverBase.h"
#include "tudat/simulation/estimation_setup/estimationInterfacesForwardDeclarations.h"
#include "tudat/simulation/estimation_setup/observationInterfacesForwardDeclarations.h"

namespace tudat
{

namespace estimatable_parameters
{
template< typename ParameterScalarType >
class EstimatableParameterSet;
}

namespace numerical_integrators
{
template< typename TimeType >
class IntegratorSettings;
}

namespace propagators
{
template< typename StateScalarType, typename TimeType >
class SimulationResults;

template< typename StateScalarType >
class PropagatorSettings;

class CombinedStateTransitionAndSensitivityMatrixInterface;

template< typename TimeType >
class DependentVariablesInterface;
}  // namespace propagators

namespace simulation_setup
{

//! Top-level class for performing orbit determination.
/*!
 *  Top-level class for performing orbit determination. All required propagation/estimation settings are provided to
 *  this class, which then creates all objects needed for the propagation and estimation process. The parameter
 *  estimation itself is performed by providing measurement data and related metadata (as EstimationInput) to the estimateParameters
 *  function.
 */
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy = 0 >
class OrbitDeterminationManager
{
public:
    //! Typedef for vector of observations.
    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;

    //! Typedef for vector of parameters.
    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ParameterVectorType;

    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > preprocessDeprecatedIntegratorSettings(
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
            const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettings,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const int integratorIndexOffset = 0 );

    OrbitDeterminationManager(
            const SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
            const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
            const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettings,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const bool propagateOnCreation = true );

    //! Constructor
    /*!
     *  Constructor
     *  \param bodies Map of body objects with names of bodies, storing all environment models used in simulation.
     *  \param parametersToEstimate Container object for all parameters that are to be estimated
     *  \param observationSettingsMap Sets of observation model settings per link ends (i.e. transmitter, receiver, etc.)
     *  per observable type for which measurement data is to be provided in orbit determination process
     *  (through estimateParameters function)
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param propagateOnCreation Boolean denoting whether initial propagatoon is to be performed upon object creation (default
     *  true)
     */
    OrbitDeterminationManager(
            const SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
            const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const bool propagateOnCreation = true );

    OrbitDeterminationManager(
            const SystemOfBodies& bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
            const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const bool propagateOnCreation = true );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > getParametersToEstimate( )
    {
        return parametersToEstimate_;
    }

    SystemOfBodies getBodies( )
    {
        return bodies_;
    }

    //! Function to retrieve map of all observation managers
    /*!
     *  Function to retrieve map of all observation managers. A single observation manager can simulate observations and
     *  calculate observation partials for all link ends involved in the given observable type.
     *  \return Map of observation managers for all observable types involved in current orbit determination.
     */
    std::map< observation_models::ObservableType,
              std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >
    getObservationManagers( ) const
    {
        return observationManagers_;
    }

    //! Function to retrieve map of all observation simulators
    /*!
     *  Function to retrieve map of all observation simulators. A single observation simulators can simulate observations all
     *  link ends involved in the given observable type. The observation simulators are retrieved from the observation manager
     *  objects (that are stored in the observationManagers_ map).
     *  \return Map of observation simulators for all observable types involved in current orbit determination.
     */
    std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >
    getObservationSimulators( ) const;

    Eigen::MatrixXd normalizeAprioriCovariance( const Eigen::MatrixXd& inverseAPrioriCovariance,
                                                const Eigen::VectorXd& normalizationValues );

    Eigen::MatrixXd normalizeCovariance( const Eigen::MatrixXd& covariance, const Eigen::VectorXd& normalizationFactors );

    //! Function to normalize the matrix of partial derivatives so that each column is in the range [-1,1]
    /*!
     * Function to normalize the matrix of partial derivatives so that each column is in the range [-1,1]
     * \param observationMatrix Matrix of partial derivatives. Matrix modified by this function, and normalized matrix is
     * returned by reference
     * \return Vector with scaling values used for normalization
     */
    Eigen::VectorXd normalizeDesignMatrix( Eigen::MatrixXd& observationMatrix );

    void getNormalizedConsiderCovariance(
            const std::shared_ptr< CovarianceAnalysisInput< ObservationScalarType, TimeType > > estimationInput,
            const Eigen::VectorXd& considerNormalizationTerms,
            Eigen::MatrixXd& normalizedConsiderCovariance );

    std::shared_ptr< CovarianceAnalysisOutput< ObservationScalarType, TimeType > > computeCovariance(
            const std::shared_ptr< CovarianceAnalysisInput< ObservationScalarType, TimeType > > estimationInput );

    //! Function to perform parameter estimation from measurement data.
    /*!
     *  Function to perform parameter estimation, including orbit determination, i.e. body initial states, from measurement data.
     *  All observable types and link ends per obsevable types that are included in the measurement data input must have been
     *  provided to the constructor by the observationSettingsMap parameter.
     *  \param estimationInput Object containing all measurement data, associated metadata, including measurement weight, and a priori
     *  estimate for covariance matrix and parameter adjustment.
     *  \param convergenceChecker Object used to check convergence/termination of algorithm
     *  \return Object containing estimated parameter value and associateed data, such as residuals and observation partials.
     */
    std::shared_ptr< EstimationOutput< ObservationScalarType, TimeType > > estimateParameters(
            std::shared_ptr< EstimationInput< ObservationScalarType, TimeType > > estimationInput );

    template< int ObservationSize >
    void computePartialsAndObservations(
            const observation_models::LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings,
            const observation_models::ObservableType observableType,
            const observation_models::LinkEndType referenceLinkEnd,
            const std::vector< TimeType >& times,
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector,
            Eigen::MatrixXd& partials );

    //! Function to reset the current parameter estimate.
    /*!
     *  Function to reset the current parameter estimate; reintegrates the variational equations and equations of motion with new estimate.
     *  \param newParameterEstimate New estimate of parameter vector.
     *  \param reintegrateVariationalEquations Boolean denoting whether the variational equations are to be reintegrated
     */
    void resetParameterEstimate( const ParameterVectorType& newParameterEstimate, const bool reintegrateVariationalEquations = 1 );

    //! Function to retrieve the object to numerical integrate and update the variational equations and equations of motion
    /*!
     *  Function to retrieve the object to numerical integrate and update the variational equations and equations of motion
     *  \return Object to numerical integrate and update the variational equations and equations of motion
     */
    std::shared_ptr< propagators::VariationalEquationsSolver< ObservationScalarType, TimeType > > getVariationalEquationsSolver( ) const
    {
        return variationalEquationsSolver_;
    }

    //! Function to retrieve an observation manager
    /*!
     *  Function to retrieve an observation manager for a single observable type. The observation manager can simulate observations and
     *  calculate observation partials for all link ends involved in the given observable type.
     *  \param observableType Type of observable for which manager is to be retrieved.
     *  \return Observation manager for given observable type.
     */
    std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > getObservationManager(
            const observation_models::ObservableType observableType ) const;

    //! Function to retrieve the current paramater estimate.
    /*!
     *  Function to retrieve the current paramater estimate.
     *  \return Current paramater estimate.
     */
    ParameterVectorType getCurrentParameterEstimate( )
    {
        return currentParameterEstimate_;
    }

    //! Function to retrieve the object used to propagate/process the solution of the variational equations/dynamics.
    /*!
     *  Function to retrieve the object used to propagate/process the numerical solution of the variational.
     *  equations/dynamics.
     *  \return Object used to propagate/process the numerical solution of the variational equations/dynamics.
     */
    std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > getStateTransitionAndSensitivityMatrixInterface( )
    {
        return stateTransitionAndSensitivityMatrixInterface_;
    }

protected:
    //! Function called by either constructor to initialize the object.
    /*!
     *  Function called by either constructor to initialize the object.
     *  \param bodies Map of body objects with names of bodies, storing all environment models used in simulation.
     *  \param observationSettingsMap Sets of observation model settings per link ends (i.e. transmitter, receiver, etc.)
     *  for which measurement data is to be provided in orbit determination process
     *  (through estimateParameters function)
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param propagateOnCreation Boolean denoting whether initial propagatoon is to be performed upon object creation (default
     *  true)
     */
    void initializeOrbitDeterminationManager(
            const SystemOfBodies& bodies,
            const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const bool propagateOnCreation = true );

    std::pair< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
    performPreEstimationSteps( std::shared_ptr< CovarianceAnalysisInput< ObservationScalarType, TimeType > > estimationInput,
                               const ParameterVectorType& newParameterEstimate,
                               const bool calculateResiduals,
                               const int numberOfIterations,
                               bool& exceptionDuringPropagation,
                               std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > >& simulationResults );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > separateEstimatedAndConsiderDesignMatrices( const Eigen::MatrixXd& designMatrix,
                                                                                              const int numberObservations );

    //! Boolean to denote whether any dynamical parameters are estimated
    bool integrateAndEstimateOrbit_;

    //! Object used to propagate/process the numerical solution of the variational equations/dynamics
    std::shared_ptr< propagators::VariationalEquationsSolver< ObservationScalarType, TimeType > > variationalEquationsSolver_;

    //! List of object that compute the values/partials of the observables
    std::map< observation_models::ObservableType,
              std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >
            observationManagers_;

    //! Container object for all parameters that are to be estimated
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate_;

    //! Container object for consider parameters (if any)
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > considerParameters_;

    SystemOfBodies bodies_;

    //! Current values of the vector of estimated parameters
    ParameterVectorType currentParameterEstimate_;

    //! Consider parameters values
    ParameterVectorType considerParametersValues_;

    //! Total number of parameters (estimated and consider parameters together)
    unsigned int totalNumberParameters_;

    //! Number of estimated parameters
    unsigned int numberEstimatedParameters_;

    //! Number of consider parameters
    unsigned int numberConsiderParameters_;

    // std::vector< int > observationLinkParameterIndices_;

    //! Object used to interpolate the numerically integrated result of the state transition/sensitivity matrices.
    std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionAndSensitivityMatrixInterface_;

    //! Object used to interpolate the numerically integrated result of the dependent variables.
    std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface_;

    //! Boolean denoting whether consider parameters are included in the orbit determination
    bool considerParametersIncluded_;
};

#if TUDAT_BUILD_ALL_TESTS
extern template class OrbitDeterminationManager< double, double >;
#endif

}  // namespace simulation_setup

}  // namespace tudat

#include "tudat/simulation/estimation_setup/orbitDeterminationManagerConstructionImplementation.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManagerUtilitiesImplementation.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManagerCovarianceImplementation.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManagerEstimationImplementation.h"

#endif  // TUDAT_ORBITDETERMINATIONMANAGER_H
