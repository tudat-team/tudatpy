/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ORBITDETERMINATIONMANAGERCOVARIANCEIMPLEMENTATION_H
#define TUDAT_ORBITDETERMINATIONMANAGERCOVARIANCEIMPLEMENTATION_H

#include <iostream>

#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::shared_ptr< CovarianceAnalysisOutput< ObservationScalarType, TimeType > >
OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::computeCovariance(
        const std::shared_ptr< CovarianceAnalysisInput< ObservationScalarType, TimeType > > estimationInput )
{
    // Get total number of observations
    int totalNumberOfObservations = estimationInput->getObservationCollection( )->getTotalObservableSize( );

    if( numberEstimatedParameters_ > static_cast< unsigned int >( totalNumberOfObservations ) &&
        estimationInput->getInverseOfAprioriCovariance( ).rows( ) == 0 )
    {
        std::cerr << "Warning when computing covariance, number of observations is smaller than number of estimated parameters, and no "
                     "a priori information is provided."
                  << std::endl;
    }

    // Define full parameters values
    ParameterVectorType parameterValues = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );

    if( estimationInput->getSaveDesignMatrix( ) )
    {
        throw std::runtime_error( "Error when computing covariance with sparse design matrix, saving the dense design matrix is disabled." );
    }

    // Compute sparse estimated-parameter design matrix through pre-estimation steps.
    bool exceptionDuringPropagation = false;
    std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > simulationResults;
    std::pair< std::pair< Eigen::SparseMatrix< double >, Eigen::MatrixXd >,
               Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
            designMatricesAndResiduals = performPreEstimationSteps< Eigen::SparseMatrix< double > >(
                    estimationInput, parameterValues, false, 0, exceptionDuringPropagation, simulationResults );
    Eigen::SparseMatrix< double > estimatedParametersDesignMatrix = designMatricesAndResiduals.first.first;
    Eigen::MatrixXd designMatrixConsiderParameters = designMatricesAndResiduals.first.second;

    // Normalise partials and inverse a priori covariance
    Eigen::VectorXd normalizationTerms = normalizeSparseDesignMatrix( estimatedParametersDesignMatrix );
    Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix =
            normalizeAprioriCovariance( estimationInput->getInverseOfAprioriCovariance( numberEstimatedParameters_ ), normalizationTerms );

    // Normalise partials w.r.t. consider parameters and consider covariance
    Eigen::VectorXd considerNormalizationTerms;
    Eigen::MatrixXd normalizedConsiderCovariance;
    if( considerParametersIncluded_ )
    {
        considerNormalizationTerms = normalizeDesignMatrix( designMatrixConsiderParameters );
        getNormalizedConsiderCovariance( estimationInput, considerNormalizationTerms, normalizedConsiderCovariance );
    }
    else
    {
        considerNormalizationTerms = Eigen::VectorXd::Zero( 0 );
        normalizedConsiderCovariance = Eigen::MatrixXd::Zero( 0, 0 );
    }

    // Retrieve constraints
    Eigen::MatrixXd constraintStateMultiplier;
    Eigen::VectorXd constraintRightHandSide;
    parametersToEstimate_->getConstraints( constraintStateMultiplier, constraintRightHandSide );

    // Compute inverse of updated covariance
    std::optional< Eigen::MatrixXd > constraintStateMultiplierOptional =
            constraintStateMultiplier.rows( ) == 0 ? std::nullopt : std::optional< Eigen::MatrixXd >( constraintStateMultiplier );
    std::optional< Eigen::VectorXd > constraintRightHandSideOptional =
            constraintStateMultiplier.rows( ) == 0 ? std::nullopt : std::optional< Eigen::VectorXd >( constraintRightHandSide );
    Eigen::MatrixXd inverseNormalizedCovariance = linear_algebra::calculateInverseOfUpdatedCovarianceMatrix(
            estimatedParametersDesignMatrix,
            estimationInput->getWeightsMatrixDiagonals( ),
            normalizedInverseAprioriCovarianceMatrix,
            constraintStateMultiplierOptional,
            constraintRightHandSideOptional,
            estimationInput->getLimitConditionNumberForWarning( ) );

    // Compute contribution consider parameters
    Eigen::MatrixXd covarianceContributionConsiderParameters;
    if( considerParametersIncluded_ )
    {
        Eigen::MatrixXd normalizedCovariance = inverseNormalizedCovariance.inverse( );
        covarianceContributionConsiderParameters =
                linear_algebra::calculateConsiderParametersCovarianceContribution( normalizedCovariance,
                                                                                   estimatedParametersDesignMatrix,
                                                                                   estimationInput->getWeightsMatrixDiagonals( ),
                                                                                   designMatrixConsiderParameters,
                                                                                   normalizedConsiderCovariance );
    }
    else
    {
        covarianceContributionConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
    }

    // Create covariance output object
    std::shared_ptr< CovarianceAnalysisOutput< ObservationScalarType, TimeType > > estimationOutput =
            std::make_shared< CovarianceAnalysisOutput< ObservationScalarType, TimeType > >(
                    Eigen::MatrixXd::Zero( 0, 0 ),
                    estimationInput->getWeightsMatrixDiagonals( ),
                    normalizationTerms,
                    inverseNormalizedCovariance,
                    Eigen::MatrixXd::Zero( 0, 0 ),
                    considerNormalizationTerms,
                    covarianceContributionConsiderParameters,
                    estimationInput->getConsiderCovariance( ),
                    exceptionDuringPropagation );

    return estimationOutput;
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ORBITDETERMINATIONMANAGERCOVARIANCEIMPLEMENTATION_H
