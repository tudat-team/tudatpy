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
    const Eigen::VectorXd weightsMatrixDiagonal = estimationInput->getWeightsMatrixDiagonals( );
    const bool hasOffDiagonalWeights = estimationInput->hasOffDiagonalWeights( );
    Eigen::SparseMatrix< double > weightsMatrix;
    if( hasOffDiagonalWeights )
    {
        weightsMatrix = estimationInput->getWeightsMatrix( );
        if( weightsMatrix.rows( ) != totalNumberOfObservations || weightsMatrix.cols( ) != totalNumberOfObservations )
        {
            throw std::runtime_error( "Error when computing covariance, size of weights matrix (" +
                                      std::to_string( weightsMatrix.rows( ) ) + ", " + std::to_string( weightsMatrix.cols( ) ) +
                                      ") is not compatible with number of observations (" + std::to_string( totalNumberOfObservations ) +
                                      ")" );
        }
    }
    else if( weightsMatrixDiagonal.rows( ) != totalNumberOfObservations )
    {
        throw std::runtime_error( "Error when computing covariance, size of weights diagonal (" +
                                  std::to_string( weightsMatrixDiagonal.rows( ) ) +
                                  ") is not compatible with number of observations (" + std::to_string( totalNumberOfObservations ) +
                                  ")" );
    }

    if( numberEstimatedParameters_ > static_cast< unsigned int >( totalNumberOfObservations ) &&
        estimationInput->getInverseOfAprioriCovariance( ).rows( ) == 0 )
    {
        std::cerr << "Warning when computing covariance, number of observations is smaller than number of estimated parameters, and no "
                     "a priori information is provided."
                  << std::endl;
    }

    // Define full parameters values
    ParameterVectorType parameterValues = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );

    // Compute design matrices (estimated and consider), and residuals (empty for covariance analysis)
    bool exceptionDuringPropagation = false;
    std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > simulationResults;
    std::pair< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
            designMatricesAndResiduals =
                    performPreEstimationSteps( estimationInput, parameterValues, false, 0, exceptionDuringPropagation, simulationResults );
    Eigen::MatrixXd designMatrixEstimatedParameters = designMatricesAndResiduals.first.first;
    Eigen::MatrixXd designMatrixConsiderParameters;
    designMatrixConsiderParameters = designMatricesAndResiduals.first.second;

    // Normalise partials and inverse a priori covariance
    Eigen::VectorXd normalizationTerms = normalizeDesignMatrix( designMatrixEstimatedParameters );
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
    Eigen::MatrixXd inverseNormalizedCovariance;
    if( hasOffDiagonalWeights )
    {
        inverseNormalizedCovariance = linear_algebra::calculateInverseOfUpdatedCovarianceMatrix(
                designMatrixEstimatedParameters.block( 0, 0, designMatrixEstimatedParameters.rows( ), numberEstimatedParameters_ ),
                weightsMatrix,
                normalizedInverseAprioriCovarianceMatrix,
                constraintStateMultiplier,
                constraintRightHandSide,
                estimationInput->getLimitConditionNumberForWarning( ) );
    }
    else
    {
        inverseNormalizedCovariance = linear_algebra::calculateInverseOfUpdatedCovarianceMatrix(
                designMatrixEstimatedParameters.block( 0, 0, designMatrixEstimatedParameters.rows( ), numberEstimatedParameters_ ),
                weightsMatrixDiagonal,
                normalizedInverseAprioriCovarianceMatrix,
                constraintStateMultiplier,
                constraintRightHandSide,
                estimationInput->getLimitConditionNumberForWarning( ) );
    }

    // Compute contribution consider parameters
    Eigen::MatrixXd covarianceContributionConsiderParameters;
    if( considerParametersIncluded_ )
    {
        if( hasOffDiagonalWeights )
        {
            covarianceContributionConsiderParameters = linear_algebra::calculateConsiderParametersCovarianceContribution(
                    inverseNormalizedCovariance.inverse( ),
                    designMatrixEstimatedParameters,
                    weightsMatrix,
                    designMatrixConsiderParameters,
                    normalizedConsiderCovariance );
        }
        else
        {
            covarianceContributionConsiderParameters = linear_algebra::calculateConsiderParametersCovarianceContribution(
                    inverseNormalizedCovariance.inverse( ),
                    designMatrixEstimatedParameters,
                    weightsMatrixDiagonal,
                    designMatrixConsiderParameters,
                    normalizedConsiderCovariance );
        }
    }
    else
    {
        covarianceContributionConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
    }

    // Create covariance output object
    std::shared_ptr< CovarianceAnalysisOutput< ObservationScalarType, TimeType > > estimationOutput =
            std::make_shared< CovarianceAnalysisOutput< ObservationScalarType, TimeType > >(
                    estimationInput->getSaveDesignMatrix( ) ? designMatrixEstimatedParameters : Eigen::MatrixXd::Zero( 0, 0 ),
                    weightsMatrixDiagonal,
                    normalizationTerms,
                    inverseNormalizedCovariance,
                    estimationInput->getSaveDesignMatrix( ) ? designMatrixConsiderParameters : Eigen::MatrixXd::Zero( 0, 0 ),
                    considerNormalizationTerms,
                    covarianceContributionConsiderParameters,
                    estimationInput->getConsiderCovariance( ),
                    exceptionDuringPropagation );

    return estimationOutput;
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ORBITDETERMINATIONMANAGERCOVARIANCEIMPLEMENTATION_H
