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

    bool exceptionDuringPropagation = false;
    std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > simulationResults;

    auto createOutputFromDesignMatrices = [&]( auto estimatedParametersDesignMatrix,
                                               Eigen::MatrixXd designMatrixConsiderParameters )
            -> std::shared_ptr< CovarianceAnalysisOutput< ObservationScalarType, TimeType > >
    {
        typedef typename std::decay< decltype( estimatedParametersDesignMatrix ) >::type EstimatedDesignMatrixType;
        typedef typename linear_algebra::from_eigen< EstimatedDesignMatrixType >::traits EstimatedDesignMatrixTraits;

        // Normalise partials and inverse a priori covariance
        Eigen::VectorXd normalizationTerms = EstimatedDesignMatrixTraits::normalize_columns( estimatedParametersDesignMatrix );
        Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix =
                normalizeAprioriCovariance( estimationInput->getInverseOfAprioriCovariance( numberEstimatedParameters_ ), normalizationTerms );

        // Normalise partials w.r.t. consider parameters and consider covariance
        Eigen::VectorXd considerNormalizationTerms;
        Eigen::MatrixXd normalizedConsiderCovariance;
        if( considerParametersIncluded_ )
        {
            considerNormalizationTerms =
                    linear_algebra::MatrixTraits< double, linear_algebra::Dense >::normalize_columns( designMatrixConsiderParameters );
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

        EstimatedDesignMatrixType designMatrixForOutput = estimationInput->getSaveDesignMatrix( ) ?
                estimatedParametersDesignMatrix : EstimatedDesignMatrixType( 0, 0 );
        Eigen::MatrixXd considerDesignMatrixForOutput = estimationInput->getSaveDesignMatrix( ) ?
                designMatrixConsiderParameters : Eigen::MatrixXd::Zero( 0, 0 );

        return std::make_shared< CovarianceAnalysisOutput< ObservationScalarType, TimeType > >(
                designMatrixForOutput,
                estimationInput->getWeightsMatrixDiagonals( ),
                normalizationTerms,
                inverseNormalizedCovariance,
                considerDesignMatrixForOutput,
                considerNormalizationTerms,
                covarianceContributionConsiderParameters,
                estimationInput->getConsiderCovariance( ),
                exceptionDuringPropagation );
    };

    if( estimationInput->getUseSparseDesignMatrix( ) )
    {
        typedef typename linear_algebra::MatrixTraits< double, linear_algebra::Sparse >::matrix_type EstimatedDesignMatrixType;
        std::pair< std::pair< EstimatedDesignMatrixType, Eigen::MatrixXd >,
                   Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                designMatricesAndResiduals = performPreEstimationSteps< EstimatedDesignMatrixType >(
                        estimationInput, parameterValues, false, 0, exceptionDuringPropagation, simulationResults );
        return createOutputFromDesignMatrices( designMatricesAndResiduals.first.first, designMatricesAndResiduals.first.second );
    }
    else
    {
        typedef Eigen::MatrixXd EstimatedDesignMatrixType;
        std::pair< std::pair< EstimatedDesignMatrixType, Eigen::MatrixXd >,
                   Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                designMatricesAndResiduals = performPreEstimationSteps< EstimatedDesignMatrixType >(
                        estimationInput, parameterValues, false, 0, exceptionDuringPropagation, simulationResults );
        return createOutputFromDesignMatrices( designMatricesAndResiduals.first.first, designMatricesAndResiduals.first.second );
    }
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ORBITDETERMINATIONMANAGERCOVARIANCEIMPLEMENTATION_H
