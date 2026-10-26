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
#include <optional>
#include <type_traits>
#include <utility>

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

    // Compute inverse of updated covariance in the storage of the given design matrix (dense or sparse). The other
    // matrices are converted to the same storage, the inverse covariance is returned as a dense matrix.
    std::optional< Eigen::VectorXd > constraintRightHandSideOptional =
            constraintStateMultiplier.rows( ) == 0 ? std::nullopt : std::optional< Eigen::VectorXd >( constraintRightHandSide );
    auto computeInverseCovariance = [ & ]( const auto& designMatrix ) {
        using MatrixType = std::decay_t< decltype( designMatrix ) >;
        using Traits = typename linear_algebra::from_eigen< MatrixType >::traits;
        std::optional< MatrixType > constraintStateMultiplierOptional = constraintStateMultiplier.rows( ) == 0
                ? std::nullopt
                : std::optional< MatrixType >( Traits::from_dense( constraintStateMultiplier ) );
        return Traits::to_dense( linear_algebra::calculateInverseOfUpdatedCovarianceMatrix< MatrixType >(
                designMatrix,
                estimationInput->getWeightsMatrixDiagonals( ),
                Traits::from_dense( normalizedInverseAprioriCovarianceMatrix ),
                constraintStateMultiplierOptional,
                constraintRightHandSideOptional,
                estimationInput->getLimitConditionNumberForWarning( ) ) );
    };

    // Design matrix to store in the output; left empty if the design matrix is not saved
    optionalAnyMatrix savedDesignMatrixEstimatedParameters;
    Eigen::MatrixXd inverseNormalizedCovariance;
    if( estimationInput->getSparseDesignMatrix( ) )
    {
        Eigen::SparseMatrix< double > sparseDesignMatrixEstimatedParameters = designMatrixEstimatedParameters.sparseView( );
        inverseNormalizedCovariance = computeInverseCovariance( sparseDesignMatrixEstimatedParameters );
        if( estimationInput->getSaveDesignMatrix( ) )
        {
            savedDesignMatrixEstimatedParameters = std::move( sparseDesignMatrixEstimatedParameters );
        }
    }
    else
    {
        inverseNormalizedCovariance = computeInverseCovariance( designMatrixEstimatedParameters );
        if( estimationInput->getSaveDesignMatrix( ) )
        {
            savedDesignMatrixEstimatedParameters = designMatrixEstimatedParameters;
        }
    }

    // Compute contribution consider parameters
    Eigen::MatrixXd covarianceContributionConsiderParameters;
    if( considerParametersIncluded_ )
    {
        Eigen::MatrixXd normalizedCovariance = inverseNormalizedCovariance.inverse( );
        covarianceContributionConsiderParameters =
                linear_algebra::calculateConsiderParametersCovarianceContribution( normalizedCovariance,
                                                                                   designMatrixEstimatedParameters,
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
                    savedDesignMatrixEstimatedParameters,
                    estimationInput->getWeightsMatrixDiagonals( ),
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
