/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ORBITDETERMINATIONMANAGERESTIMATIONIMPLEMENTATION_H
#define TUDAT_ORBITDETERMINATIONMANAGERESTIMATIONIMPLEMENTATION_H

#include <iostream>
#include <stdexcept>

#include <Eigen/SparseCore>

#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManagerHelpers.h"
#include "tudat/simulation/estimation_setup/outlierRejection.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::shared_ptr< EstimationOutput< ObservationScalarType, TimeType > >
OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::estimateParameters(
        std::shared_ptr< EstimationInput< ObservationScalarType, TimeType > > estimationInput )

{
    currentParameterEstimate_ = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );

    // Create the object that rejects and recovers outlying observations during the estimation. A null pointer is
    // returned (and no outlier rejection is performed) when the user did not provide outlier rejection settings.
    const std::shared_ptr< OutlierRejection< ObservationScalarType, TimeType > > outlierRejection =
            createOutlierRejection< ObservationScalarType, TimeType >( estimationInput->getOutlierRejectionSettings( ),
                                                                       estimationInput->getObservationDataset( ) );
    const bool applyOutlierRejection = ( outlierRejection != nullptr );

    // Flattened data of the observations that are used in the estimation (that is, all observations that are not
    // rejected). When outlier rejection is used, this data is recreated at the start of every iteration, since
    // observations may have been rejected or recovered in the previous iteration.
    observation_models::FlattenedObservationData< ObservationScalarType, TimeType > estimationData =
            estimationInput->getObservationDataset( )->createOrderedFlattenedObservationData( false );
    int totalNumberOfObservations = static_cast< int >( estimationData.getObservationVector( ).size( ) );

    // Flattened data of all observations, including the rejected ones. Outlier rejection algorithms must be able to
    // evaluate their criteria for rejected observations as well, since those observations may have to be recovered.
    // Rejecting an observation only changes whether it is active, and never the structure of the dataset, so this
    // data (and the observation covariance derived from it) is created only once.
    observation_models::FlattenedObservationData< ObservationScalarType, TimeType > computationData;
    Eigen::MatrixXd observationCovariance;
    if( applyOutlierRejection )
    {
        computationData = estimationInput->getObservationDataset( )->createOrderedFlattenedObservationData( true );

        // The observation weights represent the inverse of the observation covariance, and do not change during the
        // estimation, so the covariance is computed once here.
        observationCovariance = Eigen::MatrixXd( computationData.getSparseWeightMatrix( ) ).inverse( );
    }

    if( numberEstimatedParameters_ > static_cast< unsigned int >( totalNumberOfObservations ) &&
        estimationInput->getInverseOfAprioriCovariance( ).rows( ) == 0 )
    {
        std::cerr << "Warning when estimating parameters, number of observations is smaller than number of estimated parameters, and "
                     "no a priori information is provided."
                  << std::endl;
    }

    Eigen::VectorXd weightsMatrixDiagonals;
    Eigen::SparseMatrix< double > weightsMatrix;
    bool hasOffDiagonalWeights = false;
    retrieveObservationWeights< ObservationScalarType, TimeType >(
            estimationData, weightsMatrixDiagonals, weightsMatrix, hasOffDiagonalWeights );

    // Declare variables to be returned (i.e. results from best iteration)
    double bestCostFunction = TUDAT_NAN;
    double bestRmsResidual = TUDAT_NAN;
    ParameterVectorType bestParameterEstimate = ParameterVectorType::Constant( numberEstimatedParameters_, TUDAT_NAN );
    Eigen::VectorXd bestTransformationData = Eigen::VectorXd::Constant( numberEstimatedParameters_, TUDAT_NAN );
    Eigen::VectorXd bestResiduals = Eigen::VectorXd::Constant( totalNumberOfObservations, TUDAT_NAN );
    Eigen::MatrixXd bestDesignMatrixEstimatedParameters = Eigen::MatrixXd::Zero( 0, 0 );
    Eigen::VectorXd bestWeightsMatrixDiagonal = Eigen::VectorXd::Constant( totalNumberOfObservations, TUDAT_NAN );
    Eigen::MatrixXd bestInverseNormalizedCovarianceMatrix =
            Eigen::MatrixXd::Constant( numberEstimatedParameters_, numberEstimatedParameters_, TUDAT_NAN );

    Eigen::VectorXd bestConsiderTransformationData;
    Eigen::MatrixXd bestDesignMatrixConsiderParameters, bestConsiderCovarianceContribution;
    if( considerParametersIncluded_ )
    {
        bestConsiderTransformationData = Eigen::VectorXd::Constant( numberConsiderParameters_, TUDAT_NAN );
        bestDesignMatrixConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
        bestConsiderCovarianceContribution = Eigen::MatrixXd::Constant( numberEstimatedParameters_, numberEstimatedParameters_, TUDAT_NAN );
    }
    else
    {
        bestConsiderTransformationData = Eigen::VectorXd::Zero( 0 );
        bestDesignMatrixConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
        bestConsiderCovarianceContribution = Eigen::MatrixXd::Zero( 0, 0 );
    }

    std::vector< Eigen::VectorXd > residualHistory;
    std::vector< ParameterVectorType > parameterHistory;
    std::vector< std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > > simulationResultsPerIteration;

    // Declare residual bookkeeping variables
    std::vector< double > rmsResidualHistory;
    std::vector< double > costFunctionHistory;
    double residualRms;
    std::map< observation_models::ObservableType, double > residualRmsPerType;
    double costFunction;

    // Set current parameter estimate as both previous and current estimate
    ParameterVectorType newParameterEstimate = currentParameterEstimate_;
    ParameterVectorType oldParameterEstimate = currentParameterEstimate_;

    bool exceptionDuringPropagation = false, exceptionDuringInversion = false;

    // Iterate until convergence (at least once)
    int bestIteration = -1;
    int numberOfIterations = 0;
    while( true )
    {
        oldParameterEstimate = newParameterEstimate;

        // Update the data of the observations that are used in the estimation: observations may have been rejected or
        // recovered at the end of the previous iteration.
        if( applyOutlierRejection && numberOfIterations > 0 )
        {
            estimationData = estimationInput->getObservationDataset( )->createOrderedFlattenedObservationData( false );
            totalNumberOfObservations = static_cast< int >( estimationData.getObservationVector( ).size( ) );
            if( totalNumberOfObservations == 0 )
            {
                throw std::runtime_error( "Error during parameter estimation, all observations have been rejected by the outlier "
                                          "rejection algorithm." );
            }
            retrieveObservationWeights< ObservationScalarType, TimeType >(
                    estimationData, weightsMatrixDiagonals, weightsMatrix, hasOffDiagonalWeights );
        }

        // Compute design matrices (for estimated and consider parameters) and residuals. When outlier rejection is
        // used, these are computed for all observations (including the rejected ones), and the rows of the
        // observations that are used in the estimation are extracted from them afterwards.
        std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > simulationResults;
        std::pair< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                designMatricesAndResiduals = performPreEstimationSteps( estimationInput,
                                                                        newParameterEstimate,
                                                                        applyOutlierRejection ? computationData : estimationData,
                                                                        true,
                                                                        numberOfIterations,
                                                                        exceptionDuringPropagation,
                                                                        simulationResults );

        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > computedResiduals = std::move( designMatricesAndResiduals.second );
        Eigen::MatrixXd computedDesignMatrixEstimatedParameters = std::move( designMatricesAndResiduals.first.first );
        Eigen::MatrixXd computedDesignMatrixConsiderParameters = std::move( designMatricesAndResiduals.first.second );

        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residuals;
        Eigen::MatrixXd designMatrixEstimatedParameters;
        Eigen::MatrixXd designMatrixConsiderParameters;
        if( applyOutlierRejection )
        {
            residuals = extractEstimationObservationRows( computationData, estimationData, computedResiduals );
            designMatrixEstimatedParameters =
                    extractEstimationObservationRows( computationData, estimationData, computedDesignMatrixEstimatedParameters );
            designMatrixConsiderParameters = considerParametersIncluded_
                    ? extractEstimationObservationRows( computationData, estimationData, computedDesignMatrixConsiderParameters )
                    : computedDesignMatrixConsiderParameters;
        }
        else
        {
            residuals = std::move( computedResiduals );
            designMatrixEstimatedParameters = std::move( computedDesignMatrixEstimatedParameters );
            designMatrixConsiderParameters = std::move( computedDesignMatrixConsiderParameters );
        }

        // Set simulation results
        if( estimationInput->getSaveStateHistoryForEachIteration( ) )
        {
            simulationResultsPerIteration.push_back( simulationResults->clone( ) );
        }

        // Normalise estimated parameters partials and inverse apriori covariance
        Eigen::VectorXd normalizationTerms = normalizeDesignMatrix( designMatrixEstimatedParameters );
        Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix = normalizeAprioriCovariance(
                estimationInput->getInverseOfAprioriCovariance( numberEstimatedParameters_ ), normalizationTerms );

        // Normalise partials w.r.t. consider parameters, consider covariance and parameters deviations
        Eigen::VectorXd normalizationTermsConsider, normalizedConsiderParametersDeviation;
        Eigen::MatrixXd normalizedConsiderCovariance;
        if( considerParametersIncluded_ )
        {
            normalizationTermsConsider = normalizeDesignMatrix( designMatrixConsiderParameters );
            getNormalizedConsiderCovariance( estimationInput, normalizationTermsConsider, normalizedConsiderCovariance );
            if( estimationInput->considerParametersDeviations_.rows( ) == 0 )
            {
                normalizedConsiderParametersDeviation = Eigen::VectorXd( normalizationTermsConsider.rows( ) );
            }
            else
            {
                normalizedConsiderParametersDeviation =
                        estimationInput->considerParametersDeviations_.cwiseProduct( normalizationTermsConsider );
            }
        }
        else
        {
            normalizationTermsConsider = Eigen::VectorXd::Zero( 0 );
            normalizedConsiderCovariance = Eigen::MatrixXd::Zero( 0, 0 );
            normalizedConsiderParametersDeviation = Eigen::VectorXd::Zero( 0 );
        }

        // Perform least squares calculation for correction to parameter vector.
        std::pair< Eigen::VectorXd, Eigen::MatrixXd > leastSquaresOutput;
        try
        {
            // Get constraints
            Eigen::MatrixXd constraintStateMultiplier;
            Eigen::VectorXd constraintRightHandSide;
            parametersToEstimate_->getConstraints( constraintStateMultiplier, constraintRightHandSide );

            double conditionNumberCheck = estimationInput->getLimitConditionNumberForWarning( );
            if( numberOfIterations > 0 && estimationInput->conditionNumberWarningEachIteration_ == false )
            {
                conditionNumberCheck = TUDAT_NAN;
            }
            // Perform LSQ inversion
            if( hasOffDiagonalWeights )
            {
                leastSquaresOutput =
                        std::move( linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrixEstimatedParameters,
                                                                                                  residuals.template cast< double >( ),
                                                                                                  weightsMatrix,
                                                                                                  normalizedInverseAprioriCovarianceMatrix,
                                                                                                  conditionNumberCheck,
                                                                                                  constraintStateMultiplier,
                                                                                                  constraintRightHandSide,
                                                                                                  designMatrixConsiderParameters,
                                                                                                  normalizedConsiderParametersDeviation ) );
            }
            else
            {
                leastSquaresOutput =
                        std::move( linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrixEstimatedParameters,
                                                                                                  residuals.template cast< double >( ),
                                                                                                  weightsMatrixDiagonals,
                                                                                                  normalizedInverseAprioriCovarianceMatrix,
                                                                                                  conditionNumberCheck,
                                                                                                  constraintStateMultiplier,
                                                                                                  constraintRightHandSide,
                                                                                                  designMatrixConsiderParameters,
                                                                                                  normalizedConsiderParametersDeviation ) );
            }

            if( constraintStateMultiplier.rows( ) > 0 )
            {
                leastSquaresOutput.first.conservativeResize( numberEstimatedParameters_ );
            }
        }
        catch( std::runtime_error& error )
        {
            std::cerr << "Error when solving normal equations during parameter estimation: " << std::endl
                      << error.what( ) << std::endl
                      << "Terminating estimation" << std::endl;
            exceptionDuringInversion = true;
            break;
        }

        ParameterVectorType parameterAddition =
                ( leastSquaresOutput.first.cwiseQuotient( normalizationTerms.segment( 0, numberEstimatedParameters_ ) ) )
                        .template cast< ObservationScalarType >( );

        // Compute the (unnormalized) covariance of the estimated parameters, which is only required by outlier
        // rejection algorithms. It is computed here, since the normalized inverse covariance and the normalization
        // terms are moved into the results of the best iteration further down.
        Eigen::MatrixXd parameterCovariance;
        if( applyOutlierRejection )
        {
            parameterCovariance = normaliseUnnormaliseCovarianceMatrix( leastSquaresOutput.second.inverse( ),
                                                                        normalizationTerms.segment( 0, numberEstimatedParameters_ ),
                                                                        false );
        }

        // Compute contribution consider parameters
        Eigen::MatrixXd covarianceContributionConsiderParameters;
        if( considerParametersIncluded_ )
        {
            if( hasOffDiagonalWeights )
            {
                covarianceContributionConsiderParameters =
                        linear_algebra::calculateConsiderParametersCovarianceContribution( ( leastSquaresOutput.second ).inverse( ),
                                                                                           designMatrixEstimatedParameters,
                                                                                           weightsMatrix,
                                                                                           designMatrixConsiderParameters,
                                                                                           normalizedConsiderCovariance );
            }
            else
            {
                covarianceContributionConsiderParameters =
                        linear_algebra::calculateConsiderParametersCovarianceContribution( ( leastSquaresOutput.second ).inverse( ),
                                                                                           designMatrixEstimatedParameters,
                                                                                           weightsMatrixDiagonals,
                                                                                           designMatrixConsiderParameters,
                                                                                           normalizedConsiderCovariance );
            }
        }
        else
        {
            covarianceContributionConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
        }

        // Calculate mean residual for current iteration.
        residualRms = linear_algebra::getVectorEntryRootMeanSquare( residuals.template cast< double >( ) );
        if( hasOffDiagonalWeights )
        {
            costFunction =
                    linear_algebra::computeLeastSquaresCostFunctionFromFullWeights( weightsMatrix, residuals.template cast< double >( ) );
        }
        else
        {
            costFunction = linear_algebra::computeLeastSquaresCostFunction( weightsMatrixDiagonals, residuals.template cast< double >( ) );
        }
        rmsResidualHistory.push_back( residualRms );
        costFunctionHistory.push_back( costFunction );

        if( estimationInput->getSaveResidualsAndParametersFromEachIteration( ) )
        {
            residualHistory.push_back( residuals.template cast< double >( ) );
            if( numberOfIterations == 0 )
            {
                parameterHistory.push_back( oldParameterEstimate );
            }
        }

        if( estimationInput->getPrintOutput( ) )
        {
            std::map< observation_models::ObservableType, std::vector< double > > residualsPerObservableType;
            for( int row = 0; row < residuals.size( ); ++row )
            {
                const unsigned int setId = estimationData.getSetIds( ).at( row );
                const observation_models::ObservableType observableType =
                        estimationInput->getObservationDataset( )->getObservationSetMetadata( setId ).observableType_;
                residualsPerObservableType[ observableType ].push_back( static_cast< double >( residuals( row ) ) );
            }

            if( residualsPerObservableType.size( ) > 1 )
            {
                for( const auto& residualsForType : residualsPerObservableType )
                {
                    Eigen::VectorXd currentTypeResiduals = Eigen::VectorXd::Zero( static_cast< int >( residualsForType.second.size( ) ) );
                    for( int i = 0; i < currentTypeResiduals.size( ); ++i )
                    {
                        currentTypeResiduals( i ) = residualsForType.second.at( i );
                    }

                    const double currentTypeRms = linear_algebra::getVectorEntryRootMeanSquare( currentTypeResiduals );
                    residualRmsPerType[ residualsForType.first ] = currentTypeRms;
                    std::cout << "Current residual for observable (" << observation_models::getObservableName( residualsForType.first )
                              << "): " << currentTypeRms << std::endl;
                }
            }
            else
            {
                std::cout << "Current residual: " << residualRms << std::endl;
            }
        }

        // If current iteration is better than previous one, update 'best' data.
        if( costFunction < bestCostFunction || !( bestCostFunction == bestCostFunction ) )
        {
            bestCostFunction = costFunction;
            bestRmsResidual = residualRms;
            bestParameterEstimate = oldParameterEstimate;
            bestResiduals = std::move( residuals.template cast< double >( ) );

            if( applyOutlierRejection )
            {
                // Residuals were already computed for all observations, rejected ones included.
                estimationInput->getObservationDataset( )->setResidualVector( computationData, computedResiduals );
            }
            else
            {
                const observation_models::FlattenedObservationData< ObservationScalarType, TimeType > residualComputationData =
                        estimationInput->getObservationDataset( )->createComputationFlattenedObservationData( true );
                if( residualComputationData.getObservationVector( ).size( ) == estimationData.getObservationVector( ).size( ) )
                {
                    estimationInput->getObservationDataset( )->setResidualVector( estimationData, residuals );
                }
                else
                {
                    Eigen::MatrixXd unusedDesignMatrix;
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualComputationResiduals;
                    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >( estimationInput->getObservationDataset( ),
                                                                                          residualComputationData,
                                                                                          observationManagers_,
                                                                                          totalNumberParameters_,
                                                                                          unusedDesignMatrix,
                                                                                          residualComputationResiduals,
                                                                                          true,
                                                                                          false );
                    estimationInput->getObservationDataset( )->setResidualVector( residualComputationData, residualComputationResiduals );
                }
            }
            if( estimationInput->getSaveDesignMatrix( ) )
            {
                bestDesignMatrixEstimatedParameters = std::move( designMatrixEstimatedParameters );
                bestDesignMatrixConsiderParameters = std::move( designMatrixConsiderParameters );
            }
            bestWeightsMatrixDiagonal = weightsMatrixDiagonals;
            bestTransformationData = std::move( normalizationTerms );
            bestInverseNormalizedCovarianceMatrix = std::move( leastSquaresOutput.second );
            bestIteration = numberOfIterations;
            bestConsiderTransformationData = std::move( normalizationTermsConsider );
            bestConsiderCovarianceContribution = covarianceContributionConsiderParameters;
        }

        // Update which observations are rejected, using the data of the current iteration. The observations that are
        // rejected here are excluded from the next iteration, and rejected observations that are recovered here are
        // included again.
        if( applyOutlierRejection )
        {
            const OutlierRejectionInput< ObservationScalarType, TimeType > outlierRejectionInput( numberOfIterations,
                                                                                                  computationData,
                                                                                                  observationCovariance,
                                                                                                  computedResiduals,
                                                                                                  computedDesignMatrixEstimatedParameters,
                                                                                                  parameterCovariance,
                                                                                                  parameterAddition );
            outlierRejection->updateRejectionStatus( outlierRejectionInput );

            if( estimationInput->getPrintOutput( ) )
            {
                std::cout << "Number of rejected observations: " << outlierRejection->getNumberOfRejectedObservations( ) << " out of "
                          << estimationInput->getObservationDataset( )->getNumberOfObservations( ) << std::endl;
            }
        }

        // Increment number of iterations
        numberOfIterations++;

        // Check for convergence
        bool applyParameterCorrection = true;
        bool terminateLoop = false;
        if( estimationInput->getConvergenceChecker( )->isEstimationConverged( numberOfIterations, rmsResidualHistory ) )
        {
            terminateLoop = true;
            applyParameterCorrection = estimationInput->applyFinalParameterCorrection_;
        }

        if( applyParameterCorrection )
        {
            // Update value of parameter vector
            newParameterEstimate = oldParameterEstimate + parameterAddition;
            parametersToEstimate_->template resetParameterValues< ObservationScalarType >( newParameterEstimate );
            newParameterEstimate = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );

            if( estimationInput->getSaveResidualsAndParametersFromEachIteration( ) )
            {
                parameterHistory.push_back( newParameterEstimate );
            }

            if( estimationInput->getPrintOutput( ) )
            {
                std::cout << "Parameter update" << parameterAddition.transpose( ) << std::endl;
            }
        }

        if( terminateLoop )
        {
            break;
        }
    }

    if( estimationInput->getPrintOutput( ) )
    {
        std::cout << "Best iteration: " << bestIteration << " out of " << numberOfIterations - 1 << std::endl;
        std::cout << "Rms residual from best iteration: " << bestRmsResidual << std::endl;
        std::cout << "Cost function from best iteration: " << bestCostFunction << std::endl;
    }

    // Create estimation output object
    std::shared_ptr< EstimationOutput< ObservationScalarType, TimeType > > estimationOutput =
            std::make_shared< EstimationOutput< ObservationScalarType, TimeType > >(
                    bestParameterEstimate,
                    bestResiduals,
                    bestDesignMatrixEstimatedParameters,
                    bestWeightsMatrixDiagonal,
                    hasOffDiagonalWeights ? weightsMatrix : Eigen::SparseMatrix< double >( ),
                    bestTransformationData,
                    bestInverseNormalizedCovarianceMatrix,
                    bestRmsResidual,
                    bestIteration,
                    residualHistory,
                    parameterHistory,
                    bestDesignMatrixConsiderParameters,
                    bestConsiderTransformationData,
                    bestConsiderCovarianceContribution,
                    estimationInput->getConsiderCovariance( ),
                    exceptionDuringInversion,
                    exceptionDuringPropagation );

    if( estimationInput->getSaveStateHistoryForEachIteration( ) )
    {
        estimationOutput->setSimulationResults( simulationResultsPerIteration );
    }

    return estimationOutput;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
template< int ObservationSize >
void OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::computePartialsAndObservations(
        const observation_models::LinkEnds& linkEnds,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< TimeType >& times,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector,
        Eigen::MatrixXd& partials )
{
    std::shared_ptr< observation_models::ObservationManager< ObservationSize, ObservationScalarType, TimeType > > observationManager =
            std::dynamic_pointer_cast< observation_models::ObservationManager< ObservationSize, ObservationScalarType, TimeType > >(
                    observationManagers_.at( observableType ) );

    // Extract observation model
    std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel =
            observationManager->getObservationModel( linkEnds );

    // Compute analytical partials
    observationManager->computeObservationsWithPartials(
            times, linkEnds, referenceLinkEnd, ancillarySettings, observationsVector, partials, true, true );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::resetParameterEstimate(
        const ParameterVectorType& newParameterEstimate,
        const bool reintegrateVariationalEquations )
{
    parametersToEstimate_->template resetParameterValues< ObservationScalarType >( newParameterEstimate );
    if( integrateAndEstimateOrbit_ )
    {
        variationalEquationsSolver_->resetParameterEstimate(
                parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( ), reintegrateVariationalEquations );
    }

    currentParameterEstimate_ = newParameterEstimate;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::pair< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::performPreEstimationSteps(
        std::shared_ptr< CovarianceAnalysisInput< ObservationScalarType, TimeType > > estimationInput,
        const ParameterVectorType& newParameterEstimate,
        const observation_models::FlattenedObservationData< ObservationScalarType, TimeType >& flattenedObservationData,
        const bool calculateResiduals,
        const int numberOfIterations,
        bool& exceptionDuringPropagation,
        std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > >& simulationResults )
{
    const int totalNumberOfObservations = static_cast< int >( flattenedObservationData.getObservationVector( ).size( ) );

    // Re-integrate equations of motion and variational equations with new parameter estimate.
    try
    {
        if( ( numberOfIterations > 0 ) || ( estimationInput->getReintegrateEquationsOnFirstIteration( ) ) )
        {
            resetParameterEstimate( newParameterEstimate, estimationInput->getReintegrateVariationalEquations( ) );
        }

        if( std::dynamic_pointer_cast< EstimationInput< ObservationScalarType, TimeType > >( estimationInput ) != nullptr )
        {
            if( std::dynamic_pointer_cast< EstimationInput< ObservationScalarType, TimeType > >( estimationInput )
                        ->getSaveStateHistoryForEachIteration( ) )
            {
                simulationResults = variationalEquationsSolver_->getVariationalPropagationResults( );
            }
        }
    }
    catch( std::runtime_error& error )
    {
        std::cerr << "Error when resetting parameters during parameter estimation: " << std::endl
                  << error.what( ) << std::endl
                  << "Terminating estimation" << std::endl;
        exceptionDuringPropagation = true;
    }

    if( estimationInput->getPrintOutput( ) )
    {
        std::cout << "Calculating residuals and partials " << totalNumberOfObservations << std::endl;
    }

    // Calculate residuals and observation matrix for current parameter estimate.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residuals;
    Eigen::MatrixXd designMatrix;
    if( calculateResiduals )
    {
        calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >( estimationInput->getObservationDataset( ),
                                                                              flattenedObservationData,
                                                                              observationManagers_,
                                                                              totalNumberParameters_,
                                                                              designMatrix,
                                                                              residuals,
                                                                              true );
    }
    else
    {
        calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >( estimationInput->getObservationDataset( ),
                                                                              flattenedObservationData,
                                                                              observationManagers_,
                                                                              totalNumberParameters_,
                                                                              designMatrix,
                                                                              residuals,
                                                                              false,
                                                                              true );
    }

    // Divide partials matrix between estimated and consider parameters
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > designMatrices =
            separateEstimatedAndConsiderDesignMatrices( designMatrix, totalNumberOfObservations );

    return std::make_pair( designMatrices, residuals );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ORBITDETERMINATIONMANAGERESTIMATIONIMPLEMENTATION_H
