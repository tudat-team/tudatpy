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

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>

#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/simulation/estimation_setup/interArcContinuityConstraint.h"
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
    if( considerParametersIncluded_ && parametersToEstimate_->getConstraintSize( ) > 0 )
    {
        throw std::runtime_error(
                "Error when estimating parameters, linear parameter constraints cannot be combined with consider parameters" );
    }

    currentParameterEstimate_ = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );

    const auto observationDataset = estimationInput->getObservationDataset( );
    const std::shared_ptr< OutlierRejection< ObservationScalarType, TimeType > > outlierRejection =
            createOutlierRejection< ObservationScalarType, TimeType >( estimationInput->getOutlierRejectionSettings( ),
                                                                       observationDataset );
    const bool applyOutlierRejection = ( outlierRejection != nullptr );
    const auto carpinoSettings =
            std::dynamic_pointer_cast< CarpinoOutlierRejectionSettings >( estimationInput->getOutlierRejectionSettings( ) );

    // Get number of observations
    observation_models::ObservationVectorData< ObservationScalarType, TimeType > estimationData =
            observationDataset->createObservationVectorData( false );
    int totalNumberOfObservations = static_cast< int >( estimationData.getObservationVector( ).size( ) );
    const int initialNumberOfObservations = totalNumberOfObservations;
    for( int i = 0; i < totalNumberOfObservations; i++ )
    {
        if( !std::isfinite( static_cast< double >( estimationData.getObservationVector( )( i ) ) ) )
        {
            throw std::runtime_error( "Error in estimation, observation " + std::to_string( estimationData.getObservationIds( ).at( i ) ) +
                                      " (at time " + std::to_string( static_cast< double >( estimationData.getTimes( ).at( i ) ) ) +
                                      ") is not finite; remove or reject it before the estimation." );
        }
    }

    if( totalNumberOfObservations == 0 )
    {
        throw std::runtime_error( "Cannot run estimation or covariance analysis without active observations." );
    }

    // Residuals are computed for all observations (including rejected ones) for outlier rejection and the residual history, and
    // to store the residuals of rejected observations
    const bool computeResidualsForAllObservations = applyOutlierRejection ||
            estimationInput->getSaveResidualsAndParametersFromEachIteration( ) ||
            totalNumberOfObservations != static_cast< int >( observationDataset->getTotalScalarSize( ) );
    observation_models::ObservationVectorData< ObservationScalarType, TimeType > computationData;
    if( computeResidualsForAllObservations )
    {
        computationData = observationDataset->createObservationVectorData( true );
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
    double bestSelectionCostFunction = TUDAT_NAN;
    std::vector< bool > bestRejectionStatus;
    double bestRmsResidual = TUDAT_NAN;
    ParameterVectorType bestParameterEstimate = ParameterVectorType::Constant( numberEstimatedParameters_, TUDAT_NAN );
    Eigen::VectorXd bestTransformationData = Eigen::VectorXd::Constant( numberEstimatedParameters_, TUDAT_NAN );
    Eigen::VectorXd bestResiduals = Eigen::VectorXd::Constant( totalNumberOfObservations, TUDAT_NAN );
    std::vector< TimeType > bestResidualTimes;
    Eigen::MatrixXd bestDesignMatrixEstimatedParameters = Eigen::MatrixXd::Zero( 0, 0 );
    Eigen::VectorXd bestWeightsMatrixDiagonal = Eigen::VectorXd::Constant( totalNumberOfObservations, TUDAT_NAN );
    Eigen::SparseMatrix< double > bestWeightsMatrix;
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
    std::vector< typename EstimationOutput< ObservationScalarType, TimeType >::ActiveFlagsVector > activeFlagsPerIteration;
    std::vector< std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > > simulationResultsPerIteration;

    // Inter-arc continuity-prior setup. Empty constraint list (the default) skips the feature entirely.
    const auto& interArcConstraints = estimationInput->getInterArcContinuityConstraints( );
    double bestInterArcContinuityCost = 0.0;
    std::vector< Eigen::VectorXd > bestInterArcContinuityDiscrepancies;
    std::vector< double > interArcContinuityCostHistory;
    std::vector< std::vector< Eigen::VectorXd > > interArcContinuityDiscrepancyHistory;

    // Declare residual bookkeeping variables
    std::vector< double > rmsResidualHistory;
    std::vector< double > costFunctionHistory;
    double residualRms;
    std::map< observation_models::ObservableType, double > residualRmsPerType;
    double costFunction;

    // Set current parameter estimate as both previous and current estimate
    ParameterVectorType newParameterEstimate = currentParameterEstimate_;
    ParameterVectorType oldParameterEstimate = currentParameterEstimate_;
    const ParameterVectorType aprioriParameterEstimate = currentParameterEstimate_;

    bool exceptionDuringPropagation = false, exceptionDuringInversion = false;

    // Iterate until convergence (at least once)
    int bestIteration = -1;
    int numberOfIterations = 0;

    // If an exception leaves the iteration loop, the rejection status of the observations is reset to that of the best iteration
    // (or to that at the start of the estimation), so that the observation dataset is not left with that of an unfinished iteration
    struct RejectionStatusGuard {
        const std::shared_ptr< OutlierRejection< ObservationScalarType, TimeType > >& outlierRejection;
        const int& bestIteration;
        const std::vector< bool >& bestRejectionStatus;
        const std::vector< bool > initialRejectionStatus;
        ~RejectionStatusGuard( )
        {
            if( outlierRejection != nullptr && std::uncaught_exceptions( ) > 0 )
            {
                outlierRejection->resetRejectionStatus( bestIteration >= 0 ? bestRejectionStatus : initialRejectionStatus );
            }
        }
    };
    const RejectionStatusGuard rejectionStatusGuard{ outlierRejection,
                                                     bestIteration,
                                                     bestRejectionStatus,
                                                     applyOutlierRejection ? outlierRejection->getRejectionStatus( )
                                                                           : std::vector< bool >( ) };

    while( true )
    {
        oldParameterEstimate = newParameterEstimate;
        const bool needsParameterCovariance = carpinoSettings && numberOfIterations >= carpinoSettings->getFirstIterationWithRejection( );

        // Update the data of the observations that are used in the estimation: observations may have been rejected or
        // recovered at the end of the previous iteration.
        if( applyOutlierRejection && numberOfIterations > 0 )
        {
            estimationData = observationDataset->createObservationVectorData( false );
            computationData = observationDataset->createObservationVectorData( true );
            totalNumberOfObservations = static_cast< int >( estimationData.getObservationVector( ).size( ) );
            if( totalNumberOfObservations == 0 )
            {
                outlierRejection->resetRejectionStatus( bestRejectionStatus );
                throw std::runtime_error(
                        "Error during parameter estimation, all observations have been rejected by the outlier "
                        "rejection algorithm." );
            }
            retrieveObservationWeights< ObservationScalarType, TimeType >(
                    estimationData, weightsMatrixDiagonals, weightsMatrix, hasOffDiagonalWeights );
        }

        // Compute design matrices (for estimated and consider parameters) and residuals.
        // When saving residual history or applying outlier rejection, all observations are used and active rows extracted below.
        std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > simulationResults;
        std::pair< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                designMatricesAndResiduals =
                        performPreEstimationSteps( estimationInput,
                                                   observationDataset,
                                                   newParameterEstimate,
                                                   computeResidualsForAllObservations ? computationData : estimationData,
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
        if( computeResidualsForAllObservations && computationData.getObservationIds( ) == estimationData.getObservationIds( ) )
        {
            residuals = computedResiduals;  // Keep the original for residual history and rejection.
            designMatrixEstimatedParameters = needsParameterCovariance ? computedDesignMatrixEstimatedParameters
                                                                       : std::move( computedDesignMatrixEstimatedParameters );
            designMatrixConsiderParameters = std::move( computedDesignMatrixConsiderParameters );
        }
        else if( computeResidualsForAllObservations )
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
            if( simulationResults == nullptr )
            {
                throw std::runtime_error(
                        "Error when saving state history for each estimation iteration: no propagation results exist "
                        "in observation-only estimation." );
            }
            if( estimationInput->getSaveVariationalResultsForEachIteration( ) )
            {
                simulationResultsPerIteration.push_back( simulationResults->clone( ) );
            }
            else
            {
                simulationResultsPerIteration.push_back( simulationResults->cloneDynamicsOnly( ) );
            }
        }

        // Normalise estimated parameters partials and inverse apriori covariance
        Eigen::VectorXd normalizationTerms = normalizeDesignMatrix( designMatrixEstimatedParameters );
        const Eigen::VectorXd estimatedParameterNormalizationTerms = normalizationTerms.segment( 0, numberEstimatedParameters_ );
        Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix = normalizeAprioriCovariance(
                estimationInput->getInverseOfAprioriCovariance( numberEstimatedParameters_ ), estimatedParameterNormalizationTerms );
        Eigen::VectorXd normalizedAprioriParameterDeviation = Eigen::VectorXd::Zero( 0 );
        if( estimationInput->getApplyAprioriParameterDeviation( ) )
        {
            normalizedAprioriParameterDeviation = ( oldParameterEstimate - aprioriParameterEstimate )
                                                          .template cast< double >( )
                                                          .cwiseProduct( estimatedParameterNormalizationTerms );
        }

        InterArcConstraintContribution interArcContribution;
        if( !interArcConstraints.empty( ) )
        {
            // Assemble soft inter-arc continuity-prior contribution for this iteration. The normalisation factors are the same
            // ones just applied to the observation design matrix above.
            interArcContribution = assembleInterArcContinuityContributionFromManagerInterfaces< ObservationScalarType, TimeType >(
                    interArcConstraints,
                    parametersToEstimate_,
                    stateTransitionAndSensitivityMatrixInterface_,
                    variationalEquationsSolver_,
                    normalizationTerms,
                    static_cast< int >( numberEstimatedParameters_ ),
                    "parameter estimation",
                    static_cast< int >( designMatrixEstimatedParameters.rows( ) ) );
            interArcContinuityCostHistory.push_back( interArcContribution.totalConstraintCost );
            interArcContinuityDiscrepancyHistory.push_back( interArcContribution.perPairDiscrepancies );
        }

        // Normalise partials w.r.t. consider parameters, consider covariance and parameters deviations
        Eigen::VectorXd normalizationTermsConsider, normalizedConsiderParametersDeviation;
        Eigen::MatrixXd normalizedConsiderCovariance;
        if( considerParametersIncluded_ )
        {
            normalizationTermsConsider = normalizeDesignMatrix( designMatrixConsiderParameters );
            getNormalizedConsiderCovariance( estimationInput, normalizationTermsConsider, normalizedConsiderCovariance );
            if( estimationInput->considerParametersDeviations_.rows( ) == 0 )
            {
                normalizedConsiderParametersDeviation = Eigen::VectorXd::Zero( normalizationTermsConsider.rows( ) );
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
            normalizeLinearConstraints( constraintStateMultiplier, constraintRightHandSide, normalizationTerms );

            double conditionNumberCheck = estimationInput->getLimitConditionNumberForWarning( );
            if( numberOfIterations > 0 && estimationInput->conditionNumberWarningEachIteration_ == false )
            {
                conditionNumberCheck = TUDAT_NAN;
            }
            // Perform LSQ inversion
            if( hasOffDiagonalWeights )
            {
                leastSquaresOutput = std::move(
                        linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrixEstimatedParameters,
                                                                                       residuals.template cast< double >( ),
                                                                                       weightsMatrix,
                                                                                       normalizedInverseAprioriCovarianceMatrix,
                                                                                       conditionNumberCheck,
                                                                                       constraintStateMultiplier,
                                                                                       constraintRightHandSide,
                                                                                       designMatrixConsiderParameters,
                                                                                       normalizedConsiderParametersDeviation,
                                                                                       interArcContribution.additionalNormalMatrix,
                                                                                       interArcContribution.additionalRightHandSide,
                                                                                       normalizedAprioriParameterDeviation ) );
            }
            else
            {
                leastSquaresOutput = std::move(
                        linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrixEstimatedParameters,
                                                                                       residuals.template cast< double >( ),
                                                                                       weightsMatrixDiagonals,
                                                                                       normalizedInverseAprioriCovarianceMatrix,
                                                                                       conditionNumberCheck,
                                                                                       constraintStateMultiplier,
                                                                                       constraintRightHandSide,
                                                                                       designMatrixConsiderParameters,
                                                                                       normalizedConsiderParametersDeviation,
                                                                                       interArcContribution.additionalNormalMatrix,
                                                                                       interArcContribution.additionalRightHandSide,
                                                                                       normalizedAprioriParameterDeviation ) );
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

        ParameterVectorType parameterAddition = ( leastSquaresOutput.first.cwiseQuotient( estimatedParameterNormalizationTerms ) )
                                                        .template cast< ObservationScalarType >( );

        // Compute the (unnormalized) covariance of the estimated parameters, which is only required by outlier
        // rejection algorithms. It is computed here, since the normalized inverse covariance and the normalization
        // terms are moved into the results of the best iteration further down.
        Eigen::MatrixXd parameterCovariance;
        if( needsParameterCovariance )
        {
            // Constraints add extra rows: invert first, then keep only the estimated parameters' covariance.
            const Eigen::MatrixXd normalizedCovariance = leastSquaresOutput.second.inverse( );
            parameterCovariance = normaliseUnnormaliseCovarianceMatrix(
                    normalizedCovariance.topLeftCorner( numberEstimatedParameters_, numberEstimatedParameters_ ),
                    normalizationTerms.segment( 0, numberEstimatedParameters_ ),
                    false );
        }

        if( estimationInput->iterationCallback_ )
        {
            estimationInput->iterationCallback_(
                    numberOfIterations,
                    residuals.template cast< double >( ),
                    normaliseUnnormaliseInverseCovarianceMatrix( leastSquaresOutput.second, normalizationTerms, false ) );
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
        double aprioriCost = 0.0;
        if( estimationInput->getApplyAprioriParameterDeviation( ) )
        {
            aprioriCost = 0.5 *
                    normalizedAprioriParameterDeviation.dot( normalizedInverseAprioriCovarianceMatrix *
                                                             normalizedAprioriParameterDeviation );
        }
        // In deviation-based mode, best-iteration selection additionally includes the absolute a priori cost. Residual RMS
        // is unchanged so observation-only diagnostics remain meaningful.
        // With outlier rejection, the number of observations that is used can differ per iteration. The best iteration is then
        // selected using the observation (and inter-arc constraint) cost per used observation, scaled to the number of
        // observations used in the first iteration.
        const double selectionCostFunction = applyOutlierRejection
                ? ( costFunction + interArcContribution.totalConstraintCost ) * static_cast< double >( initialNumberOfObservations ) /
                                static_cast< double >( totalNumberOfObservations ) +
                        aprioriCost
                : costFunction + ( aprioriCost + interArcContribution.totalConstraintCost );
        costFunction += aprioriCost + interArcContribution.totalConstraintCost;
        rmsResidualHistory.push_back( residualRms );
        costFunctionHistory.push_back( costFunction );

        if( estimationInput->getSaveResidualsAndParametersFromEachIteration( ) )
        {
            residualHistory.push_back( computedResiduals.template cast< double >( ) );
            typename EstimationOutput< ObservationScalarType, TimeType >::ActiveFlagsVector activeFlags( computedResiduals.rows( ) );
            for( int row = 0; row < activeFlags.rows( ); row++ )
            {
                activeFlags( row ) = observationDataset->getObservationRow( computationData.getObservationIds( ).at( row ) ).isActive_;
            }
            activeFlagsPerIteration.push_back( activeFlags );
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
                residualsPerObservableType[ observationDataset->getObservationSetMetadata( setId ).observableType_ ].push_back(
                        static_cast< double >( residuals( row ) ) );
            }

            if( residualsPerObservableType.size( ) > 1 )
            {
                for( auto it : residualsPerObservableType )
                {
                    double currentTypeRms =
                            linear_algebra::getVectorEntryRootMeanSquare( utilities::convertStlVectorToEigenVector< double >( it.second ) );
                    residualRmsPerType[ it.first ] = currentTypeRms;
                    std::cout << "Current residual for observable (" << observation_models::getObservableName( it.first )
                              << "): " << currentTypeRms << std::endl;
                }
            }
            else
            {
                std::cout << "Current residual: " << residualRms << std::endl;
            }
        }

        // If current iteration is better than previous one, update 'best' data.
        if( selectionCostFunction < bestSelectionCostFunction || !( bestSelectionCostFunction == bestSelectionCostFunction ) )
        {
            bestSelectionCostFunction = selectionCostFunction;
            bestCostFunction = costFunction;
            if( applyOutlierRejection )
            {
                bestRejectionStatus = outlierRejection->getRejectionStatus( );
            }
            bestRmsResidual = residualRms;
            bestParameterEstimate = oldParameterEstimate;
            bestResiduals = std::move( residuals.template cast< double >( ) );
            bestResidualTimes = estimationData.getTimes( );

            if( computeResidualsForAllObservations )
            {
                observationDataset->setResidualVector( computationData, computedResiduals );
            }
            else
            {
                observationDataset->setResidualVector( estimationData, residuals );
            }
            estimationInput->synchronizeLegacyResiduals( *observationDataset );
            if( estimationInput->getSaveDesignMatrix( ) )
            {
                bestDesignMatrixEstimatedParameters = std::move( designMatrixEstimatedParameters );
                bestDesignMatrixConsiderParameters = std::move( designMatrixConsiderParameters );
            }
            bestWeightsMatrixDiagonal = weightsMatrixDiagonals;
            bestWeightsMatrix = hasOffDiagonalWeights ? weightsMatrix : Eigen::SparseMatrix< double >( );
            bestTransformationData = std::move( normalizationTerms );
            bestInverseNormalizedCovarianceMatrix = std::move( leastSquaresOutput.second );
            bestIteration = numberOfIterations;
            bestConsiderTransformationData = std::move( normalizationTermsConsider );
            bestConsiderCovarianceContribution = covarianceContributionConsiderParameters;
            if( !interArcConstraints.empty( ) )
            {
                bestInterArcContinuityCost = interArcContribution.totalConstraintCost;
                bestInterArcContinuityDiscrepancies = interArcContribution.perPairDiscrepancies;
            }
        }

        // Update which observations are rejected, using the data of the current iteration. The observations that are
        // rejected here are excluded from the next iteration, and rejected observations that are recovered here are
        // included again.
        if( applyOutlierRejection )
        {
            const OutlierRejectionInput< ObservationScalarType, TimeType > outlierRejectionInput( numberOfIterations,
                                                                                                  computationData,
                                                                                                  computedResiduals,
                                                                                                  computedDesignMatrixEstimatedParameters,
                                                                                                  parameterCovariance,
                                                                                                  parameterAddition );
            outlierRejection->updateRejectionStatus( outlierRejectionInput );

            if( estimationInput->getPrintOutput( ) )
            {
                std::cout << "Number of rejected observations: " << outlierRejection->getNumberOfRejectedObservations( ) << " out of "
                          << observationDataset->getNumberOfObservations( ) << std::endl;
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

    // Set the rejection status of the observations to that of the best iteration, to which the estimation output and the
    // residuals in the observation dataset belong (the status was last updated for the iteration after the final one).
    if( applyOutlierRejection && bestIteration >= 0 )
    {
        outlierRejection->resetRejectionStatus( bestRejectionStatus );
    }

    if( estimationInput->getPrintOutput( ) )
    {
        std::cout << "Best iteration: " << bestIteration << " out of " << numberOfIterations - 1 << std::endl;
        std::cout << "Rms residual from best iteration: " << bestRmsResidual << std::endl;
        std::cout << "Cost function from best iteration: " << bestCostFunction << std::endl;
    }

    // Reset parameters to the best iteration (residuals in the observation collection are already from the best iteration). Note
    // that the environment (e.g. ephemerides of propagated bodies) is not reset, and corresponds to the last propagation.
    if( estimationInput->resetParametersToBestIteration_ )
    {
        if( bestIteration >= 0 )
        {
            parametersToEstimate_->template resetParameterValues< ObservationScalarType >( bestParameterEstimate );
            currentParameterEstimate_ = bestParameterEstimate;
        }
    }

    // Create estimation output object
    std::shared_ptr< EstimationOutput< ObservationScalarType, TimeType > > estimationOutput =
            std::make_shared< EstimationOutput< ObservationScalarType, TimeType > >( bestParameterEstimate,
                                                                                     bestResiduals,
                                                                                     bestDesignMatrixEstimatedParameters,
                                                                                     bestWeightsMatrixDiagonal,
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
                                                                                     exceptionDuringPropagation,
                                                                                     bestInterArcContinuityCost,
                                                                                     bestInterArcContinuityDiscrepancies,
                                                                                     activeFlagsPerIteration,
                                                                                     bestWeightsMatrix );
    estimationOutput->setResidualTimes( utilities::staticCastVector< double, TimeType >( bestResidualTimes ) );
    estimationOutput->setResidualHistoryObservationIds( computationData.getObservationIds( ) );

    if( estimationInput->getSaveStateHistoryForEachIteration( ) )
    {
        estimationOutput->setSimulationResults( simulationResultsPerIteration );
    }
    if( !interArcConstraints.empty( ) )
    {
        estimationOutput->setInterArcContinuityCostHistory( interArcContinuityCostHistory );
        estimationOutput->setInterArcContinuityDiscrepancyHistory( interArcContinuityDiscrepancyHistory );
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
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > >& observationDataset,
        const ParameterVectorType& newParameterEstimate,
        const observation_models::ObservationVectorData< ObservationScalarType, TimeType >& observationVectorData,
        const bool calculateResiduals,
        const int numberOfIterations,
        bool& exceptionDuringPropagation,
        std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > >& simulationResults )
{
    const int totalNumberOfObservations = static_cast< int >( observationVectorData.getObservationVector( ).size( ) );

    std::shared_ptr< EstimationInput< ObservationScalarType, TimeType > > parameterEstimationInput =
            std::dynamic_pointer_cast< EstimationInput< ObservationScalarType, TimeType > >( estimationInput );
    if( parameterEstimationInput != nullptr && parameterEstimationInput->getSaveStateHistoryForEachIteration( ) &&
        variationalEquationsSolver_ == nullptr )
    {
        throw std::runtime_error(
                "Error when preparing estimation iteration: saving state history requires propagated dynamics and "
                "is not supported in observation-only estimation." );
    }

    // Re-integrate equations of motion and variational equations with new parameter estimate.
    try
    {
        if( ( numberOfIterations > 0 ) || ( estimationInput->getReintegrateEquationsOnFirstIteration( ) ) )
        {
            resetParameterEstimate( newParameterEstimate, estimationInput->getReintegrateVariationalEquations( ) );
        }

        if( parameterEstimationInput != nullptr )
        {
            if( parameterEstimationInput->getSaveStateHistoryForEachIteration( ) )
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
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >( observationDataset,
                                                                          observationVectorData,
                                                                          observationManagers_,
                                                                          totalNumberParameters_,
                                                                          designMatrix,
                                                                          residuals,
                                                                          calculateResiduals );

    // Divide partials matrix between estimated and consider parameters
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > designMatrices =
            separateEstimatedAndConsiderDesignMatrices( designMatrix, totalNumberOfObservations );

    return std::make_pair( designMatrices, residuals );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ORBITDETERMINATIONMANAGERESTIMATIONIMPLEMENTATION_H
