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

#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/simulation/estimation_setup/interArcContinuityConstraint.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManagerHelpers.h"

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

    // Get number of observations
    int totalNumberOfObservations = estimationInput->getObservationCollection( )->getTotalObservableSize( );

    if( numberEstimatedParameters_ > static_cast< unsigned int >( totalNumberOfObservations ) &&
        estimationInput->getInverseOfAprioriCovariance( ).rows( ) == 0 )
    {
        std::cerr << "Warning when estimating parameters, number of observations is smaller than number of estimated parameters, and "
                     "no a priori information is provided."
                  << std::endl;
    }

    const Eigen::VectorXd weightsMatrixDiagonals = estimationInput->getWeightsMatrixDiagonals( );
    if( weightsMatrixDiagonals.rows( ) != totalNumberOfObservations )
    {
        throw std::runtime_error( "Error when estimating parameters, size of weights diagonal (" +
                                  std::to_string( weightsMatrixDiagonals.rows( ) ) + ") is not compatible with number of observations (" +
                                  std::to_string( totalNumberOfObservations ) + ")" );
    }
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

    // Inter-arc continuity-prior setup. Empty constraint list (the default) skips the feature entirely.
    const auto& interArcConstraints = estimationInput->getInterArcContinuityConstraints( );
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

    bool exceptionDuringPropagation = false, exceptionDuringInversion = false;

    // Iterate until convergence (at least once)
    int bestIteration = -1;
    int numberOfIterations = 0;
    while( true )
    {
        oldParameterEstimate = newParameterEstimate;

        // Compute design matrices (for estimated and consider parameters) and residuals.
        std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > simulationResults;
        std::pair< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
                designMatricesAndResiduals = performPreEstimationSteps(
                        estimationInput, newParameterEstimate, true, numberOfIterations, exceptionDuringPropagation, simulationResults );
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residuals = designMatricesAndResiduals.second;
        Eigen::MatrixXd designMatrixEstimatedParameters = designMatricesAndResiduals.first.first;
        Eigen::MatrixXd designMatrixConsiderParameters;
        designMatrixConsiderParameters = designMatricesAndResiduals.first.second;

        // Set simulation results
        if( estimationInput->getSaveStateHistoryForEachIteration( ) )
        {
            simulationResultsPerIteration.push_back( simulationResults->clone( ) );
        }

        // Normalise estimated parameters partials and inverse apriori covariance
        Eigen::VectorXd normalizationTerms = normalizeDesignMatrix( designMatrixEstimatedParameters );
        Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix = normalizeAprioriCovariance(
                estimationInput->getInverseOfAprioriCovariance( numberEstimatedParameters_ ), normalizationTerms );

        // Assemble soft inter-arc continuity-prior contribution for this iteration. The normalisation factors are the same
        // ones just applied to the observation design matrix above.
        InterArcConstraintContribution interArcContribution =
                assembleInterArcContinuityContributionFromManagerInterfaces< ObservationScalarType, TimeType >(
                        interArcConstraints,
                        parametersToEstimate_,
                        stateTransitionAndSensitivityMatrixInterface_,
                        variationalEquationsSolver_,
                        normalizationTerms,
                        static_cast< int >( numberEstimatedParameters_ ),
                        "parameter estimation" );
        interArcContinuityCostHistory.push_back( interArcContribution.totalConstraintCost );
        interArcContinuityDiscrepancyHistory.push_back( interArcContribution.perPairDiscrepancies );

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
                                                                                   interArcContribution.additionalRightHandSide ) );

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

        // Compute contribution consider parameters
        Eigen::MatrixXd covarianceContributionConsiderParameters;
        if( considerParametersIncluded_ )
        {
            covarianceContributionConsiderParameters =
                    linear_algebra::calculateConsiderParametersCovarianceContribution( ( leastSquaresOutput.second ).inverse( ),
                                                                                       designMatrixEstimatedParameters,
                                                                                       weightsMatrixDiagonals,
                                                                                       designMatrixConsiderParameters,
                                                                                       normalizedConsiderCovariance );
        }
        else
        {
            covarianceContributionConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
        }

        // Calculate mean residual for current iteration.
        residualRms = linear_algebra::getVectorEntryRootMeanSquare( residuals.template cast< double >( ) );
        costFunction = linear_algebra::computeLeastSquaresCostFunction( weightsMatrixDiagonals, residuals.template cast< double >( ) );
        // The cost driving best-iteration selection combines the observation cost with the inter-arc continuity-prior
        // cost (zero when no continuity priors are attached). Residual RMS is unchanged so observation-only diagnostics
        // remain meaningful.
        costFunction += interArcContribution.totalConstraintCost;
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
            std::map< observation_models::ObservableType, std::pair< int, int > > indicesPerObservableType =
                    estimationInput->getObservationCollection( )->getObservableTypeStartAndEndIndices( );

            if( indicesPerObservableType.size( ) > 1 )
            {
                for( auto it : indicesPerObservableType )
                {
                    double currentTypeRms = linear_algebra::getVectorEntryRootMeanSquare(
                            residuals.segment( it.second.first, it.second.second ).template cast< double >( ) );
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
        if( costFunction < bestCostFunction || !( bestCostFunction == bestCostFunction ) )
        {
            bestCostFunction = costFunction;
            bestRmsResidual = residualRms;
            bestParameterEstimate = oldParameterEstimate;
            bestResiduals = std::move( residuals.template cast< double >( ) );
            estimationInput->getObservationCollection( )->setResiduals( residuals );
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
                                                                                     exceptionDuringPropagation );

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
        const ParameterVectorType& newParameterEstimate,
        const bool calculateResiduals,
        const int numberOfIterations,
        bool& exceptionDuringPropagation,
        std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > >& simulationResults )
{
    // Get number of observations
    int totalNumberOfObservations = estimationInput->getObservationCollection( )->getTotalObservableSize( );

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
        calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >( estimationInput->getObservationCollection( ),
                                                                              observationManagers_,
                                                                              totalNumberParameters_,
                                                                              totalNumberOfObservations,
                                                                              designMatrix,
                                                                              residuals,
                                                                              true );
    }
    else
    {
        calculateDesignMatrix< ObservationScalarType, TimeType >( estimationInput->getObservationCollection( ),
                                                                  observationManagers_,
                                                                  totalNumberParameters_,
                                                                  totalNumberOfObservations,
                                                                  designMatrix );
    }

    // Divide partials matrix between estimated and consider parameters
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > designMatrices =
            separateEstimatedAndConsiderDesignMatrices( designMatrix, totalNumberOfObservations );

    return std::make_pair( designMatrices, residuals );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ORBITDETERMINATIONMANAGERESTIMATIONIMPLEMENTATION_H
