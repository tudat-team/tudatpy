/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Sahin, M., Cross, P. A., and Sellers, P. C. (1992). Variance component estimation applied to satellite
 *          laser ranging. Bulletin Geodesique, 66(3):284-295. doi:10.1007/BF02033189
 *      Foerstner, W. (1979). Ein Verfahren zur Schaetzung von Varianz- und Kovarianzkomponenten.
 *          Allgemeine Vermessungs-Nachrichten, 86:446-453.
 */

#ifndef TUDAT_VARIANCECOMPONENTESTIMATION_H
#define TUDAT_VARIANCECOMPONENTESTIMATION_H

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/simulation/estimation_setup/observationCondition.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{

//! Variance component results for one group of observations.
struct VarianceComponentGroupResult {
    //! Number of observations in the group that have a non-zero weight.
    int numberOfObservations = 0;

    //! Redundancy of the group: number of observations minus the share absorbed by the estimated parameters.
    double redundancy = 0.0;

    //! Weighted sum of squared residuals of the group.
    double weightedSquareSum = 0.0;

    //! Estimated variance factor. The group weights should be divided by this value.
    double varianceFactor = 1.0;

    //! Variance factor from the simplified estimator (weighted square sum divided by redundancy), for comparison.
    double simplifiedVarianceFactor = 1.0;

    //! Root mean square of the (unweighted) residuals of the group.
    double rmsResidual = 0.0;

    //! Whether the group was kept fixed because its redundancy is too small.
    bool isFrozen = false;
};

//! Variance component results for all groups, from a single least-squares solution.
struct VarianceComponentEstimate {
    //! Results per group, in the order in which the groups were given.
    std::vector< VarianceComponentGroupResult > groupResults;

    //! Helmert matrix of the groups that were not frozen.
    Eigen::MatrixXd helmertMatrix;

    //! Whether the simplified estimator was used because Helmert's equations gave a non-positive variance factor.
    bool simplifiedEstimatorUsed = false;

    //! Variance factors of all groups, as a vector.
    Eigen::VectorXd getVarianceFactors( ) const
    {
        Eigen::VectorXd factors( groupResults.size( ) );
        for( unsigned int i = 0; i < groupResults.size( ); i++ )
        {
            factors( i ) = groupResults.at( i ).varianceFactor;
        }
        return factors;
    }

    //! Redundancies of all groups, as a vector.
    Eigen::VectorXd getRedundancies( ) const
    {
        Eigen::VectorXd redundancies( groupResults.size( ) );
        for( unsigned int i = 0; i < groupResults.size( ); i++ )
        {
            redundancies( i ) = groupResults.at( i ).redundancy;
        }
        return redundancies;
    }
};

//! Function to estimate variance components with Helmert's method, from a single least-squares solution.
/*!
 * Implements Helmert's equations H s = c of Sahin et al. (1992, appendix, Eqs. 15-19). For group k with normal
 * matrix N_k = A_k^T W_k A_k, total normal matrix N, n_k observations and residuals v_k:
 *
 *      h_kk = n_k - 2 tr(N^-1 N_k) + tr(N^-1 N_k N^-1 N_k)
 *      h_kl = tr(N^-1 N_k N^-1 N_l)
 *      c_k  = v_k^T W_k v_k
 *
 * The paper assumes that N is the sum of the group normal matrices. Here, N may also contain a priori information,
 * other soft constraints and observations that are in no group. These are treated as one extra component with a
 * known variance factor of 1, and moved to the right-hand side. This extension keeps the same solution as the
 * simplified estimator of Foerstner (1979), s_k = c_k / (n_k - tr(N^-1 N_k)), once the iteration has converged.
 *
 * If Helmert's equations give a non-positive or non-finite factor, the simplified estimator is used for all groups
 * instead, and this is flagged in the output. A group whose redundancy is below the given minimum is kept fixed at a
 * variance factor of 1.
 *
 * The design matrix and the inverse covariance may be normalized (as in Tudat's estimation output), as long as both
 * use the same normalization. All results are invariant under this normalization.
 *
 * \param designMatrix Design matrix of the solution (one row per observation)
 * \param weights Diagonal of the observation weight matrix
 * \param inverseCovariance Full normal matrix of the solution, including any a priori information
 * \param residuals Post-fit residuals
 * \param groupRowIndices Row indices per group. Groups must not overlap.
 * \param minimumRedundancy Groups with a smaller redundancy are kept fixed
 * \return Variance component results
 */
VarianceComponentEstimate computeHelmertVarianceComponents( const Eigen::MatrixXd& designMatrix,
                                                            const Eigen::VectorXd& weights,
                                                            const Eigen::MatrixXd& inverseCovariance,
                                                            const Eigen::VectorXd& residuals,
                                                            const std::vector< std::vector< int > >& groupRowIndices,
                                                            const double minimumRedundancy = 1.0 );

//! Function to estimate variance components with Helmert's method, for weights with correlations.
/*!
 * Same as the function above, but with a full (sparse) weight matrix. Correlations are only allowed between observations
 * of the same group: then N_k = A_k^T W_kk A_k and c_k = v_k^T W_kk v_k. The diagonal of the weight matrix must equal
 * the given weights. An empty weight matrix means diagonal weights.
 */
VarianceComponentEstimate computeHelmertVarianceComponents( const Eigen::MatrixXd& designMatrix,
                                                            const Eigen::VectorXd& weights,
                                                            const Eigen::SparseMatrix< double >& weightsMatrix,
                                                            const Eigen::MatrixXd& inverseCovariance,
                                                            const Eigen::VectorXd& residuals,
                                                            const std::vector< std::vector< int > >& groupRowIndices,
                                                            const double minimumRedundancy = 1.0 );

//! Settings for an estimation in which the weights of observation groups are estimated with Helmert's method.
template< typename ObservationScalarType = double, typename TimeType = double >
class VarianceComponentSettings
{
public:
    using ConditionType = observation_models::ObservationSelectionCondition< ObservationScalarType, TimeType >;

    /*!
     * Constructor
     * \param groups One selection condition per observation group. Each group gets its own variance factor. Observations
     * in no group keep their weights.
     * \param tolerance The iteration stops when all variance factors differ from 1 by less than this value
     * \param maximumNumberOfIterations Maximum number of estimations (each one a full least-squares estimation)
     * \param minimumRedundancy Groups with a smaller redundancy are kept fixed
     * \param warmStart Whether each estimation starts from the result of the previous one. If false, or if the a
     * priori constraint is applied to the total parameter deviation, every estimation starts from the initial
     * parameter values, so that the a priori values stay the same.
     * \param printOutput Whether to print the variance factors of each iteration
     */
    VarianceComponentSettings( const std::vector< ConditionType >& groups,
                               const double tolerance = 0.01,
                               const int maximumNumberOfIterations = 10,
                               const double minimumRedundancy = 1.0,
                               const bool warmStart = true,
                               const bool printOutput = true ):
        groups_( groups ), tolerance_( tolerance ), maximumNumberOfIterations_( maximumNumberOfIterations ),
        minimumRedundancy_( minimumRedundancy ), warmStart_( warmStart ), printOutput_( printOutput )
    {
        if( groups_.empty( ) )
        {
            throw std::runtime_error( "Error in variance component settings: no observation groups provided." );
        }
        if( maximumNumberOfIterations_ < 1 )
        {
            throw std::runtime_error( "Error in variance component settings: maximum number of iterations must be at least 1." );
        }
    }

    std::vector< ConditionType > groups_;
    double tolerance_;
    int maximumNumberOfIterations_;
    double minimumRedundancy_;
    bool warmStart_;
    bool printOutput_;
};

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< VarianceComponentSettings< ObservationScalarType, TimeType > > varianceComponentSettings(
        const std::vector< observation_models::ObservationSelectionCondition< ObservationScalarType, TimeType > >& groups,
        const double tolerance = 0.01,
        const int maximumNumberOfIterations = 10,
        const double minimumRedundancy = 1.0,
        const bool warmStart = true,
        const bool printOutput = true )
{
    return std::make_shared< VarianceComponentSettings< ObservationScalarType, TimeType > >(
            groups, tolerance, maximumNumberOfIterations, minimumRedundancy, warmStart, printOutput );
}

//! Output of an estimation with variance component estimation.
template< typename ObservationScalarType = double, typename TimeType = double >
struct VarianceComponentEstimationOutput {
    //! Output of the last least-squares estimation.
    std::shared_ptr< EstimationOutput< ObservationScalarType, TimeType > > estimationOutput;

    //! Variance component results of each iteration.
    std::vector< VarianceComponentEstimate > iterationResults;

    //! Total weight scaling per group used in each iteration's estimation (product of all earlier variance factors).
    std::vector< Eigen::VectorXd > totalVarianceFactorHistory;

    //! Best iteration of the least-squares estimation, for each iteration.
    std::vector< int > leastSquaresBestIterations;

    //! Whether all variance factors converged to 1 within the tolerance.
    bool isConverged = false;

    //! Total weight scaling per group. Each group sigma has been multiplied by the square root of this value.
    Eigen::VectorXd getTotalVarianceFactors( ) const
    {
        return totalVarianceFactorHistory.empty( ) ? Eigen::VectorXd( ) : totalVarianceFactorHistory.back( );
    }
};

//! Function to get the rows (in the order used for estimation) of the active observations matching each condition.
/*!
 * Only active observations have rows in the estimation. Observations with more than one component (e.g. position) have one row
 * per component.
 */
template< typename ObservationScalarType, typename TimeType >
std::vector< std::vector< int > > getEstimationRowIndices(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const std::vector< observation_models::ObservationSelectionCondition< ObservationScalarType, TimeType > >& groups )
{
    const std::vector< unsigned int > rowObservationIds =
            observationDataset->createOrderedObservationVectorData( false ).getObservationIds( );

    std::vector< std::vector< int > > groupRowIndices;
    for( const auto& group : groups )
    {
        const std::vector< unsigned int > selectedIds = observationDataset->getObservationIdsMatchingCondition( group );
        const std::unordered_set< unsigned int > selected( selectedIds.begin( ), selectedIds.end( ) );
        std::vector< int > rows;
        for( unsigned int row = 0; row < rowObservationIds.size( ); row++ )
        {
            if( selected.count( rowObservationIds.at( row ) ) )
            {
                rows.push_back( row );
            }
        }
        groupRowIndices.push_back( rows );
    }
    return groupRowIndices;
}

//! Function to compute the Helmert variance components of each observation group, from an estimation output.
/*!
 * The observation dataset must be in the same state as during the estimation (same active observations and weights), so that
 * its rows match those of the estimation output.
 */
template< typename ObservationScalarType, typename TimeType >
VarianceComponentEstimate computeHelmertVarianceComponents(
        const std::shared_ptr< EstimationOutput< ObservationScalarType, TimeType > > estimationOutput,
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const std::vector< observation_models::ObservationSelectionCondition< ObservationScalarType, TimeType > >& groups,
        const double minimumRedundancy = 1.0 )
{
    if( !estimationOutput->designMatrixSaved_ )
    {
        throw std::runtime_error(
                "Error when computing variance components: the design matrix was not saved in the estimation output. "
                "Enable saving of the design matrix in the estimation settings." );
    }

    // Check that the dataset still matches the estimation output, row by row.
    const Eigen::VectorXd currentWeights = observationDataset->createOrderedObservationVectorData( false ).getWeightVector( );
    if( currentWeights.rows( ) != estimationOutput->residuals_.rows( ) ||
        !( currentWeights - estimationOutput->weightsMatrixDiagonal_ ).isZero( 0.0 ) )
    {
        throw std::runtime_error(
                "Error when computing variance components: the observation dataset no longer matches the estimation output "
                "(its active observations or weights have changed)." );
    }

    return computeHelmertVarianceComponents(
            estimationOutput->getNormalizedDesignMatrix( ),
            estimationOutput->weightsMatrixDiagonal_,
            estimationOutput->hasFullWeightMatrix( ) ? estimationOutput->getWeightsMatrix( ) : Eigen::SparseMatrix< double >( ),
            estimationOutput->getNormalizedInverseCovarianceMatrix( ),
            estimationOutput->residuals_,
            getEstimationRowIndices( observationDataset, groups ),
            minimumRedundancy );
}

//! Function to estimate parameters, while estimating the weight of each observation group with Helmert's method.
/*!
 * Function to estimate parameters, while estimating the weight of each observation group with Helmert's method
 * (Sahin et al., 1992). The procedure follows Fig. 1 of that paper:
 *
 *  1. Estimate the parameters with the current weights (a full, converged least-squares estimation).
 *  2. Estimate the variance factor s_k of each group with Helmert's equations.
 *  3. Divide the weights of each group by s_k, and go back to 1 until all s_k are close to 1.
 *
 * The weights in the observation dataset of the estimation input are changed by this function. Their final values
 * are the estimated weights. Linear equality constraints on the parameters are not supported.
 *
 * \param orbitDeterminationManager Object used to perform the estimation
 * \param estimationInput Input to the estimation. The weights in its observation dataset are the starting weights. Outlier rejection is not
 * supported.
 * \param settings Settings for the variance component estimation
 * \return Output of the last estimation, and the variance factors of each iteration
 */
template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< VarianceComponentEstimationOutput< ObservationScalarType, TimeType > > estimateParametersWithVarianceComponents(
        OrbitDeterminationManager< ObservationScalarType, TimeType >& orbitDeterminationManager,
        const std::shared_ptr< EstimationInput< ObservationScalarType, TimeType > > estimationInput,
        const std::shared_ptr< VarianceComponentSettings< ObservationScalarType, TimeType > > settings )
{
    using ParameterVectorType = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >;

    if( !estimationInput->getSaveDesignMatrix( ) )
    {
        throw std::runtime_error(
                "Error in variance component estimation: the design matrix must be saved. Enable saving of the design "
                "matrix in the estimation settings." );
    }

    auto parametersToEstimate = orbitDeterminationManager.getParametersToEstimate( );
    Eigen::MatrixXd constraintMultiplier;
    Eigen::VectorXd constraintRightHandSide;
    parametersToEstimate->getConstraints( constraintMultiplier, constraintRightHandSide );
    if( constraintMultiplier.rows( ) > 0 )
    {
        throw std::runtime_error( "Error in variance component estimation: linear parameter constraints are not supported." );
    }

    if( estimationInput->getOutlierRejectionSettings( ) != nullptr )
    {
        throw std::runtime_error(
                "Error in variance component estimation: outlier rejection during the estimation is not supported, because it "
                "changes the observations of each group. Reject outliers before the variance component estimation." );
    }

    auto observationDataset = estimationInput->getObservationDataset( );

    // Check that the groups are not empty and do not overlap. Group membership is fixed during the iteration.
    const std::vector< std::vector< int > > groupRowIndices = getEstimationRowIndices( observationDataset, settings->groups_ );
    std::vector< bool > isRowInGroup( observationDataset->createOrderedObservationVectorData( false ).getObservationIds( ).size( ), false );
    for( unsigned int k = 0; k < groupRowIndices.size( ); k++ )
    {
        if( groupRowIndices.at( k ).empty( ) )
        {
            throw std::runtime_error( "Error in variance component estimation: group " + std::to_string( k ) +
                                      " contains no active observations." );
        }
        for( int row : groupRowIndices.at( k ) )
        {
            if( isRowInGroup.at( row ) )
            {
                throw std::runtime_error( "Error in variance component estimation: group " + std::to_string( k ) +
                                          " overlaps with an earlier group." );
            }
            isRowInGroup.at( row ) = true;
        }
    }

    // When the a priori constraint applies to the total deviation, it is centred on the parameter values at the start
    // of each estimation. Starting from the initial values every time keeps that centre fixed.
    const bool warmStart = settings->warmStart_ && !estimationInput->getApplyAprioriParameterDeviation( );
    const ParameterVectorType initialParameterValues = parametersToEstimate->template getFullParameterValues< ObservationScalarType >( );

    auto output = std::make_shared< VarianceComponentEstimationOutput< ObservationScalarType, TimeType > >( );
    Eigen::VectorXd totalVarianceFactors = Eigen::VectorXd::Ones( settings->groups_.size( ) );

    for( int iteration = 0; iteration < settings->maximumNumberOfIterations_; iteration++ )
    {
        // Set the parameter values at which this estimation starts. After an estimation, the parameter set holds the
        // last iterate, which is not necessarily the best one.
        ParameterVectorType startValues = initialParameterValues;
        if( warmStart && output->estimationOutput != nullptr )
        {
            startValues = output->estimationOutput->parameterEstimate_;
        }
        if( estimationInput->getReintegrateEquationsOnFirstIteration( ) )
        {
            parametersToEstimate->template resetParameterValues< ObservationScalarType >( startValues );
        }
        else
        {
            orbitDeterminationManager.resetParameterEstimate( startValues, estimationInput->getReintegrateVariationalEquations( ) );
        }

        output->totalVarianceFactorHistory.push_back( totalVarianceFactors );
        output->estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
        output->leastSquaresBestIterations.push_back( output->estimationOutput->bestIteration_ );

        VarianceComponentEstimate estimate = computeHelmertVarianceComponents(
                output->estimationOutput, observationDataset, settings->groups_, settings->minimumRedundancy_ );

        bool isConverged = true;
        for( unsigned int k = 0; k < estimate.groupResults.size( ); k++ )
        {
            const VarianceComponentGroupResult& groupResult = estimate.groupResults.at( k );
            if( !groupResult.isFrozen && std::fabs( groupResult.varianceFactor - 1.0 ) >= settings->tolerance_ )
            {
                isConverged = false;
            }
        }

        if( settings->printOutput_ )
        {
            std::cout << "Variance component iteration " << iteration
                      << ( estimate.simplifiedEstimatorUsed ? " (simplified estimator)" : "" ) << std::endl;
            for( unsigned int k = 0; k < estimate.groupResults.size( ); k++ )
            {
                const VarianceComponentGroupResult& groupResult = estimate.groupResults.at( k );
                std::cout << "  Group " << k << ": n = " << groupResult.numberOfObservations << ", redundancy = " << groupResult.redundancy
                          << ", variance factor = " << groupResult.varianceFactor << ", rms = " << groupResult.rmsResidual
                          << ( groupResult.isFrozen ? " (fixed, redundancy too small)" : "" ) << std::endl;
            }
        }

        output->iterationResults.push_back( estimate );
        output->isConverged = isConverged;

        // Stop without changing the weights, so that the final weights match the returned estimation output.
        if( isConverged || iteration == settings->maximumNumberOfIterations_ - 1 )
        {
            break;
        }

        // Rescale the weights of each group.
        for( unsigned int k = 0; k < estimate.groupResults.size( ); k++ )
        {
            const double varianceFactor = estimate.groupResults.at( k ).varianceFactor;
            if( !( varianceFactor > 0.0 ) || !std::isfinite( varianceFactor ) )
            {
                throw std::runtime_error( "Error in variance component estimation: invalid variance factor " +
                                          std::to_string( varianceFactor ) + " for group " + std::to_string( k ) + "." );
            }
            observationDataset->scaleWeights( settings->groups_.at( k ), 1.0 / varianceFactor );
            totalVarianceFactors( k ) *= varianceFactor;
        }
    }

    if( !output->isConverged )
    {
        std::cerr << "Warning in variance component estimation: variance factors did not converge within "
                  << settings->maximumNumberOfIterations_ << " iterations." << std::endl;
    }

    return output;
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_VARIANCECOMPONENTESTIMATION_H
