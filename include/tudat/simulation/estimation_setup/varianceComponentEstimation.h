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
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
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

//! Settings for an estimation in which the weights of observation groups are estimated with Helmert's method.
class VarianceComponentSettings
{
public:
    /*!
     * Constructor
     * \param groups One parser per observation group. Each group gets its own variance factor. Observations in no group
     * keep their weights.
     * \param tolerance The iteration stops when all variance factors differ from 1 by less than this value
     * \param maximumNumberOfIterations Maximum number of estimations (each one a full least-squares estimation)
     * \param minimumRedundancy Groups with a smaller redundancy are kept fixed
     * \param warmStart Whether each estimation starts from the result of the previous one. If false, or if the a
     * priori constraint is applied to the total parameter deviation, every estimation starts from the initial
     * parameter values, so that the a priori values stay the same.
     * \param printOutput Whether to print the variance factors of each iteration
     */
    VarianceComponentSettings( const std::vector< std::shared_ptr< observation_models::ObservationCollectionParser > >& groups,
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

    std::vector< std::shared_ptr< observation_models::ObservationCollectionParser > > groups_;
    double tolerance_;
    int maximumNumberOfIterations_;
    double minimumRedundancy_;
    bool warmStart_;
    bool printOutput_;
};

inline std::shared_ptr< VarianceComponentSettings > varianceComponentSettings(
        const std::vector< std::shared_ptr< observation_models::ObservationCollectionParser > >& groups,
        const double tolerance = 0.01,
        const int maximumNumberOfIterations = 10,
        const double minimumRedundancy = 1.0,
        const bool warmStart = true,
        const bool printOutput = true )
{
    return std::make_shared< VarianceComponentSettings >(
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

//! Function to get the concatenated row indices (as used in the design matrix) of the observations selected by a parser.
template< typename ObservationScalarType, typename TimeType >
std::vector< int > getObservationRowIndices(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< observation_models::ObservationCollectionParser > observationParser )
{
    std::vector< int > rowIndices;
    auto startAndSizes = observationCollection->getObservationSetStartAndSize( observationParser );
    for( auto observableIt : startAndSizes )
    {
        for( auto linkEndsIt : observableIt.second )
        {
            for( auto startAndSize : linkEndsIt.second )
            {
                for( int i = 0; i < startAndSize.second; i++ )
                {
                    rowIndices.push_back( startAndSize.first + i );
                }
            }
        }
    }
    return rowIndices;
}

//! Function to divide the weights of the observations selected by a parser by a constant factor.
/*!
 * Function to divide the weights of the observations selected by a parser by a constant factor. Relative weights
 * within the selected observations are kept.
 */
template< typename ObservationScalarType, typename TimeType >
void scaleObservationWeights(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::shared_ptr< observation_models::ObservationCollectionParser > observationParser,
        const double varianceFactor )
{
    for( auto observationSet : observationCollection->getSingleObservationSets( observationParser ) )
    {
        observationSet->setTabulatedWeights( observationSet->getWeightsVector( ) / varianceFactor );
    }
}

//! Function to compute the Helmert variance components of each observation group, from an estimation output.
template< typename ObservationScalarType, typename TimeType >
VarianceComponentEstimate computeHelmertVarianceComponents(
        const std::shared_ptr< EstimationOutput< ObservationScalarType, TimeType > > estimationOutput,
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationCollection,
        const std::vector< std::shared_ptr< observation_models::ObservationCollectionParser > >& groups,
        const double minimumRedundancy = 1.0 )
{
    if( !estimationOutput->designMatrixSaved_ )
    {
        throw std::runtime_error(
                "Error when computing variance components: the design matrix was not saved in the estimation output. "
                "Enable saving of the design matrix in the estimation settings." );
    }

    std::vector< std::vector< int > > groupRowIndices;
    for( auto group : groups )
    {
        groupRowIndices.push_back( getObservationRowIndices( observationCollection, group ) );
    }

    return computeHelmertVarianceComponents( estimationOutput->getNormalizedDesignMatrix( ),
                                             estimationOutput->weightsMatrixDiagonal_,
                                             estimationOutput->getNormalizedInverseCovarianceMatrix( ),
                                             estimationOutput->residuals_,
                                             groupRowIndices,
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
 * The weights in the observation collection of the estimation input are changed by this function. Their final values
 * are the estimated weights. Linear equality constraints on the parameters are not supported.
 *
 * \param orbitDeterminationManager Object used to perform the estimation
 * \param estimationInput Input to the estimation. The weights in its observation collection are the starting weights.
 * \param settings Settings for the variance component estimation
 * \return Output of the last estimation, and the variance factors of each iteration
 */
template< typename ObservationScalarType, typename TimeType >
std::shared_ptr< VarianceComponentEstimationOutput< ObservationScalarType, TimeType > > estimateParametersWithVarianceComponents(
        OrbitDeterminationManager< ObservationScalarType, TimeType >& orbitDeterminationManager,
        const std::shared_ptr< EstimationInput< ObservationScalarType, TimeType > > estimationInput,
        const std::shared_ptr< VarianceComponentSettings > settings )
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

    auto observationCollection = estimationInput->getObservationCollection( );

    // Check that the groups are not empty and do not overlap.
    std::vector< bool > isRowInGroup( observationCollection->getTotalObservableSize( ), false );
    for( unsigned int k = 0; k < settings->groups_.size( ); k++ )
    {
        std::vector< int > rowIndices = getObservationRowIndices( observationCollection, settings->groups_.at( k ) );
        if( rowIndices.empty( ) )
        {
            throw std::runtime_error( "Error in variance component estimation: group " + std::to_string( k ) +
                                      " contains no observations." );
        }
        for( int row : rowIndices )
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
                output->estimationOutput, observationCollection, settings->groups_, settings->minimumRedundancy_ );

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
            scaleObservationWeights( observationCollection, settings->groups_.at( k ), varianceFactor );
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
