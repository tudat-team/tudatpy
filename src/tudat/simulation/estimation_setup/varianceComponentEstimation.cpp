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
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <Eigen/Cholesky>
#include <Eigen/QR>

#include "tudat/simulation/estimation_setup/varianceComponentEstimation.h"

namespace tudat
{

namespace simulation_setup
{

namespace
{

//! Weighted normal matrix of the given rows, A_k^T W_k A_k, built in blocks of rows to limit memory use.
Eigen::MatrixXd computeGroupNormalMatrix( const Eigen::MatrixXd& designMatrix,
                                          const Eigen::VectorXd& weights,
                                          const std::vector< int >& rowIndices )
{
    const int numberOfParameters = designMatrix.cols( );
    const int blockSize = 4096;

    Eigen::MatrixXd normalMatrix = Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters );
    Eigen::MatrixXd weightedRows;
    for( unsigned int blockStart = 0; blockStart < rowIndices.size( ); blockStart += blockSize )
    {
        const int currentBlockSize = std::min( blockSize, static_cast< int >( rowIndices.size( ) - blockStart ) );
        weightedRows.resize( currentBlockSize, numberOfParameters );
        for( int i = 0; i < currentBlockSize; i++ )
        {
            const int row = rowIndices.at( blockStart + i );
            weightedRows.row( i ) = std::sqrt( weights( row ) ) * designMatrix.row( row );
        }
        normalMatrix.selfadjointView< Eigen::Lower >( ).rankUpdate( weightedRows.transpose( ) );
    }
    return normalMatrix.selfadjointView< Eigen::Lower >( );
}

//! Trace of the product of two matrices, tr(A B), without forming the product.
double traceOfProduct( const Eigen::MatrixXd& first, const Eigen::MatrixXd& second )
{
    return first.cwiseProduct( second.transpose( ) ).sum( );
}

}  // namespace

VarianceComponentEstimate computeHelmertVarianceComponents( const Eigen::MatrixXd& designMatrix,
                                                            const Eigen::VectorXd& weights,
                                                            const Eigen::MatrixXd& inverseCovariance,
                                                            const Eigen::VectorXd& residuals,
                                                            const std::vector< std::vector< int > >& groupRowIndices,
                                                            const double minimumRedundancy )
{
    const int numberOfObservations = designMatrix.rows( );
    const int numberOfParameters = designMatrix.cols( );
    const int numberOfGroups = groupRowIndices.size( );

    if( weights.rows( ) != numberOfObservations || residuals.rows( ) != numberOfObservations )
    {
        throw std::runtime_error( "Error when computing variance components: design matrix has " + std::to_string( numberOfObservations ) +
                                  " rows, but there are " + std::to_string( weights.rows( ) ) + " weights and " +
                                  std::to_string( residuals.rows( ) ) + " residuals." );
    }
    if( inverseCovariance.rows( ) != numberOfParameters || inverseCovariance.cols( ) != numberOfParameters )
    {
        throw std::runtime_error( "Error when computing variance components: inverse covariance has size " +
                                  std::to_string( inverseCovariance.rows( ) ) + "x" + std::to_string( inverseCovariance.cols( ) ) +
                                  ", expected " + std::to_string( numberOfParameters ) + "x" + std::to_string( numberOfParameters ) +
                                  ". Linear parameter constraints are not supported." );
    }
    if( numberOfGroups == 0 )
    {
        throw std::runtime_error( "Error when computing variance components: no observation groups provided." );
    }

    // Check the row indices and the weights.
    std::vector< bool > isRowInGroup( numberOfObservations, false );
    for( int k = 0; k < numberOfGroups; k++ )
    {
        for( int row : groupRowIndices.at( k ) )
        {
            if( row < 0 || row >= numberOfObservations )
            {
                throw std::runtime_error( "Error when computing variance components: row index " + std::to_string( row ) + " of group " +
                                          std::to_string( k ) + " is out of range." );
            }
            if( isRowInGroup.at( row ) )
            {
                throw std::runtime_error( "Error when computing variance components: row " + std::to_string( row ) +
                                          " is in more than one group." );
            }
            if( !( weights( row ) >= 0.0 ) )
            {
                throw std::runtime_error( "Error when computing variance components: weight of row " + std::to_string( row ) +
                                          " is negative or not a number." );
            }
            isRowInGroup.at( row ) = true;
        }
    }

    // One factorization of the full normal matrix N, used for all products N^-1 N_k.
    Eigen::LDLT< Eigen::MatrixXd > normalMatrixFactorization( inverseCovariance );
    if( normalMatrixFactorization.info( ) != Eigen::Success || !normalMatrixFactorization.isPositive( ) )
    {
        throw std::runtime_error( "Error when computing variance components: the normal matrix is not positive definite." );
    }

    // Per group: N^-1 N_k, n_k, c_k = v_k^T W_k v_k and the residual rms.
    VarianceComponentEstimate estimate;
    estimate.groupResults.resize( numberOfGroups );
    std::vector< Eigen::MatrixXd > solvedGroupNormalMatrices( numberOfGroups );
    Eigen::MatrixXd fixedNormalMatrix = inverseCovariance;
    for( int k = 0; k < numberOfGroups; k++ )
    {
        VarianceComponentGroupResult& groupResult = estimate.groupResults.at( k );
        double squaredResidualSum = 0.0;
        for( int row : groupRowIndices.at( k ) )
        {
            if( weights( row ) > 0.0 )
            {
                groupResult.numberOfObservations++;
                groupResult.weightedSquareSum += weights( row ) * residuals( row ) * residuals( row );
                squaredResidualSum += residuals( row ) * residuals( row );
            }
        }
        if( groupResult.numberOfObservations > 0 )
        {
            groupResult.rmsResidual = std::sqrt( squaredResidualSum / groupResult.numberOfObservations );
        }

        const Eigen::MatrixXd groupNormalMatrix = computeGroupNormalMatrix( designMatrix, weights, groupRowIndices.at( k ) );
        fixedNormalMatrix -= groupNormalMatrix;
        solvedGroupNormalMatrices.at( k ) = normalMatrixFactorization.solve( groupNormalMatrix );

        // Redundancy n_k - tr(N^-1 N_k), and the simplified (Foerstner) variance factor.
        groupResult.redundancy = groupResult.numberOfObservations - solvedGroupNormalMatrices.at( k ).trace( );
        groupResult.isFrozen = !( groupResult.redundancy >= minimumRedundancy );
        groupResult.simplifiedVarianceFactor = groupResult.isFrozen ? 1.0 : groupResult.weightedSquareSum / groupResult.redundancy;
    }

    // The part of N that is not in any group (a priori information, soft constraints, ungrouped observations) has a
    // known variance factor of 1.
    const Eigen::MatrixXd solvedFixedNormalMatrix = normalMatrixFactorization.solve( fixedNormalMatrix );

    // Groups whose variance factor is estimated.
    std::vector< int > activeGroups;
    for( int k = 0; k < numberOfGroups; k++ )
    {
        if( !estimate.groupResults.at( k ).isFrozen )
        {
            activeGroups.push_back( k );
        }
    }
    const int numberOfActiveGroups = activeGroups.size( );

    // Helmert's equations, Sahin et al. (1992) Eqs. 15-19. Fixed components (the part outside any group, and frozen
    // groups) have a variance factor of 1, so their columns are moved to the right-hand side.
    estimate.helmertMatrix = Eigen::MatrixXd::Zero( numberOfActiveGroups, numberOfActiveGroups );
    Eigen::VectorXd rightHandSide = Eigen::VectorXd::Zero( numberOfActiveGroups );
    for( int i = 0; i < numberOfActiveGroups; i++ )
    {
        const int k = activeGroups.at( i );
        const Eigen::MatrixXd& solvedGroupNormalMatrix = solvedGroupNormalMatrices.at( k );
        const double groupTrace = solvedGroupNormalMatrix.trace( );

        for( int j = i; j < numberOfActiveGroups; j++ )
        {
            const double coupling = traceOfProduct( solvedGroupNormalMatrix, solvedGroupNormalMatrices.at( activeGroups.at( j ) ) );
            estimate.helmertMatrix( i, j ) = coupling;
            estimate.helmertMatrix( j, i ) = coupling;
        }
        estimate.helmertMatrix( i, i ) += estimate.groupResults.at( k ).numberOfObservations - 2.0 * groupTrace;

        rightHandSide( i ) =
                estimate.groupResults.at( k ).weightedSquareSum - traceOfProduct( solvedGroupNormalMatrix, solvedFixedNormalMatrix );
        for( int l = 0; l < numberOfGroups; l++ )
        {
            if( estimate.groupResults.at( l ).isFrozen )
            {
                rightHandSide( i ) -= traceOfProduct( solvedGroupNormalMatrix, solvedGroupNormalMatrices.at( l ) );
            }
        }
    }

    Eigen::VectorXd helmertVarianceFactors = Eigen::VectorXd::Zero( numberOfActiveGroups );
    if( numberOfActiveGroups > 0 )
    {
        Eigen::ColPivHouseholderQR< Eigen::MatrixXd > helmertSolver( estimate.helmertMatrix );
        if( helmertSolver.isInvertible( ) )
        {
            helmertVarianceFactors = helmertSolver.solve( rightHandSide );
        }
        else
        {
            helmertVarianceFactors.setConstant( std::numeric_limits< double >::quiet_NaN( ) );
        }
    }

    // Helmert's method does not force the factors to be positive. Use the simplified estimator if any factor is not.
    for( int i = 0; i < numberOfActiveGroups; i++ )
    {
        if( !( helmertVarianceFactors( i ) > 0.0 ) || !std::isfinite( helmertVarianceFactors( i ) ) )
        {
            estimate.simplifiedEstimatorUsed = true;
        }
    }

    for( int i = 0; i < numberOfActiveGroups; i++ )
    {
        VarianceComponentGroupResult& groupResult = estimate.groupResults.at( activeGroups.at( i ) );
        groupResult.varianceFactor = estimate.simplifiedEstimatorUsed ? groupResult.simplifiedVarianceFactor : helmertVarianceFactors( i );
    }

    return estimate;
}

}  // namespace simulation_setup

}  // namespace tudat
