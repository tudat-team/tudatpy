/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/relativity/schwarzschildMetric.h"
namespace tudat
{

namespace relativity
{

void HarmonicSchwarzschildMetric::update( const Eigen::Matrix< double, 6, 1 >& state,
                                          const double time,
                                          const bool updateCurrentMetric,
                                          const bool updateCurrentChristoffelSymbols )
{
    currentEvaluationState_ = state;
    currentTime_ = time;
    currentCentralBodyGravitationalParameter_ = centralGravitationalParameterFunction_( );

    if ( updateCurrentMetric )
    {
        updateMetric( );
    }

    if ( updateCurrentChristoffelSymbols )
    {
        updateChristoffelSymbols( );
    }
}

void HarmonicSchwarzschildMetric::updateMetric( )
{
    const double distance = currentEvaluationState_.segment( 0, 3 ).norm( );
    const double scaledPotential = currentCentralBodyGravitationalParameter_ *
                                   physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / distance;

    currentFirstOrderCovariantMetricContributions_.setZero( );
    currentFirstOrderCovariantMetricContributions_( 0, 0 ) = 2.0 * scaledPotential;
    currentFirstOrderCovariantMetricContributions_.block( 1, 1, 3, 3 ) =
        2.0 * ppnParameterSet_->getParameterGamma( ) * scaledPotential * Eigen::Matrix3d::Identity( );

    currentSecondOrderCovariantMetricContributions_.setZero( );
    currentSecondOrderCovariantMetricContributions_( 0, 0 ) =
        -2.0 * ppnParameterSet_->getParameterBeta( ) * scaledPotential * scaledPotential;

    if ( includeSecondPostNewtonianOrder_ )
    {
        currentSecondOrderCovariantMetricContributions_.block( 1, 1, 3, 3 ) +=
            2.0 * ppnParameterSet_->getParameterEpsilon( ) * scaledPotential * scaledPotential * Eigen::Matrix3d::Identity( ) +
            4.0 * calculateExternalMatrixPotential(
                      currentCentralBodyGravitationalParameter_,
                      distance,
                      currentEvaluationState_.segment( 0, 3 ) ) * physical_constants::INVERSE_QUARTIC_SPEED_OF_LIGHT;
    }

    currentCovariantMetricContribution_ =
        currentFirstOrderCovariantMetricContributions_ + currentSecondOrderCovariantMetricContributions_;
}

void HarmonicSchwarzschildMetric::updateChristoffelSymbols( )
{
    const double gamma = ppnParameterSet_->getParameterGamma( );
    const double beta = ppnParameterSet_->getParameterBeta( );
    const double epsilon = includeSecondPostNewtonianOrder_ ? ppnParameterSet_->getParameterEpsilon( ) : 0.0;

    for ( int i = 0; i < 4; i++ )
    {
        currentChristoffelSymbols_[ i ].setZero( );
    }

    const Eigen::Vector3d position = currentEvaluationState_.segment( 0, 3 );

    const double distance = position.norm( );
    const double inverseDistanceSquared = 1.0 / ( distance * distance );
    const double scaledPotential = currentCentralBodyGravitationalParameter_ *
                                   physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / distance;

    // Placeholder: assume scalarPotentialTimeDerivative = 0 (can be replaced with dynamic model)
    const double scalarPotentialTimeDerivative = 0.0;

    // ^0_{0i} and ^0_{i0}
    currentChristoffelSymbols_[ 0 ]( 0, 0 ) = -scalarPotentialTimeDerivative;
    const Eigen::Vector3d spaceTimeCoupling = position * ( 1.0 + 2.0 * ( 1.0 - beta ) * scaledPotential );
    currentChristoffelSymbols_[ 0 ].block( 1, 0, 3, 1 ) = spaceTimeCoupling;
    currentChristoffelSymbols_[ 0 ].block( 0, 1, 1, 3 ) = spaceTimeCoupling.transpose( );
    currentChristoffelSymbols_[ 0 ].block( 1, 1, 3, 3 ) =
        scalarPotentialTimeDerivative * gamma * Eigen::Matrix3d::Identity( );

    // ^i_{00}
    const double timeTimeCommonTerm = 1.0 - 2.0 * ( gamma + beta ) * scaledPotential;
    for ( int i = 1; i < 4; i++ )
    {
        currentChristoffelSymbols_[ i ]( 0, 0 ) = position( i - 1 ) * timeTimeCommonTerm;
    }

    // ^i_{jk}
    for ( int i = 0; i < 3; i++ )
    {
        Eigen::Matrix3d spatialBlock = Eigen::Matrix3d::Zero( );

        for ( int j = 0; j < 3; j++ )
        {
            for ( int k = 0; k < 3; k++ )
            {
                spatialBlock( j, k ) =
                        gamma *
                        ( ( j == k ? position( i ) : 0.0 )
                        - ( i == j ? position( k ) : 0.0 )
                        - ( i == k ? position( j ) : 0.0 ) );
            }
        }

        if ( includeSecondPostNewtonianOrder_ )
        {
            for ( int j = 0; j < 3; j++ )
            {
                for ( int k = 0; k < 3; k++ )
                {
                    spatialBlock( j, k ) +=
                            2.0 * scaledPotential * ( epsilon - gamma * gamma ) *
                            ( j == k ? position( i ) : 0.0 );

                    spatialBlock( j, k ) +=
                            -scaledPotential * ( epsilon - 2.0 * gamma * gamma ) *
                            ( ( i == j ? position( k ) : 0.0 )
                            + ( i == k ? position( j ) : 0.0 ) );
                }
            }
        }

        currentChristoffelSymbols_[ i + 1 ].block( 1, 1, 3, 3 ) += spatialBlock;

        currentChristoffelSymbols_[ i + 1 ]( 0, i + 1 ) = gamma * scalarPotentialTimeDerivative;
        currentChristoffelSymbols_[ i + 1 ]( i + 1, 0 ) = gamma * scalarPotentialTimeDerivative;
    }

    if ( includeSecondPostNewtonianOrder_ )
    {
        Eigen::Matrix3d secondOrderCorrection =
            -position * position.transpose( ) * inverseDistanceSquared * epsilon * scaledPotential;
        for ( int i = 0; i < 3; i++ )
        {
            currentChristoffelSymbols_[ i + 1 ].block( 1, 1, 3, 3 ) += position( i ) * secondOrderCorrection;
        }
    }

    const double multiplier = scaledPotential * inverseDistanceSquared;
    for ( int i = 0; i < 4; i++ )
    {
        currentChristoffelSymbols_[ i ] *= multiplier;
    }
}

} // namespace relativity

} // namespace tudat
