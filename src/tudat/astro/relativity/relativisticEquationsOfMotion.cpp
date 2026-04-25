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

#include <tudat/astro/basic_astro/physicalConstants.h>
#include <stdexcept>

#include "tudat/astro/relativity/relativisticEquationsOfMotion.h"

namespace tudat
{

namespace relativity
{

Eigen::Vector4d evaluateFourAcceleration(
        const std::vector< Eigen::Matrix4d >& christoffelSymbols, const Eigen::Vector4d& fourVelocity )
{
    Eigen::Vector4d accelerationVector = Eigen::Vector4d::Zero( );

    for( int i = 0; i < 4; i++ )
    {
        accelerationVector( i ) = fourVelocity.transpose( ) * christoffelSymbols.at( i ) * fourVelocity;
    }

    return accelerationVector;
}

Eigen::Vector4d evaluateFourAcceleration(
         const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 8, 1 >& currentFourState )
{
    return evaluateFourAcceleration(
                spaceTimeMetric->getCurrentChristoffelSymbols( ), currentFourState.segment( 4, 4 ) );
}

Eigen::Vector4d evaluateFourAccelerationWithUpdate(
       const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 8, 1 >& currentFourState,
       const double currentTime )
{

    Eigen::Matrix< double, 6, 1 > currentState;
    currentState.segment( 0, 3 ) = currentFourState.segment( 1, 3 );
    currentState.segment( 3, 3 ) = currentFourState.segment( 5, 3 );

    spaceTimeMetric->update( currentState, currentTime, 0, 1 );

    return evaluateFourAcceleration( spaceTimeMetric, currentFourState );
}



Eigen::Vector3d evaluateRelativisticEquationsOfMotionInCoordinateTime(
        const std::vector< Eigen::Matrix4d >& christoffelSymbols, const Eigen::Vector3d& coordinateVelocity )
{
    Eigen::Vector4d spaceTimeCoordinateVelocity;
    spaceTimeCoordinateVelocity( 0 ) = physical_constants::SPEED_OF_LIGHT;
    spaceTimeCoordinateVelocity.segment( 1, 3 ) = coordinateVelocity;

    Eigen::Vector3d accelerationVector = Eigen::Vector3d::Zero( );
    const double timeChristoffelSymbolProduct =
            ( spaceTimeCoordinateVelocity.transpose( ) * christoffelSymbols.at( 0 ) * spaceTimeCoordinateVelocity )( 0, 0 ) /
            physical_constants::SPEED_OF_LIGHT;

    for( int i = 0; i < 3; i++ )
    {
        const double firstTerm = spaceTimeCoordinateVelocity( i + 1 ) * timeChristoffelSymbolProduct;
        const double secondTerm =
                ( spaceTimeCoordinateVelocity.transpose( ) * christoffelSymbols.at( i + 1 ) * spaceTimeCoordinateVelocity )( 0, 0 );
        accelerationVector( i ) = firstTerm - secondTerm;
    }

    return accelerationVector;
}

Eigen::Vector3d evaluateRelativisticEquationsOfMotionInCoordinateTime(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState )
{
    return evaluateRelativisticEquationsOfMotionInCoordinateTime(
                spaceTimeMetric->getCurrentChristoffelSymbols( ), currentState.segment( 3, 3 ) );
}

Eigen::Vector3d evaluateRelativisticEquationsOfMotionInCoordinateTimeWithUpdate(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState,
        const double currentTime )
{
    spaceTimeMetric->update( currentState, currentTime, 1, 1 );
    return evaluateRelativisticEquationsOfMotionInCoordinateTime( spaceTimeMetric, currentState );
}

double evaluateProperTimeEquation(
        const Eigen::Matrix4d& covariantMetricPerturbation, const Eigen::Vector3d& coordinateVelocity, const int squareRootOrderExpansion )
{

    Eigen::Vector4d spaceTimeCoordinateVelocity;
    spaceTimeCoordinateVelocity( 0 ) = physical_constants::SPEED_OF_LIGHT;
    spaceTimeCoordinateVelocity.segment( 1, 3 ) = coordinateVelocity;

    double coordinateSpeed = coordinateVelocity.norm( );
    double squareRootPerturbation = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
            ( spaceTimeCoordinateVelocity.transpose( ) * covariantMetricPerturbation * spaceTimeCoordinateVelocity +
            coordinateSpeed * coordinateSpeed );

    double properTimeRate = 0.0;

    for( int i = 1; i <= squareRootOrderExpansion; i++ )
    {
        switch( i )
        {
        case 1:
            properTimeRate -= squareRootPerturbation * 0.5;
            break;
        case 2:
            properTimeRate -= squareRootPerturbation * squareRootPerturbation * 0.125;
            break;
        case 3:
            properTimeRate -= squareRootPerturbation  * squareRootPerturbation * squareRootPerturbation * 0.0625;
            break;
        default:
            throw std::runtime_error(
                        "Error when calculating proper time equation, cannot expand square root at order " +
                        std::to_string( i ) );
        }
    }

    return properTimeRate;
}

double evaluateProperTimeEquation(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState,
        const int squareRootOrderExpansion )
{
    return evaluateProperTimeEquation( spaceTimeMetric->getCurrentCovariantMetricPeturbation(  ), currentState.segment( 3, 3 ),
                                       squareRootOrderExpansion );
}

double evaluateProperTimeEquationWithUpdate(
        const std::shared_ptr< relativity::Metric > spaceTimeMetric, const Eigen::Matrix< double, 6, 1 >& currentState,
        const double currentTime, const int squareRootOrderExpansion )
{
    spaceTimeMetric->update( currentState, currentTime, 1, 0 );
    return evaluateProperTimeEquation( spaceTimeMetric, currentState, squareRootOrderExpansion );
}

}

}
