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

#include <iostream>
#include <iomanip>

#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/physicalConstants.h"

#include "tudat/astro/relativity/relativisticPotentials.h"

namespace tudat
{

namespace relativity
{

std::vector< double > calculateFirstOrderExternalScalarPotentials( const std::vector< double >& gravitationalParameters,
                                                                   const std::vector< double >& distancesFromBody,
                                                                   const std::map< int, double >& sphericalHarmonicPotentials )
{
    std::vector< double > externalScalarPotentials;
    for( unsigned int i = 0; i < gravitationalParameters.size( ); i++ )
    {
        if( sphericalHarmonicPotentials.count( i ) != 0 )
        {
            externalScalarPotentials.push_back( sphericalHarmonicPotentials.at( i ) );
        }
        else
        {
            externalScalarPotentials.push_back( gravitationalParameters[ i ] / distancesFromBody[ i ] );
        }
    }
    return externalScalarPotentials;
}

double calculateFirstOrderExternalScalarPotential( const std::vector< double >& gravitationalParameters,
                                                   const std::vector< double >& distancesFromBody,
                                                   const std::map< int, double >& sphericalHarmonicPotentials )
{
    std::vector< double > externalScalarPotentials =
            calculateFirstOrderExternalScalarPotentials( gravitationalParameters, distancesFromBody, sphericalHarmonicPotentials );
    double externalScalarPotential = 0.0;
    for( unsigned int i = 0; i < externalScalarPotentials.size( ); i++ )
    {
        externalScalarPotential += externalScalarPotentials.at( i );
    }
    return externalScalarPotential;
}

double calculateSecondOrderExternalScalarPotentialCorrection( const std::vector< double >& gravitationalParameters,
                                                              const std::vector< double >& distancesFromBody,
                                                              const std::vector< double >& velocities,
                                                              const std::vector< Eigen::Vector6d >& stateVectors,
                                                              const Eigen::Vector3d centralBodyVelocity,
                                                              const std::vector< Eigen::Vector3d >& accelerations )
{
    double externalScalarPotential = 0.0;
    double singleBodyCorrectionTerm = 0.0;
    double velocityPositionProduct;

    for( unsigned int i = 0; i < gravitationalParameters.size( ); i++ )
    {
        singleBodyCorrectionTerm = 0.0;
        velocityPositionProduct = stateVectors[ i ].segment( 0, 3 ).dot( stateVectors[ i ].segment( 3, 3 ) );
        singleBodyCorrectionTerm += -2.0 * velocities[ i ] * velocities[ i ] +
                0.5 * ( velocityPositionProduct * velocityPositionProduct / ( distancesFromBody[ i ] * distancesFromBody[ i ] ) ) +
                4.0 *
                        stateVectors[ i ].segment( 3, 3 ).dot(
                                centralBodyVelocity );  // Final term included in Fienga et al 2009, but no Soffel et al. // +
        // accelerations[ i ].dot( stateVectors[ i ].segment( 0, 3 ) ) );

        for( unsigned int j = 0; j < gravitationalParameters.size( ); j++ )
        {
            if( i != j )
            {
                singleBodyCorrectionTerm += gravitationalParameters[ j ] /
                        ( ( stateVectors[ j ].segment( 0, 3 ) - stateVectors[ i ].segment( 0, 3 ) ).norm( ) );
            }
        }

        singleBodyCorrectionTerm *= gravitationalParameters[ i ] / distancesFromBody[ i ];
        externalScalarPotential += singleBodyCorrectionTerm;
    }
    return externalScalarPotential;
}

std::vector< Eigen::Vector3d > calculateExternalVectorPotentials( const std::vector< double >& gravitationalParameters,
                                                                  const std::vector< double >& distancesFromBody,
                                                                  const std::vector< Eigen::Vector6d >& stateVectors,
                                                                  const std::map< int, Eigen::Vector3d >& currentAngularMomenta )
{
    std::vector< Eigen::Vector3d > vectorPotentials;

    for( unsigned int i = 0; i < gravitationalParameters.size( ); i++ )
    {
        vectorPotentials.push_back( gravitationalParameters[ i ] / distancesFromBody[ i ] * stateVectors[ i ].segment( 3, 3 ) );
    }

    double currentDistance = 0;
    for( std::map< int, Eigen::Vector3d >::const_iterator momentaIterator = currentAngularMomenta.begin( );
         momentaIterator != currentAngularMomenta.end( );
         momentaIterator++ )
    {
        currentDistance = distancesFromBody.at( momentaIterator->first );
        vectorPotentials.at( momentaIterator->first ) -= physical_constants::GRAVITATIONAL_CONSTANT *
                ( momentaIterator->second ).cross( Eigen::Vector3d( stateVectors.at( momentaIterator->first ).segment( 0, 3 ) ) ) /
                ( 2.0 * currentDistance * currentDistance * currentDistance );
    }

    return vectorPotentials;
}

Eigen::Vector3d calculateExternalVectorPotential( const std::vector< double >& gravitationalParameters,
                                                  const std::vector< double >& distancesFromBody,
                                                  const std::vector< Eigen::Vector6d >& stateVectors,
                                                  const std::map< int, Eigen::Vector3d >& currentAngularMomenta )
{
    std::vector< Eigen::Vector3d > externalVectorPotentials =
            calculateExternalVectorPotentials( gravitationalParameters, distancesFromBody, stateVectors, currentAngularMomenta );
    Eigen::Vector3d externalVectorPotential = Eigen::Vector3d::Zero( );
    for( unsigned int i = 0; i < externalVectorPotentials.size( ); i++ )
    {
        externalVectorPotential += externalVectorPotentials.at( i );
    }
    return externalVectorPotential;
}

Eigen::Vector3d calculateExternalVectorPotential( const std::vector< double >& externalBodyScalarPotentials,
                                                  const std::vector< Eigen::Vector6d >& stateVectors,
                                                  const std::map< int, Eigen::Vector3d >& currentAngularMomenta,
                                                  const std::vector< double >& distancesFromBody )
{
    Eigen::Vector3d externalVectorPotential = Eigen::Vector3d::Zero( );
    for( unsigned int i = 0; i < externalBodyScalarPotentials.size( ); i++ )
    {
        externalVectorPotential += externalBodyScalarPotentials.at( i ) * stateVectors.at( i ).segment( 3, 3 );
    }

    for( std::map< int, Eigen::Vector3d >::const_iterator momentaIterator = currentAngularMomenta.begin( );
         momentaIterator != currentAngularMomenta.end( );
         momentaIterator++ )
    {
        double currentDistance = distancesFromBody.at( momentaIterator->first );
        externalVectorPotential -= physical_constants::GRAVITATIONAL_CONSTANT *
                ( momentaIterator->second ).cross( Eigen::Vector3d( stateVectors.at( momentaIterator->first ).segment( 0, 3 ) ) ) /
                ( 2.0 * currentDistance * currentDistance * currentDistance );
    }
    return externalVectorPotential;
}

Eigen::Matrix3d calculateExternalMatrixPotential( const double gravitationalParameter,
                                                  const double distanceFromBody,
                                                  const Eigen::Vector3d& relativePosition )
{
    Eigen::Matrix3d matrixPotential = ( gravitationalParameter * gravitationalParameter ) / ( distanceFromBody * distanceFromBody * 4.0 ) *
            ( relativePosition * relativePosition.transpose( ) / ( distanceFromBody * distanceFromBody ) - Eigen::Matrix3d::Identity( ) );
    return matrixPotential;
}

Eigen::Matrix3d calculateExternalMatrixPotential( const std::vector< double >& gravitationalParameters,
                                                  const std::vector< double >& distancesFromBody,
                                                  const std::vector< Eigen::Vector3d >& relativePositions )
{
    Eigen::Matrix3d matrixPotential = Eigen::Matrix3d::Zero( );

    for( unsigned int i = 0; i < gravitationalParameters.size( ); i++ )
    {
        matrixPotential +=
                calculateExternalMatrixPotential( gravitationalParameters.at( i ), distancesFromBody.at( i ), relativePositions.at( i ) );
    }
    return matrixPotential;
}

}  // namespace relativity

}  // namespace tudat
