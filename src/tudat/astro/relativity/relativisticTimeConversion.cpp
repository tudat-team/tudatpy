/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include <iostream>
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/relativity/relativisticTimeConversion.h"
#include <iomanip>

namespace tudat
{

namespace relativity
{

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential
double calculateFirstCentralBodyProperTimeRateDifference( const double computationPointSpeed,
                                                          double gravitationalScalarPotential,
                                                          const double equivalencePrincipleLpiViolationParameter )
{
    return ( -( 0.5 * computationPointSpeed * computationPointSpeed +
                ( 1.0 + equivalencePrincipleLpiViolationParameter ) * gravitationalScalarPotential ) ) *
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential
double calculateFirstCentralBodyProperTimeRateDifference( const Eigen::Vector6d& computationPointState,
                                                          const std::vector< Eigen::Vector6d >& perturbedInertialStates,
                                                          const std::vector< double >& centralBodyGravitationalParameters,
                                                          const double equivalencePrincipleLpiViolationParameter )
{
    double gravitationalScalarPotential = 0.0;
    for( unsigned int i = 0; i < perturbedInertialStates.size( ); i++ )
    {
        gravitationalScalarPotential += centralBodyGravitationalParameters.at( i ) /
                ( perturbedInertialStates.at( i ).segment( 0, 3 ) - computationPointState.segment( 0, 3 ) ).norm( );
    }

    return calculateFirstCentralBodyProperTimeRateDifference(
            computationPointState.segment( 3, 3 ).norm( ), gravitationalScalarPotential, equivalencePrincipleLpiViolationParameter );
}

//! Function to calculate the 1st order integrand term in Eq. (58) of Soffel et al. (2003)
double calculateFirstOrderTcbToTcgIntegrand(
        const Eigen::Vector3d velocityVector, const double gravitationalScalarPotential )
{
    return calculateFirstOrderTcbToTcgIntegrand( velocityVector.norm( ),
                                                 gravitationalScalarPotential );
}

//! Function to calculate the 1st order integrand term in Eq. (58) of Soffel et al. (2003)
double calculateFirstOrderTcbToTcgIntegrand(
        const double barycentricSpeed, const double gravitationalScalarPotential )
{
    return ( -( 0.5 * barycentricSpeed * barycentricSpeed + gravitationalScalarPotential )  ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}

double calculateFirstOrderTcbToTcgDirectCorrection( const Eigen::Vector3d& geocentricPositionVector,
                                                    const Eigen::Vector3d& barycentricVelocityVector )
{
    return -barycentricVelocityVector.dot( geocentricPositionVector ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}


double calculateSecondOrderTcbToTcgIntegrand(
        const double gravitationalScalarPotential,
        const Eigen::Vector3d& barycentricVelocity,
        const Eigen::Vector3d& gravitationalVectorPotential,
        const double currentSecondOrderExternalPotentialCorrection )
{
    double barycentricSpeed = barycentricVelocity.norm( );
    return calculateSecondOrderTcbToTcgIntegrand(
                barycentricSpeed, gravitationalScalarPotential,
                barycentricVelocity, gravitationalVectorPotential, currentSecondOrderExternalPotentialCorrection );
}

double calculateSecondOrderTcbToTcgIntegrand(
        const double barycentricSpeed,
        const double gravitationalScalarPotential,
        const Eigen::Vector3d& barycentricVelocity,
        const Eigen::Vector3d& gravitationalVectorPotential,
        const double currentSecondOrderExternalPotentialCorrection )
{
    double barycentricSpeedSquared = barycentricSpeed * barycentricSpeed;

    return ( -0.125 * barycentricSpeedSquared * barycentricSpeedSquared -
             1.5 * gravitationalScalarPotential * barycentricSpeedSquared +
             4.0 * barycentricVelocity.dot( gravitationalVectorPotential ) +
             0.5 * gravitationalScalarPotential * gravitationalScalarPotential + currentSecondOrderExternalPotentialCorrection ) *
            physical_constants::INVERSE_QUARTIC_SPEED_OF_LIGHT;
}

double calculateSecondOrderTcbToTcgDirectCorrection( const Eigen::Vector3d& geocentricPositionVector,
                                                     const Eigen::Vector3d& barycentricVelocityVector,
                                                    const  double gravitationalScalarPotential,
                                                     const Eigen::Vector3d& gravitationalVectorPotential )
{
    double barycentricSpeed = barycentricVelocityVector.norm( );
    return physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * ( (
                -barycentricVelocityVector.dot( geocentricPositionVector ) * (
                    1.0 + physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * (
                        0.5 * barycentricSpeed * barycentricSpeed + 3.0 * gravitationalScalarPotential ) ) ) +
                                                                 physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                                                                 4.0 * gravitationalVectorPotential.dot( geocentricPositionVector ) );
}

double calculateFirstOrderPlanetocentricToTopocentricConversion(
        const Eigen::Vector3d inertialPlanetCenteredPointPosition,
        const Eigen::Vector3d inertialPointVelocity,
        const double localCentralBodyPotential,
        const Eigen::Vector3d centralBodyAcceleration,
        const std::vector< Eigen::Vector3d >& relativeTidalBodyPositions,
        const std::vector< double >& tidalBodyGravitationalParameters )
{
    double firstOrderTerm = 0.0;
    firstOrderTerm += inertialPointVelocity.norm( ) * inertialPointVelocity.norm( ) / 2.0 + localCentralBodyPotential +
            centralBodyAcceleration.dot( inertialPlanetCenteredPointPosition );


    double distance, dotProduct, pointPosition;
    for( unsigned int i = 0; i < relativeTidalBodyPositions.size( ); i++ )
    {
        distance = relativeTidalBodyPositions[ i ].norm( );
        dotProduct = inertialPlanetCenteredPointPosition.dot( relativeTidalBodyPositions[ i ].normalized( ) );
        pointPosition = inertialPlanetCenteredPointPosition.norm( );
        firstOrderTerm += tidalBodyGravitationalParameters[ i ] / ( 2.0 * distance * distance * distance ) *
                ( 3.0 * dotProduct * dotProduct - pointPosition * pointPosition );
    }

    return - firstOrderTerm * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}

}  // namespace relativity

}  // namespace tudat
