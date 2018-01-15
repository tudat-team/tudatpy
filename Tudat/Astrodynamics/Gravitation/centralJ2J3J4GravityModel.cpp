/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <cmath>
#include <stdexcept>


#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2J3J4GravityModel.h"

namespace tudat
{
namespace gravitation
{

//! Compute gravitational acceleration due to J4.
Eigen::Vector3d computeGravitationalAccelerationDueToJ4(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double equatorialRadiusOfBodyExertingAcceleration,
        const double j4CoefficientOfGravityField,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    // Set constant values reused for optimal computation of acceleration components.
    const double distanceBetweenBodies = ( positionOfBodySubjectToAcceleration
                                           - positionOfBodyExertingAcceleration ).norm( );

    const double preMultiplier = gravitationalParameterOfBodyExertingAcceleration
            / std::pow( distanceBetweenBodies, 6.0 ) * 4.375 * j4CoefficientOfGravityField
            * std::pow( equatorialRadiusOfBodyExertingAcceleration, 4.0 );

    const double scaledZCoordinate = ( positionOfBodySubjectToAcceleration.z( )
                                       - positionOfBodyExertingAcceleration.z( ) )
            / distanceBetweenBodies;

    const double scaledZCoordinateSquared = scaledZCoordinate * scaledZCoordinate;

    const double scaledZCoordinateToPower4 = scaledZCoordinateSquared * scaledZCoordinateSquared;

    const double factorForXAndYDirections = ( 3.0 / 7.0 - 6.0 * scaledZCoordinateSquared
                                              + 9.0 * scaledZCoordinateToPower4 )
            / distanceBetweenBodies;

    // Compute components of acceleration due to J4-effect.
    Eigen::Vector3d gravitationalAccelerationDueToJ4 = Eigen::Vector3d::Constant( preMultiplier );

    gravitationalAccelerationDueToJ4( orbital_element_conversions::xCartesianPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.x( )
                 - positionOfBodyExertingAcceleration.x( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ4( orbital_element_conversions::yCartesianPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.y( )
                 - positionOfBodyExertingAcceleration.y( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ4( orbital_element_conversions::zCartesianPositionIndex )
            *= ( 15.0 / 7.0 - 10.0 * scaledZCoordinateSquared + 9.0 * scaledZCoordinateToPower4 )
            * scaledZCoordinate;

    return gravitationalAccelerationDueToJ4;
}

//! Compute gravitational acceleration zonal sum.
Eigen::Vector3d computeGravitationalAccelerationZonalSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double equatorialRadiusOfBodyExertingAcceleration,
        const std::map< int, double > zonalCoefficientsOfGravityField,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    // Check that only coefficients for the gravity field are given up to J2 (i.e., size of
    // vector is 3 at max), else throw an error.
    if ( zonalCoefficientsOfGravityField.size( ) > 3 )
    {
        throw std::runtime_error(
                            "Currently, accelerations can only be computed up to J4."  );
    }

    // Check if position of body subject to acceleration falls within the effective radius of
    // the central body. If so, throw an error.
    if ( ( positionOfBodySubjectToAcceleration - positionOfBodyExertingAcceleration ).norm( )
         < equatorialRadiusOfBodyExertingAcceleration )
    {
        throw std::runtime_error(
                            "Position of body subject to acceleration is within effective radius." );
    }

    // Set gravitational acceleration sum equal to central term contribution.
    Eigen::Vector3d gravitationalAccelerationSum
            = computeGravitationalAcceleration(
                positionOfBodySubjectToAcceleration,
                gravitationalParameterOfBodyExertingAcceleration,
                positionOfBodyExertingAcceleration );

    // Add contributions from zonal terms. The switch statement checks
    for ( std::map< int, double >::const_iterator mapIterator
          = zonalCoefficientsOfGravityField.begin( );
          mapIterator != zonalCoefficientsOfGravityField.end( ); mapIterator++ )
    {
        switch ( mapIterator->first )
        {
        case 2:

            gravitationalAccelerationSum += computeGravitationalAccelerationDueToJ2(
                        positionOfBodySubjectToAcceleration,
                        gravitationalParameterOfBodyExertingAcceleration,
                        equatorialRadiusOfBodyExertingAcceleration,
                        mapIterator->second,
                        positionOfBodyExertingAcceleration );

            break;

        case 3:

            gravitationalAccelerationSum += computeGravitationalAccelerationDueToJ3(
                        positionOfBodySubjectToAcceleration,
                        gravitationalParameterOfBodyExertingAcceleration,
                        equatorialRadiusOfBodyExertingAcceleration,
                        mapIterator->second,
                        positionOfBodyExertingAcceleration );

            break;

        case 4:

            gravitationalAccelerationSum += computeGravitationalAccelerationDueToJ4(
                        positionOfBodySubjectToAcceleration,
                        gravitationalParameterOfBodyExertingAcceleration,
                        equatorialRadiusOfBodyExertingAcceleration,
                        mapIterator->second,
                        positionOfBodyExertingAcceleration );

            break;

        default:

            throw std::runtime_error(
                                "Degree must be 2, 3, or 4 in current implementation." );
        };
    }

    // Return total gravitational acceleration computed.
    return gravitationalAccelerationSum;
}

//! Get gravitational acceleration.
Eigen::Vector3d CentralJ2J3J4GravitationalAccelerationModel::getAcceleration( )
{
    // Sum and return constituent acceleration terms.
    return computeGravitationalAcceleration(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->positionOfBodyExertingAcceleration )
            + computeGravitationalAccelerationDueToJ2(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->equatorialRadius,
                this->j2GravityCoefficient,
                this->positionOfBodyExertingAcceleration )
            + computeGravitationalAccelerationDueToJ3(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->equatorialRadius,
                this->j3GravityCoefficient,
                this->positionOfBodyExertingAcceleration )
            + computeGravitationalAccelerationDueToJ4(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->equatorialRadius,
                this->j4GravityCoefficient,
                this->positionOfBodyExertingAcceleration );
}

} // namespace gravitation
} // namespace tudat
