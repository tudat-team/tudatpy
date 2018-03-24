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

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h"

namespace tudat
{
namespace gravitation
{

//! Compute gravitational acceleration due to J2.
Eigen::Vector3d computeGravitationalAccelerationDueToJ2(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double equatorialRadiusOfBodyExertingAcceleration,
        const double j2CoefficientOfGravityField,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    // Set constant values reused for optimal computation of acceleration components.
    const double distanceBetweenBodies = ( positionOfBodySubjectToAcceleration
                                           - positionOfBodyExertingAcceleration ).norm( );

    const double preMultiplier = -gravitationalParameterOfBodyExertingAcceleration
            / std::pow( distanceBetweenBodies, 4.0 ) * 1.5 * j2CoefficientOfGravityField
            * equatorialRadiusOfBodyExertingAcceleration
            * equatorialRadiusOfBodyExertingAcceleration;

    const double scaledZCoordinate = ( positionOfBodySubjectToAcceleration.z( )
                                       - positionOfBodyExertingAcceleration.z( ) )
            / distanceBetweenBodies;

    const double scaledZCoordinateSquared = scaledZCoordinate * scaledZCoordinate;

    const double factorForXAndYDirections = ( 1.0 - 5.0 * scaledZCoordinateSquared )
            / distanceBetweenBodies;

    // Compute components of acceleration due to J2-effect.
    Eigen::Vector3d gravitationalAccelerationDueToJ2 = Eigen::Vector3d::Constant( preMultiplier );

    gravitationalAccelerationDueToJ2( orbital_element_conversions::xCartesianPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.x( )
                 - positionOfBodyExertingAcceleration.x( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ2( orbital_element_conversions::yCartesianPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.y( )
                 - positionOfBodyExertingAcceleration.y( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ2( orbital_element_conversions::zCartesianPositionIndex )
            *= ( 3.0 - 5.0 * scaledZCoordinateSquared ) * scaledZCoordinate;

    return gravitationalAccelerationDueToJ2;
}

//! Get gravitational acceleration.
Eigen::Vector3d CentralJ2GravitationalAccelerationModel::getAcceleration( )
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
                this->positionOfBodyExertingAcceleration );
}

} // namespace gravitation
} // namespace tudat
