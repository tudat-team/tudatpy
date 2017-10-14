/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/Gravitation/directTidalDissipationAcceleration.h"

namespace tudat
{

namespace gravitation
{

Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnPlanet(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const Eigen::Vector3d planetAngularVelocityVector,
        const double currentTidalAccelerationMultiplier, const double timeLag,
        const bool includeDirectRadialComponent )
{
    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativeVelocity = relativeStateOfBodyExertingTide.segment( 3, 3 );

    double distance = relativePosition.norm( );

    double distanceSquared = distance * distance;

    double radialComponentMultiplier = ( includeDirectRadialComponent == true ) ? 1.0 : 0.0;

    return currentTidalAccelerationMultiplier * (
                 radialComponentMultiplier * relativePosition + timeLag * (
                    2.0 * ( relativePosition.dot( relativeVelocity ) * relativePosition / distanceSquared ) +
                    ( relativePosition.cross( planetAngularVelocityVector ) + relativeVelocity ) ) );

}

Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnSatellite(
        const Eigen::Vector6d relativeStateOfBodyExertingTide,
        const double currentTidalAccelerationMultiplier,
        const double timeLag,const bool includeDirectRadialComponent )
{

    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativeVelocity = relativeStateOfBodyExertingTide.segment( 3, 3 );

    double distance = relativePosition.norm( );
    double distanceSquared = distance * distance;

    double radialComponentMultiplier = ( includeDirectRadialComponent == true ) ? 1.0 : 0.0;

    return currentTidalAccelerationMultiplier * (
                 2.0 * radialComponentMultiplier * relativePosition + timeLag * (
                   7.0 * relativePosition.dot( relativeVelocity ) * relativePosition / distanceSquared ) );

}



} // namespace gravitation

} // namespace tudat
