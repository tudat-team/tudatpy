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
        const double massOfBodyExertingTide, const double k2LoveNumber, const double timeLag, const double referenceRadius,
        const bool includeDirectRadialComponent )
{
    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativeVelocity = relativeStateOfBodyExertingTide.segment( 3, 3 );

    double distance = relativePosition.norm( );

    double distanceSquared = distance * distance;
    double distanceToEighthPower = distanceSquared * distanceSquared * distanceSquared * distanceSquared;
    double referenceRadiusToFifthPower = referenceRadius *  referenceRadius *  referenceRadius *
            referenceRadius *  referenceRadius;

    double radialComponentMultiplier = ( includeDirectRadialComponent == true ) ? 1.0 : 0.0;

    return - 3.0 * massOfBodyExertingTide * referenceRadiusToFifthPower / distanceToEighthPower * (
                 radialComponentMultiplier * k2LoveNumber * relativePosition + k2LoveNumber * timeLag * (
                    2.0 * ( relativePosition.dot( relativeVelocity ) * relativePosition / distanceSquared ) +
                    ( relativePosition.cross( planetAngularVelocityVector ) + relativeVelocity ) ) );

}

Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnSatellite(
        const Eigen::Vector6d relativeStateOfBodyExertingTide,
        const double massOfBodyExertingTide, const double massOfBodyUndergoingTide, const double k2LoveNumber,
        const double timeLag, const double referenceRadius, const bool includeDirectRadialComponent )
{

    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativeVelocity = relativeStateOfBodyExertingTide.segment( 3, 3 );

    double distance = relativePosition.norm( );

    double distanceSquared = distance * distance;
    double distanceToEighthPower = distanceSquared * distanceSquared * distanceSquared * distanceSquared;
    double referenceRadiusToFifthPower = referenceRadius *  referenceRadius *  referenceRadius *
            referenceRadius *  referenceRadius;

    // std::cout<<referenceRadius<<" "<<massOfBodyExertingTide<<std::endl;

    double radialComponentMultiplier = ( includeDirectRadialComponent == true ) ? 1.0 : 0.0;

    return -3.0 * massOfBodyExertingTide * massOfBodyExertingTide / massOfBodyUndergoingTide *
            referenceRadiusToFifthPower / distanceToEighthPower * (
                 radialComponentMultiplier * k2LoveNumber * relativePosition + k2LoveNumber * timeLag * (
                   7.0 * relativePosition.dot( relativeVelocity ) * relativePosition / distanceSquared ) );

}

} // namespace gravitation

} // namespace tudat
