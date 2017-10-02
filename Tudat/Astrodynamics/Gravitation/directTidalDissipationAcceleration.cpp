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

Eigen::Vector3d computeDirectTidalDissipationAcceleration(
        const Eigen::Vector6d satelliteRelativeState, const Eigen::Vector3d planetAngularVelocityVector,
        const double satelliteMass, const double k2LoveNumber, const double timeLag, const double planetReferenceRadius,
        const bool includeDirectRadialComponent )
{
    Eigen::Vector3d relativePosition = satelliteRelativeState.segment( 0, 3 );
    Eigen::Vector3d relativeVelocity = satelliteRelativeState.segment( 3, 3 );

    double distance = relativePosition.norm( );

    double distanceSquared = distance * distance;
    double distanceToEighthPower = distanceSquared * distanceSquared * distanceSquared * distanceSquared;
    double referenceRadiusToFifthPower = planetReferenceRadius *  planetReferenceRadius *  planetReferenceRadius *
            planetReferenceRadius *  planetReferenceRadius;

    double radialComponentMultiplier = ( includeDirectRadialComponent == true ) ? 1.0 : 0.0;

    return - 3.0 * satelliteMass * referenceRadiusToFifthPower / distanceToEighthPower * (
                 radialComponentMultiplier * k2LoveNumber * relativePosition + k2LoveNumber * timeLag * (
                    2.0 * ( relativePosition.dot( relativeVelocity ) * relativePosition / distanceSquared ) +
                    ( relativePosition.cross( planetAngularVelocityVector ) + relativeVelocity ) ) );

}


} // namespace gravitation

} // namespace tudat
