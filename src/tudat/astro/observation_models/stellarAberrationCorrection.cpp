/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/stellarAberrationCorrection.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{
namespace observation_models
{

Eigen::Vector3d calculateApparentDirectionWithStellarAberration( const Eigen::Vector3d& astrometricDirection,
                                                                 const Eigen::Vector3d& observerVelocity )
{
    if( !astrometricDirection.allFinite( ) || astrometricDirection.norm( ) <= 100.0 * std::numeric_limits< double >::epsilon( ) )
    {
        throw std::runtime_error( "Cannot apply stellar aberration to a non-finite or zero line-of-sight vector." );
    }
    if( !observerVelocity.allFinite( ) )
    {
        throw std::runtime_error( "Cannot apply stellar aberration with a non-finite observer velocity." );
    }

    const Eigen::Vector3d astrometricUnitDirection = astrometricDirection.normalized( );
    const Eigen::Vector3d observerVelocityDividedBySpeedOfLight = observerVelocity / physical_constants::SPEED_OF_LIGHT;
    const double observerSpeedSquaredDividedBySpeedOfLightSquared = observerVelocityDividedBySpeedOfLight.squaredNorm( );
    if( observerSpeedSquaredDividedBySpeedOfLightSquared >= 1.0 )
    {
        throw std::runtime_error( "Cannot apply stellar aberration for an observer moving at or above the speed of light." );
    }

    const double directionDotObserverVelocity = astrometricUnitDirection.dot( observerVelocityDividedBySpeedOfLight );
    const double astrometricDirectionScaleFactor = -directionDotObserverVelocity +
            std::sqrt( 1.0 - observerSpeedSquaredDividedBySpeedOfLightSquared +
                       directionDotObserverVelocity * directionDotObserverVelocity );
    return ( astrometricDirectionScaleFactor * astrometricUnitDirection + observerVelocityDividedBySpeedOfLight ).normalized( );
}

}  // namespace observation_models
}  // namespace tudat
