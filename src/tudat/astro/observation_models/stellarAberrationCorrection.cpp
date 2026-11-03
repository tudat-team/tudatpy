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

namespace
{

//! Check the input of the stellar-aberration functions and return the observer velocity divided by the speed of light.
Eigen::Vector3d checkStellarAberrationInputAndGetNormalizedVelocity( const Eigen::Vector3d& astrometricDirection,
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

    const Eigen::Vector3d observerVelocityDividedBySpeedOfLight = observerVelocity / physical_constants::SPEED_OF_LIGHT;
    if( observerVelocityDividedBySpeedOfLight.squaredNorm( ) >= 1.0 )
    {
        throw std::runtime_error( "Cannot apply stellar aberration for an observer moving at or above the speed of light." );
    }
    return observerVelocityDividedBySpeedOfLight;
}

}  // namespace

Eigen::Vector3d calculateApparentDirectionWithStellarAberration( const Eigen::Vector3d& astrometricDirection,
                                                                 const Eigen::Vector3d& observerVelocity )
{
    const Eigen::Vector3d observerVelocityDividedBySpeedOfLight =
            checkStellarAberrationInputAndGetNormalizedVelocity( astrometricDirection, observerVelocity );
    const Eigen::Vector3d astrometricUnitDirection = astrometricDirection.normalized( );
    const double observerSpeedSquaredDividedBySpeedOfLightSquared = observerVelocityDividedBySpeedOfLight.squaredNorm( );

    const double directionDotObserverVelocity = astrometricUnitDirection.dot( observerVelocityDividedBySpeedOfLight );
    const double astrometricDirectionScaleFactor = -directionDotObserverVelocity +
            std::sqrt( 1.0 - observerSpeedSquaredDividedBySpeedOfLightSquared +
                       directionDotObserverVelocity * directionDotObserverVelocity );
    return ( astrometricDirectionScaleFactor * astrometricUnitDirection + observerVelocityDividedBySpeedOfLight ).normalized( );
}

void calculateStellarAberrationPartials( const Eigen::Vector3d& astrometricDirection,
                                         const Eigen::Vector3d& observerVelocity,
                                         Eigen::Matrix3d& partialWrtAstrometricDirection,
                                         Eigen::Matrix3d& partialWrtObserverVelocity )
{
    const Eigen::Vector3d observerVelocityDividedBySpeedOfLight =
            checkStellarAberrationInputAndGetNormalizedVelocity( astrometricDirection, observerVelocity );
    const double astrometricDirectionNorm = astrometricDirection.norm( );
    const Eigen::Vector3d astrometricUnitDirection = astrometricDirection / astrometricDirectionNorm;

    // Unnormalized apparent direction w = s * u + beta, with s = -u.beta + sqrt( 1 - beta.beta + ( u.beta )^2 )
    const double directionDotObserverVelocity = astrometricUnitDirection.dot( observerVelocityDividedBySpeedOfLight );
    const double squareRootTerm = std::sqrt( 1.0 - observerVelocityDividedBySpeedOfLight.squaredNorm( ) +
                                             directionDotObserverVelocity * directionDotObserverVelocity );
    const double astrometricDirectionScaleFactor = -directionDotObserverVelocity + squareRootTerm;
    const Eigen::Vector3d unnormalizedApparentDirection =
            astrometricDirectionScaleFactor * astrometricUnitDirection + observerVelocityDividedBySpeedOfLight;
    const double unnormalizedApparentDirectionNorm = unnormalizedApparentDirection.norm( );
    const Eigen::Vector3d apparentDirection = unnormalizedApparentDirection / unnormalizedApparentDirectionNorm;

    // Partial of normalized w.r.t. unnormalized apparent direction
    const Eigen::Matrix3d apparentDirectionNormalizationPartial =
            ( Eigen::Matrix3d::Identity( ) - apparentDirection * apparentDirection.transpose( ) ) / unnormalizedApparentDirectionNorm;

    // Partials of s w.r.t. u and beta
    const double scaleFactorPartialWrtDirectionDotVelocity = -1.0 + directionDotObserverVelocity / squareRootTerm;
    const Eigen::Vector3d scaleFactorPartialWrtNormalizedVelocity =
            scaleFactorPartialWrtDirectionDotVelocity * astrometricUnitDirection - observerVelocityDividedBySpeedOfLight / squareRootTerm;

    // Partial w.r.t. astrometric direction, including the normalization of the input vector
    const Eigen::Matrix3d unnormalizedApparentDirectionPartialWrtUnitDirection =
            astrometricDirectionScaleFactor * Eigen::Matrix3d::Identity( ) +
            scaleFactorPartialWrtDirectionDotVelocity * astrometricUnitDirection * observerVelocityDividedBySpeedOfLight.transpose( );
    const Eigen::Matrix3d unitDirectionPartialWrtDirection =
            ( Eigen::Matrix3d::Identity( ) - astrometricUnitDirection * astrometricUnitDirection.transpose( ) ) / astrometricDirectionNorm;
    partialWrtAstrometricDirection =
            apparentDirectionNormalizationPartial * unnormalizedApparentDirectionPartialWrtUnitDirection * unitDirectionPartialWrtDirection;

    // Partial w.r.t. observer velocity
    const Eigen::Matrix3d unnormalizedApparentDirectionPartialWrtNormalizedVelocity =
            Eigen::Matrix3d::Identity( ) + astrometricUnitDirection * scaleFactorPartialWrtNormalizedVelocity.transpose( );
    partialWrtObserverVelocity = apparentDirectionNormalizationPartial * unnormalizedApparentDirectionPartialWrtNormalizedVelocity /
            physical_constants::SPEED_OF_LIGHT;
}

}  // namespace observation_models
}  // namespace tudat
