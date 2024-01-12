/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Corrected Third Printing, Springer, 2005.
 *      Bate R. Fundamentals of astro, Courier Dover Publications, 1971.
 *
 */

#include <Eigen/Core>
#include <cmath>

#include "tudat/astro/basic_astro/missionGeometry.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"

namespace tudat
{

namespace mission_geometry
{

//! Compute whether an orbit is retrograde based on inclination.
bool isOrbitRetrograde( const double inclination )
{
    bool isRetrograde = false;

    // Check which range inclination is in and return value accordingly.
    if ( inclination < 0.0 || inclination > mathematical_constants::PI )
    {
        throw std::runtime_error(
                    "The inclination is in the wrong range when determining retrogradeness" );
    }
    else if ( inclination <= mathematical_constants::PI / 2.0 )
    {
        isRetrograde = false;
    }
    else if ( inclination > mathematical_constants::PI / 2.0 )
    {
        isRetrograde = true;
    }

    return isRetrograde;
}

//! Compute whether an orbit is retrograde based on Keplerian state.
bool isOrbitRetrograde( const Eigen::Vector6d& keplerElements )
{
    // Get inclination from vector and call overloaded function.
    return isOrbitRetrograde(
                keplerElements(
                    orbital_element_conversions::inclinationIndex ) );
}

//! Compute the shadow function.
double computeShadowFunction( const Eigen::Vector3d& occultedBodyPosition,
                              const double occultedBodyRadius,
                              const Eigen::Vector3d& occultingBodyPosition,
                              const double occultingBodyRadius,
                              const Eigen::Vector3d& satellitePosition )
{
    // Calculate coordinates of the spacecraft with respect to the occulting body.
    const Eigen::Vector3d satellitePositionRelativeToOccultingBody = satellitePosition
            - occultingBodyPosition;

    // Calculate apparent radius of occulted body.
    const double occultedBodyApparentRadius
            = std::asin( occultedBodyRadius
                         / ( occultedBodyPosition - satellitePosition ).norm( ) );

    // Calculate apparent radius of occulting body.
    const double occultingBodyApparentRadius =
            std::asin( occultingBodyRadius / satellitePositionRelativeToOccultingBody.norm( ) );

    // Calculate apparent separation of the center of both bodies.
    const double apparentSeparationPartOne = -satellitePositionRelativeToOccultingBody.transpose( )
            * ( occultedBodyPosition - satellitePosition );
    const double apparentSeparationPartTwo = satellitePositionRelativeToOccultingBody.norm( )
            * ( occultedBodyPosition - satellitePosition ).norm( );
    const double apparentSeparation = std::acos( apparentSeparationPartOne
                                                 / apparentSeparationPartTwo );

    // Set initial value for the shadow function.
    double shadowFunction = 1.0;

    // Check which part of shadow satellite is in
    // Use decision schema from Fig. 5 of (Zhang et al., 2019)
    if ( std::fabs( occultedBodyApparentRadius - occultingBodyApparentRadius ) < apparentSeparation
         && apparentSeparation < occultedBodyApparentRadius + occultingBodyApparentRadius )
    {
        // Satellite is in penumbra -> partial occultation.
        // Pre-compute values for optimal computations.
        const double apparentSeparationSquared = apparentSeparation * apparentSeparation;
        const double occultedBodyApparentRadiusSquared = occultedBodyApparentRadius
                * occultedBodyApparentRadius;
        const double occultingBodyApparentRadiusSquared = occultingBodyApparentRadius
                * occultingBodyApparentRadius;

        // Partial occultation takes place, calculate the occulted area.
        const double occultedAreaPartOne
                = ( apparentSeparationSquared + occultedBodyApparentRadiusSquared
                    - occultingBodyApparentRadiusSquared ) / ( 2.0 * apparentSeparation );
        const double occultedAreaPartTwo = std::sqrt( occultedBodyApparentRadiusSquared
                                                      - occultedAreaPartOne * occultedAreaPartOne );
        const double occultedArea = occultedBodyApparentRadiusSquared
                * std::acos( occultedAreaPartOne / occultedBodyApparentRadius )
                + occultingBodyApparentRadiusSquared
                * std::acos( ( apparentSeparation - occultedAreaPartOne )
                             / occultingBodyApparentRadius )
                - apparentSeparation * occultedAreaPartTwo;
        shadowFunction = 1.0 - occultedArea / ( mathematical_constants::PI *
                                                occultedBodyApparentRadiusSquared );
    }

    else if ( apparentSeparation < occultingBodyApparentRadius - occultedBodyApparentRadius &&
         occultedBodyApparentRadius < occultingBodyApparentRadius )
    {
        // Satellite is in umbra -> total occultation.
        // Occulted circular disk is inside occulting circular disk.
        shadowFunction = 0.0;
    }

    else if ( apparentSeparation < occultedBodyApparentRadius - occultingBodyApparentRadius &&
              occultedBodyApparentRadius > occultingBodyApparentRadius )
    {
        //  Satellite is in antumbra (eclipse is annular) -> partial occultation.
        // Occulting circular disk is inside occulted circular disk.
        const double occultedBodyApparentRadiusSquared = occultedBodyApparentRadius
                                                         * occultedBodyApparentRadius;
        const double occultingBodyApparentRadiusSquared = occultingBodyApparentRadius
                                                          * occultingBodyApparentRadius;
        shadowFunction = 1.0 - occultingBodyApparentRadiusSquared / occultedBodyApparentRadiusSquared;
    }

    else if ( occultedBodyApparentRadius + occultingBodyApparentRadius <= apparentSeparation )
    {
        // No occultation.
        shadowFunction = 1.0;
    }

    // Return the shadow function
    return shadowFunction;
}

double computeSphereOfInfluence( const double distanceToCentralBody,
                                 const double ratioOfOrbitingToCentralBodyMass )
{
    // Return the radius of the sphere of influence.
    return distanceToCentralBody * std::pow( ratioOfOrbitingToCentralBodyMass, 0.4 );
}

//! Compute the sphere of influence.
double computeSphereOfInfluence( const double distanceToCentralBody,
                                 const double massOrbitingBody,
                                 const double massCentralBody )
{
    // Return the radius of the sphere of influence.
    return computeSphereOfInfluence( distanceToCentralBody, massOrbitingBody / massCentralBody );
}

} // namespace mission_geometry

} // namespace tudat
