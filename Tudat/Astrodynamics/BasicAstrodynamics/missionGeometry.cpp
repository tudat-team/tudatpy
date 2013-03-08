/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120530    M.I. Ganeff       Code created.
 *      121004    M.I. Ganeff       Input parameter types and variable-naming updated.
 *      121018    M.I. Ganeff       Added computeSphereOfInfluence().
 *      121123    D. Dirkx          Added computeSphereOfInfluence() function taking mass ratios;
 *                                  updated implementation of computeSphereOfInfluence() taking
 *                                  masses.
 *      130225    D. Dirkx          Added isOrbitRetrograde(...) functions taking inclination and
 *                                  Kepler vector.
 *      130227    D. Dirkx          Set isRetrograde at initialization to 0.
 *      130301    R.C.A. Boon       Minor textual changes, changed mathematics::PI to
 *      130305    R.C.A. Boon       Replaced Eigen::VectorXd by tudat::basic_mathematics::Vector6d
 *                                  basic_mathematics::mathematical_constants::PI.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Corrected Third Printing, Springer, 2005.
 *      Bate R. Fundamentals of Astrodynamics, Courier Dover Publications, 1971.
 *
 *    Notes
 *
 */

#include <cmath>
#include <stdexcept>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"

namespace tudat
{
namespace basic_astrodynamics
{
namespace mission_geometry
{

//! Compute whether an orbit is retrograde based on inclination.
bool isOrbitRetrograde( const double inclination )
{
    bool isRetrograde = false;

    // Check which range inclination is in and return value accordingly.
    if ( inclination < 0.0 || inclination > mathematics::PI )
    {
        throw std::runtime_error(
                    "The inclination is in the wrong range when determining retrogradeness" );
    }
    else if ( inclination <= basic_mathematics::mathematical_constants::PI / 2.0 )
    {
        isRetrograde = false;
    }
    else if ( inclination > basic_mathematics::mathematical_constants::PI / 2.0 )
    {
        isRetrograde = true;
    }

    return isRetrograde;
}

//! Compute whether an orbit is retrograde based on Keplerian state.
bool isOrbitRetrograde( const tudat::basic_mathematics::Vector6d keplerElements )
{
    // Get inclination from vector and call overloaded function.
    return isOrbitRetrograde(
                keplerElements(
                    basic_astrodynamics::orbital_element_conversions::inclinationIndex ) );
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

    // Check if partial occultation takes place
    if ( std::fabs( occultedBodyApparentRadius - occultingBodyApparentRadius ) < apparentSeparation
         && apparentSeparation < occultedBodyApparentRadius + occultingBodyApparentRadius )
    {
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
        shadowFunction = 1.0 - occultedArea / ( basic_mathematics::mathematical_constants::PI *
                                                occultedBodyApparentRadiusSquared );
    }

    else
    {
        // Full or no occultation takes place.
        // Check for type of occultation.
        if ( apparentSeparation < occultingBodyApparentRadius - occultedBodyApparentRadius &&
             occultedBodyApparentRadius < occultingBodyApparentRadius )
        {
            // Total occultation.
            shadowFunction = 0.0;
        }

        else if ( apparentSeparation < occultedBodyApparentRadius - occultingBodyApparentRadius &&
                  occultedBodyApparentRadius > occultingBodyApparentRadius )
        {
            // Maximum partial occultation.
            shadowFunction = 0.0;
        }

        else if ( occultedBodyApparentRadius + occultingBodyApparentRadius <= apparentSeparation )
        {
            // No occultation
            shadowFunction = 1.0;
        }
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
} // namespace basic_astrodynamics
} // namespace tudat
