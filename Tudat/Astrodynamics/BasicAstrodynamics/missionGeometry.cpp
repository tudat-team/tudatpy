/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      120530    M.I. Ganeff       Creation of code.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Corrected Third Printing, Springer, 2005.
 *
 *    Notes
 *
 */

#include <cmath>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace mission_geometry
{

//! Compute the shadow function.
double computeShadowFunction( const Eigen::VectorXd occultedBodyPosition,
                              const double occultedBodyRadius,
                              const Eigen::VectorXd occultingBodyPosition,
                              const double occultingBodyRadius,
                              const Eigen::VectorXd satelliteState )
{
    // Calculate coordinates of the spacecraft with respect to the occulting
    // body.
    const Eigen::Vector3d satelliteToOcultingBody = satelliteState.segment( 0, 3 )
            - occultingBodyPosition;

    // Calculate apparent radius of occulted body.
    const double occultedBodyApparentRadius
            = std::asin( occultedBodyRadius / ( occultedBodyPosition -
                                           satelliteState.segment( 0, 3 ) ).norm( ) );

    // Calculate apparent radius of occulting body.
    const double occultingBodyApparentRadius =
            std::asin( occultingBodyRadius / satelliteToOcultingBody.norm( ) );

    // Calculate apparent speration of the center of both bodies.
    const double apparentSeparationPartOne = - satelliteToOcultingBody.transpose( ) *
            ( occultedBodyPosition - satelliteState.segment( 0, 3 ) );
    const double apparentSeparationPartTwo = satelliteToOcultingBody.norm( ) *
            ( occultedBodyPosition - satelliteState.segment( 0, 3 ) ).norm( );
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
        shadowFunction = 1.0 - occultedArea / ( mathematics::PI *
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

} // namespace mission_geometry
} // namespace tudat
