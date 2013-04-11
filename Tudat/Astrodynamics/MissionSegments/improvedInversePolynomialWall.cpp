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
 *      130301    E.D. Brandon      Creation of code.
 *      130325    D.Dirkx           Split code into separate files for different functions.
 *
 *    References
 *      Novak, D.M. and M. Vasile. Improved Shaping Approach to the Preliminary Design of Low-
 *          Thrust Trajectories, Journal of Guidance, Control, and Dynamics 34(1), pp. 128–147,
 *          2011.
 *      Wall, B.J. and D. Novak. A 3D Shape-Based Approximation Method for Low-Thrust Trajectory
 *          Design, Advances in the Astronautical Sciences 142(2), pp. 1163-1176, 2012.
 *      Wall, B.J., Pols, B. and B. Lanktree. Shape-Based Approximation Method for Low-Thrust
 *          Interception and Rendezvous Trajectory Design, Advances in the Astronautical Sciences
 *          136(2), pp. 1447-1458, 2010.
 *
 *    Notes
 *      This file contains the improved inverse polynomial function, which is described by Wall et
 *      al. [2010] (See also Wall and Novak [2012]). This function is a mathematical representation
 *      of the radial position, r (in spherical coordinates), of a thrusting spacecraft. The
 *      improved inverse polynomial function can be used to approximate a continuous-thrust
 *      trajectory (also referred to as a low-thrust trajectory), when combined with a function
 *      which represents the (out-of-plane) elevation angle, phi (in spherical coordinates), of a
 *      thrusting spacecraft. The spherical coordinate system that is used for the calculations and
 *      the descriptions is taken from Novak and Vasile [2011].
 *
 */

#include <TudatCore/Basics/utilityMacros.h>

#include "Tudat/Astrodynamics/MissionSegments/improvedInversePolynomialWall.h"

namespace tudat
{
namespace mission_segments
{

//! Evaluate the function value for a given (in-plane) azimuthal angle.
double ImprovedInversePolynomialWall::evaluate( const double anAzimuthalAngle )
{
    // Inverse polynomial function (improved version from Wall et al. [2010]).
    const double radialDistance = 1.0 /
            ( boundaryParameters_(  ).first( 0 )
              + boundaryParameters_(  ).first( 1 )
              * std::cos( anAzimuthalAngle + boundaryParameters_(  ).first( 2 ) )
              + timeDependentParameter_(  )
                * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle
              + boundaryParameters_(  ).second( 0 ) * std::pow( anAzimuthalAngle , 4 )
              + boundaryParameters_(  ).second( 1 ) * std::pow( anAzimuthalAngle , 5 )
              + boundaryParameters_(  ).second( 2 ) * std::pow( anAzimuthalAngle , 6 ) );

    return radialDistance;
}

//! Compute the derivative of the function.
double ImprovedInversePolynomialWall::computeDerivative(
        const unsigned int order,
        const double anAzimuthalAngle )
{
    // Function itself.
    if ( order == 0 )
    {
        return evaluate( anAzimuthalAngle );
    }
    // First derivative.
    else if ( order == 1 )
    {
        // Radial distance.
        double radialDistance = evaluate( anAzimuthalAngle );

        // First derivative of the inverse polynomial function with respect to the (in-plane)
        // azimuthal angle.
        const double firstDerivative = -1.0 * radialDistance * radialDistance
                * ( -1.0 * boundaryParameters_(  ).first( 1 )
                    * std::sin( anAzimuthalAngle + boundaryParameters_(  ).first( 2 ) )
                    + 3.0 * timeDependentParameter_(  ) * anAzimuthalAngle * anAzimuthalAngle
                    + 4.0 * boundaryParameters_(  ).second( 0 )
                        * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle
                    + 5.0 * boundaryParameters_(  ).second( 1 ) * std::pow( anAzimuthalAngle , 4 )
                    + 6.0 * boundaryParameters_(  ).second( 2 )
                    * std::pow( anAzimuthalAngle , 5 ) );

        return firstDerivative;
    }
    // Second derivative.
    else if ( order == 2 )
    {
        // Radial distance.
        double radialDistance = evaluate( anAzimuthalAngle );

        // First derivative.
        double firstDerivative = computeDerivative( ( order - 1 ) , anAzimuthalAngle );

        // Second derivative of the inverse polynomial function with respect to the (in-plane)
        // azimuthal angle.
        const double secondDerivative = ( 2.0 / radialDistance )
                * firstDerivative * firstDerivative
                - radialDistance * radialDistance
                * ( -1.0 * boundaryParameters_(  ).first( 1 )
                    * std::cos( anAzimuthalAngle + boundaryParameters_(  ).first( 2 ) )
                    + 6.0 * timeDependentParameter_(  ) * anAzimuthalAngle
                    + 12.0 * boundaryParameters_(  ).second( 0 )
                        * anAzimuthalAngle * anAzimuthalAngle
                    + 20.0 * boundaryParameters_(  ).second( 1 )
                        * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle
                    + 30.0 * boundaryParameters_(  ).second( 2 ) * std::pow( anAzimuthalAngle , 4 ) );

        return secondDerivative;
    }
    // Throw runtime error, when order is higher that 2.
    else
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Derivatives of order higher than 2 are not supported.") ) );
    }
}

//! Compute the definite integral of the function.
double ImprovedInversePolynomialWall::computeDefiniteIntegral(
        const unsigned int order, const double lowerBound, const double upperBound )
{
    // Declare unused parameters.
    TUDAT_UNUSED_PARAMETER( order );
    TUDAT_UNUSED_PARAMETER( lowerBound );
    TUDAT_UNUSED_PARAMETER( upperBound );

    // Throw runtime error, when trying to compute the integral of the function.
    boost::throw_exception(
                boost::enable_error_info(
                    std::runtime_error(
                        "Cannot compute the integral, "
                        "this is not supported for this function class.") ) );
}

} // namespace mission_segments
} // namespace tudat
