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
 *      Novak, D.M. Methods and Tools for Preliminary Low Thrust Mission Analysis. PhD Thesis,
 *          University of Glasgow, UK, 2012.
 *      Novak, D.M. and M. Vasile. Improved Shaping Approach to the Preliminary Design of Low-
 *          Thrust Trajectories, Journal of Guidance, Control, and Dynamics 34(1), pp. 128–147,
 *          2011.
 *      Wall, B.J. and D. Novak. A 3D Shape-Based Approximation Method for Low-Thrust Trajectory
 *          Design, Advances in the Astronautical Sciences 142(2), pp. 1163-1176, 2012.
 *
 *    Notes
 *      This file contains a specific oscillating function, which is described by Novak and Vasile
 *      [2011] (See also Novak [2012] and Wall and Novak [2012]). This function is a mathematical
 *      representation of the (out-of-plane) elevation angle, phi (in spherical coordinates), of a
 *      thrusting spacecraft. This oscillating function can be used to approximate a
 *      continuous-thrust trajectory (also referred to as a low-thrust trajectory), when combined
 *      with a function which represents the radial position, r (in spherical coordinates), of a
 *      thrusting spacecraft. The spherical coordinate system that is used for the calculations and
 *      the descriptions is taken from Novak and Vasile [2011].
 */

#include <TudatCore/Basics/utilityMacros.h>

#include "Tudat/Astrodynamics/MissionSegments/oscillatingFunctionNovak.h"

namespace tudat
{
namespace mission_segments
{

//! Evaluate the function value for a given (in-plane) azimuthal angle.
double OscillatingFunctionNovak::evaluate( const double anAzimuthalAngle )
{
    // Oscillating function [Novak and Vasile, 2011].
    const double elevationAngle =
            ( boundaryParameters_(  ).first( 0 )
              + boundaryParameters_(  ).first( 1 ) * anAzimuthalAngle )
            * std::cos( anAzimuthalAngle )
            + ( boundaryParameters_(  ).second( 0 )
                + boundaryParameters_(  ).second( 1 ) * anAzimuthalAngle )
            * std::sin ( anAzimuthalAngle );

    return elevationAngle;
}

//! Compute the derivative of the function.
double OscillatingFunctionNovak::computeDerivative(
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
        // First derivative of the oscillating function with respect to the (in-plane) azimuthal angle.
        const double firstDerivative = (
                    boundaryParameters_(  ).first( 1 )
                    + boundaryParameters_(  ).second( 0 )
                    + boundaryParameters_(  ).second( 1 ) * anAzimuthalAngle )
                * std::cos( anAzimuthalAngle )
                - ( boundaryParameters_(  ).first( 0 )
                    + boundaryParameters_(  ).first( 1 ) * anAzimuthalAngle
                    - boundaryParameters_(  ).second( 1 ) ) * std::sin( anAzimuthalAngle );
        return firstDerivative;
    }
    // Second derivative.
    else if ( order == 2 )
    {
        // Second derivative of the oscillating function with respect to the (in-plane) azimuthal
        // angle.
        const double secondDerivative = (
                    - 1.0 * boundaryParameters_(  ).first( 0 )
                    - boundaryParameters_(  ).first( 1 ) * anAzimuthalAngle
                    + 2.0 * boundaryParameters_(  ).second( 1 ) )
                * std::cos( anAzimuthalAngle ) - (
                    2.0 * boundaryParameters_(  ).first( 1 )
                    + boundaryParameters_(  ).second( 0 )
                    + boundaryParameters_(  ).second( 1 ) * anAzimuthalAngle )
                * std::sin( anAzimuthalAngle );

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
double OscillatingFunctionNovak::computeDefiniteIntegral(
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
