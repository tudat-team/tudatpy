/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Novak, D.M. Methods and Tools for Preliminary Low Thrust Mission Analysis. PhD Thesis,
 *          University of Glasgow, UK, 2012.
 *      Novak, D.M. and M. Vasile. Improved Shaping Approach to the Preliminary Design of Low-
 *          Thrust Trajectories, Journal of Guidance, Control, and Dynamics 34(1), pp. 128147,
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

#include "Tudat/Basics/utilityMacros.h"

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
    // Third derivative.
    else if ( order == 3 )
    {
        // Third derivative of the oscillating function with repect to the (in-plane) azimuthal
        // angle.
        const double thirdDerivative = -1.0 * (
                    3.0 * boundaryParameters_(  ).first( 1 )
                    + boundaryParameters_(  ).second( 0 )
                    + boundaryParameters_(  ).second( 1 ) * anAzimuthalAngle )
                * std::cos( anAzimuthalAngle ) + (
                    boundaryParameters_(  ).first( 0 )
                    + boundaryParameters_(  ).first( 1 ) * anAzimuthalAngle
                    - 3.0 * boundaryParameters_(  ).second( 1 ) )
                * std::sin( anAzimuthalAngle );

        return thirdDerivative;
    }
    // Throw runtime error, when order is higher that 3.
    else
    {
        throw std::runtime_error( "Derivatives of order higher than 3 are not supported." );
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
    throw std::runtime_error( "Cannot compute the integral, this is not supported for this function class." );
}

} // namespace mission_segments
} // namespace tudat
