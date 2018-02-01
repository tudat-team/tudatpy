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
 *      Novak, D.M. and M. Vasile. Improved Shaping Approach to the Preliminary Design of Low-
 *          Thrust Trajectories, Journal of Guidance, Control, and Dynamics 34(1), pp. 128147,
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

#include "Tudat/Basics/utilityMacros.h"

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
              + boundaryParameters_(  ).second( 0 ) * anAzimuthalAngle * anAzimuthalAngle
                * anAzimuthalAngle * anAzimuthalAngle
              + boundaryParameters_(  ).second( 1 ) * anAzimuthalAngle * anAzimuthalAngle
                * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle
              + boundaryParameters_(  ).second( 2 ) * anAzimuthalAngle * anAzimuthalAngle
                * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle );

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
                    + 5.0 * boundaryParameters_(  ).second( 1 )
                        * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle
                    + 6.0 * boundaryParameters_(  ).second( 2 )
                        * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle
                        * anAzimuthalAngle );

        return firstDerivative;
    }
    // Second derivative.
    else if ( order == 2 )
    {
        // Radial distance.
        double radialDistance = evaluate( anAzimuthalAngle );

        // First derivative.
        double firstDerivative = computeDerivative( 1 , anAzimuthalAngle );

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
                    + 30.0 * boundaryParameters_(  ).second( 2 )
                        * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle
                        * anAzimuthalAngle );

        return secondDerivative;
    }
    // Third derivative.
    else if ( order == 3 )
    {
        // Radial distance.
        double radialDistance = evaluate( anAzimuthalAngle );

        // First derivative.
        double firstDerivative = computeDerivative( 1 , anAzimuthalAngle );

        // Second derivative.
        double secondDerivative = computeDerivative( 2 , anAzimuthalAngle );

        // Third derivative of the inverse polynomial function with respect to the (in-plane)
        // azimuthal angle.
        const double thirdDerivative = 6.0 *
                ( ( firstDerivative * firstDerivative * firstDerivative ) /
                  ( radialDistance * radialDistance ) )
                - 6.0 * ( firstDerivative / radialDistance ) * (
                    ( ( 2.0 * firstDerivative * firstDerivative ) / radialDistance )
                    - secondDerivative )
                - radialDistance * radialDistance * (
                    boundaryParameters_(  ).first( 1 )
                    * std::sin( anAzimuthalAngle + boundaryParameters_(  ).first( 2 ) )
                    + 6.0 * timeDependentParameter_(  )
                    + 24.0 * boundaryParameters_(  ).second( 0 ) * anAzimuthalAngle
                    + 60.0 * boundaryParameters_(  ).second( 1 )
                        * anAzimuthalAngle * anAzimuthalAngle
                    + 120.0 * boundaryParameters_(  ).second( 2 )
                        * anAzimuthalAngle * anAzimuthalAngle * anAzimuthalAngle );

        return thirdDerivative;
    }
    // Throw runtime error, when order is higher that 3.
    else
    {
        throw std::runtime_error( "Derivatives of order higher than 3 are not supported." );
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
    throw std::runtime_error( "Cannot compute the integral, this is not supported for this function class." );
}

} // namespace mission_segments
} // namespace tudat
