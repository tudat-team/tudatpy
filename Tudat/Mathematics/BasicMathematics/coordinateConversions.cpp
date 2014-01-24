/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      100903    K. Kumar          File created.
 *      100916    L. Abdulkadir     File checked.
 *      100929    K. Kumar          Checked code by D. Dirkx added.
 *      101110    K. Kumar          Added raiseToExponentPower() function.
 *      102410    D. Dirkx          Minor comment changes as code check.
 *      101213    K. Kumar          Bugfix raiseToIntegerExponent(); renamed raiseToIntegerPower().
 *                                  Added computeAbsoluteValue() functions.
 *      110202    K. Kumar          Added overload for State* for computeLinearInterpolation().
 *      110111    J. Melman         Added computeModulo() function.
 *      110411    K. Kumar          Added convertCartesianToSpherical() function.
 *      110606    J. Melman         Removed possible singularity from
 *                                  convertCartesianToSpherical.
 *      110707    K. Kumar          Added computeSampleMean(), computeSampleVariance() functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120127    D. Dirkx          First version branched from basic mathematics in Tudat Core.
 *      120127    K. Kumar          Minor comment edits.
 *      120118    D. Gondelach      Added new convertCylindricalToCartesian functions.
 *      120214    K. Kumar          Branched from old Tudat trunk for new coordinate conversions.
 *      120217    K. Kumar          Updated computeModuloForSignedValues() to computeModulo()
 *                                  from Tudat Core.
 *      120926    E. Dekens         Added spherical gradient to Cartesian conversion.
 *      131022    T. Roegiers       Added conversion from spherical state to Cartesian state.
 *                                  Added conversion from Cartesian state to spherical state.
 *      140114    E. Brandon        Reorganized includes.
 *                                  Minor changes during code check.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *      Torok, J.S. Analytical Mechanics: with an Introduction to Dynamical Systems, John Wiley and
 *          Sons, Inc., 2000.
 *      Vallado, D.A. Fundamentals of Astrodynamics and Applications. Microcosm Press, 2001.
 *
 *    Notes
 *
 */

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/math/special_functions/sign.hpp>

#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

namespace tudat
{
namespace basic_mathematics
{
namespace coordinate_conversions
{

//! Convert cylindrical to Cartesian coordinates.
Eigen::Vector3d convertCylindricalToCartesian( const double radius,
                                               const double azimuthAngle, const double z )
{
    // Create Cartesian coordinates vector.
    Eigen::Vector3d cartesianCoordinates;

    // If radius < 0, then give warning.
    if ( radius < 0.0 )
    {
        std::cerr << "Warning: cylindrical radial coordinate is negative!\n"
                  << "This could give incorrect results!\n";
    }

    // Compute and set Cartesian coordinates.
    cartesianCoordinates << radius * std::cos( azimuthAngle ),   // x-coordinate
                            radius * std::sin( azimuthAngle ),   // y-coordinate
                            z;                                   // z-coordinate

    return cartesianCoordinates;
}

//! Convert cylindrical to cartesian coordinates.
Eigen::Vector3d convertCylindricalToCartesian( const Eigen::Vector3d& cylindricalCoordinates )
{
    // Create Cartesian coordinates vector.
    Eigen::Vector3d cartesianCoordinates;

    // If radius < 0, then give warning.
    if ( cylindricalCoordinates( 0 ) < 0.0 )
    {
        std::cerr << "Warning: cylindrical radial coordinate is negative! "
                  << "This could give incorrect results!\n";
    }

    // Compute and set Cartesian coordinates.
    cartesianCoordinates
            << cylindricalCoordinates( 0 )
               * std::cos( cylindricalCoordinates( 1 ) ),    // x-coordinate
               cylindricalCoordinates( 0 )
               * std::sin( cylindricalCoordinates( 1 ) ),    // y-coordinate
               cylindricalCoordinates( 2 );                  // z-coordinate

    return cartesianCoordinates;
}

//! Convert cylindrical to Cartesian state.
Eigen::VectorXd convertCylindricalToCartesian( const Eigen::VectorXd& cylindricalState )
{
    // Create Cartesian state vector, initialized with zero entries.
    Eigen::VectorXd cartesianState = Eigen::VectorXd::Zero( 6 );

    // Get azimuth angle, theta.
    double azimuthAngle = cylindricalState( 1 );

    // Compute and set Cartesian coordinates.
    cartesianState.head( 3 ) = convertCylindricalToCartesian(
                Eigen::Vector3d( cylindricalState.head( 3 ) ) );

    // If r = 0 AND Vtheta > 0, then give warning and assume Vtheta=0.
    if ( std::fabs(cylindricalState( 0 )) <= std::numeric_limits< double >::epsilon( )
         && std::fabs(cylindricalState( 4 )) > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "Warning: cylindrical velocity Vtheta (r*thetadot) does not equal zero while "
                  << "the radius (r) is zero! Vtheta is taken equal to zero!\n";

        // Compute and set Cartesian velocities.
        cartesianState.tail( 3 )
                << cylindricalState( 3 ) * std::cos( azimuthAngle ),   // xdot
                   cylindricalState( 3 ) * std::sin( azimuthAngle ),   // ydot
                   cylindricalState( 5 );                              // zdot
    }

    else
    {
        // Compute and set Cartesian velocities.
        cartesianState.tail( 3 )
                << cylindricalState( 3 ) * std::cos( azimuthAngle )
                   - cylindricalState( 4 ) * std::sin( azimuthAngle ),   // xdot
                   cylindricalState( 3 ) * std::sin( azimuthAngle )
                   + cylindricalState( 4 ) * std::cos( azimuthAngle ),   // ydot
                   cylindricalState( 5 );                                // zdot
    }

    return cartesianState;
}

//! Convert Cartesian to cylindrical coordinates.
Eigen::Vector3d convertCartesianToCylindrical( const Eigen::Vector3d& cartesianCoordinates )
{
    // Create cylindrical coordinates vector.
    Eigen::Vector3d cylindricalCoordinates;

    // Declare new variable, the azimuth angle.
    double azimuthAngle;

    // Compute azimuth angle, theta.
    /* If x = 0, then azimuthAngle = pi/2 (y>0) or 3*pi/2 (y<0) or 0 (y=0),
       else azimuthAngle = arctan(y/x).
    */
    using tudat::basic_mathematics::mathematical_constants::PI;
    if ( std::fabs(cartesianCoordinates( 0 ) ) <= std::numeric_limits< double >::epsilon( ) )
    {
        azimuthAngle = tudat::basic_mathematics::computeModulo(
                    static_cast< double >( boost::math::sign( cartesianCoordinates( 1 ) ) )
                    * 0.5 * PI, 2.0 * PI );
    }

    else
    {
        azimuthAngle = tudat::basic_mathematics::computeModulo(
                    std::atan2( cartesianCoordinates( 1 ),
                                cartesianCoordinates( 0 ) ), 2.0 * PI );
    }

    // Compute and set cylindrical coordinates.
    cylindricalCoordinates <<
        std::sqrt( pow( cartesianCoordinates( 0 ), 2 )
                   + pow( cartesianCoordinates( 1 ), 2 ) ), // Radius
        azimuthAngle,                                       // Azimuth angle, theta
        cartesianCoordinates( 2 );                          // z-coordinate

    return cylindricalCoordinates;
}

//! Convert Cartesian to cylindrical state.
Eigen::VectorXd convertCartesianToCylindrical( const Eigen::VectorXd& cartesianState )
{
    // Create cylindrical state vector, initialized with zero entries.
    Eigen::VectorXd cylindricalState = Eigen::VectorXd::Zero( 6 );

    // Compute and set cylindrical coordinates.
    cylindricalState.head( 3 ) = convertCartesianToCylindrical(
                Eigen::Vector3d( cartesianState.head( 3 ) ) );

    // Compute and set cylindrical velocities.
    /* If radius = 0, then Vr = sqrt(xdot^2+ydot^2) and Vtheta = 0,
       else Vr = (x*xdot+y*ydot)/radius and Vtheta = (x*ydot-y*xdot)/radius.
    */
    if ( cylindricalState( 0 ) <= std::numeric_limits< double >::epsilon( ) )
    {
        cylindricalState.tail( 3 ) <<
            std::sqrt( pow( cartesianState( 3 ), 2 ) + pow( cartesianState( 4 ), 2 ) ), // Vr
            0.0,                                                                        // Vtheta
            cartesianState( 5 );                                                        // Vz
    }

    else
    {
        cylindricalState.tail( 3 ) <<
            ( cartesianState( 0 ) * cartesianState( 3 )
              + cartesianState( 1 ) * cartesianState( 4 ) ) / cylindricalState( 0 ),    // Vr
            ( cartesianState( 0 ) * cartesianState( 4 )
              - cartesianState( 1 ) * cartesianState( 3 ) ) / cylindricalState( 0 ),    // Vtheta
                cartesianState( 5 );                                                    // Vz
    }

    return cylindricalState;
}

//! Convert spherical to Cartesian gradient.
Eigen::Vector3d convertSphericalToCartesianGradient( const Eigen::Vector3d& sphericalGradient,
                                                     const Eigen::Vector3d& cartesianCoordinates )
{
    // Compute radius.
    const double radius = std::sqrt( cartesianCoordinates( 0 ) * cartesianCoordinates( 0 )
                                     + cartesianCoordinates( 1 ) * cartesianCoordinates( 1 )
                                     + cartesianCoordinates( 2 ) * cartesianCoordinates( 2 ) );

    // Compute square of distance within xy-plane.
    const double xyDistanceSquared = cartesianCoordinates( 0 ) * cartesianCoordinates( 0 )
            + cartesianCoordinates( 1 ) * cartesianCoordinates( 1 );

    // Compute distance within xy-plane.
    const double xyDistance = std::sqrt( xyDistanceSquared );

    // Compute transformation matrix.
    const Eigen::Matrix3d transformationMatrix = (
                Eigen::Matrix3d( 3, 3 ) <<
                cartesianCoordinates( 0 ) / radius,
                - cartesianCoordinates( 0 ) * cartesianCoordinates( 2 )
                / ( radius * radius * xyDistance ),
                - cartesianCoordinates( 1 ) / xyDistanceSquared,
                cartesianCoordinates( 1 ) / radius,
                - cartesianCoordinates( 1 ) * cartesianCoordinates( 2 )
                / ( radius * radius * xyDistance ),
                + cartesianCoordinates( 0 ) / xyDistanceSquared,
                cartesianCoordinates( 2 ) / radius,
                xyDistance / ( radius * radius ),
                0.0
                ).finished( );

    // Return Cartesian gradient.
    return transformationMatrix * sphericalGradient;
}

//! Convert spherical to Cartesian state.
Eigen::VectorXd convertSphericalToCartesianState( const Eigen::VectorXd& sphericalState )
{
    // Create Cartesian state vector, initialized with zero entries.
    Eigen::VectorXd convertedCartesianState = Eigen::VectorXd::Zero( 6 );

    // Create local variables.
    const double radius = sphericalState( 0 );
    const double azimuthAngle = sphericalState( 1 );
    const double elevationAngle = sphericalState( 2 );

    // Precompute sine/cosine of angles, which has multiple usages, to save computation time.
    const double cosineOfElevationAngle = std::cos( elevationAngle );
    const double sineOfElevationAngle = std::sin( elevationAngle );
    const double cosineOfAzimuthAngle = std::cos( azimuthAngle );
    const double sineOfAzimuthAngle = std::sin( azimuthAngle );

    // Set up transformation matrix for spherical to cylindrical conversion.
    Eigen::MatrixXd transformationMatrixSphericalToCylindrical = Eigen::MatrixXd::Zero( 3, 3 );
    transformationMatrixSphericalToCylindrical( 0, 0 ) = cosineOfElevationAngle;
    transformationMatrixSphericalToCylindrical( 0, 2 ) = -sineOfElevationAngle;
    transformationMatrixSphericalToCylindrical( 1, 1 ) = 1.0;
    transformationMatrixSphericalToCylindrical( 2, 0 ) = sineOfElevationAngle;
    transformationMatrixSphericalToCylindrical( 2, 2 ) = cosineOfElevationAngle;

    // Set up transformation matrix for cylindrical to Cartesian conversion.
    Eigen::MatrixXd transformationMatrixCylindricalToCartesian = Eigen::MatrixXd::Zero( 3, 3 );
    transformationMatrixCylindricalToCartesian( 0, 0 ) = cosineOfAzimuthAngle;
    transformationMatrixCylindricalToCartesian( 0, 1 ) = -sineOfAzimuthAngle;
    transformationMatrixCylindricalToCartesian( 1, 0 ) = sineOfAzimuthAngle;
    transformationMatrixCylindricalToCartesian( 1, 1 ) = cosineOfAzimuthAngle;
    transformationMatrixCylindricalToCartesian( 2, 2 ) = 1.0;

    // Compute transformation matrix for spherical to Cartesian conversion.
    const Eigen::MatrixXd transformationMatrixSphericalToCartesian
            = transformationMatrixCylindricalToCartesian
            * transformationMatrixSphericalToCylindrical;

    // Perform transformation of position coordinates.
    convertedCartesianState( 0 ) = radius * cosineOfAzimuthAngle * cosineOfElevationAngle;
    convertedCartesianState( 1 ) = radius * sineOfAzimuthAngle * cosineOfElevationAngle;
    convertedCartesianState( 2 ) = radius * sineOfElevationAngle;

    // Perform transformation of velocity vector.
    convertedCartesianState.segment( 3, 3 ) =
        transformationMatrixSphericalToCartesian * sphericalState.segment( 3, 3 );

    // Return Cartesian state vector.
    return convertedCartesianState;
}

//! Convert Cartesian to spherical state.
Eigen::VectorXd convertCartesianToSphericalState( const Eigen::VectorXd& cartesianState )
{
    // Create spherical state vector, initialized with zero entries.
    Eigen::VectorXd convertedSphericalState = Eigen::VectorXd::Zero( 6 );

    // Compute radius.
    convertedSphericalState( 0 ) = cartesianState.segment( 0, 3 ).norm( );

    // Check if radius is nonzero.
    /*
     * If r > 0, the elevation and azimuth angles are computed using trigonometric relationships.
     * If r = 0, the coordinates are at the origin, the elevation and azimuth angles equal to zero.
     * Since the state vector was initialized with zeroes, this is already the case.
     */
    if ( convertedSphericalState( 0 ) > std::numeric_limits< double >::epsilon( ) )
    {
        // Compute elevation and azimuth angles using trigonometric relationships.
        // Azimuth angle.
        convertedSphericalState( 1 ) = std::atan2( cartesianState( 1 ), cartesianState( 0 ) );
        // Elevation angle.
        convertedSphericalState( 2 ) = std::asin( cartesianState( 2 )
                                                   / convertedSphericalState( 0 ) );
    }

    // Precompute sine/cosine of angles, which has multiple usages, to save computation time.
    const double cosineOfElevationAngle = std::cos( convertedSphericalState( 2 ) );
    const double sineOfElevationAngle = std::sin( convertedSphericalState( 2 ) );
    const double cosineOfAzimuthAngle = std::cos( convertedSphericalState( 1 ) );
    const double sineOfAzimuthAngle = std::sin( convertedSphericalState( 1 ) );

    // Set up transformation matrix for cylindrical to spherical conversion.
    Eigen::MatrixXd transformationMatrixCylindricalToSpherical = Eigen::MatrixXd::Zero( 3, 3 );
    transformationMatrixCylindricalToSpherical( 0, 0 ) = cosineOfElevationAngle;
    transformationMatrixCylindricalToSpherical( 0, 2 ) = sineOfElevationAngle;
    transformationMatrixCylindricalToSpherical( 1, 1 ) = 1.0;
    transformationMatrixCylindricalToSpherical( 2, 0 ) = -sineOfElevationAngle;
    transformationMatrixCylindricalToSpherical( 2, 2 ) = cosineOfElevationAngle;

    // Set up transformation matrix for Cartesian to cylindrical conversion.
    Eigen::MatrixXd transformationMatrixCartesianToCylindrical = Eigen::MatrixXd::Zero( 3, 3 );
    transformationMatrixCartesianToCylindrical( 0, 0 ) = cosineOfAzimuthAngle;
    transformationMatrixCartesianToCylindrical( 0, 1 ) = sineOfAzimuthAngle;
    transformationMatrixCartesianToCylindrical( 1, 0 ) = -sineOfAzimuthAngle;
    transformationMatrixCartesianToCylindrical( 1, 1 ) = cosineOfAzimuthAngle;
    transformationMatrixCartesianToCylindrical( 2, 2 ) = 1.0;

    // Compute transformation matrix for Cartesian to spherical conversion.
    const Eigen::MatrixXd transformationMatrixCartesianToSpherical
            = transformationMatrixCylindricalToSpherical
            * transformationMatrixCartesianToCylindrical;

    // Perform transformation of velocity vector.
    convertedSphericalState.segment( 3, 3 )
            = transformationMatrixCartesianToSpherical * cartesianState.segment( 3, 3 );

    // Return spherical state vector.
    return convertedSphericalState;
}

} // namespace coordinate_conversions
} // namespace basic_mathematics
} // namespace tudat
