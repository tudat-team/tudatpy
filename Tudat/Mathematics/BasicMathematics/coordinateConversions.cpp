/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
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
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 */

#include <boost/math/special_functions/sign.hpp>
#include <iostream>
#include <cmath>
#include <limits>
#include <numeric>
#include <TudatCore/Mathematics/BasicMathematics/coordinateConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace mathematics
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
    if ( std::fabs(cartesianCoordinates( 0 ) ) <= std::numeric_limits< double >::epsilon( ) )
    {
        azimuthAngle = tudat::mathematics::computeModulo(
                    static_cast< double >( boost::math::sign( cartesianCoordinates( 1 ) ) )
                    * 0.5 * PI, 2.0 * PI );
    }

    else
    {
        azimuthAngle = tudat::mathematics::computeModulo(
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

} // namespace coordinate_conversions
} // namespace mathematics
} // namespace tudat
