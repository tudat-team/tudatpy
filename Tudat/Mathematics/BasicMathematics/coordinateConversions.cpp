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
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *      Torok, J.S. Analytical Mechanics: with an Introduction to Dynamical Systems, John Wiley and
 *          Sons, Inc., 2000.
 *      Vallado, D.A. Fundamentals of Astrodynamics and Applications. Microcosm Press, 2001.
 *
 */

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/math/special_functions/sign.hpp>

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

namespace tudat
{

namespace coordinate_conversions
{


//! Convert spherical (radius_, zenith, azimuth) to Cartesian (x,y,z) coordinates.
Eigen::Vector3d convertSphericalToCartesian( const Eigen::Vector3d& sphericalCoordinates )
{
    // Create local variables.
    double radius_ = sphericalCoordinates( 0 );
    double zenithAngle_ = sphericalCoordinates( 1 );
    double azimuthAngle_ = sphericalCoordinates( 2 );

    // Declaring sine which has multiple usages to save computation time.
    double sineOfZenithAngle_ = std::sin( sphericalCoordinates( 1 ) );

    // Create output Vector3d.
    Eigen::Vector3d convertedCartesianCoordinates_ = Eigen::Vector3d::Zero( 3 );

    // Perform transformation.
    convertedCartesianCoordinates_( 0 ) = radius_ * std::cos( azimuthAngle_ ) * sineOfZenithAngle_;
    convertedCartesianCoordinates_( 1 ) = radius_ * std::sin( azimuthAngle_ ) * sineOfZenithAngle_;
    convertedCartesianCoordinates_( 2 ) = radius_ * std::cos( zenithAngle_ );

    return convertedCartesianCoordinates_;
}


//! Convert cylindrical to Cartesian coordinates.
Eigen::Vector3d convertCylindricalToCartesian( const double radius,
                                               const double azimuthAngle, const double z )
{
    // Create Cartesian coordinates vector.
    Eigen::Vector3d cartesianCoordinates;

    // If radius < 0, then give warning.
    if ( radius < 0.0 )
    {
        std::cerr << "Warning: cylindrical radial coordinate is negative!, This could give incorrect results!" << std::endl;
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
        std::cerr << "Warning: cylindrical radial coordinate is negative!, This could give incorrect results!" << std::endl;
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
Eigen::Vector6d convertCylindricalToCartesianState(
        const Eigen::Vector6d& cylindricalState )
{
    // Create Cartesian state vector, initialized with zero entries.
    Eigen::Vector6d cartesianState = Eigen::Vector6d::Zero( );

    // Get azimuth angle, theta.
    double azimuthAngle = cylindricalState( 1 );

    // Compute and set Cartesian coordinates.
    cartesianState.head( 3 ) = convertCylindricalToCartesian(
                Eigen::Vector3d( cylindricalState.head( 3 ) ) );

    // If r = 0 AND Vtheta > 0, then give warning and assume Vtheta=0.
    if ( std::fabs(cylindricalState( 0 )) <= std::numeric_limits< double >::epsilon( )
         && std::fabs(cylindricalState( 4 )) > std::numeric_limits< double >::epsilon( ) )
    {
        std::cerr << "Warning: cylindrical velocity Vtheta (r*thetadot) does not equal zero while the radius (r) is zero! Vtheta is taken equal to zero!" << std::endl;

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
    using mathematical_constants::PI;
    if ( std::fabs(cartesianCoordinates( 0 ) ) <= std::numeric_limits< double >::epsilon( ) )
    {
        azimuthAngle = basic_mathematics::computeModulo(
                    static_cast< double >( boost::math::sign( cartesianCoordinates( 1 ) ) )
                    * 0.5 * PI, 2.0 * PI );
    }

    else
    {
        azimuthAngle = basic_mathematics::computeModulo(
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
Eigen::Vector6d convertCartesianToCylindricalState(
        const Eigen::Vector6d& cartesianState )
{
    // Create cylindrical state vector, initialized with zero entries.
    Eigen::Vector6d cylindricalState = Eigen::Vector6d::Zero( );

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

//! Compute matrix by which to precompute a spherical gradient vector to obtain the Cartesian gradient
Eigen::Matrix3d getSphericalToCartesianGradientMatrix( const Eigen::Vector3d& cartesianCoordinates )
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
                - cartesianCoordinates( 0 ) * cartesianCoordinates( 2 ) / ( radius * radius * xyDistance ),
                - cartesianCoordinates( 1 ) / xyDistanceSquared,
                cartesianCoordinates( 1 ) / radius,
                - cartesianCoordinates( 1 ) * cartesianCoordinates( 2 ) / ( radius * radius * xyDistance ),
                + cartesianCoordinates( 0 ) / xyDistanceSquared,
                cartesianCoordinates( 2 ) / radius,
                xyDistance / ( radius * radius ),   0.0
                ).finished( );
    return transformationMatrix;
}

//! Convert spherical to Cartesian gradient.
Eigen::Vector3d convertSphericalToCartesianGradient( const Eigen::Vector3d& sphericalGradient,
                                                     const Eigen::Vector3d& cartesianCoordinates )
{


    // Return Cartesian gradient.
    return getSphericalToCartesianGradientMatrix( cartesianCoordinates ) * sphericalGradient;
}

Eigen::Matrix3d getDerivativeOfSphericalToCartesianGradient( const Eigen::Vector3d& sphericalGradient,
                                                             const Eigen::Vector3d& cartesianCoordinates,
                                                             std::vector< Eigen::Matrix3d >& subMatrices )
{
    Eigen::Matrix3d totalPartialMatrix;
    totalPartialMatrix.setZero( );

    Eigen::Matrix3d currentPartialMatrix;

    // Precomputed quantities
    double radius = cartesianCoordinates.norm( );
    const double xyDistanceSquared = cartesianCoordinates( 0 ) * cartesianCoordinates( 0 )
            + cartesianCoordinates( 1 ) * cartesianCoordinates( 1 );
    const double xyDistance = std::sqrt( xyDistanceSquared );
    const double radiusSquaredXyDistance = xyDistance * radius * radius;

    // Precompute partials
    Eigen::Vector3d oneOverRPartial = -cartesianCoordinates / ( radius * radius * radius );
    Eigen::Vector3d oneOverRSquaredPartial = -2.0 * cartesianCoordinates / ( radius * radius * radius * radius );
    Eigen::Vector3d oneOverXyDistancePartial =
            -( Eigen::Vector3d( ) << cartesianCoordinates( 0 ), cartesianCoordinates( 1 ), 0.0 ).finished( )/
            ( xyDistanceSquared * xyDistance );
    Eigen::Vector3d oneOverXyDistanceSquaredPartial =
            -2.0 * ( Eigen::Vector3d( ) << cartesianCoordinates( 0 ), cartesianCoordinates( 1 ), 0.0 ).finished( )/
            ( xyDistanceSquared * xyDistanceSquared );
    Eigen::Vector3d oneOverRSquaredXyDistancePartial =
            oneOverRSquaredPartial / xyDistance + oneOverXyDistancePartial / ( radius * radius );


    Eigen::Vector3d xyDistancePartial =
            ( Eigen::Vector3d( ) << cartesianCoordinates( 0 ), cartesianCoordinates( 1 ), 0.0 ).finished( ) / xyDistance;

    // Compute partials w.r.t x, y and z components.
    for( unsigned int i = 0; i < 3; i++ )
    {
        currentPartialMatrix.setZero( );
        switch( i )
        {
        case 0:
        {
            currentPartialMatrix << 1.0 / radius + cartesianCoordinates( 0 ) * oneOverRPartial ( 0 ),
                    - cartesianCoordinates( 2 ) / radiusSquaredXyDistance -
                    cartesianCoordinates( 0 ) * cartesianCoordinates( 2 ) * oneOverRSquaredXyDistancePartial( 0 ),
                    - cartesianCoordinates( 1 ) * oneOverXyDistanceSquaredPartial( 0 ),
                    cartesianCoordinates( 1 ) * oneOverRPartial( 0 ),
                     - cartesianCoordinates( 1 ) * cartesianCoordinates( 2 ) * oneOverRSquaredXyDistancePartial( 0 ),
                    1.0 / ( xyDistanceSquared ) + cartesianCoordinates( 0 ) * oneOverXyDistanceSquaredPartial( 0 ),
                     cartesianCoordinates( 2 ) * oneOverRPartial( 0 ),
                    xyDistance * oneOverRSquaredPartial( 0 ) + 1.0 / ( radius * radius ) * xyDistancePartial( 0 ),
                    0.0 ;
            break;
        }
        case 1:
        {
            currentPartialMatrix << cartesianCoordinates( 0 ) * oneOverRPartial ( 1 ),
                    - cartesianCoordinates( 0 ) * cartesianCoordinates( 2 ) * oneOverRSquaredXyDistancePartial( 1 ),
                    -1.0 / ( xyDistanceSquared ) - cartesianCoordinates( 1 ) * oneOverXyDistanceSquaredPartial( 1 ),
                    1.0 / radius + cartesianCoordinates( 1 ) * oneOverRPartial ( 1 ),
                    - cartesianCoordinates( 2 ) / radiusSquaredXyDistance -
                    cartesianCoordinates( 1 ) * cartesianCoordinates( 2 ) * oneOverRSquaredXyDistancePartial( 1 ),
                    cartesianCoordinates( 0 ) * oneOverXyDistanceSquaredPartial( 1 ),
                     cartesianCoordinates( 2 ) * oneOverRPartial( 1 ),
                    xyDistance * oneOverRSquaredPartial( 1 ) + 1.0 / ( radius * radius ) * xyDistancePartial( 1 ),
                    0.0 ;
            break;
        }
        case 2:
        {
            currentPartialMatrix <<  cartesianCoordinates( 0 ) * oneOverRPartial ( 2 ),
                    - cartesianCoordinates( 0 ) / radiusSquaredXyDistance -
                    cartesianCoordinates( 0 ) * cartesianCoordinates( 2 ) * oneOverRSquaredXyDistancePartial( 2 ),
                    - cartesianCoordinates( 1 ) * oneOverXyDistanceSquaredPartial( 2 ),
                    cartesianCoordinates( 1 ) * oneOverRPartial( 2 ),
                    - cartesianCoordinates( 1 ) / radiusSquaredXyDistance -
                    cartesianCoordinates( 1 ) * cartesianCoordinates( 2 ) * oneOverRSquaredXyDistancePartial( 2 ),
                    cartesianCoordinates( 0 ) * oneOverXyDistanceSquaredPartial( 2 ),
                     1.0 / radius + cartesianCoordinates( 2 ) * oneOverRPartial( 2 ),
                    xyDistance * oneOverRSquaredPartial( 2 ) +  + 1.0 / ( radius * radius ) * xyDistancePartial( 2 ),
                    0.0 ;
            break;
        }
        }

        // Save computed matrix
        if( subMatrices.size( ) == 3 )
        {
            subMatrices[ i ] = currentPartialMatrix;
        }

        // Add current entry to results.
        totalPartialMatrix.block( 0, i, 3, 1 ) = currentPartialMatrix * sphericalGradient;
    }

    return totalPartialMatrix;
}

Eigen::Matrix3d getDerivativeOfSphericalToCartesianGradient( const Eigen::Vector3d& sphericalGradient,
                                                             const Eigen::Vector3d& cartesianCoordinates )
{
    static  std::vector< Eigen::Matrix3d > subMatrices( 3 );
    return getDerivativeOfSphericalToCartesianGradient(
                sphericalGradient, cartesianCoordinates, subMatrices );
}

//! Convert spherical to Cartesian state.
Eigen::Vector6d convertSphericalToCartesianState(
        const Eigen::Vector6d& sphericalState )
{
    // Create Cartesian state vector, initialized with zero entries.
    Eigen::Vector6d convertedCartesianState = Eigen::Vector6d::Zero( );

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
    Eigen::Matrix3d transformationMatrixSphericalToCylindrical = Eigen::Matrix3d::Zero( );
    transformationMatrixSphericalToCylindrical( 0, 0 ) = cosineOfElevationAngle;
    transformationMatrixSphericalToCylindrical( 0, 2 ) = -sineOfElevationAngle;
    transformationMatrixSphericalToCylindrical( 1, 1 ) = 1.0;
    transformationMatrixSphericalToCylindrical( 2, 0 ) = sineOfElevationAngle;
    transformationMatrixSphericalToCylindrical( 2, 2 ) = cosineOfElevationAngle;

    // Set up transformation matrix for cylindrical to Cartesian conversion.
    Eigen::Matrix3d transformationMatrixCylindricalToCartesian = Eigen::Matrix3d::Zero( );
    transformationMatrixCylindricalToCartesian( 0, 0 ) = cosineOfAzimuthAngle;
    transformationMatrixCylindricalToCartesian( 0, 1 ) = -sineOfAzimuthAngle;
    transformationMatrixCylindricalToCartesian( 1, 0 ) = sineOfAzimuthAngle;
    transformationMatrixCylindricalToCartesian( 1, 1 ) = cosineOfAzimuthAngle;
    transformationMatrixCylindricalToCartesian( 2, 2 ) = 1.0;

    // Compute transformation matrix for spherical to Cartesian conversion.
    const Eigen::Matrix3d transformationMatrixSphericalToCartesian
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
Eigen::Vector6d convertCartesianToSphericalState(
        const Eigen::Vector6d& cartesianState )
{
    // Create spherical state vector, initialized with zero entries.
    Eigen::Vector6d convertedSphericalState = Eigen::Vector6d::Zero( );

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
    Eigen::Matrix3d transformationMatrixCylindricalToSpherical = Eigen::Matrix3d::Zero( );
    transformationMatrixCylindricalToSpherical( 0, 0 ) = cosineOfElevationAngle;
    transformationMatrixCylindricalToSpherical( 0, 2 ) = sineOfElevationAngle;
    transformationMatrixCylindricalToSpherical( 1, 1 ) = 1.0;
    transformationMatrixCylindricalToSpherical( 2, 0 ) = -sineOfElevationAngle;
    transformationMatrixCylindricalToSpherical( 2, 2 ) = cosineOfElevationAngle;

    // Set up transformation matrix for Cartesian to cylindrical conversion.
    Eigen::Matrix3d transformationMatrixCartesianToCylindrical = Eigen::Matrix3d::Zero( );
    transformationMatrixCartesianToCylindrical( 0, 0 ) = cosineOfAzimuthAngle;
    transformationMatrixCartesianToCylindrical( 0, 1 ) = sineOfAzimuthAngle;
    transformationMatrixCartesianToCylindrical( 1, 0 ) = -sineOfAzimuthAngle;
    transformationMatrixCartesianToCylindrical( 1, 1 ) = cosineOfAzimuthAngle;
    transformationMatrixCartesianToCylindrical( 2, 2 ) = 1.0;

    // Compute transformation matrix for Cartesian to spherical conversion.
    const Eigen::Matrix3d transformationMatrixCartesianToSpherical
            = transformationMatrixCylindricalToSpherical
            * transformationMatrixCartesianToCylindrical;

    // Perform transformation of velocity vector.
    convertedSphericalState.segment( 3, 3 )
            = transformationMatrixCartesianToSpherical * cartesianState.segment( 3, 3 );

    // Return spherical state vector.
    return convertedSphericalState;
}

} // namespace coordinate_conversions

} // namespace tudat
