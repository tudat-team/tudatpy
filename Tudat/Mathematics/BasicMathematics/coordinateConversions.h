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

#ifndef TUDAT_COORDINATE_CONVERSIONS_H
#define TUDAT_COORDINATE_CONVERSIONS_H

#include <vector>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace coordinate_conversions
{

//! Convert spherical (radius, zenith, azimuth) to Cartesian (x,y,z) coordinates.
/*!
 * Converts spherical to cartesian coordinates. Schematic representation can be found on, e.g.,
 * http://mathworld.wolfram.com/SphericalCoordinates.html.
 * The transformation equations are the following, with \f$ r \f$ the radius,
 * \f$ \theta \f$ the azimuth angle and \f$ \phi \f$ the zenith angle:
 * \f{eqnarray*}{
 *      x &=& r\cos\theta\sin\phi \\
 *      y &=& r\sin\theta\sin\phi \\
 *      z &=& r\cos\phi \\
 * \f}
 * \param sphericalCoordinates Vector containing radius, zenith and azimuth (in that order).
 * \return Vector containing Cartesian coordinates, as calculated from sphericalCoordinates.
 */
Eigen::Vector3d convertSphericalToCartesian( const Eigen::Vector3d& sphericalCoordinates );

//! Convert Cartesian (x,y,z) to spherical (radius, zenith, azimuth) coordinates.
/*!
 * Converts Cartesian to spherical coordinates. Schematic representation can be found on, e.g.,
 * http://mathworld.wolfram.com/SphericalCoordinates.html.
 * The transformation equations are the following, with \f$ r \f$ the radius,
 * \f$ \theta \f$ the azimuth angle and \f$ \phi \f$ the zenith angle:
 * \f{eqnarray*}{
 *      r &=& \sqrt{ x^{ 2 } + y^{ 2 } + z^{ 2 } } \\
 *      \theta &=& \arctan\frac{ y }{ x } \\
 *      \phi &=& \arccos\frac{ z }{ r } \\
 * \f}
 * \param cartesianCoordinates Vector containing Cartesian coordinates.
 * \return Vector containing sphericalCoordinates radius, zenith and azimuth (in that
 *          order), as calculated from sphericalCoordinates.
*/
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 3, 1 > convertCartesianToSpherical( const Eigen::Matrix< ScalarType, 3, 1 > & cartesianCoordinates )

{
    // Create output Vector3d.
    Eigen::Matrix< ScalarType, 3, 1 > convertedSphericalCoordinates_ = Eigen::Matrix< ScalarType, 3, 1 >::Zero( );

    // Compute transformation of Cartesian coordinates to spherical coordinates.
    convertedSphericalCoordinates_( 0 ) = cartesianCoordinates.norm( );

    // Check if coordinates are at origin.
    if ( convertedSphericalCoordinates_( 0 ) < std::numeric_limits< ScalarType >::epsilon( ) )
    {
        convertedSphericalCoordinates_( 1 ) = mathematical_constants::getFloatingInteger< ScalarType >( 0 );
        convertedSphericalCoordinates_( 2 ) = mathematical_constants::getFloatingInteger< ScalarType >( 0 );
    }
    // Else compute coordinates using trigonometric relationships.
    else
    {
        convertedSphericalCoordinates_( 1 ) = std::acos( cartesianCoordinates( 2 )
                                                         / convertedSphericalCoordinates_( 0 ) );
        convertedSphericalCoordinates_( 2 ) = std::atan2( cartesianCoordinates( 1 ),
                                                          cartesianCoordinates( 0 ) );
    }

    return convertedSphericalCoordinates_;
}

//! Spherical coordinate indices.
/*!
  * Spherical coordinate indices, for position and velocity components. With r the radius, theta
  * the azimuthal angle, phi the elevational angle. Vr, Vtheta, Vphi are the velocities along the
  * corresponding base vectors of r, theta, phi.
  */
enum SphericalCoordinateIndices
{
    radiusSphericalCoordinateIndex,             // r
    azimuthSphericalCoordinateIndex,            // theta
    elevationSphericalCoordinateIndex,          // phi
    radialVelocitySphericalCoordinateIndex,     // Vr
    azimuthVelocitySphericalCoordinateIndex,    // Vtheta
    elevationVelocitySphericalCoordinateIndex   // Vphi
};

//! Cylindrical coordinate indices.
/*!
 * Cylindrical coordinate vector indices, for position and velocity components.
 */
enum CylindricalCoordinateIndices
{
    rCylindricalCoordinateIndex,
    thetaCylindricalCoordinateIndex,
    zCylindricalCoordinateIndex,
    rDotCylindricalCoordinateIndex,
    vThetaCylindricalCoordinateIndex,
    zDotCylindricalCoordinateIndex
};

//! Cartesian coordinate indices.
/*!
 * Cartesian coordinate vector indices, for position and velocity components.
 */
enum CartesianCoordinateIndices
{
    xCartesianCoordinateIndex,
    yCartesianCoordinateIndex,
    zCartesianCoordinateIndex,
    xDotCartesianCoordinateIndex,
    yDotCartesianCoordinateIndex,
    zDotCartesianCoordinateIndex
};

//! Convert cylindrical to Cartesian coordinates.
/*!
 * Converts cylindrical to Cartesian coordinates. Schematic representation can be found on, e.g.,
 * http://mathworld.wolfram.com/CylindricalCoordinates.html.
 * The transformation equations are the following, with \f$ r \f$ the radius and
 * \f$ \theta \f$ the azimuth angle [rad]:
 * \f{eqnarray*}{
 *      x &=& r\cos \theta \\
 *      y &=& r\sin \theta \\
 *      z &=& z \\
 * \f}
 * \param radius Cylindrical radial coordinate r.
 * \param azimuthAngle Cylindrical azimuthal coordinate \f$ \theta \f$ [rad].
 * \param z Cylindrical height coordinate z.
 * \return Vector of Cartesian coordinates [x,y,z].
 */
Eigen::Vector3d convertCylindricalToCartesian( const double radius,
                                               const double azimuthAngle, const double z );

//! Convert cylindrical to cartesian coordinates.
/*!
 * Converts cylindrical to cartesian coordinates. Schematic representation can be found on, e.g.,
 * http://mathworld.wolfram.com/CylindricalCoordinates.html.
 * The transformation equations are the following, with \f$ r \f$ the radius and
 * \f$ \theta \f$ the azimuth angle [rad]:
 * \f{eqnarray*}{
 *      x &=& r\cos \theta \\
 *      y &=& r\sin \theta \\
 *      z &=& z \\
 * \f}
 * \param cylindricalCoordinates Vector of cylindrical coordinates [r,theta,z].
 * \return Vector of Cartesian coordinates [x,y,z].
 */
Eigen::Vector3d convertCylindricalToCartesian( const Eigen::Vector3d& cylindricalCoordinates );

//! Convert cylindrical to Cartesian state.
/*!
 * Converts cylindrical to Cartesian state. Schematic representation can be found on, e.g.,
 * http://mathworld.wolfram.com/CylindricalCoordinates.html and
 * http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node26.html.
 * The transformation equations are the following, with \f$ r \f$ the radius and
 * \f$ \theta \f$ the azimuth angle [rad]:
 * \f{eqnarray*}{
 *      x &=& r\cos \theta \\
 *      y &=& r\sin \theta \\
 *      z &=& z \\
 *      \dot{x} &=& V_r\cos{\theta} - V_{\theta}\sin{theta} \\
 *      \dot{y} &=& V_r\sin{\theta} + V_{\theta}\cos{theta} \\
 *      \dot{z} &=& V_z
 * \f}
 * \param cylindricalState Vector of cylindrical state [r,theta,z,Vr,Vtheta,Vz],
 *           where Vtheta = r*thetadot.
 * \return Vector of Cartesian state [x,y,z,xdot,ydot,zdot].
 */
Eigen::Vector6d convertCylindricalToCartesianState(
        const Eigen::Vector6d& cylindricalState );

//! Convert Cartesian to cylindrical coordinates.
/*!
* Converts Cartesian to cylindrical coordinates. Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/CylindricalCoordinates.html.
* The transformation equations are the following, with \f$ r \f$ the radius and
* \f$ \theta \f$ the azimuth angle [rad] [0,2\f$ \pi \f$]:
* \f{eqnarray*}{
*      r &=& \sqrt{x^2+y^2} \\
*      \theta &=& \arctan{\frac{y}{x}} \\
*      z &=& z
* \f}
* \param cartesianCoordinates Vector of Cartesian coordinates [x,y,z].
* \return Vector of cylindrical coordinates [r,theta,z].
*/
Eigen::Vector3d convertCartesianToCylindrical( const Eigen::Vector3d& cartesianCoordinates );

//! Convert Cartesian to cylindrical state.
/*!
* Converts Cartesian to cylindrical state. Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/CylindricalCoordinates.html and
* http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node26.html.
* The transformation equations are the following, with \f$ r \f$  the radius,
* \f$ \theta \f$ the azimuth angle [rad] [0,2\f$\pi\f$] and \f$ V_r \f$, \f$ V_{\theta} \f$ and
* \f$ V_z \f$ the linear cylindrical velocities:
* \f{eqnarray*}{
*      r &=& \sqrt{x^2+y^2} \\
*      \theta &=& \arctan{\frac{y}{x}} \\
*      z &=& z \\
*      V_r = \dot{r} &=& \frac{x\dot{x}+y\dot{y}}{\sqrt{x^2+y^2}} \\
*      V_{\theta} = r\dot{\theta} &=& \frac{x\dot{y}-y\dot{x}}{\sqrt{x^2+y^2}} \\
*      V_z = \dot{z}
* \f}
* \param cartesianState Vector of Cartesian state [x,y,z,xdot,ydot,zdot].
* \return Vector of cylindrical state [r,theta,z,Vr,Vtheta,Vz], where Vtheta = r*thetadot.
*/
Eigen::Vector6d convertCartesianToCylindricalState(
        const Eigen::Vector6d& cartesianState );

//! Compute matrix by which to precompute a spherical gradient vector to obtain the Cartesian gradient
/*!
 * Compute matrix by which to precompute a spherical gradient vector to obtain the Cartesian gradient
 * \param cartesianCoordinates Vector with Cartesian position at which gradient is computed.
 * \return Matrix by which to precompute a spherical gradient vector to obtain the Cartesian gradient
 */
Eigen::Matrix3d getSphericalToCartesianGradientMatrix( const Eigen::Vector3d& cartesianCoordinates );

//! Convert spherical to Cartesian gradient.
/*!
* This function converts a gradient vector with respect to spherical coordinates to a gradient
* vector with respect to Cartesian coordinates. The partial derivatives are calculated according
* to Vallado [2001] as:
* \f{eqnarray*}{
*   \frac{ \partial U }{ \partial x } & = &
*   \frac{ x }{ \sqrt{ x^2 + y^2 +z^2 } } \frac{ \partial U }{ \partial r }
*   - \frac{ x z }{ ( x^2 + y^2 +z^2 ) \sqrt{ x^2 + y^2 } } \frac{ \partial U }{ \partial \phi }
*   - \frac{ y }{ x^2 + y^2 } \frac{ \partial U }{ \partial \lambda } \\
*   \frac{ \partial U }{ \partial y } & = &
*   \frac{ y }{ \sqrt{ x^2 + y^2 +z^2 } } \frac{ \partial U }{ \partial r }
*   - \frac{ y z }{ ( x^2 + y^2 +z^2 ) \sqrt{ x^2 + y^2 } } \frac{ \partial U }{ \partial \phi }
*   + \frac{ x }{ x^2 + y^2 } \frac{ \partial U }{ \partial \lambda } \\
*   \frac{ \partial U }{ \partial z } & = &
*   \frac{ z }{ \sqrt{ x^2 + y^2 +z^2 } } \frac{ \partial U }{ \partial r }
*   + \frac{ \sqrt{ x^2 + y^2 } } { x^2 + y^2 +z^2 } \frac{ \partial U }{ \partial \phi }
* \f}
* in which \f$ x \f$, \f$ y \f$ and \f$ z \f$ are the Cartesian coordinates. \f$ U \f$ is an
* arbitrary scalar. Radius \f$ r \f$ is the length of the position vector. Elevation \f$ \phi \f$
* is the the angle between the position vector and the XY-plane (with the positive direction
* towards the Z-axis). Azimuth \f$ \lambda \f$ is the dihedral angle about the Z-axis between the
* X-axis and the position vector (with the positive direction being right-handed about the Z-axis).
* \param sphericalGradient Vector with partial derivatives with respect to spherical coordinates.
*        The order is important!
*        sphericalGradient( 0 ) = partial derivative with respect to radius,
*        sphericalGradient( 1 ) = partial derivative with respect to elevation,
*        sphericalGradient( 2 ) = partial derivative with respect to azimuth.
* \param cartesianCoordinates Vector with Cartesian coordinates.
*        The order is important!
*        cartesianCoordinates( 0 ) = x coordinate,
*        cartesianCoordinates( 1 ) = y coordinate,
*        cartesianCoordinates( 2 ) = z coordinate.
* \return Vector with partial derivatives with respect to Cartesian coordinates.
*         The order is important!
*         cartesianGradient( 0 ) = partial derivative with respect to x coordinate,
*         cartesianGradient( 1 ) = partial derivative with respect to y coordinate,
*         cartesianGradient( 2 ) = partial derivative with respect to z coordinate.
*/
Eigen::Vector3d convertSphericalToCartesianGradient( const Eigen::Vector3d& sphericalGradient,
                                                     const Eigen::Vector3d& cartesianCoordinates );

//! Function to compute the derivative of the Cartesian gradient w.r.t. the  Cartesian position, keeping the
//! spherical gradient constant.
/*!
 * Function to compute the derivative of the Cartesian gradient w.r.t. the  Cartesian position, keeping the
 * spherical gradient constant.
 * \param sphericalGradient Value of spherical gradient (derivatives w.r.t. radius, latitude and longitude).
 * \param cartesianCoordinates Cartesian position at whichr result is to be computed.
 * \param subMatrices Subcomputations performed by this function: derivatives of results of
 * getSphericalToCartesianGradientMatrix, w.r.t. x, y and z component of cartesianCoordinates, respectively.
 * \return Require derivative of cartesian gradient.
 */
Eigen::Matrix3d getDerivativeOfSphericalToCartesianGradient( const Eigen::Vector3d& sphericalGradient,
                                                             const Eigen::Vector3d& cartesianCoordinates,
                                                             std::vector< Eigen::Matrix3d >& subMatrices );

//! Function to compute the derivative of the Cartesian gradient w.r.t. the  Cartesian position, keeping the
//! spherical gradient constant.
/*!
 * Function to compute the derivative of the Cartesian gradient w.r.t. the  Cartesian position, keeping the
 * spherical gradient constant.
 * \param sphericalGradient Value of spherical gradient (derivatives w.r.t. radius, latitude and longitude).
 * \param cartesianCoordinates Cartesian position at whichr result is to be computed..
 * \return Require derivative of cartesian gradient.
 */
Eigen::Matrix3d getDerivativeOfSphericalToCartesianGradient( const Eigen::Vector3d& sphericalGradient,
                                                             const Eigen::Vector3d& cartesianCoordinates );

//! Convert spherical to Cartesian state.
/*!
  * Converts a spherical state to a Cartesian state. The transformation matrices are computed
  * according to Torok [2000, pp.10-11].
  *
  * NOTE: This function is implemented separately from the other conversions due to a different
  * definition of the elevation/zenith angle. This should be consolidated in a future update.
  *
  * The transformation equations are the following, with \f$ r \f$ the radius (positive from origin
  * to the point in orbit, in meters), \f$ \theta \f$ the azimuth angle (positive from the x-axis
  * to the y-axis, in radians) and \f$ \phi \f$ the elevation angle (positive from the xy-plane
  * to the z-axis, in radians):
  * \f{eqnarray*}{
  *     x &=& r * \cos \phi * \cos \theta\\
  *     y &=& r * \cos \phi * \sin \theta\\
  *     z &=& r * \sin \phi\\
  *     CartesianVelocities = T_{cyl2cart}*T_{sph2cyl}*SphericalVelocities
  * \f}
  * with
  * \f{eqnarray*}{
  *     T_{sph2cyl} = [     \cos\phi    , 0.0   ,    -\sin \phi  ;
  *                         0.0         , 1.0   ,    0.0         ;
  *                         \sin\phi    , 0.0   ,    \cos\phi    ]\\
  *     T_{cyl2cart} = [    \cos\theta  , -\sin\theta   , 0.0   ;
  *                         \sin\theta  , \cos\theta    , 0.0   ;
  *                         0.0         , 0.0           , 1.0   ]\\
  * \f}
  *
  * \param sphericalState Vector containing the spherical coordinates and spherical velocities.
  *        The order is important!
  *        sphericalState( 0 ) = radius r [m],
  *        sphericalState( 1 ) = azimuth theta [rad],
  *        sphericalState( 2 ) = elevation phi [rad],
  *        sphericalState( 3 ) = radial velocity Vr [m/s],
  *        sphericalState( 4 ) = azimuthal velocity Vtheta [m/s],
  *        sphericalState( 5 ) = elevational velocity Vphi [m/s].
  * \return Vector containing the Cartesian state (both position and velocity, in that order).
  *         cartesianState( 0 ) = x [m],
  *         cartesianState( 1 ) = y [m],
  *         cartesianState( 2 ) = z [m],
  *         cartesianState( 3 ) = Vx [m/s],
  *         cartesianState( 4 ) = Vy [m/s],
  *         cartesianState( 5 ) = Vz [m/s].
  *
  * Take care: here the elevation is used, not the zenith angle!
  */
Eigen::Vector6d convertSphericalToCartesianState(
        const Eigen::Vector6d& sphericalState );

//! Convert Cartesian to spherical state.
/*!
  * Converts a Cartesian state to a spherical state. The transformation matrices are computed
  * according to Torok [2000, pp.10-11].
  *
  * NOTE: This function is implemented separately from the other conversions due to a different
  * definition of the elevation/zenith angle. This should be consolidated in a future update.
  *
  * The transformation equations are the following, with \f$ r \f$ the radius (positive from origin
  * to the point in orbit), \f$ \theta \f$ the azimuth angle (positive from the x-axis to the
  * y-axis) and \f$ \phi \f$ the elevation angle (positive from the xy-plane to the z-axis):
  * \f{eqnarray*}{
  *     r &=& \sqrt{ x^{ 2 } + y^{ 2 } + z^{ 2 } } \\
  *     \phi &=& \arcsin\frac{ z }{ r } \\
  *     \theta &=& \arctan\frac{ y }{ x } \\
  *     SphericalVelocities = T_{cyl2sph}*T_{cart2cyl}*CartesianVelocities
  * \f}
  * with
  * \f{eqnarray*}{
  *     T_{cyl2sph} = [     \cos\phi    , 0.0   ,    \sin \phi   ;
  *                         0.0         , 1.0   ,    0.0         ;
  *                         -\sin\phi   , 0.0   ,    \cos\phi    ]\\
  *     T_{cart2cyl} = [    \cos\theta  , \sin\theta    , 0.0   ;
  *                         -\sin\theta , \cos\theta    , 0.0   ;
  *                         0.0         , 0.0           , 1.0   ]\\
  * \f}
  *
  * \param cartesianState Vector containing the Cartesian coordinates and Cartesian velocities.
  *        The order is important!
  *        cartesianState( 0 ) = x [m],
  *        cartesianState( 1 ) = y [m],
  *        cartesianState( 2 ) = z [m],
  *        cartesianState( 3 ) = Vx [m/s],
  *        cartesianState( 4 ) = Vy [m/s],
  *        cartesianState( 5 ) = Vz [m/s].
  * \return Vector containing the spherical state (both position and velocity, in that order).
  *         The order is important!
  *         sphericalState( 0 ) = radius r [m],
  *         sphericalState( 1 ) = azimuth theta [rad],
  *         sphericalState( 2 ) = elevation phi [rad],
  *         sphericalState( 3 ) = radial velocity Vr [m/s],
  *         sphericalState( 4 ) = azimuthal velocity Vtheta [m/s],
  *         sphericalState( 5 ) = elevational velocity Vphi [m/s].
  *
  * Take care: here the elevation is used, not the zenith!
  */
Eigen::Vector6d convertCartesianToSphericalState(
        const Eigen::Vector6d& cartesianState );

} // namespace coordinate_conversions

} // namespace tudat

#endif // TUDAT_COORDINATE_CONVERSIONS_H
