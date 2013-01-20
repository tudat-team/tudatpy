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
 *      100903    K. Kumar          File header and footer added.
 *      100916    L. Abdulkadir     File checked.
 *      100929    K. Kumar          Checked code by D. Dirkx added.
 *      101110    K. Kumar          Added raiseToIntegerExponent() function.
 *      102410    D. Dirkx          Minor comment changes during code check.
 *      101213    K. Kumar          Modified raiseToIntegerExponent() function;
 *                                  renamed raiseToIntegerPower().
 *                                  Added computeAbsoluteValue() functions.
 *      110111    J. Melman         Added computeModulo() function.
 *      110202    K. Kumar          Added overload for State* for computeLinearInterpolation().
 *      110411    K. Kumar          Added convertCartesianToSpherical() function.
 *      110707    K. Kumar          Added computeSampleMean(), computeSampleVariance() functions.
 *      110810    J. Leloux         Corrected doxygen documentation (equations).
 *      110824    J. Leloux         Corrected doxygen documentation.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120127    D. Dirkx          First version branched from basic mathematics in Tudat Core.
 *      120127    K. Kumar          Minor comment edits.
 *      120118    D. Gondelach      Added new convertCylindricalToCartesian functions.
 *      120214    K. Kumar          Branched from old Tudat trunk for new coordinate conversions.
 *      120511    K. Kumar          Added enums for cylindrical and Cartesian coordinates.
 *      120926    E. Dekens         Added spherical gradient to Cartesian conversion.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *      Vallado, D.A. Fundamentals of Astrodynamics and Applications. Microcosm Press, 2001.
 *
 *    Notes
 *
 */

#ifndef TUDAT_COORDINATE_CONVERSIONS_H
#define TUDAT_COORDINATE_CONVERSIONS_H

#include <Eigen/Core>

namespace tudat
{
namespace basic_mathematics
{
namespace coordinate_conversions
{

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
    thetaDotCylindricalCoordinateIndex,
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
 * \param azimuthAngle Cylindrical azimuthal coordinate \theta [rad].
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
Eigen::VectorXd convertCylindricalToCartesian( const Eigen::VectorXd& cylindricalState );

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
Eigen::VectorXd convertCartesianToCylindrical( const Eigen::VectorXd& cartesianState );

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
Eigen::Vector3d convertSphericalToCartesianGradient( const Eigen::Vector3d sphericalGradient,
                                                     const Eigen::Vector3d cartesianCoordinates );

} // namespace coordinate_conversions
} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_COORDINATE_CONVERSIONS_H
