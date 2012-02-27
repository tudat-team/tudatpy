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
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 */

#ifndef TUDAT_COORDINATE_CONVERSIONS_H
#define TUDAT_COORDINATE_CONVERSIONS_H

#include <Eigen/Core>

namespace tudat
{
namespace mathematics
{
namespace coordinate_conversions
{

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

} // namespace coordinate_conversions
} // namespace mathematics
} // namespace tudat

#endif // TUDAT_COORDINATE_CONVERSIONS_H
