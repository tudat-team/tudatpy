/*! \file basicMathematicsFunctions.h
 *    Header file that defines the basicMathematicsFunctions namespace,
 *    containing all basic functions contained in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 13
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Author            : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Author            : D.J. Gondelach
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.J.Gondelach@student.tudelft.nl
 *
 *    Checker           : L. Abdulkadir
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : T. Secretin
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : T.A.LeitePintoSecretin@student.tudelft.nl
 *
 *    Date created      : 3 September, 2010
 *    Last modified     : 18 January, 2012
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *      Spiegel, M.R., Stephens, L.J. Statistics, Fourth Edition, Schaum's
 *          Outlines, McGraw-Hill, 2008.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      120118    D. Gondelach      Added convertCylindricalToCartesianCoordinates,
 *                                  convertCylindricalToCartesianState,
 *                                  convertCartesianToCylindricalCoordinates and
 *                                  convertCartesianToCylindricalState functions.
 *                                  Removed convertCylindricalToCartesian function.
 */

#ifndef BASICMATHEMATICSFUNCTIONS_H
#define BASICMATHEMATICSFUNCTIONS_H

// Include statements.
#include <boost/math/special_functions/sign.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <ctime>
#include <Eigen/Core>
#include <map>
#include <vector>
#include "Astrodynamics/States/state.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Mathematics namespace.
/*!
 *  Mathematics namespace.
 */
namespace mathematics
{

//! Random number generator typedef.
/*!
 * Random number generator typedef. This can be modified to any other Boost random number
 * generator type
 * (http://www.boost.org/doc/libs/1_47_0/doc/html/boost_random/reference.html.
 */
typedef boost::mt19937 globalRandomNumberGeneratorType;

//! Get global random number generator.
/*!
 * Returns global random number generator. The default seed is set to the current time.
 * \return Global random number generator.
 */
globalRandomNumberGeneratorType& getGlobalRandomNumberGenerator( );

//! Compute linear interpolation.
/*!
 * Computes linear interpolation of data provided in the form of a vector of
 * sorted indepedent variables and an associated vector of dependent variables.
 * The linear interpolation equation used is:
 * \f[
 *      y_{target} = x_{1} * ( 1 - mu ) + x_{2} * mu
 * \f]
 * where \f$ mu = \frac{ x_{target} - x_{1} } { x_{2} + x_{1} } \f$
 * and \f$ x_{2} > x_{1} \f$.
 * \param sortedIndependentVariables Vector of independent variables sorted.
 *          in ascending/descending order.
 * \param associatedDependentVariables Vector of dependent variables
 *          associated with sorted vector of independent variables.
 * \param targetIndependentVariableValue Target independent variable value
 *          in vector of sorted independent variables.
 * \return Value of dependent variable associated with target independent
 *          value in vector of sorted independent variables.
 */
double computeLinearInterpolation( Eigen::VectorXd& sortedIndependentVariables,
                                   Eigen::VectorXd& associatedDependentVariables,
                                   double targetIndependentVariableValue );
//! Compute linear interpolation.
/*!
 * Computes linear interpolation of data provided in the form of a map of
 * independent variables and associated vectors of dependent variables.
 * The linear interpolation equation used is:
 * \f[
 *      y_{target} = x_{1} * ( 1 - mu ) + x_{2} * mu
 * \f]
 * where \f$ \mu = \frac{ x_{target} - x_{1} } { x_{2} + x_{1} } \f$
 * and \f$ x_{2} > x_{1} \f$.
 * \param sortedIndepedentAndDependentVariables Map of sorted independent
 *              variables, in ascending/descending order, and associated
 *              dependent variables.
 * \param targetIndependentVariableValue Target independent variable value
 *              in vector of sorted independent variables.
 * \return Vector of dependent variable associated with target independent
 *              value in vector of sorted independent variables.
 */
Eigen::VectorXd computeLinearInterpolation(
        std::map < double, Eigen::VectorXd >& sortedIndepedentAndDependentVariables,
        double targetIndependentVariableValue );

//! Compute linear interpolation.
/*!
 * Computes linear interpolation of data provided in the form of a map of
 * sorted independent variables and associated State objects containing vectors
 * of dependent variables. The linear interpolation equation used is:
 * \f[
 *      y_{target} = x_{1} * ( 1 - mu ) + x_{2} * mu
 * \f]
 * where \f$ \mu = \frac{ x_{target} - x_{1} } { x_{2} + x_{1} } \f$
 * and \f$ x_{2} > x_{1} \f$.
 * \param sortedIndepedentAndDependentVariables Map of sorted independent
 *              variables, in ascending/descending order, and associated
 *              State objects.
 * \param targetIndependentVariableValue Target independent variable value
 *              in vector of sorted independent variables.
 * \return Vector of dependent variable associated with target independent
 *              value in vector of sorted independent variables.
 */
State* computeLinearInterpolation(
        std::map < double, State* >& sortedIndepedentAndDependentVariables,
        double targetIndependentVariableValue );

//! Convert spherical to cartesian coordinates.
/*!
* Function to convert spherical to cartesian coordinates.
* Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/SphericalCoordinates.html.
* The transformation equations are the following, with \f$ r \f$ the radius,
* \f$ \theta \f$ the azimuth angle and \f$ \phi \f$ the zenith angle:
* \f{eqnarray*}{
*      x &=& r\cos\theta\sin\phi \\
*      y &=& r\sin\theta\sin\phi \\
*      z &=& r\cos\phi \\
* \f}
*/
void convertSphericalToCartesian( double radius, double azimuthAngle, double zenithAngle,
                                  Eigen::VectorXd& cartesianCoordinates );

//! Convert cartesian to spherical coordinates.
/*!
* Converts cartesian to spherical coordinates.
* Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/SphericalCoordinates.html.
* The transformation equations are the following, with \f$ r \f$ the radius,
* \f$ \theta \f$ the azimuth angle and \f$ \phi \f$ the azimuth angle:
* \f{eqnarray*}{
*      r &=& \sqrt{ x^{ 2 } + y^{ 2 } + z^{ 2 } } \\
*      \theta &=& \arctan\frac{ y }{ x } \\
*      \phi &=& \arccos\frac{ z }{ r } \\
* \f}
*/
void convertCartesianToSpherical( const Eigen::VectorXd& cartesianCoordinates,
                                  Eigen::VectorXd& sphericalCoordinates );

//! Convert cylindrical to Cartesian coordinates.
/*!
* Function to convert cylindrical to Cartesian coordinates.
* Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/CylindricalCoordinates.html.
* The transformation equations are the following, with \f$ r \f$ the radius and
* \f$ \theta \f$ the azimuth angle [rad]:
* \f{eqnarray*}{
*      x &=& r\cos\theta \\
*      y &=& r\sin\theta \\
*      z &=& z \\
* \f}
* \param radius Cylindrical radial coordinate r.
* \param azimuthAngle Cylindrical azimuthal coordinate \theta [rad].
* \param z Cylindrical height coordinate z.
* \return Vector of Cartesian coordinates [x,y,z].
*/
Eigen::Vector3d convertCylindricalToCartesianCoordinates(
        double radius, double azimuthAngle, double z );

//! Convert cylindrical to cartesian coordinates.
/*!
* Function to convert cylindrical to cartesian coordinates.
* Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/CylindricalCoordinates.html.
* The transformation equations are the following, with \f$ r \f$ the radius and
* \f$ \theta \f$ the azimuth angle [rad]:
* \f{eqnarray*}{
*      x &=& r\cos\theta \\
*      y &=& r\sin\theta \\
*      z &=& z \\
* \f}
* \param cylindricalCoordinates Vector of cylindrical coordinates [r,theta,z].
* \return Vector of Cartesian coordinates [x,y,z].
*/
Eigen::Vector3d convertCylindricalToCartesianCoordinates( Eigen::Vector3d cylindricalCoordinates );

//! Convert cylindrical to Cartesian state.
/*!
* Function to convert cylindrical to Cartesian state.
* Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/CylindricalCoordinates.html and
* http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node26.html.
* The transformation equations are the following, with \f$ r \f$ the radius and
* \f$ \theta \f$ the azimuth angle [rad]:
* \f{eqnarray*}{
*      x &=& r\cos\theta \\
*      y &=& r\sin\theta \\
*      z &=& z \\
*      \dot{x} &=& V_r\cos{\theta} - V_{\theta}\sin{theta} \\
*      \dot{y} &=& V_r\sin{\theta} + V_{\theta}\cos{theta} \\
*      \dot{z} &=& V_z
* \f}
* \param cylindricalState Vector of cylindrical state [r,theta,z,Vr,Vtheta,Vz],
*           where Vtheta = r*thetadot.
* \return Vector of Cartesian state [x,y,z,xdot,ydot,zdot].
*/
Eigen::VectorXd convertCylindricalToCartesianState( Eigen::VectorXd cylindricalState );

//! Convert Cartesian to cylindrical coordinates.
/*!
* Function to convert Cartesian to cylindrical coordinates.
* Schematic representation can be found on, e.g.,
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
Eigen::Vector3d convertCartesianToCylindricalCoordinates( Eigen::Vector3d cartesianCoordinates );

//! Convert Cartesian to cylindrical state.
/*!
* Function to convert Cartesian to cylindrical state.
* Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/CylindricalCoordinates.html and
* http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node26.html.
* The transformation equations are the following, with \f$ r \f$  the radius,
* \f$ \theta \f$ the azimuth angle [rad] [0,2\pi] and \f$ V_r \f$, \f$ V_{\theta} \f$ and
* \f$ V_z \f$ the linear cylindrical velocities:
* \f{eqnarray*}{
*      r &=& \sqrt{x^2+y^2} \\
*      \theta &=& \arctan{\frac{y}{x}} \\
*      z &=& z \\
*      V_r = \dot{r} &=& \frac{x\dot{x}+y\dot{y}}{\sqrt{x^2+y^2}} \\
*      V_{\theta} = r\dot{\theta} &=& \frac{x\dot{y}-y\dot{x}}{\sqrt{x^2+y^2}} \\
*      V_z = \dot{z}
* \f}
* \param cylindricalState Vector of Cartesian state [x,y,z,xdot,ydot,zdot].
* \return Vector of cylindrical state [r,theta,z,Vr,Vtheta,Vz], where Vtheta = r*thetadot.
*/
Eigen::VectorXd convertCartesianToCylindricalState( Eigen::VectorXd cartesianState );

//! Compute modulo of double.
/*!
 * Computes the remainder of division of one number by another.
 * The remainder is in the range [ 0, divisor ].
 * \param dividend Number to be divided.
 * \param divisor Number that is divided by.
 * \return Remainder of division of dividend by divisor.
 */
double computeModulo( double dividend, double divisor );

//! Compute sample mean.
/*!
 * Computes sample mean based on the following unbiased estimator
 * (Spiegel and Stephens, 2008):
 * \f[
 *      \mu_{s} = \frac{ \sum_{i=1}^{N} X_{i} } { N }
 * \f]
 * where \f$\mu_{s}\f$ is the unbiased estimate of the sample mean,
 * \f$ N \f$ is the number of samples, and \f$ X \f$ is the sample value.
 * \param sampleData Sample data.
 * \return Sample mean.
 */
double computeSampleMean( const std::vector< double >& sampleData );

//! Compute sample variance.
/*!
 * Computes sample variance based on the following unbiased estimator
 * (Spiegel and Stephens, 2008):
 * \f[
 *      s^{2}_{s} = \frac{ 1 }{ N - 1 } * \sum_{i=1}^{N} X_{i}
 *                  ( X_{i} - \bar{ X } )^{ 2 } )
 * \f]
 * where \f$ s^{2}_{s} \f$ is the unbiased estimate of the sample variance,
 * \f$ N \f$ is the number of samples, \f$ X \f$ is the sample value, and
 * \f$ \bar{ X } \f$ is the sample mean.
 * \param sampleData Map containing sample data.
 * \return Sample variance.
 */
double computeSampleVariance( const std::vector< double >& sampleData );

}

}

#endif // BASICMATHEMATICSFUNCTIONS_H

// End of file.
