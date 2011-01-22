/*! \file basicMathematicsFunctions.h
 *    Header file that defines the basicMathematicsFunctions namespace,
 *    containing all basic functions contained in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 6
 *    Check status      : Unchecked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : L. Abdulkadir
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 3 September, 2010
 *    Last modified     : 13 December, 2010
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      YYMMDD    author              comment
 *      100903    K. Kumar            File header and footer added.
 *      100916    L. Abdulkadir       File checked.
 *      100929    K. Kumar            Checked code by D. Dirkx added.
 *      101110    K. Kumar            Added raiseToIntegerExponent() function.
 *      102410    D. Dirkx            Minor comment changes during code check.
 *      101213    K. Kumar            Modified raiseToIntegerExponent()
 *                                    function; renamed raiseToIntegerPower().
 *                                    Added computeAbsoluteValue() functions.
 */

#ifndef BASICMATHEMATICSFUNCTIONS_H
#define BASICMATHEMATICSFUNCTIONS_H

// Include statements.
#include <map>
#include <cmath>
#include <cfloat>
#include "basicFunctions.h"
#include "linearAlgebra.h"

//! Mathematics namespace.
/*!
 *  Mathematics namespace.
 */
namespace mathematics
{

//! Machine precision for floats.
/*!
 * Machine precision for floats from <float>.
 */
const static double MACHINE_PRECISION_FLOATS = FLT_EPSILON;

//! Machine precision for doubles.
/*!
 * Machine precision for doubles from <float>.
 */
const static double MACHINE_PRECISION_DOUBLES = DBL_EPSILON;

//! Machine precision for long doubles.
/*!
 * Machine precision for long doubles from <float>.
 */
const static double MACHINE_PRECISION_LONG_DOUBLES = LDBL_EPSILON;

//! Compute linear interpolation.
/*!
 * This function computes linear interpolation of data provided in the form of
 * a vector of sorted indepedent variables and an associated vector of
 * dependent variables. The linear interpolation equation used is:
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
double computeLinearInterpolation( VectorXd& sortedIndependentVariables,
                                   VectorXd& associatedDependentVariables,
                                   double& targetIndependentVariableValue );
//! Compute linear interpolation.
/*!
 * This function computes linear interpolation of data provided in the form of
 * a vector of sorted indepedent variables and an associated vector of
 * dependent variables. The linear interpolation equation used is:
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
VectorXd computeLinearInterpolation(
        std::map < double, VectorXd >& sortedIndepedentAndDependentVariables,
        double& targetIndependentVariableValue );

//! Convert spherical to cartesian coordinates.
/*!
* Function to convert spherical to cartesian coordinates.
* Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/SphericalCoordinates.html.
* The transformation equations are the following, with \f$ \r \f$ the radius,
* \f$ \theta \f$ the azimuth angle and \f$ \phi \f$ the azimuth angle:
* \f[
*      x=r\cos\theta\sin\phi \\
*      y=r\sin\theta\sin\phi \\
*      z=r\cos\phi \\
* \f]
*/
void convertSphericalToCartesian( const double& radius,
                                  const double& azimuthAngle,
                                  const double& zenithAngle,
                                  VectorXd& cartesianCoordinates );

//! Convert cylindrical to cartesian coordinates, z value unaffected.
/*!
* Function to convert cylindrical to cartesian coordinates.
* Schematic representation can be found on, e.g.,
* http://mathworld.wolfram.com/CylindricalCoordinates.html.
* The transformation equations are the following, with \f$ \r \f$  the radius and
* \f$ \theta \f$ the azimuth angle:
* \f[
*      x=r\cos\theta \\
*      y=r\sin\theta \\
*      z=z \\
* \f]
* Since the value of z is left unaffected by this transformation,
* it is not set or changed by this function.
*/
void convertCylindricalToCartesian( const double& radius,
                                    const double& azimuthAngle,
                                    VectorXd& cartesianCoordinates );

//! Raise double to integer power.
/*!
 * This function computes the result of raising floating-point base values to
 * integer powers. This function seems to perform faster than the STL pow()
 * function included in <cmath>. This function is based on exponentiation by
 * squares.
 */
double raiseToIntegerPower( const double& baseValue,
                            const int& integerPower );

//! Compute absolute value of integer.
/*!
 * This function computes the absolute value of an integer. This function seems
 * to perform faster than the STL abs() function included in <cmath>
 */
int computeAbsoluteValue( const int& signedInteger );

//! Compute absolute value of double.
/*!
 * This function computes the absolute value of an double. This function seems
 * to perform faster than the STL abs() function included in <cmath>
 */
double computeAbsoluteValue( const double& signedDouble );

}

#endif // BASICMATHEMATICSFUNCTIONS_H

// End of file.
