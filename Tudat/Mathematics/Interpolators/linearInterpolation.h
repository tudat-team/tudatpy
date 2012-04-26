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
 *      101110    K. Kumar          Added raiseToIntegerExponent( ) function.
 *      102410    D. Dirkx          Minor comment changes during code check.
 *      101213    K. Kumar          Modified raiseToIntegerExponent( ) function;
 *                                  renamed raiseToIntegerPower( ).
 *                                  Added computeAbsoluteValue( ) functions.
 *      110111    J. Melman         Added computeModulo( ) function.
 *      110202    K. Kumar          Added overload for State* for computeLinearInterpolation( ).
 *      110411    K. Kumar          Added convertCartesianToSpherical( ) function.
 *      110707    K. Kumar          Added computeSampleMean( ), computeSampleVariance( ) functions.
 *      110810    J. Leloux         Corrected doxygen documentation (equations).
 *      110824    J. Leloux         Corrected doxygen documentation.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120202    K. Kumar          Moved linear interpolation functions into new Interpolators
 *                                  sub-directory.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *      Spiegel, M.R., Stephens, L.J. Statistics, Fourth Edition, Schaum's
 *          Outlines, McGraw-Hill, 2008.
 *
 */

#ifndef TUDAT_LINEAR_INTERPOLATION_H
#define TUDAT_LINEAR_INTERPOLATION_H

#include <Eigen/Core>
#include <map>
#include "Tudat/Astrodynamics/States/state.h"

namespace tudat
{
namespace mathematics
{
namespace interpolators
{

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
 astrodynamics::states::State computeLinearInterpolation(
        std::map < double, astrodynamics::states::State >& sortedIndepedentAndDependentVariables,
        double targetIndependentVariableValue );

} // namespace interpolators
} // namespace mathematics
} // namespace tudat

#endif // TUDAT_LINEAR_INTERPOLATION_H
