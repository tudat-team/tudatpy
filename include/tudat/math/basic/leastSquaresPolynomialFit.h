/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Madsen, K., Nielsen, H., and Tingleff, O., Methods for Non-Linear Least Squares Problems, 2nd ed.,
 *          Technical University of Denmark, Faculty of Informatics and Mathematical Modelling, April 2004.
 */

/* These routines are currently only used in some unit tests. */

#ifndef TUDAT_LEASTSQUARESPOLYNOMIALFIT_H
#define TUDAT_LEASTSQUARESPOLYNOMIALFIT_H

#include <Eigen/Dense>
#include <vector>
#include <map>


namespace tudat {
namespace linear_algebra {

Eigen::VectorXd evaluatePolynomial( const Eigen::VectorXd& independentValues,
                                    const Eigen::VectorXd& polynomialCoefficients,
                                    const std::vector< double >& polynomialPowers );

//! Function to fit a univariate polynomial through a set of data
/*!
 *  Function to fit a univariate polynomial through a set of data. User must provide independent variables and observations
 *  (dependent variables), as well as a list of polynomial powers for which the coefficients are to be estimated.
 *  \param independentValues Independent variables of input data (e.g. time for observations as a function fo time). This
 *  variable becomes the polynomial argument.
 *  \param dependentValues Observations through which the polynomial is to be fitted, with entries defined at the
 *  corresponding entries of independentValues
 *  \param polynomialPowers List of powers of indepent variables for which coefficients are to be estimated.
 *  \return Coefficients of the polynomial powers, as estimated from the input data (in same order as polynomialPowers).
 */
Eigen::VectorXd getLeastSquaresPolynomialFit( const Eigen::VectorXd& independentValues,
                                              const Eigen::VectorXd& dependentValues,
                                              const std::vector< double >& polynomialPowers );

//! Function to fit a univariate polynomial through a set of data
/*!
 *  Function to fit a univariate polynomial through a set of data. User must provide independent variables and observations
 *  (dependent variables), as well as a list of polynomial powers for which the coefficients are to be estimated.
 *  \param independentDependentValueMap Map with key: independent variables of input data (e.g. time for observations as a
 *  function fo time), this variable becomes the polynomial argument. Map value: Observations through which the polynomial
 *  is to be fitted.
 *  \param polynomialPowers List of powers of indepent variables for which coefficients are to be estimated.
 *  \return Coefficients of the polynomial powers, as estimated from the input data (in same order as polynomialPowers).
 */
std::vector< double > getLeastSquaresPolynomialFit( const std::map< double, double >& independentDependentValueMap,
                                                    const std::vector< double >& polynomialPowers );

}} // namespace

#endif
