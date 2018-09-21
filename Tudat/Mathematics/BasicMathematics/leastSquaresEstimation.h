/*    Copyright (c) 2010-2018, Delft University of Technology
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

#ifndef TUDAT_LEASTSQUARESESTIMATION_H
#define TUDAT_LEASTSQUARESESTIMATION_H

#include <map>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <boost/function.hpp>

namespace tudat
{

namespace linear_algebra
{

//! Function to get condition number of matrix (using SVD decomposition)
/*!
 *  Function to get condition number of matrix (using SVD decomposition)
 * \param informationMatrix Matrix for which condition number is to be computed
 * \return Condition number of matrix
 */
double getConditionNumberOfInformationMatrix( const Eigen::MatrixXd informationMatrix );

//! Function to get condition number of matrix from SVD decomposition
/*!
 *  Function to get condition number of matrix from SVD decomposition
 * \param singularValueDecomposition SVD decomposition of matrix
 * \return Condition number of matrix
 */
double getConditionNumberOfDecomposedMatrix( const Eigen::JacobiSVD< Eigen::MatrixXd >& singularValueDecomposition );

//! Solve system of equations with SVD decomposition, checking condition number in the process
/*!
 * Solve system of equations with SVD decomposition, checking condition number in the process. This function solves
 * A*x = b for the vector x.
 * \param matrixToInvert Matrix A that is to be inverted to solve the equation
 * \param rightHandSideVector Vector on the righthandside of the matrix equation that is to be solved
 * \param checkConditionNumber Boolean to denote whether the condition number is checked when estimating (warning is printed
 * when value exceeds maximumAllowedConditionNumber)
 * \param maximumAllowedConditionNumber Maximum value of the condition number of the covariance matrix that is allowed
 * (warning printed when exceeded)
 * \return Solution x of matrix equation A*x=b
 */
Eigen::VectorXd solveSystemOfEquationsWithSvd( const Eigen::MatrixXd matrixToInvert,
                                               const Eigen::VectorXd rightHandSideVector,
                                               const bool checkConditionNumber = 1,
                                               const double maximumAllowedConditionNumber = 1.0E-8 );

//! Function to multiply information matrix by diagonal weights matrix
/*!
 * Function to multiply information matrix by diagonal weights matrix
 * \param informationMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \return informationMatrix, premultiplied by square matrix with diagonalOfWeightMatrix as diagonal elements
 */
Eigen::MatrixXd multiplyInformationMatrixByDiagonalWeightMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix );

//! Function to compute inverse of covariance matrix at current iteration, including influence of a priori information
/*!
 * Function to compute inverse of covariance matrix at current iteration, including influence of a priori information
 * \param informationMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \param inverseOfAPrioriCovarianceMatrix Inverse of a priori covariance matrix
 * (warning printed when exceeded)
 * \return Inverse of covariance matrix at current iteration
 */
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix );

//! Function to compute inverse of covariance matrix at current iteration
/*!
 * Function to compute inverse of covariance matrix at current iteration
 * \param informationMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \return Inverse of covariance matrix at current iteration
 */
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix );

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals and a priori
//! information
/*!
 * Function to perform an iteration least squares estimation from information matrix, weights and residuals and a priori
 * information, as is typically done in orbit determination. This function also takes an inverse if the a priori covariance
 * matrix to constrain/stabilize the inversion.
 * \param informationMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param observationResiduals Difference between measured and simulated observations
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \param inverseOfAPrioriCovarianceMatrix Inverse of a priori covariance matrix
 * (warning printed when exceeded)
 * \param checkConditionNumber Boolean to denote whether the condition number is checked when estimating (warning is printed
 * when value exceeds maximumAllowedConditionNumber)
 * \param maximumAllowedConditionNumber Maximum value of the condition number of the covariance matrix that is allowed
 * \param constraintMultiplier Multiplier for estimated parameter that defines linear constraint
 * \param constraintRightHandside Right-hand side estimation linear constraint
 * \return Pair containing: (first: parameter adjustment, second: inverse covariance)
 */
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const bool checkConditionNumber = 1,
        const double maximumAllowedConditionNumber = 1.0E8,
        const Eigen::MatrixXd& constraintMultiplier = Eigen::MatrixXd( 0, 0 ),
        const Eigen::VectorXd& constraintRightHandside = Eigen::VectorXd( 0 ) );

//! Function to perform an iteration of least squares estimation from information matrix, weights and residuals
/*!
 * Function to perform an iteration of least squares estimation from information matrix, weights and residuals, as is
 * typically done in orbit determination
 * \param informationMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param observationResiduals Difference between measured and simulated observations
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \param checkConditionNumber Boolean to denote whether the condition number is checked when estimating (warning is printed
 * when value exceeds maximumAllowedConditionNumber)
 * \param maximumAllowedConditionNumber Maximum value of the condition number of the covariance matrix that is allowed
 * (warning printed when exceeded)
 * \return Pair containing: (first: parameter adjustment, second: inverse covariance)
 */
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const bool checkConditionNumber = 1,
        const double maximumAllowedConditionNumber = 1.0E8 );

//! Function to perform an iteration of least squares estimation from information matrix and residuals
/*!
 * Function to perform an iteration of least squares estimation from information matrix and residuals, with all weights
 * fixed to 1.0.
 * \param informationMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param observationResiduals Difference between measured and simulated observations
 * \param checkConditionNumber Boolean to denote whether the condition number is checked when estimating (warning is printed
 * when value exceeds maximumAllowedConditionNumber)
 * \param maximumAllowedConditionNumber Maximum value of the condition number of the covariance matrix that is allowed
 * (warning printed when exceeded)
 * \return Pair containing: (first: parameter adjustment, second: inverse covariance)
 */
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const bool checkConditionNumber = 1,
        const double maximumAllowedConditionNumber = 1.0E8 );

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
Eigen::VectorXd getLeastSquaresPolynomialFit(
        const Eigen::VectorXd& independentValues,
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
std::vector< double > getLeastSquaresPolynomialFit(
        const std::map< double, double >& independentDependentValueMap,
        const std::vector< double >& polynomialPowers );

//! Function to perform a non-linear least squares estimation with the Levenberg-Marquardt method.
/*!
 *  Function to perform a non-linear least squares estimation. The non-linear least squares method is an iterative
 *  process, which uses the information from the actual and estimated observations, to estimate the model parameters, with
 *  the aid of a design matrix. The initial estimate of the model parameters is updated every iteration with the result of the
 *  least squares equation. The iterative process is halted whenever the norm of the update is below the user-provided
 *  threshold or when the maximum number of iterations is reached. The method used in this application is the Levenberg-Marquardt
 *  method, which uses a damping parameter \f$ \lambda \f$ to make the iterative process more stable and accurate.
 *  The reference for this implementation is (Madsen, K., et al.).
 *  \param observationAndJacobianFunctions Function returning a pair of expected observations and Jacobian of the
 *      observation function w.r.t. the model parameters (i.e., the design matrix), where the input is the current estimate
 *      of the model parameters.
 *  \param initialEstimate Initial estimate of the model parameters.
 *  \param actualObservations Vector containing the actual observations that need to be fitted by the model.
 *  \param initialScaling Double denoting the multiplicative factor to determine the damping parameter during the first iteration.
 *  \param convergenceTolerance Double denoting the convergence criterion for the norm of the update vector.
 *  \param maximumNumberOfIterations Integer denoting the maximum number of iterations.
 *  \return Optimal value of the model parameters that minimize the least squares error between expected and actual observations.
 */
Eigen::VectorXd nonLinearLeastSquaresFit(
        const std::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& observationAndJacobianFunctions,
        const Eigen::VectorXd& initialEstimate, const Eigen::VectorXd& actualObservations, const double initialScaling = 1.0e-6,
        const double convergenceTolerance = 1.0e-8, const unsigned int maximumNumberOfIterations = 25 );

} // namespace linear_algebra

} // namespace tudat

#endif // TUDAT_LEASTSQUARESESTIMATION_H
