/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_NONLINEAR_LEAST_SQUARES_ESTIMATION_H
#define TUDAT_NONLINEAR_LEAST_SQUARES_ESTIMATION_H

#include <map>

#include <Eigen/Core>
#include <boost/function.hpp>

namespace tudat
{

namespace linear_algebra
{

//! Function to perform a non-linear least squares estimation.
/*!
 *  Function to perform a non-linear least squares estimation. The non-linear least squares method is an iterative
 *  process, which uses the information from the actual and estimated observations, to estimate the model parameters, with
 *  the aid of a design matrix. The initial estimate of the model parameters is updated every iteration with the result of the
 *  least squares equation. The iterative process is halted whenever the norm of the update is below the user-provided
 *  threshold or when the maximum number of iterations is reached.
 *  \param observationAndJacobianFunctions Function returning a pair of expected observations and Jacobian of the
 *      observation function w.r.t. the model parameters (i.e., the design matrix), where the input is the current estimate
 *      of the model parameters.
 *  \param initialEstimate Initial estimate of the model parameters.
 *  \param actualObservations Vector containing the actual observations that need to be fitted by the model.
 *  \param convergenceTolerance Double denoting the convergence criterion for the norm of the update vector.
 *  \param maximumNumberOfIterations Integer denoting the maximum number of iterations.
 *  \return Optimal value of the model parameters that minimize the least squares error between expected and actual observations.
 */
Eigen::VectorXd nonLinearLeastSquaresFit(
        const boost::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& observationAndJacobianFunctions,
        const Eigen::VectorXd& initialEstimate, const Eigen::VectorXd& actualObservations,
        const double convergenceTolerance = 1.0e-8, const unsigned int maximumNumberOfIterations = 10 );

} // namespace linear_algebra

} // namespace tudat

#endif // TUDAT_NONLINEAR_LEAST_SQUARES_ESTIMATION_H
