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

#ifndef TUDAT_LEASTSQUARESESTIMATION_H
#define TUDAT_LEASTSQUARESESTIMATION_H

#include <functional>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include <Eigen/Core>
#include <Eigen/SVD>

#include "leastSquaresTraits.h"

namespace tudat
{

namespace linear_algebra
{

//! Function to get condition number of matrix from SVD decomposition
/*!
 *  Function to get condition number of matrix from SVD decomposition
 * \param singularValueDecomposition SVD decomposition of matrix
 * \tparam Options Options for the underlying Eigen JacobiSVD decomposition
 * \return Condition number of matrix
 */
template< typename M, int Options >
typename from_eigen< M >::value_type getConditionNumberOfDecomposedMatrix(
        const Eigen::JacobiSVD< M, Options >& singularValueDecomposition )
{
    typename from_eigen< M >::dense_vector_type singularValues = singularValueDecomposition.singularValues( );
    return singularValues( 0 ) / singularValues( singularValues.rows( ) - 1 );
}

//! Function to get condition number of matrix (using SVD decomposition)
/*!
 *  Function to get condition number of matrix (using SVD decomposition)
 * \param designMatrix Matrix for which condition number is to be computed
 * \return Condition number of matrix
 */
template< typename M >
typename from_eigen< M >::value_type getConditionNumberOfDesignMatrix( const M& designMatrix )
{
    return getConditionNumberOfDecomposedMatrix( designMatrix.template jacobiSvd< Eigen::ComputeThinU | Eigen::ComputeFullV >( ) );
}

//! Solve system of equations, checking condition number in the process where the solver supports it
/*!
 * Solve system of equations A*x = b for the vector x. The solver is selected by the storage type of the matrix
 * (SVD for dense matrices, sparse QR for sparse matrices), see `SolverTraits`.
 * \param matrixToInvert Matrix A that is to be inverted to solve the equation
 * \param rightHandSideVector Vector on the righthandside of the matrix equation that is to be solved
 * \param limitConditionNumberForWarning Maximum value of the condition number of the covariance matrix that is allowed
 * (warning printed when exceeded, dense matrices only)
 * \return Solution x of matrix equation A*x=b
 */
template< typename M >
typename from_eigen< M >::dense_vector_type solveSystemOfEquations(
        const M& matrixToInvert,
        const typename from_eigen< M >::dense_vector_type& rightHandSideVector,
        typename from_eigen< M >::value_type limitConditionNumberForWarning = 1.0E8 )
{
    return SolverTraits< M, typename from_eigen< M >::solver >::solve(
            matrixToInvert, rightHandSideVector, limitConditionNumberForWarning );
}

//! Solve system of equations with SVD decomposition, checking condition number in the process
/*!
 * Solve system of equations with SVD decomposition, checking condition number in the process. This function solves
 * A*x = b for the vector x. Dense matrices only; use `solveSystemOfEquations` to select the solver from the storage type.
 * \param matrixToInvert Matrix A that is to be inverted to solve the equation
 * \param rightHandSideVector Vector on the righthandside of the matrix equation that is to be solved
 * \param limitConditionNumberForWarning Maximum value of the condition number of the covariance matrix that is allowed
 * (warning printed when exceeded)
 * \return Solution x of matrix equation A*x=b
 */
template< typename M >
typename from_eigen< M >::dense_vector_type solveSystemOfEquationsWithSvd(
        const M& matrixToInvert,
        const typename from_eigen< M >::dense_vector_type& rightHandSideVector,
        typename from_eigen< M >::value_type limitConditionNumberForWarning = 1.0E8 )
{
    return SolverTraits< M, SVD >::solve( matrixToInvert, rightHandSideVector, limitConditionNumberForWarning );
}

//! Function to multiply information matrix by diagonal weights matrix
/*!
 * Function to multiply information matrix by diagonal weights matrix
 * \param designMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \return designMatrix, premultiplied by square matrix with diagonalOfWeightMatrix as diagonal elements
 */
template< typename M >
M multiplyDesignMatrixByDiagonalWeightMatrix( const M& designMatrix,
                                              const typename from_eigen< M >::dense_vector_type& diagonalOfWeightMatrix )
{
    return diagonalOfWeightMatrix.asDiagonal( ) * designMatrix;
}

// TODO: ask dominic about replacing arguments with std::optional.
//! Function to compute inverse of covariance matrix at current iteration, including influence of a priori information
/*!
 * Function to compute inverse of covariance matrix at current iteration, including influence of a priori information
 * \param designMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \param inverseOfAPrioriCovarianceMatrix Inverse of a priori covariance matrix
 * (warning printed when exceeded)
 * \return Inverse of covariance matrix at current iteration
 */
template< typename M >
M calculateInverseOfUpdatedCovarianceMatrix(
        const M& designMatrix,
        const typename from_eigen< M >::dense_vector_type& diagonalOfWeightMatrix,
        const M& inverseOfAPrioriCovarianceMatrix,
        const std::optional< M >& constraintMultiplier = std::nullopt,
        const std::optional< typename from_eigen< M >::dense_vector_type >& constraintRightHandside = std::nullopt,
        typename from_eigen< M >::value_type limitConditionNumberForWarning = 1.0E8 )
{
    using traits = typename from_eigen< M >::traits;

    // Information matrix, including a priori information
    M inverseOfCovarianceMatrix =
            inverseOfAPrioriCovarianceMatrix + designMatrix.transpose( ) * diagonalOfWeightMatrix.asDiagonal( ) * designMatrix;

    // Add constraints to inverse covariance matrix if required
    if( constraintMultiplier )
    {
        if( !constraintRightHandside )
        {
            throw std::runtime_error(
                    "Error when performing constrained least-squares, both multiplier and right-hand-side should be given." );
        }

        if( constraintMultiplier->rows( ) != constraintRightHandside->rows( ) )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible" );
        }

        if( constraintMultiplier->cols( ) != designMatrix.cols( ) )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible with partials" );
        }

        traits::extend_symmetric( inverseOfCovarianceMatrix, *constraintMultiplier );
    }

    return inverseOfCovarianceMatrix;
}

//! Function to compute the contribution of consider parameters to the covariance matrix
template< typename M >
M calculateConsiderParametersCovarianceContribution( const M& normalisedCovarianceMatrix,
                                                     const M& designMatrix,
                                                     const typename from_eigen< M >::dense_vector_type& diagonalOfWeightMatrix,
                                                     const M& considerDesignMatrix,
                                                     const M& considerCovariance )
{
    M covarianceTimesWeightedPartials =
            normalisedCovarianceMatrix * multiplyDesignMatrixByDiagonalWeightMatrix( designMatrix, diagonalOfWeightMatrix ).transpose( );
    return ( covarianceTimesWeightedPartials * considerDesignMatrix ) * considerCovariance *
            ( considerDesignMatrix.transpose( ) * covarianceTimesWeightedPartials.transpose( ) );
}

//! Function to compute the contribution of consider parameters to the covariance matrix, for a sparse design matrix
/*!
 * Same result as the dense version, (P H^T W H_c) C (P H^T W H_c)^T with P the normalised covariance, but with the products
 * ordered such that no dense matrix of the size of the design matrix is formed. The partials w.r.t. the consider parameters
 * (H_c) are dense.
 */
inline Eigen::MatrixXd calculateConsiderParametersCovarianceContribution( const Eigen::MatrixXd& normalisedCovarianceMatrix,
                                                                          const Eigen::SparseMatrix< double >& designMatrix,
                                                                          const Eigen::VectorXd& diagonalOfWeightMatrix,
                                                                          const Eigen::MatrixXd& considerDesignMatrix,
                                                                          const Eigen::MatrixXd& considerCovariance )
{
    Eigen::MatrixXd covarianceTimesConsiderNormalMatrix =
            normalisedCovarianceMatrix * ( designMatrix.transpose( ) * ( diagonalOfWeightMatrix.asDiagonal( ) * considerDesignMatrix ) );
    return covarianceTimesConsiderNormalMatrix * considerCovariance * covarianceTimesConsiderNormalMatrix.transpose( );
}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals and a priori
//! information
/*!
 * Function to perform an iteration least squares estimation from information matrix, weights and residuals and a priori
 * information, as is typically done in orbit determination. This function also takes an inverse if the a priori covariance
 * matrix to constrain/stabilize the inversion.
 * \param designMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param observationResiduals Difference between measured and simulated observations
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \param inverseOfAPrioriCovarianceMatrix Inverse of a priori covariance matrix
 * (warning printed when exceeded)
 * \param limitConditionNumberForWarning Maximum value of the condition number of the covariance matrix that is allowed
 * \param constraintMultiplier Multiplier for estimated parameter that defines linear constraint
 * \param constraintRightHandside Right-hand side estimation linear constraint
 * \param aprioriParameterDeviation Deviation of the current parameter estimate from the a priori parameter vector, expressed in
 * the same normalized coordinates as the parameter correction
 * \return Pair containing: (first: parameter adjustment, second: inverse covariance)
 */
template< typename M >
std::pair< typename from_eigen< M >::dense_vector_type, M > performLeastSquaresAdjustmentFromDesignMatrix(
        const M& designMatrix,
        const typename from_eigen< M >::dense_vector_type& observationResiduals,
        const typename from_eigen< M >::dense_vector_type& diagonalOfWeightMatrix,
        const M& inverseOfAPrioriCovarianceMatrix,
        typename from_eigen< M >::value_type limitConditionNumberForWarning = 1.0E8,
        const M& constraintMultiplier = M( 0, 0 ),
        const typename from_eigen< M >::dense_vector_type& constraintRightHandside = typename from_eigen< M >::dense_vector_type( 0 ),
        const M& designMatrixConsiderParameters = M( 0, 0 ),
        const typename from_eigen< M >::dense_vector_type& considerParametersDeviations = typename from_eigen< M >::dense_vector_type( 0 ),
        const M& additionalNormalMatrix = M( 0, 0 ),
        const typename from_eigen< M >::dense_vector_type& additionalRightHandSide = typename from_eigen< M >::dense_vector_type( 0 ),
        const typename from_eigen< M >::dense_vector_type& aprioriParameterDeviation = typename from_eigen< M >::dense_vector_type( 0 ) )
{
    using traits = typename from_eigen< M >::traits;
    using dense_vector_type = typename from_eigen< M >::dense_vector_type;

    dense_vector_type rightHandSide;
    if( considerParametersDeviations.size( ) > 0 && designMatrixConsiderParameters.size( ) > 0 )
    {
        rightHandSide = designMatrix.transpose( ) *
                ( diagonalOfWeightMatrix.cwiseProduct( observationResiduals +
                                                       designMatrixConsiderParameters * considerParametersDeviations ) );
    }
    else
    {
        rightHandSide = designMatrix.transpose( ) * ( diagonalOfWeightMatrix.cwiseProduct( observationResiduals ) );
    }

    const int nParams = static_cast< int >( designMatrix.cols( ) );
    if( inverseOfAPrioriCovarianceMatrix.rows( ) != nParams || inverseOfAPrioriCovarianceMatrix.cols( ) != nParams )
    {
        throw std::runtime_error(
                "Error in performLeastSquaresAdjustmentFromDesignMatrix: inverseOfAPrioriCovarianceMatrix has dimensions " +
                std::to_string( inverseOfAPrioriCovarianceMatrix.rows( ) ) + "x" +
                std::to_string( inverseOfAPrioriCovarianceMatrix.cols( ) ) + ", expected " + std::to_string( nParams ) + "x" +
                std::to_string( nParams ) + "." );
    }
    if( aprioriParameterDeviation.size( ) > 0 )
    {
        if( aprioriParameterDeviation.size( ) != nParams )
        {
            throw std::runtime_error( "Error in performLeastSquaresAdjustmentFromDesignMatrix: aprioriParameterDeviation has size " +
                                      std::to_string( aprioriParameterDeviation.size( ) ) + ", expected " + std::to_string( nParams ) +
                                      "." );
        }
        // The prior residual is the a priori parameter vector minus the current estimate, i.e. the negative deviation.
        rightHandSide -= inverseOfAPrioriCovarianceMatrix * aprioriParameterDeviation;
    }

    M inverseOfCovarianceMatrix;
    if( constraintMultiplier.rows( ) != 0 )
    {
        inverseOfCovarianceMatrix = calculateInverseOfUpdatedCovarianceMatrix< M >(
                designMatrix, diagonalOfWeightMatrix, inverseOfAPrioriCovarianceMatrix, constraintMultiplier, constraintRightHandside );

        int numberOfConstraints = constraintMultiplier.rows( );
        int numberOfConstrainedParameters = constraintMultiplier.cols( );

        rightHandSide.conservativeResize( numberOfConstrainedParameters + numberOfConstraints );
        rightHandSide.segment( numberOfConstrainedParameters, numberOfConstraints ) = constraintRightHandside;
    }
    else
    {
        inverseOfCovarianceMatrix =
                calculateInverseOfUpdatedCovarianceMatrix< M >( designMatrix, diagonalOfWeightMatrix, inverseOfAPrioriCovarianceMatrix );
    }

    // Soft-constraint additions (e.g. inter-arc continuity): inject into the parameter block of the LHS/RHS,
    // leaving any Lagrange-multiplier rows from the hard equality constraint above untouched.
    if( additionalNormalMatrix.size( ) > 0 )
    {
        if( additionalNormalMatrix.rows( ) != nParams || additionalNormalMatrix.cols( ) != nParams )
        {
            throw std::runtime_error( "Error in performLeastSquaresAdjustmentFromDesignMatrix: additionalNormalMatrix has shape " +
                                      std::to_string( additionalNormalMatrix.rows( ) ) + "x" +
                                      std::to_string( additionalNormalMatrix.cols( ) ) + ", expected " + std::to_string( nParams ) + "x" +
                                      std::to_string( nParams ) + "." );
        }
        traits::add_to_top_left( inverseOfCovarianceMatrix, additionalNormalMatrix );
    }
    if( additionalRightHandSide.size( ) > 0 )
    {
        if( additionalRightHandSide.size( ) != nParams )
        {
            throw std::runtime_error( "Error in performLeastSquaresAdjustmentFromDesignMatrix: additionalRightHandSide has size " +
                                      std::to_string( additionalRightHandSide.size( ) ) + ", expected " + std::to_string( nParams ) + "." );
        }
        rightHandSide.head( nParams ) += additionalRightHandSide;
    }

    return std::make_pair( solveSystemOfEquations< M >( inverseOfCovarianceMatrix, rightHandSide, limitConditionNumberForWarning ),
                           inverseOfCovarianceMatrix );
}

//! Function to perform an iteration least squares estimation, for dense matrices
/*!
 * Dense version of the function above, which also accepts arguments that are convertible to dense dynamic-size matrices and
 * vectors (e.g. fixed-size matrices or expressions), for which the template arguments of the general version cannot be deduced.
 */
inline std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const double limitConditionNumberForWarning = 1.0E8,
        const Eigen::MatrixXd& constraintMultiplier = Eigen::MatrixXd( 0, 0 ),
        const Eigen::VectorXd& constraintRightHandside = Eigen::VectorXd( 0 ),
        const Eigen::MatrixXd& designMatrixConsiderParameters = Eigen::MatrixXd( 0, 0 ),
        const Eigen::VectorXd& considerParametersDeviations = Eigen::VectorXd( 0 ),
        const Eigen::MatrixXd& additionalNormalMatrix = Eigen::MatrixXd( 0, 0 ),
        const Eigen::VectorXd& additionalRightHandSide = Eigen::VectorXd( 0 ),
        const Eigen::VectorXd& aprioriParameterDeviation = Eigen::VectorXd( 0 ) )
{
    return performLeastSquaresAdjustmentFromDesignMatrix< Eigen::MatrixXd >( designMatrix,
                                                                             observationResiduals,
                                                                             diagonalOfWeightMatrix,
                                                                             inverseOfAPrioriCovarianceMatrix,
                                                                             limitConditionNumberForWarning,
                                                                             constraintMultiplier,
                                                                             constraintRightHandside,
                                                                             designMatrixConsiderParameters,
                                                                             considerParametersDeviations,
                                                                             additionalNormalMatrix,
                                                                             additionalRightHandSide,
                                                                             aprioriParameterDeviation );
}

//! Function to perform an iteration of least squares estimation from information matrix, weights and residuals
/*!
 * Function to perform an iteration of least squares estimation from information matrix, weights and residuals, as is
 * typically done in orbit determination
 * \param designMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param observationResiduals Difference between measured and simulated observations
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \param limitConditionNumberForWarning Maximum value of the condition number of the covariance matrix that is allowed
 * (warning printed when exceeded)
 * \return Pair containing: (first: parameter adjustment, second: inverse covariance)
 */
template< typename M >
std::pair< typename from_eigen< M >::dense_vector_type, M > performLeastSquaresAdjustmentFromDesignMatrix(
        const M& designMatrix,
        const typename from_eigen< M >::dense_vector_type& observationResiduals,
        const typename from_eigen< M >::dense_vector_type& diagonalOfWeightMatrix,
        typename from_eigen< M >::value_type limitConditionNumberForWarning = 1.0E8 )
{
    return performLeastSquaresAdjustmentFromDesignMatrix< M >( designMatrix,
                                                               observationResiduals,
                                                               diagonalOfWeightMatrix,
                                                               from_eigen< M >::traits::zeros( designMatrix.cols( ), designMatrix.cols( ) ),
                                                               limitConditionNumberForWarning );
}

//! Function to perform an iteration of least squares estimation from information matrix and residuals
/*!
 * Function to perform an iteration of least squares estimation from information matrix and residuals, with all weights
 * fixed to 1.0.
 * \param designMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param observationResiduals Difference between measured and simulated observations
 * \param limitConditionNumberForWarning Maximum value of the condition number of the covariance matrix that is allowed
 * (warning printed when exceeded)
 * \return Pair containing: (first: parameter adjustment, second: inverse covariance)
 */
template< typename M >
std::pair< typename from_eigen< M >::dense_vector_type, M > performLeastSquaresAdjustmentFromDesignMatrix(
        const M& designMatrix,
        const typename from_eigen< M >::dense_vector_type& observationResiduals,
        typename from_eigen< M >::value_type limitConditionNumberForWarning = 1.0E8 )
{
    return performLeastSquaresAdjustmentFromDesignMatrix< M >(
            designMatrix,
            observationResiduals,
            from_eigen< M >::dense_vector_type::Constant( observationResiduals.size( ), typename from_eigen< M >::value_type( 1 ) ),
            limitConditionNumberForWarning );
}

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
        const Eigen::VectorXd& initialEstimate,
        const Eigen::VectorXd& actualObservations,
        const double initialScaling = 1.0e-6,
        const double convergenceTolerance = 1.0e-8,
        const unsigned int maximumNumberOfIterations = 25 );

//! Perform a non-linear least-squares estimation and return iteration diagnostics by reference.
/*!
 * This overload uses the same Levenberg-Marquardt implementation as nonLinearLeastSquaresFit, while additionally
 * returning the number of iterations and whether the update-norm convergence criterion was met.
 */
Eigen::VectorXd nonLinearLeastSquaresFit(
        const std::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& observationAndJacobianFunctions,
        const Eigen::VectorXd& initialEstimate,
        const Eigen::VectorXd& actualObservations,
        const double initialScaling,
        const double convergenceTolerance,
        const unsigned int maximumNumberOfIterations,
        unsigned int& numberOfIterations,
        bool& converged );

}  // namespace linear_algebra

}  // namespace tudat

#endif  // TUDAT_LEASTSQUARESESTIMATION_H
