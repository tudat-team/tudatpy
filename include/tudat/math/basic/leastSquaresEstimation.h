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

#include <iostream>

#include <Eigen/Core>
#include <Eigen/SVD>

#include "leastSquaresTraits.h"

namespace tudat
{

namespace linear_algebra
{

//! Function to get condition number of matrix (using SVD decomposition)
/*!
 *  Function to get condition number of matrix (using SVD decomposition)
 * \param designMatrix Matrix for which condition number is to be computed
 * \return Condition number of matrix
 */
template <typename M>
typename from_eigen<M>::value_type getConditionNumberOfDesignMatrix( const M &designMatrix ) {
    return getConditionNumberOfDesignMatrix(designMatrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeFullV));
}

// TODO: This function is only used internally, so it shouldn't be in the header file.
//! Function to get condition number of matrix from SVD decomposition
/*!
 *  Function to get condition number of matrix from SVD decomposition
 * \param singularValueDecomposition SVD decomposition of matrix
 * \tparam Options Options for the underlying Eigen JacobiSVD decomposition
 * \return Condition number of matrix
 */
template <typename M, int Options>
typename from_eigen<M>::value_type getConditionNumberOfDecomposedMatrix(
    const Eigen::JacobiSVD< M, Options >& singularValueDecomposition )
{
    typename from_eigen<M>::dense_vector_type singularValues = singularValueDecomposition.singularValues( );
    return singularValues( 0 ) / singularValues( singularValues.rows( ) - 1 );
}

//! Solve system of equations with SVD decomposition, checking condition number in the process
/*!
 * Solve system of equations with SVD decomposition, checking condition number in the process. This function solves
 * A*x = b for the vector x.
 * \param matrixToInvert Matrix A that is to be inverted to solve the equation
 * \param rightHandSideVector Vector on the righthandside of the matrix equation that is to be solved
 * \param limitConditionNumberForWarning Maximum value of the condition number of the covariance matrix that is allowed
 * (warning printed when exceeded)
 * \return Solution x of matrix equation A*x=b
 */
template <typename M>
typename from_eigen<M>::dense_vector_type solveSystemOfEquationsWithSvd(
    const M& matrixToInvert,
    const typename from_eigen<M>::dense_vector_type& rightHandSideVector,
    typename from_eigen<M>::value_type limitConditionNumberForWarning = 1.0E8 )
{
    auto svdDecomposition = matrixToInvert.jacobiSvd( Eigen::ComputeThinU | Eigen::ComputeThinV );
    if( limitConditionNumberForWarning == limitConditionNumberForWarning )
    {
        double conditionNumber = getConditionNumberOfDecomposedMatrix( svdDecomposition );

        if( conditionNumber > limitConditionNumberForWarning )
        {
            std::cerr << "Warning when performing least squares, condition number is " << conditionNumber << std::endl;
        }
    }
    return svdDecomposition.solve( rightHandSideVector );
}

//! Function to multiply information matrix by diagonal weights matrix
/*!
 * Function to multiply information matrix by diagonal weights matrix
 * \param designMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \return designMatrix, premultiplied by square matrix with diagonalOfWeightMatrix as diagonal elements
 */
template <typename M>
M multiplyDesignMatrixByDiagonalWeightMatrix( const M& designMatrix,
                                              const typename from_eigen<M>::dense_vector_type& diagonalOfWeightMatrix )
{
    return diagonalOfWeightMatrix.asDiagonal() * designMatrix;
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
template <typename M>
M calculateInverseOfUpdatedCovarianceMatrix( const M& designMatrix,
                                             const typename from_eigen<M>::dense_vector_type& diagonalOfWeightMatrix,
                                             const M& inverseOfAPrioriCovarianceMatrix,
                                             const std::optional<M> &constraintMultiplier = std::nullopt,
                                             const std::optional<typename from_eigen<M>::dense_vector_type> &constraintRightHandside =
                                                std::nullopt,
                                             typename from_eigen<M>::value_type limitConditionNumberForWarning = 1.0E8 )
{
    // Add constraints to inverse covariance matrix if required
    Eigen::MatrixXd inverseOfCovarianceMatrix = inverseOfAPrioriCovarianceMatrix +
            designMatrix.transpose( ) * diagonalOfWeightMatrix.asDiagonal( ) * designMatrix;

    if ( constraintMultiplier )
    {
        if ( !constraintRightHandside ) {
            throw std::runtime_error( "Error when performing constrained least-squares, both multiplier and right-hand-side should be given." );
        }

        if( constraintMultiplier->rows( ) != constraintRightHandside->rows( ) )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible" );
        }

        if( constraintMultiplier->cols( ) != designMatrix.cols( ) )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible with partials" );
        }

        int numberOfConstraints = constraintMultiplier->rows( );
        int numberOfParameters = constraintMultiplier->cols( );

        // TODO: Figure out what this does and replace with sparse compatible operations
        inverseOfCovarianceMatrix.conservativeResize( numberOfParameters + numberOfConstraints, numberOfParameters + numberOfConstraints );
        inverseOfCovarianceMatrix.block( numberOfParameters, 0, numberOfConstraints, numberOfParameters ) = *constraintMultiplier;
        inverseOfCovarianceMatrix.block( 0, numberOfParameters, numberOfParameters, numberOfConstraints ) =
                constraintMultiplier->transpose( );
        inverseOfCovarianceMatrix.block( numberOfParameters, numberOfParameters, numberOfConstraints, numberOfConstraints ).setZero( );
    }

    return inverseOfCovarianceMatrix;
}

//! Function to compute inverse of covariance matrix at current iteration
/*!
 * Function to compute inverse of covariance matrix at current iteration
 * \param designMatrix Matrix containing partial derivatives of observations (rows) w.r.t. estimated parameters
 * (columns)
 * \param diagonalOfWeightMatrix Diagonal of observation weights matrix (assumes all weights to be uncorrelated)
 * \return Inverse of covariance matrix at current iteration
 */
template <typename M>
M calculateInverseOfUpdatedCovarianceMatrix( const M& designMatrix,
                                             const typename from_eigen<M>::dense_vector& diagonalOfWeightMatrix,
                                             typename from_eigen<M>::value_type limitConditionNumberForWarning = 1.0E8 );

template <typename M>
M calculateConsiderParametersCovarianceContribution( const M& normalisedCovarianceMatrix,
                                                     const M& designMatrix,
                                                     const typename from_eigen<M>::dense_vector_type& diagonalOfWeightMatrix,
                                                     const M& considerDesignMatrix,
                                                     const M& considerCovariance )
{
    M covarianceTimesWeightedPartials =
            normalisedCovarianceMatrix * multiplyDesignMatrixByDiagonalWeightMatrix( designMatrix, diagonalOfWeightMatrix ).transpose( );
    return ( covarianceTimesWeightedPartials * considerDesignMatrix ) * considerCovariance *
            ( considerDesignMatrix.transpose( ) * covarianceTimesWeightedPartials.transpose( ) );
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
 * \return Pair containing: (first: parameter adjustment, second: inverse covariance)
 */
template <typename M>
std::pair< typename from_eigen<M>::dense_vector_type, M > performLeastSquaresAdjustmentFromDesignMatrix(
        const M& designMatrix,
        const typename from_eigen<M>::dense_vector_type& observationResiduals,
        const typename from_eigen<M>::dense_vector_type& diagonalOfWeightMatrix,
        const M& inverseOfAPrioriCovarianceMatrix,
        typename from_eigen<M>::value_type limitConditionNumberForWarning = 1.0E8,
        const M& constraintMultiplier = M( 0, 0 ),
        const typename from_eigen<M>::dense_vector_type& constraintRightHandside =
            typename from_eigen<M>::dense_vector_type(0),
        const M& designMatrixConsiderParameters = M( 0, 0 ),
        const typename from_eigen<M>::dense_vector_type& considerParametersDeviations =
            typename from_eigen<M>::dense_vector_type(0))
{
    typename from_eigen<M>::dense_vector_type rightHandSide =
            from_eigen<M>::dense_vector_type::Zero( observationResiduals.size( ) );
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

    M inverseOfCovarianceMatrix;
    if( constraintMultiplier.rows( ) != 0 )
    {
        inverseOfCovarianceMatrix = calculateInverseOfUpdatedCovarianceMatrix< M >(
                designMatrix,
                diagonalOfWeightMatrix,
                inverseOfAPrioriCovarianceMatrix,
                constraintMultiplier,
                constraintRightHandside );

        int numberOfConstraints = constraintMultiplier.rows( );
        int numberOfParameters = constraintMultiplier.cols( );

        rightHandSide.conservativeResize( numberOfParameters + numberOfConstraints );
        rightHandSide.segment( numberOfParameters, numberOfConstraints ) = constraintRightHandside;
    }
    else
    {
        inverseOfCovarianceMatrix = calculateInverseOfUpdatedCovarianceMatrix< M >(
                designMatrix,
                diagonalOfWeightMatrix,
                inverseOfAPrioriCovarianceMatrix );
    }

    return std::make_pair( solveSystemOfEquationsWithSvd( inverseOfCovarianceMatrix, rightHandSide, limitConditionNumberForWarning ),
                           inverseOfCovarianceMatrix );
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
template <typename M>
std::pair< typename from_eigen<M>::dense_vector_type, M > performLeastSquaresAdjustmentFromDesignMatrix(
        const M& designMatrix,
        const M& observationResiduals,
        const typename from_eigen<M>::dense_vector_type& diagonalOfWeightMatrix,
        typename from_eigen<M>::value_type limitConditionNumberForWarning = 1.0E8 );

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
template <typename M>
std::pair< typename from_eigen<M>::dense_vector_type, M > performLeastSquaresAdjustmentFromDesignMatrix(
        const M& designMatrix,
        const typename from_eigen<M>::dense_vector_type& observationResiduals,
        typename from_eigen<M>::value_type limitConditionNumberForWarning = 1.0E8 );

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
template <typename Fn, typename Vec>
Vec nonLinearLeastSquaresFit(
        // std::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >&
        Fn observationAndJacobianFunctions,
        const Vec& initialEstimate,
        const Vec& actualObservations,
        const double initialScaling = 1.0e-6,
        const double convergenceTolerance = 1.0e-8,
        const unsigned int maximumNumberOfIterations = 25 );

}  // namespace linear_algebra

}  // namespace tudat

#endif  // TUDAT_LEASTSQUARESESTIMATION_H
