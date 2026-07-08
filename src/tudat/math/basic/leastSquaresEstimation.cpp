/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <cmath>
#include <iostream>

#include <Eigen/LU>
#include "tudat/math/basic/leastSquaresEstimation.h"

namespace tudat
{

namespace linear_algebra
{

//! Solve system of equations with SVD decomposition, checking condition number in the process
Eigen::VectorXd solveSystemOfEquationsWithSvd( const Eigen::MatrixXd matrixToInvert,
                                               const Eigen::VectorXd rightHandSideVector,
                                               const double limitConditionNumberForWarning )
{
    Eigen::JacobiSVD< Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV > svdDecomposition =
            matrixToInvert.jacobiSvd< Eigen::ComputeThinU | Eigen::ComputeThinV >( );
    if( !std::isnan( limitConditionNumberForWarning ) )
    {
        double conditionNumber = getConditionNumberOfDecomposedMatrix( svdDecomposition );

        if( conditionNumber > limitConditionNumberForWarning )
        {
            std::cerr << "Warning when performing least squares, condition number is " << conditionNumber << std::endl;
        }
    }
    return svdDecomposition.solve( rightHandSideVector );
}

//! Function to compute inverse of covariance matrix at current iteration
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix( const Eigen::MatrixXd& designMatrix,
                                                           const Eigen::VectorXd& diagonalOfWeightMatrix,
                                                           const double limitConditionNumberForWarning )
{
    return calculateInverseOfUpdatedCovarianceMatrix< Eigen::MatrixXd >(
            designMatrix, diagonalOfWeightMatrix, Eigen::MatrixXd::Zero( designMatrix.cols( ), designMatrix.cols( ) ) );
}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix( const Eigen::MatrixXd& designMatrix,
                                                                                             const Eigen::VectorXd& observationResiduals,
                                                                                             const Eigen::VectorXd& diagonalOfWeightMatrix,
                                                                                             const double limitConditionNumberForWarning )
{
    return performLeastSquaresAdjustmentFromDesignMatrix< Eigen::MatrixXd >( designMatrix,
                                                          observationResiduals,
                                                          diagonalOfWeightMatrix,
                                                          Eigen::MatrixXd::Zero( designMatrix.cols( ), designMatrix.cols( ) ),
                                                          limitConditionNumberForWarning );
}

//! Function to perform an iteration of least squares estimation from information matrix and residuals
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix( const Eigen::MatrixXd& designMatrix,
                                                                                             const Eigen::VectorXd& observationResiduals,
                                                                                             const double limitConditionNumberForWarning )
{
    return performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                          observationResiduals,
                                                          Eigen::VectorXd::Constant( observationResiduals.size( ), 1, 1.0 ),
                                                          limitConditionNumberForWarning );
}


//! Function to perform a non-linear least squares estimation with the Levenberg-Marquardt method.
Eigen::VectorXd nonLinearLeastSquaresFit(
        const std::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& observationAndJacobianFunctions,
        const Eigen::VectorXd& initialEstimate,
        const Eigen::VectorXd& actualObservations,
        const double initialScaling,
        const double convergenceTolerance,
        const unsigned int maximumNumberOfIterations )
{
    // Set current estimate to initial value
    Eigen::VectorXd currentEstimate = initialEstimate;

    // Initialize variables
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > pairOfEstimatedObservationsAndDesignMatrix;
    Eigen::MatrixXd designMatrix;
    Eigen::VectorXd offsetInObservations;
    Eigen::VectorXd updateInEstimate;

    // Initial parameters for Levenberg–Marquardt method
    double levenbergMarquardtDampingParameter = 0.0;
    double scalingParameterUpdate = 2.0;
    double levenbergMarquardtGainRatio = 0.0;

    // Start iterative loop
    unsigned int iteration = 0;
    do
    {
        // Compute current system and jacobian functions
        pairOfEstimatedObservationsAndDesignMatrix = observationAndJacobianFunctions( currentEstimate );
        designMatrix = pairOfEstimatedObservationsAndDesignMatrix.second;

        // Offset in observation
        offsetInObservations = actualObservations - pairOfEstimatedObservationsAndDesignMatrix.first;

        // Compute damping parameter for first iteration
        if( iteration == 0 )
        {
            levenbergMarquardtDampingParameter = initialScaling * ( designMatrix.transpose( ) * designMatrix ).diagonal( ).maxCoeff( );
        }

        // Compute update in estimate
        Eigen::VectorXd diagonalOfWeightMatrix = Eigen::VectorXd::Ones( offsetInObservations.rows( ) );
        Eigen::MatrixXd inverseOfAPrioriCovarianceMatrix = levenbergMarquardtDampingParameter *
                //                Eigen::MatrixXd( ( designMatrix.transpose( ) * designMatrix ).diagonal( ).asDiagonal( ) ); // Marquardt’s
                //                update
                Eigen::MatrixXd::Identity( currentEstimate.rows( ), currentEstimate.rows( ) );
        updateInEstimate = linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix(
                                   designMatrix, offsetInObservations, diagonalOfWeightMatrix, inverseOfAPrioriCovarianceMatrix, false )
                                   .first;

        // Check that update is real
        if( ( !updateInEstimate.allFinite( ) ) || ( updateInEstimate.hasNaN( ) ) )
        {
            throw std::runtime_error( "Error in non-linear least squares estimation. Iterative process diverges." );
        }

        // Compute gain ratio
        levenbergMarquardtGainRatio =
                ( offsetInObservations.squaredNorm( ) -
                  ( actualObservations - observationAndJacobianFunctions( currentEstimate + updateInEstimate ).first ).squaredNorm( ) ) /
                ( updateInEstimate.transpose( ) *
                  ( levenbergMarquardtDampingParameter * updateInEstimate + designMatrix.transpose( ) * offsetInObservations ) );

        // Update damping parameter
        if( levenbergMarquardtGainRatio > 0 )
        {
            // Reduce damping parameter, since good approximation
            levenbergMarquardtDampingParameter *= std::max( 1.0 / 3.0, 1.0 - std::pow( 2.0 * levenbergMarquardtGainRatio - 1.0, 3 ) );
            scalingParameterUpdate = 2;  // reset

            // Correct estimate
            currentEstimate += updateInEstimate;
        }
        else
        {
            // Increase damping parameter, since bad approximation and reject step
            levenbergMarquardtDampingParameter *= scalingParameterUpdate;
            scalingParameterUpdate *= 2.0;
        }

        // Increase iteration counter
        iteration++;
    } while( ( updateInEstimate.norm( ) > convergenceTolerance ) && ( iteration <= maximumNumberOfIterations ) );

    // Warn user of exceeded maximum number of iterations
    if( iteration > maximumNumberOfIterations )
    {
        std::cerr << "Warning in non-linear least squares estimation. Maximum number of iterations exceeded." << std::endl;
    }

    // Give out new estimate in parameters
    return currentEstimate;
}

}  // namespace linear_algebra

}  // namespace tudat
