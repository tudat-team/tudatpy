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

#include "tudat/basics/utilities.h"
#include "tudat/math/basic/leastSquaresEstimation.h"

namespace tudat
{

namespace linear_algebra
{

//! Function to get condition number of matrix (using SVD decomposition)
double getConditionNumberOfDesignMatrix( const Eigen::MatrixXd designMatrix )
{
    return getConditionNumberOfDecomposedMatrix(
                ( designMatrix.jacobiSvd( Eigen::ComputeThinU | Eigen::ComputeFullV ) ) );
}

//! Function to get condition number of matrix from SVD decomposition
double getConditionNumberOfDecomposedMatrix( const Eigen::JacobiSVD< Eigen::MatrixXd >& singularValueDecomposition )
{
    Eigen::VectorXd singularValues = singularValueDecomposition.singularValues( );
    return singularValues( 0 ) / singularValues( singularValues.rows( ) - 1 );
}

//! Solve system of equations with SVD decomposition, checking condition number in the process
Eigen::VectorXd solveSystemOfEquationsWithSvd( const Eigen::MatrixXd matrixToInvert,
                                               const Eigen::VectorXd rightHandSideVector,
                                               const bool checkConditionNumber,
                                               const double maximumAllowedConditionNumber )
{
    Eigen::JacobiSVD< Eigen::MatrixXd > svdDecomposition = matrixToInvert.jacobiSvd(
                Eigen::ComputeThinU | Eigen::ComputeThinV );
    if( checkConditionNumber )
    {
        double conditionNumber = getConditionNumberOfDecomposedMatrix( svdDecomposition );

        if( conditionNumber > maximumAllowedConditionNumber )
        {
            std::cerr << "Warning when performing least squares, condition number is " << conditionNumber << std::endl;
        }
    }
    return svdDecomposition.solve( rightHandSideVector );
}

//! Function to multiply information matrix by diagonal weights matrix
Eigen::MatrixXd multiplyDesignMatrixByDiagonalWeightMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix )
{
    Eigen::MatrixXd weightedDesignMatrix = Eigen::MatrixXd::Zero( designMatrix.rows( ), designMatrix.cols( ) );

    for( int i = 0; i < designMatrix.cols( ); i++ )
    {
        weightedDesignMatrix.block( 0, i, designMatrix.rows( ), 1 ) =
                designMatrix.block( 0, i, designMatrix.rows( ), 1 ).cwiseProduct( diagonalOfWeightMatrix );
    }

    return weightedDesignMatrix;
}

Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const Eigen::MatrixXd& constraintMultiplier,
        const Eigen::VectorXd& constraintRightHandside )
{
    // Add constraints to inverse covariance matrix if required
    Eigen::MatrixXd inverseOfCovarianceMatrix =
            inverseOfAPrioriCovarianceMatrix + designMatrix.transpose( ) * multiplyDesignMatrixByDiagonalWeightMatrix(
                designMatrix, diagonalOfWeightMatrix );
    if( constraintMultiplier.rows( ) != 0 )
    {
        if( constraintMultiplier.rows( ) != constraintRightHandside.rows( ) )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible" );
        }

        if( constraintMultiplier.cols( ) != designMatrix.cols( ) )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible with partials" );
        }

        int numberOfConstraints = constraintMultiplier.rows( );
        int numberOfParameters = constraintMultiplier.cols( );

        inverseOfCovarianceMatrix.conservativeResize(
                    numberOfParameters + numberOfConstraints, numberOfParameters + numberOfConstraints );
        inverseOfCovarianceMatrix.block( numberOfParameters, 0, numberOfConstraints, numberOfParameters ) =
               constraintMultiplier;
        inverseOfCovarianceMatrix.block( 0, numberOfParameters, numberOfParameters, numberOfConstraints ) =
               constraintMultiplier.transpose( );
        inverseOfCovarianceMatrix.block(
                    numberOfParameters, numberOfParameters, numberOfConstraints, numberOfConstraints ).setZero( );
    }

    return inverseOfCovarianceMatrix;

}


//! Function to compute inverse of covariance matrix at current iteration
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix )
{
    return calculateInverseOfUpdatedCovarianceMatrix(
                designMatrix, diagonalOfWeightMatrix,
                Eigen::MatrixXd::Zero( designMatrix.cols( ), designMatrix.cols( ) ) );
}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals and a priori
//! information
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const bool checkConditionNumber,
        const double maximumAllowedConditionNumber,
        const Eigen::MatrixXd& constraintMultiplier,
        const Eigen::VectorXd& constraintRightHandside,
        const Eigen::MatrixXd& designMatrixConsiderParameters,
        const Eigen::VectorXd& considerParametersValues )
{
    Eigen::VectorXd rightHandSide = Eigen::VectorXd::Zero( observationResiduals.size( ) );
    if ( considerParametersValues.size( ) > 0 && designMatrixConsiderParameters.size( ) > 0 )
    {
        rightHandSide = designMatrix.transpose( ) *
                        ( diagonalOfWeightMatrix.cwiseProduct( observationResiduals + designMatrixConsiderParameters * considerParametersValues ) );
    }
    else
    {
        rightHandSide = designMatrix.transpose( ) * ( diagonalOfWeightMatrix.cwiseProduct( observationResiduals ) );
    }

    Eigen::MatrixXd inverseOfCovarianceMatrix = calculateInverseOfUpdatedCovarianceMatrix(
                designMatrix, diagonalOfWeightMatrix, inverseOfAPrioriCovarianceMatrix,
                constraintMultiplier, constraintRightHandside );

    // Add constraints to inverse covariance matrix if required
    if( constraintMultiplier.rows( ) != 0 )
    {
        int numberOfConstraints = constraintMultiplier.rows( );
        int numberOfParameters = constraintMultiplier.cols( );

        rightHandSide.conservativeResize( numberOfParameters + numberOfConstraints );
        rightHandSide.segment( numberOfParameters, numberOfConstraints ) = constraintRightHandside;
    }

    return std::make_pair( solveSystemOfEquationsWithSvd(
            inverseOfCovarianceMatrix, rightHandSide, checkConditionNumber, maximumAllowedConditionNumber ), inverseOfCovarianceMatrix );

}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const bool checkConditionNumber,
        const double maximumAllowedConditionNumber )
{
    return performLeastSquaresAdjustmentFromDesignMatrix(
                designMatrix, observationResiduals, diagonalOfWeightMatrix,
                Eigen::MatrixXd::Zero( designMatrix.cols( ), designMatrix.cols( ) ),
                checkConditionNumber, maximumAllowedConditionNumber );
}

//! Function to perform an iteration of least squares estimation from information matrix and residuals
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& observationResiduals,
        const bool checkConditionNumber,
        const double maximumAllowedConditionNumber )
{
    return performLeastSquaresAdjustmentFromDesignMatrix(
                designMatrix, observationResiduals, Eigen::VectorXd::Constant( observationResiduals.size( ), 1, 1.0 ),
                checkConditionNumber, maximumAllowedConditionNumber );
}

//! Function to fit a univariate polynomial through a set of data
Eigen::VectorXd getLeastSquaresPolynomialFit(
        const Eigen::VectorXd& independentValues,
        const Eigen::VectorXd& dependentValues,
        const std::vector< double >& polynomialPowers )
{
    if( independentValues.rows( ) != dependentValues.rows( ) )
    {
        throw std::runtime_error( "Error when doing least squares polynomial fit, size of dependent and independent "
                                  "variable vectors is not equal." );
    }

    Eigen::MatrixXd designMatrix = Eigen::MatrixXd::Zero( dependentValues.rows( ), polynomialPowers.size( ) );

    // Compute information matrix
    for( int i = 0; i < independentValues.rows( ); i++ )
    {
        for( unsigned int j = 0; j < polynomialPowers.size( ); j++ )
        {
            designMatrix( i, j ) = std::pow( independentValues( i ), polynomialPowers.at( j ) );
        }
    }

    return performLeastSquaresAdjustmentFromDesignMatrix( designMatrix, dependentValues ).first;
}

//! Function to fit a univariate polynomial through a set of data
std::vector< double > getLeastSquaresPolynomialFit(
        const std::map< double, double >& independentDependentValueMap,
        const std::vector< double >& polynomialPowers )
{
    return utilities::convertEigenVectorToStlVector(
                getLeastSquaresPolynomialFit(
                    utilities::convertStlVectorToEigenVector(
                        utilities::createVectorFromMapKeys( independentDependentValueMap ) ),
                    utilities::convertStlVectorToEigenVector(
                        utilities::createVectorFromMapValues( independentDependentValueMap ) ), polynomialPowers ) );

}

//! Function to perform a non-linear least squares estimation with the Levenberg-Marquardt method.
Eigen::VectorXd nonLinearLeastSquaresFit(
        const std::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& observationAndJacobianFunctions,
        const Eigen::VectorXd& initialEstimate, const Eigen::VectorXd& actualObservations, const double initialScaling,
        const double convergenceTolerance, const unsigned int maximumNumberOfIterations )
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
        if ( iteration == 0 )
        {
            levenbergMarquardtDampingParameter = initialScaling * ( designMatrix.transpose( ) * designMatrix ).diagonal( ).maxCoeff( );
        }

        // Compute update in estimate
        Eigen::VectorXd diagonalOfWeightMatrix = Eigen::VectorXd::Ones( offsetInObservations.rows( ) );
        Eigen::MatrixXd inverseOfAPrioriCovarianceMatrix = levenbergMarquardtDampingParameter *
//                Eigen::MatrixXd( ( designMatrix.transpose( ) * designMatrix ).diagonal( ).asDiagonal( ) ); // Marquardt’s update
                Eigen::MatrixXd::Identity( currentEstimate.rows( ), currentEstimate.rows( ) );
        updateInEstimate = linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix(
                    designMatrix, offsetInObservations, diagonalOfWeightMatrix, inverseOfAPrioriCovarianceMatrix, false ).first;

        // Check that update is real
        if ( ( !updateInEstimate.allFinite( ) ) || ( updateInEstimate.hasNaN( ) ) )
        {
            throw std::runtime_error( "Error in non-linear least squares estimation. Iterative process diverges." );
        }

        // Compute gain ratio
        levenbergMarquardtGainRatio =
                ( offsetInObservations.squaredNorm( ) -
                  ( actualObservations - observationAndJacobianFunctions( currentEstimate + updateInEstimate ).first ).squaredNorm( ) ) /
                ( updateInEstimate.transpose( ) * ( levenbergMarquardtDampingParameter * updateInEstimate +
                                                    designMatrix.transpose( ) * offsetInObservations ) );

        // Update damping parameter
        if ( levenbergMarquardtGainRatio > 0 )
        {
            // Reduce damping parameter, since good approximation
            levenbergMarquardtDampingParameter *= std::max( 1.0 / 3.0, 1.0 - std::pow( 2.0 * levenbergMarquardtGainRatio - 1.0, 3 ) );
            scalingParameterUpdate = 2; // reset

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
    }
    while ( ( updateInEstimate.norm( ) > convergenceTolerance ) && ( iteration <= maximumNumberOfIterations ) );

    // Warn user of exceeded maximum number of iterations
    if ( iteration > maximumNumberOfIterations )
    {
        std::cerr << "Warning in non-linear least squares estimation. Maximum number of iterations exceeded." << std::endl;
    }

    // Give out new estimate in parameters
    return currentEstimate;
}

} // namespace linear_algebra

} // namespace tudat
