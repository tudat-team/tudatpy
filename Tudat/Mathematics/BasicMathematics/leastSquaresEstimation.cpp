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

#include <cmath>
#include <iostream>

#include <Eigen/LU>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Mathematics/BasicMathematics/leastSquaresEstimation.h"

namespace tudat
{

namespace linear_algebra
{

//! Function to get condition number of matrix (using SVD decomposition)
double getConditionNumberOfInformationMatrix( const Eigen::MatrixXd informationMatrix )
{
    return getConditionNumberOfDecomposedMatrix(
                ( informationMatrix.jacobiSvd( Eigen::ComputeThinU | Eigen::ComputeFullV ) ) );
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
Eigen::MatrixXd multiplyInformationMatrixByDiagonalWeightMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix )
{
    Eigen::MatrixXd weightedInformationMatrix = Eigen::MatrixXd::Zero( informationMatrix.rows( ), informationMatrix.cols( ) );

    for( int i = 0; i < informationMatrix.cols( ); i++ )
    {
        weightedInformationMatrix.block( 0, i, informationMatrix.rows( ), 1 ) =
                informationMatrix.block( 0, i, informationMatrix.rows( ), 1 ).cwiseProduct( diagonalOfWeightMatrix );
    }

    return weightedInformationMatrix;
}

//! Function to compute inverse of covariance matrix at current iteration, including influence of a priori information
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix )
{
    return inverseOfAPrioriCovarianceMatrix + informationMatrix.transpose( ) * multiplyInformationMatrixByDiagonalWeightMatrix(
                informationMatrix, diagonalOfWeightMatrix );
}

//! Function to compute inverse of covariance matrix at current iteration
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix )
{
    return calculateInverseOfUpdatedCovarianceMatrix(
                informationMatrix, diagonalOfWeightMatrix,
                Eigen::MatrixXd::Zero( informationMatrix.cols( ), informationMatrix.cols( ) ) );
}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals and a priori
//! information
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const bool checkConditionNumber,
        const double maximumAllowedConditionNumber,
        const Eigen::MatrixXd& constraintMultiplier,
        const Eigen::VectorXd& constraintRightHandside )
{
//    std::cout<<"Residuals "<<observationResiduals.transpose( )<<std::endl;
//    std::cout<<"Weight diag. "<<diagonalOfWeightMatrix.transpose( )<<std::endl;
//    std::cout<<"Partials "<<informationMatrix.transpose( )<<std::endl;

    Eigen::VectorXd rightHandSide = informationMatrix.transpose( ) *
            ( diagonalOfWeightMatrix.cwiseProduct( observationResiduals ) );
    Eigen::MatrixXd inverseOfCovarianceMatrix = calculateInverseOfUpdatedCovarianceMatrix(
                informationMatrix, diagonalOfWeightMatrix, inverseOfAPrioriCovarianceMatrix );

    // Add constraints to inverse covariance matrix if required
    if( constraintMultiplier.rows( ) != 0 )
    {
        if( constraintMultiplier.rows( ) != constraintRightHandside.rows( ) )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible" );
        }

        if( constraintMultiplier.cols( ) != informationMatrix.cols( ) )
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

        rightHandSide.conservativeResize( numberOfParameters + numberOfConstraints );
        rightHandSide.segment( numberOfParameters, numberOfConstraints ) = constraintRightHandside;
    }

//    std::cout<<"RHS "<<rightHandSide.transpose( )<<std::endl;
//    std::cout<<"Inv cov "<<inverseOfCovarianceMatrix<<std::endl;

    return std::make_pair( solveSystemOfEquationsWithSvd(
                               inverseOfCovarianceMatrix, rightHandSide, checkConditionNumber, maximumAllowedConditionNumber ),
                           inverseOfCovarianceMatrix );

}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const bool checkConditionNumber,
        const double maximumAllowedConditionNumber )
{
    return performLeastSquaresAdjustmentFromInformationMatrix(
                informationMatrix, observationResiduals, diagonalOfWeightMatrix,
                Eigen::MatrixXd::Zero( informationMatrix.cols( ), informationMatrix.cols( ) ),
                checkConditionNumber, maximumAllowedConditionNumber );
}

//! Function to perform an iteration of least squares estimation from information matrix and residuals
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const bool checkConditionNumber,
        const double maximumAllowedConditionNumber )
{
    return performLeastSquaresAdjustmentFromInformationMatrix(
                informationMatrix, observationResiduals, Eigen::VectorXd::Constant( observationResiduals.size( ), 1, 1.0 ),
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

    Eigen::MatrixXd informationMatrix = Eigen::MatrixXd::Zero( dependentValues.rows( ), polynomialPowers.size( ) );

    // Compute information matrix
    for( int i = 0; i < independentValues.rows( ); i++ )
    {
        for( unsigned int j = 0; j < polynomialPowers.size( ); j++ )
        {
            informationMatrix( i, j ) = std::pow( independentValues( i ), polynomialPowers.at( j ) );
        }
    }

    return performLeastSquaresAdjustmentFromInformationMatrix( informationMatrix, dependentValues ).first;
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
        updateInEstimate = linear_algebra::performLeastSquaresAdjustmentFromInformationMatrix(
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
