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
#include <limits>

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
    return getConditionNumberOfDecomposedMatrix( ( designMatrix.jacobiSvd< Eigen::ComputeThinU | Eigen::ComputeFullV >( ) ) );
}

//! Solve system of equations with SVD decomposition, checking condition number in the process
Eigen::VectorXd solveSystemOfEquationsWithSvd( const Eigen::MatrixXd matrixToInvert,
                                               const Eigen::VectorXd rightHandSideVector,
                                               const double limitConditionNumberForWarning )
{
    Eigen::JacobiSVD< Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV > svdDecomposition =
            matrixToInvert.jacobiSvd< Eigen::ComputeThinU | Eigen::ComputeThinV >( );
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
namespace
{

//! Applies diagonal weights to each row of the design matrix.
Eigen::MatrixXd applyWeightsToDesignMatrix( const Eigen::MatrixXd& designMatrix, const Eigen::VectorXd& diagonalOfWeightMatrix )
{
    Eigen::MatrixXd weightedDesignMatrix = Eigen::MatrixXd::Zero( designMatrix.rows( ), designMatrix.cols( ) );

    for( int i = 0; i < designMatrix.cols( ); i++ )
    {
        weightedDesignMatrix.block( 0, i, designMatrix.rows( ), 1 ) =
                designMatrix.block( 0, i, designMatrix.rows( ), 1 ).cwiseProduct( diagonalOfWeightMatrix );
    }

    return weightedDesignMatrix;
}

//! Applies a full sparse weight matrix to the design matrix.
Eigen::MatrixXd applyWeightsToDesignMatrix( const Eigen::MatrixXd& designMatrix, const Eigen::SparseMatrix< double >& weightMatrix )
{
    if( weightMatrix.rows( ) != designMatrix.rows( ) || weightMatrix.cols( ) != designMatrix.rows( ) )
    {
        throw std::runtime_error( "Error when multiplying design matrix by weights matrix, sizes are incompatible." );
    }
    return ( weightMatrix * designMatrix ).eval( );
}

//! Applies diagonal weights to an observation/residual vector.
Eigen::VectorXd applyWeightsToObservationVector( const Eigen::VectorXd& observationVector, const Eigen::VectorXd& diagonalOfWeightMatrix )
{
    return diagonalOfWeightMatrix.cwiseProduct( observationVector );
}

//! Applies a full sparse weight matrix to an observation/residual vector.
Eigen::VectorXd applyWeightsToObservationVector( const Eigen::VectorXd& observationVector,
                                                 const Eigen::SparseMatrix< double >& weightMatrix )
{
    if( weightMatrix.rows( ) != observationVector.rows( ) || weightMatrix.cols( ) != observationVector.rows( ) )
    {
        throw std::runtime_error( "Error when multiplying observation vector by weights matrix, sizes are incompatible." );
    }
    return weightMatrix * observationVector;
}

//! Augments the normal matrix with linear constraint rows/columns when constraints are provided.
void addConstraintsToInverseCovarianceMatrix( Eigen::MatrixXd& inverseOfCovarianceMatrix,
                                              const Eigen::MatrixXd& designMatrix,
                                              const Eigen::MatrixXd& constraintMultiplier,
                                              const Eigen::VectorXd& constraintRightHandside )
{
    if( constraintMultiplier.rows( ) == 0 )
    {
        return;
    }

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

    inverseOfCovarianceMatrix.conservativeResize( numberOfParameters + numberOfConstraints, numberOfParameters + numberOfConstraints );
    inverseOfCovarianceMatrix.block( numberOfParameters, 0, numberOfConstraints, numberOfParameters ) = constraintMultiplier;
    inverseOfCovarianceMatrix.block( 0, numberOfParameters, numberOfParameters, numberOfConstraints ) = constraintMultiplier.transpose( );
    inverseOfCovarianceMatrix.block( numberOfParameters, numberOfParameters, numberOfConstraints, numberOfConstraints ).setZero( );
}

template< typename WeightType >
//! Shared implementation for dense-diagonal and sparse full-weight covariance updates.
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrixImplementation( const Eigen::MatrixXd& designMatrix,
                                                                         const WeightType& weightData,
                                                                         const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
                                                                         const Eigen::MatrixXd& constraintMultiplier,
                                                                         const Eigen::VectorXd& constraintRightHandside )
{
    // Build normal matrix: H^T W H + P0^-1.
    Eigen::MatrixXd inverseOfCovarianceMatrix =
            inverseOfAPrioriCovarianceMatrix + designMatrix.transpose( ) * applyWeightsToDesignMatrix( designMatrix, weightData );
    // Append constraint equations when requested.
    addConstraintsToInverseCovarianceMatrix( inverseOfCovarianceMatrix, designMatrix, constraintMultiplier, constraintRightHandside );
    return inverseOfCovarianceMatrix;
}

template< typename WeightType >
//! Shared implementation for consider-parameter covariance contribution with generic weight representation.
Eigen::MatrixXd calculateConsiderParametersCovarianceContributionImplementation( const Eigen::MatrixXd& normalisedCovarianceMatrix,
                                                                                 const Eigen::MatrixXd& designMatrix,
                                                                                 const WeightType& weightData,
                                                                                 const Eigen::MatrixXd& considerDesignMatrix,
                                                                                 const Eigen::MatrixXd& considerCovariance )
{
    Eigen::MatrixXd covarianceTimesWeightedPartials =
            normalisedCovarianceMatrix * applyWeightsToDesignMatrix( designMatrix, weightData ).transpose( );
    return ( covarianceTimesWeightedPartials * considerDesignMatrix ) * considerCovariance *
            ( considerDesignMatrix.transpose( ) * covarianceTimesWeightedPartials.transpose( ) );
}

template< typename WeightType >
//! Shared least-squares implementation for diagonal and sparse full weight inputs.
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrixImplementation(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& observationResiduals,
        const WeightType& weightData,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const double limitConditionNumberForWarning,
        const Eigen::MatrixXd& constraintMultiplier,
        const Eigen::VectorXd& constraintRightHandside,
        const Eigen::MatrixXd& designMatrixConsiderParameters,
        const Eigen::VectorXd& considerParametersDeviations,
        const Eigen::MatrixXd& additionalNormalMatrix,
        const Eigen::VectorXd& additionalRightHandSide )
{
    // Build weighted right-hand side, including consider-parameter deviations when provided.
    Eigen::VectorXd weightedRightHandSideArgument = observationResiduals;
    if( considerParametersDeviations.size( ) > 0 && designMatrixConsiderParameters.size( ) > 0 )
    {
        weightedRightHandSideArgument += designMatrixConsiderParameters * considerParametersDeviations;
    }
    Eigen::VectorXd rightHandSide =
            designMatrix.transpose( ) * applyWeightsToObservationVector( weightedRightHandSideArgument, weightData );

    Eigen::MatrixXd inverseOfCovarianceMatrix = calculateInverseOfUpdatedCovarianceMatrixImplementation(
            designMatrix, weightData, inverseOfAPrioriCovarianceMatrix, constraintMultiplier, constraintRightHandside );

    // Extend RHS with constraint values in the augmented system.
    if( constraintMultiplier.rows( ) != 0 )
    {
        int numberOfConstraints = constraintMultiplier.rows( );
        int numberOfParameters = constraintMultiplier.cols( );

        rightHandSide.conservativeResize( numberOfParameters + numberOfConstraints );
        rightHandSide.segment( numberOfParameters, numberOfConstraints ) = constraintRightHandside;
    }

    // Soft-constraint additions (e.g. inter-arc continuity) affect only the parameter block,
    // leaving any Lagrange-multiplier rows from hard equality constraints untouched.
    const int numberOfParameters = static_cast< int >( designMatrix.cols( ) );
    if( additionalNormalMatrix.size( ) > 0 )
    {
        if( additionalNormalMatrix.rows( ) != numberOfParameters || additionalNormalMatrix.cols( ) != numberOfParameters )
        {
            throw std::runtime_error( "Error in performLeastSquaresAdjustmentFromDesignMatrix: additionalNormalMatrix has shape " +
                                      std::to_string( additionalNormalMatrix.rows( ) ) + "x" +
                                      std::to_string( additionalNormalMatrix.cols( ) ) + ", expected " +
                                      std::to_string( numberOfParameters ) + "x" + std::to_string( numberOfParameters ) + "." );
        }
        inverseOfCovarianceMatrix.topLeftCorner( numberOfParameters, numberOfParameters ) += additionalNormalMatrix;
    }
    if( additionalRightHandSide.size( ) > 0 )
    {
        if( additionalRightHandSide.size( ) != numberOfParameters )
        {
            throw std::runtime_error( "Error in performLeastSquaresAdjustmentFromDesignMatrix: additionalRightHandSide has size " +
                                      std::to_string( additionalRightHandSide.size( ) ) + ", expected " +
                                      std::to_string( numberOfParameters ) + "." );
        }
        rightHandSide.head( numberOfParameters ) += additionalRightHandSide;
    }

    // Solve normal equations and return both correction and assembled normal matrix.
    return std::make_pair( solveSystemOfEquationsWithSvd( inverseOfCovarianceMatrix, rightHandSide, limitConditionNumberForWarning ),
                           inverseOfCovarianceMatrix );
}

}  // namespace

Eigen::MatrixXd multiplyDesignMatrixByDiagonalWeightMatrix( const Eigen::MatrixXd& designMatrix,
                                                            const Eigen::VectorXd& diagonalOfWeightMatrix )
{
    return applyWeightsToDesignMatrix( designMatrix, diagonalOfWeightMatrix );
}

Eigen::MatrixXd multiplyDesignMatrixByWeightMatrix( const Eigen::MatrixXd& designMatrix, const Eigen::SparseMatrix< double >& weightMatrix )
{
    return applyWeightsToDesignMatrix( designMatrix, weightMatrix );
}

Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix( const Eigen::MatrixXd& designMatrix,
                                                           const Eigen::VectorXd& diagonalOfWeightMatrix,
                                                           const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
                                                           const Eigen::MatrixXd& constraintMultiplier,
                                                           const Eigen::VectorXd& constraintRightHandside,
                                                           const double limitConditionNumberForWarning )
{
    return calculateInverseOfUpdatedCovarianceMatrixImplementation(
            designMatrix, diagonalOfWeightMatrix, inverseOfAPrioriCovarianceMatrix, constraintMultiplier, constraintRightHandside );
}

Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix( const Eigen::MatrixXd& designMatrix,
                                                           const Eigen::SparseMatrix< double >& weightMatrix,
                                                           const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
                                                           const Eigen::MatrixXd& constraintMultiplier,
                                                           const Eigen::VectorXd& constraintRightHandside,
                                                           const double limitConditionNumberForWarning )
{
    return calculateInverseOfUpdatedCovarianceMatrixImplementation(
            designMatrix, weightMatrix, inverseOfAPrioriCovarianceMatrix, constraintMultiplier, constraintRightHandside );
}

//! Function to compute inverse of covariance matrix at current iteration
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix( const Eigen::MatrixXd& designMatrix,
                                                           const Eigen::VectorXd& diagonalOfWeightMatrix,
                                                           const double limitConditionNumberForWarning )
{
    return calculateInverseOfUpdatedCovarianceMatrix(
            designMatrix, diagonalOfWeightMatrix, Eigen::MatrixXd::Zero( designMatrix.cols( ), designMatrix.cols( ) ) );
}

Eigen::MatrixXd calculateConsiderParametersCovarianceContribution( const Eigen::MatrixXd& normalisedCovarianceMatrix,
                                                                   const Eigen::MatrixXd& designMatrix,
                                                                   const Eigen::VectorXd& diagonalOfWeightMatrix,
                                                                   const Eigen::MatrixXd& considerDesignMatrix,
                                                                   const Eigen::MatrixXd& considerCovariance )
{
    return calculateConsiderParametersCovarianceContributionImplementation(
            normalisedCovarianceMatrix, designMatrix, diagonalOfWeightMatrix, considerDesignMatrix, considerCovariance );
}

Eigen::MatrixXd calculateConsiderParametersCovarianceContribution( const Eigen::MatrixXd& normalisedCovarianceMatrix,
                                                                   const Eigen::MatrixXd& designMatrix,
                                                                   const Eigen::SparseMatrix< double >& weightMatrix,
                                                                   const Eigen::MatrixXd& considerDesignMatrix,
                                                                   const Eigen::MatrixXd& considerCovariance )
{
    return calculateConsiderParametersCovarianceContributionImplementation(
            normalisedCovarianceMatrix, designMatrix, weightMatrix, considerDesignMatrix, considerCovariance );
}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals and a priori
//! information
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const double limitConditionNumberForWarning,
        const Eigen::MatrixXd& constraintMultiplier,
        const Eigen::VectorXd& constraintRightHandside,
        const Eigen::MatrixXd& designMatrixConsiderParameters,
        const Eigen::VectorXd& considerParametersDeviations,
        const Eigen::MatrixXd& additionalNormalMatrix,
        const Eigen::VectorXd& additionalRightHandSide )
{
    return performLeastSquaresAdjustmentFromDesignMatrixImplementation( designMatrix,
                                                                        observationResiduals,
                                                                        diagonalOfWeightMatrix,
                                                                        inverseOfAPrioriCovarianceMatrix,
                                                                        limitConditionNumberForWarning,
                                                                        constraintMultiplier,
                                                                        constraintRightHandside,
                                                                        designMatrixConsiderParameters,
                                                                        considerParametersDeviations,
                                                                        additionalNormalMatrix,
                                                                        additionalRightHandSide );
}

std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::SparseMatrix< double >& weightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const double limitConditionNumberForWarning,
        const Eigen::MatrixXd& constraintMultiplier,
        const Eigen::VectorXd& constraintRightHandside,
        const Eigen::MatrixXd& designMatrixConsiderParameters,
        const Eigen::VectorXd& considerParametersDeviations,
        const Eigen::MatrixXd& additionalNormalMatrix,
        const Eigen::VectorXd& additionalRightHandSide )
{
    return performLeastSquaresAdjustmentFromDesignMatrixImplementation( designMatrix,
                                                                        observationResiduals,
                                                                        weightMatrix,
                                                                        inverseOfAPrioriCovarianceMatrix,
                                                                        limitConditionNumberForWarning,
                                                                        constraintMultiplier,
                                                                        constraintRightHandside,
                                                                        designMatrixConsiderParameters,
                                                                        considerParametersDeviations,
                                                                        additionalNormalMatrix,
                                                                        additionalRightHandSide );
}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix( const Eigen::MatrixXd& designMatrix,
                                                                                             const Eigen::VectorXd& observationResiduals,
                                                                                             const Eigen::VectorXd& diagonalOfWeightMatrix,
                                                                                             const double limitConditionNumberForWarning )
{
    return performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                          observationResiduals,
                                                          diagonalOfWeightMatrix,
                                                          Eigen::MatrixXd::Zero( designMatrix.cols( ), designMatrix.cols( ) ),
                                                          limitConditionNumberForWarning );
}

std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromDesignMatrix(
        const Eigen::MatrixXd& designMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::SparseMatrix< double >& weightMatrix,
        const double limitConditionNumberForWarning )
{
    return performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                          observationResiduals,
                                                          weightMatrix,
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

Eigen::VectorXd evaluatePolynomial( const Eigen::VectorXd& independentValues,
                                    const Eigen::VectorXd& polynomialCoefficients,
                                    const std::vector< double >& polynomialPowers )
{
    Eigen::VectorXd polynomial = Eigen::VectorXd::Zero( independentValues.rows( ) );
    for( unsigned int i = 0; i < polynomialPowers.size( ); i++ )
    {
        polynomial += polynomialCoefficients( i ) * independentValues.array( ).pow( polynomialPowers.at( i ) ).matrix( );
    }
    return polynomial;
}

//! Function to fit a univariate polynomial through a set of data
Eigen::VectorXd getLeastSquaresPolynomialFit( const Eigen::VectorXd& independentValues,
                                              const Eigen::VectorXd& dependentValues,
                                              const std::vector< double >& polynomialPowers )
{
    if( independentValues.rows( ) != dependentValues.rows( ) )
    {
        throw std::runtime_error(
                "Error when doing least squares polynomial fit, size of dependent and independent "
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
std::vector< double > getLeastSquaresPolynomialFit( const std::map< double, double >& independentDependentValueMap,
                                                    const std::vector< double >& polynomialPowers )
{
    return utilities::convertEigenVectorToStlVector( getLeastSquaresPolynomialFit(
            utilities::convertStlVectorToEigenVector( utilities::createVectorFromMapKeys( independentDependentValueMap ) ),
            utilities::convertStlVectorToEigenVector( utilities::createVectorFromMapValues( independentDependentValueMap ) ),
            polynomialPowers ) );
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
    unsigned int numberOfIterations;
    bool converged;
    return nonLinearLeastSquaresFit( observationAndJacobianFunctions,
                                     initialEstimate,
                                     actualObservations,
                                     initialScaling,
                                     convergenceTolerance,
                                     maximumNumberOfIterations,
                                     numberOfIterations,
                                     converged );
}

//! Function to perform a non-linear least squares estimation and return iteration diagnostics.
Eigen::VectorXd nonLinearLeastSquaresFit(
        const std::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& observationAndJacobianFunctions,
        const Eigen::VectorXd& initialEstimate,
        const Eigen::VectorXd& actualObservations,
        const double initialScaling,
        const double convergenceTolerance,
        const unsigned int maximumNumberOfIterations,
        unsigned int& numberOfIterations,
        bool& converged )
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
        // The nonlinear solver handles poor conditioning through Levenberg-Marquardt damping. Disable the lower-level
        // SVD warning here to avoid printing once per iteration; divergence is still detected by the finite-update check.
        updateInEstimate = linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                                                          offsetInObservations,
                                                                                          diagonalOfWeightMatrix,
                                                                                          inverseOfAPrioriCovarianceMatrix,
                                                                                          std::numeric_limits< double >::infinity( ) )
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

    numberOfIterations = iteration;
    converged = updateInEstimate.norm( ) <= convergenceTolerance;

    // Give out new estimate in parameters
    return currentEstimate;
}

}  // namespace linear_algebra

}  // namespace tudat
