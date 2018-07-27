#include <iostream>

#include <Eigen/LU>

#include "Tudat/Mathematics/BasicMathematics/nonLinearLeastSquaresEstimation.h"

namespace tudat
{

namespace linear_algebra
{

//! Function to perform a non-linear least squares estimation.
Eigen::VectorXd nonLinearLeastSquaresFit(
        const boost::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& observationAndJacobianFunctions,
        const Eigen::VectorXd& initialEstimate, const Eigen::VectorXd& actualObservations,
        const double convergenceTolerance, const unsigned int maximumNumberOfIterations )
{
    // Set current estimate to initial value
    Eigen::VectorXd currentEstimate = initialEstimate;
//    std::cout << "Measured: " << actualObservations.transpose( ) << std::endl;

    // Initialize variables
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > pairOfEstimatedObservationsAndDesignMatrix;
    Eigen::MatrixXd designMatrix;
    Eigen::VectorXd offsetInObservations;
    Eigen::VectorXd updateInEstimate;

    // Start iterative loop
    unsigned int iteration = 0;
    do
    {
        // Compute current system and jacobian functions
        pairOfEstimatedObservationsAndDesignMatrix = observationAndJacobianFunctions( currentEstimate );
        designMatrix = pairOfEstimatedObservationsAndDesignMatrix.second;
//        std::cout << "Expected: " << pairOfEstimatedObservationsAndDesignMatrix.first.transpose( ) << std::endl;
//        std::cout << "Jacobian: " << designMatrix << std::endl;

        // Offset in observation
        offsetInObservations = actualObservations - pairOfEstimatedObservationsAndDesignMatrix.first;

        // Compute update in estimate
        updateInEstimate = ( designMatrix.transpose( ) * designMatrix ).inverse( ) * designMatrix.transpose( ) * offsetInObservations;
        std::cout << "Iter: " << iteration + 1 << ". Update: " << updateInEstimate.transpose( ) << std::endl;

        // Check that update is finite
        if ( ( !updateInEstimate.allFinite( ) ) || ( updateInEstimate.hasNaN( ) ) )
        {
            throw std::runtime_error( "Error in non-linear least squares estimation. Iterative process diverges." );
        }

        // Correct estimate
        currentEstimate += updateInEstimate;
        iteration++;
    }
    while ( ( updateInEstimate.norm( ) > convergenceTolerance ) && ( iteration < maximumNumberOfIterations ) );

    // Give out new estimate in parameters
    std::cout << "Final: " << currentEstimate.transpose( ) << ". Iterations: " << iteration << std::endl;
    return currentEstimate;
}

} // namespace linear_algebra

} // namespace tudat
