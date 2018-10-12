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

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

#include "Tudat/Mathematics/Filters/linearKalmanFilter.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_linear_kalman_filter )

// Test implementation of linear Kalman filter class. Tested by comparison with code by https://github.com/hmartiro/kalman-cpp.
BOOST_AUTO_TEST_CASE( testLinearKalmanFilter )
{
    using namespace tudat::filters;

    // Set initial conditions
    const double initialTime = 0;
    const double timeStep = 1.0 / 30.0;
    const unsigned int numberOfTimeSteps = 45;

    Eigen::Vector3d initialStateVector = Eigen::Vector3d::Zero( );
    initialStateVector[ 0 ] = 1.04202710058;
    initialStateVector[ 2 ] = -9.81;

    Eigen::Vector3d initialEstimatedStateVector = Eigen::Vector3d::Zero( );
    initialEstimatedStateVector[ 0 ] = 1.04202710058;
    initialEstimatedStateVector[ 2 ] = -9.81;

    Eigen::Matrix3d initialEstimatedStateCovarianceMatrix = Eigen::Matrix3d::Zero( );
    initialEstimatedStateCovarianceMatrix( 0, 0 ) = 0.1;
    initialEstimatedStateCovarianceMatrix( 0, 1 ) = 0.1;
    initialEstimatedStateCovarianceMatrix( 0, 2 ) = 0.1;
    initialEstimatedStateCovarianceMatrix( 1, 0 ) = 0.1;
    initialEstimatedStateCovarianceMatrix( 1, 1 ) = 10000.0;
    initialEstimatedStateCovarianceMatrix( 1, 2 ) = 10.0;
    initialEstimatedStateCovarianceMatrix( 2, 0 ) = 0.1;
    initialEstimatedStateCovarianceMatrix( 2, 1 ) = 10.0;
    initialEstimatedStateCovarianceMatrix( 2, 2 ) = 100.0;

    // Set system dynamics and measurement
    Eigen::Matrix3d stateTransitionMatrix = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d controlMatrix = Eigen::Matrix3d::Zero( );
    Eigen::RowVector3d measurementMatrix = Eigen::RowVector3d::Zero( );
    stateTransitionMatrix( 0, 0 ) = 1.0;
    stateTransitionMatrix( 0, 1 ) = timeStep;
    stateTransitionMatrix( 1, 1 ) = 1.0;
    stateTransitionMatrix( 1, 2 ) = timeStep;
    stateTransitionMatrix( 2, 2 ) = 1.0;
    measurementMatrix[ 0 ] = 1.0;

    // Set system and measurement uncertainty
    Eigen::Matrix3d systemUncertainty = Eigen::Matrix3d::Zero( );
    Eigen::Vector1d measurementUncertainty = Eigen::Vector1d::Zero( );
    systemUncertainty( 0, 0 ) = 0.05;
    systemUncertainty( 0, 1 ) = 0.05;
    systemUncertainty( 1, 0 ) = 0.05;
    systemUncertainty( 1, 1 ) = 0.05;
    measurementUncertainty[ 0 ] = 0.5;

    // Create linear Kalman filter object
    KalmanFilterDoublePointer linearFilter = std::make_shared< LinearKalmanFilterDouble >(
                [ & ]( const double, const Eigen::Vector3d& ){ return stateTransitionMatrix; },
                [ & ]( const double, const Eigen::Vector3d& ){ return controlMatrix; },
                [ & ]( const double, const Eigen::Vector3d& ){ return measurementMatrix; },
                systemUncertainty, measurementUncertainty, timeStep,
                initialTime, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix );

    // Load noise from file
    Eigen::MatrixXd systemNoise = input_output::readMatrixFromFile( tudat::input_output::getTudatRootPath( ) +
                                                                    "/Mathematics/Filters/UnitTests/noiseData/lkfSystemNoise1.dat" );
    Eigen::MatrixXd measurementNoise = input_output::readMatrixFromFile( tudat::input_output::getTudatRootPath( ) +
                                                                         "/Mathematics/Filters/UnitTests/noiseData/lkfMeasurementNoise1.dat" );

    // Loop over each time step
    const bool showProgress = false;
    double currentTime = initialTime;
    Eigen::Vector3d currentStateVector = initialStateVector;
    Eigen::Vector3d currentControlVector = Eigen::Vector3d::Zero( );
    Eigen::Vector1d currentMeasurementVector;
    std::map< double, Eigen::Vector3d > actualStateVectorHistory;
    std::map< double, Eigen::Vector1d > measurementVectorHistory;
    actualStateVectorHistory[ initialTime ] = initialStateVector;
    for ( unsigned int i = 0; i < numberOfTimeSteps; i++ )
    {
        // Compute actual values and perturb them
        currentStateVector = stateTransitionMatrix * currentStateVector + controlMatrix * currentControlVector +
                systemNoise.col( i );
        currentMeasurementVector = measurementMatrix * currentStateVector + measurementNoise.col( i );

        // Update filter
        linearFilter->updateFilter( currentMeasurementVector );
        currentTime = linearFilter->getCurrentTime( );

        // Store values
        actualStateVectorHistory[ currentTime ] = currentStateVector;
        measurementVectorHistory[ currentTime ] = currentMeasurementVector;

        // Print progress
        if ( showProgress )
        {
            std::cout << "Time: " << currentTime << std::endl
                      << "Measurement: " << currentMeasurementVector.transpose( ) << std::endl
                      << "Estimated State: " << linearFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
        }
    }

    // Check that final state is as expected
    Eigen::Vector3d expectedFinalState;
    expectedFinalState << -7.415393533447765, -12.421405468923618, -8.9114310345206711;
    for ( unsigned int i = 0; i < expectedFinalState.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( linearFilter->getCurrentStateEstimate( )[ i ] - expectedFinalState[ i ], 1.0e-10 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
