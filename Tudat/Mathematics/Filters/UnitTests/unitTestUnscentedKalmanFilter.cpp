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

#include "Tudat/Basics/utilities.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/Mathematics/Filters/unscentedKalmanFilter.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

#include "Tudat/Mathematics/Filters/UnitTests/controlClass.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_unscented_kalman_filter )

// Functions for unscented Kalman filter
Eigen::Vector2d stateFunction1( const double time, const Eigen::Vector2d& state,
                                const Eigen::Vector2d& control )
{
    Eigen::Vector2d stateDerivative;
    stateDerivative[ 0 ] = state[ 1 ] * std::pow( std::cos( state[ 0 ] ), 3 );
    stateDerivative[ 1 ] = std::sin( state[ 0 ] );
    return stateDerivative;
}
Eigen::Vector1d measurementFunction1( const double time, const Eigen::Vector2d& state )
{
    Eigen::Vector1d measurement;
    measurement[ 0 ] = std::pow( state[ 0 ], 3 );
    return measurement;
}

// Test implementation of unscented Kalman filter class.
BOOST_AUTO_TEST_CASE( testUnscentedKalmanFilterFirstCase )
{
    using namespace tudat::filters;

    // Set initial conditions
    const double initialTime = 0;
    const double timeStep = 0.01;
    const unsigned int numberOfTimeSteps = 1000;

    Eigen::Vector2d initialStateVector;
    initialStateVector[ 0 ] = 3.0;
    initialStateVector[ 1 ] = -0.3;

    Eigen::Vector2d initialEstimatedStateVector;
    initialEstimatedStateVector[ 0 ] = 10.0;
    initialEstimatedStateVector[ 1 ] = -3;

    Eigen::Matrix2d initialEstimatedStateCovarianceMatrix = Eigen::Matrix2d::Zero( );
    initialEstimatedStateCovarianceMatrix( 0, 0 ) = 100;
    initialEstimatedStateCovarianceMatrix( 1, 1 ) = 100;

    // Set system and measurement uncertainty
    Eigen::Matrix2d systemUncertainty = Eigen::Matrix2d::Zero( );
    Eigen::Vector1d measurementUncertainty = Eigen::Vector1d::Zero( );
    systemUncertainty( 0, 0 ) = 100;
    systemUncertainty( 1, 1 ) = 100;
    measurementUncertainty[ 0 ] = 100;

    // Set integrator settings
    boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            boost::make_shared< numerical_integrators::IntegratorSettings< > > (
                numerical_integrators::euler, initialTime, timeStep );

    // Create control class
    boost::shared_ptr< ControlWrapper< double, double, 2 > > control =
            boost::make_shared< ControlWrapper< double, double, 2 > >( boost::lambda::constant( Eigen::Vector2d::Zero( ) ) );

    // Create extended Kalman filter object
    UnscentedKalmanFilterDoublePointer unscentedFilter = boost::make_shared< UnscentedKalmanFilterDouble >(
                boost::bind( &stateFunction1, _1, _2, control->getControlVector( ) ),
                boost::bind( &measurementFunction1, _1, _2 ),
                systemUncertainty, measurementUncertainty,
                initialTime, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix,
                integratorSettings );

    // Loop over each time step
    const bool showProgress = false;
    double currentTime = initialTime;
    Eigen::Vector2d currentStateVector = initialStateVector;
    Eigen::Vector2d currentControlVector = Eigen::Vector2d::Zero( );
    Eigen::Vector1d currentMeasurementVector;
    std::map< double, Eigen::Vector2d > actualStateVectorHistory;
    std::map< double, Eigen::Vector1d > measurementVectorHistory;
    actualStateVectorHistory[ initialTime ] = initialStateVector;
    for( unsigned int i = 0; i < numberOfTimeSteps; i++ )
    {
        // Compute actual values and perturb them
        currentTime += timeStep;
        currentStateVector += ( stateFunction1( currentTime, currentStateVector, currentControlVector ) +
                                unscentedFilter->produceSystemNoise( ) ) * timeStep;
        currentMeasurementVector = measurementFunction1( currentTime, currentStateVector ) +
                unscentedFilter->produceMeasurementNoise( );
        actualStateVectorHistory[ currentTime ] = currentStateVector;
        measurementVectorHistory[ currentTime ] = currentMeasurementVector;

        // Update control class
        control->setControlVector( currentTime, unscentedFilter->getCurrentStateEstimate( ) );

        // Update filter
        unscentedFilter->updateFilter( currentTime, currentMeasurementVector );

        // Print progress
        if ( showProgress )
        {
            std::cout << "Time: " << currentTime << std::endl
                      << "Measurement: " << currentMeasurementVector.transpose( ) << std::endl
                      << "Estimated State: " << unscentedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
        }
    }

    // Check that final state is as expected
    Eigen::Vector2d expectedFinalState;
    expectedFinalState << 4.9651546003074403, -12.544916512181523;
    for ( int i = 0; i < expectedFinalState.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( unscentedFilter->getCurrentStateEstimate( )[ i ] - expectedFinalState[ i ],
                           std::numeric_limits< double >::epsilon( ) );
    }

    // Check that noise is actually normally distributed (within 5 %)
    std::pair< std::vector< Eigen::VectorXd >, std::vector< Eigen::VectorXd > > noiseHistory = unscentedFilter->getNoiseHistory( );
    Eigen::MatrixXd systemNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.first );
    Eigen::MatrixXd measurementNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.second );
    for ( unsigned int i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( statistics::computeStandardDeviationOfVectorComponents( systemNoise.row( i ) ),
                                    std::sqrt( systemUncertainty( i, i ) ), 5e-2 );
    }
    BOOST_CHECK_CLOSE_FRACTION( statistics::computeStandardDeviationOfVectorComponents( measurementNoise.row( 0 ) ),
                                std::sqrt( measurementUncertainty( 0, 0 ) ), 5e-2 );

//    // Save actual state history
//    input_output::writeDataMapToTextFile( actualStateVectorHistory, "UKFActualStateHistory.dat", "/Users/Michele/Desktop/KF" );

//    // Save estimated state history
//    input_output::writeDataMapToTextFile( unscentedFilter->getEstimatedStateHistory( ), "UKFEstimatedStateHistory.dat",
//                                          "/Users/Michele/Desktop/KF" );

//    // Save measurement history
//    input_output::writeDataMapToTextFile( measurementVectorHistory, "UKFMeasurementHistory.dat", "/Users/Michele/Desktop/KF" );

//    // Save noise histories
//    systemNoise.transposeInPlace( );
//    measurementNoise.transposeInPlace( );
//    input_output::writeMatrixToFile( systemNoise, "systemNoise.dat", 16, "/Users/Michele/Desktop/KF" );
//    input_output::writeMatrixToFile( measurementNoise, "measurementNoise.dat", 16, "/Users/Michele/Desktop/KF" );

    // Save sigma points history
//    std::map< double, Eigen::MatrixXd > sigmaPointHistory = unscentedFilter->getSigmaPointsHistory( );
//    input_output::writeDataMapToTextFile( sigmaPointHistory, "UKFSigmaPoints.dat", "/Users/Michele/Desktop/KF" );
}

//! Typedefs.
typedef Eigen::Matrix< long double, 1, 1 > Vector1ld;
typedef Eigen::Matrix< long double, 3, 1 > Vector3ld;
typedef Eigen::Matrix< long double, 3, 3 > Matrix3ld;

// Functions for unscented Kalman filter
Vector3ld stateFunction2( const long double time, const Vector3ld& state, const Vector3ld& control )
{
    Vector3ld stateFunction;
    stateFunction[ 0 ] = state[ 1 ];
    stateFunction[ 1 ] = state[ 2 ];
    stateFunction[ 2 ] = 0.05 * state[ 0 ] * ( state[ 1 ] + state[ 2 ] );
    return stateFunction;
}
Vector1ld measurementFunction2( const long double time, const Vector3ld& state )
{
    Vector1ld measurement;
    measurement[ 0 ] = state[ 0 ];
    return measurement;
}

// Test implementation of unscented Kalman filter class. Tested by comparison with code by
// https://de.mathworks.com/matlabcentral/fileexchange/18217-learning-the-unscented-kalman-filter
BOOST_AUTO_TEST_CASE( testUnscentedKalmanFilterSecondCase )
{
    using namespace tudat::filters;

    // Set initial conditions
    const long double initialTime = 0;
    const long double timeStep = 1.0;
    const unsigned int numberOfTimeSteps = 100;

    Vector3ld initialStateVector = Vector3ld::Zero( );
    initialStateVector[ 2 ] = 1.0;

    Vector3ld initialEstimatedStateVector;
    initialEstimatedStateVector[ 0 ] = 0.05376671395461;
    initialEstimatedStateVector[ 1 ] = 0.183388501459509;
    initialEstimatedStateVector[ 2 ] = 0.774115313899635;

    Matrix3ld initialEstimatedStateCovarianceMatrix = Matrix3ld::Identity( );

    // Set system and measurement uncertainty
    Matrix3ld systemUncertainty = 0.01 * Matrix3ld::Identity( );
    Vector1ld measurementUncertainty = 0.01 * Vector1ld::Identity( );

    // Set null integrator settings
    boost::shared_ptr< numerical_integrators::IntegratorSettings< long double > > integratorSettings = NULL;

    // Create control class
    boost::shared_ptr< ControlWrapper< long double, long double, 3 > > control =
            boost::make_shared< ControlWrapper< long double, long double, 3 > >(
                boost::lambda::constant( Vector3ld::Zero( ) ) );

    // Create extended Kalman filter object
    boost::shared_ptr< KalmanFilterBase< long double, long double > > unscentedFilter =
            boost::make_shared< UnscentedKalmanFilter< long double, long double > >(
                boost::bind( &stateFunction2, _1, _2, control->getControlVector( ) ),
                boost::bind( &measurementFunction2, _1, _2 ),
                systemUncertainty, measurementUncertainty,
                initialTime, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix,
                integratorSettings, custom_parameters, std::make_pair( 0.001, 0.0 ) );

    // Loop over each time step
    const bool showProgress = false;
    long double currentTime = initialTime;
    Vector3ld currentStateVector = initialStateVector;
    Vector3ld currentControlVector = Vector3ld::Zero( );
    Vector1ld currentMeasurementVector;
    std::map< long double, Vector3ld > actualStateVectorHistory;
    std::map< long double, Vector1ld > measurementVectorHistory;
    actualStateVectorHistory[ initialTime ] = initialStateVector;
    for( unsigned int i = 0; i < numberOfTimeSteps; i++ )
    {
        // Compute actual values and perturb them
        currentTime += timeStep;
        currentStateVector = stateFunction2( currentTime, currentStateVector, currentControlVector ) +
                unscentedFilter->produceSystemNoise( );
        currentMeasurementVector = measurementFunction2( currentTime, currentStateVector ) +
                unscentedFilter->produceMeasurementNoise( );
        actualStateVectorHistory[ currentTime ] = currentStateVector;
        measurementVectorHistory[ currentTime ] = currentMeasurementVector;

        // Update control class
        control->setControlVector( currentTime, unscentedFilter->getCurrentStateEstimate( ) );

        // Update filter
        unscentedFilter->updateFilter( currentTime, currentMeasurementVector );

        // Print progress
        if ( showProgress )
        {
            std::cout << "Time: " << currentTime << std::endl
                      << "Measurement: " << currentMeasurementVector.transpose( ) << std::endl
                      << "Estimated State: " << unscentedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
        }
    }

    // Check that final state is as expected
    for ( int i = 0; i < initialStateVector.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( unscentedFilter->getCurrentStateEstimate( )[ i ], 1e-16L );
    }

//    // Save actual state history
//    input_output::writeDataMapToTextFile( actualStateVectorHistory, "UKFActualStateHistory2.dat", "/Users/Michele/Desktop/KF" );

//    // Save estimated state history
//    input_output::writeDataMapToTextFile( unscentedFilter->getEstimatedStateHistory( ), "UKFEstimatedStateHistory2.dat",
//                                          "/Users/Michele/Desktop/KF" );

//    // Save measurement history
//    input_output::writeDataMapToTextFile( measurementVectorHistory, "UKFMeasurementHistory2.dat", "/Users/Michele/Desktop/KF" );

//    // Save noise histories
//    std::pair< std::vector< Eigen::VectorXld >, std::vector< Eigen::VectorXld > > noiseHistory = unscentedFilter->getNoiseHistory( );
//    Eigen::MatrixXld systemNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.first ).transpose( );
//    Eigen::MatrixXld measurementNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.second ).transpose( );
//    input_output::writeMatrixToFile( systemNoise, "systemNoise2.dat", 16, "/Users/Michele/Desktop/KF" );
//    input_output::writeMatrixToFile( measurementNoise, "measurementNoise2.dat", 16, "/Users/Michele/Desktop/KF" );
}

// Constant parameters for example
const double ballisticCoefficient = 500.0;
const double gravitationalParameter = 32.2;

// Functions for extended Kalman filter.
Eigen::Vector2d stateFunction3( const double time, const Eigen::Vector2d& state, const Eigen::Vector2d& control )
{
    Eigen::Vector2d stateDerivative;
    stateDerivative[ 0 ] = state[ 1 ];
    stateDerivative[ 1 ] = 0.0034 * gravitationalParameter * std::exp( - state[ 0 ] / 22000.0 ) *
            std::pow( state[ 1 ], 2 ) / ( 2.0 * ballisticCoefficient ) - gravitationalParameter;
    return stateDerivative;
}
Eigen::Vector1d measurementFunction3( const double time, const Eigen::Vector2d& state )
{
    Eigen::Vector1d measurement;
    measurement[ 0 ] = state[ 0 ];
    return measurement;
}

// Test implementation of extended Kalman filter class.
BOOST_AUTO_TEST_CASE( testExtendedKalmanFilterSecondCase )
{
    using namespace tudat::filters;

    // Set initial conditions
    const double initialTime = 0;
    const double timeStep = 0.1;
    const unsigned int numberOfTimeSteps = 300;

    Eigen::Vector2d initialStateVector;
    initialStateVector[ 0 ] = 200000.0;
    initialStateVector[ 1 ] = -6000.0;

    Eigen::Vector2d initialEstimatedStateVector;
    initialEstimatedStateVector[ 0 ] = 200025.0;
    initialEstimatedStateVector[ 1 ] = -6150.0;

    Eigen::Matrix2d initialEstimatedStateCovarianceMatrix = Eigen::Matrix2d::Zero( );
    initialEstimatedStateCovarianceMatrix( 0, 0 ) = std::pow( 1000.0, 2 );
    initialEstimatedStateCovarianceMatrix( 1, 1 ) = 20000.0;

    // Set system and measurement uncertainty
    Eigen::Matrix2d systemUncertainty = Eigen::Matrix2d::Zero( );
    Eigen::Vector1d measurementUncertainty = Eigen::Vector1d::Zero( );
//    systemUncertainty( 0, 0 ) = 0;
//    systemUncertainty( 1, 1 ) = 0;
    measurementUncertainty[ 0 ] = std::pow( 1000.0, 2 );

    // Set integrator settings
    boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            boost::make_shared< numerical_integrators::IntegratorSettings< > > (
                numerical_integrators::euler, initialTime, timeStep );

    // Create control class
    boost::shared_ptr< ControlWrapper< double, double, 2 > > control =
            boost::make_shared< ControlWrapper< double, double, 2 > >( boost::lambda::constant( Eigen::Vector2d::Zero( ) ) );

    // Create extended Kalman filter object
    UnscentedKalmanFilterDoublePointer unscentedFilter = boost::make_shared< UnscentedKalmanFilterDouble >(
                boost::bind( &stateFunction3, _1, _2, control->getControlVector( ) ),
                boost::bind( &measurementFunction3, _1, _2 ),
                systemUncertainty, measurementUncertainty,
                initialTime, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix,
                integratorSettings );

    // Loop over each time step
    const bool showProgress = false;
    double currentTime = initialTime;
    Eigen::Vector2d currentStateVector = initialStateVector;
    Eigen::Vector2d currentControlVector = Eigen::Vector2d::Zero( );
    Eigen::Vector1d currentMeasurementVector;
    std::map< double, Eigen::Vector2d > actualStateVectorHistory;
    std::map< double, Eigen::Vector1d > measurementVectorHistory;
    actualStateVectorHistory[ initialTime ] = initialStateVector;
    for( unsigned int i = 0; i < numberOfTimeSteps; i++ )
    {
        // Compute actual values and perturb them
        currentTime += timeStep;
        currentStateVector += ( stateFunction3( currentTime, currentStateVector, currentControlVector ) +
                                unscentedFilter->produceSystemNoise( ) ) * timeStep;
        currentMeasurementVector = measurementFunction3( currentTime, currentStateVector ) +
                unscentedFilter->produceMeasurementNoise( );
        actualStateVectorHistory[ currentTime ] = currentStateVector;
        measurementVectorHistory[ currentTime ] = currentMeasurementVector;

        // Update control class
        control->setControlVector( currentTime, unscentedFilter->getCurrentStateEstimate( ) );

        // Update filter
        unscentedFilter->updateFilter( currentTime, currentMeasurementVector );

        // Print progress
        if ( showProgress )
        {
            std::cout << "Time: " << currentTime << std::endl
                      << "Measurement: " << currentMeasurementVector.transpose( ) << std::endl
                      << "Estimated State: " << unscentedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
        }
    }

    // Check that final state is as expected
    Eigen::Vector2d expectedFinalState = Eigen::Vector2d::Zero( );
//    expectedFinalState << 4.9717991447597925, -19.824346481395978;
    for ( int i = 0; i < expectedFinalState.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( unscentedFilter->getCurrentStateEstimate( )[ i ] - expectedFinalState[ i ],
                           std::numeric_limits< double >::epsilon( ) );
    }

    // Check that noise is actually normally distributed (within 5 %)
    std::pair< std::vector< Eigen::VectorXd >, std::vector< Eigen::VectorXd > > noiseHistory = unscentedFilter->getNoiseHistory( );
    Eigen::MatrixXd systemNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.first );
    Eigen::MatrixXd measurementNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.second );
    for ( unsigned int i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( statistics::computeStandardDeviationOfVectorComponents( systemNoise.row( i ) ),
                                    std::sqrt( systemUncertainty( i, i ) ), 5e-2 );
    }
    BOOST_CHECK_CLOSE_FRACTION( statistics::computeStandardDeviationOfVectorComponents( measurementNoise.row( 0 ) ),
                                std::sqrt( measurementUncertainty( 0, 0 ) ), 5e-2 );

    // Save actual state history
    input_output::writeDataMapToTextFile( actualStateVectorHistory,
                                          "UKFActualStateHistory.dat", "/Users/Michele/Desktop/KFBook" );

    // Save estimated state history
    input_output::writeDataMapToTextFile( unscentedFilter->getEstimatedStateHistory( ),
                                          "UKFEstimatedStateHistory.dat", "/Users/Michele/Desktop/KFBook" );

    // Save measurement history
    input_output::writeDataMapToTextFile( measurementVectorHistory, "UKFMeasurementHistory.dat", "/Users/Michele/Desktop/KFBook" );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
