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

#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Mathematics/Filters/unscentedKalmanFilter.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_unscented_kalman_filter )

// Functions for extended Kalman filter
Eigen::Vector2d stateFunction1( const double time, const Eigen::Vector2d& state, const Eigen::Vector2d& control )
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
    const unsigned int numberOfTimeSteps = 5;

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

    // Create extended Kalman filter object
    KalmanFilterDoublePointer unscentedFilter = boost::make_shared< UnscentedKalmanFilterDouble >(
                boost::bind( &stateFunction1, _1, _2, _3 ),
                boost::bind( &measurementFunction1, _1, _2 ),
                systemUncertainty, measurementUncertainty,
                initialTime, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix,
                integratorSettings );

    // Loop over each time step
    const bool showProgress = true;
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
        actualStateVectorHistory[ currentTime ] = currentStateVector +
                stateFunction1( currentTime, currentStateVector, currentControlVector ) * timeStep;
        currentStateVector += ( stateFunction1( currentTime, currentStateVector, currentControlVector ) +
                                unscentedFilter->produceSystemNoise( ) ) * timeStep;
        currentMeasurementVector = measurementFunction1( currentTime, currentStateVector ) +
                unscentedFilter->produceMeasurementNoise( );
        measurementVectorHistory[ currentTime ] = currentMeasurementVector;

        // Update filter
        unscentedFilter->updateFilter( currentTime, currentControlVector, currentMeasurementVector );

        // Print progress
        if ( showProgress )
        {
            std::cout << "Time: " << currentTime << std::endl
                      << "Measurement: " << currentMeasurementVector.transpose( ) << std::endl
                      << "Estimated State: " << unscentedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
        }
    }

//    Eigen::Vector2d expectedFinalState;
//    expectedFinalState << 5.226210753689565, -11.602113592960805;
//    for ( int i = 0; i < expectedFinalState.rows( ); i++ )
//    {
//        BOOST_CHECK_SMALL( unscentedFilter->getCurrentStateEstimate( )[ i ] - expectedFinalState[ i ],
//                           std::numeric_limits< double >::epsilon( ) );
//    }

    // Save state history
    input_output::writeDataMapToTextFile( actualStateVectorHistory,
                                          "UKFActualStateHistory.dat",
                                          "/Users/Michele/Desktop" );

    // Save state history
    input_output::writeDataMapToTextFile( unscentedFilter->getEstimatedStateHistory( ),
                                          "UKFEstimatedStateHistory.dat",
                                          "/Users/Michele/Desktop" );

    // Save state history
    input_output::writeDataMapToTextFile( measurementVectorHistory,
                                          "UKFMeasurementHistory.dat",
                                          "/Users/Michele/Desktop" );
}

//// Functions for extended Kalman filter
//Eigen::Vector3d stateFunction2( const double time, const Eigen::Vector3d& state, const Eigen::Vector3d& control )
//{
//    Eigen::Vector3d stateFunction;
//    stateFunction[ 0 ] = state[ 1 ];
//    stateFunction[ 1 ] = state[ 2 ];
//    stateFunction[ 2 ] = 0.05 * state[ 0 ] * ( state[ 1 ] + state[ 2 ] );
//    return stateFunction;
//}
//Eigen::Vector1d measurementFunction2( const double time, const Eigen::Vector3d& state )
//{
//    Eigen::Vector1d measurement;
//    measurement[ 0 ] = state[ 0 ];
//    return measurement;
//}

//// Test implementation of unscented Kalman filter class. Tested by comparison with code by
//// https://de.mathworks.com/matlabcentral/fileexchange/18217-learning-the-unscented-kalman-filter
//BOOST_AUTO_TEST_CASE( testUnscentedKalmanFilterSecondCase )
//{
//    using namespace tudat::filters;

//    // Set initial conditions
//    const double initialTime = 0;
//    const double timeStep = 1.0;
//    const unsigned int numberOfTimeSteps = 100;

//    Eigen::Vector3d initialStateVector = Eigen::Vector3d::Zero( );
//    initialStateVector[ 2 ] = 1.0;

//    Eigen::Vector3d initialEstimatedStateVector;
//    initialEstimatedStateVector[ 0 ] = 0.05376671395461;
//    initialEstimatedStateVector[ 1 ] = 0.183388501459509;
//    initialEstimatedStateVector[ 2 ] = 0.774115313899635;

//    Eigen::Matrix3d initialEstimatedStateCovarianceMatrix = Eigen::Matrix3d::Identity( );

//    // Set system and measurement uncertainty
//    Eigen::Matrix3d systemUncertainty = 0.01 * Eigen::Matrix3d::Identity( );
//    Eigen::Vector1d measurementUncertainty = 0.01 * Eigen::Vector1d::Identity( );

//    // Set null integrator settings
//    boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings = NULL;

//    // Create extended Kalman filter object
//    KalmanFilterDoublePointer unscentedFilter = boost::make_shared< UnscentedKalmanFilterDouble >(
//                boost::bind( &stateFunction2, _1, _2, _3 ),
//                boost::bind( &measurementFunction2, _1, _2 ),
//                systemUncertainty, measurementUncertainty,
//                initialTime, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix,
//                integratorSettings,
//                custom_parameters, std::make_pair( 0.003, 0.0 ) );

//    // Loop over each time step
//    const bool showProgress = false;
//    double currentTime = initialTime;
//    Eigen::Vector3d currentStateVector = initialStateVector;
//    Eigen::Vector3d currentControlVector = Eigen::Vector3d::Zero( );
//    Eigen::Vector1d currentMeasurementVector;
//    std::map< double, Eigen::Vector3d > actualStateVectorHistory;
//    std::map< double, Eigen::Vector1d > measurementVectorHistory;
//    actualStateVectorHistory[ initialTime ] = initialStateVector;
//    for( unsigned int i = 0; i < numberOfTimeSteps; i++ )
//    {
//        // Compute actual values and perturb them
//        currentTime += timeStep;
//        actualStateVectorHistory[ currentTime ] = currentStateVector +
//                stateFunction2( currentTime, currentStateVector, currentControlVector ) * timeStep;
//        currentStateVector += ( stateFunction2( currentTime, currentStateVector, currentControlVector ) +
//                                unscentedFilter->produceSystemNoise( ) ) * timeStep;
//        currentMeasurementVector = measurementFunction2( currentTime, currentStateVector ) +
//                unscentedFilter->produceMeasurementNoise( );
//        measurementVectorHistory[ currentTime ] = currentMeasurementVector;

//        // Update filter
//        unscentedFilter->updateFilter( currentTime, currentControlVector, currentMeasurementVector );

//        // Print progress
//        if ( showProgress )
//        {
//            std::cout << "Time: " << currentTime << std::endl
//                      << "Measurement: " << currentMeasurementVector.transpose( ) << std::endl
//                      << "Estimated State: " << unscentedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
//        }
//    }

////    Eigen::Vector3d expectedFinalState;
////    expectedFinalState << 5.226210753689565, -11.602113592960805;
////    for ( int i = 0; i < expectedFinalState.rows( ); i++ )
////    {
////        BOOST_CHECK_SMALL( unscentedFilter->getCurrentStateEstimate( )[ i ] - expectedFinalState[ i ],
////                           std::numeric_limits< double >::epsilon( ) );
////    }

////    // Save state history
////    input_output::writeDataMapToTextFile( actualStateVectorHistory,
////                                          "UKFActualStateHistory.dat",
////                                          "/Users/Michele/Desktop" );

////    // Save state history
////    input_output::writeDataMapToTextFile( unscentedFilter->getEstimatedStateHistory( ),
////                                          "UKFEstimatedStateHistory.dat",
////                                          "/Users/Michele/Desktop" );

////    // Save state history
////    input_output::writeDataMapToTextFile( measurementVectorHistory,
////                                          "UKFMeasurementHistory.dat",
////                                          "/Users/Michele/Desktop" );
//}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
