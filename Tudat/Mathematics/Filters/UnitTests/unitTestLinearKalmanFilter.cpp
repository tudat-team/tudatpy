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

#include <vector>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/Filters/linearKalmanFilter.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_kalman_filters )

// Test implementation of linear Kalman filter class. Tested by comparison with code by https://github.com/hmartiro/kalman-cpp.
BOOST_AUTO_TEST_CASE( testLinearKalmanFilter )
{
    using namespace tudat::filters;

    int n = 3; // Number of states
    int m = 1; // Number of measurements

    double dt = 1.0/30; // Time step

    Eigen::MatrixXd A(n, n); // System dynamics matrix
    Eigen::MatrixXd C(m, n); // Output matrix
    Eigen::MatrixXd Q(n, n); // Process noise covariance
    Eigen::MatrixXd R(m, m); // Measurement noise covariance
    Eigen::MatrixXd P(n, n); // Estimate error covariance

    // Discrete LTI projectile motion, measuring position only
    A << 1, dt, 0, 0, 1, dt, 0, 0, 1;
    C << 1, 0, 0;

    // Reasonable covariance matrices
    Q << .05, .05, .0, .05, .05, .0, .0, .0, 1e-30;
    R << 0.5;
    P << .1, .1, .1, .1, 10000, 10, .1, 10, 100;

    // List of noisy position measurements ( y )
    std::vector< double > measurements =
    {
        1.04202710058, 1.10726790452, 1.2913511148, 1.48485250951, 1.72825901034,
        1.74216489744, 2.11672039768, 2.14529225112, 2.16029641405, 2.21269371128,
        2.57709350237, 2.6682215744, 2.51641839428, 2.76034056782, 2.88131780617,
        2.88373786518, 2.9448468727, 2.82866600131, 3.0006601946, 3.12920591669,
        2.858361783, 2.83808170354, 2.68975330958, 2.66533185589, 2.81613499531,
        2.81003612051, 2.88321849354, 2.69789264832, 2.4342229249, 2.23464791825,
        2.30278776224, 2.02069770395, 1.94393985809, 1.82498398739, 1.52526230354,
        1.86967808173, 1.18073207847, 1.10729605087, 0.916168349913, 0.678547664519,
        0.562381751596, 0.355468474885, -0.155607486619, -0.287198661013, -0.602973173813
    };

    // Best guess of initial states
    Eigen::VectorXd x0( n );
    double t = 0;
    x0 << measurements[ 0 ], 0, -9.81;

    // Create linear filter
    KalmanFilterPointer linearFilter = boost::make_shared< LinearKalmanFilter >( A, Eigen::MatrixXd::Zero( n, n ), C,
                                                                                 Q, R, 0.0, x0, P );

    // Feed measurements into filter, output estimated states
    Eigen::VectorXd y( m );
//    std::cout << "t = " << t << ", " << "x_hat[ 0 ]: " << linearFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
    for( unsigned int i = 0; i < measurements.size( ); i++ ) //measurements.size( )
    {
        t += dt;
        y << measurements[ i ];
        linearFilter->updateFilter( t, Eigen::Vector3d::Zero( ), y );
//        std::cout << "t = " << t << ", " << "y[" << i << "] = " << y.transpose( )
//                  << ", x_hat[" << i << "] = " << linearFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
    }

    Eigen::Vector3d expectedState;
    expectedState << -0.34094280864427917, -8.2429633777065501, -9.7238568066459514;
    for ( int i = 0; i < n; i++ )
    {
        BOOST_CHECK_SMALL( linearFilter->getCurrentStateEstimate( )[ i ] - expectedState[ i ], std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
