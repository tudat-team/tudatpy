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
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Mathematics/Filters/linearKalmanFilter.h"
#include "Tudat/Mathematics/Filters/extendedKalmanFilter.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

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
    Q << .05, .05, .0, .05, .05, .0, .0, .0, .0;
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
    KalmanFilterDoublePointer linearFilter =
            boost::make_shared< LinearKalmanFilterDouble >( A, Eigen::MatrixXd::Zero( n, n ), C, Q, R, 0.0, x0, P );

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

    Eigen::Vector3d expectedFinalState;
    expectedFinalState << -0.19336039191479262, -9.8935932082956075, -12.072317473180874;
    for ( int i = 0; i < n; i++ )
    {
        BOOST_CHECK_SMALL( linearFilter->getCurrentStateEstimate( )[ i ] - expectedFinalState[ i ], std::numeric_limits< double >::epsilon( ) );
    }
}

// Functions for extended Kalman filter
Eigen::Vector2d stateFunction( const double time, const Eigen::Vector2d& state, const Eigen::Vector2d& control )
{
    Eigen::Vector2d stateDerivative;
    stateDerivative[ 0 ] = state[ 1 ] * std::pow( std::cos( state[ 0 ] ), 3 );
    stateDerivative[ 1 ] = std::sin( state[ 0 ] );
    return stateDerivative;
}

Eigen::Vector1d measurementFunction( const double time, const Eigen::Vector2d& state )
{
    Eigen::Vector1d measurement;
    measurement[ 0 ] = std::pow( state[ 0 ], 3 );
    return measurement;
}

Eigen::Matrix2d stateJacobian( const double time, const Eigen::Vector2d& state, const Eigen::Vector2d& control )
{
    Eigen::Matrix2d stateJacobian = Eigen::Matrix2d::Zero( );
    stateJacobian( 0, 0 ) = - 3 * state[ 1 ] * std::pow( std::cos( state[ 0 ] ), 2 ) * std::sin( state[ 0 ] );
    stateJacobian( 0, 1 ) = std::pow( std::cos( state[ 0 ] ), 3 );
    stateJacobian( 1, 0 ) = std::cos( state[ 0 ] );
    return stateJacobian;
}

Eigen::Vector1d measurementJacobian( const double time, const Eigen::Vector2d& state )
{
    Eigen::Vector1d measurementJacobian;
    measurementJacobian[ 0 ] = 3.0 * std::pow( state[ 0 ], 2 );
    return measurementJacobian;
}

// Test implementation of extended Kalman filter class.
BOOST_AUTO_TEST_CASE( testExtendedKalmanFilter )
{
    using namespace tudat::filters;

    int n = 2; // Number of states
    int m = 1; // Number of measurements

    double dt = 0.01; // Time step

    // List of noisy position measurements ( y )
    Eigen::Matrix2d Q = Eigen::Matrix2d::Zero( ); // Process noise covariance
    Eigen::Vector1d R = Eigen::Vector1d::Zero( ); // Measurement noise covariance
    Eigen::Matrix2d P = Eigen::Matrix2d::Zero( ); // Estimate error covariance

    Q( 0, 0 ) = 100;
    R[ 0 ] = 100;
    P( 0, 0 ) = 100;
    P( 1, 1 ) = 100;

    // Best guess of initial states
    Eigen::Vector2d x0;
    x0[ 0 ] = 3;
    x0[ 1 ] = -0.3;
    double t = 0;

    std::vector< double > measurements =
    {
        36.9631421707894,25.2103142487346,28.4679924446531,24.5255627189394,33.9972749762844,21.2950389112620,31.1125083310715,
        34.6123799381830,55.3536702554016,46.4992735617870,21.7907204826762,46.9137540172454,71.9683675148751,47.3154618450114,
        70.8158233889969,61.4018847088800,73.8621412126349,46.8559100444617,71.5820658501241,69.1415359818752,114.108014698029,
        86.4494318613861,95.9497299298922,81.1455020984245,90.0510647224007,98.6081489545672,117.130998452190,101.340296352233,
        113.097829781596,80.3590509358706,103.224986417841,90.9689024085437,76.7194761660697,92.6873424042008,74.1243570398958,
        79.2652543902807,67.3723807956460,87.8039465277618,87.6450751276546,71.6580311114125,74.2961157093809,70.1373803830775,
        56.8934450644952,73.1726938228388,63.1207937047202,61.4433036785002,58.7812082626495,68.2213105446153,59.4000726137473,
        95.3519784403932,85.9599163312728,80.9722450217846,74.1671625537973,60.7417016608408,78.8336472573626,75.2706723602010,
        65.3409957032209,87.9189903057066,70.9359651739979,73.2820033695259,70.3416030703707,64.9397584677356,65.1312944994787,
        107.662726736862,108.028877144171,95.0696371339727,70.5988638187488,70.3311297580453,71.4831400759894,94.1840011718024,
        69.3854828860938,63.7269132767590,71.3982604895745,94.5098410692506,90.5150426497712,83.1347403397833,69.7254000282504,
        75.3836501725876,67.8344172849509,80.1768363281065,65.5104076064073,85.2858749458157,73.3491014978694,87.7854310183607,
        91.8583591176269,100.861257455023,84.4335787967761,106.694885592917,102.032931982465,87.6326534568385,79.6996429555746,
        80.0597679796415,83.3614830263125,102.678649052518,98.6428319510951,107.590681109393,114.069759914872,91.5863301007616,
        82.2569449503335,80.6627197624853
    };

    // Integrator
    boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            boost::make_shared< numerical_integrators::IntegratorSettings< > > (
                numerical_integrators::euler, t, 1.0 );

    // Create extended filter
    KalmanFilterDoublePointer extendedFilter = boost::make_shared< ExtendedKalmanFilterDouble >(
                boost::bind( &stateFunction, _1, _2, _3 ),
                boost::bind( &measurementFunction, _1, _2 ),
                boost::bind( &stateJacobian, _1, _2, _3 ),
                boost::lambda::constant( Eigen::Vector1d::Ones( ) ),
                boost::bind( &measurementJacobian, _1, _2 ),
                boost::lambda::constant( Eigen::Vector1d::Zero( ) ),
                Q, R, 0.0, x0, P, true, integratorSettings );

    // Feed measurements into filter, output estimated states
    Eigen::Vector1d y;
    std::cout << "t = " << t << ", " << "x_hat[ 0 ]: " << extendedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
    for( unsigned int i = 0; i < measurements.size( ); i++ ) //
    {
        t += dt;
        y[ 0 ] = measurements.at( i );
        extendedFilter->updateFilter( t, Eigen::Vector1d::Zero( ), y );
        std::cout << "t = " << t << ", " << "y[" << i << "] = " << y.transpose( )
                  << ", x_hat[" << i << "] = " << extendedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl;
    }

//    Eigen::Vector1d expectedFinalState;
//    expectedFinalState << ;
//    for ( int i = 0; i < n; i++ )
//    {
//        BOOST_CHECK_SMALL( extendedFilter->getCurrentStateEstimate( )[ i ] - expectedFinalState[ i ],
//                           std::numeric_limits< double >::epsilon( ) );
//    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
