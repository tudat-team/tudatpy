/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *      Montenbruck, O., Gill, E. Satellite Orbits: Models, Methods, Applications, Springer, 2005.
 *
 *    Notes
 *      There might be a problem with the RKF78 and DOPRI8 integrators, as the coefficients do not
 *      meet the required conditions to the tolerance achieved for all the other Runge-Kutta-type
 *      integrators tested (1.0e-14 versus 1.0e-15 for the other integrators). This should be
 *      looked into further to ensure that there are no bugs introduced in the coefficients
 *      used.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <atomic>
#include <thread>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/math/integrators/rungeKuttaCoefficients.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_runge_kutta_coefficients )

using numerical_integrators::CoefficientSets;
using numerical_integrators::RungeKuttaCoefficients;

BOOST_AUTO_TEST_CASE( testConcurrentDormandPrinceCoefficientLookup )
{
    // Exercise cold initialization and repeated copies without the arc executor.
    // DP87 used to rewrite a shared table on every lookup (testing order 7
    // while initializing order 8), racing with readers copying its matrices.
    constexpr int workers = 8;
    std::atomic< int > ready{ 0 };
    std::atomic< bool > start{ false }, valid{ true };
    std::vector< RungeKuttaCoefficients > snapshots( workers );
    std::vector< std::thread > threads;
    for( int worker = 0; worker < workers; ++worker )
    {
        threads.emplace_back( [ &, worker ] {
            ready.fetch_add( 1 );
            while( !start.load( ) )
            {
                std::this_thread::yield( );
            }
            for( int iteration = 0; iteration < 1000; ++iteration )
            {
                auto coefficients = RungeKuttaCoefficients::get( CoefficientSets::rungeKutta87DormandPrince );
                if( coefficients.higherOrder != 8 || coefficients.lowerOrder != 7 || coefficients.aCoefficients.rows( ) != 13 ||
                    coefficients.aCoefficients.cols( ) != 12 || !coefficients.aCoefficients.allFinite( ) ||
                    coefficients.aCoefficients( 1, 0 ) != 1.0 / 18.0 || !coefficients.bCoefficients.allFinite( ) ||
                    !coefficients.cCoefficients.allFinite( ) )
                {
                    valid.store( false );
                }
                snapshots.at( worker ) = std::move( coefficients );
            }
        } );
    }
    while( ready.load( ) != workers )
    {
        std::this_thread::yield( );
    }
    start.store( true );
    for( auto& thread : threads )
    {
        thread.join( );
    }
    BOOST_REQUIRE( valid.load( ) );
    const auto& reference = RungeKuttaCoefficients::get( CoefficientSets::rungeKutta87DormandPrince );
    for( const auto& snapshot : snapshots )
    {
        BOOST_CHECK( snapshot.aCoefficients.isApprox( reference.aCoefficients, 0.0 ) );
        BOOST_CHECK( snapshot.bCoefficients.isApprox( reference.bCoefficients, 0.0 ) );
        BOOST_CHECK( snapshot.cCoefficients.isApprox( reference.cCoefficients, 0.0 ) );
    }
}

void checkValidityOfCoefficientSet( const CoefficientSets& coefficientSet, const double tolerance )
{
    // Declare coefficient set.
    RungeKuttaCoefficients coefficients;
    coefficients = coefficients.get( coefficientSet );

    // Check that the sum of the b-coefficients for both the integrated order and the
    // error-checking order is one.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( Eigen::VectorXd::Constant( 2, 1.0 ), coefficients.bCoefficients.rowwise( ).sum( ), tolerance );

    // Check that the first c-coefficient is zero.
    BOOST_CHECK_SMALL( coefficients.cCoefficients( 0 ), tolerance );

    // Check that the c-coefficient/a-coefficient relation holds.
    for( int i = 1; i < coefficients.cCoefficients.size( ); i++ )
    {
        if( std::fabs( coefficients.cCoefficients( i ) ) < tolerance )
        {
            BOOST_CHECK_SMALL( coefficients.aCoefficients.row( i ).sum( ), tolerance );
        }

        else
        {
            BOOST_CHECK_CLOSE_FRACTION( coefficients.cCoefficients( i ), coefficients.aCoefficients.row( i ).sum( ), tolerance );
        }
    }
}

BOOST_AUTO_TEST_CASE( testRungeKuttaFehlberg45Coefficients )
{
    // Check validity of Runge-Kutta-Fehlberg 45 coefficients.
    checkValidityOfCoefficientSet( CoefficientSets::rungeKuttaFehlberg45, 1.0e-15 );
}

BOOST_AUTO_TEST_CASE( testRungeKuttaFehlberg78Coefficients )
{
    // Check validity of Runge-Kutta-Fehlberg 78 coefficients.
    // Note, for some reason, the RKF78 set fails the unit test when the tolerance in
    // checkValidityOfCoefficientSet() is set to lower than this value (rows 8 and 9 of
    // aCoefficients matrix sum does not correspond to cCoefficient counterpart with tolerance less
    // than 1.0e-14).
    checkValidityOfCoefficientSet( CoefficientSets::rungeKuttaFehlberg78, 1.0e-14 );
}

BOOST_AUTO_TEST_CASE( testRungeKutta87DormandAndPrinceCoefficients )
{
    // Check validity of Runge-Kutta 87 (Dormand and Prince) coefficients.
    // Note, for some reason, the DOPRI8 set fails the unit test when the tolerance in
    // checkValidityOfCoefficientSet() is set to lower than this value (row 10 of aCoefficients
    // matrix sum does not correspond to cCoefficient counterpart with tolerance less than
    // 1.0e-14).
    checkValidityOfCoefficientSet( CoefficientSets::rungeKutta87DormandPrince, 1.0e-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
