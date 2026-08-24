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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/math/integrators/rungeKuttaCoefficients.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_runge_kutta_coefficients )

using numerical_integrators::CoefficientSets;
using numerical_integrators::RungeKuttaCoefficients;

void checkValidityOfCoefficientSet( const CoefficientSets& coefficientSet, const double tolerance )
{
    // Declare coefficient set.
    RungeKuttaCoefficients coefficients;
    coefficients = coefficients.get( coefficientSet );

    // Check that the sum of the b-coefficients for both the integrated order and the
    // error-checking order is one.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            Eigen::VectorXhps::Constant( 2, HighPrecisionStateScalar( 1 ) ), coefficients.bCoefficients.rowwise( ).sum( ), tolerance );

    // Check that the first c-coefficient is zero.
    BOOST_CHECK_SMALL( coefficients.cCoefficients( 0 ), HighPrecisionStateScalar( tolerance ) );

    // Check that the c-coefficient/a-coefficient relation holds.
    for( int i = 1; i < coefficients.cCoefficients.size( ); i++ )
    {
        using std::abs;
        if( abs( coefficients.cCoefficients( i ) ) < tolerance )
        {
            BOOST_CHECK_SMALL( coefficients.aCoefficients.row( i ).sum( ), HighPrecisionStateScalar( tolerance ) );
        }

        else
        {
            BOOST_CHECK_CLOSE_FRACTION( coefficients.cCoefficients( i ), coefficients.aCoefficients.row( i ).sum( ), tolerance );
        }
    }
}

#if TUDAT_HIGH_PRECISION_STATE_SCALAR_IS_CPP_BIN_FLOAT_QUAD
BOOST_AUTO_TEST_CASE( testCoefficientsAreConstructedBeyondDoublePrecision )
{
    const RungeKuttaCoefficients& rk4 = RungeKuttaCoefficients::get( CoefficientSets::rungeKutta4Classic );
    BOOST_CHECK_EQUAL( rk4.bCoefficients( 0, 0 ), HighPrecisionStateScalar( 1 ) / HighPrecisionStateScalar( 6 ) );
    BOOST_CHECK_NE( rk4.bCoefficients( 0, 0 ), static_cast< HighPrecisionStateScalar >( static_cast< double >( 1.0 / 6.0 ) ) );

    const RungeKuttaCoefficients& rkf45 = RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg45 );
    BOOST_CHECK_EQUAL( rkf45.aCoefficients( 3, 0 ), HighPrecisionStateScalar( 1932 ) / HighPrecisionStateScalar( 2197 ) );
    BOOST_CHECK_NE( rkf45.aCoefficients( 3, 0 ), static_cast< HighPrecisionStateScalar >( static_cast< double >( 1932.0 / 2197.0 ) ) );

    const RungeKuttaCoefficients& feagin108 = RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFeagin108 );
    const HighPrecisionStateScalar expectedDecimal( "-0.915176561375291440520015019275342154318951387664369720564660" );
    BOOST_CHECK_EQUAL( feagin108.aCoefficients( 2, 0 ), expectedDecimal );
    BOOST_CHECK_NE( feagin108.aCoefficients( 2, 0 ), static_cast< HighPrecisionStateScalar >( static_cast< double >( expectedDecimal ) ) );
}
#endif

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
