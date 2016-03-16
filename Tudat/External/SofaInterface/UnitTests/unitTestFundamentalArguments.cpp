/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/External/SofaInterface/fundamentalArguments.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::sofa_interface;

BOOST_AUTO_TEST_SUITE( test_sofa_fundamental_arguments )

//! Test calculation of Delaunay/Doodson arguments used for calculation of Earth-tide effects
BOOST_AUTO_TEST_CASE( testSofaFundamentalArguments )
{

    // Test whether Doodson <-> Delaunay argument conversion matrices are each other's inverse.
    Eigen::Matrix< double, 5, 5 > conversionMatrixComparison =
            delaunayToDoodsonArguments - doodsonToDelaunayArguments.inverse( );
    for( unsigned int i = 0; i < 5; i++ )
    {
        for( unsigned int j = 0; j < 5; j++ )
        {
            BOOST_CHECK_EQUAL( conversionMatrixComparison( i, j ), 0.0 );
        }
    }

    // Use test case from FUNDARG.F code (provided with IERS 2010 conventions) to calculate Delaunay arguments
    double testModifiedJulianDay1 = 54465.0;
    double testSecondsSinceJ2000 = ( testModifiedJulianDay1 -
            ( - basic_astrodynamics::JULIAN_DAY_AT_0_MJD + basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) *
            physical_constants::JULIAN_DAY;
    Eigen::Matrix< double, 5, 1 > expectedFundamentalArgumentValues;
    expectedFundamentalArgumentValues <<2.291187512612069099, 6.212931111003726414, 3.658025792050572989,
            4.554139562402433228, -0.5167379217231804489;

    // Calculate Delaunay arguments.
    Eigen::Matrix< double, 5, 1 > fundamentalArgumentValues =
            calculateDelaunayFundamentalArguments( testSecondsSinceJ2000 );

    // Compare against IERS results.
    for( unsigned int i = 0; i < 5; i++ )
    {
        BOOST_CHECK_SMALL( expectedFundamentalArgumentValues( i ) - fundamentalArgumentValues( i ), 1.0E-13  );
    }

    // Calculate Doodson arguments directly and from Delaunay arguments and compare.
    Eigen::Matrix< double, 6, 1 > doodsonArguments = calculateDoodsonFundamentalArguments( testSecondsSinceJ2000 );
    Eigen::Matrix< double, 5, 1 > reconstructedFundamentalArguments =
            doodsonToDelaunayArguments * doodsonArguments.segment( 1, 5 );
    for( unsigned int i = 0; i < 5; i++ )
    {
        BOOST_CHECK_SMALL( fundamentalArgumentValues( i ) - reconstructedFundamentalArguments( i ), 1.0E-15  );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat



