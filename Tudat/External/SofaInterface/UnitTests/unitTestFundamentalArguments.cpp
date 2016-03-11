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


BOOST_AUTO_TEST_CASE( testSofaFundamentalArguments )
{

    double testModifiedJulianDay1 = 54465.0;
    double testSecondsSinceJ2000 = ( testModifiedJulianDay1 -
            ( - basic_astrodynamics::JULIAN_DAY_AT_0_MJD + basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) * physical_constants::JULIAN_DAY;

    Eigen::Matrix< double, 5, 1 > expectedFundamentalArgumentValues;// Taken from IERS code FUNDARG.F
    expectedFundamentalArgumentValues <<2.291187512612069099, 6.212931111003726414, 3.658025792050572989, 4.554139562402433228, -0.5167379217231804489;

    Eigen::Matrix< double, 5, 1 > fundamentalArgumentValues = calculateDelaunayFundamentalArguments( testSecondsSinceJ2000 );

    for( unsigned int i = 0; i < 5; i++ )
    {
        BOOST_CHECK_SMALL( expectedFundamentalArgumentValues( i ) - fundamentalArgumentValues( i ), 1.0E-13  );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat



