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

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/External/SofaInterface/earthOrientation.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::sofa_interface;

BOOST_AUTO_TEST_SUITE( test_sofa_earth_orientation )


BOOST_AUTO_TEST_CASE( testPrecessionNutation )
{

    double testJulianDay1 = 2400000.5;
    double testJulianDay2 = 54195.500754444444444 * physical_constants::JULIAN_DAY;

    std::pair< Eigen::Vector2d, double > cipInGcrs = getPositionOfCipInGcrs( testJulianDay2, testJulianDay1, iau_2006 );

    // Data from SOFA cookbook.
    Eigen::Vector2d expectedPolePosition;
    expectedPolePosition.x( ) = 0.000712264729525;
    expectedPolePosition.y( ) = 0.000044385248875;

    double expectedCioLocator = -0.002200475;

    double arcSecondToRadian = 4.848136811095359935899141E-6;
    double dXTest = 0.0001750 * arcSecondToRadian;
    double dYTest = -0.0002259 * arcSecondToRadian;

    BOOST_CHECK_CLOSE_FRACTION( expectedPolePosition.x( ), cipInGcrs.first.x( ) + dXTest, 1.0E-9 );
    BOOST_CHECK_CLOSE_FRACTION( expectedPolePosition.y( ), cipInGcrs.first.y( ) + dYTest, 1.0E-7 );
    BOOST_CHECK_CLOSE_FRACTION( expectedCioLocator * arcSecondToRadian, cipInGcrs.second, 1.0E-8 );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat


