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

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"
#include "Tudat/Mathematics/Interpolators/jumpDataLinearInterpolator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::earth_orientation;

BOOST_AUTO_TEST_SUITE( test_eop_reader )

// Compare to SOFA Cookbook
BOOST_AUTO_TEST_CASE( testDataAgainstSofaCookbookExamples )
{
    using namespace tudat::interpolators;

    std::cout<< tudat::input_output::getEarthOrientationDataFilesPath( )<<std::endl;

    EOPReader eopReader = EOPReader( tudat::input_output::getEarthOrientationDataFilesPath( ) + "eopc04_08_IAU2000.62-now",
            "C04", iau_2000_a );


    int year = 2007;

    int month = 4;
    int day = 5;
    int hour = 0;
    int minute = 0;
    double seconds = 0.0;

    double jd1, jd2;

    iauDtf2d ( ( char* )std::string( "UTC" ).c_str( ), year, month, day, hour, minute, seconds, &jd1, &jd2 );

    std::cout<<jd1<<" "<<jd2<<std::endl;
    double utcSecondsSinceJ2000 = ( ( jd1 - ( basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) + jd2 ) * physical_constants::JULIAN_DAY;

    std::cout<<utcSecondsSinceJ2000<<std::endl;

    double arcSecondToRadian = 4.848136811095359935899141E-6;

    double expectedXp = 0.033227 * arcSecondToRadian;
    double expectedYp = 0.483135 * arcSecondToRadian;
    double expectedUT1COffset = -0.0714209;
    double expecteddX = 0.000247 * arcSecondToRadian;
    double expecteddY = -0.000280 * arcSecondToRadian;

    boost::shared_ptr< CubicSplineInterpolator< double, Eigen::Vector2d > > cipInItrsInterpolator =
            boost::make_shared< CubicSplineInterpolator< double, Eigen::Vector2d > >(
                eopReader.getCipInItrsMapInSecondsSinceJ2000( ) ); // xPole, yPole

    Eigen::Vector2d cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );

    boost::shared_ptr< CubicSplineInterpolator< double, Eigen::Vector2d > > cipInGcrsCorrectionInterpolator =
            boost::make_shared< CubicSplineInterpolator< double, Eigen::Vector2d > >(
                eopReader.getCipInGcrsCorrectionMapInSecondsSinceJ2000( ) ); // dX, dY

    Eigen::Vector2d cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );

    boost::shared_ptr< JumpDataLinearInterpolator < double, double > > ut1MinusUtcInterpolator =
            boost::make_shared< JumpDataLinearInterpolator< double, double > >(
                eopReader.getUt1MinusUtcMapInSecondsSinceJ2000( ), 0.5, 1.0 ); // d(UT1-UTC)

    double utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1, expectedUT1COffset, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), expectedXp, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), expectedYp, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), expecteddX, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), expecteddY, std::numeric_limits< double >::epsilon( ) );

    day = 6;

    double expectedXp2 = 0.035739  * arcSecondToRadian;
    double expectedYp2 = 0.484209 * arcSecondToRadian;
    double expectedUT1COffset2 = -0.0727562;
    double expecteddX2 = 0.000252 * arcSecondToRadian;
    double expecteddY2 = -0.000294 * arcSecondToRadian;

    iauDtf2d ( ( char* )std::string( "UTC" ).c_str( ), year, month, day, hour, minute, seconds, &jd1, &jd2 );

    utcSecondsSinceJ2000 = ( ( jd1 - ( basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) + jd2 ) * physical_constants::JULIAN_DAY;

    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );

    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );

    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1, expectedUT1COffset2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), expectedXp2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), expectedYp2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), expecteddX2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), expecteddY2, std::numeric_limits< double >::epsilon( ) );

    day = 5;
    hour = 12;

    iauDtf2d ( ( char* )std::string( "UTC" ).c_str( ), year, month, day, hour, minute, seconds, &jd1, &jd2 );

    utcSecondsSinceJ2000 = ( ( jd1 - ( basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) + jd2 ) * physical_constants::JULIAN_DAY;

    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );

    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );

    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1, ( expectedUT1COffset + expectedUT1COffset2 ) / 2.0, 1.0E-2 );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), ( expectedXp + expectedXp2 ) / 2.0, 1.0E-2 );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), ( expectedYp + expectedYp2 ) / 2.0, 1.0E-2 );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), ( expecteddX + expecteddX2 ) / 2.0, 1.0E-2 );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), ( expecteddY + expecteddY2 ) / 2.0, 1.0E-2 );
}

BOOST_AUTO_TEST_CASE( testLeapSecondIdentification )
{

    using namespace tudat::interpolators;

    EOPReader eopReader = EOPReader( tudat::input_output::getEarthOrientationDataFilesPath( ) + "eopc04_08_IAU2000.62-now",
            "C04", iau_2000_a );

    boost::shared_ptr< JumpDataLinearInterpolator < double, double > > ut1MinusUtcInterpolator =
            boost::make_shared< JumpDataLinearInterpolator< double, double > >(
                eopReader.getUt1MinusUtcMapInSecondsSinceJ2000( ), 0.5, 1.0 );

    int year, month, day;
    int hour = 0;
    int minute = 0;
    double seconds = 0.0;

    int numberOfLeapSeconds = 24;

    Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > leapSecondDays;
    leapSecondDays.resize( numberOfLeapSeconds, 3 );
    leapSecondDays<<1, 7, 1972,
    1, 1, 1973,
    1, 1, 1974,
    1, 1, 1975,
    1, 1, 1976,
    1, 1, 1977,
    1, 1, 1978,
    1, 1, 1979,
    1, 1, 1980,
    1, 7, 1981,
    1, 7, 1982,
    1, 7, 1983,
    1, 7, 1985,
    1, 1, 1988,
    1, 1, 1990,
    1, 1, 1991,
    1, 7, 1992,
    1, 7, 1993,
    1, 7, 1994,
    1, 1, 1996,
    1, 7, 1997,
    1, 1, 1999,
    1, 1, 2006,
    1, 1, 2009;

    double differenceMinusOneDay, differenceMinusHalfDay, differenceAtEpoch, differenceAtEpochPlusSecond, differencePlusOneDay;

    double jd1, jd2;
    for( int i = 0; i < numberOfLeapSeconds; i++ )
    {
        year = leapSecondDays( i, 2 );
        month = leapSecondDays( i, 1 );
        day = leapSecondDays( i, 0 );
        iauDtf2d( ( char* )std::string( "UTC" ).c_str( ), year, month, day, hour, minute, seconds, &jd1, &jd2 );
        double utcSecondsSinceJ2000 = ( jd1 + jd2 - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;

        differenceMinusOneDay = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 - physical_constants::JULIAN_DAY );
        differenceMinusHalfDay = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 - 0.5 * physical_constants::JULIAN_DAY );
        differenceAtEpoch = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 + 1.0E-6 );
        differenceAtEpochPlusSecond = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 + 1.0 );
        differencePlusOneDay = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 + physical_constants::JULIAN_DAY );
        std::cout<<differenceMinusOneDay<<" "<<differenceMinusHalfDay<<" "<<differenceAtEpoch<<" "<<differenceAtEpochPlusSecond<<" "<<differencePlusOneDay<<std::endl;
        BOOST_CHECK_CLOSE_FRACTION( differenceMinusOneDay, differenceMinusHalfDay, 0.01 );
        BOOST_CHECK_CLOSE_FRACTION( differenceMinusOneDay + 1.0, differenceAtEpoch, 0.01);
        BOOST_CHECK_CLOSE_FRACTION( differenceAtEpochPlusSecond, differenceAtEpoch, 0.01 );
        BOOST_CHECK_CLOSE_FRACTION( differencePlusOneDay, differenceAtEpoch, 0.01 );

    }


}



BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




