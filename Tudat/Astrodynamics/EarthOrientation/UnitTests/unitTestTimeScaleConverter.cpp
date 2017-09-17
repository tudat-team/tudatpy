/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <limits>
#include <string>

#include "Tudat/Basics/testMacros.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{
namespace unit_tests
{

struct SofaTimeOutput
{
    double expectedUtcFractionalDays;
    double expectedUtcDays;
    double expectedUt1Seconds;
    double expectedUt1Days;
    double expectedTaiFractionalDays;
    double expectedTaiDays;
    double expectedTtFractionalDays;
    double expectedTtDays;
    double expectedTdbFractionalDays;
    double expectedTdbDays;
};

SofaTimeOutput getSofaDirectTimes( )
{
    SofaTimeOutput sofaTimes;

    int latnd, latnm, lonwd, lonwm, j, iy, mo, id, ih, im, ihmsf[4];
    double slatn, slonw, hm, elon, phi, xyz[3], u, v, sec,
            utc1, utc2, dut, ut11, ut12, ut, tai1, tai2, tt1, tt2,
            tcg1, tcg2, dtr, tdb1, tdb2, tcb1, tcb2;
    /* Site terrestrial coordinates (WGS84). */
    latnd = 19;
    latnm = 28;
    slatn = 52.5;
    lonwd = 155;
    lonwm = 55;
    slonw = 59.6;
    hm = 0.0;

    /* Transform to geocentric. */
    j = iauAf2a ( '+', latnd, latnm, slatn, &phi );
    j = iauAf2a ( '-', lonwd, lonwm, slonw, &elon );
    j = iauGd2gc ( 1, elon, phi, hm, xyz );
    u = sqrt ( xyz[0]*xyz[0] + xyz[1]*xyz[1] );
    v = xyz[2];

    /* UTC date and time. */
    iy = 2006;
    mo = 1;
    id = 15;
    ih = 21;
    im = 24;
    sec = 37.5;

    /* Transform into internal format. */
    j = iauDtf2d ( "UTC", iy, mo, id, ih, im, sec, &utc1, &utc2 );

    /* UT1-UTC (s, from IERS). */
    dut = 0.3340960443019867; // Value modified to coincide with value in code (difference is of order microsecond, only influences utc<->ut1).

    /* UTC -> UT1. */

    j = iauUtcut1 ( utc1, utc2, dut, &ut11, &ut12 );

    /* Extract fraction for TDB-TT calculation, later. */
    ut = fmod ( fmod(ut11,1.0) + fmod(ut12,1.0), 1.0 ) + 0.5;

    /* UTC -> TAI -> TT -> TCG. */
    j = iauUtctai ( utc1, utc2, &tai1, &tai2 );
    j = iauTaitt ( tai1, tai2, &tt1, &tt2 );
    j = iauTttcg ( tt1, tt2, &tcg1, &tcg2 );

    /* TDB-TT (using TT as a substitute for TDB). */
    dtr = iauDtdb ( tt1, tt2, ut, elon, u / 1000.0, v / 1000.0 );

    /* TT -> TDB -> TCB. */
    j = iauTttdb ( tt1, tt2, dtr, &tdb1, &tdb2 );
    j = iauTdbtcb ( tdb1, tdb2, &tcb1, &tcb2 );

    /* Report. */
    sofaTimes.expectedUtcDays = utc1;
    sofaTimes.expectedUtcFractionalDays = utc2;

    sofaTimes.expectedUt1Days = ut11;
    sofaTimes.expectedUt1Seconds = ut12;

    sofaTimes.expectedTaiDays = tai1;
    sofaTimes.expectedTaiFractionalDays = tai2;

    sofaTimes.expectedTtDays = tt1;
    sofaTimes.expectedTtFractionalDays = tt2;

    sofaTimes.expectedTdbDays = tdb1;
    sofaTimes.expectedTdbFractionalDays = tdb2;
    return sofaTimes;
}


double convertSofaOutputToSecondsSinceJ2000(
        const double fullJulianDays, const double fractionalJulianDays )
{
    return ( fullJulianDays - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY + fractionalJulianDays * physical_constants::JULIAN_DAY ;
}

using namespace tudat::earth_orientation;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_time_scale_converters )

// Compare to SOFA Cookbook
BOOST_AUTO_TEST_CASE( testDifferentTimeScaleConversions )
{


    SofaTimeOutput sofaTimes = getSofaDirectTimes( );

    boost::shared_ptr< TerrestrialTimeScaleConverter > timeScaleConverter =
            createStandardEarthOrientationCalculator( )->getTerrestrialTimeScaleConverter( );

    double year = 2006;
    double month = 1;
    double day = 15;
    double hour = 21;

    double stationLatitude, stationLongitude;

    iauAf2a( '+', 19, 28, 52.5, &stationLatitude );
    iauAf2a( '-', 155, 55, 59.6, &stationLongitude);

    double radialDistance = 6.378137E6;

    double testDays = 6.0 * 365 + 15.0;

    double g  = 6.24 + 0.017202 * testDays;

    Eigen::Vector3d stationCartesianPosition;
    stationCartesianPosition << -5492333.306498738, -2453018.508911721, 2113645.653406073;

    std::map< TimeScales, double > sofaSecondsSinceJ2000;
    sofaSecondsSinceJ2000[ tt_scale ] = convertSofaOutputToSecondsSinceJ2000( sofaTimes.expectedTtDays, sofaTimes.expectedTtFractionalDays );
    sofaSecondsSinceJ2000[ utc_scale ] = convertSofaOutputToSecondsSinceJ2000( sofaTimes.expectedUtcDays, sofaTimes.expectedUtcFractionalDays );
    sofaSecondsSinceJ2000[ ut1_scale ] = convertSofaOutputToSecondsSinceJ2000( sofaTimes.expectedUt1Days, sofaTimes.expectedUt1Seconds );
    sofaSecondsSinceJ2000[ tai_scale ] = convertSofaOutputToSecondsSinceJ2000( sofaTimes.expectedTaiDays, sofaTimes.expectedTaiFractionalDays );
    sofaSecondsSinceJ2000[ tdb_scale ] = convertSofaOutputToSecondsSinceJ2000( sofaTimes.expectedTdbDays, sofaTimes.expectedTdbFractionalDays );

    std::vector< TimeScales > originScales;
    originScales.push_back( tt_scale );
    originScales.push_back( utc_scale );
    originScales.push_back( ut1_scale );
    originScales.push_back( tai_scale );
    originScales.push_back( tdb_scale );

    for( unsigned int i = 0; i < originScales.size( ); i++ )
    {
        timeScaleConverter->resetTimes< double >( );
        timeScaleConverter->updateTimes( originScales.at( i ), sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );

        double ut1 = timeScaleConverter->getCurrentTime( originScales.at( i ), ut1_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( ut1 - sofaSecondsSinceJ2000[ ut1_scale ], std::numeric_limits< double >::epsilon( ) );

        double utc = timeScaleConverter->getCurrentTime( originScales.at( i ), utc_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( utc - sofaSecondsSinceJ2000[ utc_scale ], std::numeric_limits< double >::epsilon( ) );

        double tdb = timeScaleConverter->getCurrentTime( originScales.at( i ), tdb_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tdb - sofaSecondsSinceJ2000[ tdb_scale ], std::numeric_limits< double >::epsilon( ) );

        double tai = timeScaleConverter->getCurrentTime( originScales.at( i ), tai_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tai -sofaSecondsSinceJ2000[ tai_scale ], std::numeric_limits< double >::epsilon( ) );

        double tt = timeScaleConverter->getCurrentTime( originScales.at( i ), tt_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tt -sofaSecondsSinceJ2000[ tt_scale ], std::numeric_limits< double >::epsilon( ) );

    }

    for( unsigned int i = 0; i < originScales.size( ); i++ )
    {
        timeScaleConverter->resetTimes< Time >( );
        timeScaleConverter->updateTimes< Time >( originScales.at( i ), sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );

        double ut1 = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), ut1_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( ut1 - sofaSecondsSinceJ2000[ ut1_scale ], std::numeric_limits< double >::epsilon( ) );

        double utc = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), utc_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( utc - sofaSecondsSinceJ2000[ utc_scale ], std::numeric_limits< double >::epsilon( ) );

        double tdb = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), tdb_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tdb - sofaSecondsSinceJ2000[ tdb_scale ], std::numeric_limits< double >::epsilon( ) );

        double tai = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), tai_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tai -sofaSecondsSinceJ2000[ tai_scale ], std::numeric_limits< double >::epsilon( ) );

        double tt = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), tt_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tt -sofaSecondsSinceJ2000[ tt_scale ], std::numeric_limits< double >::epsilon( ) );

    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




