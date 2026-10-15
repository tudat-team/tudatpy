/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <algorithm>
#include <limits>
#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/math/interpolators/jumpDataLinearInterpolator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::earth_orientation;
using namespace tudat::interpolators;

BOOST_AUTO_TEST_SUITE( test_eop_reader )

//! Tests whether EOP data is properly read and processed, and whether data in correctly interpolated (linearly)
BOOST_AUTO_TEST_CASE( testEopReaderData )
{
    std::shared_ptr< EOPReader > eopReader = std::make_shared< EOPReader >(
            tudat::paths::getEarthOrientationDataFilesPath( ) + "/eopc04_14_IAU2000.62-now.txt", "C04", basic_astrodynamics::iau_2000_a );

    // Define current time/
    int year = 2007;
    int month = 4;
    int day = 5;
    int hour = 0;
    int minute = 0;
    double seconds = 0.0;
    double utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                                          year, month, day, hour, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Set EOP corrections read manually from file for 5-4-2007
    double arcSecondToRadian = 4.848136811095359935899141E-6;
    double expectedXp = 0.033194 * arcSecondToRadian;
    double expectedYp = 0.483144 * arcSecondToRadian;
    double expectedUT1COffset = -0.0714163;
    double expecteddX = 0.000250 * arcSecondToRadian;
    double expecteddY = -0.000302 * arcSecondToRadian;

    // Create interpolator for x_{p} and y_{p} (polar motion) and retrieve current value
    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector2d > > cipInItrsInterpolator =
            std::make_shared< LinearInterpolator< double, Eigen::Vector2d > >( eopReader->getCipInItrsMapInSecondsSinceJ2000( ) );
    Eigen::Vector2d cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Create interpolator for dX and dY (nutation correction) and retrieve current value
    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector2d > > cipInGcrsCorrectionInterpolator =
            std::make_shared< LinearInterpolator< double, Eigen::Vector2d > >( eopReader->getCipInGcrsCorrectionMapInSecondsSinceJ2000( ) );
    Eigen::Vector2d cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Create interpolator for d(UT1-UTC) correction  and retrieve current value
    std::shared_ptr< OneDimensionalInterpolator< double, double > > ut1MinusUtcInterpolator =
            std::make_shared< JumpDataLinearInterpolator< double, double > >(
                    eopReader->getUt1MinusUtcMapInSecondsSinceJ2000( ), 0.5, 1.0 );
    double utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Check interpolated values against test data.
    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1, expectedUT1COffset, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), expectedXp, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), expectedYp, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), expecteddX, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), expecteddY, std::numeric_limits< double >::epsilon( ) );

    // Reset time
    int day2 = 6;
    utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                                   year, month, day2, hour, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Set EOP corrections read manually from file for 6-4-2007
    double expectedXp2 = 0.035736 * arcSecondToRadian;
    double expectedYp2 = 0.484204 * arcSecondToRadian;
    double expectedUT1COffset2 = -0.0728009;
    double expecteddX2 = 0.000254 * arcSecondToRadian;
    double expecteddY2 = -0.000309 * arcSecondToRadian;

    // Read correction interpolators
    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );
    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );
    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Check interpolated values against test data.
    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1, expectedUT1COffset2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), expectedXp2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), expectedYp2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), expecteddX2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), expecteddY2, std::numeric_limits< double >::epsilon( ) );

    // Reset time to halfway between previous times
    int hour2 = 12;
    utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                                   year, month, day, hour2, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Read correction interpolators
    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );
    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );
    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    BOOST_CHECK_CLOSE_FRACTION(
            utcMinusUt1, ( expectedUT1COffset + expectedUT1COffset2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), ( expectedXp + expectedXp2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), ( expectedYp + expectedYp2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), ( expecteddX + expecteddX2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), ( expecteddY + expecteddY2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );

    // Reset time to halfway between previous times
    int hour3 = 1;
    utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                                   year, month, day, hour3, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Read correction interpolators
    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );
    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );
    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1,
                                expectedUT1COffset + ( expectedUT1COffset2 - expectedUT1COffset ) / 24.0,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            cipInItrs.x( ), expectedXp + ( expectedXp2 - expectedXp ) / 24.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            cipInItrs.y( ), expectedYp + ( expectedYp2 - expectedYp ) / 24.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            cipInGcrs.x( ), expecteddX + ( expecteddX2 - expecteddX ) / 24.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            cipInGcrs.y( ), expecteddY + ( expecteddY2 - expecteddY ) / 24.0, std::numeric_limits< double >::epsilon( ) );
}

//! Check whether UT1-UTC interpolator correctly handles leap secons
BOOST_AUTO_TEST_CASE( testLeapSecondIdentification )
{
    // Read EOP observations and construct a continuous TAI/UT1 clock.
    std::shared_ptr< EOPReader > eopReader = std::make_shared< EOPReader >(
            tudat::paths::getEarthOrientationDataFilesPath( ) + "/eopc04_14_IAU2000.62-now.txt", "C04", basic_astrodynamics::iau_2000_a );
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( eopReader->getRawEarthOrientationData( ) );
    TerrestrialTimeScaleConverter converter( evaluator, nullptr, nullptr );
    const auto rawUt1 = eopReader->getUt1MinusUtcMapInSecondsSinceJ2000( );

    // Define list of leap seconds
    Eigen::Matrix< int, Eigen::Dynamic, 3 > leapSecondDays;
    leapSecondDays.resize( 27, 3 );
    leapSecondDays << 1, 7, 1972, 1, 1, 1973, 1, 1, 1974, 1, 1, 1975, 1, 1, 1976, 1, 1, 1977, 1, 1, 1978, 1, 1, 1979, 1, 1, 1980, 1, 7,
            1981, 1, 7, 1982, 1, 7, 1983, 1, 7, 1985, 1, 1, 1988, 1, 1, 1990, 1, 1, 1991, 1, 7, 1992, 1, 7, 1993, 1, 7, 1994, 1, 1, 1996, 1,
            7, 1997, 1, 1, 1999, 1, 1, 2006, 1, 1, 2009, 1, 7, 2012, 1, 7, 2015, 1, 1, 2017;

    // Declare test variables.
    double differenceMinusOneDay, differenceMinusOneMilliSecond, differenceAtEpochPlusOneMilliSecond, differencePlusOneDay;

    // Set check time: at midnight
    int hour = 0;
    int minute = 0;
    double seconds = 0.0;

    // Iterate over all leap seconds
    for( int i = 0; i < leapSecondDays.rows( ); i++ )
    {
        // Get UTC time of current leap second
        int year, month, day;
        year = leapSecondDays( i, 2 );
        month = leapSecondDays( i, 1 );
        day = leapSecondDays( i, 0 );
        double utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                                              year, month, day, hour, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                physical_constants::JULIAN_DAY;

        // Define test times
        differenceMinusOneDay = rawUt1.at( utcSecondsSinceJ2000 - physical_constants::JULIAN_DAY );
        const Time beforeTai = sofa_interface::convertUTCtoTAI( Time( utcSecondsSinceJ2000 - 1.0E-3 ) );
        const Time afterTai = sofa_interface::convertUTCtoTAI( Time( utcSecondsSinceJ2000 + 1.0E-3 ) );
        const auto before = evaluator->evaluate( beforeTai );
        const auto after = evaluator->evaluate( afterTai );
        differenceMinusOneMilliSecond = before.ut1MinusTai + static_cast< double >( beforeTai - Time( utcSecondsSinceJ2000 - 1.0E-3 ) );
        differenceAtEpochPlusOneMilliSecond = after.ut1MinusTai + static_cast< double >( afterTai - Time( utcSecondsSinceJ2000 + 1.0E-3 ) );
        differencePlusOneDay = rawUt1.at( utcSecondsSinceJ2000 + physical_constants::JULIAN_DAY );

        // Check if leap second is properly handled
        BOOST_CHECK_SMALL( std::fabs( differenceMinusOneDay - differenceMinusOneMilliSecond ), 0.005 );
        BOOST_CHECK_SMALL( std::fabs( differencePlusOneDay - differenceAtEpochPlusOneMilliSecond ), 0.005 );

        // Scalar UTC jumps over the inserted second. The two samples are
        // 1.002 TAI seconds apart, so continuous UT1-TAI also evolves between
        // them. Account for its independently specified daily linear slopes.
        const double atMidnight = rawUt1.at( utcSecondsSinceJ2000 );
        const double beforeDuration = 1.0 + ( utcSecondsSinceJ2000 - ( utcSecondsSinceJ2000 - 1.0E-3 ) );
        const double afterDuration = ( utcSecondsSinceJ2000 + 1.0E-3 ) - utcSecondsSinceJ2000;
        const double expectedEvolution = -beforeDuration * ( atMidnight - differenceMinusOneDay - 1.0 ) / 86401.0 -
                afterDuration * ( differencePlusOneDay - atMidnight ) / 86400.0;
        const Time beforeUt1 = converter.getUt1FromTai( beforeTai );
        const Time afterUt1 = converter.getUt1FromTai( afterTai );
        const double evolution = static_cast< double >( ( beforeUt1 - afterUt1 ) - ( beforeTai - afterTai ) );

        // Four fractional-hour roundings can enter the two UT1 epochs and their differenced evolution.
        const long double evolutionTolerance = std::max(
                1.0E-14L, 2.0L * std::numeric_limits< long double >::epsilon( ) * TIME_NORMALIZATION_INTEGER_TERM );

        BOOST_CHECK_SMALL( evolution - expectedEvolution, evolutionTolerance );
        BOOST_CHECK_SMALL( std::fabs( differenceMinusOneDay - differencePlusOneDay + 1.0 ), 0.01 );
    }
}

BOOST_AUTO_TEST_CASE( test_preloaded_Eop_data_adapter_equivalence_and_provenance )
{
    const std::string sourceFile = __FILE__;
    const std::string file =
            sourceFile.substr( 0, sourceFile.find_last_of( "/\\" ) + 1 ) + "../../io/earth_orientation_fixtures/modern_c04.txt";
    const auto data = input_output::readEarthOrientationData( file );
    EOPReader fromFile( file );
    EOPReader fromData( data );
    BOOST_CHECK( fromFile.getUt1MinusUtcMapRaw( ) == fromData.getUt1MinusUtcMapRaw( ) );
    BOOST_CHECK( fromFile.getLengthOfDayMapRaw( ) == fromData.getLengthOfDayMapRaw( ) );
    const auto filePolarMotion = fromFile.getCipInItrsMapRaw( );
    const auto dataPolarMotion = fromData.getCipInItrsMapRaw( );
    BOOST_REQUIRE_EQUAL( filePolarMotion.size( ), dataPolarMotion.size( ) );
    for( const auto& entry : filePolarMotion )
    {
        BOOST_CHECK_SMALL( ( entry.second - dataPolarMotion.at( entry.first ) ).norm( ), 1e-30 );
    }
    const auto fileCpo = fromFile.getCipInGcrsCorrectionMapRaw( );
    const auto dataCpo = fromData.getCipInGcrsCorrectionMapRaw( );
    BOOST_REQUIRE_EQUAL( fileCpo.size( ), dataCpo.size( ) );
    for( const auto& entry : fileCpo )
    {
        BOOST_CHECK_SMALL( ( entry.second - dataCpo.at( entry.first ) ).norm( ), 1e-30 );
    }
    const auto& raw = fromData.getRawEarthOrientationData( );
    BOOST_CHECK_EQUAL( raw.source, data.source );
    BOOST_CHECK( raw.header == data.header );
    BOOST_CHECK_EQUAL( raw.terrestrialFrame, data.terrestrialFrame );
    BOOST_REQUIRE_EQUAL( raw.records.size( ), 3 );
    BOOST_CHECK_EQUAL( raw.records.begin( )->first, 45700.5 );
    BOOST_CHECK_EQUAL( raw.records.rbegin( )->first, 45702.5 );

    // Compatibility maps retain the same exact raw coverage.
    BOOST_CHECK_EQUAL( fromData.getUt1MinusUtcMapRaw( ).size( ), 3 );
    BOOST_CHECK_SMALL( fromData.getUt1MinusUtcMapRaw( ).begin( )->first - 45700.5, 1e-10 );
}

BOOST_AUTO_TEST_CASE( test_preloaded_Eop_data_keeps_optional_absence )
{
    const std::string sourceFile = __FILE__;
    const std::string file =
            sourceFile.substr( 0, sourceFile.find_last_of( "/\\" ) + 1 ) + "../../io/earth_orientation_fixtures/finals2000a.txt";
    auto data = input_output::readEarthOrientationData( file, input_output::EarthOrientationFormat::finals2000a );
    EOPReader reader( data );
    BOOST_CHECK_EQUAL( reader.getUt1MinusUtcMapRaw( ).count( 61396.0 ), 1 );
    BOOST_CHECK_EQUAL( reader.getLengthOfDayMapRaw( ).count( 61396.0 ), 0 );
    BOOST_CHECK_EQUAL( reader.getCipInGcrsCorrectionMapRaw( ).count( 61396.0 ), 0 );
    BOOST_CHECK( !reader.getRawEarthOrientationData( ).records.at( 61396.0 ).celestialPoleOffsets );
    BOOST_CHECK_EQUAL( reader.getRawEarthOrientationData( ).records.size( ), 4 );
    for( auto& entry : data.records )
    {
        entry.second.lengthOfDay.reset( );
        entry.second.celestialPoleOffsets.reset( );
        entry.second.lengthOfDayUncertainty.reset( );
        entry.second.celestialPoleOffsetUncertainty.reset( );
        entry.second.celestialPoleStatus.clear( );
    }
    EOPReader withoutOptionalData( data );
    BOOST_CHECK( withoutOptionalData.getLengthOfDayMapRaw( ).empty( ) );
    BOOST_CHECK( withoutOptionalData.getCipInGcrsCorrectionMapRaw( ).empty( ) );
    BOOST_CHECK_EQUAL( withoutOptionalData.getUt1MinusUtcMapRaw( ).size( ), data.records.size( ) );
    input_output::EarthOrientationData empty;
    BOOST_CHECK_THROW( static_cast< void >( EOPReader( empty ) ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( test_preloaded_Eop_data_rejects_invalid_records )
{
    input_output::EarthOrientationData data;
    input_output::EarthOrientationRecord record;
    data.records.emplace( 51544.0, record );
    BOOST_CHECK_THROW( static_cast< void >( EOPReader( data ) ), std::runtime_error );
    data.records.begin( )->second.utcModifiedJulianDate = 51544.0;
    data.records.begin( )->second.polarMotion = { 0.0, 0.0 };
    data.records.begin( )->second.ut1MinusUtc = 0.0;
    BOOST_CHECK_NO_THROW( static_cast< void >( EOPReader( data ) ) );
    data.records.begin( )->second.utcModifiedJulianDate = 51545.0;
    BOOST_CHECK_THROW( static_cast< void >( EOPReader( data ) ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
