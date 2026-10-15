#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <iostream>
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/earth_orientation/earthOrientationDataEvaluator.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "pr1092/src/astro/earth_orientation/interpolation_reference/expected.h"

using namespace tudat;
using namespace tudat::earth_orientation;
using namespace tudat::basic_astrodynamics;
namespace
{
input_output::EarthOrientationData leapData( )
{
    input_output::EarthOrientationData data;
    data.source = "controlled 2016 leap";
    for( int i = 0; i < 4; ++i )
    {
        input_output::EarthOrientationRecord record;
        record.utcModifiedJulianDate = 57752.0 + i;
        record.polarMotion = { 1e-6 + i * 1e-9, -2e-6 };
        record.ut1MinusUtc = i < 2 ? -0.4 : 0.6;
        record.celestialPoleOffsets = std::array< double, 2 >{ 1e-9, -2e-9 };
        data.records.emplace( record.utcModifiedJulianDate, record );
    }
    return data;
}
Time leapMidnightTai( )
{
    return Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0;
}
}  // namespace

//! Compare the station-dependent TT/TDB offset with SOFA using the selected UT1.
BOOST_AUTO_TEST_CASE( strict_topocentric_tdb_uses_continuous_ut1_in_sofa )
{
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    TerrestrialTimeScaleConverter converter( evaluator, nullptr, nullptr );
    const Eigen::Vector3d station( 4075578.385, 931852.890, 4801570.154 );
    const Time tai = leapMidnightTai( ) - 0.5;
    const Time tt = convertTAItoTT( tai );
    const Time ut1 = tai - 36.4;
    const double radius = std::hypot( station.x( ), station.y( ) ) / 1000.0;
    const double expected = iauDtdb( 2451545.0 + tt.fullDaysSinceEpoch( ),
                                     static_cast< double >( tt.secondsSinceNoon( ) / 86400.0L ),
                                     static_cast< double >( ut1.secondsSinceMidnight( ) / 86400.0L ),
                                     std::atan2( station.y( ), station.x( ) ),
                                     radius,
                                     station.z( ) / 1000.0 );
    const Time tdb = converter.getCurrentTime< Time >( tt_scale, tdb_scale, tt, station );

    // Forming and recovering the offset allows four fractional-hour roundings.
    const double offsetTolerance = std::max(
            1e-14, static_cast< double >( 2.0L * std::numeric_limits< long double >::epsilon( ) * TIME_NORMALIZATION_INTEGER_TERM ) );

    BOOST_CHECK_SMALL( static_cast< double >( tdb - tt ) - expected, offsetTolerance );
    BOOST_CHECK_SMALL( static_cast< double >( converter.getCurrentTime< Time >( tdb_scale, tt_scale, tdb, station ) - tt ), 1e-12 );
}

//! Check the published 1962 C04 values with the historical UTC drift.
BOOST_AUTO_TEST_CASE( pre_1972_utc_drift_uses_fractional_dat_and_continuous_ut1 )
{
    const std::string fixture = tudat::paths::getTudatTestDataPath( ) + "pr1092/src/io/earth_orientation_fixtures/legacy_c04.txt";
    const auto data = input_output::readEarthOrientationData( fixture );

    // The published polar coordinates are arcseconds; the consuming model stores radians.
    const double arcsecond = mathematical_constants::PI / 648000.0;

    BOOST_CHECK_SMALL( data.records.begin( )->second.polarMotion[ 0 ] + 0.012700 * arcsecond, 1e-20 );
    BOOST_CHECK_SMALL( data.records.begin( )->second.polarMotion[ 1 ] - 0.213000 * arcsecond, 1e-20 );
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( data );
    TerrestrialTimeScaleConverter converter( evaluator, nullptr, nullptr );

    // SOFA DAT's accepted 1962 segment: 1.845858 s at January 1 midnight,
    // increasing by 0.0011232 s per UTC day (not an integer leap offset).
    const Time utc = Time( ( 37665.5 - 51544.5 ) * 86400.0 );
    const Time tai = converter.getCurrentTime< Time >( utc_scale, tai_scale, utc );
    BOOST_CHECK_SMALL( static_cast< double >( tai - utc ) - 1.8464196, 1e-12 );
    BOOST_CHECK_SMALL( static_cast< double >( converter.getCurrentTime< Time >( tai_scale, utc_scale, tai ) - utc ), 1e-12 );

    // The daily UT1 offset and its derivative include the non-integer UTC drift.
    const auto value = evaluator->evaluate( tai );
    const double expected = 0.5 * ( ( 0.0326338 - 1.845858 ) + ( 0.0320547 - 1.8469812 ) );
    BOOST_CHECK_SMALL( value.ut1MinusTai - expected, 1e-12 );
    BOOST_CHECK_SMALL( *value.ut1MinusTaiRate - ( ( 0.0320547 - 1.8469812 ) - ( 0.0326338 - 1.845858 ) ) / ( 86400.0 + 0.0011232 ), 1e-18 );
    const Time ut1 = converter.getCurrentTime< Time >( tai_scale, ut1_scale, tai );
    BOOST_CHECK_SMALL( static_cast< double >( ut1 - tai ) - expected, 1e-12 );
    BOOST_CHECK_SMALL( static_cast< double >( converter.getCurrentTime< Time >( ut1_scale, tai_scale, ut1 ) - tai ), 2e-12 );
    const Time fractionalUtc = utc + 12345.6789;
    const Time fractionalTai = converter.getCurrentTime< Time >( utc_scale, tai_scale, fractionalUtc );
    BOOST_CHECK_SMALL( static_cast< double >( fractionalTai - fractionalUtc ) - ( 1.8464196 + 0.0011232 * 12345.6789 / 86400.0 ), 1e-12 );
    BOOST_CHECK_SMALL( static_cast< double >( converter.getCurrentTime< Time >( tai_scale, utc_scale, fractionalTai ) - fractionalUtc ),
                       1e-12 );
}

#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"

namespace
{
input_output::EarthOrientationData interpolationData( const interpolation_reference::Dataset& fixture )
{
    input_output::EarthOrientationData data;
    data.source = "controlled original LAGINT inputs";
    const double arcsecond = mathematical_constants::PI / 648000.0;
    for( const auto& sample : fixture.samples )
    {
        input_output::EarthOrientationRecord record;
        record.utcModifiedJulianDate = sample.utcMjd;
        record.polarMotion = { sample.fields[ 0 ] * arcsecond, sample.fields[ 1 ] * arcsecond };
        record.ut1MinusUtc = sample.fields[ 2 ];
        record.celestialPoleOffsets = std::array< double, 2 >{ sample.fields[ 3 ] * arcsecond, sample.fields[ 4 ] * arcsecond };
        record.lengthOfDay = sample.fields[ 5 ];
        data.records.emplace( sample.utcMjd, record );
    }
    return data;
}
EarthOrientationEvaluationSettings fourPointSettings( )
{
    EarthOrientationEvaluationSettings settings;
    settings.interpolation = EarthOrientationInterpolation::lagrange_four_point;
    return settings;
}
}  // namespace

//! Compare all interpolated EOP quantities and rates with the frozen IERS reference.
BOOST_AUTO_TEST_CASE( four_point_interpolation_matches_unmodified_original_lagint_on_tai_axis )
{
    const double arcsecond = mathematical_constants::PI / 648000.0;
    for( const auto& fixture : interpolation_reference::datasets )
    {
        EarthOrientationDataEvaluator evaluator( interpolationData( fixture ), fourPointSettings( ) );
        const Time firstTai = Time( ( fixture.samples.front( ).utcMjd - 51544.5 ) * 86400.0 ) + fixture.firstDat;
        for( const auto& query : fixture.queries )
        {
            const auto value = evaluator.evaluate( firstTai + query.taiSecondsSinceFirst, true );

            // Compare polar motion, UT1, celestial-pole offsets and LOD with the independent interpolated reference.
            BOOST_CHECK_SMALL( value.polarMotion.x( ) - query.fields[ 0 ] * arcsecond, 2e-20 );
            BOOST_CHECK_SMALL( value.polarMotion.y( ) - query.fields[ 1 ] * arcsecond, 2e-20 );
            BOOST_CHECK_SMALL( value.ut1MinusTai - query.fields[ 2 ], 1e-12 );

            // Compare the six EOP quantities and their TAI derivatives with the independent reference.
            BOOST_REQUIRE( value.celestialPoleOffsets );
            BOOST_CHECK_SMALL( value.celestialPoleOffsets->x( ) - query.fields[ 3 ] * arcsecond, 2e-23 );
            BOOST_CHECK_SMALL( value.celestialPoleOffsets->y( ) - query.fields[ 4 ] * arcsecond, 2e-23 );
            BOOST_REQUIRE( value.lengthOfDay );
            BOOST_CHECK_SMALL( *value.lengthOfDay - query.fields[ 5 ], 1e-15 );
        }

        // The original routine leaves its window index undefined at x==last.
        // This endpoint expectation comes directly from the raw input, not an
        // execution of that undefined original-reference branch.
        const auto last = evaluator.evaluate( evaluator.getEndTai( ), true );
        BOOST_CHECK_EQUAL( last.polarMotion.x( ), fixture.samples.back( ).fields[ 0 ] * arcsecond );
        BOOST_CHECK_EQUAL( last.polarMotion.y( ), fixture.samples.back( ).fields[ 1 ] * arcsecond );
        BOOST_CHECK_EQUAL( *last.lengthOfDay, fixture.samples.back( ).fields[ 5 ] );
    }
}

//! Keep raw node values available when fractional-hour roundoff crosses an interpolation gap.
BOOST_AUTO_TEST_CASE( four_point_nodes_retain_values_across_gaps_with_roundoff )
{
    auto data = interpolationData( interpolation_reference::datasets.front( ) );
    data.records.erase( 59003.0 );
    EarthOrientationDataEvaluator gaps( data, fourPointSettings( ) );
    const Time raw = Time( ( 59001.0 - 51544.5 ) * 86400.0 ) + 37.0;
    const auto exact = gaps.evaluate( raw, true );
    BOOST_CHECK_EQUAL( exact.polarMotion.x( ), data.records.at( 59001.0 ).polarMotion[ 0 ] );

    // Raw values remain available at the node, but missing neighboring data supplies no derivative.
    BOOST_CHECK( !exact.polarMotionRate );
    BOOST_CHECK( !exact.ut1MinusTaiRate );

    // Both neighboring intervals must recognize the same raw node across an interpolation gap.
    const long double roundoff = std::numeric_limits< long double >::epsilon( ) * TIME_NORMALIZATION_INTEGER_TERM;

    for( const long double shift : { -roundoff, roundoff } )
    {
        const auto nearby = gaps.evaluate( raw + shift, true );

        BOOST_CHECK_EQUAL( nearby.polarMotion.x( ), exact.polarMotion.x( ) );
        BOOST_CHECK_EQUAL( nearby.ut1MinusTai, exact.ut1MinusTai );
        BOOST_CHECK( !nearby.polarMotionRate );
        BOOST_CHECK( !nearby.ut1MinusTaiRate );
        BOOST_CHECK_THROW( gaps.evaluate( raw + 32.0L * shift ), std::runtime_error );
    }

    // Genuine missing intervals remain unavailable beyond the node roundoff allowance.
    BOOST_CHECK_THROW( gaps.evaluate( raw + 43200.0 ), std::runtime_error );
    BOOST_CHECK_NO_THROW( gaps.evaluate( gaps.getEndTai( ) - 43200.0, true ) );
}

//! Reject an interpolated clock that reverses direction between increasing raw observations.
BOOST_AUTO_TEST_CASE( four_point_ut1_must_be_monotonic_between_raw_nodes_including_derivative_extremum )
{
    auto data = interpolationData( interpolation_reference::datasets.front( ) );
    data.records.erase( std::next( data.records.begin( ), 4 ), data.records.end( ) );

    // The middle-interval dUT1/dTAI is (u-0.5)^2 - 0.01. Both
    // endpoint slopes are +0.24, but the quadratic extremum is negative.
    // All four sampled UT1 epochs nevertheless increase strictly.
    const std::array< double, 4 > relativeUt1 = { 0.0, 92736.0, 99072.0, 191808.0 };
    std::size_t i = 0;
    for( auto& entry : data.records )
    {
        entry.second.ut1MinusUtc = relativeUt1[ i ] - 86400.0 * i;
        ++i;
    }
    BOOST_CHECK_NO_THROW( static_cast< void >( EarthOrientationDataEvaluator( data ) ) );
    BOOST_CHECK_EXCEPTION(
            EarthOrientationDataEvaluator( data, fourPointSettings( ) ), std::runtime_error, []( const std::runtime_error& error ) {
                return std::string( error.what( ) ).find( "non-positive dUT1/dTAI" ) != std::string::npos;
            } );
}

//! Verify that the assembled rotation uses the selected four-point EOP interpolation.
BOOST_AUTO_TEST_CASE( public_rotation_factory_selects_four_point_reference_policy )
{
    const auto& fixture = interpolation_reference::datasets.front( );
    auto settings = std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >(
            iau_2006, "GCRS", "/must-not-read-this-file", tai_scale );
    settings->setEopData( interpolationData( fixture ) );
    settings->setEopEvaluationSettings( fourPointSettings( ) );
    auto rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( simulation_setup::createRotationModel( settings, "Earth" ) );

    // Check the selected policy and its independently referenced UT1 value.
    BOOST_REQUIRE( rotation );
    const auto evaluator = rotation->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( )->getEarthOrientationDataEvaluator( );
    BOOST_CHECK( evaluator->getSettings( ).interpolation == EarthOrientationInterpolation::lagrange_four_point );
    const auto& query = fixture.queries[ 4 ];
    const Time tai = evaluator->getStartTai( ) + query.taiSecondsSinceFirst;
    BOOST_CHECK_SMALL( evaluator->evaluate( tai ).ut1MinusTai - query.fields[ 2 ], 1e-12 );
}

//! Evaluate an interior clock epoch without requiring distant EOP endpoint coverage.
BOOST_AUTO_TEST_CASE( interior_clock_interpolator_does_not_probe_distant_raw_eop_endpoints )
{
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    const Time tai = evaluator->getStartTai( ) + 100000.0;
    const double tt = static_cast< double >( convertTAItoTT( tai ) );
    auto interpolation = interpolators::interpolatorGenerationSettings< double >(
            interpolators::linearInterpolation( interpolators::huntingAlgorithm, interpolators::throw_exception_at_boundary ),
            tt - 1000.0,
            tt + 1000.0,
            100.0 );
    TerrestrialTimeScaleConverter converter( evaluator, nullptr, interpolation );
    const Time tdb = converter.getCurrentTime< Time >( tai_scale, tdb_scale, tai );

    // An interior epoch must not query clock values outside this intentionally narrow table.
    BOOST_CHECK_SMALL( static_cast< double >( converter.getTaiTime( tdb_scale, static_cast< double >( tdb ) ) - tai ), 1e-7 );
    BOOST_CHECK_SMALL( static_cast< double >( converter.getTaiTime( tdb_scale, tdb ) - tai ), 1e-12 );

    // With interpolated TDB - TT, the rate is that of the interpolated values; it differs from the rate of the SOFA model
    // by the interpolation error only.
    TerrestrialTimeScaleConverter direct( evaluator, nullptr );
    BOOST_CHECK_SMALL( converter.getGeocentricTtRateWithRespectToTdb( tdb ) - direct.getGeocentricTtRateWithRespectToTdb( tdb ), 1e-12 );
}

//! Compare the geocentric clock derivative with independent SOFA offset differences.
BOOST_AUTO_TEST_CASE( native_geocentric_clock_rate_matches_independent_sofa_step_study )
{
    TerrestrialTimeScaleConverter converter( nullptr, nullptr );
    for( const Time& tai : { Time( 3.0e8 ), Time( 536673637.0 ), Time( 1.0e9 ) } )
    {
        const Time tt = convertTAItoTT( tai );
        const Time tdb = converter.getCurrentTime< Time >( tai_scale, tdb_scale, tai );
        const double actual = converter.getGeocentricTtRateWithRespectToTdb( tdb );

        // Propagate four fractional-hour roundings per recovered offset through the 32-second stencil's absolute weights.
        const double rateTolerance = std::max( 2e-15,
                                               static_cast< double >( 2.0L * std::numeric_limits< long double >::epsilon( ) *
                                                                      TIME_NORMALIZATION_INTEGER_TERM * 18.0L / ( 12.0L * 32.0L ) ) );

        for( double h : { 16.0, 32.0, 64.0 } )
        {
            const auto directOffset = [ &tt ]( double delta ) {
                return sofa_interface::getTDBminusTT( static_cast< double >( tt + delta ), 0.0, Eigen::Vector3d::Zero( ) );
            };
            const double derivative =
                    ( 8.0 * ( directOffset( h ) - directOffset( -h ) ) - ( directOffset( 2 * h ) - directOffset( -2 * h ) ) ) / ( 12 * h );
            const double expected = 1.0 / ( 1.0 + derivative );
            BOOST_CHECK_SMALL( actual - expected, rateTolerance );
            std::cout << std::setprecision( 17 ) << "SOFA clock epoch=" << static_cast< double >( tai ) << " h=" << h
                      << " dTT/dTDB=" << actual << " residual=" << actual - expected << '\n';
        }

        // At zero observer all fraction-of-day arguments give identical offsets.
        BOOST_CHECK_EQUAL( sofa_interface::getTDBminusTT( static_cast< double >( tt ), 0.0, Eigen::Vector3d::Zero( ) ),
                           sofa_interface::getTDBminusTT( static_cast< double >( tt ), 0.999, Eigen::Vector3d::Zero( ) ) );
    }
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    TerrestrialTimeScaleConverter withEop( evaluator, nullptr, nullptr );
    const Time tai = evaluator->getStartTai( ) + 100000.0;
    const Time tdb = withEop.getCurrentTime< Time >( tai_scale, tdb_scale, tai );
    const double before = withEop.getGeocentricTtRateWithRespectToTdb( tdb );
    BOOST_CHECK_NO_THROW( withEop.getCurrentTime< Time >( tai_scale, tdb_scale, tai, Eigen::Vector3d( 6378137.0, 0.0, 0.0 ) ) );
    BOOST_CHECK_EQUAL( withEop.getGeocentricTtRateWithRespectToTdb( tdb ), before );
}

//! Preserve the clock derivative when a Time stencil crosses a fractional-hour precision boundary.
BOOST_AUTO_TEST_CASE( precise_time_clock_stencil_accepts_roundoff_at_a_fractional_binade_boundary )
{
    const Time tdb( 0, std::nextafter( 100.0L, std::numeric_limits< long double >::infinity( ) ) );

    // Adding 32/64 s crosses the 128-s fractional binade and rounds by one
    // long-double ULP; its size depends on the platform's long-double precision.
    BOOST_CHECK_NE( tdb + 32.0 - tdb, Time( 32.0 ) );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr );
    const double actual = converter.getGeocentricTtRateWithRespectToTdb( tdb );
    const Time tt = convertTAItoTT( converter.getTaiTime( tdb_scale, tdb ) );
    const auto offset = [ &tt ]( double delta ) {
        return sofa_interface::getTDBminusTT( static_cast< double >( tt + delta ), 0.0, Eigen::Vector3d::Zero( ) );
    };
    const double derivative = ( 8.0 * ( offset( 32 ) - offset( -32 ) ) - ( offset( 64 ) - offset( -64 ) ) ) / 384.0;

    // Use the same fractional-hour offset budget and absolute stencil weights as the clock step study above.
    const double rateTolerance = std::max( 2e-15,
                                           static_cast< double >( 2.0L * std::numeric_limits< long double >::epsilon( ) *
                                                                  TIME_NORMALIZATION_INTEGER_TERM * 18.0L / ( 12.0L * 32.0L ) ) );

    BOOST_CHECK_SMALL( actual - 1.0 / ( 1.0 + derivative ), rateTolerance );
}
