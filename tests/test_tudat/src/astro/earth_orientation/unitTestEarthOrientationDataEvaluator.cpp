#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <limits>
#include <iomanip>
#include <iostream>
#include "tudat/astro/earth_orientation/earthOrientationDataEvaluator.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "interpolation_reference/expected.h"

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

BOOST_AUTO_TEST_CASE( continuous_tai_axis_values_slopes_and_exact_coverage )
{
    const auto data = leapData( );
    EarthOrientationDataEvaluator model( data );
    BOOST_CHECK_EQUAL( model.getRawData( ).records.size( ), 4 );
    BOOST_CHECK_EQUAL( model.getRawData( ).records.begin( )->first, 57752.0 );
    BOOST_CHECK_NO_THROW( model.evaluate( model.getStartTai( ), true ) );
    BOOST_CHECK_NO_THROW( model.evaluate( model.getEndTai( ), true ) );
    BOOST_CHECK_THROW( model.evaluate( model.getStartTai( ) - 1e-3 ), std::runtime_error );
    BOOST_CHECK_THROW( model.evaluate( model.getEndTai( ) + 1e-3 ), std::runtime_error );
    const auto before = model.evaluate( leapMidnightTai( ) - 0.5 );
    const auto after = model.evaluate( leapMidnightTai( ) + 0.5 );
    BOOST_CHECK_SMALL( before.ut1MinusTai + 36.4, 1e-13 );
    BOOST_CHECK_SMALL( after.ut1MinusTai + 36.4, 1e-13 );
    BOOST_CHECK_SMALL( *before.ut1MinusTaiRate, 1e-18 );
    BOOST_CHECK_SMALL( before.polarMotionRate->x( ) - 1e-9 / 86401.0, 1e-24 );
}

BOOST_AUTO_TEST_CASE( leap_second_ut1_is_continuous_and_inverse_stays_in_coverage )
{
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
    const Time midnight = leapMidnightTai( );
    for( double offset : { -2.0, -1.0, -0.5, 0.0, 0.5, 1.0 } )
    {
        const Time tai = midnight + offset;
        const Time ut1 = converter.getCurrentTime< Time >( tai_scale, ut1_scale, tai );
        BOOST_CHECK_SMALL( static_cast< double >( ut1 - tai ) + 36.4, 1e-12 );
        const Time roundTrip = converter.getCurrentTime< Time >( ut1_scale, tai_scale, ut1 );
        BOOST_CHECK_SMALL( static_cast< double >( roundTrip - tai ), 2e-12 );
        const Time fromTt = converter.getCurrentTime< Time >( tt_scale, ut1_scale, convertTAItoTT( tai ) );
        BOOST_CHECK_SMALL( static_cast< double >( fromTt - ut1 ), 1e-12 );
    }
    for( const Time& tai : { evaluator->getStartTai( ), evaluator->getEndTai( ) } )
    {
        const Time ut1 = converter.getCurrentTime< Time >( tai_scale, ut1_scale, tai );
        BOOST_CHECK_EQUAL( converter.getCurrentTime< Time >( ut1_scale, tai_scale, ut1 ), tai );
        const double roundedUt1 = converter.getCurrentTime< double >( tai_scale, ut1_scale, static_cast< double >( tai ) );
        BOOST_CHECK_EQUAL( converter.getCurrentTime< double >( ut1_scale, tai_scale, roundedUt1 ), static_cast< double >( tai ) );
    }
    const Time startUt1 = converter.getCurrentTime< Time >( tai_scale, ut1_scale, evaluator->getStartTai( ) );
    BOOST_CHECK_THROW( converter.getCurrentTime< Time >( ut1_scale, tai_scale, startUt1 - 1e-3 ), std::runtime_error );
    BOOST_CHECK_THROW( converter.getCurrentTime< Time >( tai_scale, ut1_scale, evaluator->getEndTai( ) + 1e-3 ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( atomic_conversions_do_not_require_eop_or_missing_cpo )
{
    auto data = leapData( );
    for( auto& entry : data.records )
    {
        entry.second.celestialPoleOffsets.reset( );
    }
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( data );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
    const Time outside = evaluator->getEndTai( ) + 1000.0 * 86400.0;
    const Time tt = converter.getCurrentTime< Time >( tai_scale, tt_scale, outside );
    BOOST_CHECK_SMALL( static_cast< double >( tt - outside ) - 32.184, 1e-12 );
    const Time tdb = converter.getCurrentTime< Time >( tt_scale, tdb_scale, tt );
    const Time back = converter.getCurrentTime< Time >( tdb_scale, tt_scale, tdb );
    BOOST_CHECK_SMALL( static_cast< double >( back - tt ), 1e-12 );
    BOOST_CHECK_NO_THROW( converter.getCurrentTime< Time >( tai_scale, ut1_scale, evaluator->getStartTai( ) ) );
    BOOST_CHECK_THROW( evaluator->evaluate( evaluator->getStartTai( ), true ), std::runtime_error );
    EarthOrientationEvaluationSettings allowMissing;
    allowMissing.missingCelestialPole = MissingCelestialPolePolicy::zero;
    EarthOrientationDataEvaluator explicitZero( data, allowMissing );
    BOOST_CHECK_EQUAL( explicitZero.getCelestialPoleOffsets( explicitZero.getStartTai( ) ).norm( ), 0.0 );
}

BOOST_AUTO_TEST_CASE( strict_topocentric_tdb_uses_continuous_ut1_in_sofa )
{
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
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
    BOOST_CHECK_SMALL( static_cast< double >( tdb - tt ) - expected, 1e-14 );
    BOOST_CHECK_SMALL( static_cast< double >( converter.getCurrentTime< Time >( tdb_scale, tt_scale, tdb, station ) - tt ), 1e-12 );
    BOOST_CHECK_THROW(
            converter.getCurrentTime< Time >( tt_scale, tdb_scale, convertTAItoTT( evaluator->getEndTai( ) + 86400.0 ), station ),
            std::runtime_error );
    BOOST_CHECK_NO_THROW( converter.getCurrentTime< Time >( tt_scale, tdb_scale, convertTAItoTT( evaluator->getEndTai( ) + 86400.0 ) ) );
}

BOOST_AUTO_TEST_CASE( strict_topocentric_tdb_endpoint_images_round_trip_without_extrapolation )
{
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
    const Eigen::Vector3d station( 4075578.385, 931852.890, 4801570.154 );
    for( const Time& endpoint : { evaluator->getStartTai( ), evaluator->getEndTai( ) } )
    {
        const Time tt = convertTAItoTT( endpoint );
        const Time tdb = converter.getCurrentTime< Time >( tt_scale, tdb_scale, tt, station );
        BOOST_CHECK_EQUAL( converter.getCurrentTime< Time >( tdb_scale, tt_scale, tdb, station ), tt );
        const double roundedTt = static_cast< double >( tt );
        const double roundedTdb = converter.getCurrentTime< double >( tt_scale, tdb_scale, roundedTt, station );
        BOOST_CHECK_EQUAL( converter.getCurrentTime< double >( tdb_scale, tt_scale, roundedTdb, station ), roundedTt );
    }
    for( const Time& tai : { evaluator->getStartTai( ) + 1e-9, evaluator->getEndTai( ) - 1e-9 } )
    {
        const Time tdb = converter.getCurrentTime< Time >( tai_scale, tdb_scale, tai, station );
        BOOST_CHECK_SMALL( static_cast< double >( converter.getCurrentTime< Time >( tdb_scale, tai_scale, tdb, station ) - tai ), 2e-12 );
    }
    for( const Time& tai : { evaluator->getStartTai( ) + 1.0, evaluator->getEndTai( ) - 1.0 } )
    {
        const double epoch = static_cast< double >( tai );
        const double tdb = converter.getCurrentTime< double >( tai_scale, tdb_scale, epoch, station );
        BOOST_CHECK_SMALL( converter.getCurrentTime< double >( tdb_scale, tai_scale, tdb, station ) - epoch, 2e-7 );
    }
    const Time firstTdb = converter.getCurrentTime< Time >( tai_scale, tdb_scale, evaluator->getStartTai( ), station );
    const Time lastTdb = converter.getCurrentTime< Time >( tai_scale, tdb_scale, evaluator->getEndTai( ), station );
    BOOST_CHECK_THROW( converter.getCurrentTime< Time >( tdb_scale, tai_scale, firstTdb - 1e-3, station ), std::runtime_error );
    BOOST_CHECK_THROW( converter.getCurrentTime< Time >( tdb_scale, tai_scale, lastTdb + 1e-3, station ), std::runtime_error );
    BOOST_CHECK_THROW( converter.getCurrentTime< double >( tdb_scale, tai_scale, static_cast< double >( firstTdb ) - 1.0, station ),
                       std::runtime_error );
    BOOST_CHECK_THROW( converter.getCurrentTime< double >( tdb_scale, tai_scale, static_cast< double >( lastTdb ) + 1.0, station ),
                       std::runtime_error );
}

BOOST_AUTO_TEST_CASE( pre_1972_utc_drift_uses_fractional_dat_and_continuous_ut1 )
{
    const std::string sourceFile = __FILE__;
    const std::string fixture =
            sourceFile.substr( 0, sourceFile.find_last_of( "/\\" ) + 1 ) + "../../io/earth_orientation_fixtures/legacy_c04.txt";
    const auto data = input_output::readEarthOrientationData( fixture );
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( data );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );

    // SOFA DAT's accepted 1962 segment: 1.845858 s at January 1 midnight,
    // increasing by 0.0011232 s per UTC day (not an integer leap offset).
    const Time utc = Time( ( 37665.5 - 51544.5 ) * 86400.0 );
    const Time tai = converter.getCurrentTime< Time >( utc_scale, tai_scale, utc );
    BOOST_CHECK_SMALL( static_cast< double >( tai - utc ) - 1.8464196, 1e-12 );
    BOOST_CHECK_SMALL( static_cast< double >( converter.getCurrentTime< Time >( tai_scale, utc_scale, tai ) - utc ), 1e-12 );
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

BOOST_AUTO_TEST_CASE( strict_scalar_utc_rejects_inserted_second )
{
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
    BOOST_CHECK_THROW( converter.getCurrentTime< Time >( tai_scale, utc_scale, leapMidnightTai( ) - 0.5 ), std::runtime_error );
    BOOST_CHECK_NO_THROW( converter.getCurrentTime< Time >( tai_scale, utc_scale, leapMidnightTai( ) - 2.0 ) );
    BOOST_CHECK_NO_THROW( converter.getCurrentTime< Time >( tai_scale, utc_scale, leapMidnightTai( ) ) );
    BOOST_CHECK_NO_THROW( converter.getCurrentTime< Time >( tai_scale, utc_scale, leapMidnightTai( ) - 1.0 - 1e-12 ) );
    BOOST_CHECK_THROW( converter.getCurrentTime< Time >( tai_scale, utc_scale, leapMidnightTai( ) - 1.0 + 1e-12 ), std::runtime_error );
    BOOST_CHECK_THROW( converter.getCurrentTime< Time >( tai_scale, utc_scale, leapMidnightTai( ) - 1e-12 ), std::runtime_error );
    BOOST_CHECK_NO_THROW( converter.getCurrentTime< Time >( tai_scale, utc_scale, leapMidnightTai( ) + 1e-12 ) );
}

BOOST_AUTO_TEST_CASE( strict_gap_and_optional_cpo_endpoint_policies )
{
    auto data = leapData( );
    data.records.erase( 57753.0 );
    EarthOrientationDataEvaluator model( data );
    BOOST_CHECK_NO_THROW( model.evaluate( model.getStartTai( ) ) );
    BOOST_CHECK( !model.evaluate( model.getStartTai( ) ).polarMotionRate );
    BOOST_CHECK( !model.evaluate( model.getStartTai( ) ).ut1MinusTaiRate );
    BOOST_CHECK( !model.evaluate( model.getStartTai( ) ).celestialPoleOffsetRate );
    BOOST_CHECK_THROW( model.evaluate( model.getStartTai( ) + 86400.0 ), std::runtime_error );
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( data );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
    const Time valid = evaluator->getEndTai( ) - 100.0;
    const Time ut1 = converter.getCurrentTime< Time >( tai_scale, ut1_scale, valid );
    BOOST_CHECK_SMALL( static_cast< double >( converter.getCurrentTime< Time >( ut1_scale, tai_scale, ut1 ) - valid ), 2e-12 );
}

#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"

BOOST_AUTO_TEST_CASE( public_rotation_factory_is_continuous_during_the_inserted_second )
{
    auto settings = std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >(
            iau_2006, "GCRS", "/must-not-read-this-file", tai_scale );
    auto data = leapData( );
    data.records.rbegin( )->second.polarMotionStatus = "P";
    settings->setEopData( data );
    const auto rotation = simulation_setup::createRotationModel( settings, "Earth" );
    const auto earthRotation = std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( rotation );
    BOOST_REQUIRE( earthRotation );
    const auto angles = earthRotation->getAnglesCalculator( );
    const auto converter = angles->getTerrestrialTimeScaleConverter( );
    BOOST_CHECK_EQUAL( converter->getEarthOrientationDataEvaluator( )->getRawData( ).source, data.source );
    BOOST_CHECK_EQUAL( converter->getEarthOrientationDataEvaluator( )->getRawData( ).records.rbegin( )->second.polarMotionStatus, "P" );
    const Time start = leapMidnightTai( ) - 1.5;
    const Time middle = leapMidnightTai( ) - 0.5;
    const Time end = leapMidnightTai( ) + 0.5;
    const auto a = angles->getRotationAnglesFromItrsToGcrs< Time >( start, tai_scale );
    const auto b = angles->getRotationAnglesFromItrsToGcrs< Time >( middle, tai_scale );
    const auto c = angles->getRotationAnglesFromItrsToGcrs< Time >( end, tai_scale );
    BOOST_CHECK_SMALL( static_cast< double >( b.second - a.second ) - 1.0, 2e-8 );
    BOOST_CHECK_SMALL( static_cast< double >( c.second - b.second ) - 1.0, 2e-8 );
    const Eigen::Quaterniond ra = earthRotation->getRotationToTargetFrameFromExtendedTime( start );
    const Eigen::Quaterniond rb = earthRotation->getRotationToTargetFrameFromExtendedTime( middle );
    const Eigen::Quaterniond rc = earthRotation->getRotationToTargetFrameFromExtendedTime( end );
    const double expectedSpin = 2.0 * mathematical_constants::PI * 1.00273781191135448 / 86400.0;
    BOOST_CHECK_SMALL( Eigen::AngleAxisd( rb * ra.inverse( ) ).angle( ) - expectedSpin, 1e-9 );
    BOOST_CHECK_SMALL( Eigen::AngleAxisd( rc * rb.inverse( ) ).angle( ) - expectedSpin, 1e-9 );
    BOOST_CHECK_THROW( converter->getCurrentTime< Time >( tai_scale, utc_scale, middle ), std::runtime_error );

    // The TT-only PM accessor is also used by solid/ocean pole tide.
    const Eigen::Vector2d pm =
            angles->getPolarMotionCalculator( )->getPositionOfCipInItrsAtTt( static_cast< double >( convertTAItoTT( middle ) ) );
    BOOST_CHECK_SMALL( ( pm - b.first.segment< 2 >( 3 ) ).norm( ), 1e-15 );
    const auto evaluator = converter->getEarthOrientationDataEvaluator( );
    BOOST_CHECK_NO_THROW( earthRotation->getRotationToTargetFrameFromExtendedTime( evaluator->getStartTai( ) ) );
    BOOST_CHECK_NO_THROW( earthRotation->getRotationToTargetFrameFromExtendedTime( evaluator->getEndTai( ) ) );

    // Both inverse maps include the actual production short-period UT1 terms.
    const Eigen::Vector3d station( 4075578.385, 931852.890, 4801570.154 );
    for( const Time& endpoint : { evaluator->getStartTai( ), evaluator->getEndTai( ) } )
    {
        const Time ut1 = converter->getCurrentTime< Time >( tai_scale, ut1_scale, endpoint );
        BOOST_CHECK_EQUAL( converter->getCurrentTime< Time >( ut1_scale, tai_scale, ut1 ), endpoint );
        const Time tdb = converter->getCurrentTime< Time >( tai_scale, tdb_scale, endpoint, station );
        BOOST_CHECK_EQUAL( converter->getCurrentTime< Time >( tdb_scale, tai_scale, tdb, station ), endpoint );
        BOOST_CHECK_EQUAL( converter->getCurrentTime< double >( ut1_scale, tai_scale, static_cast< double >( ut1 ) ),
                           static_cast< double >( endpoint ) );
        BOOST_CHECK_EQUAL( converter->getCurrentTime< double >( tdb_scale, tai_scale, static_cast< double >( tdb ), station ),
                           static_cast< double >( endpoint ) );
    }
    BOOST_CHECK_THROW( earthRotation->getRotationToTargetFrameFromExtendedTime( evaluator->getStartTai( ) - 1e-3 ), std::runtime_error );
    BOOST_CHECK_THROW( earthRotation->getRotationToTargetFrameFromExtendedTime( evaluator->getEndTai( ) + 1e-3 ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( public_factory_obeys_missing_cpo_and_gap_policies_only_when_needed )
{
    auto data = leapData( );
    for( auto& entry : data.records )
    {
        entry.second.celestialPoleOffsets.reset( );
    }
    auto settings = std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >(
            iau_2006, "GCRS", "/must-not-read-this-file", tai_scale );
    settings->setEopData( data );
    auto rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( simulation_setup::createRotationModel( settings, "Earth" ) );
    BOOST_REQUIRE( rotation );
    auto converter = rotation->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( );
    const Time epoch = leapMidnightTai( ) - 0.5;
    BOOST_CHECK_NO_THROW( converter->getCurrentTime< Time >( tai_scale, ut1_scale, epoch ) );
    BOOST_CHECK_NO_THROW( converter->getCurrentTime< Time >( tai_scale, tdb_scale, epoch + 1000.0 * 86400.0 ) );
    BOOST_CHECK_THROW( rotation->getRotationToTargetFrameFromExtendedTime( epoch ), std::runtime_error );
    EarthOrientationEvaluationSettings zeroCpo;
    zeroCpo.missingCelestialPole = MissingCelestialPolePolicy::zero;
    settings->setEopEvaluationSettings( zeroCpo );
    rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( simulation_setup::createRotationModel( settings, "Earth" ) );
    BOOST_CHECK_NO_THROW( rotation->getRotationToTargetFrameFromExtendedTime( epoch ) );
    data.records.erase( 57753.0 );
    settings->setEopData( data );
    rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( simulation_setup::createRotationModel( settings, "Earth" ) );
    const Time gap = Time( ( 57753.0 - 51544.5 ) * 86400.0 ) + 36.0;
    BOOST_CHECK_THROW( rotation->getRotationToTargetFrameFromExtendedTime( gap ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( factory_accepts_legacy_short_period_interpolation_settings_without_effect )
{
    auto settings = std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >(
            iau_2006, "GCRS", "/must-not-read-this-file", tai_scale );
    settings->setEopData( leapData( ) );
    auto referenceRotation = simulation_setup::createRotationModel( settings, "Earth" );
    settings->setShortTermInterpolatorSettings(
            interpolators::interpolatorGenerationSettings< double >( interpolators::linearInterpolation( ), 0.0, 86400.0, 3600.0 ) );
    auto rotation = simulation_setup::createRotationModel( settings, "Earth" );

    // Short-period corrections are evaluated at the continuous arguments, independently of these settings.
    const Time epoch = Time( ( 57753.0 - 51544.5 ) * 86400.0 ) + 1234.5;
    BOOST_CHECK_EQUAL( ( rotation->getRotationToTargetFrameFromExtendedTime( epoch ).toRotationMatrix( ) -
                         referenceRotation->getRotationToTargetFrameFromExtendedTime( epoch ).toRotationMatrix( ) )
                               .norm( ),
                       0.0 );
}

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
            BOOST_CHECK_SMALL( value.polarMotion.x( ) - query.fields[ 0 ] * arcsecond, 2e-20 );
            BOOST_CHECK_SMALL( value.polarMotion.y( ) - query.fields[ 1 ] * arcsecond, 2e-20 );
            BOOST_CHECK_SMALL( value.ut1MinusTai - query.fields[ 2 ], 1e-12 );
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

BOOST_AUTO_TEST_CASE( four_point_values_and_rates_reproduce_a_known_cubic_across_a_leap )
{
    auto data = leapData( );
    const Time firstTai = Time( ( 57752.0 - 51544.5 ) * 86400.0 ) + 36.0;
    for( auto& entry : data.records )
    {
        const double dat = entry.first < 57754.0 ? 36.0 : 37.0;
        const Time tai = Time( ( entry.first - 51544.5 ) * 86400.0 ) + dat;
        const double u = static_cast< double >( tai - firstTai ) / 86400.0;
        entry.second.polarMotion = { 1e-6 + 2e-9 * u - 3e-10 * u * u + 4e-11 * u * u * u, -2e-6 };
        entry.second.ut1MinusUtc = dat - 36.4 + 1e-3 * u - 2e-4 * u * u + 3e-5 * u * u * u;
        entry.second.celestialPoleOffsets = std::array< double, 2 >{ 1e-9 + 2e-11 * u * u * u, -2e-9 };
    }
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( data, fourPointSettings( ) );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
    for( const Time& tai : { firstTai, firstTai + 43200.0, leapMidnightTai( ) - 0.5, leapMidnightTai( ), evaluator->getEndTai( ) } )
    {
        const double u = static_cast< double >( tai - firstTai ) / 86400.0;
        const auto value = evaluator->evaluate( tai, true );
        BOOST_CHECK_SMALL( value.polarMotion.x( ) - ( 1e-6 + 2e-9 * u - 3e-10 * u * u + 4e-11 * u * u * u ), 1e-21 );
        BOOST_CHECK_SMALL( value.polarMotionRate->x( ) - ( 2e-9 - 6e-10 * u + 12e-11 * u * u ) / 86400.0, 1e-25 );
        BOOST_CHECK_SMALL( value.ut1MinusTai - ( -36.4 + 1e-3 * u - 2e-4 * u * u + 3e-5 * u * u * u ), 1e-12 );
        BOOST_CHECK_SMALL( *value.ut1MinusTaiRate - ( 1e-3 - 4e-4 * u + 9e-5 * u * u ) / 86400.0, 2e-19 );
        BOOST_CHECK_SMALL( value.celestialPoleOffsetRate->x( ) - 6e-11 * u * u / 86400.0, 1e-26 );
        const Time ut1 = converter.getCurrentTime< Time >( tai_scale, ut1_scale, tai );
        BOOST_CHECK_SMALL( static_cast< double >( converter.getCurrentTime< Time >( ut1_scale, tai_scale, ut1 ) - tai ), 2e-12 );
    }
}

BOOST_AUTO_TEST_CASE( four_point_stencils_reject_gaps_and_missing_cpo_without_reducing_order )
{
    auto data = interpolationData( interpolation_reference::datasets.front( ) );
    data.records.erase( 59003.0 );
    EarthOrientationDataEvaluator gaps( data, fourPointSettings( ) );
    const Time raw = Time( ( 59001.0 - 51544.5 ) * 86400.0 ) + 37.0;
    const auto exact = gaps.evaluate( raw, true );
    BOOST_CHECK_EQUAL( exact.polarMotion.x( ), data.records.at( 59001.0 ).polarMotion[ 0 ] );
    BOOST_CHECK( !exact.polarMotionRate );
    BOOST_CHECK( !exact.ut1MinusTaiRate );
    BOOST_CHECK_THROW( gaps.evaluate( raw + 43200.0 ), std::runtime_error );
    BOOST_CHECK_NO_THROW( gaps.evaluate( gaps.getEndTai( ) - 43200.0, true ) );
    data = interpolationData( interpolation_reference::datasets.front( ) );
    data.records.begin( )->second.celestialPoleOffsets.reset( );
    EarthOrientationDataEvaluator missing( data, fourPointSettings( ) );
    BOOST_CHECK_NO_THROW( missing.evaluate( raw + 43200.0 ) );
    BOOST_CHECK_THROW( missing.evaluate( raw + 43200.0, true ), std::runtime_error );
    BOOST_CHECK_NO_THROW( missing.evaluate( raw, true ) );
    BOOST_CHECK( !missing.evaluate( raw, true ).celestialPoleOffsetRate );
    auto zeroPolicy = fourPointSettings( );
    zeroPolicy.missingCelestialPole = MissingCelestialPolePolicy::zero;
    EarthOrientationDataEvaluator zero( data, zeroPolicy );
    BOOST_CHECK_EQUAL( zero.getCelestialPoleOffsets( raw + 43200.0 ).norm( ), 0.0 );
    data.records.erase( data.records.begin( ), std::prev( data.records.end( ), 3 ) );
    BOOST_CHECK_THROW( EarthOrientationDataEvaluator( data, fourPointSettings( ) ), std::runtime_error );
}

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

BOOST_AUTO_TEST_CASE( public_rotation_factory_selects_four_point_reference_policy )
{
    const auto& fixture = interpolation_reference::datasets.front( );
    auto settings = std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >(
            iau_2006, "GCRS", "/must-not-read-this-file", tai_scale );
    settings->setEopData( interpolationData( fixture ) );
    settings->setEopEvaluationSettings( fourPointSettings( ) );
    auto rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( simulation_setup::createRotationModel( settings, "Earth" ) );
    BOOST_REQUIRE( rotation );
    const auto evaluator = rotation->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( )->getEarthOrientationDataEvaluator( );
    BOOST_CHECK( evaluator->getSettings( ).interpolation == EarthOrientationInterpolation::lagrange_four_point );
    const auto& query = fixture.queries[ 4 ];
    const Time tai = evaluator->getStartTai( ) + query.taiSecondsSinceFirst;
    BOOST_CHECK_SMALL( evaluator->evaluate( tai ).ut1MinusTai - query.fields[ 2 ], 1e-12 );
    BOOST_CHECK_NO_THROW( rotation->getRotationToTargetFrameFromExtendedTime( tai ) );
    BOOST_CHECK_NO_THROW( rotation->getRotationToTargetFrameFromExtendedTime( evaluator->getStartTai( ) ) );
    BOOST_CHECK_NO_THROW( rotation->getRotationToTargetFrameFromExtendedTime( evaluator->getEndTai( ) ) );
}

BOOST_AUTO_TEST_CASE( rounded_tt_endpoints_use_endpoint_values_and_strict_polar_coverage )
{
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
    PolarMotionCalculator pole( nullptr, nullptr, evaluator );

    // The scalar TT CPO routes must use the same observation as polar motion.
    PrecessionNutationCalculator celestial( iau_2006, nullptr, nullptr, evaluator );

    for( const Time& boundary : { evaluator->getStartTai( ), evaluator->getEndTai( ) } )
    {
        // The rounded double TT of an end point lies within the coverage tolerance of that end point.
        const Time preciseTt = convertTAItoTT( boundary );
        const double tt = static_cast< double >( preciseTt );
        BOOST_CHECK_SMALL( static_cast< double >( converter.getTaiTime( tt_scale, tt ) - boundary ), 1e-6 );
        BOOST_CHECK_SMALL( ( pole.getPositionOfCipInItrsAtTt( tt ) - evaluator->evaluate( boundary ).polarMotion ).norm( ), 1e-20 );
        BOOST_CHECK_SMALL( ( pole.getPositionOfCipInItrsAtTt( preciseTt ) - evaluator->evaluate( boundary ).polarMotion ).norm( ), 1e-20 );
        BOOST_CHECK_SMALL( ( celestial.getPositionOfCipInGcrs( tt ) - celestial.getPositionOfCipInGcrsFromTai( tt, boundary ) ).norm( ),
                           1e-20 );
        BOOST_CHECK_SMALL(
                ( celestial.getPositionOfCipInGcrs( tt, 0.0 ) - celestial.getPositionOfCipInGcrsFromTai( tt, boundary ) ).norm( ), 1e-20 );

        // One millisecond outside the coverage is rejected.
        const double outside = boundary == evaluator->getStartTai( ) ? tt - 1.0E-3 : tt + 1.0E-3;
        BOOST_CHECK_THROW( evaluator->evaluate( convertTTtoTAI( Time( outside ) ) ), std::runtime_error );
        BOOST_CHECK_THROW( pole.getPositionOfCipInItrsAtTt( outside ), std::runtime_error );
        BOOST_CHECK_THROW( celestial.getPositionOfCipInGcrs( outside ), std::runtime_error );
        BOOST_CHECK_THROW( celestial.getPositionOfCipInGcrs( outside, 0.0 ), std::runtime_error );
    }
}

BOOST_AUTO_TEST_CASE( public_double_atomic_endpoints_reach_angles_rotation_and_convenience_methods )
{
    auto settings =
            std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >( iau_2006, "GCRS", "/must-not-read-this-file", tt_scale );
    settings->setEopData( leapData( ) );
    auto rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( simulation_setup::createRotationModel( settings, "Earth" ) );
    auto angles = rotation->getAnglesCalculator( );
    auto converter = angles->getTerrestrialTimeScaleConverter( );
    auto evaluator = converter->getEarthOrientationDataEvaluator( );
    for( const Time& boundary : { evaluator->getStartTai( ), evaluator->getEndTai( ) } )
    {
        const Time ut1 = converter->getUt1FromTai( boundary );
        const Time utc = EarthOrientationDataEvaluator::convertTaiToScalarUtc( boundary );
        for( TimeScales scale : { tt_scale, tdb_scale } )
        {
            const Time preciseInput = converter->getCurrentTime< Time >( tai_scale, scale, boundary );
            const double input = static_cast< double >( preciseInput );

            // The rounded double input lies within the coverage tolerance of the end point.
            BOOST_CHECK_SMALL( static_cast< double >( converter->getTaiTime( scale, input ) - boundary ), 1e-6 );
            BOOST_CHECK_SMALL( convertTimeScale( input, scale, tai_scale, converter ) - static_cast< double >( boundary ), 1e-6 );
            const auto result = angles->getRotationAnglesFromItrsToGcrs< double >( input, scale );
            const auto expected = angles->getRotationAnglesFromItrsToGcrs< Time >( boundary, tai_scale );
            BOOST_CHECK_SMALL( ( result.first - expected.first ).norm( ), 1e-16 );
            BOOST_CHECK_SMALL( result.second - static_cast< double >( ut1 ), 1e-6 );
            BOOST_CHECK_SMALL( converter->getUt1Correction( scale, input ) - static_cast< double >( ut1 - utc ), 1e-12 );
            BOOST_CHECK_SMALL(
                    converter->getCurrentTimeDifference( scale, ut1_scale, input ) - static_cast< double >( ut1 - Time( input ) ), 1e-6 );
            BOOST_CHECK_NO_THROW( angles->getPolarMotionCalculator( )->getPositionOfCipInItrsAtTt(
                    static_cast< double >( convertTAItoTT( boundary ) ) ) );

            // One millisecond outside the coverage is rejected by EOP consumers.
            const double outside = boundary == evaluator->getStartTai( ) ? input - 1.0E-3 : input + 1.0E-3;
            BOOST_CHECK_THROW( angles->getRotationAnglesFromItrsToGcrs< double >( outside, scale ), std::runtime_error );
            BOOST_CHECK_THROW( converter->getUt1Correction( scale, outside ), std::runtime_error );

            // Geocentric clock conversion remains available; only EOP consumers reject it.
            BOOST_CHECK_NO_THROW( converter->getTaiTime( scale, outside ) );
            if( scale == tt_scale )
            {
                BOOST_CHECK_NO_THROW( rotation->getRotationToBaseFrame( input ) );
                BOOST_CHECK_THROW( rotation->getRotationToBaseFrame( outside ), std::runtime_error );
            }
            else
            {
                auto tdbSettings = std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >(
                        iau_2006, "GCRS", "/must-not-read-this-file", tdb_scale );
                tdbSettings->setEopData( leapData( ) );
                auto tdbRotation = simulation_setup::createRotationModel( tdbSettings, "Earth" );
                BOOST_CHECK_NO_THROW( tdbRotation->getRotationToBaseFrame( input ) );
                BOOST_CHECK_THROW( tdbRotation->getRotationToBaseFrame( outside ), std::runtime_error );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( eop_free_geocentric_clock_converges_and_legacy_custom_ut1_stays_available )
{
    TerrestrialTimeScaleConverter atomic( nullptr, nullptr );
    for( const Time& tai : { Time( 536673637.0 ), Time( 1.0e9 ), Time( 3.0e8 ) } )
    {
        const Time tdb = atomic.getCurrentTime< Time >( tai_scale, tdb_scale, tai );
        BOOST_CHECK_SMALL( static_cast< double >( atomic.getTaiTime( tdb_scale, tdb ) - tai ), 1e-12 );
    }
    const double epoch = 1.0e8;
    auto daily = std::make_shared< interpolators::LinearInterpolator< double, double > >(
            std::map< double, double >{ { epoch - 86400.0, 0.25 }, { epoch + 86400.0, 0.25 } } );
    auto shortPeriod = std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< double > >(
            1.0, 0.0, std::vector< std::string >{}, std::vector< std::string >{}, []( double ) {
                return Eigen::Vector6d::Zero( ).eval( );
            } );
    TerrestrialTimeScaleConverter legacy( daily, shortPeriod );
    const double expectedUt1 = sofa_interface::convertTTtoUTC( epoch ) + 0.25;
    BOOST_CHECK_SMALL( legacy.getCurrentTime< double >( tt_scale, ut1_scale, epoch ) - expectedUt1, 1e-8 );
    BOOST_CHECK_SMALL( legacy.getUt1Correction( tt_scale, epoch ) - 0.25, 1e-15 );
    BOOST_CHECK_SMALL( legacy.getCurrentTimeDifference( tt_scale, ut1_scale, epoch ) -
                               static_cast< double >( legacy.getCurrentTime< Time >( tt_scale, ut1_scale, Time( epoch ) ) - Time( epoch ) ),
                       1e-15 );
}

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
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, interpolation, evaluator );
    const Time tdb = converter.getCurrentTime< Time >( tai_scale, tdb_scale, tai );

    // An interior epoch must not query clock values outside this intentionally narrow table.
    BOOST_CHECK_SMALL( static_cast< double >( converter.getTaiTime( tdb_scale, static_cast< double >( tdb ) ) - tai ), 1e-7 );
    BOOST_CHECK_SMALL( static_cast< double >( converter.getTaiTime( tdb_scale, tdb ) - tai ), 1e-12 );
    BOOST_CHECK_THROW( converter.getTaiTime( tdb_scale, static_cast< double >( tdb ) + 5000.0 ), std::runtime_error );
    BOOST_CHECK_EXCEPTION( converter.getGeocentricTtRateWithRespectToTdb( tdb ), std::runtime_error, []( const std::runtime_error& error ) {
        return std::string( error.what( ) ).find( "interpolator clock rates are unavailable" ) != std::string::npos;
    } );
}

BOOST_AUTO_TEST_CASE( native_geocentric_clock_rate_matches_independent_sofa_step_study )
{
    TerrestrialTimeScaleConverter converter( nullptr, nullptr );
    for( const Time& tai : { Time( 3.0e8 ), Time( 536673637.0 ), Time( 1.0e9 ) } )
    {
        const Time tt = convertTAItoTT( tai );
        const Time tdb = converter.getCurrentTime< Time >( tai_scale, tdb_scale, tai );
        const double actual = converter.getGeocentricTtRateWithRespectToTdb( tdb );
        for( double h : { 16.0, 32.0, 64.0 } )
        {
            const auto directOffset = [ &tt ]( double delta ) {
                return sofa_interface::getTDBminusTT( static_cast< double >( tt + delta ), 0.0, Eigen::Vector3d::Zero( ) );
            };
            const double derivative =
                    ( 8.0 * ( directOffset( h ) - directOffset( -h ) ) - ( directOffset( 2 * h ) - directOffset( -2 * h ) ) ) / ( 12 * h );
            const double expected = 1.0 / ( 1.0 + derivative );
            BOOST_CHECK_SMALL( actual - expected, 2e-15 );
            std::cout << std::setprecision( 17 ) << "SOFA clock epoch=" << static_cast< double >( tai ) << " h=" << h
                      << " dTT/dTDB=" << actual << " residual=" << actual - expected << '\n';
        }

        // At zero observer all fraction-of-day arguments give identical offsets.
        BOOST_CHECK_EQUAL( sofa_interface::getTDBminusTT( static_cast< double >( tt ), 0.0, Eigen::Vector3d::Zero( ) ),
                           sofa_interface::getTDBminusTT( static_cast< double >( tt ), 0.999, Eigen::Vector3d::Zero( ) ) );
    }
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    TerrestrialTimeScaleConverter withEop( nullptr, nullptr, nullptr, evaluator );
    const Time tai = evaluator->getStartTai( ) + 100000.0;
    const Time tdb = withEop.getCurrentTime< Time >( tai_scale, tdb_scale, tai );
    const double before = withEop.getGeocentricTtRateWithRespectToTdb( tdb );
    BOOST_CHECK_NO_THROW( withEop.getCurrentTime< Time >( tai_scale, tdb_scale, tai, Eigen::Vector3d( 6378137.0, 0.0, 0.0 ) ) );
    BOOST_CHECK_EQUAL( withEop.getGeocentricTtRateWithRespectToTdb( tdb ), before );
}

BOOST_AUTO_TEST_CASE( identity_time_conversions_preserve_values_without_clock_or_eop_coverage )
{
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( leapData( ) );
    const double centre = static_cast< double >( convertTAItoTT( evaluator->getStartTai( ) + 100000.0 ) );
    auto interpolation = interpolators::interpolatorGenerationSettings< double >(
            interpolators::linearInterpolation( interpolators::huntingAlgorithm, interpolators::throw_exception_at_boundary ),
            centre - 1000.0,
            centre + 1000.0,
            100.0 );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, interpolation, evaluator );
    const Time finalTt = convertTAItoTT( evaluator->getEndTai( ) );
    const Time finalTdb =
            TerrestrialTimeScaleConverter( nullptr, nullptr ).getCurrentTime< Time >( tai_scale, tdb_scale, evaluator->getEndTai( ) );
    for( double input : { static_cast< double >( finalTt ), static_cast< double >( finalTdb ), centre, 1e9 } )
    {
        for( TimeScales scale : { tai_scale, tt_scale, tdb_scale, utc_scale, ut1_scale } )
        {
            BOOST_CHECK_EQUAL( converter.getCurrentTime< double >( scale, scale, input ), input );
            BOOST_CHECK_EQUAL( converter.getCurrentTimeDifference< double >( scale, scale, input ), 0.0 );
            const Time precise = Time( input ) + 1e-12;
            BOOST_CHECK_EQUAL( converter.getCurrentTime< Time >( scale, scale, precise ), precise );
            BOOST_CHECK_EQUAL( converter.getCurrentTimeDifference< Time >( scale, scale, precise ), Time( 0.0 ) );
        }
    }
    const auto invalid = static_cast< TimeScales >( 99 );
    BOOST_CHECK_THROW( converter.getCurrentTime< double >( invalid, invalid, centre ), std::runtime_error );
    BOOST_CHECK_THROW( converter.getCurrentTimeDifference< Time >( invalid, invalid, Time( centre ) ), std::runtime_error );
    TerrestrialTimeScaleConverter legacy( nullptr, nullptr );
    BOOST_CHECK_EQUAL( legacy.getCurrentTime< double >( ut1_scale, ut1_scale, 1e9 ), 1e9 );
    BOOST_CHECK_EQUAL( legacy.getCurrentTimeDifference< double >( ut1_scale, ut1_scale, 1e9 ), 0.0 );
}

BOOST_AUTO_TEST_CASE( precise_time_clock_stencil_accepts_roundoff_at_a_fractional_binade_boundary )
{
    const Time tdb( 0, std::nextafter( 100.0L, std::numeric_limits< long double >::infinity( ) ) );

    // Adding 32/64 s crosses the 128-s fractional binade and rounds by one
    // long-double ULP (~7e-18 s), far below the precision of the returned rate.
    BOOST_CHECK_NE( tdb + 32.0 - tdb, Time( 32.0 ) );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr );
    const double actual = converter.getGeocentricTtRateWithRespectToTdb( tdb );
    const Time tt = convertTAItoTT( converter.getTaiTime( tdb_scale, tdb ) );
    const auto offset = [ &tt ]( double delta ) {
        return sofa_interface::getTDBminusTT( static_cast< double >( tt + delta ), 0.0, Eigen::Vector3d::Zero( ) );
    };
    const double derivative = ( 8.0 * ( offset( 32 ) - offset( -32 ) ) - ( offset( 64 ) - offset( -64 ) ) ) / 384.0;
    BOOST_CHECK_SMALL( actual - 1.0 / ( 1.0 + derivative ), 2e-15 );
}
