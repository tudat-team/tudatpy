#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/ground_stations/stationCoordinateHistory.h"
#include "tudat/astro/reference_frames/terrestrialFrameTransformation.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "../../../reference_files/itrf_transformations/itrf_reference.h"
#include <limits>

using namespace tudat;
using namespace tudat::ground_stations;
using namespace tudat::input_output;
namespace
{
StationHistoryMetadata metadata( const std::string& frame = "ITRF2020" )
{
    StationHistoryMetadata result;
    result.site = { "TEST", "A" };
    result.domesId = "00000M001";
    result.frame = result.release = frame;
    result.origin = "CM";
    result.tideSystem = "explicit controlled tide convention";
    result.timeConvention = "controlled affine coordinate";
    result.sources = { "controlled source history" };
    result.signals = { "linear marker coordinates" };
    return result;
}
StationHistoryTime time( double year, double yearRate )
{
    StationHistoryTime result;
    result.calendarYear = static_cast< int >( std::floor( year ) );
    result.yearFraction = year - result.calendarYear;
    result.mjdDay = 51544;

    // The MJD rate matches the year rate, so that elapsed Julian years and calendar years advance at the same rate.
    result.mjdRate = yearRate * 365.25;
    result.decimalYearRate = yearRate;
    result.convention = metadata( ).timeConvention;
    return result;
}

//! History with the given state (velocity in metres per second) at its reference epoch; the velocity is converted to
//! metres per Julian year with the year rate of the reference epoch.
std::shared_ptr< StationCoordinateHistory > history( const Eigen::Vector6d& state,
                                                     const StationHistoryTime& reference,
                                                     StationHistoryMetadata meta = metadata( ),
                                                     bool seasonal = false )
{
    StationHistoryLinearSolution linear;
    linear.source.key = { meta.site, "1" };
    linear.referenceTime = reference;
    linear.positionMetres = state.head< 3 >( );
    linear.velocityMetresPerJulianYear = state.tail< 3 >( ) * 365.25 / reference.mjdRate.value_or( 1.0 / 86400.0 );
    std::vector< StationSeasonalTerm > terms;
    if( seasonal )
    {
        terms.push_back( { "1", 1, Eigen::Vector3d( 0.01, -0.02, 0.03 ), Eigen::Vector3d::Zero( ), {} } );
    }
    return std::make_shared< StationCoordinateHistory >(
            std::vector< StationHistoryLinearSolution >{ linear }, std::vector< StationPostSeismicTerm >{}, terms, meta );
}
void closeState( const Eigen::Vector6d& actual, const Eigen::Vector6d& expected )
{
    for( int i = 0; i < 6; ++i )
    {
        BOOST_CHECK_SMALL( actual[ i ] - expected[ i ], i < 3 ? 2e-9 : 2e-17 );
    }
}
std::shared_ptr< StationCoordinateHistory > transform( std::shared_ptr< StationCoordinateHistory > source,
                                                       const std::string& target = "ITRF93" )
{
    const auto& meta = source->getMetadata( );
    return transformStationCoordinateHistory(
            source, readOfficialItrf2020TransformationData( ), { target, target }, meta.origin, meta.tideSystem );
}
}  // namespace

BOOST_AUTO_TEST_CASE( catalogueAdapterIndependentPublishedForwardAndInverseFixtures )
{
    const Eigen::Vector6d input = Eigen::Map< const Eigen::Vector6d >( itrfReferenceInputState );
    int rows = 0;
    for( const auto& fixture : itrfReferenceFixtures )
    {
        if( fixture.catalogue != 2020 || fixture.year != 2025.25 )
        {
            continue;
        }
        const auto query = time( fixture.year, fixture.yearRate );
        const auto forward = transform( history( input, query ), fixture.frame );
        closeState( forward->getState( query ), Eigen::Map< const Eigen::Vector6d >( fixture.forward ) );
        const auto inverse = transform( history( input, query, metadata( fixture.frame ) ), "ITRF2020" );
        closeState( inverse->getState( query ), Eigen::Map< const Eigen::Vector6d >( fixture.inverse ) );
        auto noRates = query;
        noRates.mjdRate.reset( );
        noRates.decimalYearRate.reset( );
        BOOST_CHECK_SMALL(
                ( forward->getPosition( noRates ) - Eigen::Map< const Eigen::Vector6d >( fixture.forward ).head< 3 >( ) ).norm( ), 2e-9 );
        BOOST_CHECK_THROW( forward->getState( noRates ), std::runtime_error );
        ++rows;
    }
    BOOST_CHECK_EQUAL( rows, 13 );
}

BOOST_AUTO_TEST_CASE( sharedNativeCalendarClockAndCompleteAffineVelocity )
{
    auto ref = getIgnNominalUtcStationTime( SinexCivilUtc{ 2020, 1, 0 } );
    auto meta = metadata( );
    meta.timeConvention = ref.convention;
    Eigen::Vector6d initial;
    initial << 4e6, -2e6, 5e6, -0.05 / 31557600.0, 0.03 / 31557600.0, 0.01 / 31557600.0;
    const auto source = history( initial, ref, meta );
    const auto transformed = transform( source );
    const Time tai = convertSinexCivilUtcToTai( { 2020, 60, 43200 } );
    const auto query = getIgnNominalUtcStationTimeFromTai( tai, 1.0 );

    // Calendar year and its TDB rate of the station time.
    BOOST_CHECK_EQUAL( static_cast< double >( query.calendarYear + query.yearFraction ), 2020 + 59.5 / 366 );
    BOOST_CHECK_EQUAL( *query.decimalYearRate, 1.0 / ( 366 * 86400.0 ) );
    const double h = 300000;
    const Eigen::Vector3d finiteDifference = ( transformed->getPosition( getIgnNominalUtcStationTimeFromTai( tai + h ) ) -
                                               transformed->getPosition( getIgnNominalUtcStationTimeFromTai( tai - h ) ) ) /
            ( 2 * h );
    BOOST_CHECK_SMALL( ( transformed->getState( query ).tail< 3 >( ) - finiteDifference ).norm( ), 1e-14 );

    // Away from the reference epoch, the transformed coefficients differ from the transformation evaluated at the calendar
    // year of the epoch by the difference between calendar and Julian years (a few micrometres at most).
    const reference_frames::HelmertTransformationModel model( readOfficialItrf2020TransformationData( ) );
    const Eigen::Vector3d exact = model.transformPosition(
            source->getPosition( query ), "ITRF2020", "ITRF93", static_cast< double >( query.calendarYear + query.yearFraction ) );
    BOOST_CHECK_SMALL( ( transformed->getPosition( query ) - exact ).norm( ), 1e-5 );
    auto missingCalendarRate = query;
    missingCalendarRate.decimalYearRate.reset( );
    BOOST_CHECK_NO_THROW( transformed->getPosition( missingCalendarRate ) );

    // The transformed coefficients are linear in elapsed Julian years, so the state needs only the MJD rate.
    BOOST_CHECK_NO_THROW( transformed->getState( missingCalendarRate ) );
    const auto nextYear = getIgnNominalUtcStationTimeFromTai( convertSinexCivilUtcToTai( { 2021, 1, 0 } ), 1.0 );
    BOOST_CHECK_EQUAL( *nextYear.decimalYearRate, 1.0 / ( 365 * 86400.0 ) );
}

BOOST_AUTO_TEST_CASE( sameFrameIdentityExplicitReleaseAlignmentAndOriginTidePreservation )
{
    const auto query = time( 2020.25, 1.0 / 31557600.0 );
    const auto source = history( Eigen::Vector6d::Zero( ), query );
    const auto catalogue = readOfficialItrf2020TransformationData( );
    const auto& meta = source->getMetadata( );
    BOOST_CHECK_THROW(
            transformStationCoordinateHistory( source, catalogue, { "ITRF2020", "ITRF2020-u2024" }, meta.origin, meta.tideSystem ),
            std::runtime_error );
    const auto aligned =
            transformStationCoordinateHistory( source, catalogue, { "ITRF2020", "ITRF2020-u2024" }, meta.origin, meta.tideSystem, true );
    auto queryWithoutCalendarRate = query;
    queryWithoutCalendarRate.decimalYearRate.reset( );
    closeState( aligned->getState( queryWithoutCalendarRate ), source->getState( queryWithoutCalendarRate ) );
    BOOST_CHECK_EQUAL( aligned->getMetadata( ).release, "ITRF2020-u2024" );
    BOOST_CHECK_EQUAL( source->getMetadata( ).release, "ITRF2020" );
    BOOST_CHECK_EQUAL( aligned->getMetadata( ).origin, meta.origin );
    BOOST_CHECK_EQUAL( aligned->getMetadata( ).tideSystem, meta.tideSystem );
    BOOST_CHECK_EQUAL( aligned->getMetadata( ).signals.back( ), "declared same-frame identity" );
    BOOST_CHECK_THROW( transformStationCoordinateHistory( source, catalogue, { "ITRF2020", "ITRF2020" }, "CF", meta.tideSystem, true ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            transformStationCoordinateHistory( source, catalogue, { "ITRF2020", "ITRF2020" }, meta.origin, "different tides", true ),
            std::runtime_error );
    BOOST_CHECK_THROW( transform( source, "unavailable realization" ), std::runtime_error );
    auto updatedMeta = meta;
    updatedMeta.release = "ITRF2020-u2024";
    const auto updated = history( Eigen::Vector6d::Zero( ), query, updatedMeta, true );
    BOOST_CHECK_THROW( transform( updated, "ITRF2014" ), std::runtime_error );
    const auto acknowledged =
            transformStationCoordinateHistory( updated, catalogue, { "ITRF2014", "ITRF2014" }, meta.origin, meta.tideSystem, true );
    BOOST_CHECK( acknowledged->hasFittedSeasonalTerms( ) );
    BOOST_CHECK_EQUAL( acknowledged->getMetadata( ).sources.front( ), meta.sources.front( ) );
    BOOST_CHECK_EQUAL( acknowledged->getMetadata( ).sources.back( ), catalogue.provenance.source );
    BOOST_CHECK_THROW( transformStationCoordinateHistory( nullptr, catalogue, { "ITRF2014", "ITRF2014" }, meta.origin, meta.tideSystem ),
                       std::runtime_error );
}

BOOST_AUTO_TEST_CASE( transformedHistoryThroughPublicMotionFactory )
{
    const auto query = time( 2020.25, 1.0 / 31557600.0 );
    Eigen::Vector6d input;
    input << 4e6, -2e6, 5e6, 0.01, 0.02, -0.03;
    const auto transformed = transform( history( input, query ) );
    const auto settings =
            simulation_setup::stationHistoryGroundStationMotionSettings( transformed, {}, [ query ]( const auto&, bool rates ) {
                auto output = query;
                if( !rates )
                {
                    output.mjdRate.reset( );
                    output.decimalYearRate.reset( );
                }
                return output;
            } );
    const auto station = simulation_setup::createGroundStationState( std::make_shared< simulation_setup::Body >( ),
                                                                     Eigen::Vector3d( 10, 20, 30 ),
                                                                     coordinate_conversions::cartesian_position,
                                                                     { settings } );
    closeState( station->getCartesianStateInTime( 0, "Earth" ), transformed->getState( query ) );
    BOOST_CHECK_SMALL( ( station->getCartesianPositionInTime( 0, "Earth" ) - transformed->getPosition( query ) ).norm( ), 1e-9 );
}
