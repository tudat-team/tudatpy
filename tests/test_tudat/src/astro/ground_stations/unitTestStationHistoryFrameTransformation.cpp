#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/ground_stations/stationCoordinateHistory.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/reference_frames/terrestrialFrameTransformation.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "pr1092/reference_files/itrf_transformations/itrf_reference.h"
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
    result.sources = { "controlled source history" };
    result.signals = { "linear marker coordinates" };
    return result;
}
//! TAI epoch with the given calendar decimal year (for years with TAI - UTC = 37 s).
Time time( double year )
{
    const int calendarYear = static_cast< int >( std::floor( year ) );
    const Time result = convertSinexCivilUtcToTai( { calendarYear, 1, 0 } ) - 37.0 +
            ( year - calendarYear ) * ( calendarYear % 4 == 0 ? 366.0 : 365.0 ) * 86400.0;
    BOOST_REQUIRE_SMALL( static_cast< double >( getCalendarDecimalYear( result ) ) - year, 1e-12 );
    return result;
}

//! History with the given state (velocity in metres per second) at its reference epoch.
std::shared_ptr< StationCoordinateHistory > history( const Eigen::Vector6d& state,
                                                     const Time& reference,
                                                     StationHistoryMetadata meta = metadata( ),
                                                     bool seasonal = false )
{
    StationHistoryLinearSolution linear;
    linear.source.key = { meta.site, "1" };
    linear.referenceTime = reference;
    linear.positionMetres = state.head< 3 >( );
    linear.velocityMetresPerJulianYear = state.tail< 3 >( ) * 31557600.0;
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

//! Compare transformed history positions with independent ITRF forward and inverse fixtures, and the velocities with
//! the transformation model (the fixtures use a different decimal-year rate).
BOOST_AUTO_TEST_CASE( catalogueAdapterIndependentPublishedForwardAndInverseFixtures )
{
    const Eigen::Vector6d input = Eigen::Map< const Eigen::Vector6d >( itrfReferenceInputState );
    const reference_frames::HelmertTransformationModel model( readOfficialItrf2020TransformationData( ) );
    int rows = 0;
    for( const auto& fixture : itrfReferenceFixtures )
    {
        if( fixture.catalogue != 2020 || fixture.year != 2025.25 )
        {
            continue;
        }
        const auto query = time( fixture.year );
        const auto forward = transform( history( input, query ), fixture.frame );
        const auto inverse = transform( history( input, query, metadata( fixture.frame ) ), "ITRF2020" );

        // Compare both directions with the frozen reference positions, and with the model states.
        const Eigen::Vector6d forwardState = forward->getState( query );
        const Eigen::Vector6d inverseState = inverse->getState( query );
        BOOST_CHECK_SMALL( ( forwardState.head< 3 >( ) - Eigen::Map< const Eigen::Vector6d >( fixture.forward ).head< 3 >( ) ).norm( ),
                           2e-9 );
        BOOST_CHECK_SMALL( ( inverseState.head< 3 >( ) - Eigen::Map< const Eigen::Vector6d >( fixture.inverse ).head< 3 >( ) ).norm( ),
                           2e-9 );
        closeState( forwardState, model.transformState( input, "ITRF2020", fixture.frame, fixture.year, 1.0 / 31557600.0 ) );
        closeState( inverseState, model.transformState( input, fixture.frame, "ITRF2020", fixture.year, 1.0 / 31557600.0 ) );
        ++rows;
    }
    BOOST_CHECK_EQUAL( rows, 13 );
}

//! Check transformed-coefficient velocity against finite differences on the native calendar clock.
BOOST_AUTO_TEST_CASE( sharedNativeCalendarClockAndCompleteVelocity )
{
    Eigen::Vector6d initial;
    initial << 4e6, -2e6, 5e6, -0.05 / 31557600.0, 0.03 / 31557600.0, 0.01 / 31557600.0;
    const auto source = history( initial, convertSinexCivilUtcToTai( { 2020, 1, 0 } ) );
    const auto transformed = transform( source );
    const Time query = convertSinexCivilUtcToTai( { 2020, 60, 43200 } );
    const double h = 300000;
    const Eigen::Vector3d finiteDifference = ( transformed->getPosition( query + h ) - transformed->getPosition( query - h ) ) / ( 2 * h );
    BOOST_CHECK_SMALL( ( transformed->getState( query ).tail< 3 >( ) - finiteDifference ).norm( ), 1e-14 );

    // Away from the reference epoch, the transformed coefficients differ from the transformation evaluated at the calendar
    // year of the epoch by the difference between calendar and Julian years (a few micrometres at most).
    const reference_frames::HelmertTransformationModel model( readOfficialItrf2020TransformationData( ) );
    const Eigen::Vector3d exact = model.transformPosition(
            source->getPosition( query ), "ITRF2020", "ITRF93", static_cast< double >( getCalendarDecimalYear( query ) ) );
    BOOST_CHECK_SMALL( ( transformed->getPosition( query ) - exact ).norm( ), 1e-5 );
}

//! Require explicit release alignment and preserve origin and tide conventions.
BOOST_AUTO_TEST_CASE( sameFrameIdentityExplicitReleaseAlignmentAndOriginTidePreservation )
{
    const auto query = time( 2020.25 );
    const auto source = history( Eigen::Vector6d::Zero( ), query );
    const auto catalogue = readOfficialItrf2020TransformationData( );
    const auto& meta = source->getMetadata( );

    // Require explicit release alignment and preserve origin and tide conventions.

    // A realization label alone does not authorize alignment of different catalogue releases.
    BOOST_CHECK_THROW(
            transformStationCoordinateHistory( source, catalogue, { "ITRF2020", "ITRF2020-u2024" }, meta.origin, meta.tideSystem ),
            std::runtime_error );
    const auto aligned =
            transformStationCoordinateHistory( source, catalogue, { "ITRF2020", "ITRF2020-u2024" }, meta.origin, meta.tideSystem, true );
    closeState( aligned->getState( query ), source->getState( query ) );

    // Identity transformations must still preserve physical origin and tide conventions.
    BOOST_CHECK_THROW( transformStationCoordinateHistory( source, catalogue, { "ITRF2020", "ITRF2020" }, "CF", meta.tideSystem, true ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            transformStationCoordinateHistory( source, catalogue, { "ITRF2020", "ITRF2020" }, meta.origin, "different tides", true ),
            std::runtime_error );
}

//! Apply transformed coordinate histories through the public motion factory.
BOOST_AUTO_TEST_CASE( transformedHistoryThroughPublicMotionFactory )
{
    const auto query = time( 2020.25 );
    Eigen::Vector6d input;
    input << 4e6, -2e6, 5e6, 0.01, 0.02, -0.03;
    const auto transformed = transform( history( input, query ) );
    const auto settings = simulation_setup::stationHistoryGroundStationMotionSettings( transformed );
    const auto station = simulation_setup::createGroundStationState( std::make_shared< simulation_setup::Body >( ),
                                                                     Eigen::Vector3d( 10, 20, 30 ),
                                                                     coordinate_conversions::cartesian_position,
                                                                     { settings } );

    // Apply transformed coordinate histories through the public motion factory, at the TAI epoch of the TDB epoch.
    const Time tai = basic_astrodynamics::getGeocentricTimeConverter( )->getTaiTime( basic_astrodynamics::tdb_scale, query );
    closeState( station->getCartesianStateInTime( query, "Earth" ), transformed->getState( tai ) );
    BOOST_CHECK_SMALL( ( station->getCartesianPositionInTime( query, "Earth" ) - transformed->getPosition( tai ) ).norm( ), 1e-9 );
}
