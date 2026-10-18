#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/ground_stations/stationHistoryMotionModel.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include <fstream>
#include <sstream>

using namespace tudat;
using namespace tudat::ground_stations;
using namespace tudat::input_output;
namespace
{
std::string fixture( const std::string& name )
{
    return paths::getTudatTestDataPath( ) + "pr1092/reference_files/station_history/" + name;
}
Time manufacturedTime( long double years )
{
    return Time( years * 31557600.0L );
}
StationHistoryMetadata metadata( )
{
    StationHistoryMetadata result;
    result.site = { "TEST", "A" };
    result.domesId = "00000M001";
    result.frame = "ITRF2020";
    result.release = "ITRF2020";
    result.origin = "CM";
    return result;
}
StationHistoryLinearSolution linear( )
{
    StationHistoryLinearSolution result;
    result.source.key = { metadata( ).site, "1" };
    result.referenceTime = manufacturedTime( 0 );
    return result;
}
StationPostSeismicTerm psd( StationPostSeismicFunction type, double amplitude, double tau, long double event = 0 )
{
    StationPostSeismicTerm result;
    result.function = type;
    result.eventTime = manufacturedTime( event );
    result.amplitudeMetres = Eigen::Vector3d( amplitude, -2 * amplitude, 3 * amplitude );
    result.relaxationJulianYears = tau;
    return result;
}
SinexStationHistory product( const std::string& name )
{
    auto result = readSinexStationHistory( fixture( name ) );
    result.frame = "ITRF2020";
    result.release = "ITRF2020";
    result.origin = "CM";
    result.tideSystem = "conventional product, supplied by caller";
    return result;
}
void closeVector( const Eigen::Vector3d& actual, const Eigen::Vector3d& expected, double tolerance )
{
    BOOST_CHECK_SMALL( ( actual - expected ).norm( ), tolerance );
}

// Controlled zero-coordinate companions exercise actual published coefficients
// without claiming an independent complete history for AB48/ANTC/ALBH.
std::pair< SinexStationHistory, SinexStationHistory > controlledCompanions( const SinexSiteRecord& site )
{
    auto coordinates = product( "ab07_solutions.snx" );
    auto intervals = product( "ab07_discontinuities.snx" );
    coordinates.sites.clear( );
    coordinates.sites.emplace( site.key, site );
    for( auto& entry : coordinates.parameters )
    {
        entry.second.key.site = site.key;
        entry.second.valueSi = 0.0;
        entry.second.originalValue = 0.0;
    }
    for( auto& span : coordinates.observationSpans )
    {
        span.key.site = site.key;
    }
    for( auto& interval : intervals.discontinuities )
    {
        interval.key.site = site.key;
    }
    return { coordinates, intervals };
}
}  // namespace

//! Compare all IGN post-seismic functions with the published scalar reference rows.
BOOST_AUTO_TEST_CASE( unchangedIgnParametricSelectorsAndSignedAmplitudes )
{
    std::ifstream stream( fixture( "ign_psd_reference.csv" ) );

    // Compare all IGN post-seismic functions with the published scalar reference rows.
    BOOST_REQUIRE( stream.good( ) );
    std::string row;
    std::getline( stream, row );
    int count = 0;
    while( std::getline( stream, row ) )
    {
        std::replace( row.begin( ), row.end( ), ',', ' ' );
        std::istringstream fields( row );
        int model;
        double dt, a1, t1, a2, t2, expected;
        BOOST_REQUIRE( fields >> model >> dt >> a1 >> t1 >> a2 >> t2 >> expected );
        std::vector< StationPostSeismicTerm > terms;
        if( model == 1 || model == 3 || model == 5 )
        {
            terms.push_back( psd( StationPostSeismicFunction::logarithmic, a1 * 0.001, t1 ) );
        }
        if( model == 2 || model == 4 )
        {
            terms.push_back( psd( StationPostSeismicFunction::exponential, a1 * 0.001, t1 ) );
        }
        if( model == 3 || model == 4 )
        {
            terms.push_back( psd( StationPostSeismicFunction::exponential, a2 * 0.001, t2 ) );
        }
        if( model == 5 )
        {
            terms.push_back( psd( StationPostSeismicFunction::logarithmic, a2 * 0.001, t2 ) );
        }
        StationCoordinateHistory history( { linear( ) }, terms, {}, metadata( ) );

        // The original LOG/EXP routine loses relative accuracy near zero. This
        // absolute bound covers its subtraction rounding, not a time-policy error.
        closeVector(
                history.getPosition( manufacturedTime( dt ) ), Eigen::Vector3d( expected, -2 * expected, 3 * expected ) * 0.001, 2e-15 );
        ++count;
    }
    BOOST_CHECK_EQUAL( count, 126 );
}

//! Check cumulative post-seismic displacements and analytic event velocities.
BOOST_AUTO_TEST_CASE( stableShortTimeCumulativeEventsAndAnalyticRates )
{
    const auto first = psd( StationPostSeismicFunction::exponential, 0.04, 0.1 );
    const auto second = psd( StationPostSeismicFunction::logarithmic, -0.02, 0.03, 0.5 );
    StationCoordinateHistory history( { linear( ) }, { first, second }, {}, metadata( ) );

    // Check cumulative post-seismic displacements and analytic event velocities.
    closeVector( history.getPosition( manufacturedTime( -0.001 ) ), Eigen::Vector3d::Zero( ), 1e-30 );
    closeVector( history.getPosition( manufacturedTime( 1e-18L ) ), first.amplitudeMetres * 1e-17, 1e-32 );

    // Independent multiplication/division order differs by a few double ULPs.
    closeVector( history.getState( manufacturedTime( 0 ) ).tail< 3 >( ), first.amplitudeMetres / ( 0.1 * 31557600.0 ), 1e-22 );
    for( const long double year : { 0.1L, 0.7L, 20.0L } )
    {
        const double h = 20.0;
        const Eigen::Vector3d difference = ( history.getPosition( manufacturedTime( year + h / 31557600.0L ) ) -
                                             history.getPosition( manufacturedTime( year - h / 31557600.0L ) ) ) /
                ( 2 * h );
        closeVector( history.getState( manufacturedTime( year ) ).tail< 3 >( ), difference, 1e-17 );
    }
    const auto eventState = history.getState( manufacturedTime( 0.5L ) );
    const Eigen::Vector3d expected = first.amplitudeMetres * ( std::exp( -5.0 ) / 0.1 ) + second.amplitudeMetres / 0.03;
    closeVector( eventState.tail< 3 >( ), expected / 31557600.0, 1e-22 );
}

//! Select the new coordinate solution and its velocity at a dated solution boundary.
BOOST_AUTO_TEST_CASE( halfOpenValidityRightHandSolutionsAndNoExtrapolation )
{
    auto first = linear( );
    first.validity = { manufacturedTime( 0 ), manufacturedTime( 1 ) };
    first.velocityMetresPerJulianYear = Eigen::Vector3d( 0.01, -0.02, 0.03 );
    auto second = first;
    second.source.key.solutionId = "2";
    second.positionMetres = Eigen::Vector3d( 0.1, -0.2, 0.3 );
    second.velocityMetresPerJulianYear *= -1;
    second.validity = { first.validity.endTai, manufacturedTime( 2 ) };
    StationCoordinateHistory history( { second, first }, {}, {}, metadata( ) );

    // Select the new coordinate solution and its velocity at a dated solution boundary.
    closeVector(
            history.getState( manufacturedTime( 1 ) ).head< 3 >( ), second.positionMetres + second.velocityMetresPerJulianYear, 1e-16 );
    closeVector( history.getState( manufacturedTime( 1 ) ).tail< 3 >( ), second.velocityMetresPerJulianYear / 31557600.0, 1e-24 );
}

//! Compare seasonal positions and velocities with calendar phases in rotated bases.
BOOST_AUTO_TEST_CASE( calendarSeasonalCoordinatesRatesAndBasisEquivalence )
{
    auto solution = linear( );
    solution.referenceTime = convertSinexCivilUtcToTai( { 2020, 1, 0 } );
    const auto meta = metadata( );
    const StationSeasonalTerm annual{ "1", 1, { 0.01, -0.02, 0.03 }, { -0.04, 0.05, -0.06 } };
    const StationSeasonalTerm semiannual{ "1", 2, { -0.002, 0.003, -0.004 }, { 0.005, -0.006, 0.007 } };
    StationCoordinateHistory history( { solution }, {}, { annual, semiannual }, meta );
    const auto time = convertSinexCivilUtcToTai( { 2020, 60, 43200 } );

    // Compare seasonal positions and velocities with calendar phases (of the TAI epoch, 37 s after UTC) in rotated bases.
    const double yearFraction = ( 59.5 + 37.0 / 86400.0 ) / 366;
    BOOST_CHECK_SMALL( static_cast< double >( getCalendarDecimalYear( time ) - 2020 ) - yearFraction, 1e-15 );
    Eigen::Vector3d expected = Eigen::Vector3d::Zero( );
    for( const auto& term : { annual, semiannual } )
    {
        const double phase = 2.0 * std::acos( -1.0 ) * term.cyclesPerCalendarYear * yearFraction;
        expected += term.cosineMetres * std::cos( phase ) + term.sineMetres * std::sin( phase );
    }
    closeVector( history.getPosition( time ), expected, 1e-16 );
    const Eigen::Vector3d difference = ( history.getPosition( time + 10.0 ) - history.getPosition( time - 10.0 ) ) / 20.0;
    closeVector( history.getState( time ).tail< 3 >( ), difference, 2e-18 );
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd( 0.73, Eigen::Vector3d( 1, 2, 3 ).normalized( ) ).toRotationMatrix( );
    auto rotatedAnnual = annual;
    auto rotatedSemiannual = semiannual;
    for( auto* term : { &rotatedAnnual, &rotatedSemiannual } )
    {
        term->cosineMetres = rotation * term->cosineMetres;
        term->sineMetres = rotation * term->sineMetres;
    }
    StationCoordinateHistory rotated( { solution }, {}, { rotatedAnnual, rotatedSemiannual }, meta );
    closeVector( rotated.getPosition( time ), rotation * history.getPosition( time ), 1e-16 );
    closeVector( rotated.getState( time ).tail< 3 >( ), rotation * history.getState( time ).tail< 3 >( ), 1e-23 );
}

//! Assemble published SINEX coordinates and their corresponding velocity records.
BOOST_AUTO_TEST_CASE( publishedHistoryPositionAndVelocity )
{
    auto solutions = product( "ab07_solutions.snx" );
    auto intervals = product( "ab07_discontinuities.snx" );
    auto deformation = product( "ab07_psd.snx" );
    const Eigen::Matrix3d basis = Eigen::Matrix3d::Identity( );
    auto history = createIgnStationCoordinateHistory( solutions, intervals, { "AB07", "A" }, &deformation, nullptr, basis );

    // Assemble published SINEX coordinates and their corresponding velocity records.
    BOOST_CHECK_EQUAL( history->getSolutions( ).size( ), 6 );
    BOOST_CHECK_EQUAL( history->getPostSeismicTerms( ).size( ), 2 );
    BOOST_CHECK_EQUAL( history->getMetadata( ).release, "ITRF2020" );
    const auto& last = history->getSolutions( ).back( );
    BOOST_CHECK_EQUAL( last.source.key.solutionId, "6" );

    // There is no V identifier 6: the SSC's complete solution 6 velocity is used.
    for( int c = 0; c < 3; ++c )
    {
        BOOST_CHECK_EQUAL( last.velocityMetresPerJulianYear[ c ], *last.source.velocityMetresPerSecond[ c ] * 31557600.0 );
    }
    const auto event = history->getPostSeismicTerms( ).front( ).eventTime;
    BOOST_CHECK( history->getState( event + 100.0 ).allFinite( ) );
}

//! Sum all published repeated post-seismic parameter pairs.
BOOST_AUTO_TEST_CASE( publishedRepeatedPsdPairs )
{
    const auto deformation = product( "repeated_psd.snx" );
    const Eigen::Matrix3d basis = Eigen::AngleAxisd( 0.7, Eigen::Vector3d( 1, 2, 3 ).normalized( ) ).toRotationMatrix( );
    for( const std::string siteCode : { "AB48", "ANTC" } )
    {
        const SinexSiteKey key{ siteCode, "A" };
        const auto companions = controlledCompanions( deformation.sites.find( key )->second );
        const auto history = createIgnStationCoordinateHistory( companions.first, companions.second, key, &deformation, nullptr, basis );

        // Sum all published repeated post-seismic parameter pairs.
        BOOST_CHECK_EQUAL( history->getPostSeismicTerms( ).size( ), siteCode == "AB48" ? 4 : 5 );
        const auto time = convertSinexCivilUtcToTai( { 2019, 160, 43200 } );
        Eigen::Vector3d expected = Eigen::Vector3d::Zero( );
        for( const auto& entry : deformation.parameters )
        {
            const auto& amplitude = entry.second;
            if( amplitude.key.site.siteCode != siteCode || amplitude.type[ 0 ] != 'A' )
            {
                continue;
            }
            const double tau = deformation.parameters.at( entry.first + 1 ).originalValue;
            const auto event = convertSinexCivilUtcToTai( *amplitude.referenceEpoch.civilUtc );
            const double dt = static_cast< double >( static_cast< long double >( time - event ) / 31557600.0L );
            const double correction = amplitude.type.substr( 1, 3 ) == "LOG" ? std::log( 1 + dt / tau ) : 1 - std::exp( -dt / tau );
            expected += basis.col( std::string( "ENU" ).find( amplitude.type.back( ) ) ) * ( amplitude.originalValue * correction );
        }
        closeVector( history->getPosition( time ), expected, 2e-15 );
    }
}

//! Build seasonal station motion from published SINEX coefficients in a fixed basis.
BOOST_AUTO_TEST_CASE( publishedSeasonalFactoryOriginIdentityAndFixedBasis )
{
    auto seasons = product( "albh_seasonal.snx" );
    const SinexSiteKey key{ "ALBH", "A" };
    const auto companions = controlledCompanions( seasons.sites.find( key )->second );
    const Eigen::Matrix3d basis = Eigen::AngleAxisd( -0.31, Eigen::Vector3d( 1, -2, 3 ).normalized( ) ).toRotationMatrix( );
    const auto history = createIgnStationCoordinateHistory( companions.first, companions.second, key, nullptr, &seasons, basis );

    // Build seasonal station motion from published SINEX coefficients in a fixed basis.
    BOOST_CHECK_EQUAL( history->getSeasonalTerms( ).size( ), 12 );
    BOOST_CHECK( history->hasFittedSeasonalTerms( ) );
    const auto time = convertSinexCivilUtcToTai( { 2020, 60, 43200 } );
    Eigen::Vector3d expected = Eigen::Vector3d::Zero( );
    for( const auto& entry : seasons.parameters )
    {
        const auto& p = entry.second;
        if( p.key.solutionId != "4" )
        {
            continue;
        }
        const double angle = 2 * std::acos( -1.0 ) * ( p.type[ 1 ] - '0' ) * ( 59.5 + 37.0 / 86400.0 ) / 366;
        const double function = p.type.substr( 2, 3 ) == "COS" ? std::cos( angle ) : std::sin( angle );
        expected += basis.col( std::string( "ENU" ).find( p.type.back( ) ) ) * ( p.originalValue * function );
    }
    closeVector( history->getPosition( time ), expected, 1e-16 );
    const auto withoutBasis =
            createIgnStationCoordinateHistory( companions.first, companions.second, key, nullptr, &seasons, Eigen::Matrix3d::Identity( ) );
    closeVector( history->getState( time ).tail< 3 >( ), basis * withoutBasis->getState( time ).tail< 3 >( ), 1e-23 );
    auto cf = seasons;
    cf.origin = "CF";

    // A seasonal product in CF cannot be combined silently with CM coordinates.
    BOOST_CHECK_THROW( createIgnStationCoordinateHistory( companions.first, companions.second, key, nullptr, &cf, basis ),
                       std::runtime_error );
    auto cfCoordinates = companions.first;
    auto cfIntervals = companions.second;
    cfCoordinates.origin = cfIntervals.origin = "CF";
    BOOST_CHECK_EQUAL( createIgnStationCoordinateHistory( cfCoordinates, cfIntervals, key, nullptr, &cf, basis )->getMetadata( ).origin,
                       "CF" );
}

//! Preserve precise solution selection through public station-motion dispatch.
BOOST_AUTO_TEST_CASE( precisePhysicalBoundaryIndependentSourceCoordinatesAndMotionIntegration )
{
    auto a = linear( );
    auto b = a;
    const Time boundary( 200000.0L );
    a.validity.endTai = boundary;
    b.validity.startTai = boundary;
    b.source.key.solutionId = "2";
    a.positionMetres = Eigen::Vector3d( 1, 2, 3 );
    b.positionMetres = Eigen::Vector3d( 4, 5, 6 );
    auto history = std::make_shared< StationCoordinateHistory >( std::vector< StationHistoryLinearSolution >{ a, b },
                                                                 std::vector< StationPostSeismicTerm >{},
                                                                 std::vector< StationSeasonalTerm >{},
                                                                 metadata( ) );

    // Select the old solution just before the boundary and the new solution at the boundary.
    closeVector( history->getPosition( boundary - 1e-12L ), a.positionMetres, 1e-30 );
    closeVector( history->getPosition( boundary ), b.positionMetres, 1e-30 );
    const auto settings = simulation_setup::stationHistoryGroundStationMotionSettings( history );
    const auto factoryStation = simulation_setup::createGroundStationState( std::make_shared< simulation_setup::Body >( ),
                                                                            Eigen::Vector3d( 100, 200, 300 ),
                                                                            coordinate_conversions::cartesian_position,
                                                                            { settings } );
    closeVector( factoryStation->getCartesianPositionInTime( 300000.0, "Earth" ), b.positionMetres, 1e-30 );
    closeVector( factoryStation->getCartesianStateInTime( 300000.0, "Earth" ).head< 3 >( ), b.positionMetres, 1e-30 );
}
