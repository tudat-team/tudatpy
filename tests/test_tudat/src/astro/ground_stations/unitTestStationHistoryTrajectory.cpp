#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
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
    const std::string source = __FILE__;
    return source.substr( 0, source.find_last_of( "/\\" ) ) + "/../../../reference_files/station_history/" + name;
}
StationHistoryTime manufacturedTime( long double years )
{
    const long double days = years * 365.25L;
    StationHistoryTime result;
    result.tai = Time( years * 31557600.0L );
    result.mjdDay = 50000 + static_cast< int >( std::floor( days ) );
    result.mjdFraction = days - std::floor( days );
    result.calendarYear = 2000 + static_cast< int >( std::floor( years ) );
    result.yearFraction = years - std::floor( years );
    result.mjdRate = 1.0 / 86400.0;
    result.decimalYearRate = 1.0 / 31557600.0;
    result.convention = "controlled source coordinate";
    return result;
}
StationHistoryMetadata metadata( )
{
    StationHistoryMetadata result;
    result.site = { "TEST", "A" };
    result.domesId = "00000M001";
    result.frame = "ITRF2020";
    result.release = "ITRF2020";
    result.origin = "CM";
    result.timeConvention = "controlled source coordinate";
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
// without claiming an independent complete trajectory for AB48/ANTC/ALBH.
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
    for( auto& span : coordinates.observationSpans ) span.key.site = site.key;
    for( auto& interval : intervals.discontinuities ) interval.key.site = site.key;
    return { coordinates, intervals };
}
}  // namespace

BOOST_AUTO_TEST_CASE( unchangedIgnParametricSelectorsAndSignedAmplitudes )
{
    std::ifstream stream( fixture( "ign_psd_reference.csv" ) );
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
        if( model == 1 || model == 3 || model == 5 ) terms.push_back( psd( StationPostSeismicFunction::logarithmic, a1 * 0.001, t1 ) );
        if( model == 2 || model == 4 ) terms.push_back( psd( StationPostSeismicFunction::exponential, a1 * 0.001, t1 ) );
        if( model == 3 || model == 4 ) terms.push_back( psd( StationPostSeismicFunction::exponential, a2 * 0.001, t2 ) );
        if( model == 5 ) terms.push_back( psd( StationPostSeismicFunction::logarithmic, a2 * 0.001, t2 ) );
        StationHistoryTrajectory trajectory( { linear( ) }, terms, {}, metadata( ) );
        // The original LOG/EXP routine loses relative accuracy near zero. This
        // absolute bound covers its subtraction rounding, not a time-policy error.
        closeVector(
                trajectory.getPosition( manufacturedTime( dt ) ), Eigen::Vector3d( expected, -2 * expected, 3 * expected ) * 0.001, 2e-15 );
        ++count;
    }
    BOOST_CHECK_EQUAL( count, 126 );
}

BOOST_AUTO_TEST_CASE( stableShortTimeCumulativeEventsAndAnalyticRates )
{
    const auto first = psd( StationPostSeismicFunction::exponential, 0.04, 0.1 );
    const auto second = psd( StationPostSeismicFunction::logarithmic, -0.02, 0.03, 0.5 );
    StationHistoryTrajectory trajectory( { linear( ) }, { first, second }, {}, metadata( ) );
    closeVector( trajectory.getPosition( manufacturedTime( -0.001 ) ), Eigen::Vector3d::Zero( ), 1e-30 );
    closeVector( trajectory.getPosition( manufacturedTime( 1e-18L ) ), first.amplitudeMetres * 1e-17, 1e-32 );
    // Independent multiplication/division order differs by a few double ULPs.
    closeVector( trajectory.getState( manufacturedTime( 0 ) ).tail< 3 >( ), first.amplitudeMetres / ( 0.1 * 31557600.0 ), 1e-22 );
    for( const long double year : { 0.1L, 0.7L, 20.0L } )
    {
        const double h = 20.0;
        const Eigen::Vector3d difference = ( trajectory.getPosition( manufacturedTime( year + h / 31557600.0L ) ) -
                                             trajectory.getPosition( manufacturedTime( year - h / 31557600.0L ) ) ) /
                ( 2 * h );
        closeVector( trajectory.getState( manufacturedTime( year ) ).tail< 3 >( ), difference, 1e-17 );
    }
    const auto eventState = trajectory.getState( manufacturedTime( 0.5L ) );
    const Eigen::Vector3d expected = first.amplitudeMetres * ( std::exp( -5.0 ) / 0.1 ) + second.amplitudeMetres / 0.03;
    closeVector( eventState.tail< 3 >( ), expected / 31557600.0, 1e-22 );
    auto noRate = manufacturedTime( 0.7L );
    noRate.mjdRate.reset( );
    BOOST_CHECK_NO_THROW( trajectory.getPosition( noRate ) );
    BOOST_CHECK_THROW( trajectory.getState( noRate ), std::runtime_error );
    auto bad = first;
    bad.relaxationJulianYears = 0.0;
    BOOST_CHECK_THROW( StationHistoryTrajectory( { linear( ) }, { bad }, {}, metadata( ) ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( halfOpenValidityRightHandSolutionsAndNoExtrapolation )
{
    auto first = linear( );
    first.validity = { manufacturedTime( 0 ).tai, manufacturedTime( 1 ).tai };
    first.velocityMetresPerJulianYear = Eigen::Vector3d( 0.01, -0.02, 0.03 );
    auto second = first;
    second.source.key.solutionId = "2";
    second.positionMetres = Eigen::Vector3d( 0.1, -0.2, 0.3 );
    second.velocityMetresPerJulianYear *= -1;
    second.validity = { first.validity.endTai, manufacturedTime( 2 ).tai };
    StationHistoryTrajectory trajectory( { second, first }, {}, {}, metadata( ) );
    closeVector(
            trajectory.getState( manufacturedTime( 1 ) ).head< 3 >( ), second.positionMetres + second.velocityMetresPerJulianYear, 1e-16 );
    closeVector( trajectory.getState( manufacturedTime( 1 ) ).tail< 3 >( ), second.velocityMetresPerJulianYear / 31557600.0, 1e-24 );
    BOOST_CHECK_THROW( trajectory.getPosition( manufacturedTime( -1e-9 ) ), std::runtime_error );
    BOOST_CHECK_THROW( trajectory.getPosition( manufacturedTime( 2 ) ), std::runtime_error );
    second.validity.startTai = manufacturedTime( 1.1 ).tai;
    BOOST_CHECK_THROW( StationHistoryTrajectory( { first, second }, {}, {}, metadata( ) ), std::runtime_error );
    second.validity.startTai = manufacturedTime( 0.9 ).tai;
    BOOST_CHECK_THROW( StationHistoryTrajectory( { first, second }, {}, {}, metadata( ) ), std::runtime_error );
    auto wrong = manufacturedTime( 0.2 );
    wrong.convention = "unannounced TT substitution";
    BOOST_CHECK_THROW( trajectory.getPosition( wrong ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( calendarSeasonalCoordinatesRatesAndBasisEquivalence )
{
    auto ref = getIgnNominalUtcStationTime( SinexCivilUtc{ 2020, 1, 0 } );
    auto solution = linear( );
    solution.referenceTime = ref;
    auto meta = metadata( );
    meta.timeConvention = ref.convention;
    const StationSeasonalTerm annual{ "1", 1, { 0.01, -0.02, 0.03 }, { -0.04, 0.05, -0.06 }, {} };
    const StationSeasonalTerm semiannual{ "1", 2, { -0.002, 0.003, -0.004 }, { 0.005, -0.006, 0.007 }, {} };
    StationHistoryTrajectory trajectory( { solution }, {}, { annual, semiannual }, meta );
    const auto tai = convertSinexCivilUtcToTai( { 2020, 60, 43200 } );
    const auto time = getIgnNominalUtcStationTimeFromTai( tai, 1.00000001 );
    BOOST_CHECK_SMALL( static_cast< double >( time.yearFraction - 59.5L / 366 ), 1e-18 );
    Eigen::Vector3d expected = Eigen::Vector3d::Zero( );
    for( const auto& term : { annual, semiannual } )
    {
        const double phase = 2.0 * std::acos( -1.0 ) * term.cyclesPerCalendarYear * ( 59.5 / 366 );
        expected += term.cosineMetres * std::cos( phase ) + term.sineMetres * std::sin( phase );
    }
    closeVector( trajectory.getPosition( time ), expected, 1e-16 );
    const Eigen::Vector3d difference = ( trajectory.getPosition( getIgnNominalUtcStationTimeFromTai( tai + 10.0 ) ) -
                                         trajectory.getPosition( getIgnNominalUtcStationTimeFromTai( tai - 10.0 ) ) ) /
            20.0 * 1.00000001;
    closeVector( trajectory.getState( time ).tail< 3 >( ), difference, 2e-18 );
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd( 0.73, Eigen::Vector3d( 1, 2, 3 ).normalized( ) ).toRotationMatrix( );
    auto rotatedAnnual = annual;
    auto rotatedSemiannual = semiannual;
    for( auto* term : { &rotatedAnnual, &rotatedSemiannual } )
    {
        term->cosineMetres = rotation * term->cosineMetres;
        term->sineMetres = rotation * term->sineMetres;
    }
    StationHistoryTrajectory rotated( { solution }, {}, { rotatedAnnual, rotatedSemiannual }, meta );
    closeVector( rotated.getPosition( time ), rotation * trajectory.getPosition( time ), 1e-16 );
    closeVector( rotated.getState( time ).tail< 3 >( ), rotation * trajectory.getState( time ).tail< 3 >( ), 1e-23 );
    const auto nextYear = getIgnNominalUtcStationTimeFromTai( convertSinexCivilUtcToTai( { 2021, 1, 0 } ), 1.0 );
    BOOST_CHECK_EQUAL( nextYear.calendarYear, 2021 );
    BOOST_CHECK_EQUAL( nextYear.yearFraction, 0.0L );
    BOOST_CHECK_EQUAL( *nextYear.decimalYearRate, 1.0 / ( 365 * 86400.0 ) );
    auto missingRate = time;
    missingRate.decimalYearRate.reset( );
    BOOST_CHECK_NO_THROW( trajectory.getPosition( missingRate ) );
    BOOST_CHECK_THROW( trajectory.getState( missingRate ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( nativeCivilUtcLeapSecondAndSofaReference )
{
    const auto leap = convertSinexCivilUtcToTai( { 2016, 366, 86400 } );
    const auto midnight = convertSinexCivilUtcToTai( { 2017, 1, 0 } );
    BOOST_CHECK_EQUAL( static_cast< double >( midnight - leap ), 1.0 );
    double u1, u2, t1, t2;
    BOOST_REQUIRE_EQUAL( iauDtf2d( "UTC", 2016, 12, 31, 23, 59, 60.0, &u1, &u2 ), 0 );
    BOOST_REQUIRE_EQUAL( iauUtctai( u1, u2, &t1, &t2 ), 0 );
    const long double directTai = ( static_cast< long double >( t1 - 2451545.0 ) + t2 ) * 86400.0L;
    BOOST_CHECK_SMALL( static_cast< double >( static_cast< long double >( leap ) - directTai ), 2e-9 );
    BOOST_CHECK_THROW( getIgnNominalUtcStationTimeFromTai( leap + 0.5 ), std::runtime_error );
    BOOST_CHECK_THROW( getIgnNominalUtcStationTime( SinexCivilUtc{ 2016, 366, 86400 } ), std::runtime_error );
    BOOST_CHECK_THROW( convertSinexCivilUtcToTai( { 2016, 365, 86400 } ), std::runtime_error );
    BOOST_CHECK_THROW( getIgnNominalUtcStationTime( SinexCivilUtc{ 1971, 365, 0 } ), std::runtime_error );
    BOOST_CHECK_NO_THROW( getIgnNominalUtcStationTimeFromTai( midnight ) );
    BOOST_CHECK_THROW( getIgnNominalUtcStationTimeFromTai( midnight, 1.0 ), std::runtime_error );
    const auto before = getIgnNominalUtcStationTimeFromTai( leap - 1e-9L, 1.0 );
    BOOST_CHECK_EQUAL( before.calendarYear, 2016 );
    BOOST_CHECK( before.mjdFraction < 1.0L );
    const auto after = getIgnNominalUtcStationTimeFromTai( midnight + 1e-9L, 1.0 );
    BOOST_CHECK_EQUAL( after.calendarYear, 2017 );
    BOOST_CHECK( after.mjdFraction > 0.0L );
    const auto converter = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );
    const Time ordinary = midnight + 12345.25;
    const Time tdb = converter->getCurrentTime< Time >( basic_astrodynamics::tai_scale, basic_astrodynamics::tdb_scale, ordinary );
    const auto native = getIgnNominalUtcStationTimeFromTdb( tdb, converter, true );
    BOOST_CHECK_SMALL( static_cast< double >( native.tai - ordinary ), 1e-12 );
    BOOST_CHECK_SMALL( static_cast< double >( native.mjdFraction * 86400 - 12345.25L ), 1e-11 );
    BOOST_REQUIRE( native.mjdRate && native.decimalYearRate );
}

BOOST_AUTO_TEST_CASE( publishedHistoryPositionVelocityIdsAndStrictMetadata )
{
    auto solutions = product( "ab07_solutions.snx" );
    auto intervals = product( "ab07_discontinuities.snx" );
    auto deformation = product( "ab07_psd.snx" );
    const Eigen::Matrix3d basis = Eigen::Matrix3d::Identity( );
    auto trajectory = createIgnStationHistoryTrajectory( solutions, intervals, { "AB07", "A" }, &deformation, nullptr, basis );
    BOOST_CHECK_EQUAL( trajectory->getSolutions( ).size( ), 6 );
    BOOST_CHECK_EQUAL( trajectory->getPostSeismicTerms( ).size( ), 2 );
    BOOST_CHECK_EQUAL( trajectory->getMetadata( ).release, "ITRF2020" );
    const auto& last = trajectory->getSolutions( ).back( );
    BOOST_CHECK_EQUAL( last.source.key.solutionId, "6" );
    BOOST_CHECK_EQUAL( last.velocitySolutionId, "4" );
    BOOST_CHECK_EQUAL( *last.velocityValidity.startTai, convertSinexCivilUtcToTai( { 2020, 204, 22365 } ) );
    // There is no V identifier 6: the SSC's complete solution 6 velocity is used.
    for( int c = 0; c < 3; ++c )
        BOOST_CHECK_EQUAL( last.velocityMetresPerJulianYear[ c ], *last.source.velocityMetresPerSecond[ c ] * 31557600.0 );
    const auto event = trajectory->getPostSeismicTerms( ).front( ).eventTime;
    auto query = getIgnNominalUtcStationTimeFromTai( event.tai + 100.0, 1.0 );
    BOOST_CHECK( trajectory->getState( query ).allFinite( ) );
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( solutions, intervals, { "AB07", "A" }, &deformation ), std::runtime_error );
    auto wrongRelease = deformation;
    wrongRelease.release = "ITRF2020-u2024";
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( solutions, intervals, { "AB07", "A" }, &wrongRelease, nullptr, basis ),
                       std::runtime_error );
    auto orphan = deformation;
    orphan.parameters.erase( 2 );
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( solutions, intervals, { "AB07", "A" }, &orphan, nullptr, basis ),
                       std::runtime_error );
    auto wrongPartner = deformation;
    wrongPartner.parameters.at( 2 ).type = "TEXP_N";
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( solutions, intervals, { "AB07", "A" }, &wrongPartner, nullptr, basis ),
                       std::runtime_error );
    auto ambiguous = deformation;
    auto duplicate = ambiguous.sites.begin( )->second;
    duplicate.domesId = "00000M002";
    ambiguous.sites.emplace( duplicate.key, duplicate );
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( solutions, intervals, { "AB07", "A" }, &ambiguous, nullptr, basis ),
                       std::runtime_error );
    BOOST_CHECK_NO_THROW( createIgnStationHistoryTrajectory(
            solutions, intervals, { "AB07", "A" }, &ambiguous, nullptr, basis, trajectory->getMetadata( ).domesId ) );
    BOOST_CHECK_THROW(
            createIgnStationHistoryTrajectory( solutions, intervals, { "AB07", "A" }, &ambiguous, nullptr, basis, "absent DOMES" ),
            std::runtime_error );
    auto contradictoryVelocity = solutions;
    *contradictoryVelocity.parameters.at( 34 ).valueSi += 0.001 / 31557600.0;
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( contradictoryVelocity, intervals, { "AB07", "A" } ), std::runtime_error );
    auto velocityGap = intervals;
    velocityGap.discontinuities.back( ).start.civilUtc->secondsOfDay += 1;
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( solutions, velocityGap, { "AB07", "A" } ), std::runtime_error );
    auto missingTide = solutions;
    missingTide.tideSystem.clear( );
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( missingTide, intervals, { "AB07", "A" } ), std::runtime_error );
    auto unresolved = intervals;
    unresolved.discontinuities.front( ).start.meaning = SinexEpochMeaning::unspecified;
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( solutions, unresolved, { "AB07", "A" } ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( publishedRepeatedPsdPairsAndInputOrderIndependence )
{
    const auto deformation = product( "repeated_psd.snx" );
    const Eigen::Matrix3d basis = Eigen::AngleAxisd( 0.7, Eigen::Vector3d( 1, 2, 3 ).normalized( ) ).toRotationMatrix( );
    for( const std::string siteCode : { "AB48", "ANTC" } )
    {
        const SinexSiteKey key{ siteCode, "A" };
        const auto companions = controlledCompanions( deformation.sites.find( key )->second );
        const auto trajectory = createIgnStationHistoryTrajectory( companions.first, companions.second, key, &deformation, nullptr, basis );
        BOOST_CHECK_EQUAL( trajectory->getPostSeismicTerms( ).size( ), siteCode == "AB48" ? 4 : 5 );
        const auto time = getIgnNominalUtcStationTimeFromTai( convertSinexCivilUtcToTai( { 2019, 160, 43200 } ), 1.0 );
        Eigen::Vector3d expected = Eigen::Vector3d::Zero( );
        for( const auto& entry : deformation.parameters )
        {
            const auto& amplitude = entry.second;
            if( amplitude.key.site.siteCode != siteCode || amplitude.type[ 0 ] != 'A' ) continue;
            const double tau = deformation.parameters.at( entry.first + 1 ).originalValue;
            const auto event = getIgnNominalUtcStationTime( *amplitude.referenceEpoch.civilUtc );
            const double dt = static_cast< double >( ( time.mjdDay - event.mjdDay + time.mjdFraction - event.mjdFraction ) / 365.25L );
            const double correction = amplitude.type.substr( 1, 3 ) == "LOG" ? std::log( 1 + dt / tau ) : 1 - std::exp( -dt / tau );
            expected += basis.col( std::string( "ENU" ).find( amplitude.type.back( ) ) ) * ( amplitude.originalValue * correction );
        }
        closeVector( trajectory->getPosition( time ), expected, 2e-15 );

        std::ifstream stream( fixture( "repeated_psd.snx" ) );
        std::string line, prefix, suffix;
        std::vector< std::string > estimates;
        bool inEstimates = false, finishedEstimates = false;
        while( std::getline( stream, line ) )
        {
            if( line == "-SOLUTION/ESTIMATE" )
            {
                inEstimates = false;
                finishedEstimates = true;
            }
            if( inEstimates )
                estimates.push_back( line );
            else if( finishedEstimates )
                suffix += line + "\n";
            else
                prefix += line + "\n";
            if( line == "+SOLUTION/ESTIMATE" ) inEstimates = true;
        }
        std::reverse( estimates.begin( ), estimates.end( ) );
        for( const auto& estimate : estimates ) prefix += estimate + "\n";
        std::istringstream reorderedStream( prefix + suffix );
        auto reordered = readSinexStationHistory( reorderedStream, "reversed published estimates" );
        reordered.frame = deformation.frame;
        reordered.release = deformation.release;
        reordered.origin = deformation.origin;
        reordered.tideSystem = deformation.tideSystem;
        const auto alternate = createIgnStationHistoryTrajectory( companions.first, companions.second, key, &reordered, nullptr, basis );
        closeVector( alternate->getPosition( time ), trajectory->getPosition( time ), 1e-30 );
        closeVector( alternate->getState( time ).tail< 3 >( ), trajectory->getState( time ).tail< 3 >( ), 1e-30 );
    }
}

BOOST_AUTO_TEST_CASE( publishedSeasonalFactoryOriginIdentityAndFixedBasis )
{
    auto seasons = product( "albh_seasonal.snx" );
    const SinexSiteKey key{ "ALBH", "A" };
    const auto companions = controlledCompanions( seasons.sites.find( key )->second );
    const Eigen::Matrix3d basis = Eigen::AngleAxisd( -0.31, Eigen::Vector3d( 1, -2, 3 ).normalized( ) ).toRotationMatrix( );
    const auto trajectory = createIgnStationHistoryTrajectory( companions.first, companions.second, key, nullptr, &seasons, basis );
    BOOST_CHECK_EQUAL( trajectory->getSeasonalTerms( ).size( ), 12 );
    BOOST_CHECK( trajectory->hasFittedSeasonalTerms( ) );
    const auto time = getIgnNominalUtcStationTimeFromTai( convertSinexCivilUtcToTai( { 2020, 60, 43200 } ), 1.0 );
    Eigen::Vector3d expected = Eigen::Vector3d::Zero( );
    for( const auto& entry : seasons.parameters )
    {
        const auto& p = entry.second;
        if( p.key.solutionId != "4" ) continue;
        const double angle = 2 * std::acos( -1.0 ) * ( p.type[ 1 ] - '0' ) * 59.5 / 366;
        const double function = p.type.substr( 2, 3 ) == "COS" ? std::cos( angle ) : std::sin( angle );
        expected += basis.col( std::string( "ENU" ).find( p.type.back( ) ) ) * ( p.originalValue * function );
    }
    closeVector( trajectory->getPosition( time ), expected, 1e-16 );
    const auto withoutBasis =
            createIgnStationHistoryTrajectory( companions.first, companions.second, key, nullptr, &seasons, Eigen::Matrix3d::Identity( ) );
    closeVector( trajectory->getState( time ).tail< 3 >( ), basis * withoutBasis->getState( time ).tail< 3 >( ), 1e-23 );
    auto cf = seasons;
    cf.origin = "CF";
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( companions.first, companions.second, key, nullptr, &cf, basis ),
                       std::runtime_error );
    auto cfCoordinates = companions.first;
    auto cfIntervals = companions.second;
    cfCoordinates.origin = cfIntervals.origin = "CF";
    BOOST_CHECK_EQUAL( createIgnStationHistoryTrajectory( cfCoordinates, cfIntervals, key, nullptr, &cf, basis )->getMetadata( ).origin,
                       "CF" );
    auto missing = seasons;
    missing.parameters.erase( 1 );
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( companions.first, companions.second, key, nullptr, &missing, basis ),
                       std::runtime_error );
    BOOST_CHECK_THROW( createIgnStationHistoryTrajectory( companions.first, companions.second, key, nullptr, &seasons, 2 * basis ),
                       std::runtime_error );
}

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
    auto trajectory = std::make_shared< StationHistoryTrajectory >( std::vector< StationHistoryLinearSolution >{ a, b },
                                                                    std::vector< StationPostSeismicTerm >{},
                                                                    std::vector< StationSeasonalTerm >{},
                                                                    metadata( ) );
    auto time = manufacturedTime( 0 );
    time.tai = boundary - 1e-12L;
    closeVector( trajectory->getPosition( time ), a.positionMetres, 1e-30 );
    time.tai = boundary;
    closeVector( trajectory->getPosition( time ), b.positionMetres, 1e-30 );
    auto motion = std::make_shared< StationHistoryMotionModel >( trajectory, [ time ]( const auto&, bool rates ) {
        auto output = time;
        if( rates ) throw std::runtime_error( "manufactured position-only source callback" );
        return output;
    } );
    auto station = std::make_shared< GroundStationState >(
            Eigen::Vector3d( 100, 200, 300 ), coordinate_conversions::cartesian_position, nullptr, motion );
    closeVector( station->getCartesianPositionInTime( 0.0, "Earth" ), b.positionMetres, 1e-30 );
    BOOST_CHECK_THROW( station->getCartesianStateInTime( 0.0, "Earth" ), std::runtime_error );
    const auto settings = simulation_setup::stationHistoryGroundStationMotionSettings( trajectory, [ time ]( const auto&, bool rates ) {
        auto output = time;
        if( !rates ) output.mjdRate.reset( );
        return output;
    } );
    const auto factoryStation = simulation_setup::createGroundStationState( std::make_shared< simulation_setup::Body >( ),
                                                                            Eigen::Vector3d( 100, 200, 300 ),
                                                                            coordinate_conversions::cartesian_position,
                                                                            { settings } );
    closeVector( factoryStation->getCartesianPositionInTime( 0.0, "Earth" ), b.positionMetres, 1e-30 );
    closeVector( factoryStation->getCartesianStateInTime( 0.0, "Earth" ).head< 3 >( ), b.positionMetres, 1e-30 );
    const auto nominalFactoryStation = simulation_setup::createGroundStationState(
            std::make_shared< simulation_setup::Body >( ), Eigen::Vector3d( 100, 200, 300 ), coordinate_conversions::cartesian_position );
    closeVector( nominalFactoryStation->getCartesianPositionInTime( 0.0, "Earth" ), Eigen::Vector3d( 100, 200, 300 ), 1e-30 );
    LinearStationMotionModel legacyLinear( Eigen::Vector3d( 1, 2, 3 ), 0.0 );
    closeVector( legacyLinear.getBodyFixedStationMotion( 0.0 ).head< 3 >( ), Eigen::Vector3d::Zero( ), 1e-30 );
}
