#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/simulation/environment_setup/createStationReferencePoint.h"
#include "tudat/astro/ground_stations/transformedStationHistoryTrajectory.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"

using namespace tudat;
using namespace tudat::ground_stations;
using namespace tudat::input_output;
namespace
{
SinexStationHistory equipment( )
{
    const std::string source = __FILE__;
    return readSinexStationHistory( source.substr( 0, source.find_last_of( "/\\" ) ) +
                                    "/../../../reference_files/station_history/abmf_equipment.snx" );
}
SinexEquipmentHeaderTimeContext context( )
{
    return { std::nullopt, convertSinexCivilUtcToTai( { 2026, 288, 0 } ), "explicit controlled coverage bound, not a source header" };
}
std::shared_ptr< StationHistoryTrajectory > marker( const Eigen::Vector3d& position = Eigen::Vector3d( 3e6, -5e6, 2e6 ),
                                                    const std::string& code = "abmf" )
{
    StationHistoryMetadata meta;
    meta.site = { code, "A" };
    meta.domesId = "97103M001";
    meta.frame = meta.release = "ITRF2020";
    meta.origin = "CM";
    meta.tideSystem = "explicit source tide convention";
    auto reference = getIgnNominalUtcStationTime( SinexCivilUtc{ 2008, 1, 0 } );
    meta.timeConvention = reference.convention;
    StationHistoryLinearSolution solution;
    solution.source.key = { meta.site, "1" };
    solution.referenceTime = reference;
    solution.positionMetres = position;
    return std::make_shared< StationHistoryTrajectory >( std::vector< StationHistoryLinearSolution >{ solution },
                                                         std::vector< StationPostSeismicTerm >{},
                                                         std::vector< StationSeasonalTerm >{},
                                                         meta );
}
std::shared_ptr< StationReferencePointTrajectory > arp( std::shared_ptr< StationHistoryTrajectory > source = marker( ) )
{
    auto geometry = createGeocentricStationReferencePointGeometry( Eigen::Vector3d( 3e6, -5e6, 2e6 ), source->getMetadata( ).frame );
    return createSinexAntennaReferencePointTrajectory( source, equipment( ), "ABMF:ARP", geometry, context( ) );
}
StationReferencePointOffset calibration( const std::shared_ptr< StationReferencePointTrajectory >& source )
{
    const auto& antenna = source->getOffsetHistory( )->getOffsets( ).front( );
    StationReferencePointOffset result;
    result.from = antenna.to;
    result.to = result.from;
    result.to.kind = StationReferencePointKind::calibratedInstrumentPoint;
    result.to.name = "ABMF:L1";
    result.validity = antenna.validity;
    result.frame = source->getMetadata( ).frame;
    result.equipmentType = antenna.equipmentType;
    result.equipmentSerialNumber = antenna.equipmentSerialNumber;
    result.calibrationId = "explicit controlled calibration";
    result.signal = "L1";
    result.source = "controlled reference-vector input";
    result.componentsMetres = Eigen::Vector3d( 0.01, -0.02, 0.03 );
    return result;
}
void close( const Eigen::Vector3d& actual, const Eigen::Vector3d& expected, double tolerance = 1e-9 )
{
    BOOST_CHECK_SMALL( ( actual - expected ).norm( ), tolerance );
}
}  // namespace

BOOST_AUTO_TEST_CASE( publishedUneDirectionEquipmentAssociationAndDatedSwitch )
{
    const auto source = marker( );
    const auto point = arp( source );
    const auto geometryForAlias = createGeocentricStationReferencePointGeometry( Eigen::Vector3d( 3e6, -5e6, 2e6 ), "ITRF2020" );
    const auto uppercaseMarker = marker( Eigen::Vector3d( 3e6, -5e6, 2e6 ), "ABMF" );
    BOOST_CHECK_THROW(
            createSinexAntennaReferencePointTrajectory( uppercaseMarker, equipment( ), "ABMF:ARP", geometryForAlias, context( ) ),
            std::runtime_error );
    const auto aliasedPoint = createSinexAntennaReferencePointTrajectory(
            uppercaseMarker, equipment( ), "ABMF:ARP", geometryForAlias, context( ), std::nullopt, SinexSiteKey{ "abmf", "A" } );
    BOOST_CHECK_EQUAL( aliasedPoint->getMetadata( ).site.siteCode, "ABMF" );
    BOOST_CHECK( aliasedPoint->getOffsetHistory( )->getOffsets( ).front( ).equipmentProvenance.text.find( "abmf" ) != std::string::npos );
    const auto query = getIgnNominalUtcStationTimeFromTai( convertSinexCivilUtcToTai( { 2008, 200, 0 } ), 1.0 );
    const auto& record = point->getOffsetHistory( )->getActiveOffset( query.tai );
    BOOST_CHECK_EQUAL( point->getOffsetHistory( )->getOffsets( ).size( ), 3 );
    BOOST_CHECK_EQUAL( record.componentsMetres.x( ), 0.05 );
    BOOST_CHECK_EQUAL( record.equipmentType, "AERAT2775_43    SPKE" );
    BOOST_CHECK_EQUAL( record.equipmentSerialNumber, "5546" );
    const Eigen::Vector3d expected = source->getPosition( query ).normalized( ) * 0.05;
    close( point->getOffsetHistory( )->getPositionCorrection( query ), expected, 2e-17 );
    close( point->getPosition( query ), source->getPosition( query ) + expected );
    BOOST_CHECK_GT( expected.dot( source->getPosition( query ) ), 0.0 );
    BOOST_CHECK( point->getMetadata( ).referencePointKind == StationReferencePointKind::antennaReferencePoint );
    BOOST_CHECK_EQUAL( point->getMetadata( ).referencePointName, "ABMF:ARP" );
    BOOST_CHECK_EQUAL( record.equipmentProvenance.source, equipment( ).source );
    BOOST_REQUIRE( record.headerTimeContext );
    BOOST_CHECK_EQUAL( record.headerTimeContext->source, context( ).source );
    const auto atSwitch = getIgnNominalUtcStationTimeFromTai( convertSinexCivilUtcToTai( { 2009, 288, 72000 } ), 1.0 );
    close( point->getPosition( atSwitch ), source->getPosition( atSwitch ) );

    // The state at a switch uses the new (right-hand) record, consistent with the position.
    BOOST_REQUIRE_NO_THROW( point->getState( atSwitch ) );
    close( point->getState( atSwitch ).head< 3 >( ), point->getPosition( atSwitch ) );
    auto outside = query;
    outside.tai = *context( ).endTai;
    BOOST_CHECK_THROW( point->getPosition( outside ), std::runtime_error );

    auto xyz = equipment( );
    xyz.frame = source->getMetadata( ).frame;
    const auto geometry = createGeocentricStationReferencePointGeometry( source->getPosition( query ), xyz.frame );
    for( auto& eccentricity : xyz.eccentricities )
    {
        const auto& une = eccentricity.componentsMetres;
        const Eigen::Vector3d converted = geometry.enuToTerrestrial * Eigen::Vector3d( une[ 2 ], une[ 1 ], une[ 0 ] );
        for( int i = 0; i < 3; ++i ) eccentricity.componentsMetres[ i ] = converted[ i ];
        eccentricity.basis = SinexEccentricityBasis::xyz;
    }
    const auto xyzPoint = createSinexAntennaReferencePointTrajectory( source, xyz, "ABMF:XYZ_ARP", std::nullopt, context( ) );
    close( xyzPoint->getPosition( query ), point->getPosition( query ), 1e-30 );
}

BOOST_AUTO_TEST_CASE( explicitOblateNormalDiffersFromRadialUp )
{
    const double radius = 6e6, flattening = 0.1, latitude = 0.7, longitude = -0.8;
    const double eccentricitySquared = flattening * ( 2 - flattening );
    const double n = radius / std::sqrt( 1 - eccentricitySquared * std::pow( std::sin( latitude ), 2 ) );
    const Eigen::Vector3d position( n * std::cos( latitude ) * std::cos( longitude ),
                                    n * std::cos( latitude ) * std::sin( longitude ),
                                    n * ( 1 - eccentricitySquared ) * std::sin( latitude ) );
    const auto normal = createGeodeticStationReferencePointGeometry( position, "ITRF2020", radius, flattening );
    const auto radial = createGeocentricStationReferencePointGeometry( position, "ITRF2020" );
    close( normal.enuToTerrestrial.col( 2 ),
           Eigen::Vector3d(
                   std::cos( latitude ) * std::cos( longitude ), std::cos( latitude ) * std::sin( longitude ), std::sin( latitude ) ),
           2e-14 );
    close( radial.enuToTerrestrial.col( 2 ), position.normalized( ), 1e-15 );
    BOOST_CHECK_GT( ( radial.enuToTerrestrial.col( 2 ) - normal.enuToTerrestrial.col( 2 ) ).norm( ), 0.05 );
    const auto source = marker( position );
    const auto point = createSinexAntennaReferencePointTrajectory( source, equipment( ), "ABMF:geodetic_ARP", normal, context( ) );
    const auto query = getIgnNominalUtcStationTime( SinexCivilUtc{ 2008, 200, 0 } );
    close( point->getOffsetHistory( )->getPositionCorrection( query ), 0.05 * normal.enuToTerrestrial.col( 2 ), 1e-16 );
}

BOOST_AUTO_TEST_CASE( actualHeaderCoverageAndEquipmentCalibrationGuards )
{
    const auto source = marker( );
    const auto geometry = createGeocentricStationReferencePointGeometry( Eigen::Vector3d( 3e6, -5e6, 2e6 ), "ITRF2020" );
    BOOST_CHECK_THROW( createSinexAntennaReferencePointTrajectory( source, equipment( ), "ARP", geometry ), std::runtime_error );
    auto gap = equipment( );
    gap.eccentricities[ 1 ].start.civilUtc->secondsOfDay += 1;
    BOOST_CHECK_THROW( createSinexAntennaReferencePointTrajectory( source, gap, "ARP", geometry, context( ) ), std::runtime_error );
    auto ambiguous = equipment( );
    ambiguous.eccentricities.push_back( ambiguous.eccentricities.front( ) );
    BOOST_CHECK_THROW( createSinexAntennaReferencePointTrajectory( source, ambiguous, "ARP", geometry, context( ) ), std::runtime_error );
    const auto point = arp( source );
    auto offset = calibration( point );
    offset.equipmentSerialNumber = "wrong antenna";
    BOOST_CHECK_THROW(
            createStationReferencePointTrajectory(
                    point, std::make_shared< StationReferencePointOffsetHistory >( std::vector< StationReferencePointOffset >{ offset } ) ),
            std::runtime_error );
    offset = calibration( point );
    offset.directionDependent = true;
    BOOST_CHECK_THROW( StationReferencePointOffsetHistory( { offset } ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( calibratedReferenceVectorsRatesAndFrameCompositionOrder )
{
    const auto point = arp( );
    auto offset = calibration( point );
    const Time reference = convertSinexCivilUtcToTai( { 2008, 200, 0 } );
    const Eigen::Vector3d velocity( 1e-6, -2e-6, 3e-6 );
    offset.positionFunction = [ reference, velocity ]( const auto& t ) -> Eigen::Vector3d {
        return Eigen::Vector3d( 0.01, -0.02, 0.03 ) + velocity * static_cast< double >( t.tai - reference );
    };
    offset.stateFunction = [ position = offset.positionFunction, velocity ]( const auto& t ) {
        Eigen::Vector6d result;
        result.head< 3 >( ) = position( t );
        result.tail< 3 >( ) = velocity * ( *t.mjdRate * 86400.0 );
        return result;
    };
    const auto calibrated = createStationReferencePointTrajectory(
            point, std::make_shared< StationReferencePointOffsetHistory >( std::vector< StationReferencePointOffset >{ offset } ) );
    const auto query = getIgnNominalUtcStationTimeFromTai( reference + 10.0, 1.0 );
    close( calibrated->getState( query ).tail< 3 >( ), velocity, 1e-20 );
    close( calibrated->getPosition( query ), point->getPosition( query ) + offset.positionFunction( query ) );
    auto noRates = query;
    noRates.mjdRate.reset( );
    noRates.decimalYearRate.reset( );
    BOOST_CHECK_NO_THROW( calibrated->getPosition( noRates ) );
    offset.stateFunction = nullptr;
    const auto positionOnly = createStationReferencePointTrajectory(
            point, std::make_shared< StationReferencePointOffsetHistory >( std::vector< StationReferencePointOffset >{ offset } ) );
    BOOST_CHECK_THROW( positionOnly->getState( query ), std::runtime_error );

    const auto catalogue = readOfficialItrf2020TransformationData( );
    const auto& meta = point->getMetadata( );
    const auto transformedPoint =
            createTransformedStationHistoryTrajectory( point, catalogue, { "ITRF2014", "ITRF2014" }, meta.origin, meta.tideSystem );

    // Both ordered wrappers expose the original record containers rather than copied histories.
    BOOST_CHECK( &transformedPoint->getSolutions( ) == &point->getSourceTrajectory( )->getSolutions( ) );
    BOOST_CHECK( &transformedPoint->getPostSeismicTerms( ) == &point->getSourceTrajectory( )->getPostSeismicTerms( ) );
    BOOST_CHECK( &transformedPoint->getSeasonalTerms( ) == &point->getSourceTrajectory( )->getSeasonalTerms( ) );
    BOOST_CHECK( getStationReferencePointIdentity( transformedPoint ).kind == StationReferencePointKind::antennaReferencePoint );
    BOOST_CHECK_THROW( arp( transformedPoint ), std::runtime_error );
    const auto transformedMarker =
            createTransformedStationHistoryTrajectory( marker( ), catalogue, { "ITRF2014", "ITRF2014" }, meta.origin, meta.tideSystem );
    const auto pointAfterFrame = arp( transformedMarker );
    BOOST_CHECK_EQUAL( pointAfterFrame->getMetadata( ).frame, "ITRF2014" );
    BOOST_CHECK_EQUAL( pointAfterFrame->getOffsetHistory( )->getGeometry( )->frame, "ITRF2014" );

    // Affine scale/rotation act on a pre-transform offset. A post-transform
    // vector must be supplied in the declared output geometry; do not silently
    // equate those distinct composition choices.
    reference_frames::TerrestrialFrameTransformationModel model( catalogue );

    // Preserve the source calendar coordinate and its TDB rate through ordered point/frame composition.
    const auto affine = model.getTransformation(
            "ITRF2020", "ITRF2014", static_cast< double >( query.calendarYear + query.yearFraction ), *query.decimalYearRate );

    close( transformedPoint->getPosition( query ), affine.transformPosition( point->getPosition( query ) ) );
    close( pointAfterFrame->getPosition( query ),
           transformedMarker->getPosition( query ) + pointAfterFrame->getOffsetHistory( )->getPositionCorrection( query ) );
    offset = calibration( point );
    offset.frame = "ITRF2014";
    BOOST_CHECK_NO_THROW( createStationReferencePointTrajectory(
            transformedPoint,
            std::make_shared< StationReferencePointOffsetHistory >( std::vector< StationReferencePointOffset >{ offset } ) ) );
}

BOOST_AUTO_TEST_CASE( precisePublicPointSwitchAndCompleteInertialOffsetVelocity )
{
    const auto point = arp( );
    const auto converter = basic_astrodynamics::getGeocentricTimeConverter( );
    const std::function< StationHistoryTime( const Time&, bool ) > sourceClock = [ converter ]( const Time& tdb, bool rates ) {
        return getIgnNominalUtcStationTimeFromTdb( tdb, converter, rates );
    };
    const auto settings = simulation_setup::stationReferencePointSettings( point, Eigen::Vector3d( 3e6, -5e6, 2e6 ), sourceClock );
    BOOST_CHECK_EQUAL( settings->getStationName( ), "ABMF:ARP" );
    const Time ordinaryTai = convertSinexCivilUtcToTai( { 2008, 200, 0 } );
    const Time tdb = converter->getCurrentTime< Time >( basic_astrodynamics::tai_scale, basic_astrodynamics::tdb_scale, ordinaryTai );
    const auto body = std::make_shared< simulation_setup::Body >( );
    body->setIsBodyGlobalFrameOrigin( 0 );
    body->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    body->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
            Eigen::Quaterniond::Identity( ), 1e-4, static_cast< double >( tdb ), "J2000", "ITRF2020" ) );
    const simulation_setup::SystemOfBodies bodies( "SSB", "J2000", { { "Earth", body } } );
    simulation_setup::createGroundStation( body, settings, bodies );
    const auto station = body->getGroundStation( "ABMF:ARP" );
    const Time switchTai = convertSinexCivilUtcToTai( { 2009, 288, 72000 } );
    const Time left =
            converter->getCurrentTime< Time >( basic_astrodynamics::tai_scale, basic_astrodynamics::tdb_scale, switchTai - 1e-10L );
    const Time right =
            converter->getCurrentTime< Time >( basic_astrodynamics::tai_scale, basic_astrodynamics::tdb_scale, switchTai + 1e-10L );
    BOOST_REQUIRE( left < right );
    BOOST_REQUIRE_EQUAL( static_cast< double >( left ), static_cast< double >( right ) );
    const auto leftSource = getIgnNominalUtcStationTimeFromTdb( left, converter );
    const auto rightSource = getIgnNominalUtcStationTimeFromTdb( right, converter );
    BOOST_REQUIRE( leftSource.tai < switchTai );
    BOOST_REQUIRE( rightSource.tai > switchTai );
    const Eigen::Vector3d leftPosition = station->getNominalStationState( )->getCartesianPositionInTime( left, "SSB" );
    const Eigen::Vector3d rightPosition = station->getNominalStationState( )->getCartesianPositionInTime( right, "SSB" );
    close( leftPosition - rightPosition, point->getOffsetHistory( )->getPositionCorrection( leftSource ), 1e-9 );
    BOOST_CHECK_GT( ( leftPosition - rightPosition ).norm( ), 0.049999999 );

    const auto ephemeris = simulation_setup::createReferencePointEphemerisFromId< Time, double >( body, "ABMF:ARP", bodies );
    BOOST_CHECK_EQUAL( ephemeris->getReferenceFrameOrigin( ), "SSB" );
    BOOST_CHECK_EQUAL( ephemeris->getReferenceFrameOrientation( ), "J2000" );
    const Eigen::Vector6d actual = ephemeris->getCartesianStateFromExtendedTime( tdb );
    const Eigen::Vector3d finiteDifference = ( ephemeris->getCartesianStateFromExtendedTime( tdb + 0.25 ).head< 3 >( ) -
                                               ephemeris->getCartesianStateFromExtendedTime( tdb - 0.25 ).head< 3 >( ) ) /
            0.5;
    close( actual.tail< 3 >( ), finiteDifference, 1e-7 );
    const auto sourceTime = getIgnNominalUtcStationTimeFromTdb( tdb, converter );
    const Eigen::Vector3d localOffset = point->getOffsetHistory( )->getPositionCorrection( sourceTime );
    const auto rotation = body->getRotationalEphemeris( );
    const Eigen::Vector3d expectedOffsetVelocity =
            rotation->getDerivativeOfRotationToBaseFrame( static_cast< double >( tdb ) ) * localOffset;
    Eigen::Vector6d markerState = Eigen::Vector6d::Zero( );
    markerState.head< 3 >( ) = point->getSourceTrajectory( )->getPosition( sourceTime );
    const auto inertialMarker = ephemerides::transformStateToInertialOrientation< double, Time >( markerState, tdb, rotation );
    close( actual.tail< 3 >( ) - inertialMarker.tail< 3 >( ), expectedOffsetVelocity, 1e-13 );
    BOOST_CHECK_GT( expectedOffsetVelocity.norm( ), 1e-6 );
}
