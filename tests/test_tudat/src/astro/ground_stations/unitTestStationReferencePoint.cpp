#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/createStationReferencePoint.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"

using namespace tudat;
using namespace tudat::ground_stations;
using namespace tudat::input_output;
namespace
{
SinexStationHistory equipment( )
{
    return readSinexStationHistory( paths::getTudatTestDataPath( ) + "pr1092/reference_files/station_history/abmf_equipment.snx" );
}
SinexEquipmentHeaderTimeContext context( )
{
    return { std::nullopt, convertSinexCivilUtcToTai( { 2026, 288, 0 } ), "explicit controlled coverage bound, not a source header" };
}
std::shared_ptr< StationCoordinateHistory > marker( const Eigen::Vector3d& position = Eigen::Vector3d( 3e6, -5e6, 2e6 ),
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
    return std::make_shared< StationCoordinateHistory >( std::vector< StationHistoryLinearSolution >{ solution },
                                                         std::vector< StationPostSeismicTerm >{},
                                                         std::vector< StationSeasonalTerm >{},
                                                         meta );
}

//! Marker-to-ARP offsets of the ABMF equipment file, starting at the marker of the given history.
std::shared_ptr< StationReferencePointOffsetHistory > arp( std::shared_ptr< StationCoordinateHistory > source = marker( ) )
{
    auto geometry = createGeocentricStationReferencePointGeometry( Eigen::Vector3d( 3e6, -5e6, 2e6 ), source->getMetadata( ).frame );
    return createSinexAntennaReferencePointOffsets( source->getMetadata( ), equipment( ), "ABMF:ARP", geometry, context( ) );
}

//! Offset from the ARP to a calibrated point, for the first antenna of the given ARP offsets.
StationReferencePointOffset calibration( const std::shared_ptr< StationReferencePointOffsetHistory >& antennaOffsets,
                                         const std::string& frame = "ITRF2020" )
{
    const auto& antenna = antennaOffsets->getOffsets( ).front( );
    StationReferencePointOffset result;
    result.from = antenna.to;
    result.to = result.from;
    result.to.kind = StationReferencePointKind::calibratedInstrumentPoint;
    result.to.name = "ABMF:L1";
    result.validity = antenna.validity;
    result.frame = frame;
    result.equipmentType = antenna.equipmentType;
    result.equipmentSerialNumber = antenna.equipmentSerialNumber;
    result.calibrationId = "explicit controlled calibration";
    result.signal = "L1";
    result.source = "controlled reference-vector input";
    result.componentsMetres = Eigen::Vector3d( 0.01, -0.02, 0.03 );
    return result;
}

//! Station motion model of the history and offsets, evaluated at a fixed station history epoch.
std::shared_ptr< StationHistoryMotionModel > motion( const std::shared_ptr< StationCoordinateHistory >& history,
                                                     const std::vector< std::shared_ptr< StationReferencePointOffsetHistory > >& offsets,
                                                     const StationHistoryTime& time )
{
    const std::function< StationHistoryTime( const Time&, bool ) > timeFunction = [ time ]( const Time&, bool rates ) {
        auto output = time;
        if( !rates )
        {
            output.mjdRate.reset( );
            output.decimalYearRate.reset( );
        }
        return output;
    };
    return std::make_shared< StationHistoryMotionModel >( history, offsets, timeFunction );
}

//! Position of the point at the end of the offsets, from the station motion model (with a station at the origin).
Eigen::Vector3d pointPosition( const std::shared_ptr< StationCoordinateHistory >& history,
                               const std::vector< std::shared_ptr< StationReferencePointOffsetHistory > >& offsets,
                               const StationHistoryTime& time )
{
    const auto station = std::make_shared< GroundStationState >( Eigen::Vector3d::Zero( ) );
    return motion( history, offsets, time )->getBodyFixedStationPositionCorrection( Time( 0.0 ), station, "" );
}

//! State of the point at the end of the offsets, from the station motion model (with a station at the origin).
Eigen::Vector6d pointState( const std::shared_ptr< StationCoordinateHistory >& history,
                            const std::vector< std::shared_ptr< StationReferencePointOffsetHistory > >& offsets,
                            const StationHistoryTime& time )
{
    const auto station = std::make_shared< GroundStationState >( Eigen::Vector3d::Zero( ) );
    return motion( history, offsets, time )->getBodyFixedStationMotion( Time( 0.0 ), station, "" );
}
void close( const Eigen::Vector3d& actual, const Eigen::Vector3d& expected, double tolerance = 1e-9 )
{
    BOOST_CHECK_SMALL( ( actual - expected ).norm( ), tolerance );
}
}  // namespace

//! Apply published antenna offsets with their equipment identities and dated changes.
BOOST_AUTO_TEST_CASE( publishedUneDirectionEquipmentAssociationAndDatedSwitch )
{
    const auto source = marker( );
    const auto point = arp( source );
    const auto query = getIgnNominalUtcStationTimeFromTai( convertSinexCivilUtcToTai( { 2008, 200, 0 } ), 1.0 );
    const auto& record = point->getActiveOffset( query.tai );

    // Apply published antenna offsets with their equipment identities and dated changes.
    BOOST_CHECK_EQUAL( point->getOffsets( ).size( ), 3 );
    BOOST_CHECK_EQUAL( record.componentsMetres.x( ), 0.05 );
    BOOST_CHECK_EQUAL( record.equipmentType, "AERAT2775_43    SPKE" );
    BOOST_CHECK_EQUAL( record.equipmentSerialNumber, "5546" );
    const Eigen::Vector3d expected = source->getPosition( query ).normalized( ) * 0.05;
    close( point->getPositionCorrection( query ), expected, 2e-17 );
    close( pointPosition( source, { point }, query ), source->getPosition( query ) + expected );
    BOOST_CHECK_GT( expected.dot( source->getPosition( query ) ), 0.0 );
    const auto atSwitch = getIgnNominalUtcStationTimeFromTai( convertSinexCivilUtcToTai( { 2009, 288, 72000 } ), 1.0 );
    close( pointPosition( source, { point }, atSwitch ), source->getPosition( atSwitch ) );

    // The state at a switch uses the new (right-hand) record, consistent with the position.
    BOOST_REQUIRE_NO_THROW( pointState( source, { point }, atSwitch ) );
    close( pointState( source, { point }, atSwitch ).head< 3 >( ), pointPosition( source, { point }, atSwitch ) );
    auto xyz = equipment( );
    xyz.frame = source->getMetadata( ).frame;
    const auto geometry = createGeocentricStationReferencePointGeometry( source->getPosition( query ), xyz.frame );
    for( auto& eccentricity : xyz.eccentricities )
    {
        const auto& une = eccentricity.componentsMetres;
        const Eigen::Vector3d converted = geometry.enuToTerrestrial * Eigen::Vector3d( une[ 2 ], une[ 1 ], une[ 0 ] );
        for( int i = 0; i < 3; ++i )
        {
            eccentricity.componentsMetres[ i ] = converted[ i ];
        }
        eccentricity.basis = SinexEccentricityBasis::xyz;
    }
    const auto xyzPoint = createSinexAntennaReferencePointOffsets( source->getMetadata( ), xyz, "ABMF:XYZ_ARP", std::nullopt, context( ) );
    close( pointPosition( source, { xyzPoint }, query ), pointPosition( source, { point }, query ), 1e-30 );
}

//! Distinguish ellipsoid-normal antenna height from radial height.
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

    // Distinguish ellipsoid-normal antenna height from radial height.
    close( normal.enuToTerrestrial.col( 2 ),
           Eigen::Vector3d(
                   std::cos( latitude ) * std::cos( longitude ), std::cos( latitude ) * std::sin( longitude ), std::sin( latitude ) ),
           2e-14 );
    close( radial.enuToTerrestrial.col( 2 ), position.normalized( ), 1e-15 );
    BOOST_CHECK_GT( ( radial.enuToTerrestrial.col( 2 ) - normal.enuToTerrestrial.col( 2 ) ).norm( ), 0.05 );
    const auto source = marker( position );
    const auto point =
            createSinexAntennaReferencePointOffsets( source->getMetadata( ), equipment( ), "ABMF:geodetic_ARP", normal, context( ) );
    const auto query = getIgnNominalUtcStationTime( SinexCivilUtc{ 2008, 200, 0 } );
    close( point->getPositionCorrection( query ), 0.05 * normal.enuToTerrestrial.col( 2 ), 1e-16 );
}

//! Reject calibration for the wrong antenna or in the wrong reference-point order.
BOOST_AUTO_TEST_CASE( equipmentCalibrationMustFollowItsAntennaReferencePoint )
{
    const auto source = marker( );
    const auto point = arp( source );
    auto offset = calibration( point );
    offset.equipmentSerialNumber = "wrong antenna";
    const auto wrongAntenna =
            std::make_shared< StationReferencePointOffsetHistory >( std::vector< StationReferencePointOffset >{ offset } );

    // Reject calibration for the wrong antenna or in the wrong reference-point order.

    // Reject calibration associated with a different physical antenna.
    BOOST_CHECK_THROW( validateStationReferencePointOffsets( source->getMetadata( ), { point, wrongAntenna } ), std::runtime_error );

    // A calibration needs the ARP offsets before it, and the first offsets must start at the point of the history.
    const auto calibrated =
            std::make_shared< StationReferencePointOffsetHistory >( std::vector< StationReferencePointOffset >{ calibration( point ) } );
    BOOST_CHECK_NO_THROW( validateStationReferencePointOffsets( source->getMetadata( ), { point, calibrated } ) );
    BOOST_CHECK_THROW( validateStationReferencePointOffsets( source->getMetadata( ), { calibrated } ), std::runtime_error );
    BOOST_CHECK_THROW( validateStationReferencePointOffsets( source->getMetadata( ), { calibrated, point } ), std::runtime_error );
}

//! Combine calibrated offsets and their rates in the transformed marker frame.
BOOST_AUTO_TEST_CASE( calibratedReferenceVectorsRatesAndFrames )
{
    const auto source = marker( );
    const auto point = arp( source );
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
    const auto calibrated = std::make_shared< StationReferencePointOffsetHistory >( std::vector< StationReferencePointOffset >{ offset } );
    const auto query = getIgnNominalUtcStationTimeFromTai( reference + 10.0, 1.0 );

    // The marker and its ARP offsets are constant, so the velocity is that of the calibration.
    close( pointState( source, { point, calibrated }, query ).tail< 3 >( ), velocity, 1e-20 );
    close( pointPosition( source, { point, calibrated }, query ),
           pointPosition( source, { point }, query ) + offset.positionFunction( query ) );

    // Offsets must be given in the frame of the history: after a frame transformation of the marker history, they are
    // created from the transformed history.
    const auto catalogue = readOfficialItrf2020TransformationData( );
    const auto& meta = source->getMetadata( );
    const auto transformedMarker =
            transformStationCoordinateHistory( source, catalogue, { "ITRF2014", "ITRF2014" }, meta.origin, meta.tideSystem );

    // Offsets in the old terrestrial frame cannot follow a transformed marker.
    BOOST_CHECK_THROW( validateStationReferencePointOffsets( transformedMarker->getMetadata( ), { point } ), std::runtime_error );
    const auto pointAfterFrame = arp( transformedMarker );
    BOOST_CHECK_EQUAL( pointAfterFrame->getOffsets( ).front( ).frame, "ITRF2014" );
    BOOST_CHECK_EQUAL( pointAfterFrame->getGeometry( )->frame, "ITRF2014" );
    close( pointPosition( transformedMarker, { pointAfterFrame }, query ),
           transformedMarker->getPosition( query ) + pointAfterFrame->getPositionCorrection( query ) );
    const auto calibratedAfterFrame = std::make_shared< StationReferencePointOffsetHistory >(
            std::vector< StationReferencePointOffset >{ calibration( pointAfterFrame, "ITRF2014" ) } );
    BOOST_CHECK_NO_THROW(
            validateStationReferencePointOffsets( transformedMarker->getMetadata( ), { pointAfterFrame, calibratedAfterFrame } ) );
}

//! Preserve a precise equipment switch and the complete inertial offset velocity.
BOOST_AUTO_TEST_CASE( precisePublicPointSwitchAndCompleteInertialOffsetVelocity )
{
    const auto source = marker( );
    const auto point = arp( source );
    const auto converter = basic_astrodynamics::getGeocentricTimeConverter( );
    const std::function< StationHistoryTime( const Time&, bool ) > sourceClock = [ converter ]( const Time& tdb, bool rates ) {
        return getIgnNominalUtcStationTimeFromTdb( tdb, converter, rates );
    };
    const auto settings =
            simulation_setup::stationReferencePointSettings( source, { point }, Eigen::Vector3d( 3e6, -5e6, 2e6 ), sourceClock );
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

    // Preserve a precise equipment switch and the complete inertial offset velocity.
    BOOST_REQUIRE( left < right );
    BOOST_REQUIRE_EQUAL( static_cast< double >( left ), static_cast< double >( right ) );
    const auto leftSource = getIgnNominalUtcStationTimeFromTdb( left, converter );
    const auto rightSource = getIgnNominalUtcStationTimeFromTdb( right, converter );
    BOOST_REQUIRE( leftSource.tai < switchTai );
    BOOST_REQUIRE( rightSource.tai > switchTai );
    const Eigen::Vector3d leftPosition = station->getNominalStationState( )->getCartesianPositionInTime( left, "SSB" );
    const Eigen::Vector3d rightPosition = station->getNominalStationState( )->getCartesianPositionInTime( right, "SSB" );
    close( leftPosition - rightPosition, point->getPositionCorrection( leftSource ), 1e-9 );
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
    const Eigen::Vector3d localOffset = point->getPositionCorrection( sourceTime );
    const auto rotation = body->getRotationalEphemeris( );
    const Eigen::Vector3d expectedOffsetVelocity =
            rotation->getDerivativeOfRotationToBaseFrame( static_cast< double >( tdb ) ) * localOffset;
    Eigen::Vector6d markerState = Eigen::Vector6d::Zero( );
    markerState.head< 3 >( ) = source->getPosition( sourceTime );
    const auto inertialMarker = ephemerides::transformStateToInertialOrientation< double, Time >( markerState, tdb, rotation );
    close( actual.tail< 3 >( ) - inertialMarker.tail< 3 >( ), expectedOffsetVelocity, 1e-13 );
    BOOST_CHECK_GT( expectedOffsetVelocity.norm( ), 1e-6 );
}
