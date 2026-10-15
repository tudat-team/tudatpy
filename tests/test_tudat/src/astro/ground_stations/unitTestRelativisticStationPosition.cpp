#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/astro/ephemerides/customRotationalEphemeris.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"

using namespace tudat;
using namespace tudat::ground_stations;
namespace
{
Eigen::Vector6d state( const Eigen::Vector3d& position, const Eigen::Vector3d& velocity )
{
    return ( Eigen::Vector6d( ) << position, velocity ).finished( );
}
const auto identity = []( double ) { return Eigen::Quaterniond::Identity( ); };
const auto noRotationRate = []( double ) -> Eigen::Matrix3d { return Eigen::Matrix3d::Zero( ); };

std::shared_ptr< BodyCentricToBarycentricRelativisticStationMotion > constantModel(
        const RelativisticStationCoordinateCompatibility compatibility,
        const bool potential = true )
{
    std::vector< RelativisticStationExternalBody > bodies;
    if( potential )
        bodies = { { "Sun",
                     []( double ) { return state( Eigen::Vector3d( 149597870700.0, 0, 0 ), Eigen::Vector3d::Zero( ) ); },
                     []( ) { return 1.32712440018e20; } },
                   { "Moon",
                     []( double ) { return state( Eigen::Vector3d( 0, 384400000.0, 0 ), Eigen::Vector3d::Zero( ) ); },
                     []( ) { return 4.9028e12; } } };
    return std::make_shared< BodyCentricToBarycentricRelativisticStationMotion >(
            []( double t ) {
                return state( ( t - 1234 ) * Eigen::Vector3d( 27000, -14000, 6200 ), Eigen::Vector3d( 27000, -14000, 6200 ) );
            },
            identity,
            noRotationRate,
            bodies,
            compatibility );
}

struct PositionOnlyMotion : StationMotionModel {
    Eigen::Vector6d getBodyFixedStationMotion( const Time&, std::shared_ptr< GroundStationState >, const std::string& ) override
    {
        throw std::runtime_error( "Intrinsic rate unavailable" );
    }
    Eigen::Vector3d getBodyFixedStationPositionCorrection( const Time&, std::shared_ptr< GroundStationState >, const std::string& ) override
    {
        return Eigen::Vector3d( 17, -21, 8 );
    }
};

struct FactoryFixture {
    FactoryFixture( )
    {
        using namespace basic_astrodynamics;
        using namespace earth_orientation;
        input_output::EarthOrientationData data;
        data.source = "Controlled native clock for relativistic station factory";
        for( int i = 0; i < 7; ++i )
        {
            input_output::EarthOrientationRecord row;
            row.utcModifiedJulianDate = 57751 + i;
            const Time utc( ( row.utcModifiedJulianDate - 51544.5 ) * 86400 );
            row.ut1MinusUtc = -36.4 + static_cast< double >( sofa_interface::convertUTCtoTAI( utc ) - utc );
            row.polarMotion = { 1e-6, -2e-6 };
            row.celestialPoleOffsets = std::array< double, 2 >{ 0, 0 };
            data.records.emplace( row.utcModifiedJulianDate, row );
        }
        auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( data );
        clock = std::make_shared< TerrestrialTimeScaleConverter >( nullptr, getDefaultUT1CorrectionCalculator( ), nullptr, evaluator );
        auto polar = std::make_shared< PolarMotionCalculator >(
                nullptr, getDefaultPolarMotionCorrectionCalculator( ), evaluator, [ this ]( const Time& tai ) {
                    return clock->getUt1FromTai( tai );
                } );
        auto pn = std::make_shared< PrecessionNutationCalculator >( iau_2006, nullptr, nullptr, evaluator );
        angles = std::make_shared< EarthOrientationAnglesCalculator >( polar, pn, clock );
        time = static_cast< double >(
                clock->getCurrentTime< Time >( tai_scale, tdb_scale, Time( ( 57754 - 51544.5 ) * 86400 ) + 10037.0 ) );
        earth = std::make_shared< simulation_setup::Body >( );
        sun = std::make_shared< simulation_setup::Body >( );
        moon = std::make_shared< simulation_setup::Body >( );
        earth->setBodyName( "Earth" );
        sun->setBodyName( "Sun" );
        moon->setBodyName( "Moon" );
        earth->setIsBodyGlobalFrameOrigin( 0 );
        sun->setIsBodyGlobalFrameOrigin( 0 );
        moon->setIsBodyGlobalFrameOrigin( 0 );
        earth->setEphemeris( std::make_shared< ephemerides::CustomEphemeris<> >(
                [ this ]( double t ) {
                    return state( Eigen::Vector3d( 1.5e11, 2e10, -1e10 ) + ( t - time ) * hostVelocity, hostVelocity );
                },
                "SSB",
                "GCRS" ) );
        sun->setEphemeris( std::make_shared< ephemerides::CustomEphemeris<> >(
                []( double ) { return state( Eigen::Vector3d( -1e8, 2e8, -3e8 ), Eigen::Vector3d::Zero( ) ); }, "SSB", "GCRS" ) );
        moon->setEphemeris( std::make_shared< ephemerides::CustomEphemeris<> >(
                [ this ]( double t ) {
                    return state( Eigen::Vector3d( 1.5e11 + 3.8e8, 2e10, -1e10 ) + ( t - time ) * hostVelocity, hostVelocity );
                },
                "SSB",
                "GCRS" ) );
        sun->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( 1.32712440018e20 ) );
        moon->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( 4.9028e12 ) );
        bodies = simulation_setup::SystemOfBodies( "SSB", "GCRS", { { "Earth", earth }, { "Sun", sun }, { "Moon", moon } } );
    }
    const Eigen::Vector3d hostVelocity{ 28000, -11000, 4000 };
    double time;
    std::shared_ptr< simulation_setup::Body > earth, sun, moon;
    simulation_setup::SystemOfBodies bodies;
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > clock;
    std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > angles;
};
}  // namespace

BOOST_AUTO_TEST_CASE( iers_equation_fixtures_use_independent_70_digit_decimal_results )
{
    const Eigen::Vector3d point( 4075578.385, 931852.890, 4801570.154 );
    // Decimal fixtures generated directly from IERS 2010 Ch. 11 Eq. 11.18/19.
    // Inputs, precision, equations and PDF SHA256 are in the accompanying README.
    const Eigen::Vector3d unscaled( -0.05927009731058432782438507239103079827,
                                    0.00067500388626950807702023348174785247,
                                    -0.05176757493587349398402470505696181319 );
    const Eigen::Vector3d scaled( -0.11962235703801889715288507239103079827,
                                  -0.01312412407358704507197976651825214753,
                                  -0.12287051583384520679542470505696181319 );
    for( const auto compatibility : { unscaled_tcg_to_tcb_station_coordinates, tt_to_tdb_station_coordinates } )
    {
        auto model = constantModel( compatibility );
        BOOST_CHECK_EQUAL( model->getExternalBodies( ).size( ), 2 );
        BOOST_CHECK_EQUAL( model->getExternalBodies( )[ 0 ].name, "Sun" );
        const auto actual = model->correctCompleteLocalPosition( 1234, point, "SSB" );
        const Eigen::Vector3d expected = compatibility == tt_to_tdb_station_coordinates ? scaled : unscaled;
        BOOST_CHECK_SMALL( ( actual - expected ).norm( ), 8e-17 );
        const Eigen::Vector3d expectedRate( -7.23953384523679299906582658942596120526e-9,
                                            -1.65526948537310939009488779777092310547e-9,
                                            -8.52912798048462520908934078398968329889e-9 );
        BOOST_CHECK_SMALL(
                ( model->correctCompleteLocalState( 1234, state( point, Eigen::Vector3d::Zero( ) ), "SSB" ).tail< 3 >( ) - expectedRate )
                        .norm( ),
                1e-22 );
    }
}

BOOST_AUTO_TEST_CASE( stationary_comoving_point_matches_exact_flat_spacetime_lorentz_contraction )
{
    const long double c = 299792458.L, speed = 30000.L;
    const long double beta2 = speed * speed / ( c * c );
    // Stable exact expression for sqrt(1-beta^2)-1, at equal barycentric time.
    // The GCRS point is stationary, so its station-event time shift changes no position.
    const long double exact = -beta2 / ( 1 + std::sqrt( 1 - beta2 ) );
    BodyCentricToBarycentricRelativisticStationMotion model(
            []( double ) { return state( Eigen::Vector3d::Zero( ), Eigen::Vector3d( 30000, 0, 0 ) ); },
            identity,
            noRotationRate,
            {},
            unscaled_tcg_to_tcb_station_coordinates );
    const Eigen::Vector3d point( 6e6, 2e6, -3e6 );
    const auto correction = model.correctCompleteLocalPosition( 0, point, "SSB" );
    BOOST_CHECK_SMALL( correction( 0 ) - static_cast< double >( exact * 6e6L ), 1e-10 );
    BOOST_CHECK_SMALL( correction( 1 ), 1e-18 );
    BOOST_CHECK_SMALL( correction( 2 ), 1e-18 );
    BOOST_CHECK_LT( correction( 0 ), 0 );
}

BOOST_AUTO_TEST_CASE( terminal_spatial_map_consumes_full_intrinsic_state_independent_of_order )
{
    const Eigen::Vector3d nominal( 4e6, 1e6, 5e6 ), drift( .02, -.01, .03 );
    const double omega = 7e-5;
    const auto host = []( double t ) {
        return state( Eigen::Vector3d( 1.5e11, 0, 0 ) + t * Eigen::Vector3d( 30000, 1000, 2000 ) +
                              .5 * t * t * Eigen::Vector3d( .006, -.002, .001 ),
                      Eigen::Vector3d( 30000, 1000, 2000 ) + t * Eigen::Vector3d( .006, -.002, .001 ) );
    };
    const auto rotation = [ = ]( double t ) { return Eigen::Quaterniond( Eigen::AngleAxisd( omega * t, Eigen::Vector3d::UnitZ( ) ) ); };
    const auto rotationRate = [ = ]( double t ) -> Eigen::Matrix3d {
        Eigen::Matrix3d spin;
        spin << 0, -omega, 0, omega, 0, 0, 0, 0, 0;
        return spin * rotation( t ).toRotationMatrix( );
    };
    auto relativistic = std::make_shared< BodyCentricToBarycentricRelativisticStationMotion >(
            host,
            rotation,
            rotationRate,
            std::vector< RelativisticStationExternalBody >{
                    { "Sun", []( double ) { return Eigen::Vector6d::Zero( ); }, []( ) { return 1.32712440018e20; } } },
            tt_to_tdb_station_coordinates );
    auto intrinsic = std::make_shared< LinearStationMotionModel >( drift, 0 );
    auto station = std::make_shared< GroundStationState >(
            nominal,
            coordinate_conversions::cartesian_position,
            nullptr,
            std::make_shared< CombinedStationMotionModel >(
                    std::vector< std::shared_ptr< StationMotionModel > >{ relativistic, intrinsic } ) );
    auto reverse = std::make_shared< GroundStationState >(
            nominal,
            coordinate_conversions::cartesian_position,
            nullptr,
            std::make_shared< CombinedStationMotionModel >(
                    std::vector< std::shared_ptr< StationMotionModel > >{ intrinsic, relativistic } ) );
    const double t = 1000;
    const auto actual = station->getCartesianStateInTime( t, "SSB" );
    const auto intrinsicState = state( nominal + t * drift, drift );
    BOOST_CHECK_SMALL( ( actual - intrinsicState - relativistic->correctCompleteLocalState( t, intrinsicState, "SSB" ) ).norm( ), 8e-10 );
    BOOST_CHECK_SMALL( ( actual - reverse->getCartesianStateInTime( t, "SSB" ) ).norm( ), 1e-18 );
    BOOST_CHECK_GT( ( actual.tail< 3 >( ) - drift ).norm( ), 1e-6 );
    BOOST_CHECK_GT( ( relativistic->correctCompleteLocalPosition( t, intrinsicState.head< 3 >( ), "SSB" ) -
                      relativistic->correctCompleteLocalPosition( t, nominal, "SSB" ) )
                            .norm( ),
                    1e-7 );
    for( const double h : { 32., 8., 2. } )
    {
        const Eigen::Vector3d derivative =
                ( 8 * ( station->getCartesianPositionInTime( t + h, "SSB" ) - station->getCartesianPositionInTime( t - h, "SSB" ) ) -
                  ( station->getCartesianPositionInTime( t + 2 * h, "SSB" ) - station->getCartesianPositionInTime( t - 2 * h, "SSB" ) ) ) /
                ( 12 * h );
        BOOST_CHECK_SMALL( ( actual.tail< 3 >( ) - derivative ).norm( ), 3e-10 );
    }
    BOOST_CHECK_SMALL( ( station->getCartesianStateInTime( t, "Earth" ) - intrinsicState ).norm( ), 1e-18 );
    auto nested = std::make_shared< CombinedStationMotionModel >(
            std::vector< std::shared_ptr< StationMotionModel > >{ intrinsic, relativistic } );
    BOOST_CHECK_THROW( CombinedStationMotionModel( ( std::vector< std::shared_ptr< StationMotionModel > >{ nested, relativistic } ) ),
                       std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( position_only_access_preserves_unavailable_intrinsic_and_external_rates )
{
    auto relativistic = std::make_shared< BodyCentricToBarycentricRelativisticStationMotion >(
            []( double ) { return state( Eigen::Vector3d::Zero( ), Eigen::Vector3d( 30000, 0, 0 ) ); },
            identity,
            []( double ) -> Eigen::Matrix3d { throw std::runtime_error( "rotation rate unavailable" ); },
            std::vector< RelativisticStationExternalBody >{
                    { "Sun",
                      []( double ) { return state( Eigen::Vector3d( 1.5e11, 0, 0 ), Eigen::Vector3d::Constant( NAN ) ); },
                      []( ) { return 1.32712440018e20; } } },
            tt_to_tdb_station_coordinates );
    auto station = std::make_shared< GroundStationState >(
            Eigen::Vector3d( 4e6, 1e6, 5e6 ),
            coordinate_conversions::cartesian_position,
            nullptr,
            std::make_shared< CombinedStationMotionModel >(
                    std::vector< std::shared_ptr< StationMotionModel > >{ relativistic, std::make_shared< PositionOnlyMotion >( ) } ) );
    BOOST_CHECK_NO_THROW( station->getCartesianPositionInTime( 0, "SSB" ) );
    BOOST_CHECK_THROW( station->getCartesianStateInTime( 0, "SSB" ), std::runtime_error );
    BOOST_CHECK_NO_THROW( relativistic->correctCompleteLocalPosition( 0, Eigen::Vector3d::Ones( ), "SSB" ) );
    BOOST_CHECK_THROW( relativistic->correctCompleteLocalState( 0, Eigen::Vector6d::Ones( ), "SSB" ), std::runtime_error );
    BodyCentricToBarycentricRelativisticStationMotion missingExternalRate(
            []( double ) { return state( Eigen::Vector3d::Zero( ), Eigen::Vector3d( 30000, 0, 0 ) ); },
            identity,
            noRotationRate,
            relativistic->getExternalBodies( ),
            tt_to_tdb_station_coordinates );
    BOOST_CHECK_NO_THROW( missingExternalRate.correctCompleteLocalPosition( 0, Eigen::Vector3d::Ones( ), "SSB" ) );
    BOOST_CHECK_THROW( missingExternalRate.correctCompleteLocalState( 0, Eigen::Vector6d::Ones( ), "SSB" ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( invalid_external_potential_and_unavailable_fallback_rates_fail_honestly )
{
    const auto host = []( double ) { return state( Eigen::Vector3d::Zero( ), Eigen::Vector3d( 30000, 0, 0 ) ); };
    const RelativisticStationExternalBody sun{ "Sun", host, []( ) { return 1e20; } };
    BOOST_CHECK_THROW( BodyCentricToBarycentricRelativisticStationMotion(
                               host, identity, noRotationRate, {}, unspecified_station_coordinate_compatibility ),
                       std::invalid_argument );
    BOOST_CHECK_THROW( BodyCentricToBarycentricRelativisticStationMotion(
                               host, identity, noRotationRate, { sun, sun }, tt_to_tdb_station_coordinates ),
                       std::invalid_argument );
    BodyCentricToBarycentricRelativisticStationMotion coincident( host, identity, noRotationRate, { sun }, tt_to_tdb_station_coordinates );
    BOOST_CHECK_THROW( coincident.correctCompleteLocalPosition( 0, Eigen::Vector3d::Ones( ), "SSB" ), std::domain_error );
    double lastEpoch = NAN;
    BodyCentricToBarycentricRelativisticStationMotion covered(
            [ & ]( double t ) {
                lastEpoch = t;
                if( t != 0 ) throw std::out_of_range( "host coverage" );
                return host( t );
            },
            identity,
            noRotationRate,
            {},
            tt_to_tdb_station_coordinates );
    BOOST_CHECK_NO_THROW( covered.correctCompleteLocalPosition( 0, Eigen::Vector3d::Ones( ), "SSB" ) );
    BOOST_CHECK_THROW( covered.correctCompleteLocalState( 0, Eigen::Vector6d::Ones( ), "SSB" ), std::out_of_range );
    const Time centre( 300000000.0 );
    const double step = std::ldexp( 1.0, -27 );
    auto scalarEphemeris = std::make_shared< ephemerides::CustomEphemeris<> >(
            []( double time ) {
                const double dt = time - 300000000.0;
                return state( Eigen::Vector3d( 30000 * dt + 0.5e6 * dt * dt, 0, 0 ), Eigen::Vector3d( 30000 + 1e6 * dt, 0, 0 ) );
            },
            "SSB",
            "GCRS" );
    double restored = NAN;
    BodyCentricToBarycentricRelativisticStationMotion unresolved(
            [ & ]( const Time& time ) {
                restored = static_cast< double >( time );
                return scalarEphemeris->getCartesianStateFromExtendedTime( time );
            },
            identity,
            noRotationRate,
            {},
            unscaled_tcg_to_tcb_station_coordinates,
            step );
    const Eigen::Vector6d intrinsic = state( Eigen::Vector3d( 6e6, 2e6, -3e6 ), Eigen::Vector3d::Zero( ) );
    BOOST_CHECK_NO_THROW( unresolved.correctCompleteLocalPosition( centre, intrinsic.head< 3 >( ), "SSB" ) );
    auto preciseEphemeris = std::make_shared< ephemerides::CustomEphemeris< Time > >(
            [ & ]( const Time time ) {
                const double dt = static_cast< double >( time - centre );
                return state( Eigen::Vector3d( 30000 * dt + 0.5e6 * dt * dt, 0, 0 ), Eigen::Vector3d( 30000 + 1e6 * dt, 0, 0 ) );
            },
            "SSB",
            "GCRS" );
    BodyCentricToBarycentricRelativisticStationMotion precise(
            [ & ]( const Time& time ) { return preciseEphemeris->getCartesianStateFromExtendedTime( time ); },
            identity,
            noRotationRate,
            {},
            unscaled_tcg_to_tcb_station_coordinates,
            step );
    const Eigen::Vector3d expectedRate( -6e6 * 30000 * 1e6 * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT, 0, 0 );
    BOOST_CHECK_SMALL( ( precise.correctCompleteLocalState( centre, intrinsic, "SSB" ).tail< 3 >( ) - expectedRate ).norm( ), 1e-15 );
}

BOOST_AUTO_TEST_CASE( native_public_factory_uses_selected_external_bodies_and_tdb_rotation )
{
    using namespace basic_astrodynamics;
    FactoryFixture fixture;
    const Eigen::Vector3d nominal( 4e6, 1e6, 5e6 ), drift( .02, -.01, .03 );
    const auto settings = simulation_setup::relativisticStationMotionSettings( { "Sun", "Moon" } );
    auto physical = std::make_shared< GroundStationState >( nominal,
                                                            coordinate_conversions::cartesian_position,
                                                            nullptr,
                                                            std::make_shared< LinearStationMotionModel >( drift, fixture.time - 1000 ) );
    auto tdbRotation = std::make_shared< ephemerides::GcrsToItrsRotationModel >( fixture.angles, tdb_scale, "GCRS" );
    const Eigen::Vector6d intrinsic = physical->getCartesianStateInTime( fixture.time, "SSB" );
    const Eigen::Vector3d inertial = tdbRotation->getRotationToBaseFrame( fixture.time ) * intrinsic.head< 3 >( );
    const Eigen::Vector3d earthPosition = fixture.earth->getPositionInBaseFrameFromEphemeris( fixture.time );
    const double potential =
            1.32712440018e20 / ( earthPosition - fixture.sun->getPositionInBaseFrameFromEphemeris( fixture.time ) ).norm( ) +
            4.9028e12 / 3.8e8;
    const Eigen::Vector3d correction = -( potential * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT + 1.48082686741e-8 ) * inertial -
            .5 * fixture.hostVelocity.dot( inertial ) * fixture.hostVelocity * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
    std::string stationName;
    // Stations use the body's rotation model as it is; inside an environment it must take TDB input.
    for( const auto scale : { tdb_scale } )
    {
        fixture.earth->setRotationalEphemeris( std::make_shared< ephemerides::GcrsToItrsRotationModel >( fixture.angles, scale, "GCRS" ) );
        stationName = "relativistic_site_" + std::to_string( scale );
        auto stationSettings = std::make_shared< simulation_setup::GroundStationSettings >(
                stationName,
                nominal,
                coordinate_conversions::cartesian_position,
                std::vector< std::shared_ptr< simulation_setup::GroundStationMotionSettings > >{
                        settings, simulation_setup::linearGroundStationMotionSettings( drift, fixture.time - 1000 ) } );
        simulation_setup::createGroundStation( fixture.earth, stationSettings, fixture.bodies );
        auto station = fixture.earth->getGroundStation( stationName )->getNominalStationState( );
        BOOST_CHECK_SMALL( ( station->getCartesianPositionInTime( fixture.time, "SSB" ) - intrinsic.head< 3 >( ) -
                             tdbRotation->getRotationToTargetFrame( fixture.time ) * correction )
                                   .norm( ),
                           8e-10 );
        BOOST_CHECK_GT( ( station->getCartesianStateInTime( fixture.time, "SSB" ).tail< 3 >( ) - drift ).norm( ), 1e-6 );
        const Eigen::Vector3d velocity = station->getCartesianStateInTime( fixture.time, "SSB" ).tail< 3 >( );
        for( const double h : { 32., 8. } )
        {
            const double t = fixture.time;
            const Eigen::Vector3d derivative =
                    ( 8 * ( station->getCartesianPositionInTime( t + h, "SSB" ) - station->getCartesianPositionInTime( t - h, "SSB" ) ) -
                      ( station->getCartesianPositionInTime( t + 2 * h, "SSB" ) -
                        station->getCartesianPositionInTime( t - 2 * h, "SSB" ) ) ) /
                    ( 12 * h );
            BOOST_CHECK_SMALL( ( velocity - derivative ).norm( ), 3e-10 );
        }

        auto ephemeris = simulation_setup::createReferencePointEphemerisFromId( fixture.earth, stationName, fixture.bodies );
        BOOST_CHECK_EQUAL( ephemeris->getReferenceFrameOrigin( ), "SSB" );
        BOOST_CHECK_EQUAL( ephemeris->getReferenceFrameOrientation( ), "GCRS" );
        BOOST_CHECK_NO_THROW( ephemeris->getCartesianState( fixture.time ) );
    }
    const auto boundStation = fixture.earth->getGroundStation( stationName )->getNominalStationState( );
    const auto eop = fixture.clock->getEarthOrientationDataEvaluator( );
    for( const Time& boundary : { eop->getStartTai( ), eop->getEndTai( ) } )
    {
        const Time tdb = fixture.clock->getCurrentTime< Time >( tai_scale, tdb_scale, boundary );
        BOOST_CHECK_NO_THROW( boundStation->getCartesianPositionInTime( static_cast< double >( tdb ), "SSB" ) );
        BOOST_CHECK_NO_THROW( boundStation->getCartesianStateInTime( static_cast< double >( tdb ), "SSB" ) );
        BOOST_CHECK_NO_THROW( boundStation->getCartesianStateInTime( tdb, "SSB" ) );
    }
    const Time outside = fixture.clock->getCurrentTime< Time >( tai_scale, tdb_scale, eop->getEndTai( ) + 1e-3 );
    BOOST_CHECK_THROW( boundStation->getCartesianPositionInTime( outside, "SSB" ), std::runtime_error );
    BOOST_CHECK_THROW( boundStation->getCartesianStateInTime( outside, "SSB" ), std::runtime_error );

    const auto makeBadStation = [ & ]( std::vector< std::string > names ) {
        auto bad = std::make_shared< simulation_setup::GroundStationSettings >(
                "bad",
                nominal,
                coordinate_conversions::cartesian_position,
                std::vector< std::shared_ptr< simulation_setup::GroundStationMotionSettings > >{
                        simulation_setup::relativisticStationMotionSettings( names ) } );
        simulation_setup::createGroundStation( fixture.earth, bad, fixture.bodies );
    };
    BOOST_CHECK_THROW( makeBadStation( { "Earth" } ), std::invalid_argument );
    BOOST_CHECK_THROW( makeBadStation( { "Sun", "Sun" } ), std::invalid_argument );
    BOOST_CHECK_THROW( makeBadStation( { "Absent" } ), std::invalid_argument );
    auto mismatched = fixture.bodies;
    mismatched = simulation_setup::SystemOfBodies( "SSB", "ECLIPJ2000", fixture.bodies.getMap( ) );
    BOOST_CHECK_THROW( simulation_setup::createGroundStationState(
                               fixture.earth,
                               nominal,
                               coordinate_conversions::cartesian_position,
                               std::vector< std::shared_ptr< simulation_setup::GroundStationMotionSettings > >{ settings },
                               mismatched ),
                       std::invalid_argument );
    const auto equivalent = simulation_setup::SystemOfBodies( "SSB", "J2000", fixture.bodies.getMap( ) );
    BOOST_CHECK_NO_THROW( simulation_setup::createGroundStationState(
            fixture.earth,
            nominal,
            coordinate_conversions::cartesian_position,
            std::vector< std::shared_ptr< simulation_setup::GroundStationMotionSettings > >{ settings },
            equivalent ) );
    auto geocentricBodies = simulation_setup::SystemOfBodies( "Earth", "GCRS" );
    fixture.earth->setIsBodyGlobalFrameOrigin( 1 );
    auto ephemeris = simulation_setup::createReferencePointEphemerisFromId( fixture.earth, stationName, geocentricBodies );
    BOOST_CHECK_EQUAL( ephemeris->getReferenceFrameOrigin( ), "Earth" );
    BOOST_CHECK_EQUAL( ephemeris->getReferenceFrameOrientation( ), "GCRS" );
}

BOOST_AUTO_TEST_CASE( body_settings_factory_retains_environment_and_converts_external_body_origin_before_station_creation )
{
    using namespace simulation_setup;
    const Eigen::Vector3d nominal( 4e6, 1e6, 5e6 ), drift( .02, -.01, .03 );
    const Eigen::Vector3d earthPosition( 1.5e11, 2e10, -1e10 ), earthVelocity( 28000, -11000, 4000 );
    const Eigen::Vector3d sunPosition( -1e8, 2e8, -3e8 );
    const double time = 1000, gm = 1.32712440018e20;
    const Eigen::Vector3d hostPosition = earthPosition + time * earthVelocity;
    const Eigen::Vector3d separation = hostPosition - sunPosition;
    const double distance = separation.norm( ), potential = gm / distance;
    const double potentialRate = -gm * separation.dot( earthVelocity ) / ( distance * distance * distance );
    const Eigen::Vector3d intrinsic = nominal + time * drift;
    const double scale = potential * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT + 1.48082686741e-8;
    const Eigen::Vector3d expectedCorrection =
            -scale * intrinsic - .5 * earthVelocity.dot( intrinsic ) * earthVelocity * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
    const Eigen::Vector3d expectedVelocityCorrection = -potentialRate * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * intrinsic -
            scale * drift - .5 * earthVelocity.dot( drift ) * earthVelocity * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
    for( const bool convertSunOrigin : { false, true } )
    {
        auto earth = std::make_shared< BodySettings >( );
        auto sun = std::make_shared< BodySettings >( );
        earth->ephemerisSettings = customEphemerisSettings(
                [ = ]( double t ) { return state( earthPosition + t * earthVelocity, earthVelocity ); }, "SSB", "GCRS" );
        earth->rotationModelSettings = customRotationModelSettings(
                "GCRS", "ITRS", []( double ) -> Eigen::Matrix3d { return Eigen::Matrix3d::Identity( ); }, 1.0 );
        sun->ephemerisSettings = customEphemerisSettings(
                [ = ]( double t ) {
                    return convertSunOrigin ? state( sunPosition - earthPosition - t * earthVelocity, -earthVelocity )
                                            : state( sunPosition, Eigen::Vector3d::Zero( ) );
                },
                convertSunOrigin ? "Earth" : "SSB",
                "GCRS" );
        sun->gravityFieldSettings = centralGravitySettings( gm );
        earth->groundStationSettings.push_back( std::make_shared< GroundStationSettings >(
                "body_settings_site",
                nominal,
                coordinate_conversions::cartesian_position,
                std::vector< std::shared_ptr< GroundStationMotionSettings > >{ relativisticStationMotionSettings( ),
                                                                               linearGroundStationMotionSettings( drift, 0 ) } ) );
        const BodyListSettings settings( { { "Earth", earth }, { "Sun", sun } }, "SSB", "GCRS" );
        const auto bodies = createSystemOfBodies( settings );
        const auto station = bodies.at( "Earth" )->getGroundStation( "body_settings_site" )->getNominalStationState( );
        BOOST_CHECK_SMALL( ( bodies.at( "Sun" )->getPositionInBaseFrameFromEphemeris( time ) - sunPosition ).norm( ), 1e-6 );
        BOOST_CHECK_SMALL( ( station->getCartesianPositionInTime( time, "SSB" ) - intrinsic - expectedCorrection ).norm( ), 2e-9 );
        const Eigen::Vector3d velocityCorrection = station->getCartesianStateInTime( time, "SSB" ).tail< 3 >( ) - drift;
        BOOST_CHECK_SMALL( ( velocityCorrection - expectedVelocityCorrection ).norm( ), 2e-17 );
        BOOST_CHECK_GT( velocityCorrection.norm( ), 1e-9 );
        const auto ephemeris = createReferencePointEphemerisFromId( bodies.at( "Earth" ), "body_settings_site", bodies );
        BOOST_CHECK_EQUAL( ephemeris->getReferenceFrameOrientation( ), "GCRS" );
        BOOST_CHECK_EQUAL( ephemeris->getReferenceFrameOrigin( ), "SSB" );
    }
}
