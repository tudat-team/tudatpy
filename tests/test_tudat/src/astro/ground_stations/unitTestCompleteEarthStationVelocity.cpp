#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/ground_stations/iers2010SolidTidalBodyDeformation.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"
#include "tudat/astro/ground_stations/poleTideDeformation.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"

using namespace tudat;
using namespace tudat::basic_astrodynamics;
using namespace tudat::earth_orientation;
using namespace tudat::ground_stations;

namespace
{
struct EarthStationFixture {
    EarthStationFixture( )
    {
        input_output::EarthOrientationData data;
        data.source = "Continuous controlled EOP for native station state closure";
        const Time reference = Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0;
        for( int i = 0; i < 7; ++i )
        {
            input_output::EarthOrientationRecord record;
            record.utcModifiedJulianDate = 57751.0 + i;
            const Time utc( ( record.utcModifiedJulianDate - 51544.5 ) * 86400.0 );
            const Time tai = sofa_interface::convertUTCtoTAI( utc );
            const double dt = static_cast< double >( tai - reference );
            record.polarMotion = { 1e-6 + 2e-12 * dt, 2e-6 - 3e-12 * dt };
            record.celestialPoleOffsets = std::array< double, 2 >{ 1e-9 + 4e-14 * dt, -2e-9 - 5e-14 * dt };
            record.ut1MinusUtc = -36.4 - 2e-8 * dt + static_cast< double >( tai - utc );
            data.records.emplace( record.utcModifiedJulianDate, record );
        }
        auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( data );
        clock = std::make_shared< TerrestrialTimeScaleConverter >( nullptr, getDefaultUT1CorrectionCalculator( ), nullptr, evaluator );
        auto polar = std::make_shared< PolarMotionCalculator >(
                nullptr, getDefaultPolarMotionCorrectionCalculator( ), evaluator, [ this ]( const Time& tai ) {
                    return clock->getUt1FromTai( tai );
                } );
        auto precession = std::make_shared< PrecessionNutationCalculator >( iau_2006, nullptr, nullptr, evaluator );
        angles = std::make_shared< EarthOrientationAnglesCalculator >( polar, precession, clock );
        rotation = std::make_shared< ephemerides::GcrsToItrsRotationModel >( angles, tdb_scale );
        referenceTdb = static_cast< double >( clock->getCurrentTime< Time >( tai_scale, tdb_scale, reference ) );
        const auto inertialState = [ this ]( const Eigen::Vector3d& position, const Eigen::Vector3d& velocity ) {
            return [ this, position, velocity ]( const double t ) -> Eigen::Vector6d {
                return ( Eigen::Vector6d( ) << position + ( t - referenceTdb ) * velocity, velocity ).finished( );
            };
        };
        moon = inertialState( Eigen::Vector3d( -179996231.920342, -312468450.131567, -169288918.592160 ),
                              Eigen::Vector3d( 850.0, -430.0, 110.0 ) );
        sun = inertialState( Eigen::Vector3d( 137859926952.015, 54228127881.435, 23509422341.696 ),
                             Eigen::Vector3d( -12000.0, 27000.0, 950.0 ) );
        zero = []( const double ) -> Eigen::Vector6d { return Eigen::Vector6d::Zero( ); };
        const auto rotate = [ this ]( const double t ) { return rotation->getRotationToTargetFrame( t ); };
        const auto contextRotate = [ this ]( const EarthDeformationEpoch& epoch ) {
            // Match native factory routing: retain the original double's exact
            // endpoint image before supplying the extended rotation epoch.
            const Time tdb = clock->getCurrentTime< Time >( tai_scale, tdb_scale, epoch.getTai( ) );
            return rotation->getRotationToTargetFrameFromExtendedTime( tdb );
        };
        const std::map< int, std::pair< double, double > > love = { { 2, { 0.6078, 0.0847 } }, { 3, { 0.292, 0.015 } } };
        solid = std::make_shared< Iers2010EarthDeformation >(
                zero,
                std::vector< std::function< Eigen::Vector6d( double ) > >{ moon, sun },
                rotate,
                []( ) { return 1.0; },
                std::vector< std::function< double( ) > >{ []( ) { return 0.0123000371; }, []( ) { return 332946.0482; } },
                6378136.6,
                love,
                std::vector< double >{ -0.0006, 0.0002 },
                std::vector< bool >( 6, true ),
                std::vector< double >{ 0.0012, 0.0024, -0.0025, -0.0007, -0.0022, -0.0007 },
                paths::getEarthDeformationDataFilesPath( ) + "/diurnalDisplacementFrequencyDependence2.txt",
                paths::getEarthDeformationDataFilesPath( ) + "/longPeriodDisplacementFrequencyDependence.txt",
                std::function< Eigen::Vector6d( double ) >( ),
                contextRotate,
                clock );
        basicSolid = std::make_shared< BasicTidalBodyDeformation >(
                zero,
                std::vector< std::function< Eigen::Vector6d( double ) > >{ moon, sun },
                rotate,
                []( ) { return 1.0; },
                std::vector< std::function< double( ) > >{ []( ) { return 0.0123000371; }, []( ) { return 332946.0482; } },
                6378136.6,
                love );
        input_output::OceanLoadingData oceanData;
        input_output::OceanLoadingStation oceanStation;
        oceanStation.stationIdentifier = "COMPLETE_VELOCITY";
        for( int row = 0; row < 3; ++row )
            for( int column = 0; column < 11; ++column )
            {
                oceanStation.amplitudesMetres( row, column ) = 0.001 + 0.0002 * row + 0.0001 * column;
                oceanStation.phaseLagsDegrees( row, column ) = 13.0 * column + 19.0 * row;
            }
        oceanData.stations.emplace( oceanStation.stationIdentifier, oceanStation );
        ocean = std::make_shared< OceanTideEarthDeformation >( oceanData );
        pole = std::make_shared< PoleTideDeformation >( polar, clock );
        models = { solid, ocean, pole };
        auto deformation = std::make_shared< BodyDeformationStationMotionModel >( [ this ]( ) -> auto& { return models; }, true, clock );
        auto secular = std::make_shared< LinearStationMotionModel >( secularVelocity, referenceTdb );
        auto combined = std::make_shared< CombinedStationMotionModel >(
                std::vector< std::shared_ptr< StationMotionModel > >{ secular, deformation } );
        auto shape = std::make_shared< OblateSpheroidBodyShapeModel >( 6378137.0, 1.0 / 298.257223563 );
        station = std::make_shared< GroundStationState >(
                Eigen::Vector3d( 4075578.385, 931852.890, 4801570.154 ), coordinate_conversions::cartesian_position, shape, combined );
        station->setSiteId( oceanStation.stationIdentifier );
    }
    Eigen::Vector3d inertialPosition( const double t, const std::shared_ptr< ephemerides::GcrsToItrsRotationModel >& model ) const
    {
        return model->getRotationToBaseFrame( t ) * station->getCartesianPositionInTime( t, "Earth" );
    }
    Eigen::Vector3d positionDerivative( const double t,
                                        const double h,
                                        const std::shared_ptr< ephemerides::GcrsToItrsRotationModel >& model ) const
    {
        return ( 8.0 * ( inertialPosition( t + h, model ) - inertialPosition( t - h, model ) ) -
                 ( inertialPosition( t + 2 * h, model ) - inertialPosition( t - 2 * h, model ) ) ) /
                ( 12.0 * h );
    }
    std::shared_ptr< TerrestrialTimeScaleConverter > clock;
    std::shared_ptr< EarthOrientationAnglesCalculator > angles;
    std::shared_ptr< ephemerides::GcrsToItrsRotationModel > rotation;
    std::shared_ptr< Iers2010EarthDeformation > solid;
    std::shared_ptr< BasicTidalBodyDeformation > basicSolid;
    std::shared_ptr< OceanTideEarthDeformation > ocean;
    std::shared_ptr< PoleTideDeformation > pole;
    std::vector< std::shared_ptr< BodyDeformationModel > > models;
    std::shared_ptr< GroundStationState > station;
    std::function< Eigen::Vector6d( double ) > moon, sun, zero;
    double referenceTdb;
    const Eigen::Vector3d secularVelocity = Eigen::Vector3d( 0.03, -0.02, 0.01 ) / physical_constants::JULIAN_YEAR;
};
}  // namespace

BOOST_AUTO_TEST_CASE( full_native_station_state_combines_each_deformation_velocity_once )
{
    EarthStationFixture fixture;
    for( const double offset : { 12345.6789, 40000.0 } )
    {
        const double t = fixture.referenceTdb + offset;
        Eigen::Vector6d deformation = Eigen::Vector6d::Zero( );
        const EarthDeformationEpoch epoch( t, fixture.clock );
        for( const auto& model : fixture.models )
        {
            const auto contribution = model->calculateDisplacementState( epoch, fixture.station );
            BOOST_CHECK_GT( contribution.head< 3 >( ).norm( ), 1e-5 );
            BOOST_CHECK_GT( contribution.tail< 3 >( ).norm( ), 1e-10 );
            deformation += contribution;
        }
        const auto local = fixture.station->getCartesianStateInTime( t, "Earth" );
        const Eigen::Vector3d expectedPosition =
                fixture.station->getNominalCartesianPosition( ) + offset * fixture.secularVelocity + deformation.head< 3 >( );
        const Eigen::Vector3d expectedVelocity = fixture.secularVelocity + deformation.tail< 3 >( );
        BOOST_CHECK_SMALL( ( local.head< 3 >( ) - expectedPosition ).norm( ), 2e-9 );
        BOOST_CHECK_SMALL( ( local.tail< 3 >( ) - expectedVelocity ).norm( ), 1e-16 );
        for( const std::string base : { "GCRS", "J2000", "EME2000", "ECLIPJ2000" } )
        {
            auto rotation = std::make_shared< ephemerides::GcrsToItrsRotationModel >( fixture.angles, tdb_scale, base );
            const auto inertial = ephemerides::transformStateToInertialOrientation< double, double >( local, t, rotation );
            const Eigen::Matrix3d r = rotation->getRotationToBaseFrame( t ).toRotationMatrix( );
            const Eigen::Matrix3d dr = rotation->getDerivativeOfRotationToBaseFrame( t );
            const Eigen::Vector3d independentComposition = dr * expectedPosition + r * expectedVelocity;
            BOOST_CHECK_SMALL( ( inertial.tail< 3 >( ) - independentComposition ).norm( ), 2e-13 );
            // These nonzero signals make omission or double application of either
            // displacement rotation or deformation velocity detectable by closure.
            BOOST_CHECK_GT( ( r * deformation.tail< 3 >( ) ).norm( ), 1e-7 );
            BOOST_CHECK_GT( ( dr * deformation.head< 3 >( ) ).norm( ), 1e-7 );
            for( const double step : { 32.0, 8.0 } )
            {
                const double error = ( inertial.tail< 3 >( ) - fixture.positionDerivative( t, step, rotation ) ).norm( );
                BOOST_TEST_MESSAGE( "Complete solid+ocean+pole station velocity " << base << " offset=" << offset << " h=" << step
                                                                                  << " s error=" << error << " m/s" );
                BOOST_CHECK_SMALL( error, 1e-7 );
            }
            const auto returned = ephemerides::transformStateToTargetFrame< double, double >( inertial, t, rotation );
            BOOST_CHECK_SMALL( ( returned.head< 3 >( ) - local.head< 3 >( ) ).norm( ), 4e-9 );
            // The inverse subtracts two ~300 m/s rotational terms to recover
            // a micrometre/s local velocity. Budget matrix/quaternion roundoff
            // against those terms rather than the small recovered result.
            const double inverseRoundoff =
                    32.0 * std::numeric_limits< double >::epsilon( ) * ( ( dr * expectedPosition ).norm( ) + expectedVelocity.norm( ) );
            BOOST_CHECK_SMALL( ( returned.tail< 3 >( ) - local.tail< 3 >( ) ).norm( ), inverseRoundoff );
        }
    }
}

BOOST_AUTO_TEST_CASE( continuous_native_solid_and_pole_states_close_through_inserted_second )
{
    EarthStationFixture fixture;
    fixture.models = { fixture.basicSolid, fixture.pole };
    const Time midnight = Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0;
    for( const double offset : { -1.5, -0.5, 0.0, 0.5 } )
    {
        const double t = static_cast< double >( fixture.clock->getCurrentTime< Time >( tai_scale, tdb_scale, midnight + offset ) );
        const auto local = fixture.station->getCartesianStateInTime( t, "Earth" );
        const auto inertial = ephemerides::transformStateToInertialOrientation< double, double >( local, t, fixture.rotation );
        BOOST_CHECK_SMALL( ( inertial.tail< 3 >( ) - fixture.positionDerivative( t, 32.0, fixture.rotation ) ).norm( ), 1e-7 );
        BOOST_CHECK_GT( local.tail< 3 >( ).norm( ), 1e-7 );
    }
}

BOOST_AUTO_TEST_CASE( native_station_state_preserves_utc_and_coverage_rate_unavailability )
{
    EarthStationFixture fixture;
    // Full continuous Earth rotation has a rate here; the native solid/ocean UTC
    // phase convention rejects a derivative interval touching this boundary.
    BOOST_CHECK_NO_THROW( fixture.rotation->getDerivativeOfRotationToBaseFrame( fixture.referenceTdb ) );
    BOOST_CHECK_NO_THROW( fixture.station->getCartesianPositionInTime( fixture.referenceTdb, "Earth" ) );
    BOOST_CHECK_THROW( fixture.station->getCartesianStateInTime( fixture.referenceTdb, "Earth" ), std::domain_error );
    const Time last = fixture.clock->getEarthOrientationDataEvaluator( )->getEndTai( );
    const double endpoint = static_cast< double >( fixture.clock->getCurrentTime< Time >( tai_scale, tdb_scale, last ) );
    BOOST_CHECK_NO_THROW( fixture.rotation->getDerivativeOfRotationToBaseFrame( endpoint ) );
    BOOST_CHECK_NO_THROW( fixture.station->getCartesianPositionInTime( endpoint, "Earth" ) );
    BOOST_CHECK_THROW( fixture.station->getCartesianStateInTime( endpoint, "Earth" ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( reference_point_ephemeris_routes_all_declared_earth_clocks_from_tdb )
{
    EarthStationFixture fixture;
    auto body = std::make_shared< simulation_setup::Body >( );
    body->setIsBodyGlobalFrameOrigin( 0 );
    body->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ) ) );
    // A constant local state isolates the factory rotation clock at the leap
    // and EOP endpoints, where native loading velocities are unavailable.
    Eigen::Vector6d local = Eigen::Vector6d::Zero( );
    local.head< 3 >( ) = fixture.station->getNominalCartesianPosition( );
    local.tail< 3 >( ) = fixture.secularVelocity;
    const auto evaluator = fixture.clock->getEarthOrientationDataEvaluator( );
    const Time leap = Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0;
    for( const auto scale : { tdb_scale, tt_scale, tai_scale, ut1_scale, utc_scale } )
    {
        body->setRotationalEphemeris( std::make_shared< ephemerides::GcrsToItrsRotationModel >( fixture.angles, scale, "ECLIPJ2000" ) );
        const std::function< Eigen::Vector6d( const double& ) > doubleState = [ local ]( const double& ) { return local; };
        const std::function< Eigen::Vector6d( const Time& ) > timeState = [ local ]( const Time& ) { return local; };
        const auto doubleEphemeris = simulation_setup::createReferencePointEphemeris< double, double >( body, doubleState );
        const auto timeEphemeris = simulation_setup::createReferencePointEphemeris< Time, double >( body, timeState );
        auto tdbRotation = std::make_shared< ephemerides::GcrsToItrsRotationModel >( fixture.angles, tdb_scale, "ECLIPJ2000" );
        for( const Time& tai : { evaluator->getStartTai( ), leap - 0.5, leap, leap + 0.5, evaluator->getEndTai( ) } )
        {
            const Time tdb = fixture.clock->getCurrentTime< Time >( tai_scale, tdb_scale, tai );
            const double t = static_cast< double >( tdb );
            const auto expectedDouble = ephemerides::transformStateToInertialOrientation< double, double >( local, t, tdbRotation );
            const auto expectedTime = ephemerides::transformStateToInertialOrientation< double, Time >( local, tdb, tdbRotation );
            const auto actualDouble = doubleEphemeris->getCartesianState( t );
            const auto actualTime = timeEphemeris->getCartesianStateFromExtendedTime( tdb );
            BOOST_CHECK_SMALL( ( actualDouble - expectedDouble ).norm( ), 1e-12 );
            BOOST_CHECK_SMALL( ( actualTime - expectedTime ).norm( ), 1e-12 );
        }
    }
}

BOOST_AUTO_TEST_CASE( public_station_ephemeris_closes_complete_native_velocity_in_tdb )
{
    EarthStationFixture fixture;
    auto body = std::make_shared< simulation_setup::Body >( );
    body->setIsBodyGlobalFrameOrigin( 0 );
    const Eigen::Vector3d bodyPosition( 1e6, -2e6, 3e6 ), bodyVelocity( 37.0, -19.0, 8.0 );
    const auto bodyState = [ & ]( const double t ) -> Eigen::Vector6d {
        return ( Eigen::Vector6d( ) << bodyPosition + ( t - fixture.referenceTdb ) * bodyVelocity, bodyVelocity ).finished( );
    };
    body->setEphemeris( std::make_shared< ephemerides::CustomEphemeris<> >( bodyState ) );
    const simulation_setup::SystemOfBodies bodies( "Earth", "ECLIPJ2000" );
    const double t = fixture.referenceTdb + 12345.6789;
    auto tdbRotation = std::make_shared< ephemerides::GcrsToItrsRotationModel >( fixture.angles, tdb_scale, "ECLIPJ2000" );
    const Eigen::Vector6d expected = bodyState( t ) +
            ephemerides::transformStateToInertialOrientation< double, double >(
                                             fixture.station->getCartesianStateInTime( t, "Earth" ), t, tdbRotation );
    for( const auto scale : { tdb_scale, tt_scale, tai_scale, ut1_scale, utc_scale } )
    {
        body->setRotationalEphemeris( std::make_shared< ephemerides::GcrsToItrsRotationModel >( fixture.angles, scale, "ECLIPJ2000" ) );
        simulation_setup::createGroundStation( body, "COMPLETE_VELOCITY", fixture.station );
        const auto ephemeris = simulation_setup::createReferencePointEphemerisFromId( body, "COMPLETE_VELOCITY", bodies );
        const Eigen::Vector6d actual = ephemeris->getCartesianState( t );
        BOOST_CHECK_SMALL( ( actual - expected ).norm( ), 1e-12 );
        for( const double h : { 32.0, 8.0 } )
        {
            const Eigen::Vector3d derivative =
                    ( 8.0 * ( ephemeris->getCartesianState( t + h ).head< 3 >( ) - ephemeris->getCartesianState( t - h ).head< 3 >( ) ) -
                      ( ephemeris->getCartesianState( t + 2 * h ).head< 3 >( ) -
                        ephemeris->getCartesianState( t - 2 * h ).head< 3 >( ) ) ) /
                    ( 12.0 * h );
            BOOST_CHECK_SMALL( ( actual.tail< 3 >( ) - derivative ).norm( ), 1e-7 );
        }
    }
}
