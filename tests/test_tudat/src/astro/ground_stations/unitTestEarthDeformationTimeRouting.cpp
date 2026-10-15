#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/ground_stations/basicTidalBodyDeformation.h"
#include "tudat/astro/ground_stations/stationHistoryTrajectory.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/earth_orientation/precessionNutationCalculator.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/astro/ephemerides/customRotationalEphemeris.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/ground_stations/poleTideDeformation.h"
#include "tudat/astro/ground_stations/oceanPoleTideDeformation.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "iers2010_reference/fixtures.h"
#include "../../../reference_files/ocean_loading/hardisp_reference_cases.h"

using namespace tudat;
using namespace tudat::basic_astrodynamics;
namespace
{
class ContextOnlyPositionModel : public BodyDeformationModel
{
public:
    explicit ContextOnlyPositionModel( std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr ):
        BodyDeformationModel( converter )
    {}
    int stateCalls = 0;

protected:
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > ) override
    {
        return Eigen::Vector3d( static_cast< double >( tdb ), 42.0, 17.0 );
    }
    Eigen::Vector6d evaluateDisplacementState( const Time&, const std::shared_ptr< ground_stations::GroundStationState > ) override
    {
        ++stateCalls;
        throw std::runtime_error( "state deliberately unavailable" );
    }
};
}  // namespace

BOOST_AUTO_TEST_CASE( station_position_context_dispatch_never_invokes_a_state_or_derivative )
{
    const auto model = std::make_shared< ContextOnlyPositionModel >( );
    std::vector< std::shared_ptr< BodyDeformationModel > > models{ model };
    const auto motion =
            std::make_shared< ground_stations::BodyDeformationStationMotionModel >( [ &models ]( ) -> auto& { return models; }, true );
    const Eigen::Vector3d nominal( 6378137.0, 0.0, 0.0 );
    const auto station =
            std::make_shared< ground_stations::GroundStationState >( nominal, coordinate_conversions::cartesian_position, nullptr, motion );
    const double tdb = 12345.25;
    const auto position = station->getCartesianPositionInTime( tdb, "Earth" );
    BOOST_CHECK_SMALL( ( position - nominal - Eigen::Vector3d( tdb, 42.0, 17.0 ) ).norm( ), 1e-30 );
    BOOST_CHECK_EQUAL( model->stateCalls, 0 );
    BOOST_CHECK_THROW( motion->getBodyFixedStationMotion( tdb, station ), std::runtime_error );
    BOOST_CHECK_EQUAL( model->stateCalls, 1 );
}

namespace
{
using earth_orientation::EarthOrientationDataEvaluator;
using earth_orientation::TerrestrialTimeScaleConverter;
input_output::EarthOrientationData routingLeapData( )
{
    input_output::EarthOrientationData data;
    data.source = "controlled leap time routing";
    for( int i = 0; i < 4; ++i )
    {
        input_output::EarthOrientationRecord row;
        row.utcModifiedJulianDate = 57752.0 + i;
        row.polarMotion = { 1e-6 + i * 1e-9, -2e-6 };
        row.ut1MinusUtc = i < 2 ? -0.4 : 0.6;
        row.celestialPoleOffsets = std::array< double, 2 >{ 0.0, 0.0 };
        data.records.emplace( row.utcModifiedJulianDate, row );
    }
    return data;
}
std::shared_ptr< TerrestrialTimeScaleConverter > routingConverter( const input_output::EarthOrientationData& data )
{
    return std::make_shared< TerrestrialTimeScaleConverter >( std::make_shared< EarthOrientationDataEvaluator >( data ), nullptr );
}
Time referenceTdb( const Time& tt )
{
    // Direct SOFA call, independently of the Tudat time wrapper/context.
    const double dtr =
            iauDtdb( 2451545.0 + tt.fullDaysSinceEpoch( ), static_cast< double >( tt.secondsSinceNoon( ) / 86400.0L ), 0, 0, 0, 0 );
    return tt + dtr;
}
double referenceTdbRate( const Time& tt, double step )
{
    const auto offset = [ & ]( double shift ) {
        const Time t = tt + shift;
        return iauDtdb( 2451545.0 + t.fullDaysSinceEpoch( ), static_cast< double >( t.secondsSinceNoon( ) / 86400.0L ), 0, 0, 0, 0 );
    };
    return 1.0 + ( 8 * ( offset( step ) - offset( -step ) ) - ( offset( 2 * step ) - offset( -2 * step ) ) ) / ( 12 * step );
}
std::shared_ptr< earth_orientation::PolarMotionCalculator > routingPolar(
        const std::shared_ptr< TerrestrialTimeScaleConverter >& converter )
{
    return std::make_shared< earth_orientation::PolarMotionCalculator >(
            converter->getEarthOrientationDataEvaluator( ), nullptr, [ converter ]( const Time& tai ) {
                return converter->getUt1FromTai( tai );
            } );
}
input_output::OceanLoadingData routingOceanData( )
{
    input_output::OceanLoadingData data;
    for( const auto& fixture : hardispReferenceStations )
    {
        input_output::OceanLoadingStation row;
        row.stationIdentifier = fixture.identifier;
        for( int i = 0; i < 3; ++i )
        {
            for( int j = 0; j < 11; ++j )
            {
                row.amplitudesMetres( i, j ) = fixture.coefficients[ i ][ j ];
                row.phaseLagsDegrees( i, j ) = fixture.coefficients[ i + 3 ][ j ];
            }
        }
        data.stations.emplace( row.stationIdentifier, row );
    }
    return data;
}
input_output::OceanPoleTideGrid routingGrid( )
{
    input_output::OceanPoleTideGrid grid;
    grid.longitudeCount = 4;
    grid.latitudeCount = 3;
    grid.firstLongitudeDegrees = 45;
    grid.longitudeStepDegrees = 90;
    grid.firstLatitudeDegrees = -89.75;
    grid.latitudeStepDegrees = 89.75;
    for( int i = 0; i < 12; ++i )
    {
        grid.coefficients.push_back( Eigen::Matrix< double, 3, 2 >::Constant( 1.0 + i ) );
    }
    return grid;
}
}  // namespace

BOOST_AUTO_TEST_CASE( context_geocentric_clocks_match_independent_sofa_at_2020_seasonal_extrema )
{
    // Independent direct-SOFA scan/refinement of both extrema in calendar 2020.
    const std::array< std::array< double, 3 >, 2 > reference = { { { { 2459128.0, 0.39173499681965041, -0.0016902241652059767 } },
                                                                   { { 2458942.0, 0.78997176452867279, 0.0016780670772794989 } } } };
    for( const auto& row : reference )
    {
        const Time tt = Time( ( row[ 0 ] - 2451545.0 ) * 86400.0 ) + static_cast< long double >( row[ 1 ] ) * 86400.0L;
        const Time tdb = tt + row[ 2 ];
        const auto converter = getGeocentricTimeConverter( );
        const Time convertedTt = converter->getCurrentTime< Time >( tdb_scale, tt_scale, tdb );
        BOOST_CHECK_SMALL( static_cast< double >( convertedTt - tt ), 2e-12 );
        BOOST_CHECK_SMALL( static_cast< double >( tdb - convertedTt ) - row[ 2 ], 2e-12 );
        BOOST_CHECK_THROW( converter->getCurrentTime< Time >( tdb_scale, ut1_scale, tdb ), std::runtime_error );
    }

    // Atomic conversion must remain independent of EOP/DAT beyond leap-table warning epochs.
    const Time future = Time( ( 62867.0 - 51544.5 ) * 86400.0 );
    BOOST_CHECK_SMALL(
            static_cast< double >( getGeocentricTimeConverter( )->getCurrentTime< Time >( tdb_scale, tt_scale, referenceTdb( future ) ) -
                                   future ),
            2e-12 );
}

BOOST_AUTO_TEST_CASE( native_pole_context_accepts_rounded_double_endpoints_and_rejects_outside_coverage )
{
    const auto converter = routingConverter( routingLeapData( ) );
    const auto polar = routingPolar( converter );
    PoleTideDeformation solid( polar, converter );
    OceanPoleTideDeformation ocean( routingGrid( ), polar, converter );
    const auto station = std::make_shared< ground_stations::GroundStationState >( Eigen::Vector3d( 4075578.385, 931852.890, 4801570.154 ) );
    const auto evaluator = converter->getEarthOrientationDataEvaluator( );
    for( const Time& endpoint : { evaluator->getStartTai( ), evaluator->getEndTai( ) } )
    {
        const Time tdb = referenceTdb( convertTAItoTT( endpoint ) );
        const Time rounded( static_cast< double >( tdb ) );
        BOOST_CHECK_SMALL( static_cast< double >( converter->getTaiTime( tdb_scale, rounded ) - endpoint ), 1e-6 );
        BOOST_CHECK_NO_THROW( solid.calculateDisplacement( rounded, station ) );
        BOOST_CHECK_NO_THROW( ocean.calculateDisplacement( rounded, station ) );
    }
    const Time later = evaluator->getEndTai( ) + 1e-3;
    const Time precise = referenceTdb( convertTAItoTT( later ) );
    BOOST_CHECK_THROW( solid.calculateDisplacement( precise, station ), std::runtime_error );
    BOOST_CHECK_THROW( ocean.calculateDisplacement( precise, station ), std::runtime_error );
    for( double dt : { -1.5, -0.5, 0.5 } )
    {
        const Time tai = Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0 + dt;
        const Time epoch = referenceTdb( convertTAItoTT( tai ) );
        BOOST_CHECK_NO_THROW( solid.calculateDisplacement( epoch, station ) );
        BOOST_CHECK_NO_THROW( ocean.calculateDisplacement( epoch, station ) );
    }
}

BOOST_AUTO_TEST_CASE( explicit_custom_ocean_arguments_retain_environment_tdb_and_native_derivatives_reject_utc_steps )
{
    const auto station = std::make_shared< ground_stations::GroundStationState >( Eigen::Vector3d( 4075578.385, 931852.890, 4801570.154 ) );
    station->setSiteId( hardispReferenceStations[ 0 ].identifier );
    double last = NAN;
    OceanTideEarthDeformation custom( routingOceanData( ), [ & ]( double time ) {
        last = time;
        return Eigen::Vector6d::Zero( ).eval( );
    } );
    const Time epoch( 639169053.56045532 );
    BOOST_CHECK_NO_THROW( custom.calculateDisplacement( epoch, station ) );
    BOOST_CHECK_EQUAL( last, static_cast< double >( epoch ) );
    BOOST_CHECK( custom.calculateDisplacementState( epoch, station ).allFinite( ) );
    OceanTideEarthDeformation native( routingOceanData( ) );
    const Time midnight = Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0;
    const Time crossing = referenceTdb( convertTAItoTT( midnight ) );
    BOOST_CHECK_THROW( native.calculateDisplacementState( crossing, station ), std::domain_error );
    BOOST_CHECK_NO_THROW( native.calculateDisplacement( referenceTdb( convertTAItoTT( midnight + 5.0 ) ), station ) );
}

namespace
{
std::shared_ptr< ephemerides::GcrsToItrsRotationModel > routingRotation( const std::shared_ptr< TerrestrialTimeScaleConverter >& converter,
                                                                         TimeScales scale )
{
    const auto cip =
            std::make_shared< earth_orientation::PrecessionNutationCalculator >( iau_2006, converter->getEarthOrientationDataEvaluator( ) );
    const auto angles =
            std::make_shared< earth_orientation::EarthOrientationAnglesCalculator >( routingPolar( converter ), cip, converter );
    return std::make_shared< ephemerides::GcrsToItrsRotationModel >( angles, scale );
}
class ContextUt1PositionModel : public ContextOnlyPositionModel
{
public:
    explicit ContextUt1PositionModel( std::shared_ptr< TerrestrialTimeScaleConverter > converter ): ContextOnlyPositionModel( converter ) {}

protected:
    Eigen::Vector3d evaluateDisplacement( const Time& epoch, const std::shared_ptr< ground_stations::GroundStationState > ) override
    {
        lastUt1 = timeConverter_->getCurrentTime< Time >( tdb_scale, ut1_scale, epoch );
        return Eigen::Vector3d( 1, 2, 3 );
    }

public:
    Time lastUt1;
};
}  // namespace

BOOST_AUTO_TEST_CASE( public_station_factory_supplies_body_owned_clock_without_requesting_derivatives )
{
    const auto converter = routingConverter( routingLeapData( ) );
    const auto body = std::make_shared< simulation_setup::Body >( );
    body->setRotationalEphemeris( routingRotation( converter, tt_scale ) );
    const auto model = std::make_shared< ContextUt1PositionModel >( converter );
    body->addBodyDeformationModel( model );
    const Eigen::Vector3d nominal( 4075578.385, 931852.890, 4801570.154 );
    simulation_setup::createGroundStation( body,
                                           "clock_site",
                                           nominal,
                                           coordinate_conversions::cartesian_position,
                                           { simulation_setup::bodyDeformationStationMotionSettings( ) } );
    const auto station = body->getGroundStation( "clock_site" )->getNominalStationState( );
    const Time tai = converter->getEarthOrientationDataEvaluator( )->getStartTai( ) + 86400.0;
    const double tdb = static_cast< double >( referenceTdb( convertTAItoTT( tai ) ) );
    BOOST_CHECK_SMALL( ( station->getCartesianPositionInTime( tdb, "Earth" ) - nominal - Eigen::Vector3d( 1, 2, 3 ) ).norm( ), 1e-20 );
    BOOST_CHECK_EQUAL( model->lastUt1, converter->getCurrentTime< Time >( tdb_scale, ut1_scale, Time( tdb ) ) );
    BOOST_CHECK_EQUAL( model->stateCalls, 0 );
}

BOOST_AUTO_TEST_CASE( solid_tide_factory_uses_registered_rotation_and_complete_mixed_tdb_rate )
{
    const Eigen::Vector3d position( 4075578.385, 931852.890, 4801570.154 );
    const auto station = std::make_shared< ground_stations::GroundStationState >( position );
    const auto zero = std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ) );
    const auto identity = std::make_shared< ephemerides::CustomRotationalEphemeris >(
            []( double ) { return Eigen::Quaterniond::Identity( ); }, "GCRS", "ITRS" );
    for( const Time& tt : { Time( 639169053.56045532 ),
                            Time( 655205045.90372527 ),
                            convertTAItoTT( Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0 - 5.0 ) } )
    {
        auto data = routingLeapData( );
        if( tt > Time( 630000000.0 ) )
        {
            data.records.clear( );
            const double first = std::floor( static_cast< double >( tt ) / 86400.0 + 51544.5 ) - 1;
            for( int i = 0; i < 4; ++i )
            {
                input_output::EarthOrientationRecord row;
                row.utcModifiedJulianDate = first + i;
                row.polarMotion = { 1e-6, -2e-6 };
                row.ut1MinusUtc = -0.1;
                row.celestialPoleOffsets = std::array< double, 2 >{ 0.0, 0.0 };
                data.records.emplace( row.utcModifiedJulianDate, row );
            }
        }
        const auto converter = routingConverter( data );

        // A non-TDB model checks compatibility with the supplied rotation; environment epochs remain TDB.
        const auto rotation = routingRotation( converter, tt_scale );
        const Time tdb = referenceTdb( tt );
        double lastMoon = NAN, lastSun = NAN;
        const auto moonState = [ & ]( double time ) {
            Eigen::Vector6d state = Eigen::Vector6d::Zero( );
            state.head< 3 >( ) = Eigen::Vector3d( 2.8e8 + 800 * ( time - static_cast< double >( tdb ) ), -2.4e8, 1.1e8 );
            return state;
        };
        const auto sunState = [ & ]( double time ) {
            Eigen::Vector6d state = Eigen::Vector6d::Zero( );
            state.head< 3 >( ) = Eigen::Vector3d( 1.1e11, 1.0e11 + 30000 * ( time - static_cast< double >( tdb ) ), -1.0e10 );
            return state;
        };
        const auto moon = std::make_shared< ephemerides::CustomEphemeris<> >( [ & ]( double time ) {
            lastMoon = time;
            return moonState( time );
        } );
        const auto sun = std::make_shared< ephemerides::CustomEphemeris<> >( [ & ]( double time ) {
            lastSun = time;
            return sunState( time );
        } );
        const auto model = simulation_setup::createDefaultEarthIers2010DeformationModel(
                zero, moon, sun, rotation, []( ) { return 1.0; }, []( ) { return 0.0123000371; }, []( ) { return 332946.0482; } );

        // The factory passes the environment epoch to this registered model without changing its declared clock.
        const Eigen::Quaterniond q = rotation->getRotationToTargetFrameFromExtendedTime( tdb );
        Eigen::Vector6d moonItrs = moonState( static_cast< double >( tdb ) ), sunItrs = sunState( static_cast< double >( tdb ) );
        moonItrs.head< 3 >( ) = q * moonItrs.head< 3 >( ).eval( );
        sunItrs.head< 3 >( ) = q * sunItrs.head< 3 >( ).eval( );
        const auto direct = simulation_setup::createDefaultEarthIers2010DeformationModel(
                zero,
                std::make_shared< ephemerides::ConstantEphemeris >( moonItrs ),
                std::make_shared< ephemerides::ConstantEphemeris >( sunItrs ),
                identity,
                []( ) { return 1.0; },
                []( ) { return 0.0123000371; },
                []( ) { return 332946.0482; } );
        const Time epoch( tdb );  // Body ephemerides and the supplied rotation receive this same environment epoch.
        const Eigen::Vector3d expected = direct->calculateDisplacement( tdb, station );
        BOOST_CHECK_SMALL( ( model->calculateDisplacement( epoch, station ) - expected ).norm( ), 2e-12 );
        BOOST_CHECK_EQUAL( lastMoon, static_cast< double >( tdb ) );
        BOOST_CHECK_EQUAL( lastSun, static_cast< double >( tdb ) );
        const Eigen::Vector6d state = model->calculateDisplacementState( epoch, station );
        BOOST_CHECK_SMALL( ( state.head< 3 >( ) - expected ).norm( ), 2e-12 );
        BOOST_CHECK_EQUAL( lastMoon, static_cast< double >( tdb ) );  // Restore central callback state.
        const auto independentRate = [ & ]( double h ) -> Eigen::Vector3d {
            // Independent wide-step fourth-order derivative of the complete routed position.
            const Eigen::Vector3d a = model->calculateDisplacement( tdb + h, station );
            const Eigen::Vector3d b = model->calculateDisplacement( tdb - h, station );
            const Eigen::Vector3d c = model->calculateDisplacement( tdb + 2 * h, station );
            const Eigen::Vector3d d = model->calculateDisplacement( tdb - 2 * h, station );
            return ( 8 * ( a - b ) - ( c - d ) ) / ( 12 * h );
        };
        const Eigen::Vector3d rate1 = independentRate( 1.0 ), rate2 = independentRate( 2.0 );
        BOOST_CHECK_SMALL( ( rate1 - rate2 ).norm( ), 1e-10 );
        BOOST_CHECK_SMALL( ( state.tail< 3 >( ) - rate2 ).norm( ), 1e-10 );
    }
}

BOOST_AUTO_TEST_CASE( custom_solid_tide_callback_retains_original_tdb )
{
    const auto& fixture = iers2010_reference::fixtures.front( );
    const auto zero = std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ) );
    Eigen::Vector6d lunar = Eigen::Vector6d::Zero( ), solar = lunar;
    lunar.head< 3 >( ) = Eigen::Vector3d( fixture.moon.data( ) );
    solar.head< 3 >( ) = Eigen::Vector3d( fixture.sun.data( ) );
    double lastArgumentTime = NAN;
    const auto model = simulation_setup::createDefaultEarthIers2010DeformationModel(
            zero,
            std::make_shared< ephemerides::ConstantEphemeris >( lunar ),
            std::make_shared< ephemerides::ConstantEphemeris >( solar ),
            std::make_shared< ephemerides::CustomRotationalEphemeris >(
                    []( double ) { return Eigen::Quaterniond::Identity( ); }, "GCRS", "ITRS" ),
            []( ) { return 1.0; },
            []( ) { return 0.0123000371; },
            []( ) { return 332946.0482; },
            [ & ]( double time ) {
                lastArgumentTime = time;
                return Eigen::Vector6d::Zero( ).eval( );
            } );
    const auto station = std::make_shared< ground_stations::GroundStationState >( Eigen::Vector3d( fixture.station.data( ) ) );
    const Time epoch( 639169053.56045532 );
    BOOST_CHECK( model->calculateDisplacementState( epoch, station ).allFinite( ) );
    BOOST_CHECK_EQUAL( lastArgumentTime, static_cast< double >( epoch ) );
}

BOOST_AUTO_TEST_CASE( native_mixed_and_tt_only_models_keep_positions_when_interpolated_clock_rates_are_unavailable )
{
    const auto& fixture = iers2010_reference::fixtures.front( );
    const auto zero = std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ) );
    Eigen::Vector6d lunar = Eigen::Vector6d::Zero( ), solar = lunar;
    lunar.head< 3 >( ) = Eigen::Vector3d( fixture.moon.data( ) );
    solar.head< 3 >( ) = Eigen::Vector3d( fixture.sun.data( ) );
    const auto makeModel = [ & ]( const std::shared_ptr< ephemerides::RotationalEphemeris >& rotation ) {
        return simulation_setup::createDefaultEarthIers2010DeformationModel(
                zero,
                std::make_shared< ephemerides::ConstantEphemeris >( lunar ),
                std::make_shared< ephemerides::ConstantEphemeris >( solar ),
                rotation,
                []( ) { return 1.0; },
                []( ) { return 0.0123000371; },
                []( ) { return 332946.0482; } );
    };
    const auto settings = interpolators::interpolatorGenerationSettings< double >(
            interpolators::linearInterpolation( interpolators::huntingAlgorithm, interpolators::throw_exception_at_boundary ),
            -1000.0,
            1000.0,
            100.0 );
    auto data = routingLeapData( );
    input_output::EarthOrientationData shifted;
    for( auto row : data.records )
    {
        row.second.utcModifiedJulianDate -= 6210.;
        shifted.records.emplace( row.second.utcModifiedJulianDate, row.second );
    }
    const auto converter = std::make_shared< TerrestrialTimeScaleConverter >(
            std::make_shared< EarthOrientationDataEvaluator >( shifted ), nullptr, settings );
    const auto model = makeModel( routingRotation( converter, tt_scale ) );
    const Time epoch( 0.0 );
    BOOST_CHECK_NO_THROW( ground_stations::getIgnNominalUtcStationTimeFromTdb( 0.0, converter ) );
    BOOST_CHECK_THROW( ground_stations::getIgnNominalUtcStationTimeFromTdb( 0.0, converter, true ), std::runtime_error );
    const auto positionOnly = std::make_shared< ContextOnlyPositionModel >( converter );
    std::vector< std::shared_ptr< BodyDeformationModel > > selected{ positionOnly };
    const auto motion =
            std::make_shared< ground_stations::BodyDeformationStationMotionModel >( [ & ]( ) -> auto& { return selected; }, true );
    const auto dispatched = std::make_shared< ground_stations::GroundStationState >(
            Eigen::Vector3d( 6378137, 0, 0 ), coordinate_conversions::cartesian_position, nullptr, motion );
    BOOST_CHECK_NO_THROW( dispatched->getCartesianPositionInTime( 0.0, "Earth" ) );
    const auto station = std::make_shared< ground_stations::GroundStationState >( Eigen::Vector3d( fixture.station.data( ) ) );
    station->setSiteId( hardispReferenceStations[ 0 ].identifier );
    OceanTideEarthDeformation loading( routingOceanData( ), OceanTideSynthesis::iersHardisp, converter );
    BOOST_CHECK_NO_THROW( model->calculateDisplacement( epoch, station ) );
    BOOST_CHECK_THROW( model->calculateDisplacementState( epoch, station ), std::runtime_error );
    BOOST_CHECK_NO_THROW( loading.calculateDisplacement( epoch, station ) );
    BOOST_CHECK_THROW( loading.calculateDisplacementState( epoch, station ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( basic_solid_factory_uses_registered_rotation )
{
    // Near the independent 2020 positive maximum of the TT/TDB clock offset.
    const Time tt( 639169053.56045532 );
    const double tdb = static_cast< double >( referenceTdb( tt ) );
    input_output::EarthOrientationData data;
    for( int i = 0; i < 4; ++i )
    {
        input_output::EarthOrientationRecord row;
        row.utcModifiedJulianDate = 58941 + i;
        row.polarMotion = { 1e-6, -2e-6 };
        row.ut1MinusUtc = -0.1;
        row.celestialPoleOffsets = std::array< double, 2 >{ 0.0, 0.0 };
        data.records.emplace( row.utcModifiedJulianDate, row );
    }
    const auto converter = routingConverter( data );
    simulation_setup::SystemOfBodies bodies( "SSB", "GCRS" );
    bodies.createEmptyBody( "Earth" );
    bodies.createEmptyBody( "Moon" );
    Eigen::Vector6d lunar = Eigen::Vector6d::Zero( );
    lunar.head< 3 >( ) = Eigen::Vector3d( 2.8e8, -2.4e8, 1.1e8 );
    bodies.at( "Earth" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "GCRS" ) );
    bodies.at( "Moon" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( lunar, "SSB", "GCRS" ) );
    bodies.at( "Earth" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( 1.0 ) );
    bodies.at( "Moon" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( 0.0123000371 ) );
    bodies.processBodyFrameDefinitions( );
    const auto settings = std::make_shared< simulation_setup::BasicSolidBodyDeformationSettings >(
            std::vector< std::string >{ "Moon" }, std::map< int, std::pair< double, double > >{ { 2, { 0.6078, 0.0847 } } }, 6378136.6 );
    const auto station = std::make_shared< ground_stations::GroundStationState >( Eigen::Vector3d( 4075578.385, 931852.890, 4801570.154 ) );
    const Time epoch( tdb );

    // TDB is the environment convention; the TT case verifies the existing supplied-rotation behavior.
    for( const auto scale : { tt_scale, tdb_scale } )
    {
        const auto registeredRotation = routingRotation( converter, scale );
        bodies.at( "Earth" )->setRotationalEphemeris( registeredRotation );
        const auto model = simulation_setup::createBodyDeformationModel( settings, "Earth", bodies );

        // The degree-two equation provides a direct reference for the factory's rotated lunar geometry.
        const Eigen::Vector3d expected = calculateDegreeTwoBasicTidalDisplacement(
                0.0123000371,
                station->getNominalCartesianPosition( ).normalized( ),
                registeredRotation->getRotationToTargetFrameFromExtendedTime( epoch ) * lunar.head< 3 >( ),
                6378136.6,
                0.6078,
                0.0847 );
        BOOST_CHECK_SMALL( ( model->calculateDisplacement( epoch, station ) - expected ).norm( ), 2e-12 );
    }
}
