#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/ground_stations/stationLoadingDeformation.h"
#include "tudat/astro/ground_stations/stationHistoryMotionModel.h"
#include "tudat/astro/ground_stations/atmosphericTideDeformation.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

using namespace tudat;
using namespace tudat::input_output;
using namespace tudat::basic_astrodynamics;
using namespace tudat::ground_stations;
using namespace tudat::earth_orientation;
namespace
{
const Eigen::Vector3d providerPoint( 4075560.075, 931580.295, 4801621.958 );
const Time originTai( ( 57754 - 51544 ) * 24 - 12, 36.0L );
std::shared_ptr< TerrestrialTimeScaleConverter > loadingClock( )
{
    return std::make_shared< TerrestrialTimeScaleConverter >( nullptr, nullptr );
}
Time tdb( const std::shared_ptr< TerrestrialTimeScaleConverter >& converter, const Time& tai )
{
    return converter->getCurrentTime< Time >( tai_scale, tdb_scale, tai );
}
Eigen::Vector3d scalarFunction( const double seconds )
{
    return { .01 + 1e-5 * seconds + 1e-7 * seconds * seconds,
             -.02 + 2e-5 * seconds - 2e-8 * seconds * seconds,
             .03 - 3e-5 * seconds + 3e-8 * seconds * seconds };
}

// Independent scalar trigonometric axes, rather than the production vector construction.
Eigen::Matrix3d referenceBasis( )
{
    const double longitude = std::atan2( providerPoint.y( ), providerPoint.x( ) );
    const double latitude = std::asin( providerPoint.z( ) / providerPoint.norm( ) );
    Eigen::Matrix3d result;
    result << std::cos( latitude ) * std::cos( longitude ), -std::sin( longitude ), -std::sin( latitude ) * std::cos( longitude ),
            std::cos( latitude ) * std::sin( longitude ), std::cos( longitude ), -std::sin( latitude ) * std::sin( longitude ),
            std::sin( latitude ), 0.0, std::cos( latitude );
    return result;
}
//! Series with four samples, 100 s apart, starting 100 s before originTai.
StationLoadingData data( )
{
    StationLoadingData result;
    result.source = "analytic sample values";
    for( int i = 0; i < 4; ++i )
    {
        result.epochs.push_back( originTai + ( i - 1 ) * 100.0 );
    }
    StationLoadingSite site;
    site.cartesianPosition = providerPoint;
    site.validityRadiusMetres = 3000;
    for( std::size_t i = 0; i < 4; ++i )
    {
        site.displacements[ i ] = scalarFunction( ( static_cast< double >( i ) - 1 ) * 100.0 );
    }
    result.sites.emplace( "arbitrary_provider_label", site );
    return result;
}
const std::map< std::string, std::string > bindings{ { "ENVIRONMENT_STATION", "arbitrary_provider_label" } };
std::shared_ptr< GroundStationState > station( const Eigen::Vector3d& point = providerPoint,
                                               const std::shared_ptr< StationMotionModel >& motion = nullptr )
{
    auto result = std::make_shared< GroundStationState >( point,
                                                          coordinate_conversions::cartesian_position,
                                                          std::make_shared< OblateSpheroidBodyShapeModel >( 6378137, 1.0 / 298.257223563 ),
                                                          motion );
    result->setSiteId( "ENVIRONMENT_STATION" );
    return result;
}
}  // namespace

//! Compare loading positions and velocities with the selected linear interpolant.
BOOST_AUTO_TEST_CASE( linear_values_rates_knot_side_and_error_convergence_use_the_selected_interpolant )
{
    const auto converter = loadingClock( );
    const auto source = data( );
    StationLoadingDeformation model( source, bindings, 86400.0, converter );
    for( const double seconds : { -100.0, -75.0, 0.0, 25.0, 99.0, 100.0, 150.0, 200.0 } )
    {
        const Time input = tdb( converter, originTai + seconds );
        const double left = seconds == 200 ? 100 : 100 * std::floor( seconds / 100 );
        const Eigen::Vector3d slope = ( scalarFunction( left + 100 ) - scalarFunction( left ) ) / 100;
        const Eigen::Vector3d expected = referenceBasis( ) * ( scalarFunction( left ) + ( seconds - left ) * slope );
        const auto state = model.calculateDisplacementState( input, station( ) );

        // The position follows linear interpolation and the velocity follows its selected slope.
        BOOST_CHECK_SMALL( ( state.head< 3 >( ) - expected ).norm( ), 2e-12 );
        BOOST_CHECK_SMALL(
                ( state.tail< 3 >( ) - referenceBasis( ) * slope * converter->getGeocentricTtRateWithRespectToTdb( input ) ).norm( ),
                2e-16 );
        BOOST_CHECK_SMALL( ( model.calculateDisplacement( input, station( ) ) - state.head< 3 >( ) ).norm( ), 1e-18 );
    }

    // A separate scalar midpoint construction verifies O(h^2) position error for
    // this explicit linear policy; it makes no claim about provider B-splines.
    const auto midpointError = [ & ]( double h ) {
        auto grid = data( );
        grid.epochs = { originTai - h / 2, originTai + h / 2 };
        auto& samples = grid.sites.begin( )->second.displacements;
        samples.clear( );
        samples[ 0 ] = scalarFunction( -h / 2 );
        samples[ 1 ] = scalarFunction( h / 2 );
        StationLoadingDeformation interpolant( grid, bindings, 86400.0, converter );
        return ( interpolant.calculateDisplacement( tdb( converter, originTai ), station( ) ) - referenceBasis( ) * scalarFunction( 0 ) )
                .norm( );
    };

    // Halving the sample spacing reduces the linear interpolation position error by four.
    BOOST_CHECK_SMALL( midpointError( 100 ) / midpointError( 50 ) - 4.0, 1e-10 );
}

//! Preserve provider displacement geometry, with explicit and radius-based site selection.
BOOST_AUTO_TEST_CASE( provider_geometry_is_fixed_and_spatial_binding_is_explicit )
{
    const auto converter = loadingClock( );
    const auto source = data( );
    StationLoadingDeformation model( source, bindings, 86400.0, converter );
    const auto original = station( );
    const Time input = tdb( converter, originTai + 25.0 );
    const auto expected = model.calculateDisplacement( input, original );
    auto geodetic = std::make_shared< GroundStationState >(
            original->getNominalGeodeticPosition( ), coordinate_conversions::geodetic_position, original->getBodySurface( ) );
    geodetic->setSiteId( "ENVIRONMENT_STATION" );

    // Equivalent coordinates and nearby bound stations retain the same provider vector.
    BOOST_CHECK_SMALL( ( model.calculateDisplacement( input, geodetic ) - expected ).norm( ), 1e-15 );
    const auto nearby = station( providerPoint + Eigen::Vector3d( 10, -20, 30 ) );
    BOOST_CHECK_SMALL( ( model.calculateDisplacement( input, nearby ) - expected ).norm( ), 1e-15 );
    const Eigen::Vector3d local = original->getConstantRotationMatrixFromBodyFixedToTopocentricFrame( ) * expected;
    const Eigen::Vector3d raw = referenceBasis( ).transpose( ) * expected;

    // Oblate geodetic up differs from the provider radial direction.
    BOOST_CHECK_GT( std::abs( local.z( ) - raw.x( ) ), 1e-6 );
    StationLoadingDeformation matched( source, {}, 86400.0, converter );
    BOOST_CHECK_SMALL( ( matched.calculateDisplacement( input, original ) - expected ).norm( ), 1e-15 );
}

//! Close the full station velocity after factory-created non-tidal loading.
BOOST_AUTO_TEST_CASE( native_factory_and_full_station_state_preserve_the_loading_velocity_once )
{
    EarthOrientationData daily;
    daily.source = "controlled continuous daily EOP";
    for( int i = 0; i < 4; ++i )
    {
        EarthOrientationRecord row;
        row.utcModifiedJulianDate = 57753.0 + i;
        const Time utc( ( row.utcModifiedJulianDate - 51544.5 ) * 86400 );
        row.ut1MinusUtc = -37.4 + static_cast< double >( sofa_interface::convertUTCtoTAI( utc ) - utc );
        row.polarMotion = { 1e-6, 2e-6 };
        row.celestialPoleOffsets = std::array< double, 2 >{ 1e-9, -2e-9 };
        daily.records.emplace( row.utcModifiedJulianDate, row );
    }
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( daily );
    const auto converter = std::make_shared< TerrestrialTimeScaleConverter >( evaluator, nullptr, nullptr );
    const auto polar = std::make_shared< PolarMotionCalculator >( evaluator, nullptr );
    const auto celestial = std::make_shared< PrecessionNutationCalculator >( iau_2006, evaluator, nullptr );
    const auto rotation = std::make_shared< ephemerides::GcrsToItrsRotationModel >(
            std::make_shared< EarthOrientationAnglesCalculator >( polar, celestial, converter ), tdb_scale );
    simulation_setup::SystemOfBodies bodies( "SSB", "GCRS" );
    bodies.createEmptyBody( "Earth" );
    const auto earth = bodies.at( "Earth" );
    earth->setRotationalEphemeris( rotation );
    earth->setShapeModel( std::make_shared< OblateSpheroidBodyShapeModel >( 6378137, 1.0 / 298.257223563 ) );
    const auto source = data( );
    const auto configured = simulation_setup::stationLoadingBodyShapeDeformation( source, bindings );
    const auto loaded = simulation_setup::createBodyDeformationModel( configured, "Earth", bodies );
    StationLoadingDeformation direct( source, bindings, 86400.0, converter );
    earth->addBodyDeformationModel( loaded );
    const auto bound = simulation_setup::createGroundStationState( earth,
                                                                   providerPoint,
                                                                   coordinate_conversions::cartesian_position,
                                                                   { simulation_setup::bodyDeformationStationMotionSettings( true ) },
                                                                   bodies );
    bound->setSiteId( "ENVIRONMENT_STATION" );
    for( const double offset : { 50.0, 150.0 } )
    {
        const Time precise = tdb( converter, originTai + offset );
        const double scalar = static_cast< double >( precise );

        // Close the full station velocity after factory-created non-tidal loading.
        BOOST_CHECK_SMALL(
                ( loaded->calculateDisplacementState( precise, bound ) - direct.calculateDisplacementState( precise, bound ) ).norm( ),
                1e-18 );
        const Eigen::Vector6d local = bound->getCartesianStateInTime( scalar, "Earth" );
        const Eigen::Vector6d deformation = direct.calculateDisplacementState( Time( scalar ), bound );
        BOOST_CHECK_SMALL( ( local.head< 3 >( ) - providerPoint - deformation.head< 3 >( ) ).norm( ), 2e-9 );
        BOOST_CHECK_SMALL( ( local.tail< 3 >( ) - deformation.tail< 3 >( ) ).norm( ), 1e-15 );
        const auto inertial = ephemerides::transformStateToInertialOrientation< double, double >( local, scalar, rotation );
        const auto position = [ & ]( double time ) -> Eigen::Vector3d {
            return rotation->getRotationToBaseFrame( time ) * bound->getCartesianPositionInTime( time, "Earth" );
        };
        for( const double h : { 16.0, 8.0 } )
        {
            const Eigen::Vector3d numerical = ( 8.0 * ( position( scalar + h ) - position( scalar - h ) ) -
                                                ( position( scalar + 2 * h ) - position( scalar - 2 * h ) ) ) /
                    ( 12 * h );
            BOOST_CHECK_SMALL( ( inertial.tail< 3 >( ) - numerical ).norm( ), 2e-7 );
        }
    }
}
