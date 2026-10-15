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
StationLoadingData data( )
{
    StationLoadingData result;
    result.metadata = imlsMerra2AtmosphericLoadingMetadata( "controlled mathematical fixture" );
    result.metadata.source = "analytic sample values on split source epochs";
    result.metadata.referenceFrame = "ITRS";
    result.basis = StationLoadingBasis::geocentricRadialEastNorth;
    result.timeScale = StationLoadingTimeScale::tai;
    result.epochs = { { 57753, 86336.0L }, { 57754, 36.0L }, { 57754, 136.0L }, { 57754, 236.0L } };
    result.nominalSampleIntervalSeconds = 100.0L;
    StationLoadingSite site;
    site.identifier = "arbitrary_provider_label";
    site.cartesianPosition = { providerPoint.x( ), providerPoint.y( ), providerPoint.z( ) };
    site.validityRadiusMetres = 3000;
    for( std::size_t i = 0; i < 4; ++i )
    {
        const auto value = scalarFunction( ( static_cast< double >( i ) - 1 ) * 100.0 );
        site.samples[ i ].displacement = { value.x( ), value.y( ), value.z( ) };
    }
    result.sites.emplace( site.identifier, site );
    return result;
}
StationLoadingEvaluationSettings settings( const StationLoadingData& source )
{
    StationLoadingEvaluationSettings result;
    result.stationBindings = { { "ENVIRONMENT_STATION", source.sites.begin( )->first } };
    result.outputOriginDefinition = source.metadata.originDefinition;
    return result;
}
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
std::shared_ptr< StationLoadingOriginTranslation > translation( const StationLoadingData& source,
                                                                const std::string& identifier = "component-specific origin shift" )
{
    auto result = std::make_shared< StationLoadingOriginTranslation >( );
    result->identifier = identifier;
    result->sourceOriginDefinition = source.metadata.originDefinition;
    result->targetOriginDefinition = "controlled supplied centre-of-figure definition";
    result->targetOrigin = StationLoadingOrigin::centreOfFigure;
    result->series = source;
    result->series.basis = StationLoadingBasis::cartesian;
    result->series.metadata.originDefinition = result->targetOriginDefinition;
    result->series.metadata.referenceOrigin = result->targetOrigin;
    for( auto& sample : result->series.sites.begin( )->second.samples )
    {
        const double seconds = ( static_cast< double >( sample.first ) - 1 ) * 100;
        sample.second.displacement = { .001 + 2e-6 * seconds, -.002 + 3e-6 * seconds, .003 - 4e-6 * seconds };
    }
    return result;
}
}  // namespace

BOOST_AUTO_TEST_CASE( linear_values_rates_knot_side_and_error_convergence_use_the_selected_interpolant )
{
    const auto converter = loadingClock( );
    const auto source = data( );
    StationLoadingDeformation model( source, settings( source ), converter );
    for( const double seconds : { -100.0, -75.0, 0.0, 25.0, 99.0, 100.0, 150.0, 200.0 } )
    {
        const Time input = tdb( converter, originTai + seconds );
        const double left = seconds == 200 ? 100 : 100 * std::floor( seconds / 100 );
        const Eigen::Vector3d slope = ( scalarFunction( left + 100 ) - scalarFunction( left ) ) / 100;
        const Eigen::Vector3d expected = referenceBasis( ) * ( scalarFunction( left ) + ( seconds - left ) * slope );
        const auto state = model.calculateDisplacementState( input, station( ) );
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
        const long double begin = 36.0L - h / 2;
        grid.epochs = { { begin < 0 ? 57753 : 57754, begin < 0 ? begin + 86400.0L : begin }, { 57754, 36.0L + h / 2 } };
        grid.nominalSampleIntervalSeconds = h;
        auto& samples = grid.sites.begin( )->second.samples;
        samples.clear( );
        const auto left = scalarFunction( -h / 2 ), right = scalarFunction( h / 2 );
        samples[ 0 ].displacement = { left.x( ), left.y( ), left.z( ) };
        samples[ 1 ].displacement = { right.x( ), right.y( ), right.z( ) };
        StationLoadingDeformation interpolant( grid, settings( grid ), converter );
        return ( interpolant.calculateDisplacement( tdb( converter, originTai ), station( ) ) - referenceBasis( ) * scalarFunction( 0 ) )
                .norm( );
    };
    BOOST_CHECK_SMALL( midpointError( 100 ) / midpointError( 50 ) - 4.0, 1e-10 );
}

BOOST_AUTO_TEST_CASE( strict_precise_coverage_singletons_gaps_and_optional_rate_failure )
{
    const auto converter = loadingClock( );
    auto source = data( );
    const auto options = settings( source );
    StationLoadingDeformation model( source, options, converter );
    const Time first = tdb( converter, originTai - 100.0 ), last = tdb( converter, originTai + 200.0 );
    BOOST_CHECK_NO_THROW( model.calculateDisplacement( first, station( ) ) );
    BOOST_CHECK_NO_THROW( model.calculateDisplacementState( last, station( ) ) );
    BOOST_CHECK_THROW( model.calculateDisplacement( first - 1e-3L, station( ) ), std::out_of_range );
    BOOST_CHECK_THROW( model.calculateDisplacement( last + 1e-3L, station( ) ), std::out_of_range );
    auto gapPolicy = options;
    gapPolicy.maximumGapSeconds = 50;
    StationLoadingDeformation gaps( source, gapPolicy, converter );
    BOOST_CHECK_NO_THROW( gaps.calculateDisplacement( first, station( ) ) );
    BOOST_CHECK_THROW( gaps.calculateDisplacementState( first, station( ) ), std::out_of_range );
    BOOST_CHECK_THROW( gaps.calculateDisplacement( tdb( converter, originTai - 50.0 ), station( ) ), std::out_of_range );
    source.epochs.resize( 1 );
    source.sites.begin( )->second.samples.erase( 1 );
    source.sites.begin( )->second.samples.erase( 2 );
    source.sites.begin( )->second.samples.erase( 3 );
    auto singleton = std::make_shared< StationLoadingDeformation >( source, options, converter );
    BOOST_CHECK_NO_THROW( singleton->calculateDisplacement( first, station( ) ) );
    BOOST_CHECK_THROW( singleton->calculateDisplacementState( first, station( ) ), std::out_of_range );
    std::vector< std::shared_ptr< BodyDeformationModel > > models{ singleton };
    BodyDeformationStationMotionModel optional( [ & ]( ) -> auto& { return models; }, false );
    // The genuine scalar adapter recognizes exactly the rounded endpoint image.
    BOOST_CHECK_GT( optional.getBodyFixedStationPositionCorrection( static_cast< double >( first ), station( ), "Earth" ).norm( ), 0.0 );
    BOOST_CHECK_THROW( optional.getBodyFixedStationMotion( static_cast< double >( first ), station( ), "Earth" ), std::out_of_range );
    BOOST_CHECK_EQUAL( optional.getBodyFixedStationPositionCorrection( static_cast< double >( first + 1.0 ), station( ), "Earth" ).norm( ),
                       0.0 );
    auto fine = data( );
    fine.epochs = { { 57754, 36.0L }, { 57754, 36.0L + 1e-9L } };
    fine.nominalSampleIntervalSeconds = 1e-9L;
    auto& fineSamples = fine.sites.begin( )->second.samples;
    fineSamples.clear( );
    fineSamples[ 0 ].displacement = { .01, 0, 0 };
    fineSamples[ 1 ].displacement = { .02, 0, 0 };
    auto preciseOnly = std::make_shared< StationLoadingDeformation >( fine, settings( fine ), converter );
    BOOST_CHECK_SMALL( preciseOnly->calculateDisplacement( tdb( converter, originTai + 1e-9L ), station( ) ).norm( ) - .02, 1e-15 );
    models = { preciseOnly };
}

BOOST_AUTO_TEST_CASE( provider_geometry_is_fixed_and_spatial_binding_is_explicit )
{
    const auto converter = loadingClock( );
    const auto source = data( );
    StationLoadingDeformation model( source, settings( source ), converter );
    const auto original = station( );
    const Time input = tdb( converter, originTai + 25.0 );
    const auto expected = model.calculateDisplacement( input, original );
    auto geodetic = std::make_shared< GroundStationState >(
            original->getNominalGeodeticPosition( ), coordinate_conversions::geodetic_position, original->getBodySurface( ) );
    geodetic->setSiteId( "ENVIRONMENT_STATION" );
    BOOST_CHECK_SMALL( ( model.calculateDisplacement( input, geodetic ) - expected ).norm( ), 1e-15 );
    const auto nearby = station( providerPoint + Eigen::Vector3d( 10, -20, 30 ) );
    BOOST_CHECK_SMALL( ( model.calculateDisplacement( input, nearby ) - expected ).norm( ), 1e-15 );
    const Eigen::Vector3d local = original->getConstantRotationMatrixFromBodyFixedToTopocentricFrame( ) * expected;
    const Eigen::Vector3d raw = referenceBasis( ).transpose( ) * expected;
    BOOST_CHECK_GT( std::abs( local.z( ) - raw.x( ) ), 1e-6 );  // Oblate geodetic up differs from source radial.
    auto wrongId = station( );
    wrongId->setSiteId( source.sites.begin( )->first );
    BOOST_CHECK_THROW( model.calculateDisplacement( input, wrongId ), std::out_of_range );
    BOOST_CHECK_THROW( model.calculateDisplacement( input, station( providerPoint + Eigen::Vector3d( 4000, 0, 0 ) ) ), std::out_of_range );
    auto radius = settings( source );
    radius.stationBindings.clear( );
    radius.spatialSelection = StationLoadingSpatialSelection::withinProviderRadius;
    StationLoadingDeformation matched( source, radius, converter );
    BOOST_CHECK_SMALL( ( matched.calculateDisplacement( input, wrongId ) - expected ).norm( ), 1e-15 );
    auto duplicate = source;
    auto other = duplicate.sites.begin( )->second;
    other.identifier = "another_internal_label";
    duplicate.sites.emplace( other.identifier, other );
    StationLoadingDeformation ambiguous( duplicate, radius, converter );
    BOOST_CHECK_THROW( ambiguous.calculateDisplacement( input, original ), std::invalid_argument );
    auto pointOnly = source;
    pointOnly.sites.begin( )->second.validityRadiusMetres = 0;
    StationLoadingDeformation exactPoint( pointOnly, settings( pointOnly ), converter );
    BOOST_CHECK_NO_THROW( exactPoint.calculateDisplacement( input, original ) );
    BOOST_CHECK_THROW( exactPoint.calculateDisplacement( input, station( providerPoint + Eigen::Vector3d( 1e-5, 0, 0 ) ) ),
                       std::out_of_range );
}

BOOST_AUTO_TEST_CASE( source_clocks_are_converted_before_interpolation_through_the_inserted_second )
{
    const auto converter = loadingClock( );
    auto source = data( );
    auto terrestrial = source;
    terrestrial.timeScale = StationLoadingTimeScale::tt;
    for( auto& epoch : terrestrial.epochs ) epoch.secondsOfDay += 32.184L;
    StationLoadingDeformation taiModel( source, settings( source ), converter ), ttModel( terrestrial, settings( terrestrial ), converter );
    for( double offset : { -.5, .5, 1.0 } )
    {
        const auto input = tdb( converter, originTai + offset );
        BOOST_CHECK_SMALL(
                ( taiModel.calculateDisplacementState( input, station( ) ) - ttModel.calculateDisplacementState( input, station( ) ) )
                        .norm( ),
                1e-14 );
    }
    // Equivalent split TT/TAI decimals can differ by long-double arithmetic at a knot.
    // Native queries retain those represented epochs; position remains continuous.
    BOOST_CHECK_SMALL( ( taiModel.calculateDisplacement( tdb( converter, originTai ), station( ) ) -
                         ttModel.calculateDisplacement( tdb( converter, originTai ), station( ) ) )
                               .norm( ),
                       1e-14 );
    auto barycentric = data( );
    barycentric.timeScale = StationLoadingTimeScale::tdb;
    barycentric.nominalSampleIntervalSeconds.reset( );  // Uniform TAI samples have nonuniform TDB spacings.
    for( auto& epoch : barycentric.epochs )
    {
        const Time day( ( epoch.modifiedJulianDay - 51544 ) * 24 - 12, 0.0L );
        epoch.secondsOfDay = static_cast< long double >( tdb( converter, day + epoch.secondsOfDay ) - day );
    }
    StationLoadingDeformation tdbModel( barycentric, settings( barycentric ), converter );
    const auto finalSource = barycentric.epochs.back( );
    const Time finalTdb( ( finalSource.modifiedJulianDay - 51544 ) * 24 - 12, finalSource.secondsOfDay );
    BOOST_CHECK_SMALL( ( tdbModel.calculateDisplacement( finalTdb, station( ) ) - referenceBasis( ) * scalarFunction( 200 ) ).norm( ),
                       1e-15 );
    BOOST_CHECK_THROW( tdbModel.calculateDisplacement( finalTdb + 1e-3L, station( ) ), std::out_of_range );
    source.timeScale = StationLoadingTimeScale::utc;
    source.basis = StationLoadingBasis::cartesian;
    source.epochs = { { 57753, 86300.0L }, { 57754, 0.0L } };
    source.nominalSampleIntervalSeconds = 100.0L;
    auto& samples = source.sites.begin( )->second.samples;
    samples.clear( );
    samples[ 0 ].displacement = { 0, 0, 0 };
    samples[ 1 ].displacement = { .0101, 0, 0 };
    StationLoadingDeformation utcModel( source, settings( source ), converter );
    const auto input = tdb( converter, originTai + .5L );
    const auto state = utcModel.calculateDisplacementState( input, station( ) );
    BOOST_CHECK_SMALL( state.x( ) - .01005, 1e-12 );
    BOOST_CHECK_SMALL( state( 3 ) - .0001 * converter->getGeocentricTtRateWithRespectToTdb( input ), 1e-16 );
}

BOOST_AUTO_TEST_CASE( matching_scoped_origin_translation_is_cartesian_and_applies_its_rate_once )
{
    const auto converter = loadingClock( );
    const auto source = data( );
    auto options = settings( source );
    options.originTranslation = translation( source );
    options.outputOriginDefinition = options.originTranslation->targetOriginDefinition;
    StationLoadingDeformation converted( source, options, converter ), original( source, settings( source ), converter );
    const Time input = tdb( converter, originTai + 25.0 );
    const Eigen::Vector6d difference =
            converted.calculateDisplacementState( input, station( ) ) - original.calculateDisplacementState( input, station( ) );
    BOOST_CHECK_SMALL( ( difference.head< 3 >( ) - Eigen::Vector3d( .00105, -.001925, .0029 ) ).norm( ), 1e-12 );
    BOOST_CHECK_SMALL(
            ( difference.tail< 3 >( ) - Eigen::Vector3d( 2e-6, 3e-6, -4e-6 ) * converter->getGeocentricTtRateWithRespectToTdb( input ) )
                    .norm( ),
            1e-16 );
    auto invalid = options;
    invalid.originTranslation.reset( );
    BOOST_CHECK_THROW( StationLoadingDeformation( source, invalid, converter ), std::invalid_argument );
    auto shift = translation( source );
    shift->sourceOriginDefinition = "CE incorrectly equated with CF";
    invalid.originTranslation = shift;
    BOOST_CHECK_THROW( StationLoadingDeformation( source, invalid, converter ), std::invalid_argument );
    shift = translation( source );
    shift->series.metadata.components = { StationLoadingComponent::ocean };
    invalid.originTranslation = shift;
    BOOST_CHECK_THROW( StationLoadingDeformation( source, invalid, converter ), std::invalid_argument );
}

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
    const auto converter = std::make_shared< TerrestrialTimeScaleConverter >( nullptr, nullptr, nullptr, evaluator );
    const auto polar = std::make_shared< PolarMotionCalculator >( nullptr, nullptr, evaluator );
    const auto celestial = std::make_shared< PrecessionNutationCalculator >( iau_2006, nullptr, nullptr, evaluator );
    const auto rotation = std::make_shared< ephemerides::GcrsToItrsRotationModel >(
            std::make_shared< EarthOrientationAnglesCalculator >( polar, celestial, converter ), tdb_scale );
    simulation_setup::SystemOfBodies bodies( "SSB", "GCRS" );
    bodies.createEmptyBody( "Earth" );
    const auto earth = bodies.at( "Earth" );
    earth->setRotationalEphemeris( rotation );
    earth->setShapeModel( std::make_shared< OblateSpheroidBodyShapeModel >( 6378137, 1.0 / 298.257223563 ) );
    const auto source = data( );
    const auto options = settings( source );
    const auto configured = simulation_setup::stationLoadingBodyShapeDeformation( source, options );
    const auto loaded = simulation_setup::createBodyDeformationModel( configured, "Earth", bodies );
    StationLoadingDeformation direct( source, options, converter );
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
    auto incompatible = source;
    incompatible.metadata.referenceFrame = "unconverted realization";
    BOOST_CHECK_THROW(
            simulation_setup::createBodyDeformationModel(
                    simulation_setup::stationLoadingBodyShapeDeformation( incompatible, settings( incompatible ) ), "Earth", bodies ),
            std::invalid_argument );
}
