#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <complex>
#include "tudat/astro/earth_orientation/eopReader.h"
#include "tudat/interface/sofa/fundamentalArguments.h"
#include <sstream>
#include "tudat/astro/ground_stations/oceanPoleTideDeformation.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"

using namespace tudat;
using namespace basic_astrodynamics;
using namespace input_output;
namespace
{
std::string fixturePath( const std::string& name )
{
    return std::string( __FILE__ ).substr( 0, std::string( __FILE__ ).find_last_of( "/\\" ) ) +
            "/../../../reference_files/ocean_pole_tide/" + name;
}
void close( const Eigen::Vector3d& actual, const Eigen::Vector3d& expected, const double tolerance = 1e-14 )
{
    for( int i = 0; i < 3; ++i )
    {
        BOOST_CHECK_SMALL( actual[ i ] - expected[ i ], tolerance );
    }
}
// Published fixtures and injected mean-pole values use TT. The model receives TDB.
Time tdbAtTt( const double tt,
              const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter = getGeocentricTimeConverter( ) )
{
    return converter->getCurrentTime< Time >( basic_astrodynamics::tt_scale, basic_astrodynamics::tdb_scale, Time( tt ) );
}
std::shared_ptr< ground_stations::GroundStationState > stationAt( const Eigen::Vector3d& position )
{
    return std::make_shared< ground_stations::GroundStationState >( position, coordinate_conversions::cartesian_position );
}
OceanPoleTideGrid syntheticGrid( )
{
    OceanPoleTideGrid grid;
    grid.longitudeCount = 4;
    grid.latitudeCount = 3;
    grid.firstLongitudeDegrees = 45;
    grid.longitudeStepDegrees = 90;
    grid.firstLatitudeDegrees = -89.75;
    grid.latitudeStepDegrees = 89.75;
    for( int row = 0; row < 3; ++row )
    {
        for( int col = 0; col < 4; ++col )
        {
            grid.coefficients.push_back( Eigen::Matrix< double, 3, 2 >::Constant( 10 * row + col ) );
        }
    }
    return grid;
}
std::shared_ptr< earth_orientation::PolarMotionCalculator > poleCalculator( const double tt, const Eigen::Vector2d& poleRadians )
{
    const double utc = sofa_interface::convertTTtoUTC( tt );
    const std::map< double, Eigen::Vector2d > samples = { { utc - 86400, poleRadians }, { utc + 86400, poleRadians } };
    auto daily = std::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector2d > >( samples );
    auto shortPeriod = std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
            1.0, 0.0, std::vector< std::string >{}, std::vector< std::string >{}, []( double ) -> Eigen::Vector6d {
                return Eigen::Vector6d::Zero( );
            } );
    return std::make_shared< earth_orientation::PolarMotionCalculator >( daily, shortPeriod );
}
}  // namespace

BOOST_AUTO_TEST_CASE( publishedIersSeriesOldTestParametersOnly )
{
    const auto grid = readOceanPoleTideGrid( fixturePath( "iers_grid_extract.txt" ) );
    const auto coefficients = interpolateOceanPoleTideCoefficients( grid, -43.75, 232.25 );
    std::ifstream input( fixturePath( "iers_published_test.txt" ) );
    BOOST_REQUIRE( input.good( ) );
    std::string line;
    int cases = 0;
    const double radiansPerArcsecond = mathematical_constants::PI / 648000.0;
    while( std::getline( input, line ) )
    {
        std::istringstream row( line );
        double mjd, xbar, ybar, xp, yp, m1, m2, radial, north, east;
        if( !( row >> mjd >> xbar >> ybar >> xp >> yp >> m1 >> m2 >> radial >> north >> east ) )
        {
            continue;
        }
        // The independently published series specifies the old 2003 mean pole and G=6.673e-11.
        // Printed m1/m2 are rounded to seven significant digits, and the output also has seven digits, giving sub-nanometre residuals.
        close( calculateOceanPoleTideLoading( coefficients, Eigen::Vector2d( m1, m2 ) * radiansPerArcsecond, 5339.4043696 ),
               Eigen::Vector3d( radial, north, east ),
               1e-9 );
        ++cases;
    }
    BOOST_CHECK_EQUAL( cases, 731 );
}

BOOST_AUTO_TEST_CASE( equationUnitsSignsScaleAndZero )
{
    Eigen::Matrix< double, 3, 2 > coefficients;
    coefficients << 1, 2, 3, 4, 5, 6;
    // Frozen via independent complex arithmetic: Re(conj(u)*(m1+i*m2)*conj(gamma))*K.
    close( calculateOceanPoleTideLoading( coefficients, Eigen::Vector2d( 1e-6, -2e-6 ) ),
           Eigen::Vector3d( -0.01108352543805147, -0.01853662753949921, -0.025989729640946948 ) );
    close( calculateOceanPoleTideLoading( coefficients, Eigen::Vector2d::Zero( ) ), Eigen::Vector3d::Zero( ) );
    close( calculateOceanPoleTideLoading( coefficients, Eigen::Vector2d( -1e-6, 2e-6 ) ),
           -calculateOceanPoleTideLoading( coefficients, Eigen::Vector2d( 1e-6, -2e-6 ) ) );
    close( calculateOceanPoleTideLoading( coefficients, Eigen::Vector2d( 1e-6, -2e-6 ), 10680.857124459354 ),
           2 * calculateOceanPoleTideLoading( coefficients, Eigen::Vector2d( 1e-6, -2e-6 ) ) );
    BOOST_CHECK_THROW( calculateOceanPoleTideLoading( coefficients, Eigen::Vector2d( NAN, 0 ) ), std::runtime_error );
    BOOST_CHECK_THROW( calculateOceanPoleTideLoading( coefficients, Eigen::Vector2d::Zero( ), -1 ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( interpolationPeriodicLongitudeBoundariesAndPolarClamp )
{
    const auto grid = syntheticGrid( );
    BOOST_CHECK_SMALL(
            ( interpolateOceanPoleTideCoefficients( grid, -44.875, 90 ) - Eigen::Matrix< double, 3, 2 >::Constant( 5.5 ) ).norm( ), 1e-14 );
    BOOST_CHECK_SMALL( ( interpolateOceanPoleTideCoefficients( grid, 0, 0 ) - Eigen::Matrix< double, 3, 2 >::Constant( 11.5 ) ).norm( ),
                       1e-14 );
    for( double lon : { -720.0, -360.0, 0.0, 360.0, 720.0 } )
    {
        BOOST_CHECK_SMALL(
                ( interpolateOceanPoleTideCoefficients( grid, 0, lon ) - Eigen::Matrix< double, 3, 2 >::Constant( 11.5 ) ).norm( ), 1e-14 );
    }
    BOOST_CHECK_SMALL( ( interpolateOceanPoleTideCoefficients( grid, 89.75, 315 ) - Eigen::Matrix< double, 3, 2 >::Constant( 23 ) ).norm( ),
                       1e-14 );
    for( double lat : { -90.0, 90.0 } )
    {
        const double edge = lat < 0 ? -89.75 : 89.75;
        BOOST_CHECK_SMALL(
                ( interpolateOceanPoleTideCoefficients( grid, lat, 45 ) - interpolateOceanPoleTideCoefficients( grid, edge, 45 ) ).norm( ),
                1e-14 );
    }
    BOOST_CHECK_THROW( interpolateOceanPoleTideCoefficients( grid, 91, 45 ), std::runtime_error );
    BOOST_CHECK_THROW( interpolateOceanPoleTideCoefficients( grid, 0, NAN ), std::runtime_error );
    const auto extract = readOceanPoleTideGrid( fixturePath( "iers_grid_extract.txt" ) );
    BOOST_CHECK_THROW( interpolateOceanPoleTideCoefficients( extract, -43.5, 240 ), std::out_of_range );
    // Latitudes more than one grid step outside a regional extract are rejected rather than clamped.
    BOOST_CHECK_THROW( interpolateOceanPoleTideCoefficients( extract, 0, 232.5 ), std::out_of_range );
    BOOST_CHECK_SMALL( ( interpolateOceanPoleTideCoefficients( extract, -43.5, 232.5 ) -
                         Eigen::Matrix< double, 3, 2 >( ( extract.coefficients[ 0 ] + extract.coefficients[ 1 ] +
                                                          extract.coefficients[ 2 ] + extract.coefficients[ 3 ] ) /
                                                        4 ) )
                               .norm( ),
                       1e-14 );
}

BOOST_AUTO_TEST_CASE( modelSecularPoleAndGeocentricBasisAcrossShapeRepresentations )
{
    const double pi = mathematical_constants::PI;
    const auto grid = syntheticGrid( );
    // At J2000 the 2018 secular pole is (55,320.5) mas. Inject exactly m1=1e-6,m2=-2e-6 radians.
    const auto pole = poleCalculator( 0, Eigen::Vector2d( 55, 320.5 ) * ( pi / 648000000.0 ) + Eigen::Vector2d( 1e-6, 2e-6 ) );
    OceanPoleTideDeformation model( grid, pole );
    const Eigen::Vector3d position( -2328813.9273189958, 3586058.9739219621, 4717062.1163217751 );
    auto sphere = std::make_shared< ground_stations::GroundStationState >( position, coordinate_conversions::cartesian_position );
    auto shape = std::make_shared< OblateSpheroidBodyShapeModel >( 6378137, 1.0 / 298.257223563 );
    auto oblate = std::make_shared< ground_stations::GroundStationState >( position, coordinate_conversions::cartesian_position, shape );
    const double lat = std::atan2( position.z( ), std::hypot( position.x( ), position.y( ) ) );
    const double lon = std::atan2( position.y( ), position.x( ) );
    const auto u = calculateOceanPoleTideLoading( interpolateOceanPoleTideCoefficients( grid, lat * 180 / pi, lon * 180 / pi ),
                                                  Eigen::Vector2d( 1e-6, -2e-6 ) );
    const Eigen::Vector3d radial( std::cos( lat ) * std::cos( lon ), std::cos( lat ) * std::sin( lon ), std::sin( lat ) );
    const Eigen::Vector3d north( -std::sin( lat ) * std::cos( lon ), -std::sin( lat ) * std::sin( lon ), std::cos( lat ) );
    const Eigen::Vector3d east( -std::sin( lon ), std::cos( lon ), 0 );
    const Eigen::Vector3d expected = u.x( ) * radial + u.y( ) * north + u.z( ) * east;
    close( model.calculateDisplacement( tdbAtTt( 0 ), sphere ), expected );
    close( model.calculateDisplacement( tdbAtTt( 0 ), oblate ), expected );
    close( model.calculateDisplacement( tdbAtTt( 0 ), stationAt( position ) ), expected );
    close( oblate->getRotationFromBodyFixedToTopocentricFrame( static_cast< double >( tdbAtTt( 0 ) ) ) *
                   model.calculateDisplacement( tdbAtTt( 0 ), oblate ),
           oblate->getRotationFromBodyFixedToTopocentricFrame( 0 ) * expected );
    BOOST_CHECK_GT( ( sphere->getRotationFromBodyFixedToTopocentricFrame( static_cast< double >( tdbAtTt( 0 ) ) ) *
                              model.calculateDisplacement( tdbAtTt( 0 ), sphere ) -
                      oblate->getRotationFromBodyFixedToTopocentricFrame( static_cast< double >( tdbAtTt( 0 ) ) ) *
                              model.calculateDisplacement( tdbAtTt( 0 ), oblate ) )
                            .norm( ),
                    1e-6 );
    OceanPoleTideDeformation zero( grid, poleCalculator( 0, Eigen::Vector2d( 55, 320.5 ) * ( pi / 648000000.0 ) ) );
    close( zero.calculateDisplacement( tdbAtTt( 0 ), stationAt( position ) ), Eigen::Vector3d::Zero( ) );
    BOOST_CHECK_THROW( model.calculateDisplacement( Time( NAN ), sphere ), std::invalid_argument );
    BOOST_CHECK_THROW( model.calculateDisplacement( tdbAtTt( 0 ), stationAt( Eigen::Vector3d::Zero( ) ) ), std::runtime_error );
    BOOST_CHECK_THROW( OceanPoleTideDeformation( grid, nullptr ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( futureSecularPoleWithIndependentMixedWobble )
{
    auto grid = syntheticGrid( );
    Eigen::Matrix< double, 3, 2 > coefficients;
    coefficients << 1, 2, 3, 4, 5, 6;
    for( auto& node : grid.coefficients )
    {
        node = coefficients;
    }
    const double tt = 20 * physical_constants::JULIAN_YEAR;
    // 2018 secular pole at 20 years: (88.54,389.70) mas. Here m1=-100 mas,m2=+200 mas.
    const auto pole = poleCalculator( tt, Eigen::Vector2d( -11.46, 189.70 ) * ( mathematical_constants::PI / 648000000.0 ) );
    OceanPoleTideDeformation model( grid, pole );
    const Eigen::Vector3d position( 6378137, 0, 0 );
    // Independent complex arithmetic yields radial,north,east. At this point Cartesian = radial,east,north.
    close( model.calculateDisplacement( tdbAtTt( tt ), stationAt( position ) ),
           Eigen::Vector3d( 0.005373444767292916, 0.012600176498269108, 0.008986810632781011 ),
           1e-14 );
}

BOOST_AUTO_TEST_CASE( nativeFactoryPreloadedFileAndConfigurationErrors )
{
    using namespace simulation_setup;
    const auto grid = readOceanPoleTideGrid( fixturePath( "iers_grid_extract.txt" ) );
    const auto pole = poleCalculator( 0, Eigen::Vector2d( 1e-6, 2e-6 ) );
    const auto angles = std::make_shared< earth_orientation::EarthOrientationAnglesCalculator >( pole, nullptr, nullptr );
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    bodies.at( "Earth" )->setRotationalEphemeris( std::make_shared< ephemerides::GcrsToItrsRotationModel >( angles ) );
    const auto preloaded = createBodyDeformationModel( oceanPoleTideBodyShapeDeformationFromData( grid ), "Earth", bodies );
    // Filename settings retain preloaded coefficients and report unavailable files immediately.
    const auto fileSettings = oceanPoleTideBodyShapeDeformation( fixturePath( "iers_grid_extract.txt" ) );

    BOOST_CHECK_EQUAL( fileSettings->getGrid( ).coefficients.size( ), grid.coefficients.size( ) );
    BOOST_CHECK_THROW( oceanPoleTideBodyShapeDeformation( fixturePath( "iers_grid_extract.txt.missing" ) ), std::runtime_error );

    const auto file = createBodyDeformationModel( fileSettings, "Earth", bodies );

    const Eigen::Vector3d spherical( 6378137, -43.75 * mathematical_constants::PI / 180, 232.25 * mathematical_constants::PI / 180 );
    const auto station = std::make_shared< ground_stations::GroundStationState >( spherical, coordinate_conversions::spherical_position );
    close( preloaded->calculateDisplacement( tdbAtTt( 0 ), station ), file->calculateDisplacement( tdbAtTt( 0 ), station ) );
    const auto typed = std::dynamic_pointer_cast< OceanPoleTideDeformation >( preloaded );
    BOOST_REQUIRE( typed );
    BOOST_CHECK( typed->getGrid( ).referenceOrigin == OceanLoadingReferenceOrigin::wholeEarthCentreOfMass );
    BOOST_CHECK_THROW( createBodyDeformationModel( oceanPoleTideBodyShapeDeformationFromData( grid ), "Mars", bodies ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            createBodyDeformationModel( std::make_shared< BodyDeformationSettings >( gravitation::ocean_pole_tide ), "Earth", bodies ),
            std::runtime_error );
    bodies.at( "Earth" )->setRotationalEphemeris( nullptr );
    BOOST_CHECK_THROW( createBodyDeformationModel( oceanPoleTideBodyShapeDeformationFromData( grid ), "Earth", bodies ),
                       std::runtime_error );
}

BOOST_AUTO_TEST_CASE( nativeFactoryContinuousEopAcrossInsertedUtcSecond )
{
    using namespace simulation_setup;
    using namespace earth_orientation;
    // Controlled daily observations are affine in TAI, including the 86401-second day.
    EarthOrientationData data;
    const Time midnightTai = Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0;
    const Eigen::Vector2d referencePole( 1.0e-6, 2.0e-6 ), slope( 2.0e-12, -3.0e-12 );
    for( int i = 0; i < 5; ++i )
    {
        EarthOrientationRecord record;
        record.utcModifiedJulianDate = 57752.0 + i;
        const Time utc( ( record.utcModifiedJulianDate - 51544.5 ) * 86400.0 );
        const Time tai = sofa_interface::convertUTCtoTAI( utc );
        const double dt = static_cast< double >( tai - midnightTai );
        const Eigen::Vector2d pole = referencePole + dt * slope;
        record.polarMotion = { pole.x( ), pole.y( ) };
        record.ut1MinusUtc = -36.4 - 2.0e-8 * dt + static_cast< double >( tai - utc );
        record.celestialPoleOffsets = std::array< double, 2 >{ 0.0, 0.0 };
        data.records.emplace( record.utcModifiedJulianDate, record );
    }
    const auto angles = createStandardEarthOrientationCalculator( std::make_shared< EOPReader >( data ) );
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    bodies.at( "Earth" )->setRotationalEphemeris( std::make_shared< ephemerides::GcrsToItrsRotationModel >( angles ) );
    auto grid = syntheticGrid( );
    Eigen::Matrix< double, 3, 2 > coefficients;
    coefficients << 1, 2, 3, 4, 5, 6;
    for( auto& node : grid.coefficients )
    {
        node = coefficients;
    }
    const auto model = createBodyDeformationModel( oceanPoleTideBodyShapeDeformationFromData( grid ), "Earth", bodies );
    const Eigen::Vector3d position( 6378137.0, 0.0, 0.0 );
    BOOST_CHECK_THROW( EarthOrientationDataEvaluator::convertTaiToScalarUtc( midnightTai - 0.5 ), std::runtime_error );
    const double midnightTt = static_cast< double >( basic_astrodynamics::convertTAItoTT( midnightTai ) );
    for( double offset : { -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0 } )
    {
        const double tt = midnightTt + offset;
        const Time tai = basic_astrodynamics::convertTTtoTAI( Time( tt ) );
        const Time ut1 = angles->getTerrestrialTimeScaleConverter( )->getCurrentTime< Time >( tai_scale, ut1_scale, tai );
        const auto arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( tt, tt, static_cast< double >( ut1 ) );
        const Eigen::Vector2d shortPeriod =
                angles->getPolarMotionCalculator( )->getShortPeriodPolarMotionCalculator( )->getCorrectionsFromFundamentalArgument(
                        arguments );
        const Eigen::Vector2d pole = referencePole + static_cast< double >( tai - midnightTai ) * slope + shortPeriod;
        const double years = tt / physical_constants::JULIAN_YEAR;
        const Eigen::Vector2d meanPole =
                Eigen::Vector2d( 55.0 + 1.677 * years, 320.5 + 3.460 * years ) * ( mathematical_constants::PI / 648000000.0 );
        const std::complex< double > wobble( pole.x( ) - meanPole.x( ), -( pole.y( ) - meanPole.y( ) ) );
        Eigen::Vector3d loading;
        for( int component = 0; component < 3; ++component )
        {
            loading[ component ] = 5340.428562229677 *
                    std::real( std::conj( std::complex< double >( coefficients( component, 0 ), coefficients( component, 1 ) ) ) * wobble *
                               std::conj( std::complex< double >( 0.687, 0.0036 ) ) );
        }
        // Independent complex loading formula and equatorial radial/east/north basis.
        close( model->calculateDisplacement( tdbAtTt( tt, angles->getTerrestrialTimeScaleConverter( ) ), stationAt( position ) ),
               Eigen::Vector3d( loading.x( ), loading.z( ), loading.y( ) ),
               1.0e-13 );
    }
}
