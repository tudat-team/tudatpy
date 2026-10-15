#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/basicInputOutput.h"
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
    return paths::getTudatTestDataPath( ) + "pr1092/reference_files/ocean_pole_tide/" + name;
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

    // Daily records with a constant pole, on the UTC days around the test epoch.
    const double utcDay = std::floor( sofa_interface::convertTTtoUTC( tt ) / 86400.0 + 51544.5 );
    input_output::EarthOrientationData data;
    for( int i = -1; i <= 2; ++i )
    {
        input_output::EarthOrientationRecord record;
        record.utcModifiedJulianDate = utcDay + i;
        record.polarMotion = { poleRadians.x( ), poleRadians.y( ) };
        record.ut1MinusUtc = 0.0;
        data.records.emplace( record.utcModifiedJulianDate, record );
    }
    auto shortPeriod = std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
            1.0, 0.0, std::vector< std::string >{}, std::vector< std::string >{}, []( double ) -> Eigen::Vector6d {
                return Eigen::Vector6d::Zero( );
            } );
    return std::make_shared< earth_orientation::PolarMotionCalculator >(
            std::make_shared< earth_orientation::EarthOrientationDataEvaluator >( data ), shortPeriod );
}
}  // namespace

//! Compare ocean-pole displacement with the published IERS series and its historical constants.
BOOST_AUTO_TEST_CASE( publishedIersSeriesOldTestParametersOnly )
{
    const auto grid = readOceanPoleTideGrid( fixturePath( "iers_grid_extract.txt" ) );
    const auto coefficients = interpolateOceanPoleTideCoefficients( grid, -43.75, 232.25 );
    std::ifstream input( fixturePath( "iers_published_test.txt" ) );

    // Read the archived published series before comparing its displacement columns.
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

//! Check ocean-pole equation signs, units and scaling with independent complex arithmetic.
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
}

//! Interpolate ocean-pole coefficients across longitude wrap and at the polar rows.
BOOST_AUTO_TEST_CASE( interpolationPeriodicLongitudeBoundariesAndPolarClamp )
{
    const auto grid = syntheticGrid( );

    // Interpolate ocean-pole coefficients across longitude wrap and at the polar rows.
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
    const auto extract = readOceanPoleTideGrid( fixturePath( "iers_grid_extract.txt" ) );

    // The centre of the published grid cell is the mean of its four corner coefficients.
    BOOST_CHECK_SMALL( ( interpolateOceanPoleTideCoefficients( extract, -43.5, 232.5 ) -
                         Eigen::Matrix< double, 3, 2 >( ( extract.coefficients[ 0 ] + extract.coefficients[ 1 ] +
                                                          extract.coefficients[ 2 ] + extract.coefficients[ 3 ] ) /
                                                        4 ) )
                               .norm( ),
                       1e-14 );
}

//! Use geocentric ocean-pole geometry independently of the station shape representation.
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

    // Cartesian displacements agree; only their local projections depend on the station shape.
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
}

//! Compare future ocean-pole displacement with independently computed secular-pole wobble.
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

//! Compare ocean-pole models constructed from a grid file and the same preloaded data.
BOOST_AUTO_TEST_CASE( nativeFactoryMatchesLoadedGridResults )
{
    using namespace simulation_setup;
    const auto grid = readOceanPoleTideGrid( fixturePath( "iers_grid_extract.txt" ) );
    const auto pole = poleCalculator( 0, Eigen::Vector2d( 1e-6, 2e-6 ) );
    const auto angles = std::make_shared< earth_orientation::EarthOrientationAnglesCalculator >( pole, nullptr, nullptr );
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    bodies.at( "Earth" )->setRotationalEphemeris( std::make_shared< ephemerides::GcrsToItrsRotationModel >( angles ) );
    const auto preloaded = createBodyDeformationModel( oceanPoleTideBodyShapeDeformationFromData( grid ), "Earth", bodies );

    const auto file =
            createBodyDeformationModel( oceanPoleTideBodyShapeDeformation( fixturePath( "iers_grid_extract.txt" ) ), "Earth", bodies );

    const Eigen::Vector3d spherical( 6378137, -43.75 * mathematical_constants::PI / 180, 232.25 * mathematical_constants::PI / 180 );
    const auto station = std::make_shared< ground_stations::GroundStationState >( spherical, coordinate_conversions::spherical_position );

    // The loaded file and preloaded grid produce the same Cartesian displacement.
    close( preloaded->calculateDisplacement( tdbAtTt( 0 ), station ), file->calculateDisplacement( tdbAtTt( 0 ), station ) );
}
