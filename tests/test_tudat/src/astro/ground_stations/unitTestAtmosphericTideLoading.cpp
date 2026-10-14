#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/ground_stations/atmosphericTideLoading.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "../../../reference_files/atmospheric_tides/origin_translation_reference.h"
using namespace tudat;
using namespace basic_astrodynamics;
using namespace input_output;
namespace
{
void close( const Eigen::Vector3d& actual, const Eigen::Vector3d& expected, double tolerance = 1e-13 )
{
    for( int i = 0; i < 3; ++i )
    {
        BOOST_CHECK_SMALL( actual[ i ] - expected[ i ], tolerance );
    }
}
AtmosphericTideGrid grid( )
{
    AtmosphericTideGrid result;
    result.longitudesDegrees = { 0, 90, 180, 270, 360 };
    result.latitudesDegrees = { 45, 0, -45 };
    for( int i = 0; i < 5; ++i )
    {
        for( int j = 0; j < 3; ++j )
        {
            result.coefficientsMetres.push_back( AtmosphericTideCoefficients::Constant( ( i == 4 ? 0 : i ) + 10 * j ) );
        }
    }
    return result;
}
}  // namespace
BOOST_AUTO_TEST_CASE( directQuarterDaySignsAndS2HalfDayPeriod )
{
    AtmosphericTideCoefficients coefficients;
    coefficients << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
    close( synthesizeAtmosphericTideLoading( coefficients, -43200 ), Eigen::Vector3d( 4, 12, 20 ) );
    close( synthesizeAtmosphericTideLoading( coefficients, -21600 ), Eigen::Vector3d( -1, -1, -1 ) );
    close( synthesizeAtmosphericTideLoading( coefficients, 0 ), Eigen::Vector3d( 2, 2, 2 ) );
    close( synthesizeAtmosphericTideLoading( coefficients, 21600 ), Eigen::Vector3d( -5, -13, -21 ) );
    auto s1 = AtmosphericTideCoefficients::Zero( ).eval( );
    s1.col( 0 ) << 1, 2, 3;
    s1.col( 1 ) << 4, 5, 6;
    close( synthesizeAtmosphericTideLoading( s1, -43200 ), Eigen::Vector3d( 1, 2, 3 ) );
    close( synthesizeAtmosphericTideLoading( s1, -21600 ), Eigen::Vector3d( 4, 5, 6 ) );
    close( synthesizeAtmosphericTideLoading( s1, 0 ), Eigen::Vector3d( -1, -2, -3 ) );
    close( synthesizeAtmosphericTideLoading( s1, 21600 ), Eigen::Vector3d( -4, -5, -6 ) );
    auto s2 = AtmosphericTideCoefficients::Zero( ).eval( );
    s2.col( 2 ) << 1, 2, 3;
    s2.col( 3 ) << 4, 5, 6;
    close( synthesizeAtmosphericTideLoading( s2, -43200 ), Eigen::Vector3d( 1, 2, 3 ) );
    close( synthesizeAtmosphericTideLoading( s2, -32400 ), Eigen::Vector3d( 4, 5, 6 ) );
    close( synthesizeAtmosphericTideLoading( s2, -21600 ), Eigen::Vector3d( -1, -2, -3 ) );
    close( synthesizeAtmosphericTideLoading( s2, -10800 ), Eigen::Vector3d( -4, -5, -6 ) );
    close( synthesizeAtmosphericTideLoading( s2, 0 ), Eigen::Vector3d( 1, 2, 3 ) );
    close( synthesizeAtmosphericTideLoading( AtmosphericTideCoefficients::Zero( ), 0 ), Eigen::Vector3d::Zero( ) );
    BOOST_CHECK_THROW( synthesizeAtmosphericTideLoading( coefficients, NAN ), std::runtime_error );
}
BOOST_AUTO_TEST_CASE( analyticDerivativeAndEpochReduction )
{
    AtmosphericTideCoefficients coefficients;
    coefficients << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
    const double omega = 2 * mathematical_constants::PI / 86400;
    close( differentiateAtmosphericTideLoadingWrtUt1( coefficients, -43200 ), omega * Eigen::Vector3d( 10, 22, 34 ) );
    close( differentiateAtmosphericTideLoadingWrtUt1( coefficients, 0 ), omega * Eigen::Vector3d( 6, 10, 14 ) );
    close( synthesizeAtmosphericTideLoading( coefficients, 86400. * 1000000 + 12345 ),
           synthesizeAtmosphericTideLoading( coefficients, 12345 ) );
    close( synthesizeAtmosphericTideLoading( coefficients, -86400. * 1000000 + 12345 ),
           synthesizeAtmosphericTideLoading( coefficients, 12345 ) );
    const double t = 12345, h = .1;
    close( differentiateAtmosphericTideLoadingWrtUt1( coefficients, t ),
           ( synthesizeAtmosphericTideLoading( coefficients, t + h ) - synthesizeAtmosphericTideLoading( coefficients, t - h ) ) /
                   ( 2 * h ),
           1e-11 );
}
BOOST_AUTO_TEST_CASE( interpolationIndependentLinearCasesAndRawEndpointRetention )
{
    const auto raw = grid( );
    AtmosphericTideGridInterpolator model( raw );
    BOOST_CHECK_SMALL( ( model.interpolate( 22.5, 45 ) - AtmosphericTideCoefficients::Constant( 5.5 ) ).norm( ), 1e-14 );
    BOOST_CHECK_SMALL( ( model.interpolate( 0, -45 ) - AtmosphericTideCoefficients::Constant( 11.5 ) ).norm( ), 1e-14 );
    for( double lon : { -720., -360., 0., 360., 720. } )
    {
        BOOST_CHECK_SMALL( ( model.interpolate( 0, lon ) - AtmosphericTideCoefficients::Constant( 10 ) ).norm( ), 1e-14 );
    }
    BOOST_CHECK_THROW( model.interpolate( 90, 0 ), std::runtime_error );
    auto inconsistent = raw;
    inconsistent.coefficientsMetres.back( )( 0, 0 ) += 8e-7;
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideGridInterpolator( inconsistent ) ), std::runtime_error );
    AtmosphericTideGridInterpolator averaged( inconsistent, AtmosphericTideEndpointPolicy::averageSeamAndPolarCartesianApproximation );
    BOOST_CHECK_CLOSE( averaged.getMaximumSeamDiscrepancyMetres( ), 8e-7, 1e-6 );
    BOOST_CHECK_EQUAL( averaged.getProviderGrid( ).longitudesDegrees.size( ), 5 );
    BOOST_CHECK_EQUAL( averaged.getProviderGrid( ).coefficientsMetres.back( )( 0, 0 ), inconsistent.coefficientsMetres.back( )( 0, 0 ) );
    BOOST_CHECK_SMALL( averaged.interpolate( -45, 0 )( 0, 0 ) - 20 - 4e-7, 1e-14 );
}
BOOST_AUTO_TEST_CASE( polarCartesianAveragingExplicitlyChangesOnlyModelCopy )
{
    auto raw = grid( );
    raw.latitudesDegrees = { 90, 0, -90 };
    // A constant radial value with zero horizontal components is Cartesian-consistent at each pole.
    for( std::size_t i = 0; i < raw.longitudesDegrees.size( ); ++i )
    {
        for( std::size_t j : { 0u, 2u } )
        {
            auto& c = raw.coefficientsMetres[ i * 3 + j ];
            c.setZero( );
            c.row( 0 ).setConstant( 1 );
            const double longitude = raw.longitudesDegrees[ i ] * mathematical_constants::PI / 180;
            const double northSign = j == 0 ? -1 : 1;
            c.row( 1 ).setConstant( northSign * ( .003 * std::cos( longitude ) - .005 * std::sin( longitude ) ) );
            c.row( 2 ).setConstant( -.003 * std::sin( longitude ) - .005 * std::cos( longitude ) );
        }
    }
    AtmosphericTideGridInterpolator exact( raw );
    raw.coefficientsMetres[ 0 ]( 0, 0 ) += 4e-7;
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideGridInterpolator( raw ) ), std::runtime_error );
    AtmosphericTideGridInterpolator averaged( raw, AtmosphericTideEndpointPolicy::averageSeamAndPolarCartesianApproximation );
    BOOST_CHECK_GT( averaged.getMaximumPolarCartesianDiscrepancyMetres( ), 3.9e-7 );
    BOOST_CHECK_EQUAL( averaged.getProviderGrid( ).coefficientsMetres[ 0 ]( 0, 0 ), 1 + 4e-7 );
    for( double lon : { 0., 25.5, 45., 90., 180., 270. } )
    {
        const auto local = averaged.interpolate( 90, lon ).col( 0 ).eval( );
        const double angle = lon * mathematical_constants::PI / 180;
        close( Eigen::Vector3d( -local.y( ) * std::cos( angle ) - local.z( ) * std::sin( angle ),
                                -local.y( ) * std::sin( angle ) + local.z( ) * std::cos( angle ),
                                local.x( ) ),
               Eigen::Vector3d( .003, -.005, 1 + 5e-8 ) );
    }
}
BOOST_AUTO_TEST_CASE( realGridNodesAndIndependentDecimalInterpolation )
{
    const std::string base = std::string( __FILE__ ).substr( 0, std::string( __FILE__ ).find_last_of( "/\\" ) );
    const auto raw = readAtmosphericTideGrid( base + "/../../../reference_files/atmospheric_tides/rp03_ce_grid_extract.txt" );
    AtmosphericTideGridInterpolator model( raw );
    BOOST_CHECK_SMALL( ( model.interpolate( -19, 47 ) - raw.coefficientsMetres.front( ) ).norm( ), 1e-18 );
    const auto polarRaw = readAtmosphericTideGrid( base + "/../../../reference_files/atmospheric_tides/rp03_ce_polar_node_extract.txt" );
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideGridInterpolator( polarRaw ) ), std::runtime_error );
    AtmosphericTideGridInterpolator polarApproximation( polarRaw,
                                                        AtmosphericTideEndpointPolicy::averageSeamAndPolarCartesianApproximation );
    BOOST_CHECK_THROW( polarApproximation.interpolate( 90, 10 ), std::runtime_error );
    // Independent Decimal interpolation of four original provider rows; this is not reference-software conformance.
    AtmosphericTideCoefficients expected;
    expected << -.229413396401398, -.140409687704627, 1.040762193571904, .276668845291711, .059365717547470, .059123892373735,
            -.033761020617976, -.014395543429494, .007941624058988, -.103894064696771, -.016544820023036, .059869142014747;
    BOOST_CHECK_SMALL( ( model.interpolate( -19.018307, 47.22921 ) - expected * 1e-3 ).norm( ), 1e-16 );
}

BOOST_AUTO_TEST_CASE( publishedGridNodeRadialQuarterDayValues )
{
    const std::string base = std::string( __FILE__ ).substr( 0, std::string( __FILE__ ).find_last_of( "/\\" ) );
    const auto grid = readAtmosphericTideGrid( base + "/../../../reference_files/atmospheric_tides/rp03_ce_polar_node_extract.txt" );
    // Direct arithmetic on the published Figure24.17 four radial coefficients, not a self-generated model fixture.
    const auto& coefficients = grid.coefficientsMetres.front( );
    BOOST_CHECK_SMALL( synthesizeAtmosphericTideLoading( coefficients, -43200 ).x( ) + .0000604, 1e-18 );
    BOOST_CHECK_SMALL( synthesizeAtmosphericTideLoading( coefficients, -21600 ).x( ) - .0001548, 1e-18 );
    BOOST_CHECK_SMALL( synthesizeAtmosphericTideLoading( coefficients, 0 ).x( ) + .0000970, 1e-18 );
    BOOST_CHECK_SMALL( synthesizeAtmosphericTideLoading( coefficients, 21600 ).x( ) - .0000026, 1e-18 );
}

BOOST_AUTO_TEST_CASE( table76CartesianSignsMidnightAndIndependentFourierReference )
{
    // Exact quarter-day arithmetic on Table7.6, which is already in metres and Cartesian X/Y/Z.
    close( calculateIers2010Rp03AtmosphericGeocentreTranslation( -43200 ), Eigen::Vector3d( .00035660, -.00105457, -.000108447 ), 1e-17 );
    close( calculateIers2010Rp03AtmosphericGeocentreTranslation( -21600 ), Eigen::Vector3d( -.00091333, .00009109, .000128514 ), 1e-17 );
    close( calculateIers2010Rp03AtmosphericGeocentreTranslation( 0 ), Eigen::Vector3d( -.00006716, .00040075, -.000084095 ), 1e-17 );
    close( calculateIers2010Rp03AtmosphericGeocentreTranslation( 21600 ), Eigen::Vector3d( .00062389, .00056273, .000064028 ), 1e-17 );
    for( const auto& reference : atmosphericOriginReferences )
    {
        close( calculateIers2010Rp03AtmosphericGeocentreTranslation( reference.ut1 ),
               Eigen::Vector3d( reference.cartesian[ 0 ], reference.cartesian[ 1 ], reference.cartesian[ 2 ] ),
               1e-17 );
    }
    BOOST_CHECK_THROW( calculateIers2010Rp03AtmosphericGeocentreTranslation( NAN ), std::runtime_error );
}
BOOST_AUTO_TEST_CASE( equation18CeToCmAdditionInCartesianAndExplicitLocalFrame )
{
    // Explicit Cartesian test position. These local projections establish no provider grid latitude convention.
    const Eigen::Vector3d position( 4500000, -3000000, 2800000 );
    const double latitude = std::atan2( position.z( ), std::hypot( position.x( ), position.y( ) ) );
    const double longitude = std::atan2( position.y( ), position.x( ) );
    Eigen::Matrix3d radialNorthEastToCartesian;
    radialNorthEastToCartesian.col( 0 ) << std::cos( latitude ) * std::cos( longitude ), std::cos( latitude ) * std::sin( longitude ),
            std::sin( latitude );
    radialNorthEastToCartesian.col( 1 ) << -std::sin( latitude ) * std::cos( longitude ), -std::sin( latitude ) * std::sin( longitude ),
            std::cos( latitude );
    radialNorthEastToCartesian.col( 2 ) << -std::sin( longitude ), std::cos( longitude ), 0;
    const Eigen::Vector3d ceLocal( .000123, -.000456, .000789 );
    const Eigen::Vector3d ceCartesian = radialNorthEastToCartesian * ceLocal;
    for( const auto& reference : atmosphericOriginReferences )
    {
        const Eigen::Vector3d d = calculateIers2010Rp03AtmosphericGeocentreTranslation( reference.ut1 );
        // Eq.18 X_CE=X_CM-d; solving for CM gives CE+d, with d added before local projection.
        const Eigen::Vector3d cmCartesian = ceCartesian + d;
        const Eigen::Vector3d independentlyProjectedTranslation(
                reference.radialNorthEast[ 0 ], reference.radialNorthEast[ 1 ], reference.radialNorthEast[ 2 ] );
        close( radialNorthEastToCartesian.transpose( ) * cmCartesian, ceLocal + independentlyProjectedTranslation, 1e-17 );
        close( cmCartesian - d, ceCartesian, 1e-17 );
    }
}
BOOST_AUTO_TEST_CASE( independentRp03TranslationRatesAndQuarterDaySigns )
{
    for( const auto& reference : atmosphericOriginReferences )
    {
        close( differentiateIers2010Rp03AtmosphericGeocentreTranslationWrtUt1( reference.ut1 ),
               Eigen::Map< const Eigen::Vector3d >( reference.cartesianRateWrtUt1 ),
               1e-22 );
    }
    const double omega = 2 * mathematical_constants::PI / 86400;
    // At midnight: omega*(sinS1+2*sinS2); at noon: omega*(-sinS1+2*sinS2).
    close( differentiateIers2010Rp03AtmosphericGeocentreTranslationWrtUt1( -43200 ),
           omega * Eigen::Vector3d( -.00112549, -.00055338, .000066195 ),
           1e-22 );
    close( differentiateIers2010Rp03AtmosphericGeocentreTranslationWrtUt1( 0 ),
           omega * Eigen::Vector3d( .00041173, -.00008174, .000001709 ),
           1e-22 );
}
