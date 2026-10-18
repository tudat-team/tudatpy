#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/ground_stations/atmosphericTideLoading.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "pr1092/reference_files/atmospheric_tides/origin_translation_reference.h"
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
}  // namespace

//! Check S1 and S2 signs and periods with exact quarter-day arithmetic.
BOOST_AUTO_TEST_CASE( directQuarterDaySignsAndS2HalfDayPeriod )
{
    Eigen::Matrix< double, 3, 4 > coefficients;
    coefficients << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;

    // The combined harmonics reproduce these direct quarter-day vectors.
    close( synthesizeAtmosphericTideLoading( coefficients, -43200 ), Eigen::Vector3d( 4, 12, 20 ) );
    close( synthesizeAtmosphericTideLoading( coefficients, -21600 ), Eigen::Vector3d( -1, -1, -1 ) );
    close( synthesizeAtmosphericTideLoading( coefficients, 0 ), Eigen::Vector3d( 2, 2, 2 ) );
    close( synthesizeAtmosphericTideLoading( coefficients, 21600 ), Eigen::Vector3d( -5, -13, -21 ) );
    auto s1 = Eigen::Matrix< double, 3, 4 >::Zero( ).eval( );
    s1.col( 0 ) << 1, 2, 3;
    s1.col( 1 ) << 4, 5, 6;
    close( synthesizeAtmosphericTideLoading( s1, -43200 ), Eigen::Vector3d( 1, 2, 3 ) );
    close( synthesizeAtmosphericTideLoading( s1, -21600 ), Eigen::Vector3d( 4, 5, 6 ) );
    close( synthesizeAtmosphericTideLoading( s1, 0 ), Eigen::Vector3d( -1, -2, -3 ) );
    close( synthesizeAtmosphericTideLoading( s1, 21600 ), Eigen::Vector3d( -4, -5, -6 ) );
    auto s2 = Eigen::Matrix< double, 3, 4 >::Zero( ).eval( );
    s2.col( 2 ) << 1, 2, 3;
    s2.col( 3 ) << 4, 5, 6;
    close( synthesizeAtmosphericTideLoading( s2, -43200 ), Eigen::Vector3d( 1, 2, 3 ) );
    close( synthesizeAtmosphericTideLoading( s2, -32400 ), Eigen::Vector3d( 4, 5, 6 ) );
    close( synthesizeAtmosphericTideLoading( s2, -21600 ), Eigen::Vector3d( -1, -2, -3 ) );
    close( synthesizeAtmosphericTideLoading( s2, -10800 ), Eigen::Vector3d( -4, -5, -6 ) );
    close( synthesizeAtmosphericTideLoading( s2, 0 ), Eigen::Vector3d( 1, 2, 3 ) );
    close( synthesizeAtmosphericTideLoading( Eigen::Matrix< double, 3, 4 >::Zero( ), 0 ), Eigen::Vector3d::Zero( ) );
}

//! Compare harmonic rates with direct derivatives over positive and negative epochs.
BOOST_AUTO_TEST_CASE( analyticDerivativeAndEpochReduction )
{
    Eigen::Matrix< double, 3, 4 > coefficients;
    coefficients << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
    const double omega = 2 * mathematical_constants::PI / 86400;

    // The quarter-day derivatives and whole-day reductions preserve their signs and phase.
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

//! Interpolate the published atmospheric grid against independent Decimal arithmetic.
BOOST_AUTO_TEST_CASE( realGridNodesAndIndependentDecimalInterpolation )
{
    const std::string base = paths::getTudatTestDataPath( ) + "pr1092/reference_files/atmospheric_tides/";
    const auto raw = readAtmosphericTideGrid( base + "rp03_ce_grid_extract.txt" );

    // At a grid node, interpolation returns the published coefficient matrix.
    BOOST_CHECK_SMALL( ( interpolateAtmosphericTideGrid( raw, -19, 47 ) - raw.coefficientsMetres.front( ) ).norm( ), 1e-18 );

    // Independent Decimal interpolation of four original provider rows; this is not reference-software conformance.
    Eigen::Matrix< double, 3, 4 > expected;
    expected << -.229413396401398, -.140409687704627, 1.040762193571904, .276668845291711, .059365717547470, .059123892373735,
            -.033761020617976, -.014395543429494, .007941624058988, -.103894064696771, -.016544820023036, .059869142014747;
    BOOST_CHECK_SMALL( ( interpolateAtmosphericTideGrid( raw, -19.018307, 47.22921 ) - expected * 1e-3 ).norm( ), 1e-16 );
}

//! Synthesize radial loading directly from published polar-node coefficients.
BOOST_AUTO_TEST_CASE( publishedGridNodeRadialQuarterDayValues )
{
    const std::string base = paths::getTudatTestDataPath( ) + "pr1092/reference_files/atmospheric_tides/";
    const auto grid = readAtmosphericTideGrid( base + "rp03_ce_polar_node_extract.txt" );

    // Direct arithmetic on the published Figure24.17 four radial coefficients, not a self-generated model fixture.
    const auto& coefficients = grid.coefficientsMetres.front( );

    // The four radial coefficients give these quarter-day displacements in metres.
    BOOST_CHECK_SMALL( synthesizeAtmosphericTideLoading( coefficients, -43200 ).x( ) + .0000604, 1e-18 );
    BOOST_CHECK_SMALL( synthesizeAtmosphericTideLoading( coefficients, -21600 ).x( ) - .0001548, 1e-18 );
    BOOST_CHECK_SMALL( synthesizeAtmosphericTideLoading( coefficients, 0 ).x( ) + .0000970, 1e-18 );
    BOOST_CHECK_SMALL( synthesizeAtmosphericTideLoading( coefficients, 21600 ).x( ) - .0000026, 1e-18 );
}

//! Compare the IERS origin translation with independent Table 7.6 reference values.
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
}

//! Apply the CE-to-CM translation before projecting into local station axes.
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

        // Project the translated vector once, and recover CE by subtracting the same translation.
        close( radialNorthEastToCartesian.transpose( ) * cmCartesian, ceLocal + independentlyProjectedTranslation, 1e-17 );
        close( cmCartesian - d, ceCartesian, 1e-17 );
    }
}

//! Compare origin-translation velocities with independent Fourier reference values.
BOOST_AUTO_TEST_CASE( independentRp03TranslationRatesAndQuarterDaySigns )
{
    for( const auto& reference : atmosphericOriginReferences )
    {

        // Each reference row supplies an independent Cartesian translation derivative.
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
