/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/gravitation/oceanPoleTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include <array>
#include <fstream>
#include <sstream>
using namespace tudat;
using namespace tudat::gravitation;
namespace
{
std::string fixture( const std::string& name )
{
    const std::string file = __FILE__;
    return file.substr( 0, file.find_last_of( "/" ) ) + "/../../../reference_files/ocean_pole_gravity/" + name;
}
OceanPoleGravityConstants constants( )
{
    return { 7.292115e-5, 6.673e-11, 1025.0, 9.7803278, 3.986004418e14, 6378136.6, { 0.687, 0.0036 } };
}
std::vector< std::vector< double > > rows( const std::string& name )
{
    std::ifstream stream( fixture( name ) );
    if( !stream ) throw std::runtime_error( "Missing gravity fixture" );
    std::vector< std::vector< double > > result;
    std::string line;
    while( std::getline( stream, line ) )
    {
        std::istringstream input( line );
        std::vector< double > row;
        double value;
        while( input >> value ) row.push_back( value );
        if( row.size( ) >= 4 && row[ 0 ] > 50000 ) result.push_back( row );
    }
    return result;
}
double printedRounding( double value )
{
    return value == 0 ? 0 : 0.5 * std::pow( 10.0, std::floor( std::log10( std::abs( value ) ) ) - 6 );
}
}  // namespace
BOOST_AUTO_TEST_CASE( all366PublishedGravityReferenceEpochs )
{
    const auto data = input_output::readOceanPoleGravityData( fixture( "desaiscopolecoef_degree180.txt" ) );
    const auto reference = rows( "desaiscopole2004.txt" );
    BOOST_REQUIRE_EQUAL( reference.size( ), 366 );
    Eigen::Vector2d wobble;
    const auto c = constants( );
    const double radians = std::acos( -1.0 ) / ( 180 * 3600 );
    OceanPoleTideGravityFieldVariations model(
            data, [ & ]( double ) { return wobble; }, { 0, 0, -0.3075 }, c, c.gravitationalParameter, c.referenceRadius, 2, 0, 2, 2 );
    // Error budget: half one printed final digit in each 7-significant-digit wobble/output,
    // propagated through the published R2/A21/B21. Remaining printed constants contribute <1e-20.
    const double r2 = 2.6871735998e-4;
    const std::array< double, 2 > cr = { r2 * radians * ( -0.24325330521 * 0.687 - 0.0054680741193 * 0.0036 ),
                                         r2 * radians * ( -0.24325330521 * 0.0036 + 0.0054680741193 * 0.687 ) };
    const std::array< double, 2 > sr = { r2 * radians * ( 0.0054680741193 * 0.687 + 0.19252111185 * 0.0036 ),
                                         r2 * radians * ( 0.0054680741193 * 0.0036 - 0.19252111185 * 0.687 ) };
    for( const auto& row : reference )
    {
        wobble = Eigen::Vector2d( row[ 5 ], row[ 6 ] ) * radians;
        const auto correction = model.calculateSphericalHarmonicsCorrections( ( row[ 0 ] - 51544.5 ) * 86400 );
        const double dc = printedRounding( row[ 7 ] ) + std::abs( cr[ 0 ] ) * printedRounding( row[ 5 ] ) +
                std::abs( cr[ 1 ] ) * printedRounding( row[ 6 ] ) + 1e-20;
        const double ds = printedRounding( row[ 8 ] ) + std::abs( sr[ 0 ] ) * printedRounding( row[ 5 ] ) +
                std::abs( sr[ 1 ] ) * printedRounding( row[ 6 ] ) + 1e-20;
        BOOST_CHECK_SMALL( correction.first( 0, 1 ) - row[ 7 ], dc );
        BOOST_CHECK_SMALL( correction.second( 0, 1 ) - row[ 8 ], ds );
    }
}
BOOST_AUTO_TEST_CASE( independentCostgFullDegree180Forces )
{
    const auto data = input_output::readOceanPoleGravityData( fixture( "desaiscopolecoef_degree180.txt" ) );
    const auto position = rows( "costg_00orbit_itrf.txt" ), rotation = rows( "costg_01earthRotation_rotaryMatrix.txt" ),
               eop = rows( "costg_01earthRotation_interpolatedEOP.txt" ), force = rows( "costg_06oceanPoleTide_icrf.txt" );
    BOOST_REQUIRE_EQUAL( force.size( ), 4 );
    BOOST_REQUIRE_EQUAL( position.size( ), 4 );
    BOOST_REQUIRE_EQUAL( rotation.size( ), 4 );
    BOOST_REQUIRE_EQUAL( eop.size( ), 4 );
    std::ifstream loveStream( fixture( "loadLoveNumbers_Gegout97_degree180.txt" ) );
    std::vector< double > love;
    double value;
    while( loveStream >> value ) love.push_back( value );
    BOOST_REQUIRE_EQUAL( love.size( ), 181 );
    Eigen::Vector2d wobble;
    const auto c = constants( );
    const double radians = std::acos( -1.0 ) / ( 180 * 3600 );
    // COST-G nominal field's GM/R, distinct from IERS forcing constants.
    const double gm = 3.986004415e14, radius = 6378136.460;
    OceanPoleTideGravityFieldVariations model( data, [ & ]( double ) { return wobble; }, love, c, gm, radius, 2, 0, 180, 180 );
    basic_mathematics::SphericalHarmonicsCache cache( 180, 180 );
    std::map< std::pair< int, int >, Eigen::Vector3d > perTerm;
    for( std::size_t i = 0; i < force.size( ); ++i )
    {
        // Inject the benchmark's historical secular pole and its supplied polar motion.
        const double years = ( eop[ i ][ 0 ] - 51544.5 ) / 365.25;
        wobble = { eop[ i ][ 1 ] - ( 0.055 + 0.001677 * years ) * radians, -eop[ i ][ 2 ] + ( 0.3205 + 0.003460 * years ) * radians };
        const auto block = model.calculateSphericalHarmonicsCorrections( 0.0 );
        Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( 181, 181 ), sine = cosine;
        cosine.block( 2, 0, 179, 181 ) = block.first;
        sine.block( 2, 0, 179, 181 ) = block.second;
        Eigen::Matrix3d matrix;
        for( int j = 0; j < 9; ++j ) matrix( j / 3, j % 3 ) = rotation[ i ][ j + 1 ];
        const Eigen::Vector3d acceleration = computeGeodesyNormalizedGravitationalAccelerationSum(
                Eigen::Vector3d( position[ i ][ 1 ], position[ i ][ 2 ], position[ i ][ 3 ] ), gm, radius, cosine, sine, cache, perTerm );
        const Eigen::Vector3d inertial = matrix.transpose( ) * acceleration;
        const Eigen::Vector3d expected( force[ i ][ 1 ], force[ i ][ 2 ], force[ i ][ 3 ] );
        BOOST_TEST_MESSAGE( "COST-G MJD " << force[ i ][ 0 ] << ": force residual " << ( inertial - expected ).norm( ) << " m/s^2" );
        BOOST_CHECK_SMALL( ( inertial - expected ).norm( ), 1e-11 );
    }
}
BOOST_AUTO_TEST_CASE( selectionScalingMissingCoverageAndSnapshotOwnership )
{
    auto data = input_output::readOceanPoleGravityData( fixture( "desaiscopolecoef_degree180.txt" ) );
    auto c = constants( );
    Eigen::Vector2d wobble( 1e-6, -2e-6 );
    auto callback = [ & ]( double ) { return wobble; };
    std::vector< double > love( 11, -0.1 );
    OceanPoleTideGravityFieldVariations full( data, callback, love, c, c.gravitationalParameter, c.referenceRadius, 2, 0, 10, 10 );
    OceanPoleTideGravityFieldVariations scaled( data, callback, love, c, 2 * c.gravitationalParameter, 2 * c.referenceRadius, 3, 1, 10, 6 );
    const auto a = full.calculateSphericalHarmonicsCorrections( 0 ), b = scaled.calculateSphericalHarmonicsCorrections( 0 );
    for( int n = 3; n <= 10; ++n )
        for( int m = 1; m <= std::min( n, 6 ); ++m )
        {
            BOOST_CHECK_SMALL( b.first( n - 3, m - 1 ) - a.first( n - 2, m ) * std::pow( 0.5, n + 1 ), 1e-25 );
            BOOST_CHECK_SMALL( b.second( n - 3, m - 1 ) - a.second( n - 2, m ) * std::pow( 0.5, n + 1 ), 1e-25 );
        }
    data.coefficients.erase( { 3, 1 } );
    love[ 3 ] = 1e99;
    c.gamma2 = 0;
    const auto frozen = scaled.calculateSphericalHarmonicsCorrections( 0 );
    BOOST_CHECK_SMALL( ( frozen.first - b.first ).norm( ), 1e-30 );
    BOOST_CHECK_THROW(
            OceanPoleTideGravityFieldVariations( data, callback, std::vector< double >( 11, -0.1 ), constants( ), 1e14, 6e6, 2, 0, 10, 10 ),
            std::invalid_argument );
    BOOST_CHECK_THROW( OceanPoleTideGravityFieldVariations( data, callback, { 0, 0, -0.3 }, constants( ), 1e14, 6e6, 2, 0, 10, 10 ),
                       std::invalid_argument );
    wobble.x( ) = std::numeric_limits< double >::quiet_NaN( );
    BOOST_CHECK_THROW( full.calculateSphericalHarmonicsCorrections( 0 ), std::invalid_argument );
}
