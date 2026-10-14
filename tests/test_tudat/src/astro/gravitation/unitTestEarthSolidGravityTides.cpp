/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <fstream>
#include <iostream>

#include "tudat/astro/earth_orientation/polarMotionCalculator.h"
#include "tudat/astro/gravitation/earthSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{
namespace
{
constexpr double sourceGM = 3.986004415e14, sourceRadius = 6378136.3;
constexpr double fieldGM = 3.986004415e14, fieldRadius = 6378136.46;
constexpr double sunGM = 1.32712442076e20, moonGM = 0.49028010560e13;
using basic_astrodynamics::TideSystem;

std::shared_ptr< gravitation::EarthSolidBodyTideGravityFieldVariations > makeSolid( const Eigen::Vector3d& sun,
                                                                                    const Eigen::Vector3d& moon,
                                                                                    const Eigen::Matrix3d& rotation,
                                                                                    const Eigen::Vector6d& arguments,
                                                                                    const TideSystem tideSystem = TideSystem::tideFree,
                                                                                    const double outputGM = fieldGM,
                                                                                    const double outputRadius = fieldRadius )
{
    const auto state = []( const Eigen::Vector3d& position ) {
        Eigen::Vector6d result = Eigen::Vector6d::Zero( );
        result.head< 3 >( ) = position;
        return result;
    };
    return std::make_shared< gravitation::EarthSolidBodyTideGravityFieldVariations >(
            []( double ) -> Eigen::Vector6d { return Eigen::Vector6d::Zero( ); },
            [ rotation ]( double ) { return Eigen::Quaterniond( rotation ); },
            std::vector< std::function< Eigen::Vector6d( double ) > >{ [ sun, state ]( double ) { return state( sun ); },
                                                                       [ moon, state ]( double ) { return state( moon ); } },
            sourceRadius,
            []( ) { return sourceGM; },
            std::vector< std::function< double( ) > >{ []( ) { return sunGM; }, []( ) { return moonGM; } },
            std::vector< std::string >{ "Sun", "Moon" },
            [ arguments ]( double ) { return arguments; },
            outputRadius,
            outputGM,
            tideSystem );
}

Eigen::Vector3d acceleration( gravitation::GravityFieldVariations& model,
                              const double time,
                              const Eigen::Vector3d& position,
                              const Eigen::Matrix3d& rotation,
                              const double gm = fieldGM,
                              const double radius = fieldRadius )
{
    Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( 5, 5 ), sine = cosine;
    model.addSphericalHarmonicsCorrections( time, sine, cosine );
    basic_mathematics::SphericalHarmonicsCache cache( 5, 5 );
    std::map< std::pair< int, int >, Eigen::Vector3d > terms;
    const Eigen::Vector3d fixed = gravitation::computeGeodesyNormalizedGravitationalAccelerationSum(
            rotation * position, gm, radius, cosine, sine, cache, terms );
    return rotation.transpose( ) * fixed;
}
}  // namespace

BOOST_AUTO_TEST_SUITE( earth_solid_gravity_tides )

BOOST_AUTO_TEST_CASE( costg_published_solid_and_pole_forces )
{
    const std::string source = __FILE__;
    std::ifstream data( source.substr( 0, source.find_last_of( "/\\" ) ) +
                        "/../../../reference_files/earth_gravity_solid/costg_samples.txt" );
    BOOST_REQUIRE( data.good( ) );
    std::string header;
    std::getline( data, header );
    double largestSolidError = 0.0, largestPoleError = 0.0;
    int count = 0, index;
    while( data >> index )
    {
        double mjdGps, tdb, tt;
        data >> mjdGps >> tdb >> tt;
        Eigen::Matrix3d rotation;
        Eigen::Vector3d position, sun, moon, expectedSolid, expectedPole;
        Eigen::Vector6d arguments;
        Eigen::Vector2d pole;
        for( int i = 0; i < 3; ++i )
            for( int j = 0; j < 3; ++j ) data >> rotation( i, j );
        for( int i = 0; i < 3; ++i ) data >> position[ i ];
        for( int i = 0; i < 6; ++i ) data >> arguments[ i ];
        data >> pole.x( ) >> pole.y( );
        for( auto* vector : { &sun, &moon, &expectedSolid, &expectedPole } )
            for( int i = 0; i < 3; ++i ) data >> ( *vector )[ i ];
        BOOST_REQUIRE( data.good( ) );
        auto solid = makeSolid( sun, moon, rotation, arguments );
        const Eigen::Vector2d secularArcsec = earth_orientation::getSecularPolePositionInMas( tt ) / 1000.0;
        const Eigen::Vector2d poleArcsec = pole * ( 180.0 * 3600.0 / mathematical_constants::PI );
        const Eigen::Vector2d wobble( poleArcsec.x( ) - secularArcsec.x( ), -( poleArcsec.y( ) - secularArcsec.y( ) ) );
        gravitation::SolidPoleTideGravityFieldVariations poleModel(
                [ wobble ]( double ) { return wobble; }, sourceGM, sourceRadius, fieldGM, fieldRadius );
        const double solidError = ( acceleration( *solid, tdb, position, rotation ) - expectedSolid ).norm( );
        const double poleError = ( acceleration( poleModel, tdb, position, rotation ) - expectedPole ).norm( );
        BOOST_CHECK_SMALL( solidError, 1e-11 );
        BOOST_CHECK_SMALL( poleError, 1e-11 );
        largestSolidError = std::max( largestSolidError, solidError );
        largestPoleError = std::max( largestPoleError, poleError );
        ++count;
    }
    BOOST_CHECK_EQUAL( count, 4 );
    std::cout << "COST-G max solid/pole force errors " << largestSolidError << "/" << largestPoleError << " m/s^2\n";
}

BOOST_AUTO_TEST_CASE( permanent_tide_and_output_rescaling )
{
    const Eigen::Vector3d sun( 1.4e11, 3e10, -2e10 ), moon( 3e8, -1.2e8, 1.8e8 );
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd( 0.7, Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );
    const Eigen::Vector6d arguments = ( Eigen::Vector6d( ) << 0.3, -0.4, 0.9, 0.6, -0.8, 0.1 ).finished( );
    auto tideFree = makeSolid( sun, moon, rotation, arguments );
    auto zeroTide = makeSolid( sun, moon, rotation, arguments, TideSystem::zeroTide );
    const auto free = tideFree->calculateSphericalHarmonicsCorrections( 0.0 );
    const auto zero = zeroTide->calculateSphericalHarmonicsCorrections( 0.0 );
    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero( 3, 4 );
    expected( 0, 0 ) = -4.4228e-8 * ( -0.31460 ) * 0.30190 * std::pow( sourceRadius / fieldRadius, 2 );
    BOOST_CHECK_SMALL( ( zero.first - free.first - expected ).norm( ), 1e-23 );
    BOOST_CHECK_SMALL( ( zero.second - free.second ).norm( ), 1e-23 );
    auto scaled = makeSolid( sun, moon, rotation, arguments, TideSystem::tideFree, fieldGM * 1.7, fieldRadius * 1.2 );
    const auto coefficients = scaled->calculateSphericalHarmonicsCorrections( 0.0 );
    for( int n = 2; n <= 4; ++n )
    {
        const double scale = 1.0 / ( 1.7 * std::pow( 1.2, n ) );
        BOOST_CHECK_SMALL( ( coefficients.first.row( n - 2 ) - scale * free.first.row( n - 2 ) ).norm( ), 1e-22 );
        BOOST_CHECK_SMALL( ( coefficients.second.row( n - 2 ) - scale * free.second.row( n - 2 ) ).norm( ), 1e-22 );
    }
    BOOST_CHECK_SMALL( ( acceleration( *scaled, 0.0, Eigen::Vector3d( 7e6, 1e6, -2e6 ), rotation, fieldGM * 1.7, fieldRadius * 1.2 ) -
                         acceleration( *tideFree, 0.0, Eigen::Vector3d( 7e6, 1e6, -2e6 ), rotation ) )
                               .norm( ),
                       1e-21 );
    BOOST_CHECK_THROW( makeSolid( sun, moon, rotation, arguments, TideSystem::unknown ), std::invalid_argument );
    BOOST_CHECK_THROW( makeSolid( sun, moon, rotation, arguments, TideSystem::meanTide ), std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( degree_two_forces_degree_four_and_pole_block_signs )
{
    const double sunDistance = 1.5e11, moonDistance = 4e8;
    auto solid = makeSolid( Eigen::Vector3d( sunDistance, 0, 0 ),
                            Eigen::Vector3d( moonDistance, 0, 0 ),
                            Eigen::Matrix3d::Identity( ),
                            Eigen::Vector6d::Zero( ),
                            TideSystem::tideFree,
                            sourceGM,
                            sourceRadius );
    const auto coefficients = solid->calculateSphericalHarmonicsCorrections( 0.0 );
    const double forcing2 = ( sunGM * std::pow( sourceRadius / sunDistance, 3 ) + moonGM * std::pow( sourceRadius / moonDistance, 3 ) ) /
            ( 5.0 * sourceGM );
    BOOST_CHECK_SMALL( coefficients.first( 2, 0 ) - ( -0.00089 ) * ( -std::sqrt( 5.0 ) / 2 ) * forcing2, 1e-25 );
    BOOST_CHECK_SMALL( coefficients.first( 2, 2 ) - ( -0.00057 ) * ( std::sqrt( 15.0 ) / 2 ) * forcing2, 1e-25 );
    BOOST_CHECK_EQUAL( coefficients.first( 2, 3 ), 0.0 );
    BOOST_CHECK_EQUAL( coefficients.second.col( 0 ).norm( ), 0.0 );
    const double forcing3 = ( sunGM * std::pow( sourceRadius / sunDistance, 4 ) + moonGM * std::pow( sourceRadius / moonDistance, 4 ) ) /
            ( 7.0 * sourceGM );
    BOOST_CHECK_SMALL( coefficients.first( 1, 3 ) - 0.094 * std::sqrt( 70.0 ) / 4 * forcing3, 1e-25 );
    gravitation::SolidPoleTideGravityFieldVariations pole(
            []( double ) { return Eigen::Vector2d( 0.2, -0.3 ); }, sourceGM, sourceRadius, sourceGM, sourceRadius );
    const auto poleCoefficients = pole.calculateSphericalHarmonicsCorrections( 0.0 );
    BOOST_CHECK_EQUAL( poleCoefficients.first.rows( ), 1 );
    BOOST_CHECK_EQUAL( poleCoefficients.first.cols( ), 1 );
    BOOST_CHECK_SMALL( poleCoefficients.first( 0, 0 ) - ( -1.333e-9 * ( 0.2 - 0.0115 * 0.3 ) ), 1e-25 );
    BOOST_CHECK_SMALL( poleCoefficients.second( 0, 0 ) - ( -1.333e-9 * ( -0.3 - 0.0115 * 0.2 ) ), 1e-25 );
}

BOOST_AUTO_TEST_SUITE_END( )
}  // namespace unit_tests
}  // namespace tudat
