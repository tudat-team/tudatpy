#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <atomic>
#include <future>
#include <sstream>
#include "tudat/math/integrators/rungeKutta4Integrator.h"
#include <cereal/archives/binary.hpp>
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#if defined( TUDAT_BUILD_WITH_OPENMP ) && TUDAT_BUILD_WITH_OPENMP
#include <omp.h>
#endif

using namespace tudat;
using namespace basic_mathematics;

#if defined( TUDAT_BUILD_WITH_OPENMP ) && TUDAT_BUILD_WITH_OPENMP
constexpr int workers = 4;
#else
constexpr int workers = 1;
#endif

BOOST_AUTO_TEST_CASE( policyValidation )
{
    LegendreCache cache( 20, 20 );
    BOOST_CHECK_EQUAL( cache.getNumberOfThreads( ), 1 );
    BOOST_CHECK_THROW( cache.setNumberOfThreads( 0 ), std::invalid_argument );
    BOOST_CHECK_THROW( simulation_setup::sphericalHarmonicAcceleration( 20, 20, false, -1 ), std::invalid_argument );
#if !defined( TUDAT_BUILD_WITH_OPENMP ) || !TUDAT_BUILD_WITH_OPENMP
    BOOST_CHECK_THROW( cache.setNumberOfThreads( 2 ), std::runtime_error );
    BOOST_CHECK_THROW( simulation_setup::sphericalHarmonicAcceleration( 20, 20, false, 2 ), std::runtime_error );
#endif
    auto settings = std::dynamic_pointer_cast< simulation_setup::SphericalHarmonicAccelerationSettings >(
            simulation_setup::sphericalHarmonicAcceleration( 20, 20, true, workers ) );
    BOOST_CHECK_EQUAL( settings->numberOfThreads_, workers );
    BOOST_CHECK( settings->removePointMass_ );
}

void compareCaches( const int count )
{
    for( bool normalized : { false, true } )
        for( int degree : { 1, 20, 100 } )
            for( int order : { 0, 1, 5, 100 } )
            {
                if( order > degree ) continue;
                for( bool second : { false, true } )
                {
                    LegendreCache serial( degree, order, normalized ), parallel( degree, order, normalized );
                    serial.setComputeFirstDerivatives( true );
                    parallel.setComputeFirstDerivatives( true );
                    serial.setComputeSecondDerivatives( second );
                    parallel.setComputeSecondDerivatives( second );
                    parallel.setNumberOfThreads( count );
                    for( bool consistency : { false, true } )
                        for( double x : { -.99999, -.3, 0., .5, .99999, .5 } )
                        {
                            serial.update( x, consistency );
                            parallel.update( x, consistency );
                            // Same recurrence and operation order: exact parity is expected on this build.
                            for( int n = 0; n <= degree; ++n )
                                for( int m = 0; m <= std::min( n, order ); ++m )
                                {
                                    BOOST_CHECK_EQUAL( serial.getLegendrePolynomial( n, m ), parallel.getLegendrePolynomial( n, m ) );
                                    BOOST_CHECK_EQUAL( serial.getLegendrePolynomialDerivative( n, m ),
                                                       parallel.getLegendrePolynomialDerivative( n, m ) );
                                    if( second )
                                        BOOST_CHECK_EQUAL( serial.getLegendrePolynomialSecondDerivative( n, m ),
                                                           parallel.getLegendrePolynomialSecondDerivative( n, m ) );
                                }
                        }
                    auto copy = parallel;
                    copy.resetMaximumDegreeAndOrder( 3, 1 );
                    copy.update( -.2 );
                    BOOST_CHECK_EQUAL( parallel.getMaximumDegree( ), degree );
                    BOOST_CHECK_EQUAL( copy.getNumberOfThreads( ), count );
                }
            }
}

BOOST_AUTO_TEST_CASE( cacheParityAndCopies )
{
    compareCaches( workers );
}

BOOST_AUTO_TEST_CASE( joinedExceptionAndRetry )
{
    std::atomic< int > completed{ 0 };
    BOOST_CHECK_THROW( utilities::parallelFor( 20,
                                               workers,
                                               1,
                                               [ & ]( int i ) {
                                                   ++completed;
                                                   if( i == 7 ) throw std::runtime_error( "injected worker failure" );
                                               } ),
                       std::runtime_error );
    // Parallel workers finish even after one fails; serial execution stops at the failing index.
    BOOST_CHECK_EQUAL( completed.load( ), workers > 1 ? 20 : 8 );
    LegendreCache cache( 20, 20, true );
    cache.setNumberOfThreads( workers );
    cache.setComputeFirstDerivatives( true );
    cache.setthrowExceptionForSingularDerivative( true );
    BOOST_CHECK_THROW( cache.update( 1.0 ), std::runtime_error );
    BOOST_CHECK( std::isnan( cache.getCurrentPolynomialParameter( ) ) );
    BOOST_CHECK_THROW( cache.update( 1.0 ), std::runtime_error );
    cache.setthrowExceptionForSingularDerivative( false );
    BOOST_CHECK_NO_THROW( cache.update( 1.0 ) );
    BOOST_CHECK_EQUAL( cache.getCurrentPolynomialParameter( ), 1.0 );
    cache.update( .5 );
    LegendreCache reference( 20, 20, true );
    reference.setComputeFirstDerivatives( true );
    reference.update( .5 );
    BOOST_CHECK_EQUAL( cache.getLegendrePolynomialDerivative( 15, 5 ), reference.getLegendrePolynomialDerivative( 15, 5 ) );
}

BOOST_AUTO_TEST_CASE( modelBuffersAndFallbacks )
{
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero( 21, 21 ), s = c;
    c( 0, 0 ) = 1.;
    c( 2, 0 ) = -.001;
    c( 20, 5 ) = 1.e-8;
    Eigen::Vector3d position( 7.e6, 2.e6, 1.e6 );
    const auto caller = std::this_thread::get_id( );
    auto getPosition = [ & ]( Eigen::Vector3d& out ) {
        BOOST_CHECK( std::this_thread::get_id( ) == caller );
        out = position;
    };
    gravitation::SphericalHarmonicsGravitationalAccelerationModel serial( getPosition, 3.986e14, 6.378e6, c, s );
    auto parallel = serial;
    parallel.setNumberOfThreads( workers );
    BOOST_CHECK_EQUAL( serial.getNumberOfThreads( ), 1 );
    for( auto p : { Eigen::Vector3d( 7.e6, 2.e6, 1.e6 ), Eigen::Vector3d( 7., 0., 7.e6 ), Eigen::Vector3d( 7., 0., -7.e6 ) } )
    {
        position = p;
        serial.updateMembers( );
        parallel.updateMembers( );
        BOOST_CHECK( serial.getAcceleration( ).allFinite( ) );
        BOOST_CHECK_EQUAL( ( serial.getAcceleration( ) - parallel.getAcceleration( ) ).norm( ), 0. );
    }
    parallel.setSaveSphericalHarmonicTermsSeparately( true );
    position = Eigen::Vector3d( 7.e6, 2.e6, 1.e6 );
    serial.updateMembers( );
    parallel.updateMembers( );
    BOOST_CHECK_EQUAL( ( serial.getAcceleration( ) - parallel.getAcceleration( ) ).norm( ), 0. );
}

BOOST_AUTO_TEST_CASE( independentConcurrentCachesAndNestedTeam )
{
    // Separate caches on std::threads share neither policy nor mutable scratch.
    auto evaluate = []( int count ) {
        LegendreCache cache( 100, 100, true );
        cache.setNumberOfThreads( count );
        cache.setComputeFirstDerivatives( true );
        for( int i = 0; i < 20; ++i ) cache.update( .1 + i * .01 );
        return cache.getLegendrePolynomialDerivative( 90, 10 );
    };
    auto a = std::async( std::launch::async, evaluate, workers );
    auto b = std::async( std::launch::async, evaluate, 1 );
    BOOST_CHECK_EQUAL( a.get( ), b.get( ) );
#if defined( TUDAT_BUILD_WITH_OPENMP ) && TUDAT_BUILD_WITH_OPENMP
    std::atomic< int > wrongLevel{ 0 };
#pragma omp parallel num_threads( 2 )
    {
        utilities::parallelFor( 10, 4, 1, [ & ]( int ) {
            if( omp_get_level( ) != 1 ) ++wrongLevel;
        } );
    }
    BOOST_CHECK_EQUAL( wrongLevel.load( ), 0 );
#endif
}

BOOST_AUTO_TEST_CASE( executionPolicyExcludedFromArchive )
{
    simulation_setup::SphericalHarmonicAccelerationSettings source( 100, 50, true, workers );
    std::stringstream bytes;
    {
        cereal::BinaryOutputArchive archive( bytes );
        archive( source );
    }
    simulation_setup::SphericalHarmonicAccelerationSettings restored( 1, 1, false, workers );
    {
        cereal::BinaryInputArchive archive( bytes );
        archive( restored );
    }
    BOOST_CHECK_EQUAL( restored.maximumDegree_, 100 );
    BOOST_CHECK_EQUAL( restored.maximumOrder_, 50 );
    BOOST_CHECK( restored.removePointMass_ );
    BOOST_CHECK_EQUAL( restored.numberOfThreads_, 1 );
    std::stringstream serialBytes;
    source.numberOfThreads_ = 1;
    {
        cereal::BinaryOutputArchive archive( serialBytes );
        archive( source );
    }
    BOOST_CHECK_EQUAL( bytes.str( ), serialBytes.str( ) );
}

BOOST_AUTO_TEST_CASE( reusableTermStorageResizesAndOverwrites )
{
    std::vector< Eigen::Vector3d > terms;
    for( int degree : { 3, 20, 3 } )
    {
        Eigen::MatrixXd c = Eigen::MatrixXd::Zero( degree + 1, degree + 1 ), s = c;
        c( 0, 0 ) = 1.;
        c( 2, 0 ) = -.001;
        SphericalHarmonicsCache serialCache( degree + 1, degree + 1 ), parallelCache = serialCache;
        parallelCache.getLegendreCache( ).setNumberOfThreads( workers );
        std::map< std::pair< int, int >, Eigen::Vector3d > dummy;
        for( const Eigen::Vector3d position : { Eigen::Vector3d( 7.e6, 2.e6, 1.e6 ), Eigen::Vector3d( 6.e6, -3.e6, 2.e6 ) } )
        {
            for( auto& term : terms ) term.setConstant( TUDAT_NAN );
            const auto reference = gravitation::computeGeodesyNormalizedGravitationalAccelerationSum(
                    position, 3.986e14, 6.378e6, c, s, serialCache, dummy );
            const auto actual = gravitation::computeGeodesyNormalizedGravitationalAccelerationSum( position,
                                                                                                   3.986e14,
                                                                                                   6.378e6,
                                                                                                   c,
                                                                                                   s,
                                                                                                   parallelCache,
                                                                                                   dummy,
                                                                                                   false,
                                                                                                   Eigen::Matrix3d::Identity( ),
                                                                                                   true,
                                                                                                   1.e-6,
                                                                                                   workers,
                                                                                                   &terms );
            BOOST_CHECK( actual.allFinite( ) );
            BOOST_CHECK_EQUAL( ( reference - actual ).norm( ), 0. );
            if( workers > 1 ) BOOST_CHECK_EQUAL( terms.size( ), ( degree + 1 ) * ( degree + 1 ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( nativePropagationParity )
{
    const auto propagate = []( const int count ) {
        Eigen::MatrixXd c = Eigen::MatrixXd::Zero( 21, 21 ), s = c;
        c( 0, 0 ) = 1.;
        c( 2, 0 ) = -.001;
        c( 20, 5 ) = 1.e-8;
        Eigen::Vector3d position = Eigen::Vector3d::Zero( );
        gravitation::SphericalHarmonicsGravitationalAccelerationModel model(
                [ & ]( Eigen::Vector3d& out ) { out = position; }, 3.986e14, 6.378e6, c, s );
        model.setNumberOfThreads( count );
        const auto derivative = [ & ]( double, const Eigen::VectorXd& state ) -> Eigen::VectorXd {
            position = state.head< 3 >( );
            model.updateMembers( );
            Eigen::VectorXd result( 6 );
            result.head< 3 >( ) = state.tail< 3 >( );
            result.tail< 3 >( ) = model.getAcceleration( );
            return result;
        };
        Eigen::VectorXd initial( 6 );
        initial << 7.e6, 0., 0., 0., 6000., 4500.;
        numerical_integrators::RungeKutta4Integrator<> integrator( derivative, 0., initial, 30. );
        std::vector< Eigen::VectorXd > history;
        for( int step = 0; step < 500; ++step ) history.push_back( integrator.performIntegrationStep( 30. ) );
        return history;
    };
    const auto serial = propagate( 1 ), parallel = propagate( workers );
    BOOST_REQUIRE_EQUAL( serial.size( ), parallel.size( ) );
    for( std::size_t i = 0; i < serial.size( ); ++i )
    {
        BOOST_CHECK( parallel[ i ].allFinite( ) );
        BOOST_CHECK_EQUAL( ( serial[ i ] - parallel[ i ] ).norm( ), 0. );
    }
}
