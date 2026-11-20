#include <cmath>
#include <iostream>
#include <stdexcept>
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"

#ifndef TUDAT_BUILD_WITH_OPENMP
#error "Installed target does not export its OpenMP feature definition"
#endif
static_assert( ( TUDAT_BUILD_WITH_OPENMP != 0 ) == EXPECT_OPENMP, "Installed config and target disagree on OpenMP availability" );

int main( )
{
    using namespace tudat;
#if TUDAT_BUILD_WITH_OPENMP
    constexpr int threads = 4;
#ifndef _OPENMP
#error "The exported OpenMP compile requirements did not reach this consumer"
#endif
#else
    constexpr int threads = 1;
#endif
    int actualThreads = 1;
#if TUDAT_BUILD_WITH_OPENMP
    utilities::parallelFor( 20, threads, 1, [ & ]( int i ) {
        if( i == 0 ) actualThreads = omp_get_num_threads( );
    } );
#endif
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero( 21, 21 ), s = c;
    c( 0, 0 ) = 1.;
    c( 2, 0 ) = -.001;
    c( 20, 5 ) = 1.e-8;
    Eigen::Vector3d position( 7.e6, 2.e6, 1.e6 );
    auto getPosition = [ & ]( Eigen::Vector3d& out ) { out = position; };
    gravitation::SphericalHarmonicsGravitationalAccelerationModel serial( getPosition, 3.986e14, 6.378e6, c, s );
    auto parallel = serial;
    parallel.setNumberOfThreads( threads );
    for( int step = 0; step < 20; ++step )
    {
        position[ 1 ] += 100.;
        serial.updateMembers( );
        parallel.updateMembers( );
        if( !parallel.getAcceleration( ).allFinite( ) || ( serial.getAcceleration( ) - parallel.getAcceleration( ) ).norm( ) != 0. )
            throw std::runtime_error( "Installed model differs from sequential reference" );
    }
#if !TUDAT_BUILD_WITH_OPENMP
    bool rejected = false;
    try
    {
        parallel.setNumberOfThreads( 2 );
    }
    catch( const std::runtime_error& )
    {
        rejected = true;
    }
    if( !rejected ) throw std::runtime_error( "Unavailable backend was not rejected" );
#endif
    std::cout << "Installed gravity model passed; requested threads=" << threads << "; observed helper team=" << actualThreads << '\n';
}
