#ifndef TUDAT_PARALLEL_FOR_H
#define TUDAT_PARALLEL_FOR_H

#include <exception>
#include <stdexcept>
#if defined( TUDAT_BUILD_WITH_OPENMP ) && TUDAT_BUILD_WITH_OPENMP
#include <omp.h>
#endif

namespace tudat
{
namespace utilities
{
//! Validate a local execution policy, including availability of the optional backend.
inline void validateNumberOfThreads( const int numberOfThreads )
{
    if( numberOfThreads < 1 ) throw std::invalid_argument( "Number of threads must be positive" );
#if !defined( TUDAT_BUILD_WITH_OPENMP ) || !TUDAT_BUILD_WITH_OPENMP
    if( numberOfThreads > 1 ) throw std::runtime_error( "Multiple gravity threads require TUDAT_BUILD_WITH_OPENMP=ON" );
#endif
}

//! True only when a new inner team may be created. Does not detect outer std::threads.
inline bool useParallelTeam( const int numberOfThreads )
{
#if defined( TUDAT_BUILD_WITH_OPENMP ) && TUDAT_BUILD_WITH_OPENMP
    return numberOfThreads > 1 && !omp_in_parallel( );
#else
    return false;
#endif
}

//! Execute independent indices; join all workers before rethrowing an exception on the caller.
//! The callable must not invoke user callbacks or mutate another index's output.
template< typename Function >
void parallelFor( const int end, const int numberOfThreads, const int chunk, const Function& function )
{
#if defined( TUDAT_BUILD_WITH_OPENMP ) && TUDAT_BUILD_WITH_OPENMP
    if( useParallelTeam( numberOfThreads ) )
    {
        std::exception_ptr failure;
#pragma omp parallel for num_threads( numberOfThreads ) schedule( static, chunk )
        for( int i = 0; i < end; ++i )
        {
            try
            {
                function( i );
            }
            catch( ... )
            {
#pragma omp critical( tudat_parallel_failure )
                {
                    if( !failure ) failure = std::current_exception( );
                }
            }
        }
        if( failure ) std::rethrow_exception( failure );
        return;
    }
#endif
    for( int i = 0; i < end; ++i ) function( i );
}
}  // namespace utilities
}  // namespace tudat
#endif
