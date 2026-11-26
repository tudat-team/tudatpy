#ifndef TUDAT_PARALLEL_FOR_H
#define TUDAT_PARALLEL_FOR_H

#include <exception>

namespace tudat
{
namespace utilities
{
//! Validate a local execution policy, including availability of the optional backend.
void validateNumberOfThreads( const int numberOfThreads );

//! True only when a new inner team may be created. Does not detect outer std::threads.
bool useParallelTeam( const int numberOfThreads );

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
