#include <stdexcept>
#if defined( TUDAT_BUILD_WITH_OPENMP ) && TUDAT_BUILD_WITH_OPENMP
#include <omp.h>
#endif

#include "tudat/basics/parallelFor.h"

namespace tudat
{
namespace utilities
{

void validateNumberOfThreads( const int numberOfThreads )
{
    if( numberOfThreads < 1 ) throw std::invalid_argument( "Number of threads must be positive" );
#if !defined( TUDAT_BUILD_WITH_OPENMP ) || !TUDAT_BUILD_WITH_OPENMP
    if( numberOfThreads > 1 ) throw std::runtime_error( "Multiple gravity threads require TUDAT_BUILD_WITH_OPENMP=ON" );
#endif
}

bool useParallelTeam( const int numberOfThreads )
{
#if defined( TUDAT_BUILD_WITH_OPENMP ) && TUDAT_BUILD_WITH_OPENMP
    return numberOfThreads > 1 && !omp_in_parallel( );
#else
    return false;
#endif
}

}  // namespace utilities
}  // namespace tudat
