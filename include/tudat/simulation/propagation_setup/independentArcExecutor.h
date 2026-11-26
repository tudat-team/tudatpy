#ifndef TUDAT_INDEPENDENT_ARC_EXECUTOR_H
#define TUDAT_INDEPENDENT_ARC_EXECUTOR_H

#include <algorithm>
#include <atomic>
#include <exception>
#include <functional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace tudat
{
namespace propagators
{

// Task executor. The caller must supply independent mutable
// environments; this executor cannot make shared Tudat models thread-safe.
inline void executeIsolatedArcs( const std::size_t arcs, const std::size_t workers, const std::function< void( std::size_t ) >& integrate )
{
    if( workers == 0 )
    {
        throw std::invalid_argument( "workers must be positive" );
    }
    std::atomic< std::size_t > next{ 0 };
    std::vector< std::exception_ptr > failures( arcs );
    auto work = [ & ] {
        for( ;; )
        {
            const auto arc = next.fetch_add( 1 );
            if( arc >= arcs )
            {
                return;
            }
            try
            {
                integrate( arc );
            }
            catch( ... )
            {
                failures[ arc ] = std::current_exception( );
            }
        }
    };
    std::vector< std::thread > threads;
    try
    {
        for( std::size_t i = 0; i < std::min( workers, arcs ); ++i )
        {
            threads.emplace_back( work );
        }
    }
    catch( ... )
    {
        for( auto& thread : threads )
        {
            thread.join( );
        }
        throw;
    }
    for( auto& thread : threads )
    {
        thread.join( );
    }
    // Every task completes before failure is reported. Arc index determines
    // which error is reported, independently of thread scheduling.
    for( std::size_t arc = 0; arc < arcs; ++arc )
    {
        if( failures[ arc ] )
        {
            try
            {
                std::rethrow_exception( failures[ arc ] );
            }
            catch( const std::exception& error )
            {
                throw std::runtime_error( "arc " + std::to_string( arc ) + ": " + error.what( ) );
            }
            catch( ... )
            {
                throw std::runtime_error( "arc " + std::to_string( arc ) + ": unknown exception" );
            }
        }
    }
}
}  // namespace propagators
}  // namespace tudat
#endif
