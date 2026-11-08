#include "isolatedArcExecutor.h"
#include <iostream>

int main( )
{
    for( const std::size_t workers : { 1, 2, 4, 16 } )
    {
        std::vector< int > counts( 1000, 0 );
        executeIsolatedArcs( counts.size( ), workers, [ & ]( std::size_t arc ) { ++counts.at( arc ); } );
        for( int value : counts )
        {
            if( value != 1 )
            {
                throw std::runtime_error( "task skipped or repeated" );
            }
        }
    }
    executeIsolatedArcs( 0, 4, []( std::size_t ) { throw std::runtime_error( "unexpected task" ); } );
    int one = 0;
    executeIsolatedArcs( 1, 16, [ & ]( std::size_t ) { ++one; } );
    if( one != 1 )
    {
        throw std::runtime_error( "single arc failed" );
    }
    bool invalid = false;
    try
    {
        executeIsolatedArcs( 1, 0, []( std::size_t ) {} );
    }
    catch( const std::invalid_argument& )
    {
        invalid = true;
    }
    if( !invalid )
    {
        throw std::runtime_error( "accepted zero workers" );
    }
    std::vector< int > completed( 9, 0 );
    bool reported = false;
    try
    {
        executeIsolatedArcs( completed.size( ), 4, [ & ]( std::size_t arc ) {
            completed.at( arc ) = 1;
            if( arc == 3 || arc == 7 )
            {
                throw std::runtime_error( "injected" );
            }
        } );
    }
    catch( const std::runtime_error& error )
    {
        reported = std::string( error.what( ) ) == "arc 3: injected";
    }
    if( !reported )
    {
        throw std::runtime_error( "wrong error attribution" );
    }
    for( int value : completed )
    {
        if( value != 1 )
        {
            throw std::runtime_error( "returned before all tasks finished" );
        }
    }
    std::cout << "executor tests passed\n";
}
