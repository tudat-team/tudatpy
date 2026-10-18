#include "tudat/io/readOceanPoleTideData.h"
#include "tudat/io/textParsing.h"

#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>

namespace tudat
{
namespace input_output
{

OceanPoleTideGrid readOceanPoleTideGrid( std::istream& stream, const std::string& source )
{
    OceanPoleTideGrid grid;
    std::map< std::string, double > geometry;
    std::string line;
    std::size_t lineNumber = 0;
    const auto fail = [ & ]( const std::string& message ) {
        throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": invalid ocean pole-tide grid: " + message );
    };
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        const auto tokens = text_parsing::splitIntoTokens( line );
        if( tokens.empty( ) )
        {
            continue;
        }

        // Geometry header lines have the form "name = value"; data lines contain eight numbers.
        if( tokens.size( ) == 3 && tokens[ 1 ] == "=" )
        {
            geometry[ tokens[ 0 ] ] = text_parsing::parseFiniteNumber( tokens[ 2 ] ).value_or( NAN );
            continue;
        }
        std::vector< double > row;
        for( const auto& token : tokens )
        {
            const auto value = text_parsing::parseFiniteNumber( token );
            if( !value )
            {
                break;
            }
            row.push_back( *value );
        }
        if( row.size( ) != tokens.size( ) )
        {
            continue;
        }
        if( row.size( ) != 8 )
        {
            fail( "expected eight columns" );
        }
        if( grid.coefficients.empty( ) )
        {
            const auto get = [ & ]( const std::string& key ) {
                if( !geometry.count( key ) || !std::isfinite( geometry.at( key ) ) )
                {
                    fail( "missing geometry field " + key );
                }
                return geometry.at( key );
            };
            grid.longitudeCount = static_cast< std::size_t >( get( "Number_longitude_Grid_Points" ) );
            grid.latitudeCount = static_cast< std::size_t >( get( "Number_latitude_grid_points" ) );
            grid.firstLongitudeDegrees = get( "First_longitude_degrees" );
            grid.longitudeStepDegrees = get( "Longitude_step_degrees" );
            grid.firstLatitudeDegrees = get( "First_latitude_degrees" );
            grid.latitudeStepDegrees = get( "Latitude_step_degrees" );
        }
        const std::size_t node = grid.coefficients.size( );
        if( std::abs( row[ 0 ] - grid.firstLongitudeDegrees - ( node % grid.longitudeCount ) * grid.longitudeStepDegrees ) > 1e-9 ||
            std::abs( row[ 1 ] - grid.firstLatitudeDegrees - ( node / grid.longitudeCount ) * grid.latitudeStepDegrees ) > 1e-9 )
        {
            fail( "missing, reordered or off-grid node" );
        }
        grid.coefficients.push_back( Eigen::Map< const Eigen::Matrix< double, 2, 3 > >( row.data( ) + 2 ).transpose( ) );
    }
    if( grid.longitudeCount < 2 || grid.latitudeCount < 2 || grid.coefficients.size( ) != grid.longitudeCount * grid.latitudeCount )
    {
        fail( "incomplete grid" );
    }
    return grid;
}

OceanPoleTideGrid readOceanPoleTideGrid( const std::string& filePath )
{
    std::ifstream stream( filePath );
    if( !stream )
    {
        throw std::runtime_error( "Cannot open ocean pole-tide grid: " + filePath );
    }
    return readOceanPoleTideGrid( stream, filePath );
}
}  // namespace input_output
}  // namespace tudat
