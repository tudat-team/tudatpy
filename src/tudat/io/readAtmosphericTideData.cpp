#include "tudat/io/readAtmosphericTideData.h"
#include "tudat/io/textParsing.h"
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
//! Parse the numbers of a line, which must contain the given number of values.
std::vector< double > parseNumbers( const std::string& line,
                                    const std::size_t count,
                                    const std::string& source,
                                    const std::size_t lineNumber )
{
    std::vector< double > values;
    for( const auto& token : text_parsing::splitIntoTokens( line ) )
    {
        const auto value = text_parsing::parseFiniteNumber( token );
        values.push_back( value ? *value : std::nan( "" ) );
    }
    if( values.size( ) != count || !Eigen::Map< Eigen::VectorXd >( values.data( ), values.size( ) ).allFinite( ) )
    {
        throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": expected " + std::to_string( count ) + " numbers" );
    }
    return values;
}
}  // namespace

AtmosphericTideGrid readAtmosphericTideGrid( std::istream& stream, const std::string& source )
{
    AtmosphericTideGrid grid;
    std::string line;
    std::size_t lineNumber = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        line = text_parsing::trim( line );
        if( line.empty( ) || line.front( ) == '$' )
        {
            continue;
        }
        const auto values = parseNumbers( line, 14, source, lineNumber );
        if( grid.longitudesDegrees.empty( ) || values[ 0 ] != grid.longitudesDegrees.back( ) )
        {
            grid.longitudesDegrees.push_back( values[ 0 ] );
        }
        if( grid.longitudesDegrees.size( ) == 1 )
        {
            grid.latitudesDegrees.push_back( values[ 1 ] );
        }
        const std::size_t index = grid.coefficientsMetres.size( ) % grid.latitudesDegrees.size( );
        if( grid.longitudesDegrees.size( ) > 1 && values[ 1 ] != grid.latitudesDegrees.at( index ) )
        {
            throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": missing or reordered grid point" );
        }
        grid.coefficientsMetres.push_back( Eigen::Map< const Eigen::Matrix< double, 4, 3 > >( values.data( ) + 2 ).transpose( ) * 1.0E-3 );
    }
    if( grid.latitudesDegrees.size( ) < 2 || grid.longitudesDegrees.size( ) < 2 ||
        grid.coefficientsMetres.size( ) != grid.longitudesDegrees.size( ) * grid.latitudesDegrees.size( ) )
    {
        throw std::runtime_error( source + ": incomplete atmospheric tide grid" );
    }
    return grid;
}

AtmosphericTideGrid readAtmosphericTideGrid( const std::string& path )
{
    std::ifstream stream( path );
    if( !stream )
    {
        throw std::runtime_error( "Cannot open atmospheric tide grid " + path );
    }
    return readAtmosphericTideGrid( stream, path );
}

std::map< std::string, AtmosphericTideCoefficients > readAtmosphericTideStationData( std::istream& stream, const std::string& source )
{
    std::map< std::string, AtmosphericTideCoefficients > stations;
    std::string line, station;
    int row = -1;
    std::size_t lineNumber = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        line = text_parsing::trim( line );
        if( line.empty( ) || line.front( ) == '$' )
        {
            continue;
        }
        if( row < 0 )
        {
            station = text_parsing::splitIntoTokens( line ).front( );
            if( !stations.emplace( station, AtmosphericTideCoefficients::Zero( ) ).second )
            {
                throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": duplicate station " + station );
            }
            row = 0;
            continue;
        }
        const auto values = parseNumbers( line, 4, source, lineNumber );
        for( int column = 0; column < 4; ++column )
        {
            stations.at( station )( row, column ) = values[ column ] * 1.0E-3;
        }
        row = row == 2 ? -1 : row + 1;
    }
    if( row >= 0 || stations.empty( ) )
    {
        throw std::runtime_error( source + ": truncated or empty atmospheric tide station file" );
    }
    return stations;
}

std::map< std::string, AtmosphericTideCoefficients > readAtmosphericTideStationData( const std::string& path )
{
    std::ifstream stream( path );
    if( !stream )
    {
        throw std::runtime_error( "Cannot open atmospheric tide station file " + path );
    }
    return readAtmosphericTideStationData( stream, path );
}

}  // namespace input_output
}  // namespace tudat
