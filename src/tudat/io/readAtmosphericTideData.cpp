#include "tudat/io/readAtmosphericTideData.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
void require( bool value, const std::string& message )
{
    if( !value )
    {
        throw std::runtime_error( "Invalid atmospheric tide data: " + message );
    }
}
std::string trim( std::string text )
{
    const auto begin = text.find_first_not_of( " \t\r" );
    if( begin == std::string::npos )
    {
        return "";
    }
    return text.substr( begin, text.find_last_not_of( " \t\r" ) - begin + 1 );
}
double number( std::string text )
{
    std::replace( text.begin( ), text.end( ), 'D', 'E' );
    std::replace( text.begin( ), text.end( ), 'd', 'e' );
    std::size_t end;
    const double value = std::stod( text, &end );
    require( end == text.size( ) && std::isfinite( value ), "malformed or nonfinite number" );
    return value;
}
std::vector< double > numbers( const std::string& line, const std::size_t count )
{
    std::istringstream input( line );
    std::string token;
    std::vector< double > values;
    while( input >> token )
    {
        values.push_back( number( token ) );
    }
    require( values.size( ) == count, "wrong number of coefficient columns" );
    return values;
}
void validateMetadata( const AtmosphericTideMetadata& metadata )
{
    std::string declarations;
    for( const auto& h : metadata.sourceHeaders )
    {
        declarations += h + "\n";
    }
    std::transform(
            declarations.begin( ), declarations.end( ), declarations.begin( ), []( unsigned char c ) { return std::tolower( c ); } );
    declarations.erase( std::remove_if( declarations.begin( ),
                                        declarations.end( ),
                                        []( char c ) { return c == ' ' || c == '\t' || c == ':' || c == '$' || c == '\n' || c == ','; } ),
                        declarations.end( ) );
    std::size_t rowOrder = declarations.find( "roworder" );
    while( rowOrder != std::string::npos )
    {
        require( declarations.compare(
                         rowOrder + 8, std::string( "radialtangentialnstangentialew" ).size( ), "radialtangentialnstangentialew" ) == 0,
                 "incompatible declared component row order" );
        rowOrder = declarations.find( "roworder", rowOrder + 8 );
    }
    require( metadata.originalCoefficientUnits == "mm", "this format requires millimetres" );
    require( metadata.originalCoefficientOrder == "radial,north,east; cosS1,sinS1,cosS2,sinS2", "unsupported coefficient order" );
    require( metadata.referenceOrigin == AtmosphericLoadingReferenceOrigin::unspecified ||
                     metadata.referenceOrigin == AtmosphericLoadingReferenceOrigin::solidEarthCentre ||
                     metadata.referenceOrigin == AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass,
             "unknown origin" );
    require( metadata.displacementBasis == AtmosphericTideDisplacementBasis::unspecified ||
                     metadata.displacementBasis == AtmosphericTideDisplacementBasis::geocentricRadialNorthEast ||
                     metadata.displacementBasis == AtmosphericTideDisplacementBasis::geodeticUpNorthEast,
             "unknown displacement basis" );
    require( metadata.sourceProfile == AtmosphericTideSourceProfile::unspecified ||
                     metadata.sourceProfile == AtmosphericTideSourceProfile::iers2010Rp03NoIbFarrell,
             "unknown source profile" );
    require( metadata.latitudeConvention == AtmosphericTideLatitudeConvention::unspecified ||
                     metadata.latitudeConvention == AtmosphericTideLatitudeConvention::geocentric ||
                     metadata.latitudeConvention == AtmosphericTideLatitudeConvention::geodetic,
             "unknown latitude convention" );
}
void header( AtmosphericTideMetadata& metadata, const std::string& line )
{
    metadata.sourceHeaders.push_back( line );
    std::string lower = line;
    std::transform( lower.begin( ), lower.end( ), lower.begin( ), []( unsigned char c ) { return std::tolower( c ); } );
    if( lower.find( "units" ) != std::string::npos )
    {
        require( lower.find( "millimet" ) != std::string::npos || lower.find( " mm" ) != std::string::npos, "incompatible declared units" );
    }
    if( lower.find( "column order" ) != std::string::npos )
    {
        std::string order = lower.substr( lower.find( "column order" ) + 12 );
        order.erase( std::remove_if( order.begin( ), order.end( ), []( char c ) { return c == ' ' || c == '\t' || c == ':' || c == ','; } ),
                     order.end( ) );
        require( order == "coss1sins1coss2sins2", "incompatible declared harmonic order" );
    }
    if( lower.find( "positive" ) != std::string::npos &&
        ( lower.find( "west" ) != std::string::npos || lower.find( "south" ) != std::string::npos ) )
    {
        require( false, "this format requires positive north and east" );
    }
    if( lower.find( "reference origin" ) != std::string::npos )
    {
        std::string declaration = lower.substr( lower.find( "reference origin" ) + 16 );
        declaration.erase(
                std::remove_if( declaration.begin( ), declaration.end( ), []( char c ) { return c == ' ' || c == '\t' || c == ':'; } ),
                declaration.end( ) );
        require( declaration == "ce" || declaration == "cm", "unknown declared origin" );
        const auto origin = declaration == "ce" ? AtmosphericLoadingReferenceOrigin::solidEarthCentre
                                                : AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass;
        require( metadata.referenceOrigin == AtmosphericLoadingReferenceOrigin::unspecified || metadata.referenceOrigin == origin,
                 "conflicting origin declarations" );
        metadata.referenceOrigin = origin;
    }
    if( lower.find( "farrell" ) != std::string::npos )
    {
        require( metadata.greensFunction.empty( ) || metadata.greensFunction == "Farrell elastic Green functions",
                 "Green-function declaration conflicts with source header" );
        metadata.greensFunction = "Farrell elastic Green functions";
    }
    const std::string joined = [ & ] {
        std::string all;
        for( const auto& h : metadata.sourceHeaders )
        {
            all += h + "\n";
        }
        return all;
    }( );
    if( joined.find( "Ray" ) != std::string::npos && joined.find( "Ponte" ) != std::string::npos &&
        joined.find( "2003" ) != std::string::npos )
    {
        require( metadata.modelName.empty( ) || metadata.modelName == "Ray-Ponte 2003", "model declaration conflicts with source header" );
        metadata.modelName = "Ray-Ponte 2003";
    }
}
void axis( const std::vector< double >& values, const bool latitude )
{
    require( values.size( ) >= 2, "each grid axis needs at least two nodes" );
    const double step = values[ 1 ] - values[ 0 ];
    require( std::isfinite( step ) && step != 0, "zero or nonfinite grid spacing" );
    for( std::size_t i = 0; i < values.size( ); ++i )
    {
        require( std::isfinite( values[ i ] ), "nonfinite axis node" );
        require( std::abs( values[ i ] - values[ 0 ] - i * step ) < 1e-9, "irregular or reordered axis" );
        require( std::abs( values[ i ] ) <= ( latitude ? 90.0 : 360.0 ), "axis coordinate outside degree bounds" );
    }
    if( !latitude )
    {
        require( std::abs( values.back( ) - values.front( ) ) <= 360.0, "longitude span exceeds one period" );
    }
}
Eigen::Matrix< double, 3, 4 > matrix( const std::vector< double >& values, const std::size_t offset )
{
    Eigen::Matrix< double, 3, 4 > result;
    for( int row = 0; row < 3; ++row )
    {
        for( int column = 0; column < 4; ++column )
        {
            result( row, column ) = values[ offset + 4 * row + column ] * 1e-3;
        }
    }
    return result;
}
}  // namespace

void validateAtmosphericTideGrid( const AtmosphericTideGrid& grid )
{
    validateMetadata( grid.metadata );
    axis( grid.longitudesDegrees, false );
    axis( grid.latitudesDegrees, true );
    require( grid.metadata.expectedLongitudeCount == 0 || grid.metadata.expectedLongitudeCount == grid.longitudesDegrees.size( ),
             "unexpected longitude count or terminal truncation" );
    require( grid.metadata.expectedLatitudeCount == 0 || grid.metadata.expectedLatitudeCount == grid.latitudesDegrees.size( ),
             "unexpected latitude count or terminal truncation" );
    require( grid.longitudesDegrees.size( ) <= std::numeric_limits< std::size_t >::max( ) / grid.latitudesDegrees.size( ),
             "node count overflow" );
    require( grid.coefficientsMetres.size( ) == grid.longitudesDegrees.size( ) * grid.latitudesDegrees.size( ), "missing or excess nodes" );
    for( const auto& node : grid.coefficientsMetres )
    {
        require( node.allFinite( ), "nonfinite coefficients" );
    }
}

//! Select a provider station, reporting explicit absence separately from malformed data.
const AtmosphericTideStation& AtmosphericTideStationData::getStation( const std::string& identifier ) const
{
    const auto station = stations.find( identifier );

    if( station == stations.end( ) )
    {
        throw std::out_of_range( "Atmospheric tide station not found: " + identifier );
    }
    return station->second;
}

void validateAtmosphericTideStationData( const AtmosphericTideStationData& data )
{
    validateMetadata( data.metadata );
    require( !data.stations.empty( ), "no stations" );
    for( const auto& entry : data.stations )
    {
        const auto& station = entry.second;
        require( !entry.first.empty( ) && entry.first == station.stationIdentifier, "inconsistent station identity" );
        validateMetadata( station.metadata );
        require( station.coefficientsMetres.allFinite( ), "nonfinite station coefficients" );
        if( station.hasProviderCoordinates )
        {
            require( station.providerLongitudeLatitudeDegrees.allFinite( ) &&
                             std::abs( station.providerLongitudeLatitudeDegrees.y( ) ) <= 90,
                     "invalid provider station coordinates" );
        }
    }
}

AtmosphericTideGrid readAtmosphericTideGrid( std::istream& stream, const AtmosphericTideMetadata& metadata )
{
    std::size_t lineNumber = 0;
    try
    {
        AtmosphericTideGrid grid;
        grid.metadata = metadata;
        validateMetadata( metadata );
        std::vector< std::vector< double > > rows;
        std::string line;
        while( std::getline( stream, line ) )
        {
            ++lineNumber;
            const std::string providerLine = !line.empty( ) && line.back( ) == '\r' ? line.substr( 0, line.size( ) - 1 ) : line;
            line = trim( line );
            if( line.empty( ) )
            {
                continue;
            }
            if( line.front( ) == '$' )
            {
                header( grid.metadata, providerLine );
                continue;
            }
            rows.push_back( numbers( line, 14 ) );
        }
        require( !stream.bad( ) && !( stream.fail( ) && !stream.eof( ) ) && !rows.empty( ), "empty grid or stream failure" );
        const double firstLongitude = rows.front( )[ 0 ];
        std::size_t latitudeCount = 0;
        while( latitudeCount < rows.size( ) && rows[ latitudeCount ][ 0 ] == firstLongitude )
        {
            grid.latitudesDegrees.push_back( rows[ latitudeCount++ ][ 1 ] );
        }
        require( latitudeCount >= 2 && rows.size( ) % latitudeCount == 0, "truncated grid longitude block" );
        for( std::size_t i = 0; i < rows.size( ); ++i )
        {
            if( i % latitudeCount == 0 )
            {
                grid.longitudesDegrees.push_back( rows[ i ][ 0 ] );
            }
            require( rows[ i ][ 0 ] == grid.longitudesDegrees.back( ) && rows[ i ][ 1 ] == grid.latitudesDegrees[ i % latitudeCount ],
                     "duplicate, missing, reordered or off-grid node" );
            grid.coefficientsMetres.push_back( matrix( rows[ i ], 2 ) );
        }
        validateAtmosphericTideGrid( grid );
        return grid;
    }
    catch( const std::exception& error )
    {
        throw std::runtime_error( ( metadata.source.empty( ) ? "<stream>" : metadata.source ) + ":" + std::to_string( lineNumber ) + ": " +
                                  error.what( ) );
    }
}

AtmosphericTideStationData readAtmosphericTideStationData( std::istream& stream, const AtmosphericTideMetadata& metadata )
{
    std::size_t lineNumber = 0;
    try
    {
        AtmosphericTideStationData data;
        data.metadata = metadata;
        validateMetadata( metadata );
        AtmosphericTideStation current;
        int row = -1;
        bool hasPendingCoordinates = false;
        Eigen::Vector2d pendingCoordinates = Eigen::Vector2d::Zero( );
        std::string line;
        while( std::getline( stream, line ) )
        {
            ++lineNumber;
            const std::string providerLine = !line.empty( ) && line.back( ) == '\r' ? line.substr( 0, line.size( ) - 1 ) : line;
            line = trim( line );
            if( line.empty( ) )
            {
                continue;
            }
            if( line.front( ) == '$' )
            {
                header( data.metadata, providerLine );
                const auto coordinates = line.find( "coord." );
                if( coordinates != std::string::npos )
                {
                    const auto end = line.find( ')', coordinates );
                    require( end != std::string::npos, "malformed station coordinates header" );
                    std::istringstream values( line.substr( end + 1 ) );
                    std::string lon, lat, extra;
                    require( static_cast< bool >( values >> lon >> lat ) && !( values >> extra ), "malformed station coordinates" );
                    const Eigen::Vector2d declaredCoordinates( number( lon ), number( lat ) );
                    require( std::abs( declaredCoordinates.y( ) ) <= 90, "provider latitude outside [-90,90]" );
                    if( row >= 0 )
                    {
                        require( !current.hasProviderCoordinates ||
                                         ( current.providerLongitudeLatitudeDegrees - declaredCoordinates ).norm( ) < 1e-9,
                                 "conflicting current station coordinates" );
                        current.providerLongitudeLatitudeDegrees = declaredCoordinates;
                        current.hasProviderCoordinates = true;
                    }
                    else
                    {
                        require( !hasPendingCoordinates || ( pendingCoordinates - declaredCoordinates ).norm( ) < 1e-9,
                                 "conflicting pending station coordinates" );
                        pendingCoordinates = declaredCoordinates;
                        hasPendingCoordinates = true;
                    }
                }
                if( row >= 0 )
                {
                    header( current.metadata, providerLine );
                }
                continue;
            }
            if( row < 0 )
            {
                std::istringstream identity( line );
                std::string identifier;
                identity >> identifier;
                // Numeric identifiers are legitimate (e.g. SLR7839). Distinguish an orphan FOUR-number coefficient row structurally.
                std::istringstream fields( line );
                std::vector< std::string > tokens;
                std::string token;
                while( fields >> token )
                {
                    tokens.push_back( token );
                }
                bool orphanRow = tokens.size( ) == 4;
                for( auto value : tokens )
                {
                    std::replace( value.begin( ), value.end( ), 'D', 'E' );
                    std::replace( value.begin( ), value.end( ), 'd', 'e' );
                    char* end = nullptr;
                    std::strtod( value.c_str( ), &end );
                    orphanRow = orphanRow && end != value.c_str( ) && *end == '\0';
                }
                require( !identifier.empty( ) && !orphanRow, "missing station identifier before coefficient row" );
                require( data.stations.count( identifier ) == 0, "duplicate station identifier" );
                current = AtmosphericTideStation( );
                current.stationIdentifier = identifier;
                current.providerStationIdentifierLine = providerLine;
                current.metadata = data.metadata;
                current.hasProviderCoordinates = hasPendingCoordinates;
                current.providerLongitudeLatitudeDegrees = pendingCoordinates;
                hasPendingCoordinates = false;
                row = 0;
            }
            else
            {
                const auto values = numbers( line, 4 );
                for( int column = 0; column < 4; ++column )
                {
                    current.coefficientsMetres( row, column ) = values[ column ] * 1e-3;
                }
                if( ++row == 3 )
                {
                    data.stations.emplace( current.stationIdentifier, current );
                    row = -1;
                }
            }
        }
        require( !stream.bad( ) && !( stream.fail( ) && !stream.eof( ) ) && row < 0, "truncated station block or stream failure" );
        validateAtmosphericTideStationData( data );
        return data;
    }
    catch( const std::exception& error )
    {
        throw std::runtime_error( ( metadata.source.empty( ) ? "<stream>" : metadata.source ) + ":" + std::to_string( lineNumber ) + ": " +
                                  error.what( ) );
    }
}
AtmosphericTideGrid readAtmosphericTideGrid( const std::string& path, const AtmosphericTideMetadata& metadata )
{
    std::ifstream stream( path );
    require( stream.good( ), "cannot open " + path );
    auto supplied = metadata;
    supplied.source = path;
    return readAtmosphericTideGrid( stream, supplied );
}
AtmosphericTideStationData readAtmosphericTideStationData( const std::string& path, const AtmosphericTideMetadata& metadata )
{
    std::ifstream stream( path );
    require( stream.good( ), "cannot open " + path );
    auto supplied = metadata;
    supplied.source = path;
    return readAtmosphericTideStationData( stream, supplied );
}
}  // namespace input_output
}  // namespace tudat
