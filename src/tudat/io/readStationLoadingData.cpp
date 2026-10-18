#include "tudat/io/readStationLoadingData.h"
#include "tudat/io/textParsing.h"

#include <fstream>
#include <optional>
#include <stdexcept>

namespace tudat
{
namespace input_output
{

StationLoadingData readEphedispStationLoadingData( std::istream& stream, const std::string& source )
{
    const std::string signature = "EPHEDISP  Format version of 2005.06.30";
    StationLoadingData data;
    data.source = source;
    std::size_t lineNumber = 0, epochCount = 0;
    std::optional< Time > begin;
    std::optional< long double > interval;
    double radius = 0.0;
    bool started = false;
    std::string text;
    const auto fail = [ & ]( const std::string& message ) {
        throw std::runtime_error( "Station loading " + source + " line " + std::to_string( lineNumber ) + ": " + message );
    };
    const auto field = [ & ]( std::size_t start, std::size_t width ) {
        const auto value = text.size( ) < start + width
                ? std::nullopt
                : text_parsing::parseFiniteNumber( text_parsing::trim( text.substr( start, width ) ) );
        if( !value )
        {
            fail( "invalid or truncated numerical field" );
        }
        return *value;
    };
    const auto label = [ & ]( std::size_t start ) {
        const auto value = text.size( ) < start + 8 ? "" : text_parsing::trim( text.substr( start, 8 ) );
        if( value.empty( ) )
        {
            fail( "missing site identifier" );
        }
        return value;
    };
    while( std::getline( stream, text ) )
    {
        ++lineNumber;
        if( !text.empty( ) && text.back( ) == '\r' )
        {
            text.pop_back( );
        }
        if( text == signature )
        {
            if( started )
            {
                break;
            }
            started = true;
            continue;
        }
        if( !started )
        {
            fail( "missing EPHEDISP signature" );
        }
        if( text.empty( ) || text.front( ) == '#' )
        {
            continue;
        }
        if( text.rfind( "P T 3 S ", 0 ) == 0 )
        {
            epochCount = static_cast< std::size_t >( field( 21, 6 ) );
        }
        else if( text.front( ) == 'A' )
        {
            radius = field( 2, 14 );
        }
        else if( text.rfind( "T sample", 0 ) == 0 )
        {
            interval = static_cast< long double >( field( 10, 16 ) ) * 86400.0L;
        }
        else if( text.rfind( "T begin ", 0 ) == 0 )
        {
            begin = Time( ( static_cast< int >( field( 10, 5 ) ) - 51544 ) * 24 - 12, static_cast< long double >( field( 16, 7 ) ) );
        }
        else if( text.front( ) == 'S' )
        {
            StationLoadingSite site;
            site.cartesianPosition << field( 13, 13 ), field( 27, 13 ), field( 41, 13 );
            site.validityRadiusMetres = radius;
            if( !data.sites.emplace( label( 3 ), site ).second )
            {
                fail( "duplicate site" );
            }
        }
        else if( text.front( ) == 'D' )
        {
            const int index = static_cast< int >( field( 2, 5 ) );
            const auto site = data.sites.find( label( 45 ) );
            if( index < 1 || static_cast< std::size_t >( index ) > epochCount || site == data.sites.end( ) )
            {
                fail( "displacement outside the declared epochs or sites" );
            }
            auto& displacements = site->second.displacements;
            if( !displacements.empty( ) && static_cast< std::size_t >( index - 1 ) != displacements.rbegin( )->first + 1 )
            {
                fail( "displacements of a site must be consecutive" );
            }
            displacements[ index - 1 ] << field( 54, 8 ), field( 63, 8 ), field( 72, 8 );
        }
        else if( text.front( ) != 'T' )
        {
            fail( "unknown record" );
        }
    }
    if( !begin || !interval || epochCount == 0 || data.sites.empty( ) )
    {
        throw std::runtime_error( "Station loading " + source + ": incomplete EPHEDISP header or no sites" );
    }
    for( std::size_t i = 0; i < epochCount; ++i )
    {
        data.epochs.push_back( *begin + static_cast< long double >( i ) * *interval );
    }
    return data;
}

StationLoadingData readEphedispStationLoadingData( const std::string& path )
{
    std::ifstream stream( path );
    if( !stream )
    {
        throw std::runtime_error( "Station loading " + path + ": cannot open file" );
    }
    return readEphedispStationLoadingData( stream, path );
}

}  // namespace input_output
}  // namespace tudat
