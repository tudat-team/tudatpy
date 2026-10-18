#include "tudat/io/readEarthOrientationData.h"
#include "tudat/io/textParsing.h"

#include <cctype>
#include <fstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
constexpr double arcsecondToRadian = 3.141592653589793238462643383279502884 / ( 180.0 * 3600.0 );

double parseEarthOrientationNumber( const std::string& text )
{
    const auto value = text_parsing::parseFiniteNumber( text );
    if( !value )
    {
        throw std::runtime_error( "invalid finite number: " + text );
    }
    return *value;
}

EarthOrientationRecord c04Record( const std::vector< std::string >& tokens )
{
    const bool modern = tokens.size( ) == 21;
    if( !modern && tokens.size( ) != 16 )
    {
        throw std::runtime_error( "C04 requires 16 or 21 columns" );
    }
    std::vector< double > v;
    for( const auto& token : tokens )
    {
        v.push_back( parseEarthOrientationNumber( token ) );
    }

    // The 21-column format has an additional hour column, and the excess length of day after the pole offset rates.
    const int offset = modern ? 1 : 0;
    EarthOrientationRecord record;
    record.utcModifiedJulianDate = v[ 3 + offset ];
    record.polarMotion = { v[ 4 + offset ] * arcsecondToRadian, v[ 5 + offset ] * arcsecondToRadian };
    record.ut1MinusUtc = v[ 6 + offset ];
    record.celestialPoleOffsets = std::array< double, 2 >{ v[ 8 ] * arcsecondToRadian, v[ 9 ] * arcsecondToRadian };
    record.lengthOfDay = modern ? v[ 12 ] : v[ 7 ];
    return record;
}

std::optional< double > optionalNumber( const std::string& line, const int first, const int count, const double scale = 1.0 )
{
    // USNO specification columns are one-based.
    const std::string text = text_parsing::trim( line.substr( first - 1, count ) );
    return text.empty( ) ? std::nullopt : std::optional< double >( parseEarthOrientationNumber( text ) * scale );
}

std::optional< EarthOrientationRecord > finalsRecord( const std::string& line )
{
    if( line.size( ) < 134 )
    {
        throw std::runtime_error( "truncated finals2000A record (requires columns 1-134)" );
    }
    EarthOrientationRecord record;
    record.utcModifiedJulianDate = *optionalNumber( line, 8, 8 );
    const auto x = optionalNumber( line, 19, 9, arcsecondToRadian );
    const auto y = optionalNumber( line, 38, 9, arcsecondToRadian );
    const auto ut1 = optionalNumber( line, 59, 10 );
    if( !x && !y && !ut1 )
    {
        return std::nullopt;  // date-only future placeholders
    }
    if( !x || !y || !ut1 )
    {
        throw std::runtime_error( "missing mandatory Bulletin A polar motion or UT1" );
    }
    record.polarMotion = { *x, *y };
    record.ut1MinusUtc = *ut1;
    const auto dX = optionalNumber( line, 98, 9, arcsecondToRadian / 1000.0 );
    const auto dY = optionalNumber( line, 117, 9, arcsecondToRadian / 1000.0 );
    if( dX && dY )
    {
        record.celestialPoleOffsets = std::array< double, 2 >{ *dX, *dY };
    }
    record.lengthOfDay = optionalNumber( line, 80, 7, 0.001 );
    return record;  // Columns 135-185 are Bulletin B and are intentionally not mixed in.
}
}  // namespace

EarthOrientationData readEarthOrientationData( std::istream& stream, const EarthOrientationFormat format, const std::string& source )
{
    EarthOrientationData data;
    data.format = format;
    data.source = source;
    std::string line;
    std::size_t lineNumber = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        const std::string content = text_parsing::trim( line );

        // C04 header lines do not start with a digit.
        if( content.empty( ) ||
            ( format == EarthOrientationFormat::c04 && !std::isdigit( static_cast< unsigned char >( content.front( ) ) ) ) )
        {
            continue;
        }
        try
        {
            const auto record = format == EarthOrientationFormat::c04
                    ? std::optional< EarthOrientationRecord >( c04Record( text_parsing::splitIntoTokens( content ) ) )
                    : finalsRecord( line );
            if( !record )
            {
                continue;
            }
            if( !data.records.empty( ) && record->utcModifiedJulianDate <= data.records.rbegin( )->first )
            {
                throw std::runtime_error( "duplicate or unsorted UTC epoch" );
            }
            data.records.emplace( record->utcModifiedJulianDate, *record );
        }
        catch( const std::exception& exception )
        {
            throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": " + exception.what( ) );
        }
    }
    if( data.records.empty( ) )
    {
        throw std::runtime_error( source + ": no Earth orientation records" );
    }
    return data;
}

EarthOrientationData readEarthOrientationData( const std::string& fileName, const EarthOrientationFormat format )
{
    std::ifstream stream( fileName );
    if( !stream )
    {
        throw std::runtime_error( "Cannot open Earth orientation file: " + fileName );
    }
    return readEarthOrientationData( stream, format, fileName );
}

}  // namespace input_output
}  // namespace tudat
