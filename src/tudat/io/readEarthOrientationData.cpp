#include "tudat/io/readEarthOrientationData.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
constexpr double arcsecondToRadian = 3.141592653589793238462643383279502884 / ( 180.0 * 3600.0 );
constexpr double secondsPerDay = 86400.0;

std::string trim( const std::string& text )
{
    const auto first = text.find_first_not_of( " \t\r\n" );
    return first == std::string::npos ? "" : text.substr( first, text.find_last_not_of( " \t\r\n" ) - first + 1 );
}

bool unsupportedCelestialPoleConvention( std::string text )
{
    std::transform( text.begin( ), text.end( ), text.begin( ), []( unsigned char character ) { return std::tolower( character ); } );
    return text.find( "dpsi" ) != std::string::npos || text.find( "deps" ) != std::string::npos;
}

double number( const std::string& text )
{
    std::size_t consumed = 0;
    const double value = std::stod( text, &consumed );
    if( consumed != text.size( ) || !std::isfinite( value ) )
    {
        throw std::runtime_error( "invalid finite number: " + text );
    }
    return value;
}

int integer( const std::string& text )
{
    const double value = number( text );
    if( std::floor( value ) != value || value < 0 || value > 9999 )
    {
        throw std::runtime_error( "invalid calendar integer: " + text );
    }
    return static_cast< int >( value );
}

// Gregorian midnight MJD, independently of any dynamics or SOFA library.
double calendarMjd( int year, const int month, const int day )
{
    const bool leap = year % 4 == 0 && ( year % 100 != 0 || year % 400 == 0 );
    const int monthDays[] = { 31, leap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    if( year < 1600 || month < 1 || month > 12 || day < 1 || day > monthDays[ month - 1 ] )
    {
        throw std::runtime_error( "invalid Gregorian calendar date" );
    }
    const int a = ( 14 - month ) / 12;
    year += 4800 - a;
    const int m = month + 12 * a - 3;
    const int julianDayNumber = day + ( 153 * m + 2 ) / 5 + 365 * year + year / 4 - year / 100 + year / 400 - 32045;
    return julianDayNumber - 2400001.0;
}

void validateDate( const int year, const int month, const int day, const double hour, const double mjd )
{
    if( hour < 0 || hour >= 24 || std::abs( mjd - calendarMjd( year, month, day ) - hour / 24.0 ) > 1.0e-7 )
    {
        throw std::runtime_error( "calendar date/hour and UTC MJD disagree" );
    }
}

void nonnegative( const double value )
{
    if( value < 0 ) throw std::runtime_error( "negative uncertainty" );
}

EarthOrientationRecord c04Record( const std::vector< std::string >& tokens )
{
    const bool modern = tokens.size( ) == 21;
    if( !modern && tokens.size( ) != 16 ) throw std::runtime_error( "C04 requires 16 or 21 columns" );
    std::vector< double > v;
    for( const auto& token : tokens ) v.push_back( number( token ) );
    EarthOrientationRecord record;
    const int offset = modern ? 1 : 0;
    record.utcModifiedJulianDate = v[ 3 + offset ];
    validateDate(
            integer( tokens[ 0 ] ), integer( tokens[ 1 ] ), integer( tokens[ 2 ] ), modern ? v[ 3 ] : 0.0, record.utcModifiedJulianDate );
    record.polarMotion = { v[ 4 + offset ] * arcsecondToRadian, v[ 5 + offset ] * arcsecondToRadian };
    record.ut1MinusUtc = v[ 6 + offset ];
    if( modern )
    {
        record.celestialPoleOffsets = std::array< double, 2 >{ v[ 8 ] * arcsecondToRadian, v[ 9 ] * arcsecondToRadian };
        record.polarMotionRates =
                std::array< double, 2 >{ v[ 10 ] * arcsecondToRadian / secondsPerDay, v[ 11 ] * arcsecondToRadian / secondsPerDay };
        record.lengthOfDay = v[ 12 ];
        record.polarMotionUncertainty = std::array< double, 2 >{ v[ 13 ] * arcsecondToRadian, v[ 14 ] * arcsecondToRadian };
        record.ut1MinusUtcUncertainty = v[ 15 ];
        record.celestialPoleOffsetUncertainty = std::array< double, 2 >{ v[ 16 ] * arcsecondToRadian, v[ 17 ] * arcsecondToRadian };
        record.polarMotionRateUncertainty =
                std::array< double, 2 >{ v[ 18 ] * arcsecondToRadian / secondsPerDay, v[ 19 ] * arcsecondToRadian / secondsPerDay };
        record.lengthOfDayUncertainty = v[ 20 ];
        for( int i = 13; i < 21; ++i ) nonnegative( v[ i ] );
    }
    else
    {
        record.lengthOfDay = v[ 7 ];
        record.celestialPoleOffsets = std::array< double, 2 >{ v[ 8 ] * arcsecondToRadian, v[ 9 ] * arcsecondToRadian };
        record.polarMotionUncertainty = std::array< double, 2 >{ v[ 10 ] * arcsecondToRadian, v[ 11 ] * arcsecondToRadian };
        record.ut1MinusUtcUncertainty = v[ 12 ];
        record.lengthOfDayUncertainty = v[ 13 ];
        record.celestialPoleOffsetUncertainty = std::array< double, 2 >{ v[ 14 ] * arcsecondToRadian, v[ 15 ] * arcsecondToRadian };
        for( int i = 10; i < 16; ++i ) nonnegative( v[ i ] );
    }
    return record;
}

std::string field( const std::string& line, const std::size_t first, const std::size_t count )
{
    return trim( line.substr( first - 1, count ) );  // USNO specification columns are one-based.
}

std::optional< double > optionalNumber( const std::string& line, const int first, const int count, const double scale = 1.0 )
{
    const std::string text = field( line, first, count );
    if( text.empty( ) ) return std::nullopt;
    return number( text ) * scale;
}

std::optional< std::array< double, 2 > > optionalPair( const std::optional< double >& x, const std::optional< double >& y )
{
    if( x.has_value( ) != y.has_value( ) ) throw std::runtime_error( "incomplete optional coordinate pair" );
    return x ? std::optional< std::array< double, 2 > >( std::array< double, 2 >{ *x, *y } ) : std::nullopt;
}

std::string statusFlag( const std::string& line, const int column, const bool required )
{
    const std::string flag = field( line, column, 1 );
    if( flag != "I" && flag != "P" && !( !required && flag.empty( ) ) )
        throw std::runtime_error( "invalid or missing Bulletin A I/P flag" );
    return flag;
}

std::optional< EarthOrientationRecord > finalsRecord( const std::string& line )
{
    // Optional fields must still occupy their fixed columns, even when blank.
    if( line.size( ) < 134 ) throw std::runtime_error( "truncated finals2000A record (requires columns 1-134)" );
    EarthOrientationRecord record;
    record.utcModifiedJulianDate = number( field( line, 8, 8 ) );
    const int year = integer( field( line, 1, 2 ) ) + ( record.utcModifiedJulianDate <= 51543 ? 1900 : 2000 );
    validateDate( year,
                  integer( field( line, 3, 2 ) ),
                  integer( field( line, 5, 2 ) ),
                  24.0 * ( record.utcModifiedJulianDate - std::floor( record.utcModifiedJulianDate ) ),
                  record.utcModifiedJulianDate );
    const auto x = optionalNumber( line, 19, 9, arcsecondToRadian );
    const auto y = optionalNumber( line, 38, 9, arcsecondToRadian );
    const auto ut1 = optionalNumber( line, 59, 10 );
    if( !x && !y && !ut1 && trim( line.substr( 16 ) ).empty( ) ) return std::nullopt;  // date-only future placeholders
    if( !x || !y || !ut1 ) throw std::runtime_error( "missing mandatory Bulletin A polar motion or UT1" );
    record.polarMotion = { *x, *y };
    record.ut1MinusUtc = *ut1;
    record.polarMotionStatus = statusFlag( line, 17, true );
    record.ut1Status = statusFlag( line, 58, true );
    record.celestialPoleOffsets = optionalPair( optionalNumber( line, 98, 9, arcsecondToRadian / 1000.0 ),
                                                optionalNumber( line, 117, 9, arcsecondToRadian / 1000.0 ) );
    record.celestialPoleStatus = statusFlag( line, 96, record.celestialPoleOffsets.has_value( ) );
    if( !record.celestialPoleOffsets && !record.celestialPoleStatus.empty( ) )
        throw std::runtime_error( "CPO flag present without offsets" );
    record.lengthOfDay = optionalNumber( line, 80, 7, 0.001 );
    record.polarMotionUncertainty =
            optionalPair( optionalNumber( line, 28, 9, arcsecondToRadian ), optionalNumber( line, 47, 9, arcsecondToRadian ) );
    record.ut1MinusUtcUncertainty = optionalNumber( line, 69, 10 );
    record.lengthOfDayUncertainty = optionalNumber( line, 87, 7, 0.001 );
    record.celestialPoleOffsetUncertainty = optionalPair( optionalNumber( line, 107, 9, arcsecondToRadian / 1000.0 ),
                                                          optionalNumber( line, 126, 9, arcsecondToRadian / 1000.0 ) );
    if( record.lengthOfDayUncertainty && !record.lengthOfDay ) throw std::runtime_error( "LOD uncertainty without LOD" );
    if( record.celestialPoleOffsetUncertainty && !record.celestialPoleOffsets )
        throw std::runtime_error( "CPO uncertainty without offsets" );
    if( record.polarMotionUncertainty )
        for( auto v : *record.polarMotionUncertainty ) nonnegative( v );
    if( record.ut1MinusUtcUncertainty ) nonnegative( *record.ut1MinusUtcUncertainty );
    if( record.lengthOfDayUncertainty ) nonnegative( *record.lengthOfDayUncertainty );
    if( record.celestialPoleOffsetUncertainty )
        for( auto v : *record.celestialPoleOffsetUncertainty ) nonnegative( v );
    return record;  // Columns 135-185 are Bulletin B and are intentionally not mixed in.
}
}  // namespace

void validateEarthOrientationData( const EarthOrientationData& data )
{
    const std::string source = data.source.empty( ) ? "<data>" : data.source;
    if( data.records.empty( ) ) throw std::runtime_error( source + ": no Earth orientation records" );
    if( data.format != EarthOrientationFormat::c04 && data.format != EarthOrientationFormat::finals2000a )
        throw std::runtime_error( source + ": unsupported Earth orientation format" );
    if( unsupportedCelestialPoleConvention( data.celestialPoleConvention ) )
        throw std::runtime_error( source + ": unsupported dPsi/dEps convention" );
    for( const auto& line : data.header )
        if( unsupportedCelestialPoleConvention( line ) ) throw std::runtime_error( source + ": unsupported dPsi/dEps header" );
    const auto finite = []( const double value ) {
        if( !std::isfinite( value ) ) throw std::runtime_error( "non-finite Earth orientation value" );
    };
    const auto pair = [ &finite ]( const std::optional< std::array< double, 2 > >& values, const bool uncertainty ) {
        if( values )
            for( const double value : *values )
            {
                finite( value );
                if( uncertainty ) nonnegative( value );
            }
    };
    const auto scalar = [ &finite ]( const std::optional< double >& value, const bool uncertainty ) {
        if( value )
        {
            finite( *value );
            if( uncertainty ) nonnegative( *value );
        }
    };
    for( const auto& entry : data.records )
    {
        try
        {
            const auto& record = entry.second;
            finite( entry.first );
            finite( record.utcModifiedJulianDate );
            if( entry.first != record.utcModifiedJulianDate ) throw std::runtime_error( "map key differs from record UTC MJD" );
            finite( record.polarMotion[ 0 ] );
            finite( record.polarMotion[ 1 ] );
            finite( record.ut1MinusUtc );
            scalar( record.lengthOfDay, false );
            pair( record.celestialPoleOffsets, false );
            pair( record.polarMotionRates, false );
            pair( record.polarMotionUncertainty, true );
            scalar( record.ut1MinusUtcUncertainty, true );
            scalar( record.lengthOfDayUncertainty, true );
            pair( record.celestialPoleOffsetUncertainty, true );
            pair( record.polarMotionRateUncertainty, true );
            if( record.lengthOfDayUncertainty && !record.lengthOfDay ) throw std::runtime_error( "LOD uncertainty without LOD" );
            if( record.celestialPoleOffsetUncertainty && !record.celestialPoleOffsets )
                throw std::runtime_error( "CPO uncertainty without offsets" );
            if( record.polarMotionRateUncertainty && !record.polarMotionRates )
                throw std::runtime_error( "rate uncertainty without rates" );
            for( const auto& status : { record.polarMotionStatus, record.ut1Status, record.celestialPoleStatus } )
                if( !status.empty( ) && status != "I" && status != "P" ) throw std::runtime_error( "invalid I/P status" );
            if( !record.celestialPoleStatus.empty( ) && !record.celestialPoleOffsets )
                throw std::runtime_error( "CPO status without offsets" );
        }
        catch( const std::exception& error )
        {
            throw std::runtime_error( source + ": UTC MJD " + std::to_string( entry.first ) + ": " + error.what( ) );
        }
    }
}

EarthOrientationData readEarthOrientationData( std::istream& stream, const EarthOrientationFormat format, const std::string& source )
{
    EarthOrientationData data;
    data.format = format;
    data.source = source;
    data.formatDescription = format == EarthOrientationFormat::finals2000a ? "finals2000A Bulletin A" : "C04";
    data.celestialPoleConvention = format == EarthOrientationFormat::finals2000a ? "dX,dY relative to IAU2000A; Bulletin A" : "dX,dY; C04";
    std::string line;
    std::size_t lineNumber = 0, c04Columns = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        const std::string content = trim( line );
        if( content.empty( ) ) continue;
        try
        {
            if( format == EarthOrientationFormat::c04 )
            {
                if( !std::isdigit( static_cast< unsigned char >( content.front( ) ) ) )
                {
                    std::string lowerContent = content;
                    std::transform( lowerContent.begin( ), lowerContent.end( ), lowerContent.begin( ), []( unsigned char character ) {
                        return std::tolower( character );
                    } );
                    if( lowerContent.find( "dpsi" ) != std::string::npos || lowerContent.find( "deps" ) != std::string::npos )
                        throw std::runtime_error( "unsupported dPsi/dEps celestial pole convention; require dX/dY" );
                    if( !data.records.empty( ) && content.front( ) != '#' )
                        throw std::runtime_error( "unexpected text inside C04 records" );
                    data.header.push_back( line );
                    if( content.find( "ITRF" ) != std::string::npos )
                    {
                        const auto first = content.find( "ITRF" );
                        const auto last = content.find( " -", first );
                        data.terrestrialFrame = trim( content.substr( first, last == std::string::npos ? last : last - first ) );
                    }
                    if( content.find( "YR" ) != std::string::npos && content.find( "HH" ) != std::string::npos &&
                        content.find( "MJD" ) != std::string::npos )
                        c04Columns = 21;
                    continue;
                }
                std::istringstream tokensStream( content );
                std::vector< std::string > tokens;
                for( std::string token; tokensStream >> token; ) tokens.push_back( token );
                if( c04Columns && tokens.size( ) != c04Columns ) throw std::runtime_error( "mixed C04 record formats" );
                auto record = c04Record( tokens );
                c04Columns = tokens.size( );
                data.formatDescription = c04Columns == 21 ? "C04 21-column (hour/rates)" : "C04 16-column";
                if( !data.records.empty( ) && record.utcModifiedJulianDate <= data.records.rbegin( )->first )
                    throw std::runtime_error( "duplicate or unsorted UTC epoch" );
                data.records.emplace( record.utcModifiedJulianDate, std::move( record ) );
            }
            else if( format == EarthOrientationFormat::finals2000a )
            {
                auto record = finalsRecord( line );
                if( !record ) continue;
                if( !data.records.empty( ) && record->utcModifiedJulianDate <= data.records.rbegin( )->first )
                    throw std::runtime_error( "duplicate or unsorted UTC epoch" );
                data.records.emplace( record->utcModifiedJulianDate, std::move( *record ) );
            }
            else
                throw std::runtime_error( "unsupported Earth orientation format" );
        }
        catch( const std::exception& exception )
        {
            throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": " + exception.what( ) );
        }
    }
    if( stream.bad( ) || ( stream.fail( ) && !stream.eof( ) ) ) throw std::runtime_error( source + ": input stream read failure" );
    validateEarthOrientationData( data );
    return data;
}

EarthOrientationData readEarthOrientationData( const std::string& fileName, const EarthOrientationFormat format )
{
    std::ifstream stream( fileName );
    if( !stream ) throw std::runtime_error( "Cannot open Earth orientation file: " + fileName );
    return readEarthOrientationData( stream, format, fileName );
}

}  // namespace input_output
}  // namespace tudat
