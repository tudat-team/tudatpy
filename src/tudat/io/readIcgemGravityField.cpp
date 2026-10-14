/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/io/readIcgemGravityField.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
[[noreturn]] void fail( const std::string& source, std::size_t line, const std::string& message )
{
    throw std::runtime_error( source + ":" + std::to_string( line ) + ": " + message );
}
std::vector< std::string > tokens( const std::string& line )
{
    std::istringstream stream( line );
    std::vector< std::string > result;
    std::string token;
    while( stream >> token ) result.push_back( token );
    return result;
}
double number( std::string token, const std::string& source, std::size_t line )
{
    std::replace( token.begin( ), token.end( ), 'D', 'E' );
    std::replace( token.begin( ), token.end( ), 'd', 'e' );
    try
    {
        std::size_t used;
        const double value = std::stod( token, &used );
        if( used == token.size( ) && std::isfinite( value ) ) return value;
    }
    catch( const std::exception& )
    {}
    fail( source, line, "Invalid finite number: " + token );
}
int integer( const std::string& token, const std::string& source, std::size_t line )
{
    if( !token.empty( ) && std::all_of( token.begin( ), token.end( ), []( char c ) { return c >= '0' && c <= '9'; } ) ) try
        {
            return std::stoi( token );
        }
        catch( const std::exception& )
        {}
    fail( source, line, "Invalid nonnegative integer: " + token );
}
std::size_t errorColumns( IcgemCoefficientErrors errors )
{
    switch( errors )
    {
        case IcgemCoefficientErrors::none:
            return 0;
        case IcgemCoefficientErrors::calibrated:
        case IcgemCoefficientErrors::formal:
            return 2;
        case IcgemCoefficientErrors::calibratedAndFormal:
            return 4;
        default:
            throw std::runtime_error( "Invalid ICGEM error-column convention" );
    }
}
void readHeader( IcgemGravityFieldData& data )
{
    const std::set< std::string > mandatory{ "product_type", "modelname", "earth_gravity_constant", "radius", "max_degree", "errors" };
    auto supported = mandatory;
    supported.insert( { "format", "norm", "tide_system" } );
    std::set< std::string > seen;
    std::size_t start = 0;
    bool beginSeen = false;
    for( std::size_t i = 0; i < data.headerLines.size( ); ++i )
    {
        const auto row = tokens( data.headerLines[ i ] );
        if( !row.empty( ) && row.front( ) == "begin_of_head" )
        {
            if( beginSeen ) fail( data.source, i + 1, "Duplicate begin_of_head" );
            beginSeen = true;
            start = i + 1;
        }
    }
    for( std::size_t i = start; i < data.headerLines.size( ); ++i )
    {
        const auto row = tokens( data.headerLines[ i ] );
        if( row.empty( ) || !supported.count( row.front( ) ) ) continue;
        if( row.size( ) < 2 || !seen.insert( row.front( ) ).second )
            fail( data.source, i + 1, "Missing value or duplicate header keyword: " + row.front( ) );
        const auto& key = row.front( );
        const auto& value = row[ 1 ];
        if( key == "product_type" )
            data.productType = value;
        else if( key == "modelname" )
            data.modelName = value;
        else if( key == "earth_gravity_constant" )
            data.gravitationalParameter = number( value, data.source, i + 1 );
        else if( key == "radius" )
            data.referenceRadius = number( value, data.source, i + 1 );
        else if( key == "max_degree" )
            data.maximumDegree = integer( value, data.source, i + 1 );
        else if( key == "format" )
            data.format = value;
        else if( key == "errors" )
        {
            if( value == "no" )
                data.errors = IcgemCoefficientErrors::none;
            else if( value == "calibrated" )
                data.errors = IcgemCoefficientErrors::calibrated;
            else if( value == "formal" )
                data.errors = IcgemCoefficientErrors::formal;
            else if( value == "calibrated_and_formal" )
                data.errors = IcgemCoefficientErrors::calibratedAndFormal;
            else
                fail( data.source, i + 1, "Unsupported errors convention: " + value );
        }
        else if( key == "norm" )
        {
            if( value == "fully_normalized" )
                data.normalization = IcgemNormalization::fullyNormalized;
            else if( value == "unnormalized" )
                data.normalization = IcgemNormalization::unnormalized;
            else
                fail( data.source, i + 1, "Unsupported normalization: " + value );
        }
        else if( key == "tide_system" )
        {
            using basic_astrodynamics::TideSystem;
            if( value == "unknown" )
                data.tideSystem = TideSystem::unknown;
            else if( value == "tide_free" )
                data.tideSystem = TideSystem::tideFree;
            else if( value == "zero_tide" )
                data.tideSystem = TideSystem::zeroTide;
            else if( value == "mean_tide" )
                data.tideSystem = TideSystem::meanTide;
            else
                fail( data.source, i + 1, "Unsupported tide system: " + value );
        }
    }
    for( const auto& key : mandatory )
        if( !seen.count( key ) ) fail( data.source, data.headerLines.size( ), "Missing mandatory header keyword: " + key );
}
}  // namespace

void validateIcgemGravityFieldData( const IcgemGravityFieldData& data )
{
    using basic_astrodynamics::TideSystem;
    if( data.productType != "gravity_field" || data.modelName.empty( ) || data.maximumDegree < 0 ||
        !std::isfinite( data.gravitationalParameter ) || data.gravitationalParameter <= 0 || !std::isfinite( data.referenceRadius ) ||
        data.referenceRadius <= 0 || ( data.format != "icgem1.0" && data.format != "icgem2.0" ) ||
        ( data.normalization != IcgemNormalization::fullyNormalized && data.normalization != IcgemNormalization::unnormalized ) ||
        ( data.tideSystem != TideSystem::unknown && data.tideSystem != TideSystem::tideFree && data.tideSystem != TideSystem::zeroTide &&
          data.tideSystem != TideSystem::meanTide ) )
        fail( data.source, 0, "Invalid static ICGEM metadata" );
    errorColumns( data.errors );
    for( const auto& entry : data.coefficients )
    {
        const auto& key = entry.first;
        const auto& coefficient = entry.second;
        if( key.first < 0 || key.second < 0 || key.second > key.first || key.first > data.maximumDegree ||
            !std::isfinite( coefficient.cosine ) || !std::isfinite( coefficient.sine ) ||
            coefficient.calibratedErrors.has_value( ) !=
                    ( data.errors == IcgemCoefficientErrors::calibrated || data.errors == IcgemCoefficientErrors::calibratedAndFormal ) ||
            coefficient.formalErrors.has_value( ) !=
                    ( data.errors == IcgemCoefficientErrors::formal || data.errors == IcgemCoefficientErrors::calibratedAndFormal ) )
            fail( data.source, coefficient.sourceLine, "Invalid coefficient indices, values or declared error columns" );
        for( const auto& sigma : { coefficient.calibratedErrors, coefficient.formalErrors } )
            if( sigma )
                for( double value : *sigma )
                    if( !std::isfinite( value ) || value < 0 )
                        fail( data.source, coefficient.sourceLine, "Invalid coefficient standard deviation" );
    }
}
IcgemGravityFieldData readIcgemGravityField( std::istream& stream, const std::string& source )
{
    IcgemGravityFieldData result;
    result.source = source;
    std::string line;
    std::size_t lineNumber = 0;
    bool headerEnded = false;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        if( !line.empty( ) && line.back( ) == '\r' ) line.pop_back( );
        const auto row = tokens( line );
        if( !headerEnded )
        {
            result.headerLines.push_back( line );
            if( !row.empty( ) && row.front( ) == "end_of_head" )
            {
                headerEnded = true;
                readHeader( result );
                validateIcgemGravityFieldData( result );
            }
            continue;
        }
        if( row.empty( ) ) continue;
        const auto& key = row.front( );
        if( key == "gfct" || key == "trnd" || key == "dot" || key == "asin" || key == "acos" )
            fail( source, lineNumber, "Temporal ICGEM records are unsupported by the static reader: " + key );
        if( key != "gfc" ) continue;
        if( row.size( ) < 5 + errorColumns( result.errors ) ) fail( source, lineNumber, "Truncated gfc record" );
        const auto indices = std::make_pair( integer( row[ 1 ], source, lineNumber ), integer( row[ 2 ], source, lineNumber ) );
        IcgemGravityCoefficient coefficient;
        coefficient.cosine = number( row[ 3 ], source, lineNumber );
        coefficient.sine = number( row[ 4 ], source, lineNumber );
        std::size_t column = 5;
        if( result.errors == IcgemCoefficientErrors::calibrated || result.errors == IcgemCoefficientErrors::calibratedAndFormal )
        {
            coefficient.calibratedErrors =
                    std::array< double, 2 >{ number( row[ column ], source, lineNumber ), number( row[ column + 1 ], source, lineNumber ) };
            column += 2;
        }
        if( result.errors == IcgemCoefficientErrors::formal || result.errors == IcgemCoefficientErrors::calibratedAndFormal )
            coefficient.formalErrors =
                    std::array< double, 2 >{ number( row[ column ], source, lineNumber ), number( row[ column + 1 ], source, lineNumber ) };
        coefficient.sourceLine = lineNumber;
        coefficient.originalRecord = line;
        if( !result.coefficients.emplace( indices, coefficient ).second ) fail( source, lineNumber, "Duplicate gfc degree/order" );
    }
    if( stream.bad( ) ) fail( source, lineNumber, "Input stream read failure" );
    if( !headerEnded ) fail( source, lineNumber, "Missing end_of_head" );
    validateIcgemGravityFieldData( result );
    return result;
}
IcgemGravityFieldData readIcgemGravityField( const std::string& fileName )
{
    std::ifstream stream( fileName );
    if( !stream ) throw std::runtime_error( fileName + ": Cannot open ICGEM gravity field" );
    return readIcgemGravityField( stream, fileName );
}
}  // namespace input_output
}  // namespace tudat
