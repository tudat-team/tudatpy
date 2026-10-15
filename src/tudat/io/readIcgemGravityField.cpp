/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/io/readIcgemGravityField.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <boost/date_time/gregorian/gregorian.hpp>

namespace tudat
{
namespace input_output
{
namespace
{
//! Attach source and one-based line context to format errors.
[[noreturn]] void fail( const std::string& source, std::size_t line, const std::string& message )
{
    throw std::runtime_error( source + ":" + std::to_string( line ) + ": " + message );
}

//! Split the format's whitespace-separated fields while retaining the original record elsewhere.
std::vector< std::string > tokens( const std::string& line )
{
    // The original line remains in the coefficient record; these tokens are temporary parse fields.
    std::istringstream stream( line );

    std::vector< std::string > result;

    std::string token;

    while( stream >> token )
    {
        result.push_back( token );
    }
    return result;
}

//! Parse a complete finite coefficient/error value, including Fortran D exponents.
double number( std::string token, const std::string& source, std::size_t line )
{
    std::replace( token.begin( ), token.end( ), 'D', 'E' );
    std::replace( token.begin( ), token.end( ), 'd', 'e' );
    try
    {
        // Require the complete field, rather than accepting a finite numeric prefix.
        std::size_t used = 0;

        const double value = std::stod( token, &used );

        if( used == token.size( ) && std::isfinite( value ) )
        {
            return value;
        }
    }
    catch( const std::exception& )
    {}
    fail( source, line, "Invalid finite number: " + token );
}

//! Reject signed, fractional or truncated degree/order fields.
int integer( const std::string& token, const std::string& source, std::size_t line )
{
    if( !token.empty( ) && std::all_of( token.begin( ), token.end( ), []( char c ) { return c >= '0' && c <= '9'; } ) )
    {
        try
        {
            return std::stoi( token );
        }
        catch( const std::exception& )
        {}
    }
    fail( source, line, "Invalid nonnegative integer: " + token );
}

//! Locate the date/period fields after the declared uncertainty columns.
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

//! Validate a published date without choosing its time scale or converting subday digits to seconds.
void validateDate( const std::string& date, const IcgemGravityFieldData& data, const IcgemGravityCoefficient& coefficient )
{
    // Fixed-width date tokens can be ordered exactly before a consuming time convention is selected.
    const bool piecewise = data.format == "icgem2.0";

    const bool digits = std::all_of( date.begin( ), date.end( ), []( const char c ) { return ( c >= '0' && c <= '9' ) || c == '.'; } );

    if( date.size( ) != ( piecewise ? 13 : 8 ) || !digits || ( piecewise ? date[ 8 ] != '.' : date.find( '.' ) != std::string::npos ) ||
        ( piecewise && ( date.find( '.' ) != 8 || date.find( '.', 9 ) != std::string::npos ) ) )
    {
        fail( data.source, coefficient.sourceLine, "Invalid temporal date token: " + date );
    }
    try
    {
        // Boost validates leap years and month lengths; no parallel calendar arithmetic is needed.
        const boost::gregorian::date calendar(
                std::stoi( date.substr( 0, 4 ) ), std::stoi( date.substr( 4, 2 ) ), std::stoi( date.substr( 6, 2 ) ) );

        (void)calendar;
    }
    catch( const std::exception& )
    {
        fail( data.source, coefficient.sourceLine, "Invalid Gregorian temporal date: " + date );
    }
}

//! Apply the same degree/order, finiteness and uncertainty checks to every record type.
void validateCoefficient( const IcgemGravityFieldData& data, const int degree, const int order, const IcgemGravityCoefficient& coefficient )
{
    if( degree < 0 || order < 0 || order > degree || degree > data.maximumDegree || !std::isfinite( coefficient.cosine ) ||
        !std::isfinite( coefficient.sine ) ||
        coefficient.calibratedErrors.has_value( ) !=
                ( data.errors == IcgemCoefficientErrors::calibrated || data.errors == IcgemCoefficientErrors::calibratedAndFormal ) ||
        coefficient.formalErrors.has_value( ) !=
                ( data.errors == IcgemCoefficientErrors::formal || data.errors == IcgemCoefficientErrors::calibratedAndFormal ) )
    {
        fail( data.source, coefficient.sourceLine, "Invalid coefficient indices, values or declared error columns" );
    }
    for( const auto& sigma : { coefficient.calibratedErrors, coefficient.formalErrors } )
    {
        if( sigma )
        {
            for( const double value : *sigma )
            {
                if( !std::isfinite( value ) || value < 0 )
                {
                    fail( data.source, coefficient.sourceLine, "Invalid coefficient standard deviation" );
                }
            }
        }
    }
}

//! Preserve independently dated terms while requiring their coefficient's reference model to cover them.
void validateTemporalCoefficients( const IcgemGravityFieldData& data )
{
    // Group only reference intervals for coverage checks; periodic records keep their own intervals.
    std::map< std::pair< int, int >, std::vector< const IcgemTemporalCoefficient* > > bases;

    std::set< std::tuple< std::string, int, int, std::string, std::string, double > > records;

    bool oldTrend = false, periodic = false;

    for( const auto& record : data.temporalCoefficients )
    {
        validateCoefficient( data, record.degree, record.order, record.coefficient );

        const bool base = record.type == "gfct";

        const bool harmonic = record.type == "asin" || record.type == "acos";

        const bool piecewise = data.format == "icgem2.0";

        if( ( !base && !harmonic && record.type != "trnd" && record.type != "dot" ) || ( record.type == "dot" && piecewise ) ||
            harmonic != record.periodYears.has_value( ) ||
            ( record.periodYears && ( !std::isfinite( *record.periodYears ) || *record.periodYears <= 0.0 ) ) )
        {
            fail( data.source, record.coefficient.sourceLine, "Invalid temporal keyword, format or period" );
        }
        if( base || piecewise )
        {
            validateDate( record.referenceDate, data, record.coefficient );
        }
        else if( !record.referenceDate.empty( ) )
        {
            fail( data.source, record.coefficient.sourceLine, "ICGEM 1.0 trend/amplitude has no printed reference date" );
        }
        if( piecewise )
        {
            validateDate( record.validityEndDate, data, record.coefficient );
            if( record.validityEndDate <= record.referenceDate )
            {
                fail( data.source, record.coefficient.sourceLine, "Temporal validity interval is empty or reversed" );
            }
        }
        else if( !record.validityEndDate.empty( ) )
        {
            fail( data.source, record.coefficient.sourceLine, "ICGEM 1.0 has no printed validity end date" );
        }
        if( !records.emplace( record.type == "dot" ? "trnd" : record.type,
                              record.degree,
                              record.order,
                              record.referenceDate,
                              record.validityEndDate,
                              record.periodYears.value_or( 0.0 ) )
                     .second )
        {
            fail( data.source, record.coefficient.sourceLine, "Duplicate temporal coefficient record" );
        }
        if( base )
        {
            bases[ { record.degree, record.order } ].push_back( &record );
        }
        oldTrend = oldTrend || record.type == "dot";
        periodic = periodic || harmonic;
    }
    if( oldTrend && periodic )
    {
        fail( data.source, 0, "Old dot records cannot be combined with asin/acos" );
    }
    for( auto& entry : bases )
    {
        std::sort( entry.second.begin( ), entry.second.end( ), []( const auto* a, const auto* b ) {
            return a->referenceDate < b->referenceDate;
        } );
        for( std::size_t index = 1; index < entry.second.size( ); ++index )
        {
            if( data.format == "icgem1.0" || entry.second[ index - 1 ]->validityEndDate > entry.second[ index ]->referenceDate )
            {
                fail( data.source, entry.second[ index ]->coefficient.sourceLine, "Multiple or overlapping gfct baselines" );
            }
        }
    }
    for( const auto& record : data.temporalCoefficients )
    {
        if( record.type == "gfct" )
        {
            continue;
        }

        const auto found = bases.find( { record.degree, record.order } );

        if( found == bases.end( ) )
        {
            fail( data.source, record.coefficient.sourceLine, "Temporal term requires a corresponding gfct coefficient" );
        }
        if( data.format == "icgem2.0" )
        {
            // A term may span several adjacent gfct pieces, as real EIGEN C20 annual terms do.
            std::string coveredUntil = record.referenceDate;

            for( const auto* base : found->second )
            {
                if( base->validityEndDate <= coveredUntil )
                {
                    continue;
                }
                if( base->referenceDate > coveredUntil )
                {
                    break;
                }
                coveredUntil = base->validityEndDate;
                if( coveredUntil >= record.validityEndDate )
                {
                    break;
                }
            }
            if( coveredUntil < record.validityEndDate )
            {
                fail( data.source, record.coefficient.sourceLine, "Temporal term extends beyond supplied gfct coverage" );
            }
        }
    }
}

//! Interpret the header after any explicit begin_of_head separator, preserving its preceding comments.
void readHeader( IcgemGravityFieldData& data )
{
    // Only recognized header keywords are interpreted; the complete original header remains available.
    const std::set< std::string > mandatory{ "product_type", "modelname", "earth_gravity_constant", "radius", "max_degree", "errors" };

    auto supported = mandatory;

    supported.insert( { "format", "norm", "tide_system" } );

    // The optional separator distinguishes a preamble from the effective header.
    std::set< std::string > seen;

    std::size_t start = 0;

    bool beginSeen = false;

    for( std::size_t i = 0; i < data.headerLines.size( ); ++i )
    {
        const auto row = tokens( data.headerLines[ i ] );

        if( !row.empty( ) && row.front( ) == "begin_of_head" )
        {
            if( beginSeen )
            {
                fail( data.source, i + 1, "Duplicate begin_of_head" );
            }
            beginSeen = true;
            start = i + 1;
        }
    }
    for( std::size_t i = start; i < data.headerLines.size( ); ++i )
    {
        const auto row = tokens( data.headerLines[ i ] );

        if( row.empty( ) || !supported.count( row.front( ) ) )
        {
            continue;
        }
        if( row.size( ) < 2 || !seen.insert( row.front( ) ).second )
        {
            fail( data.source, i + 1, "Missing value or duplicate header keyword: " + row.front( ) );
        }

        const auto& key = row.front( );

        const auto& value = row[ 1 ];

        if( key == "product_type" )
        {
            data.productType = value;
        }
        else if( key == "modelname" )
        {
            data.modelName = value;
        }
        else if( key == "earth_gravity_constant" )
        {
            data.gravitationalParameter = number( value, data.source, i + 1 );
        }
        else if( key == "radius" )
        {
            data.referenceRadius = number( value, data.source, i + 1 );
        }
        else if( key == "max_degree" )
        {
            data.maximumDegree = integer( value, data.source, i + 1 );
        }
        else if( key == "format" )
        {
            data.format = value;
        }
        else if( key == "errors" )
        {
            if( value == "no" )
            {
                data.errors = IcgemCoefficientErrors::none;
            }
            else if( value == "calibrated" )
            {
                data.errors = IcgemCoefficientErrors::calibrated;
            }
            else if( value == "formal" )
            {
                data.errors = IcgemCoefficientErrors::formal;
            }
            else if( value == "calibrated_and_formal" )
            {
                data.errors = IcgemCoefficientErrors::calibratedAndFormal;
            }
            else
            {
                fail( data.source, i + 1, "Unsupported errors convention: " + value );
            }
        }
        else if( key == "norm" )
        {
            if( value == "fully_normalized" )
            {
                data.normalization = IcgemNormalization::fullyNormalized;
            }
            else if( value == "unnormalized" )
            {
                data.normalization = IcgemNormalization::unnormalized;
            }
            else
            {
                fail( data.source, i + 1, "Unsupported normalization: " + value );
            }
        }
        else if( key == "tide_system" )
        {
            using basic_astrodynamics::TideSystem;

            if( value == "unknown" )
            {
                data.tideSystem = TideSystem::unknown;
            }
            else if( value == "tide_free" )
            {
                data.tideSystem = TideSystem::tideFree;
            }
            else if( value == "zero_tide" )
            {
                data.tideSystem = TideSystem::zeroTide;
            }
            else if( value == "mean_tide" )
            {
                data.tideSystem = TideSystem::meanTide;
            }
            else
            {
                fail( data.source, i + 1, "Unsupported tide system: " + value );
            }
        }
    }
    for( const auto& key : mandatory )
    {
        if( !seen.count( key ) )
        {
            fail( data.source, data.headerLines.size( ), "Missing mandatory header keyword: " + key );
        }
    }
}
}  // namespace

//! Validate the shared header and all supplied records without evaluating a temporal field.
void validateIcgemGravityFieldData( const IcgemGravityFieldData& data )
{
    using basic_astrodynamics::TideSystem;

    if( data.productType != "gravity_field" || data.modelName.empty( ) || data.maximumDegree < 0 ||
        !std::isfinite( data.gravitationalParameter ) || data.gravitationalParameter <= 0 || !std::isfinite( data.referenceRadius ) ||
        data.referenceRadius <= 0 || ( data.format != "icgem1.0" && data.format != "icgem2.0" ) ||
        ( data.normalization != IcgemNormalization::fullyNormalized && data.normalization != IcgemNormalization::unnormalized ) ||
        ( data.tideSystem != TideSystem::unknown && data.tideSystem != TideSystem::tideFree && data.tideSystem != TideSystem::zeroTide &&
          data.tideSystem != TideSystem::meanTide ) )
    {
        fail( data.source, 0, "Invalid ICGEM metadata" );
    }
    errorColumns( data.errors );
    for( const auto& entry : data.coefficients )
    {
        validateCoefficient( data, entry.first.first, entry.first.second, entry.second );
    }
    validateTemporalCoefficients( data );
}

//! Read every supported record through the same coefficient/error parser and preserve source order.
IcgemGravityFieldData readIcgemGravityField( std::istream& stream, const std::string& source )
{
    // Header and data records share one source line counter for exact diagnostics/provenance.
    IcgemGravityFieldData result;

    result.source = source;

    std::string line;

    std::size_t lineNumber = 0;

    bool headerEnded = false;

    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        if( !line.empty( ) && line.back( ) == '\r' )
        {
            line.pop_back( );
        }

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
        if( row.empty( ) )
        {
            continue;
        }

        const auto& key = row.front( );

        const bool temporal = key == "gfct" || key == "trnd" || key == "dot" || key == "asin" || key == "acos";

        if( key != "gfc" && !temporal )
        {
            continue;
        }

        const bool harmonic = key == "asin" || key == "acos";

        const std::size_t dateColumns = temporal ? ( result.format == "icgem2.0" ? 2 : ( key == "gfct" ? 1 : 0 ) ) : 0;

        if( row.size( ) < 5 + errorColumns( result.errors ) + dateColumns + ( harmonic ? 1 : 0 ) )
        {
            fail( source, lineNumber, "Truncated " + key + " record" );
        }
        const auto indices = std::make_pair( integer( row[ 1 ], source, lineNumber ), integer( row[ 2 ], source, lineNumber ) );

        // Static values, references, trends and amplitudes use the identical uncertainty columns.
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
        {
            coefficient.formalErrors =
                    std::array< double, 2 >{ number( row[ column ], source, lineNumber ), number( row[ column + 1 ], source, lineNumber ) };
            column += 2;
        }
        coefficient.sourceLine = lineNumber;
        coefficient.originalRecord = line;
        if( temporal )
        {
            IcgemTemporalCoefficient record;

            record.type = key;
            record.degree = indices.first;
            record.order = indices.second;
            record.coefficient = std::move( coefficient );
            if( dateColumns > 0 )
            {
                record.referenceDate = row[ column++ ];
            }
            if( dateColumns == 2 )
            {
                record.validityEndDate = row[ column++ ];
            }
            if( harmonic )
            {
                record.periodYears = number( row[ column ], source, lineNumber );
            }
            result.temporalCoefficients.push_back( std::move( record ) );
        }
        else if( !result.coefficients.emplace( indices, std::move( coefficient ) ).second )
        {
            fail( source, lineNumber, "Duplicate gfc degree/order" );
        }
    }
    if( stream.bad( ) )
    {
        fail( source, lineNumber, "Input stream read failure" );
    }
    if( !headerEnded )
    {
        fail( source, lineNumber, "Missing end_of_head" );
    }
    validateIcgemGravityFieldData( result );
    return result;
}

//! Read an existing file through the same stream implementation.
IcgemGravityFieldData readIcgemGravityField( const std::string& fileName )
{
    // All format handling stays in the stream reader, including temporal records and comments.
    std::ifstream stream( fileName );

    if( !stream )
    {
        throw std::runtime_error( fileName + ": Cannot open ICGEM gravity field" );
    }
    return readIcgemGravityField( stream, fileName );
}
}  // namespace input_output
}  // namespace tudat
