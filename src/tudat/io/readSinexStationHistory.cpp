/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/io/readSinexStationHistory.h"
#include "tudat/io/textParsing.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
constexpr double julianYearSeconds = 365.25 * 86400.0;
constexpr double radiansPerDegree = 0.017453292519943295769;

//! Line of a SINEX file that is being parsed, used for field extraction and error messages.
struct SinexLine {
    const std::string& source;
    std::size_t lineNumber;
    const std::string& text;

    [[noreturn]] void fail( const std::string& message ) const
    {
        throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": " + message );
    }

    std::string field( std::size_t start, std::size_t width, bool required = true ) const
    {
        if( text.size( ) < start + width )
        {
            fail( "Truncated fixed-width SINEX record" );
        }
        const auto result = text_parsing::trim( text.substr( start, width ) );
        if( required && result.empty( ) )
        {
            fail( "Missing field at column " + std::to_string( start + 1 ) );
        }
        return result;
    }

    int integer( const std::string& value ) const
    {
        const auto result = text_parsing::parseNonNegativeInteger( value );
        if( !result )
        {
            fail( "Invalid integer: " + value );
        }
        return *result;
    }

    double number( const std::string& value ) const
    {
        const auto result = text_parsing::parseFiniteNumber( value );
        if( !result )
        {
            fail( "Invalid/nonfinite numerical value: " + value );
        }
        return *result;
    }

    //! Epoch (YY:DDD:SSSSS) at the given column; an all-zero epoch has the given meaning.
    SinexEpoch epoch( std::size_t start, SinexEpochMeaning zeroMeaning ) const
    {
        const std::string token = field( start, 12 );
        if( token.size( ) != 12 || token[ 2 ] != ':' || token[ 6 ] != ':' )
        {
            fail( "Invalid SINEX epoch: " + token );
        }
        SinexEpoch result;
        result.original = token;
        if( token == "00:000:00000" )
        {
            result.meaning = zeroMeaning;
            return result;
        }
        const int year = integer( token.substr( 0, 2 ) );
        result.meaning = SinexEpochMeaning::specified;
        result.civilUtc =
                SinexCivilUtc{ year <= 50 ? 2000 + year : 1900 + year, integer( token.substr( 3, 3 ) ), integer( token.substr( 7, 5 ) ) };
        if( result.civilUtc->dayOfYear < 1 || result.civilUtc->dayOfYear > 366 || result.civilUtc->secondsOfDay > 86400 )
        {
            fail( "Invalid civil UTC epoch" );
        }
        return result;
    }

    SinexSolutionKey solutionKey( std::size_t offset = 0 ) const
    {
        return { { field( offset + 1, 4 ), field( offset + 6, 2 ) }, field( offset + 9, 4 ) };
    }
};

int cartesianComponent( const std::string& type )
{
    const std::array< std::string, 6 > types{ { "STAX", "STAY", "STAZ", "VELX", "VELY", "VELZ" } };
    const auto found = std::find( types.begin( ), types.end( ), type );
    return found == types.end( ) ? -1 : static_cast< int >( found - types.begin( ) );
}

//! Convert the value and standard deviation of a parameter to SI units, if its unit is supported.
void convertToSi( SinexParameterRecord& parameter, const SinexLine& line )
{
    // IGN publishes relaxation-time extension records with the literal unit "m"; their unit is not converted.
    if( parameter.type.rfind( "TEXP", 0 ) == 0 || parameter.type.rfind( "TLOG", 0 ) == 0 )
    {
        return;
    }
    static const std::map< std::string, std::pair< std::string, double > > units{ { "m", { "m", 1.0 } },
                                                                                  { "mm", { "m", 0.001 } },
                                                                                  { "m/y", { "m/s", 1.0 / julianYearSeconds } },
                                                                                  { "m/s", { "m/s", 1.0 } },
                                                                                  { "rad", { "rad", 1.0 } },
                                                                                  { "mas", { "rad", radiansPerDegree / 3600000.0 } },
                                                                                  { "y", { "s", julianYearSeconds } } };
    const auto unit = units.find( parameter.originalUnit );
    const int component = cartesianComponent( parameter.type );
    if( component >= 0 && ( unit == units.end( ) || unit->second.first != ( component < 3 ? "m" : "m/s" ) ) )
    {
        line.fail( "Invalid unit for Cartesian parameter " + parameter.type + ": " + parameter.originalUnit );
    }
    if( unit != units.end( ) )
    {
        parameter.siUnit = unit->second.first;
        parameter.valueSi = parameter.originalValue * unit->second.second;
        parameter.standardDeviationSi = parameter.originalStandardDeviation * unit->second.second;
    }
}

SinexSiteRecord siteRecord( const SinexLine& line )
{
    SinexSiteRecord site;
    site.key = { line.field( 1, 4 ), line.field( 6, 2 ) };
    site.domesId = line.field( 9, 9 );
    site.observationTechnique = line.field( 19, 1, false );
    site.description = line.field( 21, 22, false );

    // Approximate coordinates (degrees, minutes, seconds) are optional. F4.1 seconds can round to 60.0 (for example CEDU).
    if( line.text.size( ) >= 75 && !text_parsing::trim( line.text.substr( 44, 31 ) ).empty( ) )
    {
        // SLRF2020 repeats the sign of a southern latitude in the minutes and seconds.
        const auto coordinate = [ & ]( std::size_t start ) {
            const std::string degrees = line.field( start, 3 );
            const double value = std::fabs( line.number( degrees ) ) + std::fabs( line.number( line.field( start + 3, 3 ) ) ) / 60 +
                    std::fabs( line.number( line.field( start + 6, 5 ) ) ) / 3600;
            return ( degrees.front( ) == '-' ? -value : value ) * radiansPerDegree;
        };
        site.approximateLongitudeRadians = coordinate( 44 );
        site.approximateLatitudeRadians = coordinate( 56 );
        site.approximateHeightMetres = line.number( line.field( 68, 7 ) );
    }
    return site;
}

SinexParameterRecord parameterRecord( const SinexLine& line )
{
    SinexParameterRecord parameter;
    parameter.index = line.integer( line.field( 1, 5 ) );
    parameter.type = line.field( 7, 6 );
    parameter.key = line.solutionKey( 13 );
    parameter.referenceEpoch = line.epoch( 27, SinexEpochMeaning::unspecified );
    parameter.originalUnit = line.field( 40, 4 );
    parameter.constraintCode = line.field( 45, 1 );
    parameter.originalValue = line.number( line.field( 47, 21 ) );
    parameter.originalStandardDeviation = line.number( line.field( 69, 11 ) );
    convertToSi( parameter, line );
    return parameter;
}

std::string blockBase( const std::string& block )
{
    return block.substr( 0, block.find_first_of( " \t" ) );
}
}  // namespace

SinexStationHistory readSinexStationHistory( std::istream& stream, const std::string& source )
{
    SinexStationHistory history;
    history.source = source;
    std::string text, block;
    bool header = false, estimates = false;
    std::size_t lineNumber = 0;
    while( std::getline( stream, text ) )
    {
        ++lineNumber;
        if( !text.empty( ) && text.back( ) == '\r' )
        {
            text.pop_back( );
        }
        const SinexLine line{ source, lineNumber, text };
        const std::string stripped = text_parsing::trim( text );
        if( stripped.empty( ) || stripped.front( ) == '*' )
        {
            continue;
        }
        if( text.rfind( "%=SNX", 0 ) == 0 )
        {
            history.version = line.field( 6, 4 );
            history.agency = line.field( 11, 3 );
            history.creationEpoch = line.epoch( 15, SinexEpochMeaning::unspecified );
            history.headerStart = line.epoch( 32, SinexEpochMeaning::unspecified );
            history.headerEnd = line.epoch( 45, SinexEpochMeaning::unspecified );
            history.observationTechnique = line.field( 58, 1 );
            history.declaredParameterCount = line.integer( line.field( 60, 5 ) );
            header = true;
            continue;
        }
        if( !header )
        {
            line.fail( "Station history requires a SINEX header" );
        }
        if( stripped == "%ENDSNX" )
        {
            break;
        }
        // Block delimiters occupy column one; indented INPUT/HISTORY entries are data.
        if( text.front( ) == '+' || text.front( ) == '-' )
        {
            const std::string name = text_parsing::trim( stripped.substr( 1 ) );
            if( ( stripped.front( ) == '+' ) != block.empty( ) || ( stripped.front( ) == '-' && blockBase( name ) != blockBase( block ) ) )
            {
                line.fail( "Nested or mismatched SINEX block " + name );
            }
            block = stripped.front( ) == '+' ? name : "";
            estimates = estimates || blockBase( block ) == "SOLUTION/ESTIMATE";
            continue;
        }
        const auto name = blockBase( block );
        if( name == "SITE/ID" )
        {
            const auto site = siteRecord( line );
            history.sites.emplace( site.key, site );
        }
        else if( name == "SOLUTION/ESTIMATE" )
        {
            const auto parameter = parameterRecord( line );
            if( !history.parameters.emplace( parameter.index, parameter ).second )
            {
                line.fail( "Duplicate parameter index" );
            }
        }
        else if( name == "SOLUTION/EPOCHS" )
        {
            history.observationSpans.push_back( { line.solutionKey( ),
                                                  line.field( 14, 1 ),
                                                  line.epoch( 16, SinexEpochMeaning::headerStartDefault ),
                                                  line.epoch( 29, SinexEpochMeaning::headerEndDefault ),
                                                  line.epoch( 42, SinexEpochMeaning::unspecified ) } );
        }
        else if( name == "SOLUTION/DISCONTINUITY" )
        {
            SinexDiscontinuityRecord discontinuity;
            discontinuity.key = line.solutionKey( );
            discontinuity.observationTechnique = line.field( 14, 1 );
            discontinuity.start = line.epoch( 16, SinexEpochMeaning::openStart );
            discontinuity.end = line.epoch( 29, SinexEpochMeaning::openEnd );
            discontinuity.kind = line.field( 42, 1 );
            if( discontinuity.kind != "P" && discontinuity.kind != "V" )
            {
                line.fail( "Unsupported discontinuity kind" );
            }
            discontinuity.description = text.size( ) > 44 ? text_parsing::trim( text.substr( 44 ) ) : "";
            history.discontinuities.push_back( discontinuity );
        }
        else if( name == "SITE/ECCENTRICITY" )
        {
            SinexEccentricityRecord eccentricity;
            eccentricity.key = line.solutionKey( );
            eccentricity.observationTechnique = line.field( 14, 1 );
            eccentricity.start = line.epoch( 16, SinexEpochMeaning::headerStartDefault );
            eccentricity.end = line.epoch( 29, SinexEpochMeaning::headerEndDefault );
            const auto basis = line.field( 42, 3 );
            if( basis != "UNE" && basis != "XYZ" )
            {
                line.fail( "Unsupported eccentricity basis " + basis );
            }
            eccentricity.basis = basis == "UNE" ? SinexEccentricityBasis::une : SinexEccentricityBasis::xyz;
            for( int i = 0; i < 3; ++i )
            {
                eccentricity.componentsMetres[ i ] = line.number( line.field( 46 + i * 9, 8 ) );
            }
            history.eccentricities.push_back( eccentricity );
        }
        else if( name == "SITE/RECEIVER" || name == "SITE/ANTENNA" )
        {
            SinexEquipmentRecord equipment;
            equipment.kind = name == "SITE/RECEIVER" ? SinexEquipmentKind::receiver : SinexEquipmentKind::antenna;
            equipment.key = line.solutionKey( );
            equipment.observationTechnique = line.field( 14, 1 );
            equipment.start = line.epoch( 16, SinexEpochMeaning::headerStartDefault );
            equipment.end = line.epoch( 29, SinexEpochMeaning::headerEndDefault );
            equipment.equipmentType = line.field( 42, 20 );
            equipment.serialNumber = line.field( 63, 5, false );
            if( equipment.kind == SinexEquipmentKind::receiver )
            {
                equipment.firmwareVersion = line.field( 69, 11, false );
            }
            history.equipmentHistories.push_back( equipment );
        }
        else if( block.empty( ) )
        {
            line.fail( "Record outside a SINEX block" );
        }
    }
    if( !header || !block.empty( ) )
    {
        throw std::runtime_error( source + ": Empty or truncated SINEX file" );
    }
    if( estimates && history.parameters.size( ) != history.declaredParameterCount )
    {
        throw std::runtime_error( source + ": SOLUTION/ESTIMATE count differs from the declared parameter count" );
    }
    return history;
}

SinexStationHistory readSinexStationHistory( const std::string& fileName )
{
    std::ifstream stream( fileName );
    if( !stream )
    {
        throw std::runtime_error( fileName + ": Cannot open SINEX file" );
    }
    return readSinexStationHistory( stream, fileName );
}

std::vector< SinexSolutionKey > getSinexStationCoordinateSolutionKeys( const SinexStationHistory& history )
{
    std::set< SinexSolutionKey > keys;
    for( const auto& entry : history.parameters )
    {
        if( cartesianComponent( entry.second.type ) >= 0 )
        {
            keys.insert( entry.second.key );
        }
    }
    return { keys.begin( ), keys.end( ) };
}

SinexCartesianSolution getSinexStationCoordinateSolution( const SinexStationHistory& history,
                                                          const SinexSolutionKey& key,
                                                          bool requireCompleteVelocity )
{
    const auto fail = [ & ]( const std::string& message ) {
        throw std::runtime_error( history.source + ": " + message + " (" + key.site.siteCode + "/" + key.site.pointCode + "/" +
                                  key.solutionId + ")" );
    };
    SinexCartesianSolution solution;
    solution.key = key;
    const auto site = history.sites.find( key.site );
    if( site == history.sites.end( ) )
    {
        fail( "Requested solution has no SITE/ID identity" );
    }
    if( history.sites.count( key.site ) > 1 )
    {
        const auto matchingSites = history.sites.equal_range( key.site );
        for( auto candidate = matchingSites.first; candidate != matchingSites.second; ++candidate )
        {
            if( candidate->second.domesId != site->second.domesId )
            {
                fail( "Ambiguous SITE/ID monument for selected site/point" );
            }
        }
    }
    solution.site = site->second;
    std::array< bool, 6 > present{ { false, false, false, false, false, false } };
    for( const auto& entry : history.parameters )
    {
        const auto& parameter = entry.second;
        const int component = cartesianComponent( parameter.type );
        if( component < 0 || !( parameter.key == key ) )
        {
            continue;
        }
        if( present[ component ] )
        {
            fail( "Multiple reference epochs for Cartesian component " + parameter.type );
        }
        present[ component ] = true;
        if( !parameter.referenceEpoch.civilUtc )
        {
            fail( "Cartesian parameter reference epoch is unspecified" );
        }
        if( component < 3 )
        {
            if( solution.positionReferenceEpoch.civilUtc &&
                !( *solution.positionReferenceEpoch.civilUtc == *parameter.referenceEpoch.civilUtc ) )
            {
                fail( "Cartesian position components have different reference epochs" );
            }
            solution.positionReferenceEpoch = parameter.referenceEpoch;
            solution.positionMetres[ component ] = *parameter.valueSi;
            solution.positionParameterIndices[ component ] = parameter.index;
            solution.positionStandardDeviationsMetres[ component ] = *parameter.standardDeviationSi;
        }
        else
        {
            const int index = component - 3;
            solution.velocityMetresPerSecond[ index ] = parameter.valueSi;
            solution.velocityParameterIndices[ index ] = parameter.index;
            solution.velocityStandardDeviationsMetresPerSecond[ index ] = parameter.standardDeviationSi;
            solution.velocityReferenceEpochs[ index ] = parameter.referenceEpoch;
        }
    }
    for( int i = 0; i < ( requireCompleteVelocity ? 6 : 3 ); ++i )
    {
        if( !present[ i ] )
        {
            fail( "Selected Cartesian solution is incomplete" );
        }
    }
    return solution;
}
}  // namespace input_output
}  // namespace tudat
