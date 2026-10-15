/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/io/readSinexStationHistory.h"
#include "tudat/io/textParsing.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
constexpr double julianYearSeconds = 365.25 * 86400.0;
constexpr double radiansPerDegree = 0.017453292519943295769;
[[noreturn]] void failSinexRecord( const SinexSourceRecord& record, const std::string& message )
{
    throw std::runtime_error( record.source + ":" + std::to_string( record.line ) + ": " + message );
}
std::string sinexField( const SinexSourceRecord& record, std::size_t start, std::size_t width, bool required = true )
{
    if( record.text.size( ) < start + width ) failSinexRecord( record, "Truncated fixed-width SINEX record" );
    const auto result = text_parsing::trim( record.text.substr( start, width ) );
    if( required && result.empty( ) ) failSinexRecord( record, "Missing field at column " + std::to_string( start + 1 ) );
    return result;
}
int parseSinexInteger( const std::string& value, const SinexSourceRecord& record )
{
    const auto result = text_parsing::parseNonNegativeInteger( value );
    if( !result ) failSinexRecord( record, "Invalid integer: " + value );
    return *result;
}
double parseSinexNumber( const std::string& value, const SinexSourceRecord& record )
{
    const auto result = text_parsing::parseFiniteNumber( value );
    if( !result ) failSinexRecord( record, "Invalid/nonfinite numerical value: " + value );
    return *result;
}
bool leapYear( int year )
{
    return year % 4 == 0 && ( year % 100 != 0 || year % 400 == 0 );
}
void validateEpoch( const SinexEpoch& value, const SinexSourceRecord& record )
{
    switch( value.meaning )
    {
        case SinexEpochMeaning::specified:
        case SinexEpochMeaning::unspecified:
        case SinexEpochMeaning::headerStartDefault:
        case SinexEpochMeaning::headerEndDefault:
        case SinexEpochMeaning::openStart:
        case SinexEpochMeaning::openEnd:
            break;
        default:
            failSinexRecord( record, "Invalid epoch meaning" );
    }
    if( value.meaning == SinexEpochMeaning::specified )
    {
        if( !value.civilUtc ) failSinexRecord( record, "Specified epoch has no civil UTC fields" );
        const auto& date = *value.civilUtc;

        // Preserve a possible inserted-second label; physical leap-table validation belongs to a model.
        if( date.year < 1951 || date.year > 2050 || date.dayOfYear < 1 || date.dayOfYear > ( leapYear( date.year ) ? 366 : 365 ) ||
            date.secondsOfDay < 0 || date.secondsOfDay > 86400 )
            failSinexRecord( record, "Invalid civil UTC epoch" );
    }
    else if( value.civilUtc )
        failSinexRecord( record, "Unspecified/default/open epoch must not contain a civil date" );
    if( !value.original.empty( ) )
    {
        if( value.meaning == SinexEpochMeaning::specified )
        {
            if( value.original.size( ) != 12 || value.original[ 2 ] != ':' || value.original[ 6 ] != ':' )
                failSinexRecord( record, "Invalid original epoch token" );
            const int shortYear = parseSinexInteger( value.original.substr( 0, 2 ), record );
            const SinexCivilUtc original{ shortYear <= 50 ? 2000 + shortYear : 1900 + shortYear,
                                          parseSinexInteger( value.original.substr( 3, 3 ), record ),
                                          parseSinexInteger( value.original.substr( 7, 5 ), record ) };
            if( !( original == *value.civilUtc ) ) failSinexRecord( record, "Civil fields disagree with original epoch token" );
        }
        else if( value.original != "00:000:00000" )
            failSinexRecord( record, "Default/open epoch has a nonzero original token" );
    }
}
SinexEpoch epoch( const std::string& token, SinexEpochMeaning zeroMeaning, const SinexSourceRecord& record )
{
    SinexEpoch result;
    result.original = token;
    if( token.size( ) != 12 || token[ 2 ] != ':' || token[ 6 ] != ':' ) failSinexRecord( record, "Invalid SINEX epoch: " + token );
    if( token == "00:000:00000" )
        result.meaning = zeroMeaning;
    else
    {
        const int year = parseSinexInteger( token.substr( 0, 2 ), record );
        result.meaning = SinexEpochMeaning::specified;
        result.civilUtc = SinexCivilUtc{ year <= 50 ? 2000 + year : 1900 + year,
                                         parseSinexInteger( token.substr( 3, 3 ), record ),
                                         parseSinexInteger( token.substr( 7, 5 ), record ) };
    }
    validateEpoch( result, record );
    return result;
}
void validateInterval( const SinexEpoch& start, const SinexEpoch& end, const SinexSourceRecord& record, bool allowInstant = false )
{
    validateEpoch( start, record );
    validateEpoch( end, record );
    if( start.civilUtc && end.civilUtc && ( *end.civilUtc < *start.civilUtc || ( !allowInstant && *start.civilUtc == *end.civilUtc ) ) )
        failSinexRecord( record, "Empty/reversed station interval" );
}
SinexSolutionKey solutionKey( const SinexSourceRecord& record, bool estimate )
{
    const std::size_t offset = estimate ? 13 : 0;
    return { { sinexField( record, offset + 1, 4 ), sinexField( record, offset + 6, 2 ) }, sinexField( record, offset + 9, 4 ) };
}
int cartesianComponent( const std::string& type )
{
    const std::array< std::string, 6 > types{ { "STAX", "STAY", "STAZ", "VELX", "VELY", "VELZ" } };
    const auto found = std::find( types.begin( ), types.end( ), type );
    return found == types.end( ) ? -1 : static_cast< int >( found - types.begin( ) );
}
void normalize( SinexParameterRecord& parameter )
{
    const auto& unit = parameter.originalUnit;
    const auto& type = parameter.type;

    // IGN publishes relaxation-time extension records with the literal unit "m".
    // Only an explicitly declared provider profile may resolve that conflict.
    if( type.rfind( "TEXP", 0 ) == 0 || type.rfind( "TLOG", 0 ) == 0 )
    {
        parameter.normalization = SinexNormalization::unresolvedExtension;
        return;
    }
    double factor = 1.0;
    std::string siUnit;
    if( unit == "m" )
        siUnit = "m";
    else if( unit == "mm" )
    {
        siUnit = "m";
        factor = 0.001;
    }
    else if( unit == "m/y" )
    {
        siUnit = "m/s";
        factor = 1.0 / julianYearSeconds;
    }
    else if( unit == "m/s" )
        siUnit = "m/s";
    else if( unit == "rad" )
        siUnit = "rad";
    else if( unit == "mas" )
    {
        siUnit = "rad";
        factor = radiansPerDegree / 3600000.0;
    }
    else if( unit == "y" )
    {
        siUnit = "s";
        factor = julianYearSeconds;
    }
    const int component = cartesianComponent( type );
    if( component >= 0 && ( siUnit != ( component < 3 ? "m" : "m/s" ) ) )
        failSinexRecord( parameter.provenance, "Invalid unit for Cartesian parameter " + type + ": " + unit );
    if( !siUnit.empty( ) )
    {
        if( !std::isfinite( parameter.originalValue * factor ) ||
            !std::isfinite( parameter.originalStandardDeviation * std::fabs( factor ) ) )
            failSinexRecord( parameter.provenance, "SI normalization overflow" );
        parameter.normalization = SinexNormalization::si;
        parameter.siUnit = siUnit;
        parameter.valueSi = parameter.originalValue * factor;
        parameter.standardDeviationSi = parameter.originalStandardDeviation * std::fabs( factor );
    }
}
SinexSiteRecord siteRecord( const SinexSourceRecord& record )
{
    SinexSiteRecord site;
    site.key = { sinexField( record, 1, 4 ), sinexField( record, 6, 2 ) };
    site.domesId = sinexField( record, 9, 9 );
    site.observationTechnique = sinexField( record, 19, 1, false );
    site.description = sinexField( record, 21, 22, false );
    site.provenance = record;
    const auto coordinate = [ & ]( std::size_t start, bool latitude ) {
        const std::string degrees = sinexField( record, start, 3 );
        const double deg = parseSinexNumber( degrees, record );
        const double minutes = parseSinexNumber( sinexField( record, start + 4, 2 ), record );
        const double seconds = parseSinexNumber( sinexField( record, start + 7, 4 ), record );

        // F4.1 approximate coordinates in the IGS template can round to 60.0
        // seconds (for example CEDU). Addition carries that exact value into
        // the next minute while the original fields remain in provenance.
        if( minutes < 0 || minutes >= 60 || seconds < 0 || seconds > 60 ) failSinexRecord( record, "Invalid approximate DMS coordinate" );
        const double sign = degrees.front( ) == '-' ? -1 : 1;
        const double value = sign * ( std::fabs( deg ) + minutes / 60 + seconds / 3600 );
        if( latitude ? std::fabs( value ) > 90 : value < 0 || value >= 360 )
            failSinexRecord( record, "Approximate coordinate outside SINEX range" );
        return value * radiansPerDegree;
    };

    // SITE/ID may omit all approximate coordinates; partial fields are invalid.
    if( record.text.size( ) >= 75 && !text_parsing::trim( record.text.substr( 44, 31 ) ).empty( ) )
    {
        site.approximateLongitudeRadians = coordinate( 44, false );
        site.approximateLatitudeRadians = coordinate( 56, true );
        site.approximateHeightMetres = parseSinexNumber( sinexField( record, 68, 7 ), record );
    }
    else if( record.text.size( ) > 43 && !text_parsing::trim( record.text.substr( 43 ) ).empty( ) )
        failSinexRecord( record, "Truncated SITE/ID coordinates" );
    return site;
}
SinexParameterRecord parameterRecord( const SinexSourceRecord& record )
{
    SinexParameterRecord parameter;
    parameter.index = parseSinexInteger( sinexField( record, 1, 5 ), record );
    if( parameter.index == 0 ) failSinexRecord( record, "Parameter index must be positive" );
    parameter.type = sinexField( record, 7, 6 );
    parameter.key = solutionKey( record, true );
    parameter.referenceEpoch = epoch( sinexField( record, 27, 12 ), SinexEpochMeaning::unspecified, record );
    parameter.originalUnit = sinexField( record, 40, 4 );
    parameter.constraintCode = sinexField( record, 45, 1 );
    parameter.originalValue = parseSinexNumber( sinexField( record, 47, 21 ), record );
    parameter.originalStandardDeviation = parseSinexNumber( sinexField( record, 69, 11 ), record );
    parameter.provenance = record;
    if( parameter.constraintCode != "0" && parameter.constraintCode != "1" && parameter.constraintCode != "2" )
        failSinexRecord( record, "Invalid constraint code" );
    if( parameter.originalStandardDeviation < 0 ) failSinexRecord( record, "Negative parameter standard deviation" );
    normalize( parameter );
    return parameter;
}
std::string blockBase( const std::string& block )
{
    return block.substr( 0, block.find_first_of( " \t" ) );
}
bool knownBlock( const std::string& name )
{
    return name == "SITE/ID" || name == "SOLUTION/ESTIMATE" || name == "SOLUTION/EPOCHS" || name == "SOLUTION/DISCONTINUITY" ||
            name == "SITE/ECCENTRICITY" || name == "SITE/RECEIVER" || name == "SITE/ANTENNA";
}
}  // namespace

SinexStationHistory readSinexStationHistory( std::istream& stream, const std::string& source )
{
    SinexStationHistory history;
    history.source = source;
    std::string line, block;
    bool header = false, ended = false, estimates = false;
    std::size_t lineNumber = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        if( !line.empty( ) && line.back( ) == '\r' ) line.pop_back( );
        const SinexSourceRecord record{ source, lineNumber, line };
        const std::string stripped = text_parsing::trim( line );
        if( stripped.empty( ) ) continue;
        if( ended ) failSinexRecord( record, "Content after %ENDSNX" );
        if( stripped.front( ) == '*' )
        {
            history.comments.push_back( record );
            if( !block.empty( ) && !knownBlock( blockBase( block ) ) ) history.opaqueBlocks.back( ).records.push_back( record );
            continue;
        }
        if( line.rfind( "%=SNX", 0 ) == 0 )
        {
            if( header || !block.empty( ) ) failSinexRecord( record, "Duplicate/misplaced SINEX header" );
            history.originalHeader = line;
            history.version = sinexField( record, 6, 4 );
            if( history.version != "2.02" && history.version != "2.01" && history.version != "2.00" )
                failSinexRecord( record, "Unsupported SINEX version " + history.version );
            history.agency = sinexField( record, 11, 3 );
            history.creationEpoch = epoch( sinexField( record, 15, 12 ), SinexEpochMeaning::unspecified, record );
            history.headerStart = epoch( sinexField( record, 32, 12 ), SinexEpochMeaning::unspecified, record );
            history.headerEnd = epoch( sinexField( record, 45, 12 ), SinexEpochMeaning::unspecified, record );
            history.observationTechnique = sinexField( record, 58, 1 );
            history.declaredParameterCount = parseSinexInteger( sinexField( record, 60, 5 ), record );
            header = true;
            continue;
        }
        if( !header ) failSinexRecord( record, "Station history requires a SINEX header" );
        if( stripped == "%ENDSNX" )
        {
            if( !block.empty( ) ) failSinexRecord( record, "Unclosed SINEX block " + block );
            ended = true;
            continue;
        }
        if( stripped.front( ) == '+' )
        {
            if( !block.empty( ) ) failSinexRecord( record, "Nested SINEX blocks" );
            block = text_parsing::trim( stripped.substr( 1 ) );
            if( block.empty( ) ) failSinexRecord( record, "Missing SINEX block name" );
            if( blockBase( block ) == "SOLUTION/ESTIMATE" ) estimates = true;
            if( !knownBlock( blockBase( block ) ) ) history.opaqueBlocks.push_back( { block, {} } );
            continue;
        }
        if( stripped.front( ) == '-' )
        {
            if( block.empty( ) || blockBase( text_parsing::trim( stripped.substr( 1 ) ) ) != blockBase( block ) )
                failSinexRecord( record, "Mismatched SINEX block end" );
            block.clear( );
            continue;
        }
        if( block.empty( ) ) failSinexRecord( record, "Record outside a SINEX block" );
        const auto name = blockBase( block );
        if( name == "SITE/ID" )
        {
            const auto site = siteRecord( record );
            history.sites.emplace( site.key, site );
        }
        else if( name == "SOLUTION/ESTIMATE" )
        {
            const auto parameter = parameterRecord( record );
            if( !history.parameters.emplace( parameter.index, parameter ).second ) failSinexRecord( record, "Duplicate parameter index" );
        }
        else if( name == "SOLUTION/EPOCHS" )
        {
            SinexObservationSpan span;
            span.key = solutionKey( record, false );
            span.observationTechnique = sinexField( record, 14, 1 );
            span.start = epoch( sinexField( record, 16, 12 ), SinexEpochMeaning::headerStartDefault, record );
            span.end = epoch( sinexField( record, 29, 12 ), SinexEpochMeaning::headerEndDefault, record );
            span.mean = epoch( sinexField( record, 42, 12 ), SinexEpochMeaning::unspecified, record );
            span.provenance = record;
            history.observationSpans.push_back( span );
        }
        else if( name == "SOLUTION/DISCONTINUITY" )
        {
            SinexDiscontinuityRecord discontinuity;
            discontinuity.key = solutionKey( record, false );
            discontinuity.observationTechnique = sinexField( record, 14, 1 );
            discontinuity.start = epoch( sinexField( record, 16, 12 ), SinexEpochMeaning::openStart, record );
            discontinuity.end = epoch( sinexField( record, 29, 12 ), SinexEpochMeaning::openEnd, record );
            discontinuity.kind = sinexField( record, 42, 1 );
            if( discontinuity.kind != "P" && discontinuity.kind != "V" ) failSinexRecord( record, "Unsupported discontinuity kind" );
            discontinuity.description = line.size( ) > 44 ? text_parsing::trim( line.substr( 44 ) ) : "";
            discontinuity.provenance = record;
            history.discontinuities.push_back( discontinuity );
        }
        else if( name == "SITE/ECCENTRICITY" )
        {
            SinexEccentricityRecord eccentricity;
            eccentricity.key = solutionKey( record, false );
            eccentricity.observationTechnique = sinexField( record, 14, 1 );
            eccentricity.start = epoch( sinexField( record, 16, 12 ), SinexEpochMeaning::headerStartDefault, record );
            eccentricity.end = epoch( sinexField( record, 29, 12 ), SinexEpochMeaning::headerEndDefault, record );
            const auto basis = sinexField( record, 42, 3 );
            if( basis != "UNE" && basis != "XYZ" ) failSinexRecord( record, "Unsupported eccentricity basis " + basis );
            eccentricity.basis = basis == "UNE" ? SinexEccentricityBasis::une : SinexEccentricityBasis::xyz;
            for( int i = 0; i < 3; ++i )
                eccentricity.componentsMetres[ i ] = parseSinexNumber( sinexField( record, 46 + i * 9, 8 ), record );
            eccentricity.provenance = record;
            history.eccentricities.push_back( eccentricity );
        }
        else if( name == "SITE/RECEIVER" || name == "SITE/ANTENNA" )
        {
            SinexEquipmentRecord equipment;
            equipment.kind = name == "SITE/RECEIVER" ? SinexEquipmentKind::receiver : SinexEquipmentKind::antenna;
            equipment.key = solutionKey( record, false );
            equipment.observationTechnique = sinexField( record, 14, 1 );
            equipment.start = epoch( sinexField( record, 16, 12 ), SinexEpochMeaning::headerStartDefault, record );
            equipment.end = epoch( sinexField( record, 29, 12 ), SinexEpochMeaning::headerEndDefault, record );
            equipment.equipmentType = sinexField( record, 42, 20 );
            equipment.serialNumber = sinexField( record, 63, 5, false );
            if( equipment.kind == SinexEquipmentKind::receiver ) equipment.firmwareVersion = sinexField( record, 69, 11, false );
            equipment.provenance = record;
            history.equipmentHistories.push_back( equipment );
        }
        else
            history.opaqueBlocks.back( ).records.push_back( record );
    }
    if( stream.bad( ) ) throw std::runtime_error( source + ": I/O failure while reading SINEX" );
    if( !header || !ended || !block.empty( ) )
        throw std::runtime_error( source + ":" + std::to_string( lineNumber ) + ": Empty/truncated SINEX file" );
    if( estimates && history.parameters.size( ) != history.declaredParameterCount )
        failSinexRecord( { source, 1, history.originalHeader }, "SOLUTION/ESTIMATE count differs from the declared parameter count" );
    validateSinexStationHistory( history );
    return history;
}

SinexStationHistory readSinexStationHistory( const std::string& fileName )
{
    std::ifstream stream( fileName );
    if( !stream ) throw std::runtime_error( fileName + ": Cannot open SINEX file" );
    return readSinexStationHistory( stream, fileName );
}

void validateSinexStationHistory( const SinexStationHistory& history )
{
    if( history.originalHeader.empty( ) ) throw std::runtime_error( history.source + ": Missing source SINEX header" );
    if( history.version != "2.00" && history.version != "2.01" && history.version != "2.02" )
        throw std::runtime_error( history.source + ": Unsupported preloaded SINEX version" );
    validateEpoch( history.creationEpoch, { history.source, 1, history.originalHeader } );
    validateEpoch( history.headerStart, { history.source, 1, history.originalHeader } );
    validateEpoch( history.headerEnd, { history.source, 1, history.originalHeader } );
    validateInterval( history.headerStart, history.headerEnd, { history.source, 1, history.originalHeader }, true );
    std::set< std::tuple< std::string, std::string, std::string, std::string, std::string > > components;
    const auto validateRole = [ & ]( const SinexEpoch& value, SinexEpochMeaning zeroRole, const SinexSourceRecord& record ) {
        if( value.meaning != SinexEpochMeaning::specified && value.meaning != zeroRole )
            failSinexRecord( record, "Epoch meaning does not match its field role" );
    };
    const SinexSourceRecord headerRecord{ history.source, 1, history.originalHeader };
    validateRole( history.creationEpoch, SinexEpochMeaning::unspecified, headerRecord );
    validateRole( history.headerStart, SinexEpochMeaning::unspecified, headerRecord );
    validateRole( history.headerEnd, SinexEpochMeaning::unspecified, headerRecord );
    const auto validateKey = [ & ]( const SinexSolutionKey& key, const SinexSourceRecord& record ) {
        if( key.site.siteCode.empty( ) || key.site.siteCode.size( ) > 4 || key.site.pointCode.empty( ) || key.site.pointCode.size( ) > 2 ||
            key.solutionId.empty( ) || key.solutionId.size( ) > 4 )
            failSinexRecord( record, "Invalid station/point/solution identity" );
    };
    for( const auto& entry : history.sites )
    {
        const auto& site = entry.second;
        if( entry.first.siteCode != site.key.siteCode || entry.first.pointCode != site.key.pointCode || site.key.siteCode.empty( ) ||
            site.key.pointCode.empty( ) || site.domesId.empty( ) )
            failSinexRecord( site.provenance, "Invalid SITE/ID key/identity" );
        for( const auto& value : { site.approximateLongitudeRadians, site.approximateLatitudeRadians, site.approximateHeightMetres } )
            if( value && !std::isfinite( *value ) ) failSinexRecord( site.provenance, "Nonfinite site coordinate" );
        if( static_cast< bool >( site.approximateLongitudeRadians ) != static_cast< bool >( site.approximateLatitudeRadians ) ||
            static_cast< bool >( site.approximateLatitudeRadians ) != static_cast< bool >( site.approximateHeightMetres ) )
            failSinexRecord( site.provenance, "Partial approximate SITE/ID coordinate" );
        if( ( site.approximateLongitudeRadians &&
              ( *site.approximateLongitudeRadians < 0 || *site.approximateLongitudeRadians >= 360 * radiansPerDegree ) ) ||
            ( site.approximateLatitudeRadians && std::fabs( *site.approximateLatitudeRadians ) > 90 * radiansPerDegree ) )
            failSinexRecord( site.provenance, "SITE/ID coordinate outside allowed range" );
    }
    for( const auto& entry : history.parameters )
    {
        const auto& parameter = entry.second;
        if( entry.first != parameter.index || parameter.index == 0 || parameter.type.empty( ) || parameter.key.site.siteCode.empty( ) ||
            parameter.key.site.pointCode.empty( ) || parameter.key.solutionId.empty( ) )
            failSinexRecord( parameter.provenance, "Invalid parameter key/identity" );
        validateKey( parameter.key, parameter.provenance );
        if( parameter.constraintCode != "0" && parameter.constraintCode != "1" && parameter.constraintCode != "2" )
            failSinexRecord( parameter.provenance, "Invalid constraint code" );
        validateEpoch( parameter.referenceEpoch, parameter.provenance );
        validateRole( parameter.referenceEpoch, SinexEpochMeaning::unspecified, parameter.provenance );
        if( !std::isfinite( parameter.originalValue ) || !std::isfinite( parameter.originalStandardDeviation ) ||
            parameter.originalStandardDeviation < 0 )
            failSinexRecord( parameter.provenance, "Invalid parameter value/standard deviation" );
        SinexParameterRecord expected = parameter;
        expected.valueSi.reset( );
        expected.standardDeviationSi.reset( );
        expected.siUnit.clear( );
        expected.normalization = SinexNormalization::unsupportedUnit;
        normalize( expected );
        if( expected.normalization != parameter.normalization || expected.valueSi != parameter.valueSi ||
            expected.standardDeviationSi != parameter.standardDeviationSi || expected.siUnit != parameter.siUnit )
            failSinexRecord( parameter.provenance, "Invalid parameter normalization" );
        const auto& civil = parameter.referenceEpoch.civilUtc;
        const std::string epochIdentity = civil
                ? std::to_string( civil->year ) + ":" + std::to_string( civil->dayOfYear ) + ":" + std::to_string( civil->secondsOfDay )
                : "unspecified";
        const auto identity = std::make_tuple(
                parameter.key.site.siteCode, parameter.key.site.pointCode, parameter.key.solutionId, parameter.type, epochIdentity );
        if( cartesianComponent( parameter.type ) >= 0 && !components.insert( identity ).second )
            failSinexRecord( parameter.provenance, "Duplicate/conflicting Cartesian parameter component" );

        // SITE/ID is optional in ITRF discontinuity-only products, but Cartesian parameters need an identity.
        if( cartesianComponent( parameter.type ) >= 0 && history.sites.count( parameter.key.site ) == 0 )
            failSinexRecord( parameter.provenance, "Cartesian parameter has no matching SITE/ID point identity" );
    }
    for( const auto& span : history.observationSpans )
    {
        validateKey( span.key, span.provenance );
        validateInterval( span.start, span.end, span.provenance, true );
        validateRole( span.start, SinexEpochMeaning::headerStartDefault, span.provenance );
        validateRole( span.end, SinexEpochMeaning::headerEndDefault, span.provenance );
        validateRole( span.mean, SinexEpochMeaning::unspecified, span.provenance );
        validateEpoch( span.mean, span.provenance );
        if( span.mean.civilUtc &&
            ( ( span.start.civilUtc && *span.mean.civilUtc < *span.start.civilUtc ) ||
              ( span.end.civilUtc && *span.end.civilUtc < *span.mean.civilUtc ) ) )
            failSinexRecord( span.provenance, "Mean observation epoch outside observation span" );
    }
    for( const auto& discontinuity : history.discontinuities )
    {
        validateKey( discontinuity.key, discontinuity.provenance );
        validateInterval( discontinuity.start, discontinuity.end, discontinuity.provenance );
        validateRole( discontinuity.start, SinexEpochMeaning::openStart, discontinuity.provenance );
        validateRole( discontinuity.end, SinexEpochMeaning::openEnd, discontinuity.provenance );
        if( discontinuity.kind != "P" && discontinuity.kind != "V" )
            failSinexRecord( discontinuity.provenance, "Invalid discontinuity kind" );
    }
    for( const auto& eccentricity : history.eccentricities )
    {
        validateKey( eccentricity.key, eccentricity.provenance );
        if( eccentricity.basis != SinexEccentricityBasis::une && eccentricity.basis != SinexEccentricityBasis::xyz )
            failSinexRecord( eccentricity.provenance, "Invalid eccentricity component basis" );
        validateInterval( eccentricity.start, eccentricity.end, eccentricity.provenance );
        validateRole( eccentricity.start, SinexEpochMeaning::headerStartDefault, eccentricity.provenance );
        validateRole( eccentricity.end, SinexEpochMeaning::headerEndDefault, eccentricity.provenance );
        for( const double value : eccentricity.componentsMetres )
            if( !std::isfinite( value ) ) failSinexRecord( eccentricity.provenance, "Nonfinite eccentricity" );
    }
    for( const auto& equipment : history.equipmentHistories )
    {
        validateKey( equipment.key, equipment.provenance );
        if( equipment.kind != SinexEquipmentKind::receiver && equipment.kind != SinexEquipmentKind::antenna )
            failSinexRecord( equipment.provenance, "Invalid equipment kind" );
        if( equipment.equipmentType.empty( ) ) failSinexRecord( equipment.provenance, "Missing equipment type" );
        validateInterval( equipment.start, equipment.end, equipment.provenance );
        validateRole( equipment.start, SinexEpochMeaning::headerStartDefault, equipment.provenance );
        validateRole( equipment.end, SinexEpochMeaning::headerEndDefault, equipment.provenance );
    }
}

std::vector< SinexSolutionKey > getSinexCartesianSolutionKeys( const SinexStationHistory& history )
{
    std::set< SinexSolutionKey > keys;
    for( const auto& entry : history.parameters )
        if( cartesianComponent( entry.second.type ) >= 0 ) keys.insert( entry.second.key );
    return { keys.begin( ), keys.end( ) };
}

SinexCartesianSolution getSinexCartesianSolution( const SinexStationHistory& history,
                                                  const SinexSolutionKey& key,
                                                  bool requireCompleteVelocity )
{
    validateSinexStationHistory( history );
    SinexCartesianSolution solution;
    solution.key = key;
    const auto site = history.sites.find( key.site );
    if( site == history.sites.end( ) ) throw std::runtime_error( history.source + ": Requested solution has no SITE/ID identity" );
    const auto matchingSites = history.sites.equal_range( key.site );
    for( auto candidate = matchingSites.first; candidate != matchingSites.second; ++candidate )
    {
        const auto& first = site->second;
        const auto& other = candidate->second;
        if( first.domesId != other.domesId || first.observationTechnique != other.observationTechnique ||
            first.description != other.description || first.approximateLongitudeRadians != other.approximateLongitudeRadians ||
            first.approximateLatitudeRadians != other.approximateLatitudeRadians ||
            first.approximateHeightMetres != other.approximateHeightMetres )
            failSinexRecord(
                    other.provenance,
                    "Ambiguous SITE/ID monument metadata for selected site/point; supply a source product with an unambiguous DOMES "
                    "identity" );
    }
    solution.site = site->second;
    std::array< bool, 6 > present{ { false, false, false, false, false, false } };
    for( const auto& entry : history.parameters )
    {
        const auto& parameter = entry.second;
        if( parameter.key < key || key < parameter.key ) continue;
        const int component = cartesianComponent( parameter.type );
        if( component < 0 ) continue;
        if( present[ component ] )
            failSinexRecord( parameter.provenance,
                             "Multiple reference epochs for selected Cartesian component; select an explicit source solution" );
        present[ component ] = true;
        if( !parameter.referenceEpoch.civilUtc )
            failSinexRecord( parameter.provenance, "Cartesian parameter reference epoch is unspecified" );
        if( component < 3 )
        {
            // Parameter order does not determine the common epoch.
            if( solution.positionReferenceEpoch.civilUtc &&
                !( *solution.positionReferenceEpoch.civilUtc == *parameter.referenceEpoch.civilUtc ) )
                failSinexRecord( parameter.provenance, "Cartesian position components have different reference epochs" );
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
        if( !present[ i ] )
            throw std::runtime_error( history.source + ": Selected Cartesian solution is incomplete: " + key.site.siteCode + "/" +
                                      key.site.pointCode + "/" + key.solutionId );
    return solution;
}
}  // namespace input_output
}  // namespace tudat
