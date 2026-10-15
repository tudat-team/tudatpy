#include "tudat/io/readTerrestrialFrameTransformationData.h"
#include "tudat/paths.hpp"
#include <cmath>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{
constexpr double julianYear = 31557600.0;
constexpr double masToRadians = 3.141592653589793238462643383279502884 / ( 180.0 * 3600.0 * 1000.0 );
std::vector< std::string > tokens( const std::string& line )
{
    std::istringstream input( line );
    std::vector< std::string > result;
    std::string token;
    while( input >> token ) result.push_back( token );
    return result;
}
std::runtime_error error( const std::string& source, std::size_t line, const std::string& message )
{
    return std::runtime_error( source + ":" + std::to_string( line ) + ": " + message );
}
double number( std::string value, const std::string& source, std::size_t line )
{
    for( auto& c : value )
        if( c == 'd' || c == 'D' ) c = 'E';
    try
    {
        std::size_t used = 0;
        const double result = std::stod( value, &used );
        if( used != value.size( ) || !std::isfinite( result ) ) throw std::invalid_argument( "nonfinite or trailing characters" );
        return result;
    }
    catch( const std::exception& )
    {
        throw error( source, line, "Invalid finite parameter: " + value );
    }
}
}  // namespace
const TerrestrialFrameTransformationParameters& TerrestrialFrameTransformationData::getParameters( const std::string& target ) const
{
    const auto found = records.find( target );
    if( found == records.end( ) ) throw std::runtime_error( "No terrestrial-frame parameters for " + target );
    return found->second;
}
void validateTerrestrialFrameTransformationParameters( const TerrestrialFrameTransformationParameters& p )
{
    if( p.sourceFrame.realization.empty( ) || p.targetFrame.realization.empty( ) || !std::isfinite( p.referenceDecimalYear ) ||
        !p.translationMetres.allFinite( ) || !p.rotationRadians.allFinite( ) || !std::isfinite( p.scale ) ||
        !p.translationRateMetresPerSecond.allFinite( ) || !p.rotationRateRadiansPerSecond.allFinite( ) ||
        !std::isfinite( p.scaleRatePerSecond ) || !std::isfinite( p.parameterRateYearSeconds ) || p.parameterRateYearSeconds <= 0 )
        throw std::runtime_error( "Invalid terrestrial-frame identities or finite SI parameters" );
    if( !p.provenance.sourceRelease.empty( ) && p.provenance.sourceRelease != p.sourceFrame.release )
        throw std::runtime_error( "Conflicting terrestrial-frame record release provenance" );
    if( p.rotationConvention != HelmertRotationConvention::unspecified &&
        p.rotationConvention != HelmertRotationConvention::positionVector &&
        p.rotationConvention != HelmertRotationConvention::coordinateFrame )
        throw std::runtime_error( "Unknown Helmert rotation convention" );
}
void validateTerrestrialFrameTransformationData( const TerrestrialFrameTransformationData& data )
{
    if( data.referenceFrame.realization.empty( ) || data.records.empty( ) )
        throw std::runtime_error( "Empty terrestrial-frame parameter catalogue" );
    if( !data.provenance.sourceRelease.empty( ) && data.provenance.sourceRelease != data.referenceFrame.release )
        throw std::runtime_error( "Conflicting terrestrial-frame catalogue release provenance" );
    for( const auto& entry : data.records )
    {
        validateTerrestrialFrameTransformationParameters( entry.second );
        if( entry.first != entry.second.targetFrame.realization ||
            entry.second.sourceFrame.realization != data.referenceFrame.realization ||
            entry.second.sourceFrame.release != data.referenceFrame.release )
            throw std::runtime_error( "Conflicting terrestrial-frame catalogue hub/identity" );
    }
}
TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( std::istream& stream,
                                                                           const TerrestrialFrameTransformationSource& supplied )
{
    TerrestrialFrameTransformationData result;
    result.provenance = supplied;
    const auto source = supplied.source.empty( ) ? std::string( "<stream>" ) : supplied.source;
    result.provenance.source = source;
    bool headline = false, columns = false, rateColumns = false, units = false, rateUnits = false, closed = false;
    bool xEquation = false, yEquation = false, zEquation = false;
    bool awaitingRates = false;
    TerrestrialFrameTransformationParameters pending;
    std::string line;
    std::size_t count = 0;
    while( std::getline( stream, line ) )
    {
        ++count;
        if( !line.empty( ) && line.back( ) == '\r' ) line.pop_back( );
        const auto row = tokens( line );
        if( row.empty( ) )
        {
            result.sourceHeaders.push_back( line );
            continue;
        }
        if( !headline )
        {
            std::smatch match;
            if( !std::regex_match( line, match, std::regex( "\\s*Transformation parameters from (ITRF[0-9]+) to past ITRFs\\.\\s*" ) ) )
                throw error( source, count, "Expected IGN transformation-table identity header" );
            result.referenceFrame = { match[ 1 ], supplied.sourceRelease.empty( ) ? match[ 1 ].str( ) : supplied.sourceRelease };
            headline = true;
            result.sourceHeaders.push_back( line );
            continue;
        }
        if( row[ 0 ] == "SOLUTION" )
        {
            if( columns || row != std::vector< std::string >{ "SOLUTION", "Tx", "Ty", "Tz", "D", "Rx", "Ry", "Rz", "EPOCH" } )
                throw error( source, count, "Conflicting or unsupported parameter column order" );
            columns = true;
        }
        else if( row[ 0 ] == "RATES" )
        {
            if( rateColumns || row != std::vector< std::string >{ "RATES", "Tx", "Ty", "Tz", "D", "Rx", "Ry", "Rz" } )
                throw error( source, count, "Conflicting or unsupported rate column order" );
            rateColumns = true;
        }
        else if( row[ 0 ] == "UNITS---------->" )
        {
            if( !units )
            {
                if( row != std::vector< std::string >{ "UNITS---------->", "mm", "mm", "mm", "ppb", ".001\"", ".001\"", ".001\"" } )
                    throw error( source, count, "Expected mm, ppb and mas parameter units" );
                units = true;
            }
            else
            {
                if( rateUnits ||
                    row !=
                            std::vector< std::string >{
                                    "UNITS---------->", "mm/y", "mm/y", "mm/y", "ppb/y", ".001\"/y", ".001\"/y", ".001\"/y" } )
                    throw error( source, count, "Expected published per-year parameter-rate units" );
                rateUnits = true;
            }
        }
        else if( row[ 0 ].rfind( "ITRF", 0 ) == 0 && !closed )
        {
            if( !columns || !rateColumns || !units || !rateUnits || awaitingRates || row.size( ) != 9 )
                throw error( source, count, "Incomplete declarations, malformed parameter row or missing paired rates" );
            pending = TerrestrialFrameTransformationParameters( );
            pending.sourceFrame = result.referenceFrame;
            pending.targetFrame.realization = row[ 0 ];
            pending.provenance = result.provenance;
            pending.rotationConvention = HelmertRotationConvention::positionVector;
            for( int i = 0; i < 3; ++i )
            {
                pending.translationMetres[ i ] = number( row[ i + 1 ], source, count ) * 1e-3;
                pending.rotationRadians[ i ] = number( row[ i + 5 ], source, count ) * masToRadians;
            }
            pending.scale = number( row[ 4 ], source, count ) * 1e-9;
            pending.referenceDecimalYear = number( row[ 8 ], source, count );
            pending.originalParameterLine = line;
            pending.parameterSourceLine = count;
            awaitingRates = true;
            continue;
        }
        else if( row[ 0 ] == "rates" && !closed )
        {
            if( !awaitingRates || row.size( ) != 8 ) throw error( source, count, "Orphan or malformed parameter-rate row" );
            for( int i = 0; i < 3; ++i )
            {
                pending.translationRateMetresPerSecond[ i ] = number( row[ i + 1 ], source, count ) * 1e-3 / julianYear;
                pending.rotationRateRadiansPerSecond[ i ] = number( row[ i + 5 ], source, count ) * masToRadians / julianYear;
            }
            pending.scaleRatePerSecond = number( row[ 4 ], source, count ) * 1e-9 / julianYear;
            pending.originalRateLine = line;
            pending.rateSourceLine = count;
            if( !result.records.emplace( pending.targetFrame.realization, pending ).second )
                throw error( source, count, "Duplicate target realization " + pending.targetFrame.realization );
            awaitingRates = false;
            continue;
        }
        else if( row[ 0 ][ 0 ] == '_' )
        {
            if( awaitingRates ) throw error( source, count, "Missing paired parameter rates" );
            closed = true;
        }
        else if( !closed && row[ 0 ][ 0 ] != '-' && row[ 0 ] != "." )
        {
            throw error( source, count, "Unexpected table record" );
        }
        // Confirm the published position-vector sign convention without using file-name inference.
        if( closed )
        {
            xEquation = xEquation || std::regex_search( line, std::regex( ":\\s*D\\s+-Rz\\s+Ry\\s*:" ) );
            yEquation = yEquation || std::regex_search( line, std::regex( ":\\s*Rz\\s+D\\s+-Rx\\s*:" ) );
            zEquation = zEquation || std::regex_search( line, std::regex( ":\\s*-Ry\\s+Rx\\s+D\\s*:" ) );
        }
        result.sourceHeaders.push_back( line );
    }
    if( stream.fail( ) && !stream.eof( ) ) throw error( source, count, "Input stream read failure" );
    if( !headline || !closed || awaitingRates || !xEquation || !yEquation || !zEquation )
        throw error( source, count, "Truncated table or missing explicit position-vector equation" );
    try
    {
        validateTerrestrialFrameTransformationData( result );
    }
    catch( const std::exception& e )
    {
        throw error( source, count, e.what( ) );
    }
    return result;
}
TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( const std::string& filePath,
                                                                           const TerrestrialFrameTransformationSource& supplied )
{
    std::ifstream stream( filePath );
    if( !stream ) throw std::runtime_error( filePath + ": Cannot open terrestrial-frame parameter table" );
    auto source = supplied;
    if( source.source.empty( ) || source.source == "<stream>" ) source.source = filePath;
    return readTerrestrialFrameTransformationData( stream, source );
}
//! Read an official IGN table from the station_locations resource directory, retaining its published source and checksum.
TerrestrialFrameTransformationData readOfficialItrfTransformationData( const std::string& fileName,
                                                                       const TerrestrialFrameTransformationSource& source )
{
    const std::string filePath = paths::get_resources_path( ) + "/station_locations/" + fileName;
    std::ifstream stream( filePath );
    if( !stream )
    {
        throw std::runtime_error( filePath + ": Cannot open official terrestrial-frame parameter table (tudat resources)" );
    }
    return readTerrestrialFrameTransformationData( stream, source );
}

TerrestrialFrameTransformationData readOfficialItrf2014TransformationData( )
{
    return readOfficialItrfTransformationData( "Transfo-ITRF2014_ITRFs.txt",
                                               { "https://itrf.ign.fr/docs/solutions/itrf2014/Transfo-ITRF2014_ITRFs.txt",
                                                 "ITRF2014",
                                                 "6cc1c068d0d7e9ce53277752795dc2ef3ea5f8890c4b5d1878e1e9ccc05a0543" } );
}

TerrestrialFrameTransformationData readOfficialItrf2020TransformationData( )
{
    return readOfficialItrfTransformationData( "Transfo-ITRF2020_TRFs.txt",
                                               { "https://itrf.ign.fr/docs/solutions/itrf2020/Transfo-ITRF2020_TRFs.txt",
                                                 "ITRF2020",
                                                 "0d8dc872e7ef994966063e614adb5b1b0fe20582359c6e77f1d83259d1ff1a68" } );
}
}  // namespace input_output
}  // namespace tudat
