#include "tudat/io/readTerrestrialFrameTransformationData.h"
#include "tudat/io/textParsing.h"
#include "tudat/paths.hpp"
#include <fstream>
#include <regex>
#include <stdexcept>

namespace tudat
{
namespace input_output
{

const TerrestrialFrameTransformationParameters& TerrestrialFrameTransformationData::getParameters( const std::string& target ) const
{
    const auto found = records.find( target );
    if( found == records.end( ) )
    {
        throw std::runtime_error( "No terrestrial-frame parameters for " + target );
    }
    return found->second;
}

TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( std::istream& stream, const std::string& source )
{
    constexpr double julianYear = 31557600.0;
    constexpr double masToRadians = 3.141592653589793238462643383279502884 / ( 180.0 * 3600.0 * 1000.0 );
    TerrestrialFrameTransformationData result;
    result.source = source;
    TerrestrialFrameTransformationParameters pending;
    std::string line;
    std::size_t count = 0;
    const auto number = [ & ]( const std::string& value ) {
        const auto parsed = text_parsing::parseFiniteNumber( value );
        if( !parsed )
        {
            throw std::runtime_error( source + ":" + std::to_string( count ) + ": invalid parameter " + value );
        }
        return *parsed;
    };

    // The table has a header naming the reference frame, then for each target frame a row of parameters (Tx, Ty, Tz in mm,
    // D in ppb, Rx, Ry, Rz in mas, epoch) followed by a row of their rates (per year).
    std::smatch match;
    while( std::getline( stream, line ) )
    {
        ++count;
        const auto row = text_parsing::splitIntoTokens( line );
        if( std::regex_search( line, match, std::regex( "Transformation parameters from (ITRF[0-9]+)" ) ) )
        {
            result.referenceFrame = { match[ 1 ], match[ 1 ] };
        }
        else if( row.size( ) == 9 && row[ 0 ].rfind( "ITRF", 0 ) == 0 )
        {
            pending = TerrestrialFrameTransformationParameters( );
            pending.sourceFrame = result.referenceFrame;
            pending.targetFrame.realization = row[ 0 ];
            for( int i = 0; i < 3; ++i )
            {
                pending.translationMetres[ i ] = number( row[ i + 1 ] ) * 1e-3;
                pending.rotationRadians[ i ] = number( row[ i + 5 ] ) * masToRadians;
            }
            pending.scale = number( row[ 4 ] ) * 1e-9;
            pending.referenceDecimalYear = number( row[ 8 ] );
        }
        else if( row.size( ) == 8 && row[ 0 ] == "rates" )
        {
            for( int i = 0; i < 3; ++i )
            {
                pending.translationRateMetresPerSecond[ i ] = number( row[ i + 1 ] ) * 1e-3 / julianYear;
                pending.rotationRateRadiansPerSecond[ i ] = number( row[ i + 5 ] ) * masToRadians / julianYear;
            }
            pending.scaleRatePerSecond = number( row[ 4 ] ) * 1e-9 / julianYear;
            if( pending.targetFrame.realization.empty( ) || !result.records.emplace( pending.targetFrame.realization, pending ).second )
            {
                throw std::runtime_error( source + ":" + std::to_string( count ) + ": orphan or duplicate parameter rates" );
            }
            pending.targetFrame.realization.clear( );
        }
    }
    if( result.referenceFrame.realization.empty( ) || result.records.empty( ) )
    {
        throw std::runtime_error( source + ": incomplete terrestrial-frame parameter table" );
    }
    return result;
}

TerrestrialFrameTransformationData readTerrestrialFrameTransformationData( const std::string& filePath )
{
    std::ifstream stream( filePath );
    if( !stream )
    {
        throw std::runtime_error( filePath + ": Cannot open terrestrial-frame parameter table" );
    }
    return readTerrestrialFrameTransformationData( stream, filePath );
}

TerrestrialFrameTransformationData readOfficialItrf2014TransformationData( )
{
    return readTerrestrialFrameTransformationData( paths::get_resources_path( ) + "/station_locations/Transfo-ITRF2014_ITRFs.txt" );
}

TerrestrialFrameTransformationData readOfficialItrf2020TransformationData( )
{
    return readTerrestrialFrameTransformationData( paths::get_resources_path( ) + "/station_locations/Transfo-ITRF2020_TRFs.txt" );
}
}  // namespace input_output
}  // namespace tudat
