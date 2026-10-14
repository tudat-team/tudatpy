#include "tudat/io/readOceanLoadingData.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace tudat
{
namespace input_output
{
namespace
{

std::string trimOceanLoadingText( const std::string& text )
{
    const auto first = text.find_first_not_of( " \t\r\n" );
    return first == std::string::npos ? "" : text.substr( first, text.find_last_not_of( " \t\r\n" ) - first + 1 );
}

std::string uppercaseOceanLoadingText( std::string text )
{
    std::transform( text.begin( ), text.end( ), text.begin( ), []( const unsigned char character ) { return std::toupper( character ); } );
    return text;
}

std::vector< std::string > tokenizeOceanLoadingText( const std::string& text )
{
    std::istringstream stream( text );
    stream.imbue( std::locale::classic( ) );
    std::vector< std::string > tokens;
    std::string token;
    while( stream >> token )
    {
        tokens.push_back( token );
    }
    return tokens;
}

[[noreturn]] void throwOceanLoadingParseError( const std::string& source, const std::size_t line, const std::string& message )
{
    throw std::runtime_error( "Invalid BLQ data in " + source + " at line " + std::to_string( line ) + ": " + message );
}

double parseOceanLoadingNumber( std::string token, const std::string& source, const std::size_t line )
{
    std::replace( token.begin( ), token.end( ), 'D', 'E' );
    std::replace( token.begin( ), token.end( ), 'd', 'e' );
    std::size_t consumed = 0;
    double value;
    try
    {
        value = std::stod( token, &consumed );
    }
    catch( const std::exception& )
    {
        throwOceanLoadingParseError( source, line, "invalid numerical coefficient '" + token + "'" );
    }
    if( consumed != token.size( ) || !std::isfinite( value ) )
    {
        throwOceanLoadingParseError( source, line, "nonfinite or invalid numerical coefficient '" + token + "'" );
    }
    return value;
}

void validateOceanLoadingHeaderRowOrder( const std::vector< std::string >& headers, const std::string& source )
{
    std::vector< char > declared;
    std::vector< char > sections;
    bool reading = false;
    bool inlineOrder = false;
    for( std::size_t line = 0; line < headers.size( ); ++line )
    {
        std::string text = uppercaseOceanLoadingText( trimOceanLoadingText( headers[ line ] ) );
        if( text.compare( 0, 2, "$$" ) == 0 )
        {
            text = trimOceanLoadingText( text.substr( 2 ) );
        }
        const auto order = text.find( "ROW ORDER:" );
        if( order != std::string::npos )
        {
            reading = true;
            text = trimOceanLoadingText( text.substr( order + 10 ) );
            if( !text.empty( ) )
            {
                inlineOrder = true;
                const auto tokens = tokenizeOceanLoadingText( text );
                for( const auto& token : tokens )
                {
                    if( token == "U" || token == "UP" || token == "RADIAL" )
                        declared.push_back( 'U' );
                    else if( token == "W" || token == "WEST" )
                        declared.push_back( 'W' );
                    else if( token == "S" || token == "SOUTH" )
                        declared.push_back( 'S' );
                    else
                        throwOceanLoadingParseError( source, line + 1, "unsupported declared BLQ row order" );
                }
            }
            continue;
        }
        if( !reading || inlineOrder || declared.size( ) >= 6 ) continue;
        if( text.find( "AMPLITUDES" ) == 0 )
            sections.push_back( 'A' );
        else if( text.find( "PHASES" ) == 0 )
            sections.push_back( 'P' );
        else if( text == "U" || text == "UP" || text == "VERTICAL" || text == "RADIAL" )
            declared.push_back( 'U' );
        else if( text == "W" || text == "WEST" || text.find( "TANGENTL" ) == 0 || text.find( "TANGENTIAL" ) == 0 )
        {
            if( text == "W" || text == "WEST" || text.find( "EW" ) != std::string::npos || text.find( "E-W" ) != std::string::npos )
                declared.push_back( 'W' );
            else if( text.find( "NS" ) != std::string::npos || text.find( "N-S" ) != std::string::npos )
                declared.push_back( 'S' );
            else
                throwOceanLoadingParseError( source, line + 1, "unsupported declared BLQ tangent row" );
        }
        else if( text == "S" || text == "SOUTH" )
            declared.push_back( 'S' );
        else if( text == "E" || text == "EAST" || text == "N" || text == "NORTH" )
            throwOceanLoadingParseError( source, line + 1, "BLQ rows must use up, west and south directions" );
    }
    if( reading )
    {
        const std::vector< char > full = { 'U', 'W', 'S', 'U', 'W', 'S' };
        const std::vector< char > shortOrder = { 'U', 'W', 'S' };
        if( declared != full && !( inlineOrder && declared == shortOrder ) )
            throwOceanLoadingParseError( source, 1, "BLQ row order must be up/west/south amplitudes followed by phases" );
        if( !sections.empty( ) && sections != std::vector< char >{ 'A', 'P' } )
            throwOceanLoadingParseError( source, 1, "BLQ amplitude rows must precede phase rows" );
    }
}

void readOceanLoadingHeaderMetadata( const std::string& header,
                                     OceanLoadingStation& metadata,
                                     const std::string& source,
                                     const std::size_t line )
{
    const std::string text = trimOceanLoadingText( header.substr( 2 ) );
    const std::string upper = uppercaseOceanLoadingText( text );
    const auto columnOrder = upper.find( "COLUMN ORDER:" );
    if( columnOrder != std::string::npos )
    {
        const auto names = tokenizeOceanLoadingText( upper.substr( columnOrder + 13 ) );
        const auto& expectedNames = getOceanLoadingConstituentNames( );
        if( names.size( ) != expectedNames.size( ) )
        {
            throwOceanLoadingParseError( source, line, "column order must contain the eleven standard BLQ constituents" );
        }
        for( std::size_t column = 0; column < names.size( ); ++column )
        {
            if( names[ column ] != uppercaseOceanLoadingText( expectedNames[ column ] ) )
            {
                throwOceanLoadingParseError( source, line, "unsupported BLQ constituent ordering" );
            }
        }
    }
    if( upper.find( "AMPLITUDES (" ) != std::string::npos && upper.find( "AMPLITUDES (M)" ) == std::string::npos )
    {
        throwOceanLoadingParseError( source, line, "BLQ amplitudes must be declared in metres" );
    }
    if( upper.find( "PHASES (" ) != std::string::npos && upper.find( "PHASES (DEGREES)" ) == std::string::npos )
    {
        throwOceanLoadingParseError( source, line, "BLQ phases must be declared in degrees" );
    }
    const auto cmc = upper.find( "CMC:" );
    const auto frame = upper.find( "FRAME:" );
    if( cmc != std::string::npos || frame != std::string::npos )
    {
        const auto values = tokenizeOceanLoadingText( upper.substr( cmc != std::string::npos ? cmc + 4 : frame + 6 ) );
        if( !values.empty( ) )
        {
            if( values.front( ) == "NO" || values.front( ) == "CE" )
            {
                metadata.referenceOrigin = OceanLoadingReferenceOrigin::solidEarthCentre;
            }
            else if( values.front( ) == "YES" || values.front( ) == "CM" )
            {
                metadata.referenceOrigin = OceanLoadingReferenceOrigin::wholeEarthCentreOfMass;
            }
            else
            {
                metadata.referenceOrigin = OceanLoadingReferenceOrigin::unspecified;
            }
        }
    }
    const auto model = upper.find( "OCEAN TIDE MODEL:" );
    if( model != std::string::npos )
    {
        metadata.oceanTideModel = trimOceanLoadingText( text.substr( model + 17 ) );
    }
    else if( upper.find( "COMPLETE " ) == 0 )
    {
        metadata.oceanTideModel = trimOceanLoadingText( text.substr( 9 ) );
    }
    else if( upper.find( " ID:" ) != std::string::npos )
    {
        metadata.oceanTideModel = trimOceanLoadingText( text.substr( 0, upper.find( " ID:" ) ) );
    }
    if( upper.find( "GREEN" ) != std::string::npos )
    {
        // Keep the complete declaration rather than guessing a provider-specific model identifier.
        metadata.greensFunction = text;
    }
    const auto coordinates = upper.find( "LON/LAT:" );
    if( coordinates != std::string::npos )
    {
        const auto values = tokenizeOceanLoadingText( text.substr( coordinates + 8 ) );
        if( values.size( ) >= 2 )
        {
            metadata.providerLongitudeLatitudeHeight( 0 ) = parseOceanLoadingNumber( values[ 0 ], source, line );
            metadata.providerLongitudeLatitudeHeight( 1 ) = parseOceanLoadingNumber( values[ 1 ], source, line );
            metadata.providerLongitudeLatitudeHeight( 2 ) = values.size( ) > 2 ? parseOceanLoadingNumber( values[ 2 ], source, line ) : 0.0;
            if( std::abs( metadata.providerLongitudeLatitudeHeight( 1 ) ) > 90.0 )
            {
                throwOceanLoadingParseError( source, line, "provider latitude is outside [-90,90] degrees" );
            }
            metadata.hasProviderCoordinates = true;
        }
    }
}

}  // namespace

const std::array< std::string, 11 >& getOceanLoadingConstituentNames( )
{
    static const std::array< std::string, 11 > names = { "M2", "S2", "N2", "K2", "K1", "O1", "P1", "Q1", "Mf", "Mm", "Ssa" };
    return names;
}

const OceanLoadingStation& OceanLoadingData::getStation( const std::string& stationIdentifier ) const
{
    const auto station = stations.find( stationIdentifier );
    if( station == stations.end( ) )
    {
        throw std::runtime_error( "Ocean loading station '" + stationIdentifier + "' is absent from " + source );
    }
    return station->second;
}

//! Check provider metadata and station coefficients, optionally permitting an initial empty station set.
void validateOceanLoadingData( const OceanLoadingData& data, const bool allowEmpty )
{
    if( data.stations.empty( ) && !allowEmpty )
    {
        throw std::runtime_error( "Ocean loading data contains no stations: " + data.source );
    }
    validateOceanLoadingHeaderRowOrder( data.sourceHeaders, data.source );
    for( const auto& item : data.stations )
    {
        const auto& station = item.second;
        validateOceanLoadingHeaderRowOrder( station.sourceHeaders, station.source );
        if( item.first.empty( ) || station.stationIdentifier != item.first )
        {
            throw std::runtime_error( "Ocean loading station identifier does not match its map key: " + item.first );
        }
        if( !station.amplitudesMetres.allFinite( ) || !station.phaseLagsDegrees.allFinite( ) ||
            ( station.amplitudesMetres.array( ) < 0.0 ).any( ) )
        {
            throw std::runtime_error( "Ocean loading station '" + item.first + "' has nonfinite coefficients or negative amplitudes" );
        }
        if( station.hasProviderCoordinates &&
            ( !station.providerLongitudeLatitudeHeight.allFinite( ) || std::abs( station.providerLongitudeLatitudeHeight( 1 ) ) > 90.0 ) )
        {
            throw std::runtime_error( "Ocean loading station '" + item.first + "' has invalid provider coordinates" );
        }
    }
}

//! Combine neutral BLQ records without changing station provenance or accepting duplicate site identities.
void appendOceanLoadingData( OceanLoadingData& destination, const OceanLoadingData& source )
{
    validateOceanLoadingData( source );
    for( const auto& item : source.stations )
    {
        if( destination.stations.count( item.first ) != 0 )
        {
            throw std::runtime_error( "Duplicate ocean loading station: " + item.first );
        }
    }

    // Preserve the existing model-append convention and every provider station's original source/header records.
    destination.source = destination.source.empty( ) ? source.source : "multiple BLQ inputs";
    destination.sourceHeaders.insert( destination.sourceHeaders.end( ), source.sourceHeaders.begin( ), source.sourceHeaders.end( ) );
    destination.stations.insert( source.stations.begin( ), source.stations.end( ) );
}

OceanLoadingData readOceanLoadingData( std::istream& stream, const std::string& source )
{
    OceanLoadingData data;
    data.source = source;
    OceanLoadingStation globalMetadata;
    globalMetadata.source = source;
    OceanLoadingStation currentStation;
    int row = 0;
    std::size_t lineNumber = 0;
    auto finishStation = [ & ] {
        if( currentStation.stationIdentifier.empty( ) )
        {
            return;
        }
        if( row != 6 )
        {
            throwOceanLoadingParseError(
                    source, lineNumber, "station '" + currentStation.stationIdentifier + "' requires six coefficient rows" );
        }
        if( !data.stations.emplace( currentStation.stationIdentifier, currentStation ).second )
        {
            throwOceanLoadingParseError( source, lineNumber, "duplicate station '" + currentStation.stationIdentifier + "'" );
        }
    };

    std::string line;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        const std::string text = trimOceanLoadingText( line );
        if( text.empty( ) )
        {
            continue;
        }
        if( text.compare( 0, 2, "$$" ) == 0 )
        {
            if( currentStation.stationIdentifier.empty( ) )
            {
                data.sourceHeaders.push_back( line );
                readOceanLoadingHeaderMetadata( text, globalMetadata, source, lineNumber );
            }
            else
            {
                currentStation.sourceHeaders.push_back( line );
                readOceanLoadingHeaderMetadata( text, currentStation, source, lineNumber );
            }
            continue;
        }
        const auto tokens = tokenizeOceanLoadingText( text );
        if( tokens.size( ) == 1 && ( currentStation.stationIdentifier.empty( ) || row == 6 ) )
        {
            finishStation( );
            currentStation = globalMetadata;
            currentStation.stationIdentifier = tokens.front( );
            currentStation.sourceHeaders.clear( );
            row = 0;
            continue;
        }
        if( currentStation.stationIdentifier.empty( ) )
        {
            throwOceanLoadingParseError( source, lineNumber, "coefficient row has no preceding station identifier" );
        }
        if( row >= 6 || tokens.size( ) != 11 )
        {
            throwOceanLoadingParseError(
                    source, lineNumber, "station '" + currentStation.stationIdentifier + "' requires exactly eleven coefficients per row" );
        }
        for( int column = 0; column < 11; ++column )
        {
            const double value = parseOceanLoadingNumber( tokens[ column ], source, lineNumber );
            if( row < 3 )
            {
                if( value < 0.0 )
                {
                    throwOceanLoadingParseError( source, lineNumber, "negative displacement amplitude" );
                }
                currentStation.amplitudesMetres( row, column ) = value;
            }
            else
            {
                currentStation.phaseLagsDegrees( row - 3, column ) = value;
            }
        }
        ++row;
    }
    if( stream.bad( ) || ( stream.fail( ) && !stream.eof( ) ) )
    {
        throwOceanLoadingParseError( source, lineNumber, "input stream failed" );
    }
    finishStation( );
    validateOceanLoadingData( data );
    return data;
}

OceanLoadingData readOceanLoadingData( const std::string& fileName )
{
    std::ifstream stream( fileName );
    if( !stream )
    {
        throw std::runtime_error( "Ocean loading BLQ file could not be opened: " + fileName );
    }
    return readOceanLoadingData( stream, fileName );
}

}  // namespace input_output
}  // namespace tudat
