/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readPsfFile.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

namespace tudat
{

namespace input_output
{

namespace psf
{

Eigen::Vector2d RawPsfMeasurement::getEffectivePixelLine( ) const
{
    return observedPixelLine_ - localCorrection_;
}

namespace
{

std::string trimCopy( const std::string& s )
{
    return boost::algorithm::trim_copy( s );
}

bool startsWithCaseSensitive( const std::string& s, const std::string& prefix )
{
    return ( s.size( ) >= prefix.size( ) && s.compare( 0, prefix.size( ), prefix ) == 0 );
}

std::string stripSingleQuotesIfPresent( const std::string& s )
{
    std::string t = trimCopy( s );
    if( t.size( ) >= 2 && t.front( ) == '\'' && t.back( ) == '\'' )
    {
        return t.substr( 1, t.size( ) - 2 );
    }
    return t;
}

std::string stripTrailingCommaAndTrim( const std::string& s )
{
    std::string t = trimCopy( s );
    while( !t.empty( ) && t.back( ) == ',' )
    {
        t.pop_back( );
        t = trimCopy( t );
    }
    return t;
}

OpticalImageType parseOpticalImageTypeFromPsfToken( const std::string& raw )
{
    std::string t = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( raw ) );

    if( t == "STAR" )
    {
        return OpticalImageType::star;
    }
    if( t == "PLAN" )
    {
        return OpticalImageType::planet;
    }
    if( t == "SAT" )
    {
        return OpticalImageType::satellite;
    }
    if( t == "ROCK" )
    {
        return OpticalImageType::rock;
    }
    if( t == "END" )
    {
        return OpticalImageType::end_marker;
    }

    return OpticalImageType::unknown;
}

double toDouble( const std::string& s )
{
    return std::stod( trimCopy( s ) );
}
int toInt( const std::string& s )
{
    return std::stoi( trimCopy( s ) );
}
std::int64_t toInt64( const std::string& s )
{
    return static_cast< std::int64_t >( std::stoll( trimCopy( s ) ) );
}

std::string readNamelistBlockContents( std::istream& in, const std::string& firstLineWithDollarKeyword )
{
    std::string blockText;
    bool hasEndMarker = false;

    // Keep remainder of first line after "$XXX"
    {
        std::string remainder = firstLineWithDollarKeyword;
        const std::size_t firstSpace = remainder.find_first_of( " \t" );
        if( firstSpace != std::string::npos )
        {
            remainder = remainder.substr( firstSpace + 1 );
            remainder = trimCopy( remainder );
            if( !remainder.empty( ) )
            {
                blockText += " " + remainder;
            }
        }
    }

    std::string line;
    while( std::getline( in, line ) )
    {
        const std::string trimmed = trimCopy( line );
        if( trimmed.empty( ) )
        {
            continue;
        }

        const std::size_t endPos = trimmed.find( "$END" );
        if( endPos != std::string::npos )
        {
            const std::string beforeEnd = trimCopy( trimmed.substr( 0, endPos ) );
            if( !beforeEnd.empty( ) )
            {
                blockText += " " + beforeEnd;
            }
            hasEndMarker = true;
            break;
        }

        blockText += " " + trimmed;
    }

    if( !hasEndMarker )
    {
        throw std::runtime_error( "Error when reading PSF file: namelist block starting with '" + firstLineWithDollarKeyword +
                                  "' has no $END marker." );
    }

    return blockText;
}

std::vector< std::pair< std::string, std::string > > parseFortranNamelistAssignments( const std::string& block )
{
    std::vector< std::pair< std::string, std::string > > assignments;

    const std::size_t n = block.size( );
    std::size_t pos = 0;

    auto skipDelims = [ & ]( ) {
        while( pos < n )
        {
            const char c = block[ pos ];
            if( std::isspace( static_cast< unsigned char >( c ) ) || c == ',' )
            {
                ++pos;
            }
            else
            {
                break;
            }
        }
    };

    auto isKeyStart = [ & ]( std::size_t p ) -> bool {
        while( p < n && std::isspace( static_cast< unsigned char >( block[ p ] ) ) )
        {
            ++p;
        }
        if( p >= n )
        {
            return false;
        }
        const char c = block[ p ];
        return ( std::isalpha( static_cast< unsigned char >( c ) ) || c == '_' );
    };

    auto looksLikeKeyEqualsAt = [ & ]( std::size_t p ) -> bool {
        if( !isKeyStart( p ) )
        {
            return false;
        }

        bool inQuotes = false;
        int parenDepth = 0;

        // Scan forward to see if we hit an '=' before hitting a top-level comma
        for( std::size_t i = p; i < n; ++i )
        {
            const char c = block[ i ];

            if( c == '\'' )
            {
                inQuotes = !inQuotes;
                continue;
            }

            if( inQuotes )
            {
                continue;
            }

            if( c == '(' )
            {
                ++parenDepth;
                continue;
            }
            if( c == ')' )
            {
                if( parenDepth > 0 )
                {
                    --parenDepth;
                }
                continue;
            }

            if( c == '=' )
            {
                return true;
            }

            // Only a comma at top-level (not inside parentheses) terminates the key candidate
            if( c == ',' && parenDepth == 0 )
            {
                return false;
            }
        }

        return false;
    };

    while( true )
    {
        skipDelims( );
        if( pos >= n )
        {
            break;
        }

        // Find '=' for current assignment (must be outside quotes)
        bool inQuotes = false;
        std::size_t eqPos = std::string::npos;
        for( std::size_t i = pos; i < n; ++i )
        {
            const char c = block[ i ];
            if( c == '\'' )
            {
                inQuotes = !inQuotes;
                continue;
            }
            if( !inQuotes && c == '=' )
            {
                eqPos = i;
                break;
            }
        }

        if( eqPos == std::string::npos )
        {
            break;
        }

        const std::string key = trimCopy( block.substr( pos, eqPos - pos ) );
        pos = eqPos + 1;

        // Value runs until a comma that starts the next KEY= (outside quotes)
        const std::size_t valueStart = pos;

        inQuotes = false;
        std::size_t valueEnd = n;

        for( std::size_t i = pos; i < n; ++i )
        {
            const char c = block[ i ];
            if( c == '\'' )
            {
                inQuotes = !inQuotes;
                continue;
            }

            if( !inQuotes && c == ',' )
            {
                if( looksLikeKeyEqualsAt( i + 1 ) )
                {
                    valueEnd = i;
                    pos = i + 1;
                    break;
                }
            }
        }

        if( valueEnd == n )
        {
            pos = n;
        }

        const std::string value = trimCopy( block.substr( valueStart, valueEnd - valueStart ) );
        if( !key.empty( ) )
        {
            assignments.push_back( std::make_pair( key, value ) );
        }

        if( pos >= n )
        {
            break;
        }
    }

    return assignments;
}

bool parseIndexedKey( const std::string& keyWithOptionalIndices, std::string& baseKeyOut, std::vector< int >& indicesOut )
{
    baseKeyOut.clear( );
    indicesOut.clear( );

    const std::size_t parenOpen = keyWithOptionalIndices.find( '(' );
    if( parenOpen == std::string::npos )
    {
        baseKeyOut = trimCopy( keyWithOptionalIndices );
        return true;
    }

    const std::size_t parenClose = keyWithOptionalIndices.find( ')', parenOpen );
    if( parenClose == std::string::npos )
    {
        return false;
    }

    baseKeyOut = trimCopy( keyWithOptionalIndices.substr( 0, parenOpen ) );
    std::string inside = keyWithOptionalIndices.substr( parenOpen + 1, parenClose - parenOpen - 1 );
    boost::algorithm::trim( inside );

    std::vector< std::string > parts;
    boost::algorithm::split( parts, inside, boost::algorithm::is_any_of( "," ), boost::algorithm::token_compress_on );

    for( std::size_t i = 0; i < parts.size( ); ++i )
    {
        const std::string p = trimCopy( parts[ i ] );
        if( !p.empty( ) )
        {
            indicesOut.push_back( std::stoi( p ) );
        }
    }

    return true;
}

std::vector< std::string > splitCommaSeparatedRespectingQuotes( const std::string& text )
{
    std::vector< std::string > tokens;
    std::string current;
    bool inSingleQuotes = false;

    for( std::size_t i = 0; i < text.size( ); ++i )
    {
        const char c = text[ i ];
        if( c == '\'' )
        {
            inSingleQuotes = !inSingleQuotes;
            current.push_back( c );
        }
        else if( c == ',' && !inSingleQuotes )
        {
            const std::string t = trimCopy( current );
            if( !t.empty( ) )
            {
                tokens.push_back( t );
            }
            current.clear( );
        }
        else
        {
            current.push_back( c );
        }
    }

    const std::string t = trimCopy( current );
    if( !t.empty( ) )
    {
        tokens.push_back( t );
    }
    return tokens;
}

Eigen::Vector2d parseVector2FromValueText( const std::string& valueText, const std::string& context )
{
    std::vector< std::string > vals = splitCommaSeparatedRespectingQuotes( valueText );
    if( vals.size( ) < 2 )
    {
        throw std::runtime_error( "Error when reading PSF file: malformed 2-vector for " + context + ", value='" + valueText + "'" );
    }

    Eigen::Vector2d v;
    v( 0 ) = toDouble( vals.at( 0 ) );
    v( 1 ) = toDouble( vals.at( 1 ) );
    return v;
}

void throwMissingPsfFieldError( const std::string& blockName, const std::string& fieldName, const std::string& context = "" )
{
    throw std::runtime_error( "Error when reading PSF file: missing required " + blockName + " field " + fieldName +
                              ( context.empty( ) ? "." : " for " + context + "." ) );
}

}  // namespace

RawPsfFileContents readPsfFile( const std::string& psfFile )
{
    std::ifstream dataFile( psfFile.c_str( ) );
    if( !dataFile.good( ) )
    {
        throw std::runtime_error( "Error when opening PSF file, file " + psfFile + " could not be opened." );
    }

    RawPsfFileContents fileContents;

    bool idParsed = false;
    bool camParsed = false;

    bool hasPendingLine = false;
    std::string pendingLine;

    std::string line;
    while( true )
    {
        if( hasPendingLine )
        {
            line = pendingLine;
            hasPendingLine = false;
        }
        else
        {
            if( !std::getline( dataFile, line ) )
            {
                break;
            }
        }

        const std::string trimmedLine = trimCopy( line );
        if( trimmedLine.empty( ) )
        {
            continue;
        }

        // -----------------------------
        // $ID
        // -----------------------------
        if( startsWithCaseSensitive( trimmedLine, "$ID" ) )
        {
            const std::string block = readNamelistBlockContents( dataFile, trimmedLine );
            const std::vector< std::pair< std::string, std::string > > assigns = parseFortranNamelistAssignments( block );

            for( std::size_t i = 0; i < assigns.size( ); ++i )
            {
                const std::string key = trimCopy( assigns.at( i ).first );
                const std::string value = trimCopy( assigns.at( i ).second );

                if( key == "SCID" )
                {
                    fileContents.spacecraftId_ = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( value ) );
                }
                else if( key == "NCAM" )
                {
                    fileContents.numberOfCameras_ = toInt( value );
                }
                else if( key == "EQUNOX" )
                {
                    fileContents.equinox_ = toInt( value );
                }
                else if( key == "PSFID" )
                {
                    fileContents.psfId_ = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( value ) );
                }
                else if( key == "PSFPRG" )
                {
                    fileContents.psfProgram_ = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( value ) );
                }
                else if( key == "PSFTIM" )
                {
                    fileContents.psfGenerationTimeUtcString_ = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( value ) );
                }
                else if( key == "PSFCOM" )
                {
                    // Value is a comma-separated list of quoted strings
                    const std::vector< std::string > parts = splitCommaSeparatedRespectingQuotes( value );
                    for( std::size_t j = 0; j < parts.size( ); ++j )
                    {
                        const std::string comment = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( parts.at( j ) ) );
                        if( !comment.empty( ) )
                        {
                            fileContents.psfComments_.push_back( comment );
                        }
                    }
                }
            }

            idParsed = true;
            continue;
        }

        // -----------------------------
        // $CAM
        // -----------------------------
        if( startsWithCaseSensitive( trimmedLine, "$CAM" ) )
        {
            const std::string block = readNamelistBlockContents( dataFile, trimmedLine );
            const std::vector< std::pair< std::string, std::string > > assigns = parseFortranNamelistAssignments( block );

            std::map< int, RawPsfCameraProperties > camerasByIndex;

            for( std::size_t i = 0; i < assigns.size( ); ++i )
            {
                const std::string rawKey = trimCopy( assigns.at( i ).first );
                const std::string rawValue = trimCopy( assigns.at( i ).second );

                std::string baseKey;
                std::vector< int > indices;
                if( !parseIndexedKey( rawKey, baseKey, indices ) )
                {
                    continue;
                }

                if( baseKey == "CAMID" && indices.size( ) == 1 )
                {
                    const int icam = indices.at( 0 );
                    camerasByIndex[ icam ].cameraId_ = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( rawValue ) );
                }
                else if( baseKey == "FL" && indices.size( ) == 1 )
                {
                    const int icam = indices.at( 0 );
                    camerasByIndex[ icam ].focalLengthMm_ = toDouble( rawValue );
                }
                else if( baseKey == "PLCTR" && indices.size( ) == 2 )
                {
                    const int icam = indices.at( 1 );
                    const Eigen::Vector2d v = parseVector2FromValueText( rawValue, "PLCTR" );
                    camerasByIndex[ icam ].principalPoint_ = v;
                }
                else if( baseKey == "PLSIZ" && indices.size( ) == 2 )
                {
                    const int icam = indices.at( 1 );
                    std::vector< std::string > vals = splitCommaSeparatedRespectingQuotes( rawValue );
                    if( vals.size( ) < 4 )
                    {
                        throw std::runtime_error( "Error when reading PSF file: malformed PLSIZ entry, value='" + rawValue + "'" );
                    }
                    camerasByIndex[ icam ].fieldOfViewBounds_( 0 ) = toDouble( vals.at( 0 ) );
                    camerasByIndex[ icam ].fieldOfViewBounds_( 1 ) = toDouble( vals.at( 1 ) );
                    camerasByIndex[ icam ].fieldOfViewBounds_( 2 ) = toDouble( vals.at( 2 ) );
                    camerasByIndex[ icam ].fieldOfViewBounds_( 3 ) = toDouble( vals.at( 3 ) );
                }
                else if( baseKey == "KMAT" && indices.size( ) == 3 )
                {
                    // KMAT(1,j,icam) = v1, v2  -> K(1,j)=v1 and K(2,j)=v2 (1-based -> 0-based)
                    const int iIndex = indices.at( 0 );
                    const int jIndex = indices.at( 1 );
                    const int icam = indices.at( 2 );

                    if( iIndex != 1 )
                    {
                        continue;
                    }

                    const int col = jIndex - 1;
                    if( col < 0 || col >= 3 )
                    {
                        continue;
                    }

                    const Eigen::Vector2d v = parseVector2FromValueText( rawValue, "KMAT" );
                    camerasByIndex[ icam ].kMatrix_( 0, col ) = v( 0 );
                    camerasByIndex[ icam ].kMatrix_( 1, col ) = v( 1 );
                }
                else if( baseKey == "EM" && indices.size( ) == 2 )
                {
                    const int startIndex = indices.at( 0 );  // 1 or 4
                    const int icam = indices.at( 1 );
                    const int start0 = startIndex - 1;

                    std::vector< std::string > vals = splitCommaSeparatedRespectingQuotes( rawValue );
                    for( std::size_t j = 0; j < vals.size( ); ++j )
                    {
                        const int idx = start0 + static_cast< int >( j );
                        if( idx >= 0 && idx < 6 )
                        {
                            camerasByIndex[ icam ].distortionCoefficients_( idx ) = toDouble( vals.at( j ) );
                        }
                    }
                }
                else if( baseKey == "OFFSET" && indices.size( ) == 2 )
                {
                    const int icam = indices.at( 1 );
                    std::vector< std::string > vals = splitCommaSeparatedRespectingQuotes( rawValue );
                    if( vals.size( ) < 3 )
                    {
                        throw std::runtime_error( "Error when reading PSF file: malformed OFFSET entry, value='" + rawValue + "'" );
                    }
                    camerasByIndex[ icam ].mountingOffsetsDegrees_( 0 ) = toDouble( vals.at( 0 ) );
                    camerasByIndex[ icam ].mountingOffsetsDegrees_( 1 ) = toDouble( vals.at( 1 ) );
                    camerasByIndex[ icam ].mountingOffsetsDegrees_( 2 ) = toDouble( vals.at( 2 ) );
                }
            }

            for( std::map< int, RawPsfCameraProperties >::const_iterator it = camerasByIndex.begin( ); it != camerasByIndex.end( ); ++it )
            {
                const RawPsfCameraProperties& cam = it->second;
                if( !cam.cameraId_.empty( ) )
                {
                    fileContents.cameraProperties_[ cam.cameraId_ ] = cam;
                }
            }

            camParsed = true;
            continue;
        }

        // -----------------------------
        // $PIC
        // -----------------------------
        if( startsWithCaseSensitive( trimmedLine, "$PIC" ) )
        {
            if( !idParsed || !camParsed )
            {
                throw std::runtime_error( "Error when reading PSF file: encountered $PIC before parsing $ID/$CAM." );
            }

            RawPsfFileImageContents imageContents;
            bool hasPictureName = false;
            bool hasPictureNumber = false;
            bool hasEndOfExposureTime = false;
            bool hasCameraId = false;
            bool hasExposureTime = false;
            bool hasPictureDeletionFlag = false;
            bool hasRightAscension = false;
            bool hasDeclination = false;
            bool hasTwist = false;

            const std::string block = readNamelistBlockContents( dataFile, trimmedLine );
            const std::vector< std::pair< std::string, std::string > > assigns = parseFortranNamelistAssignments( block );

            for( std::size_t i = 0; i < assigns.size( ); ++i )
            {
                const std::string key = trimCopy( assigns.at( i ).first );
                const std::string value = trimCopy( assigns.at( i ).second );

                if( key == "PICNM" )
                {
                    imageContents.pictureName_ = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( value ) );
                    hasPictureName = true;
                }
                else if( key == "PICNO" )
                {
                    imageContents.pictureNumber_ = toInt( value );
                    hasPictureNumber = true;
                }
                else if( key == "TOB" )
                {
                    imageContents.endOfExposureTimeUtcString_ = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( value ) );
                    hasEndOfExposureTime = true;
                }
                else if( key == "CAMERA" )
                {
                    imageContents.cameraId_ = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( value ) );
                    hasCameraId = true;
                }
                else if( key == "EXPTIM" )
                {
                    imageContents.exposureTimeSeconds_ = toDouble( value );
                    hasExposureTime = true;
                }
                else if( key == "PICDEL" )
                {
                    imageContents.pictureDeletionFlag_ = toInt( value );
                    hasPictureDeletionFlag = true;
                }
                else if( key == "RA" )
                {
                    imageContents.rightAscensionDegrees_ = toDouble( value );
                    hasRightAscension = true;
                }
                else if( key == "DEC" )
                {
                    imageContents.declinationDegrees_ = toDouble( value );
                    hasDeclination = true;
                }
                else if( key == "TWIST" )
                {
                    imageContents.twistDegrees_ = toDouble( value );
                    hasTwist = true;
                }
            }

            if( !hasPictureName )
            {
                throwMissingPsfFieldError( "$PIC", "PICNM" );
            }

            // End-of-file sentinel picture
            if( imageContents.pictureName_ == "END" )
            {
                break;
            }

            const std::string pictureContext = "picture " + imageContents.pictureName_;
            if( !hasPictureNumber )
            {
                throwMissingPsfFieldError( "$PIC", "PICNO", pictureContext );
            }
            if( !hasEndOfExposureTime )
            {
                throwMissingPsfFieldError( "$PIC", "TOB", pictureContext );
            }
            if( !hasCameraId )
            {
                throwMissingPsfFieldError( "$PIC", "CAMERA", pictureContext );
            }
            if( !hasExposureTime )
            {
                throwMissingPsfFieldError( "$PIC", "EXPTIM", pictureContext );
            }
            if( !hasPictureDeletionFlag )
            {
                throwMissingPsfFieldError( "$PIC", "PICDEL", pictureContext );
            }
            if( !hasRightAscension )
            {
                throwMissingPsfFieldError( "$PIC", "RA", pictureContext );
            }
            if( !hasDeclination )
            {
                throwMissingPsfFieldError( "$PIC", "DEC", pictureContext );
            }
            if( !hasTwist )
            {
                throwMissingPsfFieldError( "$PIC", "TWIST", pictureContext );
            }
            if( fileContents.cameraProperties_.count( imageContents.cameraId_ ) == 0 )
            {
                throw std::runtime_error( "Error when reading PSF file: picture " + imageContents.pictureName_ +
                                          " references unknown camera " + imageContents.cameraId_ + "." );
            }

            // -----------------------------
            // Read $IM blocks until IMGTYP='END' or IMG='END'
            // -----------------------------
            while( std::getline( dataFile, line ) )
            {
                const std::string imLineTrimmed = trimCopy( line );
                if( imLineTrimmed.empty( ) )
                {
                    continue;
                }

                if( startsWithCaseSensitive( imLineTrimmed, "$PIC" ) )
                {
                    // Missing $IM END marker; carry $PIC back to outer loop
                    pendingLine = imLineTrimmed;
                    hasPendingLine = true;
                    break;
                }

                if( !startsWithCaseSensitive( imLineTrimmed, "$IM" ) )
                {
                    continue;
                }

                const std::string imBlock = readNamelistBlockContents( dataFile, imLineTrimmed );
                const std::vector< std::pair< std::string, std::string > > imAssigns = parseFortranNamelistAssignments( imBlock );

                // Convert assignments to map for easy lookup
                std::map< std::string, std::string > fields;
                for( std::size_t i = 0; i < imAssigns.size( ); ++i )
                {
                    fields[ trimCopy( imAssigns.at( i ).first ) ] = trimCopy( imAssigns.at( i ).second );
                }

                OpticalImageType opticalType = OpticalImageType::unknown;
                if( fields.count( "IMGTYP" ) )
                {
                    opticalType = parseOpticalImageTypeFromPsfToken( fields[ "IMGTYP" ] );
                }

                std::string imgNameForSentinel;
                if( fields.count( "IMG" ) )
                {
                    imgNameForSentinel = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( fields[ "IMG" ] ) );
                }

                if( opticalType == OpticalImageType::end_marker || imgNameForSentinel == "END" )
                {
                    break;  // end of this picture
                }

                std::string imageName;
                std::int64_t imageId = 0;
                int useFlag = 0;

                if( fields.count( "IMG" ) )
                {
                    imageName = stripSingleQuotesIfPresent( stripTrailingCommaAndTrim( fields[ "IMG" ] ) );
                }
                if( fields.count( "IMGID" ) )
                {
                    imageId = toInt64( fields[ "IMGID" ] );
                }
                if( fields.count( "USE" ) )
                {
                    useFlag = toInt( fields[ "USE" ] );
                }

                const bool hasZ = fields.count( "Z" ) > 0;
                const bool hasZc = fields.count( "ZC" ) > 0;
                const bool hasSig = fields.count( "SIG" ) > 0;

                if( !hasZ || !hasZc || !hasSig )
                {
                    std::ostringstream errorStream;
                    errorStream << "Error when reading PSF file: incomplete $IM entry.\n"
                                << "  Picture name     : '" << imageContents.pictureName_ << "'\n"
                                << "  Picture number   : " << imageContents.pictureNumber_ << "\n"
                                << "  Camera           : '" << imageContents.cameraId_ << "'\n"
                                << "  IMGTYP           : " << static_cast< int >( opticalType ) << "\n"
                                << "  IMG              : '" << imageName << "'\n"
                                << "  IMGID            : " << imageId << "\n"
                                << "  Missing fields   : " << ( hasZ ? "" : "Z " ) << ( hasZc ? "" : "ZC " ) << ( hasSig ? "" : "SIG " )
                                << "\n"
                                << "  Raw IM block     : " << imBlock << "\n";
                    throw std::runtime_error( errorStream.str( ) );
                }

                const Eigen::Vector2d z = parseVector2FromValueText( fields[ "Z" ], "Z" );
                const Eigen::Vector2d zc = parseVector2FromValueText( fields[ "ZC" ], "ZC" );
                const Eigen::Vector2d sig = parseVector2FromValueText( fields[ "SIG" ], "SIG" );

                std::shared_ptr< RawPsfMeasurement > measurement;

                if( opticalType == OpticalImageType::star )
                {
                    if( fields.count( "STRA" ) == 0 || fields.count( "STDEC" ) == 0 )
                    {
                        std::ostringstream errorStream;
                        errorStream << "Error when reading PSF file: STAR entry missing STRA/STDEC.\n"
                                    << "  Picture name     : '" << imageContents.pictureName_ << "'\n"
                                    << "  Picture number   : " << imageContents.pictureNumber_ << "\n"
                                    << "  IMG              : '" << imageName << "'\n"
                                    << "  Raw IM block     : " << imBlock << "\n";
                        throw std::runtime_error( errorStream.str( ) );
                    }

                    std::shared_ptr< RawPsfStarMeasurement > starMeasurement( new RawPsfStarMeasurement( ) );
                    starMeasurement->opticalImageType_ = opticalType;
                    starMeasurement->imageName_ = imageName;
                    starMeasurement->imageId_ = imageId;
                    starMeasurement->useFlag_ = useFlag;
                    starMeasurement->observedPixelLine_ = z;
                    starMeasurement->localCorrection_ = zc;
                    starMeasurement->sigmaPixelLine_ = sig;
                    starMeasurement->rightAscensionDegrees_ = toDouble( fields[ "STRA" ] );
                    starMeasurement->declinationDegrees_ = toDouble( fields[ "STDEC" ] );

                    measurement = starMeasurement;
                }
                else
                {
                    std::shared_ptr< RawPsfMeasurement > genericMeasurement( new RawPsfMeasurement( ) );
                    genericMeasurement->opticalImageType_ = opticalType;
                    genericMeasurement->imageName_ = imageName;
                    genericMeasurement->imageId_ = imageId;
                    genericMeasurement->useFlag_ = useFlag;
                    genericMeasurement->observedPixelLine_ = z;
                    genericMeasurement->localCorrection_ = zc;
                    genericMeasurement->sigmaPixelLine_ = sig;

                    measurement = genericMeasurement;
                }

                imageContents.measurements_.push_back( measurement );
            }

            fileContents.images_.push_back( imageContents );
            continue;
        }

        // Ignore other lines
    }

    if( !idParsed )
    {
        throw std::runtime_error( "Error when reading PSF file: no $ID block found." );
    }
    if( !camParsed )
    {
        throw std::runtime_error( "Error when reading PSF file: no $CAM block found." );
    }
    if( fileContents.numberOfCameras_ != static_cast< int >( fileContents.cameraProperties_.size( ) ) )
    {
        throw std::runtime_error( "Error when reading PSF file: NCAM is " + std::to_string( fileContents.numberOfCameras_ ) + " but " +
                                  std::to_string( fileContents.cameraProperties_.size( ) ) + " camera properties were parsed." );
    }

    return fileContents;
}

}  // namespace psf

}  // namespace input_output

}  // namespace tudat
