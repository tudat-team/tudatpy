/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readSumLmkFiles.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

namespace tudat
{

namespace input_output
{

namespace sum_lmk
{

namespace
{

constexpr double kilometersToMeters = 1000.0;
constexpr double parserTolerance = 1.0E-9;

enum class SumSection { none, landmarks, limbFits };
enum class LmkSection { none, pictures, mapOverlaps, limbFits };

std::string trimCopy( const std::string& input )
{
    return boost::algorithm::trim_copy( input );
}

std::string toUpperCopy( const std::string& input )
{
    std::string output = input;
    std::transform( output.begin( ), output.end( ), output.begin( ), []( const unsigned char character ) {
        return static_cast< char >( std::toupper( character ) );
    } );
    return output;
}

std::string stripTrailingComma( const std::string& input )
{
    std::string output = trimCopy( input );
    while( !output.empty( ) && output.back( ) == ',' )
    {
        output.pop_back( );
        output = trimCopy( output );
    }
    return output;
}

std::vector< std::string > tokenize( const std::string& input )
{
    std::stringstream stream( input );
    std::vector< std::string > tokens;
    std::string token;
    while( stream >> token )
    {
        tokens.push_back( stripTrailingComma( token ) );
    }
    return tokens;
}

std::string replaceFortranExponent( std::string token )
{
    std::replace( token.begin( ), token.end( ), 'D', 'E' );
    std::replace( token.begin( ), token.end( ), 'd', 'E' );
    return token;
}

std::string fileLineContext( const std::string& file, const int lineNumber, const std::string& id = "" )
{
    std::string context = file + ":" + std::to_string( lineNumber );
    if( !id.empty( ) )
    {
        context += " (" + id + ")";
    }
    return context;
}

double parseDoubleToken( const std::string& token, const std::string& context )
{
    const std::string normalizedToken = replaceFortranExponent( stripTrailingComma( token ) );
    try
    {
        std::size_t processedCharacters = 0;
        const double value = std::stod( normalizedToken, &processedCharacters );
        if( processedCharacters != normalizedToken.size( ) )
        {
            throw std::runtime_error( "" );
        }
        return value;
    }
    catch( const std::exception& )
    {
        throw std::runtime_error( "Error when parsing floating-point value '" + token + "' in " + context + "." );
    }
}

int parseIntToken( const std::string& token, const std::string& context )
{
    const std::string normalizedToken = stripTrailingComma( token );
    try
    {
        std::size_t processedCharacters = 0;
        const int value = std::stoi( normalizedToken, &processedCharacters );
        if( processedCharacters != normalizedToken.size( ) )
        {
            throw std::runtime_error( "" );
        }
        return value;
    }
    catch( const std::exception& )
    {
        throw std::runtime_error( "Error when parsing integer value '" + token + "' in " + context + "." );
    }
}

bool isStandaloneMarker( const std::string& trimmedLine, const std::string& marker )
{
    return toUpperCopy( trimmedLine ) == marker;
}

bool getLeadingTextForTrailingLabel( const std::string& trimmedLine, const std::string& label, std::string& leadingText )
{
    const std::string upperLine = toUpperCopy( trimmedLine );
    const std::string upperLabel = toUpperCopy( label );
    if( upperLine.size( ) < upperLabel.size( ) ||
        upperLine.compare( upperLine.size( ) - upperLabel.size( ), upperLabel.size( ), upperLabel ) != 0 )
    {
        return false;
    }

    if( upperLine.size( ) > upperLabel.size( ) )
    {
        const char precedingCharacter = upperLine.at( upperLine.size( ) - upperLabel.size( ) - 1 );
        if( !std::isspace( static_cast< unsigned char >( precedingCharacter ) ) )
        {
            return false;
        }
    }

    leadingText = trimCopy( trimmedLine.substr( 0, trimmedLine.size( ) - label.size( ) ) );
    return true;
}

std::vector< std::string > getLeadingTokensForLabel( const std::string& trimmedLine,
                                                     const std::string& label,
                                                     const std::string& file,
                                                     const int lineNumber,
                                                     const std::string& id,
                                                     const int expectedNumberOfTokens )
{
    std::string leadingText;
    if( !getLeadingTextForTrailingLabel( trimmedLine, label, leadingText ) )
    {
        return std::vector< std::string >( );
    }

    const std::vector< std::string > tokens = tokenize( leadingText );
    if( static_cast< int >( tokens.size( ) ) != expectedNumberOfTokens )
    {
        throw std::runtime_error( "Error when parsing label '" + label + "' in " + fileLineContext( file, lineNumber, id ) + ": expected " +
                                  std::to_string( expectedNumberOfTokens ) + " leading fields, found " + std::to_string( tokens.size( ) ) +
                                  "." );
    }
    return tokens;
}

std::vector< std::string > getLeadingTokensForLabelRange( const std::string& trimmedLine,
                                                          const std::string& label,
                                                          const std::string& file,
                                                          const int lineNumber,
                                                          const std::string& id,
                                                          const int minimumNumberOfTokens,
                                                          const int maximumNumberOfTokens )
{
    std::string leadingText;
    if( !getLeadingTextForTrailingLabel( trimmedLine, label, leadingText ) )
    {
        return std::vector< std::string >( );
    }

    const std::vector< std::string > tokens = tokenize( leadingText );
    if( static_cast< int >( tokens.size( ) ) < minimumNumberOfTokens || static_cast< int >( tokens.size( ) ) > maximumNumberOfTokens )
    {
        throw std::runtime_error( "Error when parsing label '" + label + "' in " + fileLineContext( file, lineNumber, id ) +
                                  ": expected between " + std::to_string( minimumNumberOfTokens ) + " and " +
                                  std::to_string( maximumNumberOfTokens ) + " leading fields, found " + std::to_string( tokens.size( ) ) +
                                  "." );
    }
    return tokens;
}

Eigen::Vector3d parseVector3( const std::vector< std::string >& tokens, const std::string& context, const double scale = 1.0 )
{
    return ( Eigen::Vector3d( ) << scale * parseDoubleToken( tokens.at( 0 ), context ),
             scale * parseDoubleToken( tokens.at( 1 ), context ),
             scale * parseDoubleToken( tokens.at( 2 ), context ) )
            .finished( );
}

bool isFinite( const Eigen::Vector3d& vector )
{
    return vector.array( ).isFinite( ).all( );
}

bool isFinite( const Eigen::Matrix3d& matrix )
{
    return matrix.array( ).isFinite( ).all( );
}

bool isFinite( const Eigen::Matrix< double, 2, 3 >& matrix )
{
    return matrix.array( ).isFinite( ).all( );
}

bool isClose( const double lhs, const double rhs, const double tolerance = parserTolerance )
{
    return std::fabs( lhs - rhs ) <= tolerance * std::max( 1.0, std::max( std::fabs( lhs ), std::fabs( rhs ) ) );
}

template< typename Derived1, typename Derived2 >
bool isClose( const Eigen::MatrixBase< Derived1 >& lhs, const Eigen::MatrixBase< Derived2 >& rhs, const double tolerance = parserTolerance )
{
    return ( lhs - rhs ).norm( ) <= tolerance * std::max( 1.0, std::max( lhs.norm( ), rhs.norm( ) ) );
}

void require( const bool condition, const std::string& message )
{
    if( !condition )
    {
        throw std::runtime_error( message );
    }
}

void validateSumImage( const SumImageData& image,
                       const std::string& file,
                       const bool hasImageSize,
                       const bool hasCalibration,
                       const bool hasCameraAxes,
                       const bool hasKMatrix )
{
    require( !image.imageId_.empty( ), "Error when reading SUM file " + file + ": missing image ID." );
    require( !image.utcEpochString_.empty( ), "Error when reading SUM file " + file + " (" + image.imageId_ + "): missing UTC epoch." );
    require( hasImageSize, "Error when reading SUM file " + file + " (" + image.imageId_ + "): missing NPX/NLN/THRSH record." );
    require( hasCalibration, "Error when reading SUM file " + file + " (" + image.imageId_ + "): missing MMFL/CTR record." );
    require( hasCameraAxes && isFinite( image.cameraAxes_ ),
             "Error when reading SUM file " + file + " (" + image.imageId_ + "): missing one or more CX/CY/CZ records." );
    require( hasKMatrix && isFinite( image.kMatrix_ ),
             "Error when reading SUM file " + file + " (" + image.imageId_ + "): missing K-MATRIX record." );
}

void validateLmkLandmark( const LmkLandmarkData& landmark, const std::string& file )
{
    require( !landmark.landmarkId_.empty( ), "Error when reading LMK file " + file + ": missing landmark ID." );
    require( isFinite( landmark.bodyFixedPosition_ ),
             "Error when reading LMK file " + file + " (" + landmark.landmarkId_ + "): missing VLM record." );
}

bool haveSameImageMetadata( const SumImageData& lhs, const SumImageData& rhs )
{
    return lhs.utcEpochString_ == rhs.utcEpochString_ && lhs.imageSize_ == rhs.imageSize_ && lhs.threshold_ == rhs.threshold_ &&
            lhs.maxDn_ == rhs.maxDn_ && isClose( lhs.focalLengthMm_, rhs.focalLengthMm_ ) &&
            isClose( lhs.opticalCenter_, rhs.opticalCenter_ ) && isClose( lhs.cameraAxes_, rhs.cameraAxes_ ) &&
            isClose( lhs.kMatrix_, rhs.kMatrix_ );
}

bool haveSameLandmarkDefinition( const LmkLandmarkData& lhs, const LmkLandmarkData& rhs )
{
    return lhs.landmarkId_ == rhs.landmarkId_ && lhs.typeFlag_ == rhs.typeFlag_ &&
            isClose( lhs.bodyFixedPosition_, rhs.bodyFixedPosition_ ) && isClose( lhs.localXAxis_, rhs.localXAxis_ ) &&
            isClose( lhs.localYAxis_, rhs.localYAxis_ ) && isClose( lhs.localZAxis_, rhs.localZAxis_ );
}

void addSumLandmarkObservation( SumImageData& image,
                                std::map< std::string, Eigen::Vector2d >& landmarkRows,
                                const std::string& file,
                                const int lineNumber,
                                const std::vector< std::string >& tokens )
{
    if( tokens.size( ) != 3 && !( tokens.size( ) == 4 && tokens.at( 3 ) == "-" ) )
    {
        throw std::runtime_error( "Error when parsing SUM LANDMARKS row in " + fileLineContext( file, lineNumber, image.imageId_ ) +
                                  ": expected 3 fields plus optional '-' flag." );
    }

    SumLandmarkObservation observation;
    observation.landmarkId_ = tokens.at( 0 );
    const std::string context = fileLineContext( file, lineNumber, image.imageId_ ) + " LANDMARKS";
    observation.pixelCoordinates_ =
            ( Eigen::Vector2d( ) << parseDoubleToken( tokens.at( 1 ), context ), parseDoubleToken( tokens.at( 2 ), context ) ).finished( );

    if( landmarkRows.count( observation.landmarkId_ ) != 0 )
    {
        throw std::runtime_error( "Error when parsing SUM file " + file + " (" + image.imageId_ +
                                  "): duplicate LANDMARKS row for landmark '" + observation.landmarkId_ + "'." );
    }

    landmarkRows[ observation.landmarkId_ ] = observation.pixelCoordinates_;
    image.landmarkObservations_.push_back( observation );
}

void addSumLimbFitObservation( SumImageData& image,
                               const std::string& file,
                               const int lineNumber,
                               const std::vector< std::string >& tokens )
{
    if( tokens.size( ) != 4 )
    {
        throw std::runtime_error( "Error when parsing SUM LIMB FITS row in " + fileLineContext( file, lineNumber, image.imageId_ ) +
                                  ": expected 4 fields." );
    }

    SumLimbFitObservation observation;
    observation.featureId_ = tokens.at( 0 );
    const std::string context = fileLineContext( file, lineNumber, image.imageId_ ) + " LIMB FITS";
    observation.pixelCoordinates_ =
            ( Eigen::Vector2d( ) << parseDoubleToken( tokens.at( 1 ), context ), parseDoubleToken( tokens.at( 2 ), context ) ).finished( );
    observation.sigma_ = parseDoubleToken( tokens.at( 3 ), context );
    image.limbFitObservations_.push_back( observation );
}

}  // namespace

SumImageData readSumFile( const std::string& sumFile )
{
    std::ifstream stream( sumFile );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when reading SUM file " + sumFile + ": file could not be opened." );
    }

    SumImageData image;
    image.sourceFile_ = sumFile;
    SumSection activeSection = SumSection::none;
    std::map< std::string, Eigen::Vector2d > landmarkRows;
    bool hasImageSize = false;
    bool hasCalibration = false;
    bool hasCx = false;
    bool hasCy = false;
    bool hasCz = false;
    bool hasKMatrix = false;

    std::string line;
    int lineNumber = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        const std::string trimmedLine = trimCopy( line );
        if( trimmedLine.empty( ) )
        {
            continue;
        }

        if( image.imageId_.empty( ) )
        {
            image.imageId_ = trimmedLine;
            continue;
        }
        if( image.utcEpochString_.empty( ) )
        {
            image.utcEpochString_ = trimmedLine;
            continue;
        }

        if( isStandaloneMarker( trimmedLine, "END FILE" ) )
        {
            break;
        }
        if( isStandaloneMarker( trimmedLine, "LANDMARKS" ) )
        {
            activeSection = SumSection::landmarks;
            continue;
        }
        if( isStandaloneMarker( trimmedLine, "LIMB FITS" ) )
        {
            activeSection = SumSection::limbFits;
            continue;
        }

        if( activeSection == SumSection::landmarks )
        {
            addSumLandmarkObservation( image, landmarkRows, sumFile, lineNumber, tokenize( trimmedLine ) );
            continue;
        }
        if( activeSection == SumSection::limbFits )
        {
            addSumLimbFitObservation( image, sumFile, lineNumber, tokenize( trimmedLine ) );
            continue;
        }

        std::vector< std::string > tokens;
        const std::string context = fileLineContext( sumFile, lineNumber, image.imageId_ );
        if( !( tokens = getLeadingTokensForLabel( trimmedLine, "NPX, NLN, THRSH", sumFile, lineNumber, image.imageId_, 4 ) ).empty( ) )
        {
            image.imageSize_ = Eigen::Vector2i( parseIntToken( tokens.at( 0 ), context ), parseIntToken( tokens.at( 1 ), context ) );
            image.threshold_ = parseIntToken( tokens.at( 2 ), context );
            image.maxDn_ = parseIntToken( tokens.at( 3 ), context );
            hasImageSize = true;
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "MMFL, CTR", sumFile, lineNumber, image.imageId_, 3 ) ).empty( ) )
        {
            image.focalLengthMm_ = parseDoubleToken( tokens.at( 0 ), context );
            image.opticalCenter_ =
                    ( Eigen::Vector2d( ) << parseDoubleToken( tokens.at( 1 ), context ), parseDoubleToken( tokens.at( 2 ), context ) )
                            .finished( );
            hasCalibration = true;
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "SCOBJ", sumFile, lineNumber, image.imageId_, 3 ) ).empty( ) )
        {
            image.spacecraftObjectVector_ = parseVector3( tokens, context, kilometersToMeters );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "CX", sumFile, lineNumber, image.imageId_, 3 ) ).empty( ) )
        {
            image.cameraAxes_.row( 0 ) = parseVector3( tokens, context ).transpose( );
            hasCx = true;
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "CY", sumFile, lineNumber, image.imageId_, 3 ) ).empty( ) )
        {
            image.cameraAxes_.row( 1 ) = parseVector3( tokens, context ).transpose( );
            hasCy = true;
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "CZ", sumFile, lineNumber, image.imageId_, 3 ) ).empty( ) )
        {
            image.cameraAxes_.row( 2 ) = parseVector3( tokens, context ).transpose( );
            hasCz = true;
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "SZ", sumFile, lineNumber, image.imageId_, 3 ) ).empty( ) )
        {
            image.sunDirectionBodyFixed_ = parseVector3( tokens, context );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "K-MATRIX", sumFile, lineNumber, image.imageId_, 6 ) ).empty( ) )
        {
            image.kMatrix_( 0, 0 ) = parseDoubleToken( tokens.at( 0 ), context );
            image.kMatrix_( 0, 1 ) = parseDoubleToken( tokens.at( 1 ), context );
            image.kMatrix_( 0, 2 ) = parseDoubleToken( tokens.at( 2 ), context );
            image.kMatrix_( 1, 0 ) = parseDoubleToken( tokens.at( 3 ), context );
            image.kMatrix_( 1, 1 ) = parseDoubleToken( tokens.at( 4 ), context );
            image.kMatrix_( 1, 2 ) = parseDoubleToken( tokens.at( 5 ), context );
            hasKMatrix = true;
        }
        else if( !( tokens = getLeadingTokensForLabelRange( trimmedLine, "SIGMA_PTG", sumFile, lineNumber, image.imageId_, 1, 3 ) )
                          .empty( ) )
        {
            if( tokens.size( ) == 1 )
            {
                image.pointingSigma_ = Eigen::Vector3d::Constant( parseDoubleToken( tokens.at( 0 ), context ) );
            }
            else
            {
                image.pointingSigma_ = parseVector3( tokens, context );
            }
        }
        else if( !( tokens = getLeadingTokensForLabelRange( trimmedLine, "SIGMA_VSO", sumFile, lineNumber, image.imageId_, 1, 3 ) )
                          .empty( ) )
        {
            if( tokens.size( ) == 1 )
            {
                image.spacecraftObjectSigma_ =
                        Eigen::Vector3d::Constant( kilometersToMeters * parseDoubleToken( tokens.at( 0 ), context ) );
            }
            else
            {
                image.spacecraftObjectSigma_ = parseVector3( tokens, context, kilometersToMeters );
            }
        }
        else if( !( tokens = getLeadingTokensForLabelRange( trimmedLine, "DISTORTION", sumFile, lineNumber, image.imageId_, 4, 4 ) )
                          .empty( ) ||
                 !( tokens = getLeadingTokensForLabelRange(
                            trimmedLine, "DISTORTION COEFFICIENTS", sumFile, lineNumber, image.imageId_, 4, 4 ) )
                          .empty( ) ||
                 !( tokens = getLeadingTokensForLabelRange( trimmedLine, "DCOEF", sumFile, lineNumber, image.imageId_, 4, 4 ) ).empty( ) )
        {
            for( int i = 0; i < 4; ++i )
            {
                image.distortionCoefficients_( i ) = parseDoubleToken( tokens.at( i ), context );
            }
        }
    }

    validateSumImage( image, sumFile, hasImageSize, hasCalibration, hasCx && hasCy && hasCz, hasKMatrix );
    return image;
}

std::vector< SumImageData > readSumFiles( const std::vector< std::string >& sumFiles )
{
    std::vector< SumImageData > images;
    std::map< std::string, std::size_t > imageIndexById;

    for( const std::string& sumFile : sumFiles )
    {
        const SumImageData image = readSumFile( sumFile );
        if( image.landmarkObservations_.empty( ) )
        {
            continue;
        }
        if( imageIndexById.count( image.imageId_ ) == 0 )
        {
            imageIndexById[ image.imageId_ ] = images.size( );
            images.push_back( image );
            continue;
        }

        SumImageData& existingImage = images.at( imageIndexById.at( image.imageId_ ) );
        if( !haveSameImageMetadata( existingImage, image ) )
        {
            throw std::runtime_error( "Error when merging SUM files: image '" + image.imageId_ + "' has conflicting metadata in '" +
                                      existingImage.sourceFile_ + "' and '" + image.sourceFile_ + "'." );
        }

        for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
        {
            auto existingObservationIterator = std::find_if( existingImage.landmarkObservations_.begin( ),
                                                             existingImage.landmarkObservations_.end( ),
                                                             [ & ]( const SumLandmarkObservation& existingObservation ) {
                                                                 return existingObservation.landmarkId_ == observation.landmarkId_;
                                                             } );
            if( existingObservationIterator == existingImage.landmarkObservations_.end( ) )
            {
                existingImage.landmarkObservations_.push_back( observation );
            }
            else if( !isClose( existingObservationIterator->pixelCoordinates_, observation.pixelCoordinates_ ) )
            {
                throw std::runtime_error( "Error when merging SUM files: image/landmark pair ('" + image.imageId_ + "', '" +
                                          observation.landmarkId_ + "') has conflicting pixel coordinates in '" +
                                          existingImage.sourceFile_ + "' and '" + image.sourceFile_ + "'." );
            }
        }

        existingImage.limbFitObservations_.insert(
                existingImage.limbFitObservations_.end( ), image.limbFitObservations_.begin( ), image.limbFitObservations_.end( ) );
    }

    return images;
}

LmkLandmarkData readLmkFile( const std::string& lmkFile )
{
    std::ifstream stream( lmkFile );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when reading LMK file " + lmkFile + ": file could not be opened." );
    }

    LmkLandmarkData landmark;
    landmark.sourceFile_ = lmkFile;
    LmkSection activeSection = LmkSection::none;

    std::string line;
    int lineNumber = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        const std::string trimmedLine = trimCopy( line );
        if( trimmedLine.empty( ) )
        {
            continue;
        }

        if( landmark.landmarkId_.empty( ) )
        {
            const std::vector< std::string > headerTokens = tokenize( trimmedLine );
            if( headerTokens.empty( ) )
            {
                continue;
            }
            landmark.landmarkId_ = headerTokens.at( 0 );
            if( headerTokens.size( ) > 1 && !headerTokens.at( 1 ).empty( ) )
            {
                landmark.typeFlag_ = headerTokens.at( 1 ).at( 0 );
            }
            continue;
        }

        if( isStandaloneMarker( trimmedLine, "END FILE" ) )
        {
            break;
        }
        if( isStandaloneMarker( trimmedLine, "PICTURES" ) )
        {
            activeSection = LmkSection::pictures;
            continue;
        }
        if( isStandaloneMarker( trimmedLine, "MAP OVERLAPS" ) )
        {
            activeSection = LmkSection::mapOverlaps;
            continue;
        }
        if( isStandaloneMarker( trimmedLine, "LIMB FITS" ) )
        {
            activeSection = LmkSection::limbFits;
            continue;
        }

        const std::vector< std::string > sectionTokens = tokenize( trimmedLine );
        if( activeSection == LmkSection::pictures )
        {
            if( sectionTokens.size( ) != 3 && sectionTokens.size( ) != 4 )
            {
                throw std::runtime_error( "Error when parsing LMK PICTURES row in " +
                                          fileLineContext( lmkFile, lineNumber, landmark.landmarkId_ ) +
                                          ": expected 3 fields plus optional '*' flag." );
            }
            LmkPictureObservation picture;
            picture.imageId_ = sectionTokens.at( 0 );
            const std::string context = fileLineContext( lmkFile, lineNumber, landmark.landmarkId_ ) + " PICTURES";
            picture.pixelCoordinates_ = ( Eigen::Vector2d( ) << parseDoubleToken( sectionTokens.at( 1 ), context ),
                                          parseDoubleToken( sectionTokens.at( 2 ), context ) )
                                                .finished( );
            picture.flagged_ = ( sectionTokens.size( ) == 4 && sectionTokens.at( 3 ) == "*" );
            landmark.pictures_.push_back( picture );
            continue;
        }
        if( activeSection == LmkSection::mapOverlaps )
        {
            if( sectionTokens.size( ) != 4 )
            {
                throw std::runtime_error( "Error when parsing LMK MAP OVERLAPS row in " +
                                          fileLineContext( lmkFile, lineNumber, landmark.landmarkId_ ) + ": expected 4 fields." );
            }
            continue;
        }
        if( activeSection == LmkSection::limbFits )
        {
            continue;
        }

        std::vector< std::string > tokens;
        const std::string context = fileLineContext( lmkFile, lineNumber, landmark.landmarkId_ );
        // Production SPC LMK files combine parsed-only fields under multi-token trailing labels
        // ("SIZE, SCALE(KM)" and "SIGKM, RMSLMK"). These must be matched before the single-token
        // labels below, because e.g. a "... SIGKM, RMSLMK" line also ends with "RMSLMK".
        if( !( tokens = getLeadingTokensForLabel( trimmedLine, "SIZE, SCALE(KM)", lmkFile, lineNumber, landmark.landmarkId_, 2 ) )
                     .empty( ) )
        {
            landmark.patchSize_ = parseIntToken( tokens.at( 0 ), context );
            landmark.patchScaleKm_ = parseDoubleToken( tokens.at( 1 ), context );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "SIGKM, RMSLMK", lmkFile, lineNumber, landmark.landmarkId_, 2 ) )
                          .empty( ) )
        {
            landmark.sigKm_ = parseDoubleToken( tokens.at( 0 ), context );
            landmark.rmsLmk_ = parseDoubleToken( tokens.at( 1 ), context );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "SIZE", lmkFile, lineNumber, landmark.landmarkId_, 1 ) ).empty( ) )
        {
            landmark.patchSize_ = parseIntToken( tokens.at( 0 ), context );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "SCALE", lmkFile, lineNumber, landmark.landmarkId_, 1 ) ).empty( ) )
        {
            landmark.patchScaleKm_ = parseDoubleToken( tokens.at( 0 ), context );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "HORIZON", lmkFile, lineNumber, landmark.landmarkId_, 4 ) ).empty( ) )
        {
            for( int i = 0; i < 4; ++i )
            {
                landmark.horizonFlags_( i ) = parseIntToken( tokens.at( i ), context );
            }
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "SIGKM", lmkFile, lineNumber, landmark.landmarkId_, 1 ) ).empty( ) )
        {
            landmark.sigKm_ = parseDoubleToken( tokens.at( 0 ), context );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "RMSLMK", lmkFile, lineNumber, landmark.landmarkId_, 1 ) ).empty( ) )
        {
            landmark.rmsLmk_ = parseDoubleToken( tokens.at( 0 ), context );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "VLM", lmkFile, lineNumber, landmark.landmarkId_, 3 ) ).empty( ) )
        {
            landmark.bodyFixedPosition_ = parseVector3( tokens, context, kilometersToMeters );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "UX", lmkFile, lineNumber, landmark.landmarkId_, 3 ) ).empty( ) )
        {
            landmark.localXAxis_ = parseVector3( tokens, context );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "UY", lmkFile, lineNumber, landmark.landmarkId_, 3 ) ).empty( ) )
        {
            landmark.localYAxis_ = parseVector3( tokens, context );
        }
        else if( !( tokens = getLeadingTokensForLabel( trimmedLine, "UZ", lmkFile, lineNumber, landmark.landmarkId_, 3 ) ).empty( ) )
        {
            landmark.localZAxis_ = parseVector3( tokens, context );
        }
        else if( !( tokens = getLeadingTokensForLabelRange( trimmedLine, "SIGMA_LMK", lmkFile, lineNumber, landmark.landmarkId_, 1, 3 ) )
                          .empty( ) )
        {
            if( tokens.size( ) == 1 )
            {
                landmark.landmarkPositionSigma_ =
                        Eigen::Vector3d::Constant( kilometersToMeters * parseDoubleToken( tokens.at( 0 ), context ) );
            }
            else
            {
                landmark.landmarkPositionSigma_ = parseVector3( tokens, context, kilometersToMeters );
            }
        }
    }

    validateLmkLandmark( landmark, lmkFile );
    return landmark;
}

std::map< std::string, LmkLandmarkData > readLmkFiles( const std::vector< std::string >& lmkFiles )
{
    std::map< std::string, LmkLandmarkData > landmarks;
    for( const std::string& lmkFile : lmkFiles )
    {
        const LmkLandmarkData landmark = readLmkFile( lmkFile );
        if( landmarks.count( landmark.landmarkId_ ) == 0 )
        {
            landmarks[ landmark.landmarkId_ ] = landmark;
        }
        else if( !haveSameLandmarkDefinition( landmarks.at( landmark.landmarkId_ ), landmark ) )
        {
            throw std::runtime_error( "Error when merging LMK files: landmark '" + landmark.landmarkId_ +
                                      "' has conflicting definitions in '" + landmarks.at( landmark.landmarkId_ ).sourceFile_ + "' and '" +
                                      landmark.sourceFile_ + "'." );
        }
    }
    return landmarks;
}

}  // namespace sum_lmk

}  // namespace input_output

}  // namespace tudat
