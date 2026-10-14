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

#include <Eigen/LU>

#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/io/readTrackingTxtFile.h"

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

std::vector< std::string > getLeadingTokensForLabel( const std::string& trimmedLine,
                                                     const std::string& label,
                                                     const std::string& file,
                                                     const int lineNumber,
                                                     const std::string& id,
                                                     const int expectedNumberOfTokens )
{
    return getLeadingTokensForLabelRange( trimmedLine, label, file, lineNumber, id, expectedNumberOfTokens, expectedNumberOfTokens );
}

Eigen::Vector3d parseVector3( const std::vector< std::string >& tokens, const std::string& context, const double scale = 1.0 )
{
    return ( Eigen::Vector3d( ) << scale * parseDoubleToken( tokens.at( 0 ), context ),
             scale * parseDoubleToken( tokens.at( 1 ), context ),
             scale * parseDoubleToken( tokens.at( 2 ), context ) )
            .finished( );
}

template< typename Derived >
bool isFinite( const Eigen::MatrixBase< Derived >& matrix )
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

template< typename Derived1, typename Derived2 >
bool isSameOptionalVector( const Eigen::MatrixBase< Derived1 >& lhs, const Eigen::MatrixBase< Derived2 >& rhs )
{
    const bool lhsIsFinite = lhs.array( ).isFinite( ).all( );
    const bool rhsIsFinite = rhs.array( ).isFinite( ).all( );
    if( lhsIsFinite || rhsIsFinite )
    {
        return lhsIsFinite && rhsIsFinite && isClose( lhs, rhs );
    }
    return lhs.array( ).isNaN( ).all( ) && rhs.array( ).isNaN( ).all( );
}

void require( const bool condition, const std::string& message )
{
    if( !condition )
    {
        throw std::runtime_error( message );
    }
}

void validateSumLandmarkObservation( const SumLandmarkObservation& observation, const std::string& errorPrefix )
{
    require( observation.pixelCoordinates_.array( ).isFinite( ).all( ),
             errorPrefix + "LANDMARKS pixel coordinates for '" + observation.landmarkId_ + "' must be finite." );
}

//! Shared value validation for both file-loaded and programmatically constructed SUM images.
void validateSumImageValues( const SumImageData& image, const std::string& errorPrefix )
{
    require( image.imageSize_( 0 ) > 0 && image.imageSize_( 1 ) > 0 && image.threshold_ >= 0 && image.maxDn_ >= image.threshold_,
             errorPrefix + "NPX and NLN must be positive and 0 <= THRSH <= MAXDN." );
    require( std::isfinite( image.focalLengthMm_ ) && image.focalLengthMm_ > std::numeric_limits< double >::epsilon( ),
             errorPrefix + "MMFL must be finite and positive." );
    require( isFinite( image.opticalCenter_ ), errorPrefix + "CTR coordinates must be finite." );
    require( isFinite( image.cameraAxes_ ), errorPrefix + "CX/CY/CZ entries must be finite." );
    const bool hasProperCameraAttitude = ( image.cameraAxes_ * image.cameraAxes_.transpose( ) - Eigen::Matrix3d::Identity( ) ).norm( ) <=
                    sumCameraRotationMatrixTolerance &&
            std::fabs( image.cameraAxes_.determinant( ) - 1.0 ) <= sumCameraRotationMatrixTolerance;
    require( hasProperCameraAttitude, errorPrefix + "CX/CY/CZ must form an orthonormal right-handed rotation matrix." );
    require( isFinite( image.kMatrix_ ) && image.kMatrix_.block< 2, 2 >( 0, 0 ).fullPivLu( ).isInvertible( ),
             errorPrefix + "leading 2x2 K-MATRIX block must be finite and nonsingular." );
    if( !image.pointingSigma_.array( ).isNaN( ).all( ) )
    {
        require( image.pointingSigma_.array( ).isFinite( ).all( ) && ( image.pointingSigma_.array( ) > 0.0 ).all( ),
                 errorPrefix + "SIGMA_PTG must contain three finite positive values, or all NaNs when absent." );
    }
    for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
    {
        validateSumLandmarkObservation( observation, errorPrefix );
    }
}

void validateSumImage( const SumImageData& image,
                       const std::string& file,
                       const bool hasImageSize,
                       const bool hasCalibration,
                       const bool hasCameraAxes,
                       const bool hasKMatrix )
{
    const std::string errorPrefix = "Error when reading SUM file " + file + " (" + image.imageId_ + "): ";
    require( !image.imageId_.empty( ), errorPrefix + "missing image ID." );
    require( !image.utcEpochString_.empty( ), errorPrefix + "missing UTC epoch." );
    require( hasImageSize, errorPrefix + "missing NPX/NLN/THRSH record." );
    require( hasCalibration, errorPrefix + "missing MMFL/CTR record." );
    require( hasCameraAxes, errorPrefix + "missing one or more CX/CY/CZ records." );
    require( hasKMatrix, errorPrefix + "missing K-MATRIX record." );
    validateSumImageValues( image, errorPrefix );
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
            isClose( lhs.kMatrix_, rhs.kMatrix_ ) && isSameOptionalVector( lhs.spacecraftObjectVector_, rhs.spacecraftObjectVector_ ) &&
            isSameOptionalVector( lhs.sunDirectionBodyFixed_, rhs.sunDirectionBodyFixed_ ) &&
            isClose( lhs.distortionCoefficients_, rhs.distortionCoefficients_ ) &&
            isSameOptionalVector( lhs.spacecraftObjectSigma_, rhs.spacecraftObjectSigma_ ) &&
            isSameOptionalVector( lhs.pointingSigma_, rhs.pointingSigma_ );
}

bool haveSameLandmarkDefinition( const LmkLandmarkData& lhs, const LmkLandmarkData& rhs )
{
    return lhs.landmarkId_ == rhs.landmarkId_ && lhs.typeFlag_ == rhs.typeFlag_ &&
            isClose( lhs.bodyFixedPosition_, rhs.bodyFixedPosition_ ) && isSameOptionalVector( lhs.localXAxis_, rhs.localXAxis_ ) &&
            isSameOptionalVector( lhs.localYAxis_, rhs.localYAxis_ ) && isSameOptionalVector( lhs.localZAxis_, rhs.localZAxis_ );
}

//! Merge repeated images and retain each image/landmark measurement once, preserving input order.
std::vector< SumImageData > mergeSumImages( const std::vector< SumImageData >& inputImages )
{
    std::vector< SumImageData > images;
    std::map< std::string, std::size_t > imageIndexById;
    for( const SumImageData& image : inputImages )
    {
        if( image.landmarkObservations_.empty( ) )
        {
            continue;
        }
        const auto inserted = imageIndexById.emplace( image.imageId_, images.size( ) );
        if( inserted.second )
        {
            images.push_back( image );
            images.back( ).landmarkObservations_.clear( );
            images.back( ).limbFitObservations_.clear( );
        }
        SumImageData& existingImage = images.at( inserted.first->second );
        if( !inserted.second && !haveSameImageMetadata( existingImage, image ) )
        {
            throw std::runtime_error( "Error when merging SUM data: image '" + image.imageId_ + "' has conflicting metadata in '" +
                                      existingImage.sourceFile_ + "' and '" + image.sourceFile_ + "'." );
        }
        for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
        {
            const auto existing =
                    std::find_if( existingImage.landmarkObservations_.begin( ),
                                  existingImage.landmarkObservations_.end( ),
                                  [ & ]( const SumLandmarkObservation& entry ) { return entry.landmarkId_ == observation.landmarkId_; } );
            if( existing == existingImage.landmarkObservations_.end( ) )
            {
                existingImage.landmarkObservations_.push_back( observation );
            }
            else if( !isClose( existing->pixelCoordinates_, observation.pixelCoordinates_ ) )
            {
                throw std::runtime_error( "Error when merging SUM data: image/landmark pair ('" + image.imageId_ + "', '" +
                                          observation.landmarkId_ + "') has conflicting pixel coordinates in '" +
                                          existingImage.sourceFile_ + "' and '" + image.sourceFile_ + "'." );
            }
        }
        for( const SumLimbFitObservation& observation : image.limbFitObservations_ )
        {
            const auto existing = std::find_if( existingImage.limbFitObservations_.begin( ),
                                                existingImage.limbFitObservations_.end( ),
                                                [ & ]( const SumLimbFitObservation& entry ) {
                                                    return entry.featureId_ == observation.featureId_ &&
                                                            isClose( entry.pixelCoordinates_, observation.pixelCoordinates_ ) &&
                                                            isClose( entry.sigma_, observation.sigma_ );
                                                } );
            if( existing == existingImage.limbFitObservations_.end( ) )
            {
                existingImage.limbFitObservations_.push_back( observation );
            }
        }
    }
    return images;
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
    validateSumLandmarkObservation( observation, "Error when reading " + context + ": " );

    const auto inserted = landmarkRows.emplace( observation.landmarkId_, observation.pixelCoordinates_ );
    if( !inserted.second )
    {
        if( !isClose( inserted.first->second, observation.pixelCoordinates_ ) )
        {
            throw std::runtime_error( "Error when parsing SUM file " + file + " (" + image.imageId_ +
                                      "): conflicting LANDMARKS rows for landmark '" + observation.landmarkId_ + "'." );
        }
        return;
    }

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
    for( const std::string& sumFile : sumFiles )
    {
        images.push_back( readSumFile( sumFile ) );
    }
    return mergeSumImages( images );
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

namespace
{

//! Check the SUM camera calibration and attitude, also for SumImageData that was built in code.
void validateSumImageCalibrationAndAttitude( const SumImageData& image )
{
    const std::string errorPrefix = "Error when converting SUM image '" + image.imageId_ + "': ";
    validateSumImageValues( image, errorPrefix );

    // The SUM DISTORTION row is not a full distortion model. SPC keeps it zero and gives a non-zero Owen model
    // separately in INIT_LITHOS. Rejecting non-zero values avoids silently simulating a different camera.
    if( !image.distortionCoefficients_.array( ).isFinite( ).all( ) )
    {
        throw std::runtime_error( errorPrefix + "DISTORTION coefficients must be finite." );
    }
    if( !image.distortionCoefficients_.isZero( 0.0 ) )
    {
        throw std::runtime_error( errorPrefix +
                                  "non-zero SUM DISTORTION coefficients are unsupported. SPC supplies a non-zero Owen distortion "
                                  "model separately through INIT_LITHOS. It cannot be reconstructed from this SUM row." );
    }
}

void validateReferencedLandmarksHaveDefinitions( const std::vector< SumImageData >& sumImages,
                                                 const std::map< std::string, LmkLandmarkData >& landmarks )
{
    std::map< std::string, std::set< std::string > > missingLandmarkImages;
    for( const SumImageData& image : sumImages )
    {
        for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
        {
            if( landmarks.count( observation.landmarkId_ ) == 0 )
            {
                missingLandmarkImages[ observation.landmarkId_ ].insert( image.imageId_ );
            }
        }
    }

    if( missingLandmarkImages.empty( ) )
    {
        return;
    }

    std::ostringstream errorMessage;
    errorMessage << "Error when converting SUM/LMK observations: missing LMK data for " << missingLandmarkImages.size( )
                 << " landmark(s): ";
    bool firstLandmark = true;
    for( const auto& missingEntry : missingLandmarkImages )
    {
        errorMessage << ( firstLandmark ? "" : "; " ) << missingEntry.first << " referenced by image(s) ";
        firstLandmark = false;

        bool firstImage = true;
        for( const std::string& imageId : missingEntry.second )
        {
            errorMessage << ( firstImage ? "" : ", " ) << imageId;
            firstImage = false;
        }
    }
    errorMessage << ".";
    throw std::runtime_error( errorMessage.str( ) );
}

std::vector< SumImageData > removeObservationsWithMissingLandmarks( const std::vector< SumImageData >& sumImages,
                                                                    const std::map< std::string, LmkLandmarkData >& landmarks )
{
    std::vector< SumImageData > filteredImages;
    for( const SumImageData& image : sumImages )
    {
        SumImageData filteredImage = image;
        filteredImage.landmarkObservations_.clear( );
        for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
        {
            if( landmarks.count( observation.landmarkId_ ) != 0 )
            {
                filteredImage.landmarkObservations_.push_back( observation );
            }
        }

        if( !filteredImage.landmarkObservations_.empty( ) )
        {
            filteredImages.push_back( filteredImage );
        }
    }
    return filteredImages;
}

//! Check that every landmark of an image lies in front of the camera, using the SUM SCOBJ vector.
void validateSumImageGeometryWithScObj( const SumImageData& image, const std::map< std::string, LmkLandmarkData >& landmarks )
{
    if( !image.spacecraftObjectVector_.array( ).isFinite( ).all( ) )
    {
        return;
    }

    for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
    {
        // SCOBJ points from the spacecraft to the target centre, so the spacecraft sits at -SCOBJ. The landmark
        // position relative to the spacecraft is therefore VLM + SCOBJ.
        const Eigen::Vector3d relativeCameraFramePosition =
                image.cameraAxes_ * ( landmarks.at( observation.landmarkId_ ).bodyFixedPosition_ + image.spacecraftObjectVector_ );
        if( !relativeCameraFramePosition.array( ).isFinite( ).all( ) ||
            relativeCameraFramePosition.z( ) <= std::numeric_limits< double >::epsilon( ) )
        {
            throw std::runtime_error( "Error when converting SUM image '" + image.imageId_ + "': landmark '" + observation.landmarkId_ +
                                      "' projects to non-positive camera-frame z using CX/CY/CZ and SCOBJ." );
        }
    }
}

}  // namespace

std::string getSumImageCameraName( const std::string& imageId )
{
    const std::string trimmedImageId = trimCopy( imageId );
    std::string sanitizedId;
    sanitizedId.reserve( trimmedImageId.size( ) );
    for( const char character : trimmedImageId )
    {
        const bool isAllowed = std::isalnum( static_cast< unsigned char >( character ) ) || character == '_' || character == '-';
        sanitizedId.push_back( isAllowed ? character : '_' );
    }

    if( sanitizedId.empty( ) )
    {
        throw std::runtime_error( "Error when generating SUM camera name: image ID is empty after trimming." );
    }
    return "Camera_" + sanitizedId;
}

double getSumImageObservationTime( const SumImageData& image )
{
    const std::string errorPrefix = "Error when converting SUM UTC epoch for image '" + image.imageId_ + "': ";
    try
    {
        return parseTrackingDateTime( image.utcEpochString_ ).epoch< double >( );
    }
    catch( const std::exception& exception )
    {
        throw std::runtime_error( errorPrefix + exception.what( ) );
    }
}

Eigen::Quaterniond getSumImageRotationFromBodyFixedToCameraFrame( const SumImageData& image )
{
    return Eigen::Quaterniond( image.cameraAxes_ ).normalized( );
}

std::vector< SumImageData > getSumImagesToConvert( const std::vector< SumImageData >& sumImages,
                                                   const std::map< std::string, LmkLandmarkData >& landmarks,
                                                   const bool skipObservationsWithMissingLandmarks,
                                                   const bool validateSpacecraftObjectGeometry )
{
    std::vector< SumImageData > imagesToConvert = sumImages;
    if( skipObservationsWithMissingLandmarks )
    {
        imagesToConvert = removeObservationsWithMissingLandmarks( sumImages, landmarks );
    }
    else
    {
        validateReferencedLandmarksHaveDefinitions( imagesToConvert, landmarks );
    }

    std::map< std::string, std::string > cameraNameToImageId;
    for( const SumImageData& image : imagesToConvert )
    {
        validateSumImageCalibrationAndAttitude( image );
        getSumImageObservationTime( image );
        if( validateSpacecraftObjectGeometry )
        {
            validateSumImageGeometryWithScObj( image, landmarks );
        }

        const std::string cameraName = getSumImageCameraName( image.imageId_ );
        const auto inserted = cameraNameToImageId.emplace( cameraName, image.imageId_ );
        if( !inserted.second && inserted.first->second != image.imageId_ )
        {
            throw std::runtime_error( "Error when converting SUM/LMK observations: image IDs '" + inserted.first->second + "' and '" +
                                      image.imageId_ + "' both map to camera name '" + cameraName + "'." );
        }
    }
    return mergeSumImages( imagesToConvert );
}

std::vector< std::shared_ptr< data::TrackingSupplementaryData > > getSumLmkTrackingSupplementaryData(
        const std::vector< SumImageData >& sumImages,
        const std::map< std::string, LmkLandmarkData >& landmarks,
        const std::string& targetBodyName,
        const std::string& receiverBodyName )
{
    std::vector< std::shared_ptr< data::TrackingSupplementaryData > > supplementaryData;
    for( const auto& landmarkEntry : landmarks )
    {
        auto landmarkData = std::make_shared< data::TrackingSupplementaryData >( targetBodyName, landmarkEntry.second.landmarkId_ );
        landmarkData->setReferencePointSupplementaryData(
                std::make_shared< data::ReferencePointSupplementaryData >( landmarkEntry.second.bodyFixedPosition_ ) );
        supplementaryData.push_back( landmarkData );
    }

    for( const SumImageData& image : sumImages )
    {
        const std::string cameraName = getSumImageCameraName( image.imageId_ );
        auto camera = std::make_shared< data::CameraInstrumentSupplementaryData >(
                cameraName,
                image.focalLengthMm_,
                image.opticalCenter_,
                Eigen::Vector4d( 0.0, static_cast< double >( image.imageSize_( 0 ) ), 0.0, static_cast< double >( image.imageSize_( 1 ) ) ),
                image.kMatrix_,
                Eigen::Matrix< double, 6, 1 >::Zero( ),
                Eigen::Vector3d::Zero( ),
                std::map< double, Eigen::Quaterniond >(
                        { { getSumImageObservationTime( image ), getSumImageRotationFromBodyFixedToCameraFrame( image ) } } ) );
        camera->setPointingReferenceBodyName( targetBodyName );
        camera->setPointingUncertainty( image.pointingSigma_ );

        auto cameraData = std::make_shared< data::TrackingSupplementaryData >( receiverBodyName, cameraName );
        cameraData->setInstrumentSupplementaryData( { camera } );
        supplementaryData.push_back( cameraData );
    }
    return supplementaryData;
}

}  // namespace sum_lmk

}  // namespace input_output

}  // namespace tudat
