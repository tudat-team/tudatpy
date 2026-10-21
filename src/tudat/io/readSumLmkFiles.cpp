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
#include <stdexcept>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/LU>

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

//! Tolerance for the orthonormality and determinant checks on a SUM CX/CY/CZ attitude matrix. Valid matrices in the SPC
//! archive have a Frobenius orthogonality error just below 8e-8, due to their stored numerical precision.
constexpr double sumCameraRotationMatrixTolerance = 1.0E-7;

enum class SumSection { none, landmarks, limbFits };
enum class LmkSection { none, pictures, mapOverlaps, limbFits };

//! Split a line into whitespace/comma separated tokens.
std::vector< std::string > tokenize( const std::string& line )
{
    std::vector< std::string > tokens;
    const std::string trimmedLine = boost::algorithm::trim_copy_if( line, boost::algorithm::is_any_of( " \t," ) );
    if( !trimmedLine.empty( ) )
    {
        boost::algorithm::split( tokens, trimmedLine, boost::algorithm::is_any_of( " \t," ), boost::algorithm::token_compress_on );
    }
    return tokens;
}

std::string getContext( const std::string& file, const int lineNumber, const std::string& id )
{
    return file + ":" + std::to_string( lineNumber ) + " (" + id + ")";
}

//! Parse a floating-point value, which may use the Fortran D exponent notation.
double parseDouble( std::string token, const std::string& context )
{
    std::replace( token.begin( ), token.end( ), 'D', 'E' );
    std::replace( token.begin( ), token.end( ), 'd', 'E' );
    try
    {
        return boost::lexical_cast< double >( token );
    }
    catch( const boost::bad_lexical_cast& )
    {
        throw std::runtime_error( "Error when parsing floating-point value '" + token + "' in " + context + "." );
    }
}

int parseInt( const std::string& token, const std::string& context )
{
    try
    {
        return boost::lexical_cast< int >( token );
    }
    catch( const boost::bad_lexical_cast& )
    {
        throw std::runtime_error( "Error when parsing integer value '" + token + "' in " + context + "." );
    }
}

Eigen::Vector3d parseVector3( const std::vector< std::string >& tokens, const std::string& context, const double scale = 1.0 )
{
    return scale *
            Eigen::Vector3d( parseDouble( tokens.at( 0 ), context ),
                             parseDouble( tokens.at( 1 ), context ),
                             parseDouble( tokens.at( 2 ), context ) );
}

//! Parse a one- or three-valued uncertainty record; a single value applies to all three components.
Eigen::Vector3d parseUncertainty( const std::vector< std::string >& tokens, const std::string& context, const double scale = 1.0 )
{
    if( tokens.size( ) == 1 )
    {
        return Eigen::Vector3d::Constant( scale * parseDouble( tokens.at( 0 ), context ) );
    }
    return parseVector3( tokens, context, scale );
}

//! Check whether a record ends with the given (case-insensitive) label and, if so, retrieve the values preceding it.
bool getRecordValues( const std::string& trimmedLine,
                      const std::string& label,
                      const std::string& context,
                      const unsigned int minimumNumberOfValues,
                      const unsigned int maximumNumberOfValues,
                      std::vector< std::string >& values )
{
    if( trimmedLine.size( ) <= label.size( ) || !boost::algorithm::iends_with( trimmedLine, label ) ||
        !std::isspace( static_cast< unsigned char >( trimmedLine.at( trimmedLine.size( ) - label.size( ) - 1 ) ) ) )
    {
        return false;
    }

    values = tokenize( trimmedLine.substr( 0, trimmedLine.size( ) - label.size( ) ) );
    if( values.size( ) < minimumNumberOfValues || values.size( ) > maximumNumberOfValues )
    {
        throw std::runtime_error( "Error when parsing record '" + label + "' in " + context + ": found " +
                                  std::to_string( values.size( ) ) + " values." );
    }
    return true;
}

SumLandmarkObservation parseSumLandmarkRow( const std::vector< std::string >& tokens, const std::string& context )
{
    // Production SUM files may carry an optional trailing '-' flag.
    if( tokens.size( ) != 3 && !( tokens.size( ) == 4 && tokens.at( 3 ) == "-" ) )
    {
        throw std::runtime_error( "Error when parsing SUM LANDMARKS row in " + context + ": expected 3 fields plus optional '-' flag." );
    }
    SumLandmarkObservation observation;
    observation.landmarkId_ = tokens.at( 0 );
    observation.pixelCoordinates_ = Eigen::Vector2d( parseDouble( tokens.at( 1 ), context ), parseDouble( tokens.at( 2 ), context ) );
    return observation;
}

SumLimbFitObservation parseSumLimbFitRow( const std::vector< std::string >& tokens, const std::string& context )
{
    if( tokens.size( ) != 4 )
    {
        throw std::runtime_error( "Error when parsing SUM LIMB FITS row in " + context + ": expected 4 fields." );
    }
    SumLimbFitObservation observation;
    observation.featureId_ = tokens.at( 0 );
    observation.pixelCoordinates_ = Eigen::Vector2d( parseDouble( tokens.at( 1 ), context ), parseDouble( tokens.at( 2 ), context ) );
    observation.sigma_ = parseDouble( tokens.at( 3 ), context );
    return observation;
}

//! Check the camera calibration, attitude, pointing uncertainty and measurements of a SUM image before conversion.
void validateSumImage( const SumImageData& image )
{
    const std::string errorPrefix = "Error when converting SUM image " + image.imageId_ + ": ";
    if( !( image.focalLengthMm_ > 0.0 ) || !image.opticalCenter_.allFinite( ) || !image.kMatrix_.allFinite( ) ||
        !image.kMatrix_.block< 2, 2 >( 0, 0 ).fullPivLu( ).isInvertible( ) )
    {
        throw std::runtime_error( errorPrefix + "MMFL must be positive, CTR finite, and the K-MATRIX finite and invertible." );
    }
    if( !image.cameraAxes_.allFinite( ) ||
        ( image.cameraAxes_ * image.cameraAxes_.transpose( ) - Eigen::Matrix3d::Identity( ) ).norm( ) > sumCameraRotationMatrixTolerance ||
        std::fabs( image.cameraAxes_.determinant( ) - 1.0 ) > sumCameraRotationMatrixTolerance )
    {
        throw std::runtime_error( errorPrefix + "CX/CY/CZ must form a rotation matrix." );
    }
    const bool hasPointingSigma = !image.pointingSigma_.array( ).isNaN( ).all( );
    if( hasPointingSigma && !( image.pointingSigma_.allFinite( ) && ( image.pointingSigma_.array( ) > 0.0 ).all( ) ) )
    {
        throw std::runtime_error( errorPrefix + "SIGMA_PTG must contain three finite positive values, or be absent." );
    }
    // The SUM DISTORTION row is not a full distortion model: SPC keeps it zero and supplies a non-zero (Owen) model separately
    // in INIT_LITHOS. A non-zero row is rejected to avoid silently simulating a different camera.
    if( !image.distortionCoefficients_.isZero( 0.0 ) )
    {
        throw std::runtime_error( errorPrefix + "non-zero SUM DISTORTION coefficients are not supported." );
    }
    for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
    {
        if( !observation.pixelCoordinates_.allFinite( ) )
        {
            throw std::runtime_error( errorPrefix + "pixel coordinates of landmark " + observation.landmarkId_ + " must be finite." );
        }
    }
}

//! Check that every landmark of an image lies in front of the camera, using the SUM SCOBJ vector.
void validateSumImageGeometry( const SumImageData& image, const std::map< std::string, LmkLandmarkData >& landmarks )
{
    if( !image.spacecraftObjectVector_.allFinite( ) )
    {
        return;
    }

    for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
    {
        // SCOBJ points from the spacecraft to the target centre, so the landmark position relative to the spacecraft is VLM + SCOBJ.
        const Eigen::Vector3d cameraFramePosition =
                image.cameraAxes_ * ( landmarks.at( observation.landmarkId_ ).bodyFixedPosition_ + image.spacecraftObjectVector_ );
        if( !( cameraFramePosition.z( ) > 0.0 ) )
        {
            throw std::runtime_error( "Error when converting SUM image " + image.imageId_ + ": landmark " + observation.landmarkId_ +
                                      " is behind the camera, using CX/CY/CZ and SCOBJ." );
        }
    }
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
    bool hasImageSize = false, hasCalibration = false, hasCx = false, hasCy = false, hasCz = false, hasKMatrix = false;

    std::string line;
    int lineNumber = 0;
    while( std::getline( stream, line ) )
    {
        ++lineNumber;
        const std::string trimmedLine = boost::algorithm::trim_copy( line );
        if( trimmedLine.empty( ) )
        {
            continue;
        }

        // The first two lines hold the image ID and the UTC epoch.
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

        if( boost::algorithm::iequals( trimmedLine, "END FILE" ) )
        {
            break;
        }
        else if( boost::algorithm::iequals( trimmedLine, "LANDMARKS" ) )
        {
            activeSection = SumSection::landmarks;
            continue;
        }
        else if( boost::algorithm::iequals( trimmedLine, "LIMB FITS" ) )
        {
            activeSection = SumSection::limbFits;
            continue;
        }

        const std::string context = getContext( sumFile, lineNumber, image.imageId_ );
        if( activeSection == SumSection::landmarks )
        {
            image.landmarkObservations_.push_back( parseSumLandmarkRow( tokenize( trimmedLine ), context ) );
            continue;
        }
        else if( activeSection == SumSection::limbFits )
        {
            image.limbFitObservations_.push_back( parseSumLimbFitRow( tokenize( trimmedLine ), context ) );
            continue;
        }

        std::vector< std::string > values;
        auto isRecord = [ & ]( const std::string& label, const unsigned int minimumSize, const unsigned int maximumSize ) {
            return getRecordValues( trimmedLine, label, context, minimumSize, maximumSize, values );
        };
        if( isRecord( "NPX, NLN, THRSH", 4, 4 ) )
        {
            image.imageSize_ = Eigen::Vector2i( parseInt( values.at( 0 ), context ), parseInt( values.at( 1 ), context ) );
            image.threshold_ = parseInt( values.at( 2 ), context );
            image.maxDn_ = parseInt( values.at( 3 ), context );
            hasImageSize = true;
        }
        else if( isRecord( "MMFL, CTR", 3, 3 ) )
        {
            image.focalLengthMm_ = parseDouble( values.at( 0 ), context );
            image.opticalCenter_ = Eigen::Vector2d( parseDouble( values.at( 1 ), context ), parseDouble( values.at( 2 ), context ) );
            hasCalibration = true;
        }
        else if( isRecord( "SCOBJ", 3, 3 ) )
        {
            image.spacecraftObjectVector_ = parseVector3( values, context, kilometersToMeters );
        }
        else if( isRecord( "CX", 3, 3 ) )
        {
            image.cameraAxes_.row( 0 ) = parseVector3( values, context ).transpose( );
            hasCx = true;
        }
        else if( isRecord( "CY", 3, 3 ) )
        {
            image.cameraAxes_.row( 1 ) = parseVector3( values, context ).transpose( );
            hasCy = true;
        }
        else if( isRecord( "CZ", 3, 3 ) )
        {
            image.cameraAxes_.row( 2 ) = parseVector3( values, context ).transpose( );
            hasCz = true;
        }
        else if( isRecord( "SZ", 3, 3 ) )
        {
            image.sunDirectionBodyFixed_ = parseVector3( values, context );
        }
        else if( isRecord( "K-MATRIX", 6, 6 ) )
        {
            for( int i = 0; i < 6; ++i )
            {
                image.kMatrix_( i / 3, i % 3 ) = parseDouble( values.at( i ), context );
            }
            hasKMatrix = true;
        }
        else if( isRecord( "SIGMA_PTG", 1, 3 ) )
        {
            image.pointingSigma_ = parseUncertainty( values, context );
        }
        else if( isRecord( "SIGMA_VSO", 1, 3 ) )
        {
            image.spacecraftObjectSigma_ = parseUncertainty( values, context, kilometersToMeters );
        }
        else if( isRecord( "DISTORTION", 4, 4 ) || isRecord( "DISTORTION COEFFICIENTS", 4, 4 ) || isRecord( "DCOEF", 4, 4 ) )
        {
            for( int i = 0; i < 4; ++i )
            {
                image.distortionCoefficients_( i ) = parseDouble( values.at( i ), context );
            }
        }
    }

    if( image.utcEpochString_.empty( ) || !hasImageSize || !hasCalibration || !hasCx || !hasCy || !hasCz || !hasKMatrix )
    {
        throw std::runtime_error( "Error when reading SUM file " + sumFile +
                                  ": image ID, epoch, NPX/NLN/THRSH, MMFL/CTR, CX/CY/CZ and K-MATRIX records are required." );
    }
    return image;
}

std::vector< SumImageData > readSumFiles( const std::vector< std::string >& sumFiles )
{
    std::vector< SumImageData > images;
    for( const std::string& sumFile : sumFiles )
    {
        images.push_back( readSumFile( sumFile ) );
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
        const std::string trimmedLine = boost::algorithm::trim_copy( line );
        if( trimmedLine.empty( ) )
        {
            continue;
        }

        // The first line holds the landmark ID, optionally followed by a type flag.
        if( landmark.landmarkId_.empty( ) )
        {
            const std::vector< std::string > headerTokens = tokenize( trimmedLine );
            landmark.landmarkId_ = headerTokens.at( 0 );
            if( headerTokens.size( ) > 1 )
            {
                landmark.typeFlag_ = headerTokens.at( 1 ).at( 0 );
            }
            continue;
        }

        if( boost::algorithm::iequals( trimmedLine, "END FILE" ) )
        {
            break;
        }
        else if( boost::algorithm::iequals( trimmedLine, "PICTURES" ) )
        {
            activeSection = LmkSection::pictures;
            continue;
        }
        else if( boost::algorithm::iequals( trimmedLine, "MAP OVERLAPS" ) )
        {
            activeSection = LmkSection::mapOverlaps;
            continue;
        }
        else if( boost::algorithm::iequals( trimmedLine, "LIMB FITS" ) )
        {
            activeSection = LmkSection::limbFits;
            continue;
        }

        const std::string context = getContext( lmkFile, lineNumber, landmark.landmarkId_ );
        if( activeSection == LmkSection::pictures )
        {
            const std::vector< std::string > tokens = tokenize( trimmedLine );
            if( tokens.size( ) != 3 && tokens.size( ) != 4 )
            {
                throw std::runtime_error( "Error when parsing LMK PICTURES row in " + context +
                                          ": expected 3 fields plus optional '*' flag." );
            }
            LmkPictureObservation picture;
            picture.imageId_ = tokens.at( 0 );
            picture.pixelCoordinates_ = Eigen::Vector2d( parseDouble( tokens.at( 1 ), context ), parseDouble( tokens.at( 2 ), context ) );
            picture.flagged_ = ( tokens.size( ) == 4 && tokens.at( 3 ) == "*" );
            landmark.pictures_.push_back( picture );
            continue;
        }
        else if( activeSection != LmkSection::none )
        {
            // MAP OVERLAPS and LIMB FITS rows are not used.
            continue;
        }

        std::vector< std::string > values;
        auto isRecord = [ & ]( const std::string& label, const unsigned int numberOfValues ) {
            return getRecordValues( trimmedLine, label, context, numberOfValues, numberOfValues, values );
        };
        // Combined labels must be matched before the single ones, since e.g. "... SIGKM, RMSLMK" also ends with "RMSLMK".
        if( isRecord( "SIZE, SCALE(KM)", 2 ) )
        {
            landmark.patchSize_ = parseInt( values.at( 0 ), context );
            landmark.patchScaleKm_ = parseDouble( values.at( 1 ), context );
        }
        else if( isRecord( "SIGKM, RMSLMK", 2 ) )
        {
            landmark.sigKm_ = parseDouble( values.at( 0 ), context );
            landmark.rmsLmk_ = parseDouble( values.at( 1 ), context );
        }
        else if( isRecord( "SIZE", 1 ) )
        {
            landmark.patchSize_ = parseInt( values.at( 0 ), context );
        }
        else if( isRecord( "SCALE", 1 ) )
        {
            landmark.patchScaleKm_ = parseDouble( values.at( 0 ), context );
        }
        else if( isRecord( "HORIZON", 4 ) )
        {
            for( int i = 0; i < 4; ++i )
            {
                landmark.horizonFlags_( i ) = parseInt( values.at( i ), context );
            }
        }
        else if( isRecord( "SIGKM", 1 ) )
        {
            landmark.sigKm_ = parseDouble( values.at( 0 ), context );
        }
        else if( isRecord( "RMSLMK", 1 ) )
        {
            landmark.rmsLmk_ = parseDouble( values.at( 0 ), context );
        }
        else if( isRecord( "VLM", 3 ) )
        {
            landmark.bodyFixedPosition_ = parseVector3( values, context, kilometersToMeters );
        }
        else if( isRecord( "UX", 3 ) )
        {
            landmark.localXAxis_ = parseVector3( values, context );
        }
        else if( isRecord( "UY", 3 ) )
        {
            landmark.localYAxis_ = parseVector3( values, context );
        }
        else if( isRecord( "UZ", 3 ) )
        {
            landmark.localZAxis_ = parseVector3( values, context );
        }
        else if( getRecordValues( trimmedLine, "SIGMA_LMK", context, 1, 3, values ) )
        {
            landmark.landmarkPositionSigma_ = parseUncertainty( values, context, kilometersToMeters );
        }
    }

    if( !landmark.bodyFixedPosition_.allFinite( ) )
    {
        throw std::runtime_error( "Error when reading LMK file " + lmkFile + ": landmark ID and VLM record are required." );
    }
    return landmark;
}

std::map< std::string, LmkLandmarkData > readLmkFiles( const std::vector< std::string >& lmkFiles )
{
    std::map< std::string, LmkLandmarkData > landmarks;
    for( const std::string& lmkFile : lmkFiles )
    {
        const LmkLandmarkData landmark = readLmkFile( lmkFile );
        if( !landmarks.emplace( landmark.landmarkId_, landmark ).second )
        {
            throw std::runtime_error( "Error when reading LMK files: landmark " + landmark.landmarkId_ + " is provided more than once." );
        }
    }
    return landmarks;
}

double getSumImageObservationTime( const SumImageData& image )
{
    return parseTrackingDateTime( image.utcEpochString_ ).epoch< double >( );
}

Eigen::Quaterniond getSumImageRotationFromBodyFixedToCameraFrame( const SumImageData& image )
{
    return Eigen::Quaterniond( image.cameraAxes_ ).normalized( );
}

data::CameraPictureCalibration getSumImageCalibration( const SumImageData& image )
{
    data::CameraPictureCalibration calibration;
    calibration.focalLength_ = image.focalLengthMm_;
    calibration.principalPoint_ = image.opticalCenter_;
    calibration.fieldOfViewBounds_ =
            Eigen::Vector4d( 0.0, static_cast< double >( image.imageSize_( 0 ) ), 0.0, static_cast< double >( image.imageSize_( 1 ) ) );
    calibration.kMatrix_ = image.kMatrix_;
    return calibration;
}

std::vector< SumImageData > getSumImagesToConvert( const std::vector< SumImageData >& sumImages,
                                                   const std::map< std::string, LmkLandmarkData >& landmarks,
                                                   const bool skipObservationsWithMissingLandmarks,
                                                   const bool validateSpacecraftObjectGeometry )
{
    std::vector< SumImageData > imagesToConvert;
    std::set< std::string > imageIds;
    for( SumImageData image : sumImages )
    {
        if( !imageIds.insert( image.imageId_ ).second )
        {
            throw std::runtime_error( "Error when converting SUM/LMK data: image " + image.imageId_ + " is provided more than once." );
        }
        validateSumImage( image );

        std::vector< SumLandmarkObservation > landmarkObservations;
        for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
        {
            if( landmarks.count( observation.landmarkId_ ) != 0 )
            {
                landmarkObservations.push_back( observation );
            }
            else if( !skipObservationsWithMissingLandmarks )
            {
                throw std::runtime_error( "Error when converting SUM/LMK data: no LMK data for landmark " + observation.landmarkId_ +
                                          ", observed in image " + image.imageId_ + "." );
            }
        }
        image.landmarkObservations_ = landmarkObservations;

        if( validateSpacecraftObjectGeometry )
        {
            validateSumImageGeometry( image, landmarks );
        }
        if( !image.landmarkObservations_.empty( ) )
        {
            imagesToConvert.push_back( image );
        }
    }

    // All images are taken by one camera, so they must have distinct epochs. Their calibration may differ (for instance by pixel
    // binning), and is kept per image.
    std::sort( imagesToConvert.begin( ), imagesToConvert.end( ), []( const SumImageData& first, const SumImageData& second ) {
        return getSumImageObservationTime( first ) < getSumImageObservationTime( second );
    } );
    for( unsigned int i = 1; i < imagesToConvert.size( ); ++i )
    {
        const SumImageData& firstImage = imagesToConvert.front( );
        const SumImageData& image = imagesToConvert.at( i );
        if( image.pointingSigma_.array( ).isNaN( ).all( ) != firstImage.pointingSigma_.array( ).isNaN( ).all( ) ||
            !( getSumImageObservationTime( image ) > getSumImageObservationTime( imagesToConvert.at( i - 1 ) ) ) )
        {
            throw std::runtime_error( "Error when converting SUM image " + image.imageId_ +
                                      ": all images must have distinct epochs, and either all or none must give SIGMA_PTG." );
        }
    }
    return imagesToConvert;
}

std::vector< std::shared_ptr< data::TrackingSupplementaryData > > getSumLmkTrackingSupplementaryData(
        const std::vector< SumImageData >& sumImages,
        const std::map< std::string, LmkLandmarkData >& landmarks,
        const std::string& targetBodyName,
        const std::string& receiverBodyName,
        const std::string& cameraName )
{
    std::vector< std::shared_ptr< data::TrackingSupplementaryData > > supplementaryData;
    for( const auto& landmarkEntry : landmarks )
    {
        auto landmarkData = std::make_shared< data::TrackingSupplementaryData >( targetBodyName, landmarkEntry.second.landmarkId_ );
        landmarkData->setReferencePointSupplementaryData(
                std::make_shared< data::ReferencePointSupplementaryData >( landmarkEntry.second.bodyFixedPosition_ ) );
        supplementaryData.push_back( landmarkData );
    }

    if( sumImages.empty( ) )
    {
        return supplementaryData;
    }

    std::map< double, Eigen::Quaterniond > pointingHistory;
    std::map< double, data::CameraPictureCalibration > pictureCalibrationHistory;
    std::map< double, std::string > pictureIdHistory;
    std::map< double, Eigen::Vector3d > pointingUncertaintyHistory;
    for( const SumImageData& image : sumImages )
    {
        pointingHistory[ getSumImageObservationTime( image ) ] = getSumImageRotationFromBodyFixedToCameraFrame( image );
        pictureCalibrationHistory[ getSumImageObservationTime( image ) ] = getSumImageCalibration( image );
        pictureIdHistory[ getSumImageObservationTime( image ) ] = image.imageId_;
        if( !image.pointingSigma_.array( ).isNaN( ).all( ) )
        {
            pointingUncertaintyHistory[ getSumImageObservationTime( image ) ] = image.pointingSigma_;
        }
    }

    // Each image has its own calibration. The camera-level calibration is that of the first image.
    const data::CameraPictureCalibration& firstCalibration = pictureCalibrationHistory.begin( )->second;
    auto camera = std::make_shared< data::CameraInstrumentSupplementaryData >( cameraName,
                                                                               firstCalibration.focalLength_,
                                                                               firstCalibration.principalPoint_,
                                                                               firstCalibration.fieldOfViewBounds_,
                                                                               firstCalibration.kMatrix_,
                                                                               Eigen::Matrix< double, 6, 1 >::Zero( ),
                                                                               Eigen::Vector3d::Zero( ),
                                                                               pointingHistory );
    camera->setPictureCalibrationHistory( pictureCalibrationHistory );
    camera->setPointingReferenceBodyName( targetBodyName );
    camera->setPictureIdHistory( pictureIdHistory );
    camera->setPointingUncertaintyHistory( pointingUncertaintyHistory );

    auto cameraData = std::make_shared< data::TrackingSupplementaryData >( receiverBodyName, cameraName );
    cameraData->setInstrumentSupplementaryData( { camera } );
    supplementaryData.push_back( cameraData );
    return supplementaryData;
}

}  // namespace sum_lmk

}  // namespace input_output

}  // namespace tudat
