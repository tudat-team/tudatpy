/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROCESS_SUM_LMK_FILES_H
#define TUDAT_PROCESS_SUM_LMK_FILES_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <memory>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/io/readSumLmkFiles.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/astro/system_models/camera.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"
#include "tudat/simulation/estimation_setup/estimatableParameterSettings.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/processPsfFile.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/singleObservationSet.h"

namespace tudat
{

namespace observation_models
{

class SumLmkObservationConversionSettings
{
public:
    SumLmkObservationConversionSettings( const std::string& targetBodyName, const std::string& receiverBodyName ):
        targetBodyName_( targetBodyName ), receiverBodyName_( receiverBodyName )
    {}

    std::string targetBodyName_;
    std::string receiverBodyName_;
    Eigen::Vector3d bodyFixedCameraPosition_ = Eigen::Vector3d::Zero( );
    bool validateSpacecraftObjectGeometry_ = true;
    bool skipObservationsWithMissingLandmarks_ = false;
};

template< typename ObservationScalarType = double, typename TimeType = double >
struct SumLmkObservationConversionResult {
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection_;
    std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
            inverseAprioriCovarianceDiagonalEntries_;
    std::map< std::string, std::string > imageIdToCameraName_;
    // Body carrying the per-image cameras, i.e. the body the pointing-correction parameters belong to.
    std::string receiverBodyName_;
    // Observation model settings matching observationCollection_ (one pixel_coordinates setting per
    // (image, landmark) link end), ready for createObservationSimulators / residual computation.
    std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings_;
};

namespace detail
{

inline std::string trimCopyForSumLmk( const std::string& input )
{
    const std::string whitespace = " \t\r\n";
    const std::size_t first = input.find_first_not_of( whitespace );
    if( first == std::string::npos )
    {
        return "";
    }
    const std::size_t last = input.find_last_not_of( whitespace );
    return input.substr( first, last - first + 1 );
}

inline std::string sanitizeSumImageIdToCameraName( const std::string& imageId )
{
    const std::string trimmedImageId = trimCopyForSumLmk( imageId );
    std::string sanitizedId;
    sanitizedId.reserve( trimmedImageId.size( ) );
    for( const char character : trimmedImageId )
    {
        const unsigned char unsignedCharacter = static_cast< unsigned char >( character );
        if( std::isalnum( unsignedCharacter ) || character == '_' || character == '-' )
        {
            sanitizedId.push_back( character );
        }
        else
        {
            sanitizedId.push_back( '_' );
        }
    }

    if( sanitizedId.empty( ) )
    {
        throw std::runtime_error( "Error when generating SUM camera name: image ID is empty after trimming." );
    }
    return "Camera_" + sanitizedId;
}

inline bool isFiniteVector( const Eigen::Vector3d& vector )
{
    return vector.array( ).isFinite( ).all( );
}

inline bool isCloseVector( const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs, const double tolerance = 1.0E-9 )
{
    return ( lhs - rhs ).norm( ) <= tolerance * std::max( 1.0, std::max( lhs.norm( ), rhs.norm( ) ) );
}

template< typename TimeType >
TimeType convertSumUtcStringToSecondsSinceJ2000( const input_output::sum_lmk::SumImageData& image )
{
    try
    {
        return convertPsfUtcStringToSecondsSinceJ2000< TimeType >( image.utcEpochString_ );
    }
    catch( const std::exception& exception )
    {
        throw std::runtime_error( "Error when converting SUM UTC epoch for image '" + image.imageId_ +
                                  "': " + std::string( exception.what( ) ) );
    }
}

inline void validateReceiverAndTargetBodies( const simulation_setup::SystemOfBodies& bodies,
                                             const SumLmkObservationConversionSettings& conversionSettings )
{
    if( bodies.count( conversionSettings.targetBodyName_ ) == 0 )
    {
        throw std::runtime_error( "Error when converting SUM/LMK observations: target body '" + conversionSettings.targetBodyName_ +
                                  "' not found." );
    }
    if( bodies.count( conversionSettings.receiverBodyName_ ) == 0 )
    {
        throw std::runtime_error( "Error when converting SUM/LMK observations: receiver body '" + conversionSettings.receiverBodyName_ +
                                  "' not found." );
    }
    if( bodies.at( conversionSettings.targetBodyName_ )->getRotationalEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when converting SUM/LMK observations: target body '" + conversionSettings.targetBodyName_ +
                                  "' does not have a rotational ephemeris." );
    }
    if( bodies.at( conversionSettings.receiverBodyName_ )->getRotationalEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when converting SUM/LMK observations: receiver body '" + conversionSettings.receiverBodyName_ +
                                  "' does not have a rotational ephemeris." );
    }
}

inline void validateSanitizedCameraNames( const std::vector< input_output::sum_lmk::SumImageData >& sumImages,
                                          std::map< std::string, std::string >& imageIdToCameraName )
{
    std::map< std::string, std::string > cameraNameToImageId;
    for( const input_output::sum_lmk::SumImageData& image : sumImages )
    {
        const std::string cameraName = sanitizeSumImageIdToCameraName( image.imageId_ );
        if( cameraNameToImageId.count( cameraName ) != 0 && cameraNameToImageId.at( cameraName ) != image.imageId_ )
        {
            throw std::runtime_error( "Error when converting SUM/LMK observations: image IDs '" + cameraNameToImageId.at( cameraName ) +
                                      "' and '" + image.imageId_ + "' both map to camera name '" + cameraName + "'." );
        }
        cameraNameToImageId[ cameraName ] = image.imageId_;
        imageIdToCameraName[ image.imageId_ ] = cameraName;
    }
}

//! SUM's four-value DISTORTION row is not a complete camera-distortion model. In SPC, the row is
//! conventionally zero and a non-zero Owen distortion model is defined separately in INIT_LITHOS.
//! Rejecting non-zero values here avoids silently simulating a different camera from the one used
//! to produce the SPC measurements.
inline void validateSumImageDistortionCoefficients( const std::vector< input_output::sum_lmk::SumImageData >& sumImages )
{
    for( const input_output::sum_lmk::SumImageData& image : sumImages )
    {
        if( !image.distortionCoefficients_.array( ).isFinite( ).all( ) )
        {
            throw std::runtime_error( "Error when converting SUM image '" + image.imageId_ + "': DISTORTION coefficients must be finite." );
        }
        if( !image.distortionCoefficients_.isZero( 0.0 ) )
        {
            throw std::runtime_error( "Error when converting SUM image '" + image.imageId_ +
                                      "': non-zero SUM DISTORTION coefficients are unsupported. SPC supplies a non-zero Owen distortion "
                                      "model separately through INIT_LITHOS; it cannot be reconstructed from this SUM row." );
        }
    }
}

inline std::shared_ptr< system_models::PsfCameraProjectionModel > createSumCameraProjectionModel(
        const input_output::sum_lmk::SumImageData& image )
{
    // validateSumImageDistortionCoefficients has established that the SUM row is zero. The
    // PsfCameraProjectionModel implements the SUM focal-length/principal-point/K-MATRIX mapping.
    return std::make_shared< system_models::PsfCameraProjectionModel >(
            image.focalLengthMm_,
            image.opticalCenter_,
            image.kMatrix_,
            Eigen::Matrix< double, 6, 1 >::Zero( ),
            Eigen::Vector3d::Zero( ),
            Eigen::Vector4d( 0.0, static_cast< double >( image.imageSize_( 0 ) ), 0.0, static_cast< double >( image.imageSize_( 1 ) ) ) );
}

inline void addSumLmkLandmarksToBody( const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks,
                                      const std::shared_ptr< simulation_setup::Body >& targetBody )
{
    for( const auto& landmarkEntry : landmarks )
    {
        const input_output::sum_lmk::LmkLandmarkData& landmark = landmarkEntry.second;
        if( targetBody->getGroundStationMap( ).count( landmark.landmarkId_ ) == 0 )
        {
            simulation_setup::createGroundStation( targetBody,
                                                   landmark.landmarkId_,
                                                   landmark.bodyFixedPosition_,
                                                   coordinate_conversions::cartesian_position,
                                                   std::vector< std::shared_ptr< simulation_setup::GroundStationMotionSettings > >( ) );
        }
        else
        {
            const Eigen::Vector3d existingPosition =
                    targetBody->getGroundStation( landmark.landmarkId_ )->getNominalStationState( )->getNominalCartesianPosition( );
            if( !isCloseVector( existingPosition, landmark.bodyFixedPosition_ ) )
            {
                throw std::runtime_error( "Error when converting SUM/LMK observations: target body already has landmark/station '" +
                                          landmark.landmarkId_ + "' at a conflicting position." );
            }
        }
    }
}

inline void validateReferencedLandmarksHaveDefinitions( const std::vector< input_output::sum_lmk::SumImageData >& sumImages,
                                                        const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks )
{
    std::map< std::string, std::set< std::string > > missingLandmarkImages;
    for( const input_output::sum_lmk::SumImageData& image : sumImages )
    {
        for( const input_output::sum_lmk::SumLandmarkObservation& observation : image.landmarkObservations_ )
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
        if( !firstLandmark )
        {
            errorMessage << "; ";
        }
        firstLandmark = false;

        errorMessage << missingEntry.first << " referenced by image(s) ";
        bool firstImage = true;
        for( const std::string& imageId : missingEntry.second )
        {
            if( !firstImage )
            {
                errorMessage << ", ";
            }
            firstImage = false;
            errorMessage << imageId;
        }
    }
    errorMessage << ".";
    throw std::runtime_error( errorMessage.str( ) );
}

inline std::vector< input_output::sum_lmk::SumImageData > filterSumImagesForAvailableLandmarks(
        const std::vector< input_output::sum_lmk::SumImageData >& sumImages,
        const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks )
{
    std::vector< input_output::sum_lmk::SumImageData > filteredImages;
    for( const input_output::sum_lmk::SumImageData& image : sumImages )
    {
        input_output::sum_lmk::SumImageData filteredImage = image;
        filteredImage.landmarkObservations_.clear( );

        for( const input_output::sum_lmk::SumLandmarkObservation& observation : image.landmarkObservations_ )
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

inline void validateSumImageGeometryWithScoBj( const input_output::sum_lmk::SumImageData& image,
                                               const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks,
                                               const std::shared_ptr< system_models::PsfCameraProjectionModel >& projectionModel )
{
    if( !isFiniteVector( image.spacecraftObjectVector_ ) )
    {
        return;
    }

    for( const input_output::sum_lmk::SumLandmarkObservation& observation : image.landmarkObservations_ )
    {
        if( landmarks.count( observation.landmarkId_ ) == 0 )
        {
            throw std::runtime_error( "Error when converting SUM image '" + image.imageId_ + "': missing LMK data for landmark '" +
                                      observation.landmarkId_ + "'." );
        }
        // SCOBJ is the spacecraft-to-object (target-centre) vector in the target body-fixed frame,
        // so the spacecraft position relative to the target centre is -SCOBJ. The landmark position
        // relative to the spacecraft is therefore VLM - (-SCOBJ) = VLM + SCOBJ.
        const Eigen::Vector3d spacecraftBodyFixedPosition = -image.spacecraftObjectVector_;
        const Eigen::Vector3d relativeBodyFixedPosition =
                landmarks.at( observation.landmarkId_ ).bodyFixedPosition_ - spacecraftBodyFixedPosition;
        const Eigen::Vector3d relativeCameraFramePosition = image.cameraAxes_ * relativeBodyFixedPosition;
        if( relativeCameraFramePosition.z( ) <= std::numeric_limits< double >::epsilon( ) )
        {
            throw std::runtime_error( "Error when converting SUM image '" + image.imageId_ + "': landmark '" + observation.landmarkId_ +
                                      "' projects to non-positive camera-frame z using CX/CY/CZ and SCOBJ." );
        }
        projectionModel->projectUnitVectorToPixelLine( relativeCameraFramePosition );
    }
}

inline void addSumLmkCamerasToBody( const std::vector< input_output::sum_lmk::SumImageData >& sumImages,
                                    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks,
                                    const std::shared_ptr< simulation_setup::Body >& receiverBody,
                                    const std::shared_ptr< simulation_setup::Body >& targetBody,
                                    const SumLmkObservationConversionSettings& conversionSettings,
                                    const std::map< std::string, std::string >& imageIdToCameraName )
{
    for( const input_output::sum_lmk::SumImageData& image : sumImages )
    {
        const std::string cameraName = imageIdToCameraName.at( image.imageId_ );
        if( receiverBody->getVehicleSystems( )->getCameraMap( ).count( cameraName ) != 0 )
        {
            throw std::runtime_error( "Error when converting SUM/LMK observations: receiver body '" + receiverBody->getBodyName( ) +
                                      "' already has camera '" + cameraName + "'." );
        }

        std::shared_ptr< system_models::PsfCameraProjectionModel > projectionModel = createSumCameraProjectionModel( image );
        if( conversionSettings.validateSpacecraftObjectGeometry_ )
        {
            validateSumImageGeometryWithScoBj( image, landmarks, projectionModel );
        }

        const Eigen::Matrix3d rotationFromTargetBodyFixedToCamera = image.cameraAxes_;
        const std::shared_ptr< ephemerides::RotationalEphemeris > targetRotationalEphemeris = targetBody->getRotationalEphemeris( );
        // Nominal (uncorrected) picture-specific pointing; the estimated pointing correction is applied on top of
        // this by the Camera itself, so that it is applied identically for every way the camera frame is reached.
        std::function< Eigen::Quaterniond( const double ) > rotationFromInertialToCameraFrameFunction =
                [ rotationFromTargetBodyFixedToCamera, targetRotationalEphemeris ]( const double time ) {
                    return Eigen::Quaterniond( rotationFromTargetBodyFixedToCamera *
                                               targetRotationalEphemeris->getRotationToTargetFrame( time ).toRotationMatrix( ) )
                            .normalized( );
                };

        std::shared_ptr< system_models::Camera > camera =
                std::make_shared< system_models::Camera >( cameraName,
                                                           Eigen::Quaterniond( rotationFromTargetBodyFixedToCamera ),
                                                           projectionModel,
                                                           rotationFromInertialToCameraFrameFunction );
        receiverBody->getVehicleSystems( )->addCamera( cameraName, camera, conversionSettings.bodyFixedCameraPosition_ );
    }
}

inline void addPointingAprioriEntryIfAvailable(
        const input_output::sum_lmk::SumImageData& image,
        const std::string& receiverBodyName,
        const std::string& cameraName,
        std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                inverseAprioriCovarianceDiagonalEntries )
{
    const bool hasAnyFiniteSigma = image.pointingSigma_.array( ).isFinite( ).any( );
    if( !hasAnyFiniteSigma )
    {
        return;
    }
    if( !image.pointingSigma_.array( ).isFinite( ).all( ) )
    {
        throw std::runtime_error( "Error when converting SUM image '" + image.imageId_ +
                                  "': SIGMA_PTG must provide either all three components or no components." );
    }

    Eigen::VectorXd inverseVariance = Eigen::VectorXd::Zero( 3 );
    for( int i = 0; i < 3; ++i )
    {
        if( image.pointingSigma_( i ) <= 0.0 )
        {
            throw std::runtime_error( "Error when converting SUM image '" + image.imageId_ + "': SIGMA_PTG entries must be positive." );
        }
        inverseVariance( i ) = 1.0 / ( image.pointingSigma_( i ) * image.pointingSigma_( i ) );
    }

    inverseAprioriCovarianceDiagonalEntries.push_back(
            std::make_pair( estimatable_parameters::EstimatebleParameterIdentifier( estimatable_parameters::camera_pointing_correction,
                                                                                    std::make_pair( receiverBodyName, cameraName ) ),
                            inverseVariance ) );
}

}  // namespace detail

template< typename ObservationScalarType = double, typename TimeType = double >
std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > >
createSumLmkObservationSets( const std::vector< input_output::sum_lmk::SumImageData >& sumImages,
                             const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks,
                             const SumLmkObservationConversionSettings& conversionSettings,
                             const std::map< std::string, std::string >& imageIdToCameraName )
{
    std::map< ObservableType,
              std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > > > >
            observationSets;

    for( const input_output::sum_lmk::SumImageData& image : sumImages )
    {
        const TimeType observationTime = detail::convertSumUtcStringToSecondsSinceJ2000< TimeType >( image );
        const std::string cameraName = imageIdToCameraName.at( image.imageId_ );

        for( const input_output::sum_lmk::SumLandmarkObservation& observation : image.landmarkObservations_ )
        {
            if( landmarks.count( observation.landmarkId_ ) == 0 )
            {
                throw std::runtime_error( "Error when converting SUM image '" + image.imageId_ + "': missing LMK data for landmark '" +
                                          observation.landmarkId_ + "'." );
            }

            LinkEnds currentLinkEnds;
            currentLinkEnds[ transmitter ] = LinkEndId( conversionSettings.targetBodyName_, observation.landmarkId_ );
            currentLinkEnds[ receiver ] = LinkEndId( conversionSettings.receiverBodyName_, cameraName );

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentObservable( 2 );
            currentObservable( 0 ) = static_cast< ObservationScalarType >( observation.pixelCoordinates_( 0 ) );
            currentObservable( 1 ) = static_cast< ObservationScalarType >( observation.pixelCoordinates_( 1 ) );

            observationSets[ pixel_coordinates ][ currentLinkEnds ].push_back(
                    std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                            pixel_coordinates,
                            currentLinkEnds,
                            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( 1, currentObservable ),
                            std::vector< TimeType >( 1, observationTime ),
                            receiver ) );
        }
    }

    return observationSets;
}

//! Build the pixel-coordinate observation model settings matching every (image, landmark) link end
//! in a SUM/LMK observation collection. Light-time is geometric and single-leg (an empty correction
//! list still triggers the model's light-time iteration); stellar aberration is off, per the SUM/LMK
//! conventions. No bodies are needed - the camera is resolved from VehicleSystems when simulators are
//! created. The resulting settings plug directly into createObservationSimulators and the generic
//! residual/estimation machinery, exactly like any other observable type.
template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< ObservationModelSettings > > createSumLmkObservationModelSettings(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& observationCollection,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ) )
{
    if( observationCollection == nullptr )
    {
        throw std::runtime_error( "Error when creating SUM/LMK observation model settings: observation collection is null." );
    }
    std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings;
    for( const LinkDefinition& linkDefinition : observationCollection->getLinkDefinitionsForSingleObservable( pixel_coordinates ) )
    {
        observationModelSettings.push_back( pixelCoordinatesSettings( linkDefinition, lightTimeCorrections ) );
    }
    return observationModelSettings;
}

template< typename ObservationScalarType = double, typename TimeType = double >
SumLmkObservationConversionResult< ObservationScalarType, TimeType > createSumLmkObservationCollection(
        const std::vector< input_output::sum_lmk::SumImageData >& sumImages,
        const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks,
        const simulation_setup::SystemOfBodies& bodies,
        const SumLmkObservationConversionSettings& conversionSettings )
{
    detail::validateReceiverAndTargetBodies( bodies, conversionSettings );

    std::vector< input_output::sum_lmk::SumImageData > sumImagesToConvert = sumImages;
    if( conversionSettings.skipObservationsWithMissingLandmarks_ )
    {
        sumImagesToConvert = detail::filterSumImagesForAvailableLandmarks( sumImages, landmarks );
    }
    else
    {
        detail::validateReferencedLandmarksHaveDefinitions( sumImagesToConvert, landmarks );
    }
    // Do this before adding stations/cameras, so unsupported calibration input leaves the caller's
    // environment unchanged.
    detail::validateSumImageDistortionCoefficients( sumImagesToConvert );

    SumLmkObservationConversionResult< ObservationScalarType, TimeType > result;
    result.receiverBodyName_ = conversionSettings.receiverBodyName_;
    detail::validateSanitizedCameraNames( sumImagesToConvert, result.imageIdToCameraName_ );

    std::shared_ptr< simulation_setup::Body > targetBody = bodies.at( conversionSettings.targetBodyName_ );
    std::shared_ptr< simulation_setup::Body > receiverBody = bodies.at( conversionSettings.receiverBodyName_ );

    detail::addSumLmkLandmarksToBody( landmarks, targetBody );
    detail::addSumLmkCamerasToBody(
            sumImagesToConvert, landmarks, receiverBody, targetBody, conversionSettings, result.imageIdToCameraName_ );

    result.observationCollection_ = std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >(
            createSumLmkObservationSets< ObservationScalarType, TimeType >(
                    sumImagesToConvert, landmarks, conversionSettings, result.imageIdToCameraName_ ) );

    for( const input_output::sum_lmk::SumImageData& image : sumImagesToConvert )
    {
        detail::addPointingAprioriEntryIfAvailable( image,
                                                    conversionSettings.receiverBodyName_,
                                                    result.imageIdToCameraName_.at( image.imageId_ ),
                                                    result.inverseAprioriCovarianceDiagonalEntries_ );
    }

    result.observationModelSettings_ =
            createSumLmkObservationModelSettings< ObservationScalarType, TimeType >( result.observationCollection_ );

    return result;
}

template< typename ObservationScalarType = double, typename TimeType = double >
SumLmkObservationConversionResult< ObservationScalarType, TimeType > createSumLmkObservationCollection(
        const std::vector< std::string >& sumFiles,
        const std::vector< std::string >& lmkFiles,
        const simulation_setup::SystemOfBodies& bodies,
        const SumLmkObservationConversionSettings& conversionSettings )
{
    return createSumLmkObservationCollection< ObservationScalarType, TimeType >(
            input_output::sum_lmk::readSumFiles( sumFiles ), input_output::sum_lmk::readLmkFiles( lmkFiles ), bodies, conversionSettings );
}

//! Build the camera_pointing_correction parameter settings for the per-image cameras registered by a
//! SUM/LMK conversion: one 3-vector pointing parameter per image, on the receiver body that carries the
//! cameras. The settings are returned ordered by camera name, so that the resulting parameter vector has a
//! reproducible layout. Pass a non-empty imageIds to restrict the settings to a subset of the images (for
//! instance to estimate pointing only for the images that have a SIGMA_PTG a-priori); an unknown image ID
//! is an error rather than being silently ignored.
template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > createSumLmkPointingParameterSettings(
        const SumLmkObservationConversionResult< ObservationScalarType, TimeType >& conversionResult,
        const std::vector< std::string >& imageIds = std::vector< std::string >( ) )
{
    if( conversionResult.receiverBodyName_.empty( ) )
    {
        throw std::runtime_error( "Error when creating SUM/LMK pointing parameter settings: conversion result has no receiver body name." );
    }

    std::set< std::string > cameraNames;
    if( imageIds.empty( ) )
    {
        for( const auto& imageEntry : conversionResult.imageIdToCameraName_ )
        {
            cameraNames.insert( imageEntry.second );
        }
    }
    else
    {
        for( const std::string& imageId : imageIds )
        {
            if( conversionResult.imageIdToCameraName_.count( imageId ) == 0 )
            {
                throw std::runtime_error( "Error when creating SUM/LMK pointing parameter settings: image '" + imageId +
                                          "' is not part of the converted observations." );
            }
            cameraNames.insert( conversionResult.imageIdToCameraName_.at( imageId ) );
        }
    }

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterSettings;
    for( const std::string& cameraName : cameraNames )
    {
        parameterSettings.push_back( estimatable_parameters::cameraPointingCorrection( conversionResult.receiverBodyName_, cameraName ) );
    }
    return parameterSettings;
}

//! Assemble an inverse a-priori covariance matrix for an EstimationInput from per-parameter inverse-variance
//! diagonal entries, such as the per-image SIGMA_PTG entries produced by a SUM/LMK conversion. Each entry
//! identifies a parameter by its (type, body, reference point) identifier; the entry's values are written onto
//! the diagonal of that parameter's block in the full parameter vector.
//!
//! Entries whose parameter is not part of parametersToEstimate are skipped: a conversion typically produces an
//! a-priori for every image, while only a subset of the images may have their pointing estimated. Pass
//! baseInverseAprioriCovariance to add these entries on top of an existing a-priori (for instance one already
//! constraining the initial state); it must then be square with the estimated parameter-set size.
template< typename InitialStateParameterType = double >
Eigen::MatrixXd createInverseAprioriCovarianceFromDiagonalEntries(
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                inverseAprioriCovarianceDiagonalEntries,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parametersToEstimate,
        const Eigen::MatrixXd& baseInverseAprioriCovariance = Eigen::MatrixXd::Zero( 0, 0 ) )
{
    if( parametersToEstimate == nullptr )
    {
        throw std::runtime_error( "Error when creating inverse a-priori covariance: parameter set is null." );
    }

    const int numberOfParameters = parametersToEstimate->getEstimatedParameterSetSize( );
    Eigen::MatrixXd inverseAprioriCovariance;
    if( baseInverseAprioriCovariance.rows( ) == 0 && baseInverseAprioriCovariance.cols( ) == 0 )
    {
        inverseAprioriCovariance = Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters );
    }
    else if( baseInverseAprioriCovariance.rows( ) != numberOfParameters || baseInverseAprioriCovariance.cols( ) != numberOfParameters )
    {
        throw std::runtime_error( "Error when creating inverse a-priori covariance: base matrix is " +
                                  std::to_string( baseInverseAprioriCovariance.rows( ) ) + "x" +
                                  std::to_string( baseInverseAprioriCovariance.cols( ) ) + ", but the estimated parameter set has size " +
                                  std::to_string( numberOfParameters ) + "." );
    }
    else
    {
        inverseAprioriCovariance = baseInverseAprioriCovariance;
    }

    std::set< int > alreadySetStartIndices;
    for( const auto& entry : inverseAprioriCovarianceDiagonalEntries )
    {
        const std::vector< std::pair< int, int > > parameterIndices = parametersToEstimate->getIndicesForParameterType( entry.first );
        if( parameterIndices.empty( ) )
        {
            // Parameter is not estimated; its a-priori simply does not apply to this parameter set.
            continue;
        }
        if( parameterIndices.size( ) > 1 )
        {
            throw std::runtime_error( "Error when creating inverse a-priori covariance: parameter " +
                                      estimatable_parameters::getParameterTypeString( entry.first.first ) + "of (" +
                                      entry.first.second.first + ", " + entry.first.second.second +
                                      ") matches more than one parameter block." );
        }

        const int startIndex = parameterIndices.at( 0 ).first;
        const int parameterSize = parameterIndices.at( 0 ).second;
        if( entry.second.size( ) != parameterSize )
        {
            throw std::runtime_error( "Error when creating inverse a-priori covariance: a-priori for parameter " +
                                      estimatable_parameters::getParameterTypeString( entry.first.first ) + "of (" +
                                      entry.first.second.first + ", " + entry.first.second.second + ") has size " +
                                      std::to_string( entry.second.size( ) ) + ", but the parameter has size " +
                                      std::to_string( parameterSize ) + "." );
        }
        if( !alreadySetStartIndices.insert( startIndex ).second )
        {
            throw std::runtime_error( "Error when creating inverse a-priori covariance: parameter " +
                                      estimatable_parameters::getParameterTypeString( entry.first.first ) + "of (" +
                                      entry.first.second.first + ", " + entry.first.second.second +
                                      ") is given an a-priori more than once." );
        }

        for( int i = 0; i < parameterSize; ++i )
        {
            inverseAprioriCovariance( startIndex + i, startIndex + i ) += entry.second( i );
        }
    }

    return inverseAprioriCovariance;
}

//! Convenience overload assembling the inverse a-priori covariance directly from a SUM/LMK conversion result,
//! i.e. from the per-image SIGMA_PTG pointing a-priori it carries.
template< typename ObservationScalarType = double, typename TimeType = double, typename InitialStateParameterType = double >
Eigen::MatrixXd createSumLmkInverseAprioriCovariance(
        const SumLmkObservationConversionResult< ObservationScalarType, TimeType >& conversionResult,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parametersToEstimate,
        const Eigen::MatrixXd& baseInverseAprioriCovariance = Eigen::MatrixXd::Zero( 0, 0 ) )
{
    return createInverseAprioriCovarianceFromDiagonalEntries< InitialStateParameterType >(
            conversionResult.inverseAprioriCovarianceDiagonalEntries_, parametersToEstimate, baseInverseAprioriCovariance );
}

//! Compute observed-minus-computed (O-C) pixel residuals for a SUM/LMK observation collection given a
//! fixed environment (e.g. a spacecraft SPICE trajectory). Builds the matching pixel-coordinate
//! observation simulators, evaluates the model at each observation's link end/time, and stores the
//! residuals into the collection (readable via getResiduals / getResidualStatistics). Returns the
//! concatenated residual vector for convenience. This is the pixel-landmark analogue of computing
//! Doppler residuals against a reference trajectory.
template< typename ObservationScalarType = double, typename TimeType = double >
Eigen::VectorXd computeSumLmkResiduals(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& observationCollection,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ) )
{
    if( observationCollection == nullptr )
    {
        throw std::runtime_error( "Error when computing SUM/LMK residuals: observation collection is null." );
    }
    const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings =
            createSumLmkObservationModelSettings< ObservationScalarType, TimeType >( observationCollection, lightTimeCorrections );
    const std::vector< std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators =
            createObservationSimulators< ObservationScalarType, TimeType >( observationModelSettings, bodies );
    simulation_setup::computeResidualsAndDependentVariables< ObservationScalarType, TimeType >(
            observationCollection, observationSimulators, bodies );
    return observationCollection->getConcatenatedResiduals( ).template cast< double >( );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_PROCESS_SUM_LMK_FILES_H
