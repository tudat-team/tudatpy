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
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/io/readSumLmkFiles.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/system_models/camera.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"
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
};

template< typename ObservationScalarType = double, typename TimeType = double >
struct SumLmkObservationConversionResult {
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection_;
    std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
            inverseAprioriCovarianceDiagonalEntries_;
    std::map< std::string, std::string > imageIdToCameraName_;
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

inline Eigen::Quaterniond getQuaternionFromRotationVector( const Eigen::Vector3d& rotationVector )
{
    const double angle = rotationVector.norm( );
    if( angle <= std::numeric_limits< double >::epsilon( ) )
    {
        return Eigen::Quaterniond::Identity( );
    }
    return Eigen::Quaterniond( Eigen::AngleAxisd( angle, rotationVector / angle ) );
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

inline std::shared_ptr< system_models::PsfCameraProjectionModel > createSumCameraProjectionModel(
        const input_output::sum_lmk::SumImageData& image )
{
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
        const std::shared_ptr< Eigen::Vector3d > pointingCorrection = std::make_shared< Eigen::Vector3d >( Eigen::Vector3d::Zero( ) );
        const std::shared_ptr< ephemerides::RotationalEphemeris > targetRotationalEphemeris = targetBody->getRotationalEphemeris( );
        std::function< Eigen::Quaterniond( const double ) > rotationFromInertialToCameraFrameFunction =
                [ rotationFromTargetBodyFixedToCamera, pointingCorrection, targetRotationalEphemeris ]( const double time ) {
                    const Eigen::Quaterniond correctionRotation = getQuaternionFromRotationVector( *pointingCorrection );
                    const Eigen::Quaterniond nominalRotation(
                            rotationFromTargetBodyFixedToCamera *
                            targetRotationalEphemeris->getRotationToTargetFrame( time ).toRotationMatrix( ) );
                    return Eigen::Quaterniond( correctionRotation * nominalRotation ).normalized( );
                };

        std::shared_ptr< system_models::Camera > camera =
                std::make_shared< system_models::Camera >( cameraName,
                                                           Eigen::Quaterniond( rotationFromTargetBodyFixedToCamera ),
                                                           projectionModel,
                                                           rotationFromInertialToCameraFrameFunction,
                                                           pointingCorrection );
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

    SumLmkObservationConversionResult< ObservationScalarType, TimeType > result;
    detail::validateSanitizedCameraNames( sumImages, result.imageIdToCameraName_ );

    std::shared_ptr< simulation_setup::Body > targetBody = bodies.at( conversionSettings.targetBodyName_ );
    std::shared_ptr< simulation_setup::Body > receiverBody = bodies.at( conversionSettings.receiverBodyName_ );

    detail::addSumLmkLandmarksToBody( landmarks, targetBody );
    detail::addSumLmkCamerasToBody( sumImages, landmarks, receiverBody, targetBody, conversionSettings, result.imageIdToCameraName_ );

    result.observationCollection_ = std::make_shared< ObservationCollection< ObservationScalarType, TimeType > >(
            createSumLmkObservationSets< ObservationScalarType, TimeType >(
                    sumImages, landmarks, conversionSettings, result.imageIdToCameraName_ ) );

    for( const input_output::sum_lmk::SumImageData& image : sumImages )
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
