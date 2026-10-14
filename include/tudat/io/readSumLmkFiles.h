/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_READ_SUM_LMK_FILES_H
#define TUDAT_READ_SUM_LMK_FILES_H

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/basics/timeType.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace input_output
{

namespace sum_lmk
{

//! Tolerance for the orthonormality and determinant checks on a SUM CX/CY/CZ attitude matrix.
//!
//! The SPC archive contains valid matrices with a Frobenius orthogonality error just below 8e-8,
//! due to its stored numerical precision. This tolerance admits those matrices while still
//! rejecting matrices that are materially non-rotational.
constexpr double sumCameraRotationMatrixTolerance = 1.0E-7;

struct SumLandmarkObservation {
    std::string landmarkId_;
    Eigen::Vector2d pixelCoordinates_ = Eigen::Vector2d::Zero( );
};

struct SumLimbFitObservation {
    std::string featureId_;
    Eigen::Vector2d pixelCoordinates_ = Eigen::Vector2d::Zero( );
    double sigma_ = TUDAT_NAN;
};

struct SumImageData {
    std::string imageId_;
    std::string utcEpochString_;
    Eigen::Vector2i imageSize_ = Eigen::Vector2i::Zero( );
    int threshold_ = 0;
    int maxDn_ = 0;
    double focalLengthMm_ = TUDAT_NAN;
    Eigen::Vector2d opticalCenter_ = Eigen::Vector2d::Zero( );
    // SCOBJ: spacecraft-to-object (target-centre) vector in the target body-fixed frame (m).
    // Spacecraft position relative to the target centre is therefore -spacecraftObjectVector_.
    Eigen::Vector3d spacecraftObjectVector_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Matrix3d cameraAxes_ = Eigen::Matrix3d::Constant( TUDAT_NAN );
    Eigen::Vector3d sunDirectionBodyFixed_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Matrix< double, 2, 3 > kMatrix_ = Eigen::Matrix< double, 2, 3 >::Constant( TUDAT_NAN );
    // Raw four-value DISTORTION row from the SUM file. SPC SUM files conventionally keep this
    // zero; an actual Owen distortion model is supplied separately in INIT_LITHOS. The SUM/LMK
    // observation converter therefore rejects non-zero values rather than silently ignoring them.
    Eigen::Vector4d distortionCoefficients_ = Eigen::Vector4d::Zero( );
    Eigen::Vector3d spacecraftObjectSigma_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Vector3d pointingSigma_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    std::vector< SumLandmarkObservation > landmarkObservations_;
    std::vector< SumLimbFitObservation > limbFitObservations_;
    std::string sourceFile_;
};

struct LmkPictureObservation {
    std::string imageId_;
    Eigen::Vector2d pixelCoordinates_ = Eigen::Vector2d::Zero( );
    bool flagged_ = false;
};

struct LmkLandmarkData {
    std::string landmarkId_;
    char typeFlag_ = '\0';
    int patchSize_ = 0;
    double patchScaleKm_ = TUDAT_NAN;
    Eigen::Vector4i horizonFlags_ = Eigen::Vector4i::Zero( );
    double sigKm_ = TUDAT_NAN;
    double rmsLmk_ = TUDAT_NAN;
    Eigen::Vector3d bodyFixedPosition_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Vector3d localXAxis_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Vector3d localYAxis_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Vector3d localZAxis_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Vector3d landmarkPositionSigma_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    std::vector< LmkPictureObservation > pictures_;
    std::string sourceFile_;
};

SumImageData readSumFile( const std::string& sumFile );

std::vector< SumImageData > readSumFiles( const std::vector< std::string >& sumFiles );

LmkLandmarkData readLmkFile( const std::string& lmkFile );

std::map< std::string, LmkLandmarkData > readLmkFiles( const std::vector< std::string >& lmkFiles );

//! Name of the camera created for one SUM image: "Camera_" followed by the image ID.
std::string getSumImageCameraName( const std::string& imageId );

//! Observation epoch of a SUM image, in seconds since J2000 in the UTC time scale.
double getSumImageObservationTime( const SumImageData& image );

//! Rotation from the target body-fixed frame to the camera frame, given by the SUM CX/CY/CZ rows.
Eigen::Quaterniond getSumImageRotationFromBodyFixedToCameraFrame( const SumImageData& image );

//! Check SUM/LMK data before conversion and return the images to convert.
/*!
 *  Throws if an image has unusable camera data, or if two image IDs give the same camera name. If
 *  skipObservationsWithMissingLandmarks is true, landmarks without an LMK file are dropped from the images (and images
 *  left without landmarks are dropped). Otherwise a missing landmark is an error. If validateSpacecraftObjectGeometry
 *  is true, every landmark must lie in front of the camera when placed with the SUM SCOBJ vector.
 */
std::vector< SumImageData > getSumImagesToConvert( const std::vector< SumImageData >& sumImages,
                                                   const std::map< std::string, LmkLandmarkData >& landmarks,
                                                   const bool skipObservationsWithMissingLandmarks = false,
                                                   const bool validateSpacecraftObjectGeometry = true );

//! Create the supplementary data for SUM/LMK observations.
/*!
 *  Each landmark becomes a reference point on the target body. Each image becomes a camera on the receiver body,
 *  pointed relative to the target body-fixed frame and carrying the SUM SIGMA_PTG pointing uncertainty.
 */
std::vector< std::shared_ptr< data::TrackingSupplementaryData > > getSumLmkTrackingSupplementaryData(
        const std::vector< SumImageData >& sumImages,
        const std::map< std::string, LmkLandmarkData >& landmarks,
        const std::string& targetBodyName,
        const std::string& receiverBodyName );

inline data::PlainLinkDefinition getSumLmkLinkEnds( const std::string& targetBodyName,
                                                    const std::string& landmarkId,
                                                    const std::string& receiverBodyName,
                                                    const std::string& cameraName )
{
    return data::PlainLinkDefinition( { std::make_pair( std::make_pair( targetBodyName, landmarkId ), "transmitter" ),
                                        std::make_pair( std::make_pair( receiverBodyName, cameraName ), "receiver" ) } );
}

//! Convert parsed SUM/LMK data to pixel-coordinate tracking data and supplementary data.
/*!
 *  Every (image, landmark) pair becomes one pixel-coordinate observation, sent by the landmark on the target body and
 *  received by the image's camera on the receiver body. Epochs are kept in UTC.
 */
template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
convertSumLmkData( const std::vector< SumImageData >& sumImages,
                   const std::map< std::string, LmkLandmarkData >& landmarks,
                   const std::string& targetBodyName,
                   const std::string& receiverBodyName,
                   const bool skipObservationsWithMissingLandmarks = false,
                   const bool validateSpacecraftObjectGeometry = true )
{
    const std::vector< SumImageData > imagesToConvert =
            getSumImagesToConvert( sumImages, landmarks, skipObservationsWithMissingLandmarks, validateSpacecraftObjectGeometry );

    std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataSets;
    for( const SumImageData& image : imagesToConvert )
    {
        const TimeType observationTime = static_cast< TimeType >( getSumImageObservationTime( image ) );
        const std::string cameraName = getSumImageCameraName( image.imageId_ );
        for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
        {
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observable( 2 );
            observable( 0 ) = static_cast< ObservationScalarType >( observation.pixelCoordinates_( 0 ) );
            observable( 1 ) = static_cast< ObservationScalarType >( observation.pixelCoordinates_( 1 ) );

            trackingDataSets.push_back( std::make_shared< data::TrackingData< ObservationScalarType, TimeType > >(
                    "PixelCoordinates",
                    getSumLmkLinkEnds( targetBodyName, observation.landmarkId_, receiverBodyName, cameraName ),
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( 1, observable ),
                    std::vector< TimeType >( 1, observationTime ),
                    "receiver",
                    "UTC" ) );
        }
    }

    return std::make_pair( trackingDataSets,
                           getSumLmkTrackingSupplementaryData( imagesToConvert, landmarks, targetBodyName, receiverBodyName ) );
}

//! Read SUM and LMK files and convert them to pixel-coordinate tracking data and supplementary data.
template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
readSumLmkData( const std::vector< std::string >& sumFiles,
                const std::vector< std::string >& lmkFiles,
                const std::string& targetBodyName,
                const std::string& receiverBodyName,
                const bool skipObservationsWithMissingLandmarks = false,
                const bool validateSpacecraftObjectGeometry = true )
{
    return convertSumLmkData< ObservationScalarType, TimeType >( readSumFiles( sumFiles ),
                                                                 readLmkFiles( lmkFiles ),
                                                                 targetBodyName,
                                                                 receiverBodyName,
                                                                 skipObservationsWithMissingLandmarks,
                                                                 validateSpacecraftObjectGeometry );
}

}  // namespace sum_lmk

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READ_SUM_LMK_FILES_H
