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

//! Pixel coordinates of one landmark, as measured in one SUM image.
struct SumLandmarkObservation {
    std::string landmarkId_;
    Eigen::Vector2d pixelCoordinates_ = Eigen::Vector2d::Zero( );
};

//! Limb-fit measurement in one SUM image (read, but not used for observations).
struct SumLimbFitObservation {
    std::string featureId_;
    Eigen::Vector2d pixelCoordinates_ = Eigen::Vector2d::Zero( );
    double sigma_ = TUDAT_NAN;
};

//! Raw contents of one SUM file, which describes one image. Distances are converted from km to m.
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
    // Raw four-value DISTORTION row; SPC keeps this zero, and conversion rejects non-zero values.
    Eigen::Vector4d distortionCoefficients_ = Eigen::Vector4d::Zero( );
    Eigen::Vector3d spacecraftObjectSigma_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Vector3d pointingSigma_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    std::vector< SumLandmarkObservation > landmarkObservations_;
    std::vector< SumLimbFitObservation > limbFitObservations_;
    std::string sourceFile_;
};

//! Image in which a landmark is seen, as listed in the PICTURES section of an LMK file.
struct LmkPictureObservation {
    std::string imageId_;
    Eigen::Vector2d pixelCoordinates_ = Eigen::Vector2d::Zero( );
    bool flagged_ = false;
};

//! Raw contents of one LMK file, which describes one landmark. Distances are converted from km to m.
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

//! Read a single SUM file.
/*!
 *  Read a single SUM file. Only the presence of the records required for conversion is checked here; the values are checked
 *  when converting the data (see getSumImagesToConvert).
 *  \param sumFile Path to the SUM file.
 *  \return Raw contents of the file.
 */
SumImageData readSumFile( const std::string& sumFile );

//! Read a list of SUM files.
/*!
 *  Read a list of SUM files.
 *  \param sumFiles Paths to the SUM files.
 *  \return Raw contents of each file, in the order of the input.
 */
std::vector< SumImageData > readSumFiles( const std::vector< std::string >& sumFiles );

//! Read a single LMK file.
/*!
 *  Read a single LMK file. Only the presence of the landmark ID and VLM record is checked.
 *  \param lmkFile Path to the LMK file.
 *  \return Raw contents of the file.
 */
LmkLandmarkData readLmkFile( const std::string& lmkFile );

//! Read a list of LMK files.
/*!
 *  Read a list of LMK files. Each landmark may be provided only once.
 *  \param lmkFiles Paths to the LMK files.
 *  \return Raw contents of each file, with the landmark ID as key.
 */
std::map< std::string, LmkLandmarkData > readLmkFiles( const std::vector< std::string >& lmkFiles );

//! Get the observation epoch of a SUM image, in seconds since J2000 in the UTC time scale.
double getSumImageObservationTime( const SumImageData& image );

//! Get the rotation from the target body-fixed frame to the camera frame, given by the SUM CX/CY/CZ rows.
Eigen::Quaterniond getSumImageRotationFromBodyFixedToCameraFrame( const SumImageData& image );

//! Get the camera calibration of a SUM image: MMFL, CTR, K-MATRIX and the image size as field-of-view bounds.
data::CameraPictureCalibration getSumImageCalibration( const SumImageData& image );

//! Check SUM/LMK data before conversion, and return the images to convert.
/*!
 *  Check SUM/LMK data before conversion, and return the images to convert, sorted by epoch. Each image may be provided only once,
 *  and its calibration, attitude, pointing uncertainty and pixel coordinates must be valid. All images are taken by one camera,
 *  so they must have distinct epochs. Images without landmark observations are not converted.
 *  \param sumImages Raw SUM image contents.
 *  \param landmarks Raw LMK landmark contents, with the landmark ID as key.
 *  \param skipObservationsWithMissingLandmarks If true, observations of landmarks without LMK data are skipped; if false,
 *  they are an error.
 *  \param validateSpacecraftObjectGeometry If true, every landmark must lie in front of the camera when placed with the
 *  SUM SCOBJ vector.
 *  \return Images to convert, containing only the observations of landmarks with LMK data.
 */
std::vector< SumImageData > getSumImagesToConvert( const std::vector< SumImageData >& sumImages,
                                                   const std::map< std::string, LmkLandmarkData >& landmarks,
                                                   const bool skipObservationsWithMissingLandmarks = false,
                                                   const bool validateSpacecraftObjectGeometry = true );

//! Create the supplementary data for SUM/LMK observations.
/*!
 *  Create the supplementary data for SUM/LMK observations. Each landmark becomes a reference point on the target body. The
 *  images become one camera on the receiver body, with the ID and calibration of each image, its pointing relative to the
 *  target body-fixed frame and its SUM SIGMA_PTG pointing uncertainty.
 *  \param sumImages SUM images to convert (output of getSumImagesToConvert).
 *  \param landmarks Raw LMK landmark contents, with the landmark ID as key.
 *  \param targetBodyName Name of the body on which the landmarks are located.
 *  \param receiverBodyName Name of the body carrying the camera.
 *  \param cameraName Name of the camera that took the images.
 *  \return Supplementary data for the landmarks and the camera.
 */
std::vector< std::shared_ptr< data::TrackingSupplementaryData > > getSumLmkTrackingSupplementaryData(
        const std::vector< SumImageData >& sumImages,
        const std::map< std::string, LmkLandmarkData >& landmarks,
        const std::string& targetBodyName,
        const std::string& receiverBodyName,
        const std::string& cameraName );

//! Get the link ends of a SUM/LMK observation: transmitted by a landmark on the target, received by the camera.
inline data::PlainLinkDefinition getSumLmkLinkEnds( const std::string& targetBodyName,
                                                    const std::string& landmarkId,
                                                    const std::string& receiverBodyName,
                                                    const std::string& cameraName )
{
    return data::PlainLinkDefinition( { std::make_pair( std::make_pair( targetBodyName, landmarkId ), "transmitter" ),
                                        std::make_pair( std::make_pair( receiverBodyName, cameraName ), "receiver" ) } );
}

//! Convert raw SUM/LMK data to pixel-coordinate tracking data and supplementary data.
/*!
 *  Convert raw SUM/LMK data to pixel-coordinate tracking data and supplementary data. Every (image, landmark) pair becomes one
 *  pixel-coordinate observation, sent by the landmark on the target body and received by the camera on the receiver body. The
 *  observations of each landmark are collected in one tracking data object. Epochs are kept in UTC.
 *  \param sumImages Raw SUM image contents.
 *  \param landmarks Raw LMK landmark contents, with the landmark ID as key.
 *  \param targetBodyName Name of the body on which the landmarks are located.
 *  \param receiverBodyName Name of the body carrying the camera.
 *  \param cameraName Name of the camera that took the images.
 *  \param skipObservationsWithMissingLandmarks If true, observations of landmarks without LMK data are skipped; if false,
 *  they are an error.
 *  \param validateSpacecraftObjectGeometry If true, every landmark must lie in front of the camera when placed with the
 *  SUM SCOBJ vector.
 *  \return Pair of tracking data (one object per landmark) and supplementary data.
 */
template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
convertSumLmkData( const std::vector< SumImageData >& sumImages,
                   const std::map< std::string, LmkLandmarkData >& landmarks,
                   const std::string& targetBodyName,
                   const std::string& receiverBodyName,
                   const std::string& cameraName,
                   const bool skipObservationsWithMissingLandmarks = false,
                   const bool validateSpacecraftObjectGeometry = true )
{
    const std::vector< SumImageData > imagesToConvert =
            getSumImagesToConvert( sumImages, landmarks, skipObservationsWithMissingLandmarks, validateSpacecraftObjectGeometry );

    // Images are sorted by epoch, so the observations of each landmark are as well.
    std::map< std::string, std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > observablesPerLandmark;
    std::map< std::string, std::vector< TimeType > > timesPerLandmark;
    for( const SumImageData& image : imagesToConvert )
    {
        const TimeType observationTime = static_cast< TimeType >( getSumImageObservationTime( image ) );
        for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
        {
            observablesPerLandmark[ observation.landmarkId_ ].push_back( observation.pixelCoordinates_.cast< ObservationScalarType >( ) );
            timesPerLandmark[ observation.landmarkId_ ].push_back( observationTime );
        }
    }

    std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataSets;
    for( const auto& landmarkEntry : observablesPerLandmark )
    {
        trackingDataSets.push_back( std::make_shared< data::TrackingData< ObservationScalarType, TimeType > >(
                "PixelCoordinates",
                getSumLmkLinkEnds( targetBodyName, landmarkEntry.first, receiverBodyName, cameraName ),
                landmarkEntry.second,
                timesPerLandmark.at( landmarkEntry.first ),
                "receiver",
                "UTC" ) );
    }

    return std::make_pair( trackingDataSets,
                           getSumLmkTrackingSupplementaryData( imagesToConvert, landmarks, targetBodyName, receiverBodyName, cameraName ) );
}

//! Read SUM and LMK files and convert them to pixel-coordinate tracking data and supplementary data.
/*!
 *  Read SUM and LMK files and convert them to pixel-coordinate tracking data and supplementary data, see convertSumLmkData.
 *  \param sumFiles Paths to the SUM files.
 *  \param lmkFiles Paths to the LMK files.
 *  \param targetBodyName Name of the body on which the landmarks are located.
 *  \param receiverBodyName Name of the body carrying the camera.
 *  \param cameraName Name of the camera that took the images.
 *  \param skipObservationsWithMissingLandmarks If true, observations of landmarks without LMK data are skipped; if false,
 *  they are an error.
 *  \param validateSpacecraftObjectGeometry If true, every landmark must lie in front of the camera when placed with the
 *  SUM SCOBJ vector.
 *  \return Pair of tracking data (one object per landmark) and supplementary data.
 */
template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
readSumLmkData( const std::vector< std::string >& sumFiles,
                const std::vector< std::string >& lmkFiles,
                const std::string& targetBodyName,
                const std::string& receiverBodyName,
                const std::string& cameraName,
                const bool skipObservationsWithMissingLandmarks = false,
                const bool validateSpacecraftObjectGeometry = true )
{
    return convertSumLmkData< ObservationScalarType, TimeType >( readSumFiles( sumFiles ),
                                                                 readLmkFiles( lmkFiles ),
                                                                 targetBodyName,
                                                                 receiverBodyName,
                                                                 cameraName,
                                                                 skipObservationsWithMissingLandmarks,
                                                                 validateSpacecraftObjectGeometry );
}

}  // namespace sum_lmk

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READ_SUM_LMK_FILES_H
