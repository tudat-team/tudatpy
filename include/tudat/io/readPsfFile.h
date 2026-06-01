/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: 820-013, TRK-2-18 Tracking System Interfaces Orbit Data File
 * Interface, Revision E, 2008, JPL/DSN
 *
 */

#ifndef TUDAT_READ_PSF_FILE_H
#define TUDAT_READ_PSF_FILE_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cstdint>

#include <Eigen/Core>

namespace tudat
{

namespace input_output
{

namespace psf
{

enum class OpticalImageType { star, planet, satellite, rock, end_marker, unknown };

class CameraModel
{
public:
    CameraModel( ) = default;

    std::string cameraId_;
    double focalLengthMm_ = 0.0;

    Eigen::Vector2d principalPoint_ = Eigen::Vector2d::Zero( );     // [pixel; line]
    Eigen::Vector4d fieldOfViewBounds_ = Eigen::Vector4d::Zero( );  // [min_p; max_p; min_l; max_l]

    // Per PSF spec: KMAT is 2x3 per camera
    Eigen::Matrix< double, 2, 3 > kMatrix_ = Eigen::Matrix< double, 2, 3 >::Zero( );

    Eigen::Matrix< double, 6, 1 > distortionCoefficients_ = Eigen::Matrix< double, 6, 1 >::Zero( );
    Eigen::Vector3d mountingOffsetsDegrees_ = Eigen::Vector3d::Zero( );  // [elev; cross-elev; twist]
};

class RawPsfMeasurement
{
public:
    virtual ~RawPsfMeasurement( ) {}

    OpticalImageType opticalImageType_ = OpticalImageType::unknown;

    std::string imageName_;     // raw IMG (string, without quotes)
    std::int64_t imageId_ = 0;  // raw IMGID
    int useFlag_ = 0;           // raw USE

    Eigen::Vector2d observedPixelLine_ = Eigen::Vector2d::Zero( );  // Z
    Eigen::Vector2d localCorrection_ = Eigen::Vector2d::Zero( );    // ZC
    Eigen::Vector2d sigmaPixelLine_ = Eigen::Vector2d::Zero( );     // SIG

    Eigen::Vector2d getEffectivePixelLine( ) const;
};

class RawPsfStarMeasurement : public RawPsfMeasurement
{
public:
    double rightAscensionDegrees_ = 0.0;  // STRA
    double declinationDegrees_ = 0.0;     // STDEC
};

class RawPsfFileImageContents
{
public:
    RawPsfFileImageContents( ) = default;

    std::string pictureName_;                 // PICNM
    int pictureNumber_ = 0;                   // PICNO
    std::string endOfExposureTimeUtcString_;  // TOB
    std::string cameraId_;                    // CAMERA

    double exposureTimeSeconds_ = 0.0;  // EXPTIM
    int pictureDeletionFlag_ = 0;       // PICDEL

    double rightAscensionDegrees_ = 0.0;  // RA
    double declinationDegrees_ = 0.0;     // DEC
    double twistDegrees_ = 0.0;           // TWIST

    std::vector< std::shared_ptr< RawPsfMeasurement > > measurements_;
};

class RawPsfFileContents
{
public:
    RawPsfFileContents( ) = default;

    // $ID
    std::string spacecraftId_;  // SCID
    int numberOfCameras_ = 0;   // NCAM
    int equinox_ = 0;           // EQUNOX

    std::string psfId_;                       // PSFID
    std::string psfProgram_;                  // PSFPRG
    std::string psfGenerationTimeUtcString_;  // PSFTIM
    std::vector< std::string > psfComments_;  // PSFCOM list

    // $CAM
    std::map< std::string, CameraModel > cameraModels_;  // keyed by CAMID ('A','B',...)

    // $PIC/$IM
    std::vector< RawPsfFileImageContents > images_;
};

RawPsfFileContents readPsfFile( const std::string& psfFile );

}  // namespace psf

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READ_PSF_FILE_H
