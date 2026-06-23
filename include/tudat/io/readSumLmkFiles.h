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
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace input_output
{

namespace sum_lmk
{

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
    double tdbSecondsSinceJ2000_ = TUDAT_NAN;
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

}  // namespace sum_lmk

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READ_SUM_LMK_FILES_H
