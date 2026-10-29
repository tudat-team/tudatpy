/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_CAMERA_H
#define TUDAT_CAMERA_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/nearestNeighbourSearch.h"
#include "tudat/math/basic/rotationRepresentations.h"

namespace tudat
{

namespace system_models
{

//! Base class for camera projection models.
/*!
 *  This class defines the mapping between a direction in the camera frame and image-plane pixel/line coordinates.
 *  The derivative is with respect to the Cartesian position vector in the camera frame.
 */
class CameraProjectionModel
{
public:
    virtual ~CameraProjectionModel( ) {}

    virtual Eigen::Vector2d projectUnitVectorToPixelLine( const Eigen::Vector3d& cameraFrameUnitVector ) const = 0;

    virtual Eigen::Vector3d pixelLineToUnitVector( const Eigen::Vector2d& pixelLine ) const = 0;

    virtual Eigen::Matrix< double, 2, 3 > getPixelLinePartialWrtCameraFramePosition( const Eigen::Vector3d& cameraFramePosition ) const = 0;

    virtual Eigen::DiagonalMatrix< double, 2 > getFocalLengthsMatrix( ) const
    {
        throw std::runtime_error( "Error when requesting diagonal focal lengths from a non-pinhole camera projection model." );
    }

    virtual Eigen::Vector2d getOpticalCenter( ) const
    {
        throw std::runtime_error( "Error when requesting optical center from a camera projection model that does not define one." );
    }
};

//! Pinhole camera projection model.
/*!
 *  Implements the existing Tudat camera projection:
 *      pixelLine = opticalCenter + diag(fx, fy) * [ X/Z, Y/Z ].
 */
class PinholeCameraProjectionModel : public CameraProjectionModel
{
public:
    PinholeCameraProjectionModel( const std::pair< double, double >& focalLengths = std::make_pair( 1.0, 1.0 ),
                                  const std::pair< double, double >& opticalCenter = std::make_pair( 0.0, 0.0 ) ):
        focalLengthsMatrix_( focalLengths.first, focalLengths.second ), opticalCenter_( opticalCenter.first, opticalCenter.second )
    {}

    PinholeCameraProjectionModel( const Eigen::DiagonalMatrix< double, 2 >& focalLengthsMatrix, const Eigen::Vector2d& opticalCenter ):
        focalLengthsMatrix_( focalLengthsMatrix ), opticalCenter_( opticalCenter )
    {}

    Eigen::Vector2d projectUnitVectorToPixelLine( const Eigen::Vector3d& cameraFrameUnitVector ) const override
    {
        validateProjectionDirection( cameraFrameUnitVector );
        return opticalCenter_ +
                focalLengthsMatrix_ *
                Eigen::Vector2d( cameraFrameUnitVector.x( ) / cameraFrameUnitVector.z( ),
                                 cameraFrameUnitVector.y( ) / cameraFrameUnitVector.z( ) );
    }

    Eigen::Vector3d pixelLineToUnitVector( const Eigen::Vector2d& pixelLine ) const override
    {
        Eigen::Vector2d focalPlaneCoordinates = focalLengthsMatrix_.inverse( ) * ( pixelLine - opticalCenter_ );
        return Eigen::Vector3d( focalPlaneCoordinates.x( ), focalPlaneCoordinates.y( ), 1.0 ).normalized( );
    }

    Eigen::Matrix< double, 2, 3 > getPixelLinePartialWrtCameraFramePosition( const Eigen::Vector3d& cameraFramePosition ) const override
    {
        validateProjectionDirection( cameraFramePosition );
        Eigen::Matrix< double, 2, 3 > partial = Eigen::Matrix< double, 2, 3 >::Zero( );
        partial( 0, 0 ) = focalLengthsMatrix_.diagonal( )( 0 ) / cameraFramePosition.z( );
        partial( 0, 2 ) = -focalLengthsMatrix_.diagonal( )( 0 ) * cameraFramePosition.x( ) / std::pow( cameraFramePosition.z( ), 2 );
        partial( 1, 1 ) = focalLengthsMatrix_.diagonal( )( 1 ) / cameraFramePosition.z( );
        partial( 1, 2 ) = -focalLengthsMatrix_.diagonal( )( 1 ) * cameraFramePosition.y( ) / std::pow( cameraFramePosition.z( ), 2 );
        return partial;
    }

    Eigen::DiagonalMatrix< double, 2 > getFocalLengthsMatrix( ) const override
    {
        return focalLengthsMatrix_;
    }

    Eigen::Vector2d getOpticalCenter( ) const override
    {
        return opticalCenter_;
    }

private:
    void validateProjectionDirection( const Eigen::Vector3d& cameraFrameVector ) const
    {
        if( cameraFrameVector.z( ) <= std::numeric_limits< double >::epsilon( ) )
        {
            throw std::runtime_error(
                    "Error when projecting camera-frame vector: observed point is behind the camera or too close to the focal plane." );
        }
    }

    Eigen::DiagonalMatrix< double, 2 > focalLengthsMatrix_;
    Eigen::Vector2d opticalCenter_;
};

//! PSF/Jacobson camera projection model.
/*!
 *  Implements the projection used by Jacobson's PSF files:
 *      [p, l] = KMAT * [x, y, x*y] + PLCTR,
 *      x = f * X / Z, y = f * Y / Z.
 *
 *  Distortion coefficients, field-of-view bounds, and mounting offsets are stored so that PSF camera properties are kept together.
 *  The present projection implements the KMAT/PLCTR/FL mapping; distortion and mounting-offset handling can be added without changing
 *  the Camera interface.
 */
class PsfCameraProjectionModel : public CameraProjectionModel
{
public:
    PsfCameraProjectionModel( const double focalLength,
                              const Eigen::Vector2d& principalPoint,
                              const Eigen::Matrix< double, 2, 3 >& kMatrix,
                              const Eigen::Matrix< double, 6, 1 >& distortionCoefficients = Eigen::Matrix< double, 6, 1 >::Zero( ),
                              const Eigen::Vector3d& mountingOffsetsDegrees = Eigen::Vector3d::Zero( ),
                              const Eigen::Vector4d& fieldOfViewBounds = Eigen::Vector4d::Zero( ) ):
        focalLength_( focalLength ), principalPoint_( principalPoint ), kMatrix_( kMatrix ),
        distortionCoefficients_( distortionCoefficients ), mountingOffsetsDegrees_( mountingOffsetsDegrees ),
        fieldOfViewBounds_( fieldOfViewBounds )
    {
        validateCalibration( );
    }

    Eigen::Vector2d projectUnitVectorToPixelLine( const Eigen::Vector3d& cameraFrameUnitVector ) const override
    {
        validateProjectionDirection( cameraFrameUnitVector );
        const double x = focalLength_ * cameraFrameUnitVector.x( ) / cameraFrameUnitVector.z( );
        const double y = focalLength_ * cameraFrameUnitVector.y( ) / cameraFrameUnitVector.z( );
        return focalPlaneCoordinatesToPixelLine( x, y );
    }

    Eigen::Vector3d pixelLineToUnitVector( const Eigen::Vector2d& pixelLine ) const override
    {
        Eigen::Vector2d focalPlaneCoordinates = pixelLineToFocalPlaneCoordinates( pixelLine );
        return Eigen::Vector3d( focalPlaneCoordinates.x( ) / focalLength_, focalPlaneCoordinates.y( ) / focalLength_, 1.0 ).normalized( );
    }

    Eigen::Matrix< double, 2, 3 > getPixelLinePartialWrtCameraFramePosition( const Eigen::Vector3d& cameraFramePosition ) const override
    {
        validateProjectionDirection( cameraFramePosition );
        const double x = focalLength_ * cameraFramePosition.x( ) / cameraFramePosition.z( );
        const double y = focalLength_ * cameraFramePosition.y( ) / cameraFramePosition.z( );

        Eigen::Matrix2d partialPixelLineWrtFocalPlaneCoordinates;
        partialPixelLineWrtFocalPlaneCoordinates( 0, 0 ) = kMatrix_( 0, 0 ) + kMatrix_( 0, 2 ) * y;
        partialPixelLineWrtFocalPlaneCoordinates( 0, 1 ) = kMatrix_( 0, 1 ) + kMatrix_( 0, 2 ) * x;
        partialPixelLineWrtFocalPlaneCoordinates( 1, 0 ) = kMatrix_( 1, 0 ) + kMatrix_( 1, 2 ) * y;
        partialPixelLineWrtFocalPlaneCoordinates( 1, 1 ) = kMatrix_( 1, 1 ) + kMatrix_( 1, 2 ) * x;

        Eigen::Matrix< double, 2, 3 > partialFocalPlaneCoordinatesWrtCameraFramePosition = Eigen::Matrix< double, 2, 3 >::Zero( );
        partialFocalPlaneCoordinatesWrtCameraFramePosition( 0, 0 ) = focalLength_ / cameraFramePosition.z( );
        partialFocalPlaneCoordinatesWrtCameraFramePosition( 0, 2 ) =
                -focalLength_ * cameraFramePosition.x( ) / std::pow( cameraFramePosition.z( ), 2 );
        partialFocalPlaneCoordinatesWrtCameraFramePosition( 1, 1 ) = focalLength_ / cameraFramePosition.z( );
        partialFocalPlaneCoordinatesWrtCameraFramePosition( 1, 2 ) =
                -focalLength_ * cameraFramePosition.y( ) / std::pow( cameraFramePosition.z( ), 2 );

        return partialPixelLineWrtFocalPlaneCoordinates * partialFocalPlaneCoordinatesWrtCameraFramePosition;
    }

    Eigen::Vector2d getOpticalCenter( ) const override
    {
        return principalPoint_;
    }

    double getFocalLength( ) const
    {
        return focalLength_;
    }

    Eigen::Matrix< double, 2, 3 > getKMatrix( ) const
    {
        return kMatrix_;
    }

    Eigen::Matrix< double, 6, 1 > getDistortionCoefficients( ) const
    {
        return distortionCoefficients_;
    }

    Eigen::Vector3d getMountingOffsetsDegrees( ) const
    {
        return mountingOffsetsDegrees_;
    }

    Eigen::Vector4d getFieldOfViewBounds( ) const
    {
        return fieldOfViewBounds_;
    }

private:
    void validateCalibration( ) const
    {
        if( !( focalLength_ > 0.0 ) || !principalPoint_.allFinite( ) || !kMatrix_.allFinite( ) )
        {
            throw std::runtime_error(
                    "Error when creating PSF camera model: focal length must be positive, and principal point and "
                    "K-MATRIX must be finite." );
        }
    }

    void validateProjectionDirection( const Eigen::Vector3d& cameraFrameVector ) const
    {
        if( cameraFrameVector.z( ) <= std::numeric_limits< double >::epsilon( ) )
        {
            throw std::runtime_error(
                    "Error when projecting PSF camera-frame vector: observed point is behind the camera or too close to the focal plane." );
        }
    }

    Eigen::Vector2d focalPlaneCoordinatesToPixelLine( const double x, const double y ) const
    {
        return principalPoint_ + kMatrix_ * Eigen::Vector3d( x, y, x * y );
    }

    Eigen::Vector2d pixelLineToFocalPlaneCoordinates( const Eigen::Vector2d& pixelLine ) const
    {
        const Eigen::Vector2d target = pixelLine - principalPoint_;
        Eigen::Vector2d focalPlaneCoordinates = kMatrix_.block< 2, 2 >( 0, 0 ).fullPivLu( ).solve( target );

        for( unsigned int i = 0; i < 20; ++i )
        {
            const double x = focalPlaneCoordinates.x( );
            const double y = focalPlaneCoordinates.y( );
            const Eigen::Vector2d residual = focalPlaneCoordinatesToPixelLine( x, y ) - pixelLine;
            if( residual.norm( ) < 1.0E-13 )
            {
                return focalPlaneCoordinates;
            }

            Eigen::Matrix2d jacobian;
            jacobian( 0, 0 ) = kMatrix_( 0, 0 ) + kMatrix_( 0, 2 ) * y;
            jacobian( 0, 1 ) = kMatrix_( 0, 1 ) + kMatrix_( 0, 2 ) * x;
            jacobian( 1, 0 ) = kMatrix_( 1, 0 ) + kMatrix_( 1, 2 ) * y;
            jacobian( 1, 1 ) = kMatrix_( 1, 1 ) + kMatrix_( 1, 2 ) * x;

            focalPlaneCoordinates -= jacobian.fullPivLu( ).solve( residual );
        }

        throw std::runtime_error( "Error when converting PSF pixel/line coordinates to camera-frame unit vector: Newton solve failed." );
    }

    double focalLength_;
    Eigen::Vector2d principalPoint_;
    Eigen::Matrix< double, 2, 3 > kMatrix_;
    Eigen::Matrix< double, 6, 1 > distortionCoefficients_;
    Eigen::Vector3d mountingOffsetsDegrees_;
    Eigen::Vector4d fieldOfViewBounds_;
};

class Camera
{
public:
    //! Constructor of the Camera class.
    /*! \brief sets the camera name, rotation from body-fixed to camera frame, and focal lengths in x and y directions.
     *  Camera is fixed on the body and points along the z-axis of its own frame.
     *  \param cameraId Name of the camera.
     *  \param rotationFromBodyFixedToCameraFrame Rotation from body-fixed frame to camera frame.
     *  \param focal_lengths Focal lengths in x and y directions [px].
     *  \param optical_center Optical center of the camera [px].
     */
    Camera( const std::string& cameraId,
            const Eigen::Quaterniond& rotationFromBodyFixedToCameraFrame,
            std::pair< double, double > focal_lengths = std::make_pair( 1.0, 1.0 ),
            std::pair< double, double > optical_center = std::make_pair( 0.0, 0.0 ),
            std::function< Eigen::Quaterniond( const double ) > rotationFromInertialToCameraFrameFunction = nullptr ):
        cameraId_( cameraId ), rotationFromBodyFixedToCameraFrame_( rotationFromBodyFixedToCameraFrame ),
        projectionModels_( { std::make_shared< PinholeCameraProjectionModel >( focal_lengths, optical_center ) } ),
        rotationFromInertialToCameraFrameFunction_( rotationFromInertialToCameraFrameFunction )
    {}

    Camera( const std::string& cameraId,
            const Eigen::Quaterniond& rotationFromBodyFixedToCameraFrame,
            const std::shared_ptr< CameraProjectionModel > projectionModel,
            std::function< Eigen::Quaterniond( const double ) > rotationFromInertialToCameraFrameFunction = nullptr ):
        cameraId_( cameraId ), rotationFromBodyFixedToCameraFrame_( rotationFromBodyFixedToCameraFrame ),
        projectionModels_( { projectionModel } ), rotationFromInertialToCameraFrameFunction_( rotationFromInertialToCameraFrameFunction )
    {
        if( projectionModel == nullptr )
        {
            throw std::runtime_error( "Error when creating camera: projection model is nullptr." );
        }
    }

    /*! \brief Get the camera ID.
     *  \return The camera ID.
     */
    std::string getCameraId( ) const
    {
        return cameraId_;
    }

    Eigen::DiagonalMatrix< double, 2 > getFocalLengthsMatrix( const double secondsSinceEpoch = TUDAT_NAN ) const
    {
        return getProjectionModel( secondsSinceEpoch )->getFocalLengthsMatrix( );
    }

    Eigen::Vector2d getOpticalCenter( const double secondsSinceEpoch = TUDAT_NAN ) const
    {
        return getProjectionModel( secondsSinceEpoch )->getOpticalCenter( );
    }

    //! Function to get the projection model at the given time: the constant one, or that of the picture closest in time. A time is
    //! only required when the calibration changes between pictures.
    std::shared_ptr< CameraProjectionModel > getProjectionModel( const double secondsSinceEpoch = TUDAT_NAN ) const
    {
        return projectionModels_.size( ) == 1 ? projectionModels_.at( 0 ) : projectionModels_.at( getPictureIndex( secondsSinceEpoch ) );
    }

    bool hasRotationFromInertialToCameraFrameFunction( ) const
    {
        return static_cast< bool >( rotationFromInertialToCameraFrameFunction_ );
    }

    //! Body whose rotation is composed into the direct pointing function; empty for inertial pointing.
    const std::string& getPointingReferenceBodyName( ) const
    {
        return pointingReferenceBodyName_;
    }

    void setPointingReferenceBodyName( const std::string& pointingReferenceBodyName )
    {
        pointingReferenceBodyName_ = pointingReferenceBodyName;
    }

    /*! \brief Get quaternion representing active rotation from body-fixed to camera frame.
     *  \return The rotation from body-fixed to camera frame.
     */
    Eigen::Quaterniond getRotationFromBodyFixedToCameraFrame( ) const
    {
        return rotationFromBodyFixedToCameraFrame_;
    }

    //! Function to get the rotation from inertial to camera frame of a camera fixed to its host body, including the pointing
    //! correction of the picture closest to the given time.
    /*!
     *  \param rotationFromInertialToBodyFixedFrame Rotation from inertial frame to the body-fixed frame of the camera host body.
     *  \param secondsSinceEpoch Time at which the rotation is evaluated.
     *  \return The rotation from inertial to camera frame.
     */
    Eigen::Quaterniond getRotationFromInertialToCameraFrame( const Eigen::Quaterniond& rotationFromInertialToBodyFixedFrame,
                                                             const double secondsSinceEpoch ) const
    {
        return applyPointingCorrection_( rotationFromBodyFixedToCameraFrame_ * rotationFromInertialToBodyFixedFrame, secondsSinceEpoch );
    }

    Eigen::Quaterniond getRotationFromInertialToCameraFrame( const double secondsSinceEpoch ) const
    {
        if( !rotationFromInertialToCameraFrameFunction_ )
        {
            throw std::runtime_error(
                    "Error when calculating rotation from inertial to camera frame: no direct inertial-to-camera rotation function was "
                    "provided for the camera." );
        }
        if( std::isnan( secondsSinceEpoch ) )
        {
            throw std::runtime_error( "Error when calculating rotation from inertial to camera frame of camera " + cameraId_ +
                                      ": a time is required, since the camera pointing is time-dependent." );
        }
        return applyPointingCorrection_( rotationFromInertialToCameraFrameFunction_( secondsSinceEpoch ), secondsSinceEpoch );
    }

    /*! \brief Calculate the pixel coordinates of a point, for a camera fixed to its host body.
     *  The position is rotated to the camera frame with the camera orientation w.r.t. its host body and the pointing correction at
     *  the given time, and projected with the calibration at that time. A time is only required for a camera with several
     *  pictures. Not available for a camera whose orientation is given directly w.r.t. the inertial frame.
     * \param positionOfObservedBodyInBodyFrame Position of the observed point relative to the camera, in the body-fixed frame of
     * the camera host body.
     * \param secondsSinceEpoch Time at which the pixel coordinates are computed.
     * \return The pixel coordinates of the observed point.
     */
    Eigen::Vector2d calculateObservableFromBodyFixed( const Eigen::Vector3d& positionOfObservedBodyInBodyFrame,
                                                      const double secondsSinceEpoch = TUDAT_NAN ) const
    {
        if( hasRotationFromInertialToCameraFrameFunction( ) )
        {
            throw std::runtime_error( "Error when calculating pixel coordinates of camera " + cameraId_ +
                                      " from a body-fixed position: the camera orientation is given w.r.t. the inertial frame." );
        }
        return calculateObservableFromInertial( positionOfObservedBodyInBodyFrame, Eigen::Quaterniond::Identity( ), secondsSinceEpoch );
    }

    Eigen::Vector3d pixelLineToCameraFrameUnitVector( const Eigen::Vector2d& pixelLine, const double secondsSinceEpoch = TUDAT_NAN ) const
    {
        return getProjectionModel( secondsSinceEpoch )->pixelLineToUnitVector( pixelLine );
    }

    Eigen::Matrix< double, 2, 3 > getPixelLinePartialWrtCameraFramePosition( const Eigen::Vector3d& positionOfObservedBodyInCameraFrame,
                                                                           const double secondsSinceEpoch = TUDAT_NAN ) const
    {
        return getProjectionModel( secondsSinceEpoch )->getPixelLinePartialWrtCameraFramePosition( positionOfObservedBodyInCameraFrame );
    }

    /*! \brief Calculate the pixel coordinates of a point, given its position relative to the camera in the inertial frame.
     *  The position is rotated to the camera frame with the camera orientation and pointing correction at the given time, and
     *  projected with the calibration at that time. For a camera fixed to its host body, the orientation is composed from the
     *  given inertial-to-body-fixed rotation; for a camera whose orientation is given directly w.r.t. the inertial frame, that
     *  rotation is not used. A time is only required for a camera with several pictures, or with a direct orientation.
     *  \param positionOfObservedBodyInInertialFrame Position of the observed point relative to the camera, in the inertial frame.
     *  \param rotationFromInertialToBodyFixedFrame Rotation from inertial frame to the body-fixed frame of the camera host body.
     *  \param secondsSinceEpoch Time at which the pixel coordinates are computed.
     *  \return The pixel coordinates of the observed point.
     */
    Eigen::Vector2d calculateObservableFromInertial( const Eigen::Vector3d& positionOfObservedBodyInInertialFrame,
                                                     const Eigen::Quaterniond& rotationFromInertialToBodyFixedFrame,
                                                     const double secondsSinceEpoch = TUDAT_NAN ) const
    {
        const Eigen::Quaterniond rotationFromInertialToCameraFrame = hasRotationFromInertialToCameraFrameFunction( )
                ? getRotationFromInertialToCameraFrame( secondsSinceEpoch )
                : getRotationFromInertialToCameraFrame( rotationFromInertialToBodyFixedFrame, secondsSinceEpoch );
        return getProjectionModel( secondsSinceEpoch )
                ->projectUnitVectorToPixelLine( rotationFromInertialToCameraFrame * positionOfObservedBodyInInertialFrame );
    }

    Eigen::Vector2d calculateObservableFromInertial( const Eigen::Vector3d& positionOfObservedBodyInInertialFrame,
                                                     const double secondsSinceEpoch ) const
    {
        return getProjectionModel( secondsSinceEpoch )
                ->projectUnitVectorToPixelLine( getRotationFromInertialToCameraFrame( secondsSinceEpoch ) *
                                                positionOfObservedBodyInInertialFrame );
    }

    //! Function to set the pictures taken by the camera, and reset all pointing corrections to zero.
    /*!
     *  The camera has one pointing correction per picture. Without pictures, a single correction applies at all times.
     *  \param pictureEpochs Epoch of each picture, in increasing order.
     *  \param pictureIds Identifier of each picture, or empty if the pictures have no identifier.
     *  \param pictureProjectionModels Projection model of each picture, for a camera whose calibration changes between pictures
     *  (for instance by pixel binning). If empty, the constant projection model of the camera is kept.
     */
    void setPictures( const std::vector< double >& pictureEpochs,
                      const std::vector< std::string >& pictureIds,
                      const std::vector< std::shared_ptr< CameraProjectionModel > >& pictureProjectionModels = {} )
    {
        if( ( !pictureIds.empty( ) && pictureIds.size( ) != pictureEpochs.size( ) ) ||
            ( !pictureProjectionModels.empty( ) && pictureProjectionModels.size( ) != pictureEpochs.size( ) ) )
        {
            throw std::runtime_error( "Error when setting pictures of camera " + cameraId_ +
                                      ": the number of picture identifiers and projection models must match the number of picture "
                                      "epochs." );
        }
        if( !pictureProjectionModels.empty( ) )
        {
            projectionModels_ = pictureProjectionModels;
        }
        pictureEpochs_ = Eigen::Map< const Eigen::VectorXd >( pictureEpochs.data( ), pictureEpochs.size( ) );
        pictureIds_ = pictureIds;
        pointingCorrection_ = Eigen::VectorXd::Zero( 3 * getNumberOfPointingCorrections( ) );
    }

    const Eigen::VectorXd& getPictureEpochs( ) const
    {
        return pictureEpochs_;
    }

    const std::vector< std::string >& getPictureIds( ) const
    {
        return pictureIds_;
    }

    int getNumberOfPointingCorrections( ) const
    {
        return std::max( 1, static_cast< int >( pictureEpochs_.size( ) ) );
    }

    //! Function to get the index of the picture closest in time to the given time (0 if the camera has at most one picture, in which
    //! case the time is not used).
    int getPictureIndex( const double secondsSinceEpoch ) const
    {
        if( pictureEpochs_.size( ) <= 1 )
        {
            return 0;
        }
        if( std::isnan( secondsSinceEpoch ) )
        {
            throw std::runtime_error( "Error when getting picture of camera " + cameraId_ +
                                      ": a time is required, since the camera has several pictures." );
        }
        return basic_mathematics::computeNearestNeighborUsingBinarySearch( pictureEpochs_, secondsSinceEpoch );
    }

    //! Function to get the pointing corrections of all pictures, as consecutive rotation vectors (axis times angle, in radians).
    Eigen::VectorXd getPointingCorrection( ) const
    {
        return pointingCorrection_;
    }

    void setPointingCorrection( const Eigen::VectorXd& pointingCorrection )
    {
        if( pointingCorrection.size( ) != 3 * getNumberOfPointingCorrections( ) )
        {
            throw std::runtime_error( "Error when setting pointing correction of camera " + cameraId_ + ": size must be " +
                                      std::to_string( 3 * getNumberOfPointingCorrections( ) ) + "." );
        }
        pointingCorrection_ = pointingCorrection;
    }

    //! Function to get the pointing correction of the picture closest in time to the given time.
    Eigen::Vector3d getPointingCorrection( const double secondsSinceEpoch ) const
    {
        return pointingCorrection_.segment< 3 >( 3 * getPictureIndex( secondsSinceEpoch ) );
    }

private:
    /*! \brief Apply the pointing correction of the picture closest to the given time to a nominal rotation to the camera frame.
     *  The correction is a small rotation expressed in the camera frame, and is therefore applied on the left:
     *  a direction that nominally maps to r_cam maps to Exp(theta) * r_cam once the correction is applied. This
     *  convention is the one differentiated by PixelCoordinatesScaling, and must be kept consistent with it.
     *  \param nominalRotationToCameraFrame Rotation to the camera frame without pointing correction.
     *  \param secondsSinceEpoch Time at which the rotation is evaluated.
     *  \return The rotation to the camera frame including the pointing correction.
     */
    Eigen::Quaterniond applyPointingCorrection_( const Eigen::Quaterniond& nominalRotationToCameraFrame,
                                                 const double secondsSinceEpoch ) const
    {
        const Eigen::Vector3d pointingCorrection = getPointingCorrection( secondsSinceEpoch );
        if( pointingCorrection.isZero( 0.0 ) )
        {
            return nominalRotationToCameraFrame;
        }
        return Eigen::Quaterniond( basic_mathematics::getQuaternionFromRotationVector( pointingCorrection ) * nominalRotationToCameraFrame )
                .normalized( );
    }

    //! Name of the camera
    std::string cameraId_;

    //! Rotation from body-fixed frame to camera frame
    Eigen::Quaterniond rotationFromBodyFixedToCameraFrame_;

    //! Projection model(s) mapping camera-frame directions to pixel/line coordinates: a single one for a constant calibration,
    //! or one per picture for a calibration that changes between pictures.
    std::vector< std::shared_ptr< CameraProjectionModel > > projectionModels_;

    //! Optional direct rotation from inertial to camera frame, used for picture-specific camera pointing data.
    std::function< Eigen::Quaterniond( const double ) > rotationFromInertialToCameraFrameFunction_;

    //! Reference body used by a body-relative direct pointing function.
    std::string pointingReferenceBodyName_;

    //! Epoch of each picture taken by the camera, in increasing order. Empty for a camera without pictures.
    Eigen::VectorXd pictureEpochs_;

    //! Identifier of each picture, or empty if the pictures have no identifier.
    std::vector< std::string > pictureIds_;

    //! Small rotation-vector correction of each picture, applied by camera pointing estimation.
    Eigen::VectorXd pointingCorrection_ = Eigen::VectorXd::Zero( 3 );
};
}  // namespace system_models
}  // namespace tudat

#endif
