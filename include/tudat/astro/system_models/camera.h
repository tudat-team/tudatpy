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

#include <cmath>
#include <string>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <tudat/astro/ephemerides/rotationalEphemeris.h>

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
        if( std::fabs( cameraFrameVector.z( ) ) <= std::numeric_limits< double >::epsilon( ) )
        {
            throw std::runtime_error( "Error when projecting camera-frame vector: z-component is zero." );
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
    {}

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
    void validateProjectionDirection( const Eigen::Vector3d& cameraFrameVector ) const
    {
        if( std::fabs( focalLength_ ) <= std::numeric_limits< double >::epsilon( ) )
        {
            throw std::runtime_error( "Error when projecting with PSF camera model: focal length is zero." );
        }
        if( std::fabs( cameraFrameVector.z( ) ) <= std::numeric_limits< double >::epsilon( ) )
        {
            throw std::runtime_error( "Error when projecting PSF camera-frame vector: z-component is zero." );
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
     *  \param bodyRotationalEphemeris Rotational ephemeris of the body to which the camera is fixed.
     *  which can be used to determine the rotation from inertial to body fixed frame at a given time.
     */
    Camera( const std::string& cameraId,
            const Eigen::Quaterniond& rotationFromBodyFixedToCameraFrame,
            std::pair< double, double > focal_lengths = std::make_pair( 1.0, 1.0 ),
            std::pair< double, double > optical_center = std::make_pair( 0.0, 0.0 ),
            std::shared_ptr< tudat::ephemerides::RotationalEphemeris > bodyRotationalEphemeris = nullptr,
            std::function< Eigen::Quaterniond( const double ) > rotationFromInertialToCameraFrameFunction = nullptr ):
        cameraId_( cameraId ), rotationFromBodyFixedToCameraFrame_( rotationFromBodyFixedToCameraFrame ),
        projectionModel_( std::make_shared< PinholeCameraProjectionModel >( focal_lengths, optical_center ) ),
        bodyRotationalEphemeris_( bodyRotationalEphemeris ),
        rotationFromInertialToCameraFrameFunction_( rotationFromInertialToCameraFrameFunction )
    {}

    Camera( const std::string& cameraId,
            const Eigen::Quaterniond& rotationFromBodyFixedToCameraFrame,
            const std::shared_ptr< CameraProjectionModel > projectionModel,
            std::shared_ptr< tudat::ephemerides::RotationalEphemeris > bodyRotationalEphemeris = nullptr,
            std::function< Eigen::Quaterniond( const double ) > rotationFromInertialToCameraFrameFunction = nullptr ):
        cameraId_( cameraId ), rotationFromBodyFixedToCameraFrame_( rotationFromBodyFixedToCameraFrame ),
        projectionModel_( projectionModel ), bodyRotationalEphemeris_( bodyRotationalEphemeris ),
        rotationFromInertialToCameraFrameFunction_( rotationFromInertialToCameraFrameFunction )
    {
        if( projectionModel_ == nullptr )
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

    Eigen::DiagonalMatrix< double, 2 > getFocalLengthsMatrix( )
    {
        return projectionModel_->getFocalLengthsMatrix( );
    }

    Eigen::Vector2d getOpticalCenter( )
    {
        return projectionModel_->getOpticalCenter( );
    }

    std::shared_ptr< CameraProjectionModel > getProjectionModel( )
    {
        return projectionModel_;
    }

    /*! \brief Get quaternion representing active rotation from body-fixed to camera frame.
     *  \return The rotation from body-fixed to camera frame.
     */
    Eigen::Quaterniond getRotationFromBodyFixedToCameraFrame( ) const
    {
        return rotationFromBodyFixedToCameraFrame_;
    }

    Eigen::Quaterniond getRotationFromInertialToCameraFrame( const double secondsSinceEpoch ) const
    {
        if( rotationFromInertialToCameraFrameFunction_ )
        {
            return rotationFromInertialToCameraFrameFunction_( secondsSinceEpoch );
        }
        if( bodyRotationalEphemeris_ == nullptr )
        {
            throw std::runtime_error(
                    "Error when calculating rotation from inertial to camera frame: no rotational ephemeris provided for the body "
                    "using the camera." );
        }
        Eigen::Quaterniond rotationFromInertialToBodyFixed = bodyRotationalEphemeris_->getRotationToTargetFrame( secondsSinceEpoch );
        return rotationFromBodyFixedToCameraFrame_ * rotationFromInertialToBodyFixed;
    }

    /*! \brief Calculate the observable position of a body in the camera frame.
     *  Calculates the pixel coordinates of the observed body. Uses gnomonic projection, i.e the camera boresight points along the z-axis of
     * its own frame, and the focal plane is located at the camera frame x-y plane. The optical center and focal lengths are then used to
     * calculate the pixel coordinates from the position in the camera frame.
     * \param positionOfObservedBodyInBodyFrame Position of the
     * observed body in the body-fixed frame.
     * \return The observable position in the camera frame.
     */
    Eigen::Vector2d calculateObservableFromBodyFixed( const Eigen::Vector3d& positionOfObservedBodyInBodyFrame ) const
    {
        Eigen::Vector3d positionInCameraFrame = bodyFixedToCameraFrame_( positionOfObservedBodyInBodyFrame );
        return projectionModel_->projectUnitVectorToPixelLine( positionInCameraFrame );
    }

    Eigen::Vector3d pixelLineToCameraFrameUnitVector( const Eigen::Vector2d& pixelLine ) const
    {
        return projectionModel_->pixelLineToUnitVector( pixelLine );
    }

    Eigen::Matrix< double, 2, 3 > getPixelLinePartialWrtCameraFramePosition(
            const Eigen::Vector3d& positionOfObservedBodyInCameraFrame ) const
    {
        return projectionModel_->getPixelLinePartialWrtCameraFramePosition( positionOfObservedBodyInCameraFrame );
    }

    /*! \brief Calculate the observable position of a body in the camera frame, given its position in the observer centered inertial frame and the time at which this position is valid.
     *  \param positionOfObservedBodyInInertialFrame Position of the observed body in the observer centered inertial frame.
     *  \param secondsSinceEpoch Seconds since Julian day reference epoch at which the position of the observed body is valid.
     *  \return The observable position in the camera frame.
     */
    Eigen::Vector2d calculateObservableFromInertial( const Eigen::Vector3d& positionOfObservedBodyInInertialFrame,
                                                     const double secondsSinceEpoch ) const
    {
        if( bodyRotationalEphemeris_ == nullptr )
        {
            if( !rotationFromInertialToCameraFrameFunction_ )
            {
                throw std::runtime_error(
                        "Error when calculating observable: no rotational ephemeris or direct inertial-to-camera rotation function "
                        "provided "
                        "for the body using the camera." );
            }
        }
        Eigen::Vector3d positionOfObservedBodyInCameraFrame =
                getRotationFromInertialToCameraFrame( secondsSinceEpoch ) * positionOfObservedBodyInInertialFrame;
        return projectionModel_->projectUnitVectorToPixelLine( positionOfObservedBodyInCameraFrame );
    }

private:
    /*! \brief Calculate the position of a body in the camera frame from its position in the body-fixed frame.
     *  \param positionOfObservedBodyInBodyFrame Position of the observed body in the body-fixed frame.
     *  \return The position in the camera frame.
     */
    Eigen::Vector3d bodyFixedToCameraFrame_( const Eigen::Vector3d& positionOfObservedBodyInBodyFrame ) const
    {
        return rotationFromBodyFixedToCameraFrame_ * positionOfObservedBodyInBodyFrame;
    }

    //! Name of the camera
    std::string cameraId_;

    //! Rotation from body-fixed frame to camera frame
    Eigen::Quaterniond rotationFromBodyFixedToCameraFrame_;

    //! Projection model mapping camera-frame directions to pixel/line coordinates.
    std::shared_ptr< CameraProjectionModel > projectionModel_;

    //! Shared pointer to a rotational ephemeris of the body, which can be used to determine the rotation from inertial to body fixed frame at a given time.
    std::shared_ptr< tudat::ephemerides::RotationalEphemeris > bodyRotationalEphemeris_;

    //! Optional direct rotation from inertial to camera frame, used for picture-specific camera pointing data.
    std::function< Eigen::Quaterniond( const double ) > rotationFromInertialToCameraFrameFunction_;
};
}  // namespace system_models
}  // namespace tudat

#endif
