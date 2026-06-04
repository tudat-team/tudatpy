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

#include <string>
#include <functional>
#include <memory>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <tudat/astro/ephemerides/rotationalEphemeris.h>

namespace tudat
{

namespace system_models
{

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
            std::shared_ptr< tudat::ephemerides::RotationalEphemeris > bodyRotationalEphemeris = nullptr ):
        cameraId_( cameraId ), rotationFromBodyFixedToCameraFrame_( rotationFromBodyFixedToCameraFrame ),
        K_( focal_lengths.first, focal_lengths.second ), opticalCenter_( optical_center.first, optical_center.second ),
        bodyRotationalEphemeris_( bodyRotationalEphemeris )
    {}

    /*! \brief Get the camera ID.
     *  \return The camera ID.
     */
    std::string getCameraId( ) const
    {
        return cameraId_;
    }

    Eigen::DiagonalMatrix< double, 2 > getFocalLengthsMatrix( )
    {
        return K_;
    }

    Eigen::Vector2d getOpticalCenter( )
    {
        return opticalCenter_;
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
        return opticalCenter_ +
                K_ *
                Eigen::Vector2d( positionInCameraFrame.x( ) / positionInCameraFrame.z( ),
                                 positionInCameraFrame.y( ) / positionInCameraFrame.z( ) );
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
            throw std::runtime_error(
                    "Error when calculating observable: no rotational ephemeris provided for the body using the camera." );
        }
        Eigen::Quaterniond rotationFromInertialToBodyFixed = bodyRotationalEphemeris_->getRotationToTargetFrame( secondsSinceEpoch );
        Eigen::Vector3d positionOfObservedBodyInBodyFrame = rotationFromInertialToBodyFixed * positionOfObservedBodyInInertialFrame;
        return calculateObservableFromBodyFixed( positionOfObservedBodyInBodyFrame );
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

    //! Focal lengths diagonal matrix of the camera, containing the focal lengths in the x and y directions, respectively.
    Eigen::DiagonalMatrix< double, 2 > K_;

    //! Optical center of the camera
    Eigen::Vector2d opticalCenter_;

    //! Shared pointer to a rotational ephemeris of the body, which can be used to determine the rotation from inertial to body fixed frame at a given time.
    std::shared_ptr< tudat::ephemerides::RotationalEphemeris > bodyRotationalEphemeris_;
};
}  // namespace system_models
}  // namespace tudat

#endif