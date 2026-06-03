/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATECAMERAS_H
#define TUDAT_CREATECAMERAS_H

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/astro/system_models/camera.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace simulation_setup
{

//! Class to hold settings for creating a camera
/*!
 *  Class to hold settings for creating a camera, which can then be passed to createCamera function.
 *  This class allows setting the camera name, Euler angles with respect to the body-fixed frame, focal lengths, optical center, and
 *  body-fixed camera position.
 *
 *  Note that the camera is assumed to be fixed on the body and to point along the z-axis of its own frame, i.e. the focal plane is located
 * at z = 1 in the camera frame.
 */
class CameraSettings
{
public:
    /*! \brief Constructor, sets the camera name, rotation from body-fixed to camera frame, and focal lengths in x and y directions.
     *  Euler angles in 3-2-3 convention representing the right ascension, declination and twist angle of camera boresight.
     *  \param cameraName Name of the camera.
     *  \param boresightEulerAngles Euler angles of camera boresight with respect to body-fixed frame [rad].
     *  \param focalLengths Focal lengths in x and y directions [px].
     *  \param opticalCenter Optical center of the camera [px].
     *  \param bodyFixedCameraPosition Position of the camera in the body-fixed frame [m].
     */
    CameraSettings( const std::string& cameraName,
                    const Eigen::Vector3d& boresightEulerAngles,
                    const std::pair< double, double > focalLengths = std::make_pair( 1.0, 1.0 ),
                    const std::pair< double, double > opticalCenter = std::make_pair( 0.0, 0.0 ),
                    const Eigen::Vector3d& bodyFixedCameraPosition = Eigen::Vector3d::Zero( ) ):
        cameraName_( cameraName ), boresightEulerAngles_( boresightEulerAngles ), focalLengths_( focalLengths ),
        opticalCenter_( opticalCenter ), bodyFixedCameraPosition_( bodyFixedCameraPosition )
    {}

    std::string getCameraName( )
    {
        return cameraName_;
    }

    void setBoresightEulerAngles( const Eigen::Vector3d& boresightEulerAngles )
    {
        boresightEulerAngles_ = boresightEulerAngles;
    }

    Eigen::Vector3d getBoresightEulerAngles( )
    {
        return boresightEulerAngles_;
    }

    Eigen::Vector3d getCamera313EulerAngles( )
    {
        Eigen::Vector3d cameraEulerAngles = boresightEulerAngles_;
        cameraEulerAngles( 0 ) = boresightEulerAngles_( 2 );
        cameraEulerAngles( 1 ) = -boresightEulerAngles_( 1 ) + mathematical_constants::PI / 2.0;
        cameraEulerAngles( 2 ) = boresightEulerAngles_( 0 ) + mathematical_constants::PI / 2.0;
        return cameraEulerAngles;
    }

    std::pair< double, double > getFocalLengths( )
    {
        return focalLengths_;
    }

    void setFocalLengths( const std::pair< double, double >& focalLengths )
    {
        focalLengths_ = focalLengths;
    }

    std::pair< double, double > getOpticalCenter( )
    {
        return opticalCenter_;
    }

    void setOpticalCenter( const std::pair< double, double >& opticalCenter )
    {
        opticalCenter_ = opticalCenter;
    }

    Eigen::Vector3d getBodyFixedCameraPosition( )
    {
        return bodyFixedCameraPosition_;
    }

    void setBodyFixedCameraPosition( const Eigen::Vector3d& bodyFixedCameraPosition )
    {
        bodyFixedCameraPosition_ = bodyFixedCameraPosition;
    }

protected:
    std::string cameraName_;

    Eigen::Vector3d boresightEulerAngles_;

    std::pair< double, double > focalLengths_;

    std::pair< double, double > opticalCenter_;

    Eigen::Vector3d bodyFixedCameraPosition_;
};

//! Function to create a shared pointer to CameraSettings object
/*!
 *  Function to create a shared pointer to CameraSettings object, which can then be passed to createCamera function.
 *  This function allows setting the camera name, Euler angles with respect to the body-fixed frame, focal lengths, optical center, and
 *  body-fixed camera position.
 *  Note that the camera is assumed to be fixed on the body and to point along the z-axis of its own frame, i.e. the focal plane is located
 * at z = 1 in the camera frame.
 *  Euler angles are in 3-2-3 convention representing the right ascension, declination and twist angle of camera boresight.
 *
 *  \param cameraName Name of the camera.
 *  \param boresightEulerAngles Euler angles of camera boresight with respect to
 * body-fixed frame.
 *  \param focalLengths Focal lengths in x and y directions [px].
 *  \param opticalCenter Optical center of the camera [px].
 *  \param bodyFixedCameraPosition Position of the camera in the body-fixed frame [m].
 *  \return Shared pointer to CameraSettings object created from input parameters.
 */
inline std::shared_ptr< CameraSettings > cameraSettings( const std::string& cameraName,
                                                         const Eigen::Vector3d& boresightEulerAngles,
                                                         const std::pair< double, double > focalLengths = std::make_pair( 1.0, 1.0 ),
                                                         const std::pair< double, double > opticalCenter = std::make_pair( 0.0, 0.0 ),
                                                         const Eigen::Vector3d& bodyFixedCameraPosition = Eigen::Vector3d::Zero( ) )
{
    return std::make_shared< CameraSettings >( cameraName, boresightEulerAngles, focalLengths, opticalCenter, bodyFixedCameraPosition );
}

//! Function to create a camera and add it to a Body object
/*!
 * Function to create a camera and add it to a Body object
 * \param body Body object in which the newly created camera is to be added.
 * \param cameraName Name of camera that is to be created.
 * \param boresightEulerAngles Euler angles of camera boresight with respect to body-fixed frame.
 * \param focalLengths Focal lengths in x and y directions [px].
 * \param opticalCenter Optical center of the camera [px].
 * \param bodyFixedCameraPosition Position of the camera in the body-fixed frame [m].
 */
void createCamera( const std::shared_ptr< Body > body,
                   const std::string& cameraName,
                   const Eigen::Vector3d& boresightEulerAngles,
                   const std::pair< double, double > focalLengths = std::make_pair( 1.0, 1.0 ),
                   const std::pair< double, double > opticalCenter = std::make_pair( 0.0, 0.0 ),
                   const Eigen::Vector3d& bodyFixedCameraPosition = Eigen::Vector3d::Zero( ) );

//! Function to create a camera and add it to a Body object
/*!
 * Function to create a camera and add it to a Body object
 * \param body Body object in which the newly created camera is to be added.
 * \param cameraSettings Camera settings that are to be used for creating the camera
 */
void createCamera( const std::shared_ptr< Body > body, const std::shared_ptr< CameraSettings > cameraSettings );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATECAMERAS_H
