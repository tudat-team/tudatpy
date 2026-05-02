/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createCameras.h"
#include "tudat/math/basic/rotationRepresentations.h"
#include "tudat/astro/system_models/vehicleSystems.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a camera from pre-defined camera settings object, and add it to a Body object
void createCamera( const std::shared_ptr< Body > body, const std::shared_ptr< CameraSettings > cameraSettings )
{
    Eigen::Quaterniond rotationQuaternion =
            basic_mathematics::getQuaternionFrom313EulerAngles( cameraSettings->getCamera313EulerAngles( ) );
    std::pair< double, double > focalLengths = cameraSettings->getFocalLengths( );
    std::pair< double, double > opticalCenter = cameraSettings->getOpticalCenter( );
    std::shared_ptr< tudat::ephemerides::RotationalEphemeris > rotationalEphemeris = body->getRotationalEphemeris( );
    if( rotationalEphemeris == nullptr )
    {
        throw std::runtime_error(
                "Error when creating camera: body does not have a rotational ephemeris."
                "If rotational states are being propagated, make sure to add an empty tabulated rotational ephemeris to the body before "
                "creating the camera." );
    }
    std::shared_ptr< system_models::Camera > camera = std::make_shared< system_models::Camera >(
            cameraSettings->getCameraName( ), rotationQuaternion, focalLengths, opticalCenter, rotationalEphemeris );
    body->getVehicleSystems( )->addCamera( cameraSettings->getCameraName( ), camera );
}

//! Function to create a camera and add it to a Body object
void createCamera( const std::shared_ptr< Body > body,
                   const std::string& cameraName,
                   const Eigen::Vector3d& boresightEulerAngles,
                   const std::pair< double, double > focalLengths,
                   const std::pair< double, double > opticalCenter )
{
    std::shared_ptr< CameraSettings > cameraSettings =
            std::make_shared< CameraSettings >( cameraName, boresightEulerAngles, focalLengths, opticalCenter );
    createCamera( body, cameraSettings );
}

}  // namespace simulation_setup

}  // namespace tudat
