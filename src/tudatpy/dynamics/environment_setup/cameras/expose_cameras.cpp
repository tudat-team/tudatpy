/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_cameras.h"

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include <pybind11/native_enum.h>
#include <tudat/simulation/environment_setup/createCameras.h>
#include <tudat/simulation/environment_setup/defaultBodies.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tcc = tudat::coordinate_conversions;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace cameras
{

void expose_cameras_setup( py::module& m )
{
    py::class_< tss::CameraSettings, std::shared_ptr< tss::CameraSettings > >( m,
                                                                               "CameraSettings",
                                                                               R"doc(

         Base class for providing settings for the creation of a camera.


      )doc" )
            .def_property(
                    "boresight_euler_angles", &tss::CameraSettings::getBoresightEulerAngles, &tss::CameraSettings::setBoresightEulerAngles )
            .def_property( "focal_lengths", &tss::CameraSettings::getFocalLengths, &tss::CameraSettings::setFocalLengths )
            .def_property( "optical_center", &tss::CameraSettings::getOpticalCenter, &tss::CameraSettings::setOpticalCenter )

            .def_property_readonly( "camera_name", &tss::CameraSettings::getCameraName );

    m.def( "pinhole_camera",
           &tss::cameraSettings,
           py::arg( "camera_name" ),
           py::arg( "boresight_euler_angles" ),
           py::arg( "focal_lengths" ) = std::make_pair( 1.0, 1.0 ),
           py::arg( "optical_center" ) = std::make_pair( 0.0, 0.0 ),
           R"doc(

 Function for creating settings for a camera

 Function for creating settings for a camera, defining only its name and orientation.
 Camera is fixed to the body and points along the z-axis of its own frame.
 Represents a pinhole camera model without any distortions added.
 The orientation of the camera is defined by the Euler angles of the camera boresight (z-axis) with respect to the body-fixed frame, in a 3-2-3 rotation sequence.
 A zero twist angle will make the body-fixed frame x-axis aligned with the pixels u direction (positive horizontal direction).
   The focal lengths and optical center of the camera can also be defined, but default to (1.0, 1.0) and (0.0, 0.0), respectively.

 Parameters
 ----------
 camera_name : string
    Name (unique identifier) by which the camera is to be known.
 boresight_euler_angles : numpy.ndarray([3,1])
   [RA, DEC, Twist] angles in radians, following a 3-2-3 rotation sequence.
 focal_lengths : tuple[float, float], optional
    Tuple of the focal lengths of the camera in the x and y directions, respectively. Default is (1.0, 1.0)
 optical_center : tuple[float, float], optional
    Tuple of the optical center coordinates of the camera in the x and y directions, respectively. Default is (0.0, 0.0)
 Returns
 -------
 CameraSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.camera.CameraSettings` defining settings of the to be created camera





 Examples
 --------
 In this example, we create a basic camera settings aligned with y axis:

    .. code-block:: python
    
        from tudatpy.dynamics.environment_setup.cameras import pinhole_camera
    
        camera_settings = pinhole_camera( camera_name = "Camera1", boresight_euler_angles = [ 0.0, np.pi / 2.0, 0.0 ] )

 )doc" );
}

}  // namespace cameras
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
