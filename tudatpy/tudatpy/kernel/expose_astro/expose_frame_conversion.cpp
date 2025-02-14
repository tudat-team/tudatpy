/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_frame_conversion.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/ephemerides/rotationalEphemeris.h>
#include <tudat/astro/reference_frames.h>

namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;

namespace py = pybind11;

namespace tudatpy
{

namespace astro
{
namespace frame_conversion
{

void expose_frame_conversion( py::module &m )
{
    m.def( "inertial_to_rsw_rotation_matrix",
           &trf::getInertialToRswSatelliteCenteredFrameRotationMatrix,
           py::arg( "inertial_cartesian_state" ),
           R"doc(

Computes the rotation matrix from inertial to RSW frame.


Function to compute the rotation matrix from inertial to RSW frame.
The RSW frame is defined  by the state of a body w.r.t. to some
central body. The x-axis of the RSW frame points away from the
origin, and the y-axis lies in the orbital plane, and is positive
for in the direction of the velocity vector (but is not colinear
with the velocity vector, except for circular orbits). The z-axis
is perpendicular to the orbital plane, and completes the
right-handed coordinate system.


Parameters
----------
inertial_cartesian_state : numpy.ndarray
    Cartesian state, in an inertial frame, for which the rotation
    matrix is to be calculated. Note that the RSW frame is defined
    w.r.t. some central body, and this Cartesian state must be
    defined w.r.t. that central body (e.g. central body at the
    origin).

Returns
-------
numpy.ndarray
    Rotation matrix from inertial to RSW frame.






    )doc" );

    m.def( "rsw_to_inertial_rotation_matrix",
           &trf::getRswSatelliteCenteredToInertialFrameRotationMatrix,
           py::arg( "inertial_cartesian_state" ),
           R"doc(

Computes the rotation matrix from RSW to inertial frame.


Function to compute the rotation matrix from RSW to inertial. The
RSW frame is defined  by the state of a body w.r.t. to some central
body. The x-axis of the RSW frame points away from the origin, and
the y-axis lies in the orbital plane, and is positive for in the
direction of the velocity vector (but is not colinear with the
velocity vector, except for circular orbits). The z-axis is
perpendicular to the orbital plane, and completes the right-handed
coordinate system.


Parameters
----------
inertial_cartesian_state : numpy.ndarray
    Cartesian state, in an inertial frame, for which the rotation
    matrix is to be calculated. Note that the RSW frame is defined
    w.r.t. some central body, and this Cartesian state must be
    defined w.r.t. that central body (e.g. central body at the
    origin).

Returns
-------
numpy.ndarray
    Rotation matrix from RSW to inertial frame.






    )doc" );

    m.def( "tnw_to_inertial_rotation_matrix",
           py::overload_cast< const Eigen::Vector6d &, const bool >( &trf::getTnwToInertialRotation ),
           py::arg( "inertial_cartesian_state" ),
           py::arg( "n_axis_points_away_from_central_body" ) = true,
           R"doc(

Computes the rotation matrix from TNW to inertial frame.


Function to compute the rotation matrix from TNW to inertial frame.
The TNW frame is defined by the state of a body w.r.t. to some
central body. The x-axis of the TNW frame points along the velocity
vector, and the y-axis lies in the orbital plane, and is positive
in the direction away from the central body (or positive **towards**
the central body if the ``n_axis_points_away_from_central_body``
variable is set to false, see below). The z-axis is perpendicular
to the orbital plane, and completes the right-handed coordinate
system.


Parameters
----------
inertial_cartesian_state : numpy.ndarray
    Cartesian state, in an inertial frame, for which the rotation
    matrix is to be calculated. Note that the TNW frame is defined
    w.r.t. some central body, and this Cartesian state must be
    defined w.r.t. that central body (e.g. central body at the
    origin).

n_axis_points_away_from_central_body : bool
    Boolean (default=``True``) defining whether the N axis of the
    TNW frame points away from the central body (if ``True``) or
    towards the central body (if ``False``).

Returns
-------
numpy.ndarray
    Rotation matrix from TNW to inertial frame






    )doc" );

    m.def( "inertial_to_tnw_rotation_matrix",
           py::overload_cast< const Eigen::Vector6d &, const bool >( &trf::getInertialToTnwRotation ),
           py::arg( "inertial_cartesian_state" ),
           py::arg( "n_axis_points_away_from_central_body" ) = true,
           R"doc(

Computes the rotation matrix from inertial to TNW frame.


Function to compute the rotation matrix from inertial to TNW frame.
The TNW frame is defined by the state of a body w.r.t. to some
central body. The x-axis of the TNW frame points along the velocity
vector, and the y-axis lies in the orbital plane, and is positive
in the direction away from the central body (or positive **towards**
the central body if the ``n_axis_points_away_from_central_body``
variable is set to false, see below). The z-axis is perpendicular
to the orbital plane, and completes the right-handed coordinate
system.


Parameters
----------
inertial_cartesian_state : numpy.ndarray
    Cartesian state, in an inertial frame, for which the rotation
    matrix is to be calculated. Note that the RSW frame is defined
    w.r.t. some central body, and this Cartesian state must be
    defined w.r.t. that central body (e.g. central body at the
    origin).

n_axis_points_away_from_central_body : Boolean
    Boolean (default is ``True``) defining whether the N axis of the
    TNW frame points away from the central body (if ``True``) or
    towards the central body (if ``False``).

Returns
-------
numpy.ndarray
    Rotation matrix from inertial to TNW frame.






    )doc" );

    m.def( "inertial_to_body_fixed_rotation_matrix",
           py::overload_cast< const double, const double, const double >( &trf::getInertialToPlanetocentricFrameTransformationMatrix ),
           py::arg( "pole_declination" ),
           py::arg( "pole_right_ascension" ),
           py::arg( "prime_meridian_longitude" ),
           R"doc(

Computes the rotation matrix from inertial to body-fixed frame.


Function to compute the rotation matrix from inertial to body-fixed
frame, using typical pole right ascension (:math:`\alpha`), pole
declination (:math:`\delta`), and prime meridian longitude
(:math:`W`) angles.


Parameters
----------
pole_declination : float
    Declination of body pole in inertial frame (:math:`\delta`).

pole_right_ascension : float
    Right ascension of body pole in inertial frame (:math:`\alpha`).

prime_meridian_longitude : float
    Longitude of prime meridian w.r.t. intermediate frame
    (:math:`W`).

Returns
-------
numpy.ndarray
    Rotation matrix from inertial to body-fixed frame



Notes
-----
This definition of a body-fixed orientation is used by, for
instance, the IAU Working Group on Cartographic Coordinates and
Rotational Elements. Rotation is performed by a successive z-x-z
Euler angle rotation (see Archinal et al. [1]_).




    )doc" );

    m.def( "body_fixed_to_inertial_rotation_matrix",
           py::overload_cast< const double, const double, const double >(
                   &trf::getRotatingPlanetocentricToInertialFrameTransformationMatrix ),
           py::arg( "pole_declination" ),
           py::arg( "pole_right_ascension" ),
           py::arg( "pole_meridian" ),
           R"doc(

Computes the rotation matrix from body-fixed to inertial frame.


Function to compute the rotation matrix from body-fixed to inertial
frame, using typical pole right ascension (:math:`\alpha`), pole
declination (:math:`\delta`), and prime meridian longitude
(:math:`W`) angles.


Parameters
----------
pole_declination : float
    Declination of body pole in inertial frame (:math:`\delta`).

pole_right_ascension : float
    Right ascension of body pole in inertial frame (:math:`\alpha`).

prime_meridian_longitude : float
    Longitude of prime meridian w.r.t. intermediate frame
    (:math:`W`).

Returns
-------
numpy.ndarray
    Rotation matrix from body-fixed to inertial frame.




Notes
-----
This definition of a body-fixed orientation is used by,
for instance, the IAU Working Group on Cartographic Coordinates
and Rotational Elements. Rotation is performed by a successive z-x-z
Euler angle rotation (see Archinal et al. [1]_).




    )doc" );

    m.def( "rotate_state_to_frame",
           py::overload_cast< const Eigen::Vector6d &, const Eigen::Matrix3d &, const Eigen::Matrix3d & >(
                   &te::transformStateToFrameFromRotations< double > ),
           py::arg( "original_state" ),
           py::arg( "rotation_matrix" ),
           py::arg( "rotation_matrix_time_derivative" ) = Eigen::Matrix3d::Zero( ),
           R"doc(

Rotates a Cartesian state (position and velocity) from one frame :math:`B` to another frame :math:`A`, using the rotation matrix :math:`\mathbf{R}^{(A/B)}` from frame :math:`B` to :math:`A`, and its time derivative
:math:`\dot{\mathbf{R}}^{(A/B)}`.
This function computes:

.. math::
   \mathbf{r}^{(A)}=\mathbf{R}^{(A/B)}\mathbf{r}^{(B)}+\dot{\mathbf{R}}^{(A/B)}\mathbf{v}^{(B)}\\
   \mathbf{v}^{(A)}=\mathbf{R}^{(A/B)}\mathbf{v}^{(B)}

Parameters
----------
original_state : ndarray[numpy.float64[6, 1]]
    Cartesian state vector :math:`\mathbf{x}^{(B)}=[\mathbf{r}^{(B)};\mathbf{v}^{(B)}]` in frame :math:`B`
rotation_matrix: ndarray[numpy.float64[3, 3]]
    Rotation matrix :math:`\mathbf{R}^{(A/B)}` from frame :math:`B` to :math:`A`
rotation_matrix_time_derivative: ndarray[numpy.float64[3, 3]], default = numpy.zeros((3, 3))
    Time derivative of rotation matrix (:math:`\dot{\mathbf{R}}^{(A/B)})` from frame :math:`B` to :math:`A`; default zero indicates that frames :math:`A` and :math:`B` have a constant orientation w.r.t. one another.
Returns
-------
numpy.ndarray
    Input state in frame :math:`B`, rotated to frame :math:`A`


    )doc" );

    m.attr( "transform_cartesian_state_to_frame" ) = m.attr( "rotate_state_to_frame" );
}

}  // namespace frame_conversion
}  // namespace astro
}  // namespace tudatpy
