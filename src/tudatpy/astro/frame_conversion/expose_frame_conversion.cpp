/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/reference_frames.h"

#include "tudatpy/docstrings.h"

namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;

namespace py = pybind11;

namespace tudatpy {
    namespace astro {
        namespace frame_conversion {

            PYBIND11_MODULE(expose_frame_conversion, m) {
                m.def(
                    "inertial_to_rsw_rotation_matrix",
                    &trf::getInertialToRswSatelliteCenteredFrameRotationMatrix,
                    py::arg("inertial_cartesian_state"),
R"doc(Computes the rotation matrix from inertial to RSW frame.


	Function to compute the rotation matrix from inertial to RSW frame.
	The RSW frame is defined  by the state of a body w.r.t. to some
	central body. The x-axis of the RSW frame points away from the
	origin, and the y-axis lies in the orbital plane, and is positive
	for in the direction of the velocity vector (but is not colinear
	with the velocity vector, except for circular orbits). The z-axis
	is perpendicular to the orbital plane, and completes the
	right-handed coordinate system.
	

	:param inertial_cartesian_state:
		Cartesian state, in an inertial frame, for which the rotation
		matrix is to be calculated. Note that the RSW frame is defined
		w.r.t. some central body, and this Cartesian state must be
		defined w.r.t. that central body (e.g. central body at the
		origin).
		
	:return:
		Rotation matrix from inertial to RSW frame.
)doc");

                m.def(
                    "rsw_to_inertial_rotation_matrix",
                    &trf::getRswSatelliteCenteredToInertialFrameRotationMatrix,
                    py::arg("inertial_cartesian_state"),
R"doc(Computes the rotation matrix from RSW to inertial frame.


	Function to compute the rotation matrix from RSW to inertial. The
	RSW frame is defined  by the state of a body w.r.t. to some central
	body. The x-axis of the RSW frame points away from the origin, and
	the y-axis lies in the orbital plane, and is positive for in the
	direction of the velocity vector (but is not colinear with the
	velocity vector, except for circular orbits). The z-axis is
	perpendicular to the orbital plane, and completes the right-handed
	coordinate system.
	

	:param inertial_cartesian_state:
		Cartesian state, in an inertial frame, for which the rotation
		matrix is to be calculated. Note that the RSW frame is defined
		w.r.t. some central body, and this Cartesian state must be
		defined w.r.t. that central body (e.g. central body at the
		origin).
		
	:return:
		Rotation matrix from RSW to inertial frame.
)doc");

                m.def("tnw_to_inertial_rotation_matrix",
                      py::overload_cast<const Eigen::Vector6d &, const bool>(
                          &trf::getTnwToInertialRotation),
                      py::arg("inertial_cartesian_state"),
                      py::arg("n_axis_points_away_from_central_body") = true,
R"doc(Computes the rotation matrix from TNW to inertial frame.


	Function to compute the rotation matrix from TNW to inertial frame.
	The TNW frame is defined by the state of a body w.r.t. to some
	central body. The x-axis of the TNW frame points along the velocity
	vector, and the y-axis lies in the orbital plane, and is positive
	in the direction away from the central body (or positive **towards**
	the central body if the ``n_axis_points_away_from_central_body``
	variable is set to false, see below). The z-axis is perpendicular
	to the orbital plane, and completes the right-handed coordinate
	system.
	

	:param inertial_cartesian_state:
		Cartesian state, in an inertial frame, for which the rotation
		matrix is to be calculated. Note that the TNW frame is defined
		w.r.t. some central body, and this Cartesian state must be
		defined w.r.t. that central body (e.g. central body at the
		origin).
		
	:param n_axis_points_away_from_central_body:
		Boolean (default=``True``) defining whether the N axis of the
		TNW frame points away from the central body (if ``True``) or
		towards the central body (if ``False``).
		
	:return:
		Rotation matrix from TNW to inertial frame
)doc");

                m.def("inertial_to_tnw_rotation_matrix",
                      py::overload_cast<const Eigen::Vector6d &, const bool>(
                          &trf::getInertialToTnwRotation),
                      py::arg("inertial_cartesian_state"),
                      py::arg("n_axis_points_away_from_central_body") = true,
R"doc(Computes the rotation matrix from inertial to TNW frame.


	Function to compute the rotation matrix from inertial to TNW frame.
	The TNW frame is defined by the state of a body w.r.t. to some
	central body. The x-axis of the TNW frame points along the velocity
	vector, and the y-axis lies in the orbital plane, and is positive
	in the direction away from the central body (or positive **towards**
	the central body if the ``n_axis_points_away_from_central_body``
	variable is set to false, see below). The z-axis is perpendicular
	to the orbital plane, and completes the right-handed coordinate
	system.
	

	:param inertial_cartesian_state:
		Cartesian state, in an inertial frame, for which the rotation
		matrix is to be calculated. Note that the RSW frame is defined
		w.r.t. some central body, and this Cartesian state must be
		defined w.r.t. that central body (e.g. central body at the
		origin).
		
	:param n_axis_points_away_from_central_body:
		Boolean (default is ``True``) defining whether the N axis of the
		TNW frame points away from the central body (if ``True``) or
		towards the central body (if ``False``).
		
	:return:
		Rotation matrix from inertial to TNW frame.
)doc");


                m.def(
                    "inertial_to_body_fixed_rotation_matrix",
                    py::overload_cast<const double, const double, const double>(
                        &trf::
                            getInertialToPlanetocentricFrameTransformationMatrix),
                    py::arg("pole_declination"),
                    py::arg("pole_right_ascension"),
                    py::arg("prime_meridian_longitude"),
R"doc(Computes the rotation matrix from inertial to body-fixed frame.


	Function to compute the rotation matrix from inertial to body-fixed
	frame, using typical pole right ascension (:math:`\alpha`), pole
	declination (:math:`\delta`), and prime meridian longitude
	(:math:`W`) angles.
	

	:param pole_declination:
		Declination of body pole in inertial frame (:math:`\delta`).
		
	:param pole_right_ascension:
		Right ascension of body pole in inertial frame (:math:`\alpha`).
		
	:param prime_meridian_longitude:
		Longitude of prime meridian w.r.t. intermediate frame
		(:math:`W`).
		
	:return:
		Rotation matrix from inertial to body-fixed frame
)doc");

                m.def(
                    "body_fixed_to_inertial_rotation_matrix",
                    py::overload_cast<const double, const double, const double>(
                        &trf::
                            getRotatingPlanetocentricToInertialFrameTransformationMatrix),
                    py::arg("pole_declination"),
                    py::arg("pole_right_ascension"), py::arg("pole_meridian"),
R"doc(Computes the rotation matrix from body-fixed to inertial frame.


	Function to compute the rotation matrix from body-fixed to inertial
	frame, using typical pole right ascension (:math:`\alpha`), pole
	declination (:math:`\delta`), and prime meridian longitude
	(:math:`W`) angles.
	

	:param pole_declination:
		Declination of body pole in inertial frame (:math:`\delta`).
		
	:param pole_right_ascension:
		Right ascension of body pole in inertial frame (:math:`\alpha`).
		
	:param prime_meridian_longitude:
		Longitude of prime meridian w.r.t. intermediate frame
		(:math:`W`).
		
	:return:
		Rotation matrix from body-fixed to inertial frame.
		
)doc");


                m.def("transform_cartesian_state_to_frame",
                      py::overload_cast<const Eigen::Vector6d &,
                                        const Eigen::Matrix3d &,
                                        const Eigen::Matrix3d &>(
                          &te::transformStateToFrameFromRotations<double>),
                      py::arg("original_state"), py::arg("rotation_matrix"),
                      py::arg("transform_cartesian_state_to_frame"));
            }

        }  // namespace frame_conversion
    }  // namespace astro
}  // namespace tudatpy
