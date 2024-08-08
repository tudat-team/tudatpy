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
#include <tudat/astro/ephemerides/rotationalEphemeris.h>
#include <tudat/astro/reference_frames.h>

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
                    get_docstring("inertial_to_rsw_rotation_matrix").c_str());

                m.def(
                    "rsw_to_inertial_rotation_matrix",
                    &trf::getRswSatelliteCenteredToInertialFrameRotationMatrix,
                    py::arg("inertial_cartesian_state"),
                    get_docstring("rsw_to_inertial_rotation_matrix").c_str());

                m.def("tnw_to_inertial_rotation_matrix",
                      py::overload_cast<const Eigen::Vector6d &, const bool>(
                          &trf::getTnwToInertialRotation),
                      py::arg("inertial_cartesian_state"),
                      py::arg("n_axis_points_away_from_central_body") = true,
                      get_docstring("tnw_to_inertial_rotation_matrix").c_str());

                m.def("inertial_to_tnw_rotation_matrix",
                      py::overload_cast<const Eigen::Vector6d &, const bool>(
                          &trf::getInertialToTnwRotation),
                      py::arg("inertial_cartesian_state"),
                      py::arg("n_axis_points_away_from_central_body") = true,
                      get_docstring("inertial_to_tnw_rotation_matrix").c_str());


                m.def(
                    "inertial_to_body_fixed_rotation_matrix",
                    py::overload_cast<const double, const double, const double>(
                        &trf::
                            getInertialToPlanetocentricFrameTransformationMatrix),
                    py::arg("pole_declination"),
                    py::arg("pole_right_ascension"),
                    py::arg("prime_meridian_longitude"),
                    get_docstring("inertial_to_body_fixed_rotation_matrix")
                        .c_str());

                m.def(
                    "body_fixed_to_inertial_rotation_matrix",
                    py::overload_cast<const double, const double, const double>(
                        &trf::
                            getRotatingPlanetocentricToInertialFrameTransformationMatrix),
                    py::arg("pole_declination"),
                    py::arg("pole_right_ascension"), py::arg("pole_meridian"),
                    get_docstring("body_fixed_to_inertial_rotation_matrix")
                        .c_str());


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
