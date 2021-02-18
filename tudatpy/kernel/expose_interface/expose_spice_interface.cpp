/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_spice_interface.h"

#include "../docstrings.h"
#include "tudat/interface/spice.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tsi = tudat::spice_interface;

namespace tudatpy {

void expose_spice_interface(py::module &m) {

  // time related
  m.def("convert_julian_date_to_ephemeris_time",
        &tudat::spice_interface::convertJulianDateToEphemerisTime,
        py::arg("julian_date"),
        "Convert a Julian date to ephemeris time (equivalent to TDB in Spice).");

//  m.def("jd2tdb", m.attr("convert_julian_date_to_ephemeris_time"));

  m.def("convert_ephemeris_time_to_julian_date",
        &tudat::spice_interface::convertEphemerisTimeToJulianDate,
        py::arg("ephemeris_time"),
        "Convert ephemeris time (equivalent to TDB) to a Julian date.");

//  m.def("tdb2jd", m.attr("convert_ephemeris_time_to_julian_date"));

  m.def("convert_date_string_to_ephemeris_time",
        &tudat::spice_interface::convertDateStringToEphemerisTime,
        py::arg("date_string"),
        "Converts a date string to ephemeris time.");

//  m.def("dstr2jd", m.attr("convert_date_string_to_ephemeris_time"));

  // positional state related
  m.def("get_body_cartesian_position_at_epoch",
        &tudat::spice_interface::getBodyCartesianPositionAtEpoch,
        py::arg("target_body_name"),
        py::arg("observer_body_name"),
        py::arg("reference_frame_name"),
        py::arg("aberration_corrections"),
        py::arg("ephemeris_time"),
        "Get Cartesian position of a body, as observed from another body.");

  m.def("get_body_cartesian_state_at_epoch",
        &tudat::spice_interface::getBodyCartesianStateAtEpoch,
        py::arg("target_body_name"),
        py::arg("observer_body_name"),
        py::arg("reference_frame_name"),
        py::arg("aberration_corrections"),
        py::arg("ephemeris_time"),
        "Get Cartesian position of a body, as observed from another body.");

  m.def("get_cartesian_state_from_tle_at_epoch",
        &tudat::spice_interface::getCartesianStateFromTleAtEpoch,
        py::arg("epoch"),
        py::arg("tle"),
        "Get Cartesian state of a satellite from its two-line element set at a specified epoch.");

  // rotational state related
  m.def("compute_rotation_matrix_between_frames",
        &tudat::spice_interface::computeRotationMatrixBetweenFrames,
        py::arg("original_frame"),
        py::arg("new_frame"),
        py::arg("ephemeris_time"),
        "Compute matrix of rotation between two frames.");

  m.def("compute_rotation_matrix_derivative_between_frames",
        &tudat::spice_interface::computeRotationMatrixDerivativeBetweenFrames,
        py::arg("original_frame"),
        py::arg("new_frame"),
        py::arg("ephemeris_time"),
        "Computes time derivative of rotation matrix between two frames.");

  m.def("get_angular_velocity_vector_of_frame_in_original_frame",
        &tudat::spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame,
        py::arg("original_frame"),
        py::arg("new_frame"),
        py::arg("ephemeris_time"),
        "Computes the angular velocity of one frame w.r.t. to another frame.");

  m.def("compute_rotation_quaternion_and_rotation_matrix_derivative_between_frames",
        &tudat::spice_interface::computeRotationQuaternionAndRotationMatrixDerivativeBetweenFrames,
        py::arg("original_frame"),
        py::arg("new_frame"),
        py::arg("ephemeris_time"),
        "<no doc>");

  m.def("get_body_properties",
        &tudat::spice_interface::getBodyProperties,
        py::arg("body_name"),
        py::arg("property"),
        py::arg("max_n_val"),
        "Get property of a body from Spice. (Wrapper for bodvrd_c routine)");

  m.def("get_body_gravitational_parameter",
        &tudat::spice_interface::getBodyGravitationalParameter,
        py::arg("body_name"),
        "Get gravitational parameter of a body.");

  m.def("get_average_radius",
        &tudat::spice_interface::getAverageRadius,
        py::arg("body_name"),
        "Get the (arithmetic) mean of the three principal axes of the tri-axial ellipsoid shape.");

  m.def("convert_body_name_to_naif_id",
        &tudat::spice_interface::convertBodyNameToNaifId,
        py::arg("body_name"),
        "Convert a body name to its NAIF identification number.");

  // kernel pool related
  m.def("get_standard_kernels",
        &tudat::spice_interface::getStandardSpiceKernels,
        "Get all standard Spice kernels used in tudat.");

  m.def("load_standard_kernels",
        &tudat::spice_interface::loadStandardSpiceKernels,
        py::arg("alternative_kernels") = std::vector<std::string>(),// <pybind11/stl.h>
        tudatpy::load_standard_kernels_docstring().c_str());

  m.def("get_total_count_of_kernels_loaded",
        &tudat::spice_interface::getTotalCountOfKernelsLoaded,
        "Get the amount of loaded Spice kernels.");

  m.def("load_kernel",
        &tudat::spice_interface::loadSpiceKernelInTudat,
        "<no_doc>");

  m.def("clear_kernels",
        &tudat::spice_interface::clearSpiceKernels,
        tudatpy::clear_spice_kernels_docstring().c_str());
};

}// namespace tudatpy
