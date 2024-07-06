/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_spice.h"
#include <tudat/astro/basic_astro.h>

#include "tudatpy/docstrings.h"
#include "tudat/interface/spice.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tsi = tudat::spice_interface;
namespace tba = tudat::basic_astrodynamics;

namespace tudat
{

namespace spice_interface
{

void loadStandardDepracatedSpiceKernels(const std::vector<std::string> alternativeEphemerisKernels)
{
    std::string kernelPath = paths::getSpiceKernelPath();
    loadSpiceKernelInTudat(kernelPath + "/pck00010.tpc");
    loadSpiceKernelInTudat(kernelPath + "/gm_de431.tpc");

    if (alternativeEphemerisKernels.size() == 0) {
        loadSpiceKernelInTudat(kernelPath + "/tudat_merged_spk_kernel.bsp");
    } else {
        for (unsigned int i = 0; i < alternativeEphemerisKernels.size(); i++) {
            loadSpiceKernelInTudat(alternativeEphemerisKernels.at(i));
        }
    }
    loadSpiceKernelInTudat(kernelPath + "/naif0012.tls");
}

}

}

namespace tudatpy {
namespace interface {
namespace spice {

    void expose_spice(py::module &m) {

        // time related
        m.def("convert_julian_date_to_ephemeris_time",
              &tudat::spice_interface::convertJulianDateToEphemerisTime,
              py::arg("julian_date"),
              get_docstring("convert_julian_date_to_ephemeris_time").c_str());

//  m.def("jd2tdb", m.attr("convert_julian_date_to_ephemeris_time"));

        m.def("convert_ephemeris_time_to_julian_date",
              &tudat::spice_interface::convertEphemerisTimeToJulianDate,
              py::arg("ephemeris_time"),
              get_docstring("convert_ephemeris_time_to_julian_date").c_str());

//  m.def("tdb2jd", m.attr("convert_ephemeris_time_to_julian_date"));

        m.def("convert_date_string_to_ephemeris_time",
              &tudat::spice_interface::convertDateStringToEphemerisTime,
              py::arg("date_string"),
              get_docstring("convert_date_string_to_ephemeris_time").c_str());

//  m.def("dstr2jd", m.attr("convert_date_string_to_ephemeris_time"));

        // positional state related
        m.def("get_body_cartesian_position_at_epoch",
              &tudat::spice_interface::getBodyCartesianPositionAtEpoch,
              py::arg("target_body_name"),
              py::arg("observer_body_name"),
              py::arg("reference_frame_name"),
              py::arg("aberration_corrections"),
              py::arg("ephemeris_time"),
              get_docstring("get_body_cartesian_position_at_epoch").c_str());

        m.def("get_body_cartesian_state_at_epoch",
              &tudat::spice_interface::getBodyCartesianStateAtEpoch,
              py::arg("target_body_name"),
              py::arg("observer_body_name"),
              py::arg("reference_frame_name"),
              py::arg("aberration_corrections"),
              py::arg("ephemeris_time"),
              get_docstring("get_body_cartesian_state_at_epoch").c_str());

        m.def("get_cartesian_state_from_tle_at_epoch",
              &tudat::spice_interface::getCartesianStateFromTleAtEpoch,
              py::arg("epoch"),
              py::arg("tle"),
              get_docstring("get_cartesian_state_from_tle_at_epoch").c_str());

        // rotational state related
        m.def("compute_rotation_matrix_between_frames",
              &tudat::spice_interface::computeRotationMatrixBetweenFrames,
              py::arg("original_frame"),
              py::arg("new_frame"),
              py::arg("ephemeris_time"),
              get_docstring("compute_rotation_matrix_between_frames").c_str());

      //   m.def("compute_rotation_quaternion_between_frames",
      //         &tudat::spice_interface::computeRotationQuaternionBetweenFrames,
      //         py::arg("original_frame"),
      //         py::arg("new_frame"),
      //         py::arg("ephemeris_time"),
      //         get_docstring("compute_rotation_quaternion_between_frames").c_str());

        m.def("compute_rotation_matrix_derivative_between_frames",
              &tudat::spice_interface::computeRotationMatrixDerivativeBetweenFrames,
              py::arg("original_frame"),
              py::arg("new_frame"),
              py::arg("ephemeris_time"),
              get_docstring("compute_rotation_matrix_derivative_between_frames").c_str());

        m.def("get_angular_velocity_vector_of_frame_in_original_frame",
              &tudat::spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame,
              py::arg("original_frame"),
              py::arg("new_frame"),
              py::arg("ephemeris_time"),
              get_docstring("get_angular_velocity_vector_of_frame_in_original_frame").c_str());

        m.def("compute_rotation_quaternion_and_rotation_matrix_derivative_between_frames",
              &tudat::spice_interface::computeRotationQuaternionAndRotationMatrixDerivativeBetweenFrames,
              py::arg("original_frame"),
              py::arg("new_frame"),
              py::arg("ephemeris_time"),
              get_docstring("compute_rotation_quaternion_and_rotation_matrix_derivative_between_frames").c_str());

        m.def("get_body_properties",
              &tudat::spice_interface::getBodyProperties,
              py::arg("body_name"),
              py::arg("property"),
              py::arg("max_n_val"),
              get_docstring("get_body_properties").c_str());

        m.def("get_body_gravitational_parameter",
              &tudat::spice_interface::getBodyGravitationalParameter,
              py::arg("body_name"),
              get_docstring("get_body_gravitational_parameter").c_str());

        m.def("get_average_radius",
              &tudat::spice_interface::getAverageRadius,
              py::arg("body_name"),
              get_docstring("get_average_radius").c_str());

        m.def("convert_body_name_to_naif_id",
              &tudat::spice_interface::convertBodyNameToNaifId,
              py::arg("body_name"),
              get_docstring("convert_body_name_to_naif_id").c_str());

//        // kernel pool related
//        m.def("get_standard_kernels",
//              &tudat::spice_interface::getStandardSpiceKernels,
//              get_docstring("get_standard_kernels").c_str());

        m.def("load_standard_kernels",
              &tudat::spice_interface::loadStandardSpiceKernels,
              py::arg("alternative_kernels") = std::vector<std::string>(),// <pybind11/stl.h>
              get_docstring("load_standard_kernels").c_str());

        m.def("load_standard_deprecated_kernels",
              &tudat::spice_interface::loadStandardDepracatedSpiceKernels,
              py::arg("alternative_kernels") = std::vector<std::string>(),// <pybind11/stl.h>
              get_docstring("load_standard_kernels").c_str());

        m.def("get_total_count_of_kernels_loaded",
              &tudat::spice_interface::getTotalCountOfKernelsLoaded,
              get_docstring("get_total_count_of_kernels_loaded").c_str());

        m.def("check_body_property_in_kernel_pool",
              &tudat::spice_interface::checkBodyPropertyInKernelPool,
              py::arg("body_name"),
              py::arg("body_property"),
              get_docstring("check_body_property_in_kernel_pool").c_str());

        m.def("load_kernel",
              &tudat::spice_interface::loadSpiceKernelInTudat,
              py::arg("kernel_file"),
              get_docstring("load_kernel").c_str());

        m.def("clear_kernels",
              &tudat::spice_interface::clearSpiceKernels,
              get_docstring("clear_kernels").c_str());

//      py::class_<tudat::ephemerides::SpiceEphemeris,
//            std::shared_ptr<tudat::ephemerides::SpiceEphemeris>>(m, "SpiceEphemeris",
//                                    get_docstring("SpiceEphemeris").c_str())
//            .def(py::init<
//                  const std::string &,
//                  const std::string &,
//                  const bool,
//                  const bool,
//                  const bool,
//                  const std::string &,
//                  const double>(),
//                  py::arg("target_body_name"),
//                  py::arg("observer_body_name"),
//                  py::arg("correct_for_stellar_aberration") = false,
//                  py::arg("correct_for_light_time_aberration") = true,
//                  py::arg("converge_light_time_aberration") = false,
//                  py::arg("reference_frame_name") = "ECLIPJ2000",
//                  py::arg("reference_julian_day") = tba::JULIAN_DAY_ON_J2000,
//                  get_docstring("SpiceEphemeris.ctor").c_str())
//            .def("get_cartesian_state",
//                  &tudat::ephemerides::SpiceEphemeris::getCartesianState,
//                  py::arg("seconds_since_epoch"),
//                  get_docstring("SpiceEphemeris.get_cartesian_state").c_str());

    };

}// namespace spice
}// namespace interface
}// namespace tudatpy
