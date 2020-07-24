/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_ephemerides.h"

#include <tudat/astro/ephemerides.h>
#include <tudat/simulation/simulation.h>// TODO: EphemerisType should be in <tudat/astro/ephemerides.h>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace te = tudat::ephemerides;
namespace tss = tudat::simulation_setup;

namespace tudatpy {
void expose_ephemerides(py::module &m) {
  py::class_<te::Ephemeris, std::shared_ptr<te::Ephemeris>> ephemeris(m, "Ephemeris");

  py::enum_<tss::EphemerisType>(ephemeris, "EphemerisType")
      .value("approximate_planet_positions", tss::approximate_planet_positions)
      .value("direct_spice_ephemeris", tss::direct_spice_ephemeris)
      .value("interpolated_spice", tss::interpolated_spice)
      .value("constant_ephemeris", tss::constant_ephemeris)
      .value("kepler_ephemeris", tss::kepler_ephemeris)
      .value("custom_ephemeris", tss::custom_ephemeris);

  py::class_<te::RotationalEphemeris,
             std::shared_ptr<te::RotationalEphemeris>>
      RotationalEphemeris_(m, "RotationalEphemeris");

  m.def("transform_state_to_global_frame",
        &te::transformStateToGlobalFrame<double, double>,
        py::arg("state_in_local_frame"),
        py::arg("current_time"),
        py::arg("rotational_ephemeris"));

  py::class_<te::ConstantEphemeris,
             std::shared_ptr<te::ConstantEphemeris>,
             te::Ephemeris>(
      m, "ConstantEphemeris")
      .def(py::init<
               const std::function<Eigen::Vector6d()>,//<pybind11/functional.h>,<pybind11/eigen.h>
               const std::string &,
               const std::string &>(),
           py::arg("constant_state_function"),
           py::arg("reference_frame_origin") = "SSB",
           py::arg("reference_frame_orientation") = "ECLIPJ2000")
      .def(py::init<
               const Eigen::Vector6d,//<pybind11/eigen.h>
               const std::string &,
               const std::string &>(),
           py::arg("constant_state"),
           py::arg("reference_frame_origin") = "SSB",
           py::arg("reference_frame_orientation") = "ECLIPJ2000")
      .def("get_cartesian_state", &te::ConstantEphemeris::getCartesianState,
           py::arg("seconds_since_epoch") = 0.0)
      .def("update_constant_state", &te::ConstantEphemeris::updateConstantState,
           py::arg("new_state"))
      .def("get_cartesian_position", &te::ConstantEphemeris::getCartesianPosition,
           py::arg("seconds_since_epoch"))
      .def("get_cartesian_velocity", &te::ConstantEphemeris::getCartesianVelocity,
           py::arg("seconds_since_epoch"))
      .def("get_cartesian_long_state", &te::ConstantEphemeris::getCartesianLongState,
           py::arg("seconds_since_epoch"))
      .def("get_cartesian_state_from_extended_time",
           &te::ConstantEphemeris::getCartesianStateFromExtendedTime,
           py::arg("current_time"))
      .def("get_cartesian_long_state_from_extended_time",
           &te::ConstantEphemeris::getCartesianLongStateFromExtendedTime,
           py::arg("current_time"));
};
}// namespace tudatpy