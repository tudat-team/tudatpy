/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_frames.h"
#include "prototype/frames.h"

#include <pybind11/pybind11.h>

#include <tudat/astro/reference_frames.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

namespace trf = tudat::reference_frames;

namespace py = pybind11;

namespace tudatpy {

void expose_frames(py::module &m) {
  //////////////////////////////////////////////////////////////////////
  // prototype [submodule] / frames
  //////////////////////////////////////////////////////////////////////
  py::class_<ReferenceFrame>(m, "ReferenceFrame")
      .def(py::init<std::string, std::string>(),
           py::arg("origin"),
           py::arg("orientation"))
      .def(py::init<>())
      .def("__str__", &ReferenceFrame::getString)
      .def_property_readonly("origin", &ReferenceFrame::getOrigin)
      .def_property_readonly("orientation", &ReferenceFrame::getOrientation);

  m.attr("SSB_J2000") = SSB_J2000();
  m.attr("SSB_ECLIPJ2000") = SSB_ECLIPJ2000();


  py::class_<trf::AerodynamicAngleCalculator,
      std::shared_ptr<trf::AerodynamicAngleCalculator>>(m, "AerodynamicAngleCalculator")
      .def("set_orientation_angle_functions",
           py::overload_cast<
               const std::function<double()>,
               const std::function<double()>,
               const std::function<double()>,
               const std::function<void(const double)>>(&trf::AerodynamicAngleCalculator::setOrientationAngleFunctions),
           py::arg("angle_of_attack_function") = std::function<double()>(),       // <pybind11/functional.h>
           py::arg("angle_of_sideslip_function") = std::function<double()>(),     // <pybind11/functional.h>
           py::arg("bank_angle_function") = std::function<double()>(),            // <pybind11/functional.h>
           py::arg("angle_update_function") = std::function<void(const double)>(),// <pybind11/functional.h>
           "<no_doc>")
      .def("set_orientation_angle_functions",
           py::overload_cast<
               const double,
               const double,
               const double>(&trf::AerodynamicAngleCalculator::setOrientationAngleFunctions),
           py::arg("angle_of_attack") = TUDAT_NAN,
           py::arg("angle_of_sideslip") = TUDAT_NAN,
           py::arg("bank_angle") = TUDAT_NAN,
           "<no_doc>");

}

}// namespace tudatpy
