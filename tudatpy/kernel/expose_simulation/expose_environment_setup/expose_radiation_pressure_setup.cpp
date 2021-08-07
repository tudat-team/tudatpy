/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_radiation_pressure_setup.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/environment_setup.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

//#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
//#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy {
namespace simulation {
namespace environment_setup {


    void expose_radiation_pressure_setup(py::module &m) {

        /////////////////////////////////////////////////////////////////////////////
        // createRadiationPressureInterface.h
        /////////////////////////////////////////////////////////////////////////////
        py::enum_<tss::RadiationPressureType>(m, "RadiationPressureType", "<no_doc>")
                .value(
                        "cannonball_radiation_pressure_interface",
                        tss::RadiationPressureType::cannon_ball_radiation_pressure_interface)
                .value("panelled_radiation_pressure_interface",
                       tss::RadiationPressureType::panelled_radiation_pressure_interface)
                .value("solar_sailing_radiation_pressure_interface",
                       tss::RadiationPressureType::
                       solar_sailing_radiation_pressure_interface)
                .export_values();


        py::class_<tss::RadiationPressureInterfaceSettings,
                std::shared_ptr<tss::RadiationPressureInterfaceSettings>>(
                m, "RadiationPressureInterfaceSettings", "<no_doc>");
//            .def(py::init<const tss::RadiationPressureType, const std::string &,
//                 const std::vector<std::string>>(),
//                 py::arg("radiation_pressure_type"), py::arg("source_body"),
//                 py::arg("occulting_bodies") = std::vector<std::string>());

        py::class_<tss::CannonBallRadiationPressureInterfaceSettings,
                std::shared_ptr<tss::CannonBallRadiationPressureInterfaceSettings>,
                tss::RadiationPressureInterfaceSettings>(
                m, "CannonBallRadiationPressureInterfaceSettings", "<no_doc>");
//            .def(py::init<const std::string &, const double, const double,
//                 const std::vector<std::string> &>(),
//                 py::arg("source_body"), py::arg("area"),
//                 py::arg("radiation_pressure_coefficient"),
//                 py::arg("occulting_bodies") = std::vector<std::string>())

        m.def("cannonball",
              py::overload_cast<const std::string &, const double, const double, const std::vector<std::string> &>(
                      &tss::cannonBallRadiationPressureSettings),
              py::arg("source_body"), py::arg("reference_area"),
              py::arg("radiation_pressure_coefficient"),
              py::arg("occulting_bodies") = std::vector<std::string>());

        m.def("panelled",
              &tss::panelledRadiationPressureInterfaceSettings,
              py::arg("source_body"),
              py::arg("emissivities"),
              py::arg("areas"),
              py::arg("diffusion_coefficients"),
              py::arg("surface_normals_in_body_fixed_frame"),
              py::arg("occulting_bodies") = std::vector<std::string>());

    }


}// namespace environment_setup
}// namespace simulation
}// namespace tudatpy