/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_vehicle_systems_setup.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/environment_setup.h>

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
namespace numerical_simulation {
namespace environment_setup {
namespace vehicle_systems {

    void expose_vehicle_systems_setup(py::module &m)
    {

        py::class_<tss::BodyPanelGeometrySettings,
            std::shared_ptr<tss::BodyPanelGeometrySettings>>(m, "BodyPanelGeometrySettings",
                                                     get_docstring("BodyPanelGeometrySettings").c_str());

        m.def("frame_fixed_panel_geometry",
              &tss::frameFixedPanelGeometry,
              py::arg("surface_normal"),
              py::arg("area"),
              py::arg("frame_orientation"),
              get_docstring("frame_fixed_geometry").c_str());

        m.def("body_tracking_panel_geometry",
              &tss::bodyTrackingPanelGeometry,
              py::arg("body_to_track"),
              py::arg("towards_tracked_body"),
              py::arg("area"),
              py::arg("frame_orientation"),
              get_docstring("frame_fixed_geometry").c_str());

        m.def("time_varying_panel_geometry",
              &tss::timeVaryingPanelGeometry,
              py::arg("surface_normal_function"),
              py::arg("area"),
              py::arg("frame_orientation"),
              get_docstring("frame_fixed_geometry").c_str());

        py::class_<tss::BodyPanelSettings,
            std::shared_ptr<tss::BodyPanelSettings>>(m, "BodyPanelSettings",
                                                     get_docstring("BodyPanelSettings").c_str())
            .def_readwrite("reflection_law_settings", &tss::BodyPanelSettings::reflectionLawSettings_ );

        m.def("body_panel_settings",
              &tss::bodyPanelSettings,
              py::arg("panel_geometry"),
              py::arg("panel_reflection_law") = nullptr,
              py::arg("panel_type_id") = "",
              get_docstring("frame_fixed_geometry").c_str());

        py::class_<tss::FullPanelledBodySettings,
            std::shared_ptr<tss::FullPanelledBodySettings>>(m, "FullPanelledBodySettings",
                                                     get_docstring("FullPanelledBodySettings").c_str());


        m.def("full_panelled_body_settings",
              &tss::fullPanelledBodySettings,
              py::arg("panel_settings"),
              py::arg("part_rotation_model_settings") = std::map< std::string, std::shared_ptr< tss::RotationModelSettings > >( ),
              get_docstring("frame_fixed_geometry").c_str() );
    }

}// namespace vehicle_systems
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
