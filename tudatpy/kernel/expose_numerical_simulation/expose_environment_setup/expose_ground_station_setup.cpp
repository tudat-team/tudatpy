/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_ground_station_setup.h"

#include "docstrings.h"
#include <tudat/simulation/environment_setup/createGroundStations.h>
#include <tudat/simulation/environment_setup/defaultBodies.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tcc = tudat::coordinate_conversions;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy {
namespace numerical_simulation {
namespace environment_setup {
namespace ground_station{

    void expose_ground_station_setup(py::module &m) {

        py::class_<tss::GroundStationSettings,
                std::shared_ptr<tss::GroundStationSettings>>(m, "GroundStationSettings",
                                                         get_docstring("GroundStationSettings").c_str());

        py::class_<tss::GroundStationMotionSettings,
                std::shared_ptr<tss::GroundStationMotionSettings>>(m, "GroundStationMotionSettings",
                                        get_docstring("GroundStationMotionSettings").c_str());

        py::class_<tss::LinearGroundStationMotionSettings,
                std::shared_ptr<tss::LinearGroundStationMotionSettings>,
                tss::GroundStationMotionSettings>(m, "LinearGroundStationMotionSettings",
                                        get_docstring("LinearGroundStationMotionSettings").c_str());

        py::class_<tss::PiecewiseConstantGroundStationMotionSettings,
                std::shared_ptr<tss::PiecewiseConstantGroundStationMotionSettings>,
                tss::GroundStationMotionSettings>(m, "PiecewiseConstantGroundStationMotionSettings",
                                        get_docstring("PiecewiseConstantGroundStationMotionSettings").c_str());

        py::class_<tss::CustomGroundStationMotionSettings,
                std::shared_ptr<tss::CustomGroundStationMotionSettings>,
                tss::GroundStationMotionSettings>(m, "CustomGroundStationMotionSettings",
                                        get_docstring("CustomGroundStationMotionSettings").c_str());


        m.def("add_motion_model_to_each_groun_station",
              &tss::addStationMotionModelToEachGroundStation,
              py::arg("ground_station_settings_list"),
              py::arg("station_motion_setting") );

        m.def("basic_station",
              &tss::groundStationSettings,
              py::arg("station_name"),
              py::arg("station_nominal_position"),
              py::arg("station_position_element_type") = tcc::cartesian_position,
              py::arg("station_motion_settings") = std::vector< std::shared_ptr< tss::GroundStationMotionSettings    > >( ),
              get_docstring("basic_station").c_str());

        m.def("dsn_stations",
              &tss::getDsnStationSettings,
              get_docstring("dsn_stations").c_str());

        m.def("linear_station_motion",
              &tss::linearGroundStationMotionSettings,
              py::arg("linear_velocity"),
              py::arg("reference_epoch") = 0.0,
              get_docstring("linear_station_motion").c_str());

        m.def("piecewise_constant_station_motion",
              &tss::piecewiseConstantGroundStationMotionSettings,
              py::arg("displacement_list"),
              get_docstring("piecewise_constant_station_motion").c_str());


        m.def("custom_station_motion",
              &tss::customGroundStationMotionSettings,
              py::arg("custom_displacement_function"),
              get_docstring("custom_station_motion").c_str());

    }

}// namespace ground_station
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
