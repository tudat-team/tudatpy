/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_astro_setup.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <tudat/astro/mission_segments/createTransferTrajectory.h>
#include <tudat/simulation/propagation_setup/accelerationSettings.h>


namespace py = pybind11;
namespace tms = tudat::mission_segments;
namespace tss = tudat::simulation_setup;


namespace tudatpy {
namespace simulation {
namespace astro_setup {

void expose_astro_setup(py::module &m)
{

//    m.def("defaut_mga_no_dsm_settings",
//          &tms::createTransferTrajectory,
//          py::arg( "bodies" ),
//          py::arg( "leg_settings" ),
//          py::arg( "node_settings" ),
//          py::arg( "node_names" ),
//          py::arg( "central_body" ),
//          py::arg( "node_times" ),
//          py::arg( "leg_free_parameters" ),
//          py::arg( "node_free_parameters" ) );


    m.def("create_transfer_trajectory",
          &tms::createTransferTrajectory,
          py::arg( "bodies" ),
          py::arg( "leg_settings" ),
          py::arg( "node_settings" ),
          py::arg( "node_names" ),
          py::arg( "central_body" ) );


    m.def("get_low_thrust_acceleration_settings",
          &tss::getLowThrustLegAccelerationSettings,
          py::arg("low_thrust_leg"),
          py::arg("bodies"),
          py::arg("body_to_propagate"),
          py::arg("specific_impulse_function"),
          py::arg("low_thrust_leg_initial_time") );


}

}
}
}// namespace tudatpy
