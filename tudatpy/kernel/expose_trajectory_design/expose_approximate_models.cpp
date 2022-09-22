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
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "expose_approximate_models.h"

#include "tudat/simulation/propagation_setup/propagationCR3BPFullProblem.h"
#include "tudat/astro/basic_astro/clohessyWiltshirePropagator.h"

#include "tudatpy/docstrings.h"

namespace py = pybind11;
//namespace tms = tudat::mission_segments;
namespace tg = tudat::gravitation;
namespace tba = tudat::basic_astrodynamics;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tpc = tudat::physical_constants;
namespace trf = tudat::root_finders;

namespace tudatpy {
namespace trajectory_design {
namespace approximate_models {

void expose_approximate_models(py::module &m)
{


    m.def("propagate_clohessy_wiltshire",
          &tba::propagateClohessyWiltshire,
          py::arg("initial_state"),
          py::arg("propagation_duration"),
          py::arg("central_body_graviational_parameter"),
          py::arg("reference_orbital_radius"),
          get_docstring("propagate_clohessy_wiltshire").c_str() );


    m.def("create_cr3bp_equivalent_body_settings",
          &tp::setupBodySettingsCR3BP,
          py::arg("primary_secondary_distance"),
          py::arg("primary_name"),
          py::arg("secondary_name"),
          py::arg("primary_gravitational_parameter"),
          py::arg("secondary_gravitational_parameter"),
          py::arg("frame_orientation"),
          get_docstring("create_cr3bp_equivalent_body_settings").c_str( ) );



}

} // namespace approximate_models
} // namespace trajectory_design
} // namespace tudatpy
