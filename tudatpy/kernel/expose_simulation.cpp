/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_simulation.h"

#include "expose_simulation/expose_environment_setup.h"
#include "expose_simulation/expose_estimation_setup.h"
#include "expose_simulation/expose_propagation_setup.h"
#include "expose_simulation/expose_astro_setup.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy {
namespace simulation {

void expose_simulation(py::module &m) {

  auto environment_setup_submodule = m.def_submodule("environment_setup");
  environment_setup::expose_environment_setup(environment_setup_submodule);

  auto propagation_setup_submodule = m.def_submodule("propagation_setup");
  propagation_setup::expose_propagation_setup(propagation_setup_submodule);

  auto estimation_setup_submodule = m.def_submodule("estimation_setup");
  estimation_setup::expose_estimation_setup(estimation_setup_submodule);

  auto astro_setup = m.def_submodule("astro_setup");
  astro_setup::expose_astro_setup(astro_setup);

};

}// namespace simulation
}// namespace tudatpy
