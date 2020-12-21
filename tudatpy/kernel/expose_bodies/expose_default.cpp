/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_default.h"

// TODO(Geoffrey): This needs to be modified accordingly in the tudat source.
//#include <tudat/simulation/environment_setup/defaultBodies.h>
#include "prototype/default.h"

#include <pybind11/pybind11.h>

namespace tb = tudat::bodies;

namespace tbss = tudat::bodies::solar_system;

namespace py = pybind11;

namespace tudatpy {

void expose_default(py::module &m) {

  // solar_system
  auto solar_system = m.def_submodule("solar_system");
  solar_system.attr("SolarSystem") = tbss::SolarSystem;
  solar_system.attr("Sun") = tbss::Sun;
  solar_system.attr("Mercury") = tbss::Mercury;
//  solar_system.attr("Venus") = tbss::Venus;
//  solar_system.attr("Earth") = tbss::Earth;
//  solar_system.attr("Mars") = tbss::Mars;
//  solar_system.attr("name_map") = tbss::nameMap;

  //    solar_system.attr("Mercury", tbss::Mercury());
  //    solar_system.attr("Venus", tbss::Venus());
  //    solar_system.attr("Earth", tbss::Earth());
  //    solar_system.attr("Mars", tbss::Mars());
  //    solar_system.attr("Jupiter", tbss::Jupiter());
  //    solar_system.attr("Saturn", tbss::Saturn());
  //    solar_system.attr("Uranus", tbss::Uranus());
  //    solar_system.attr("Neptune", tbss::Neptune());
  //    solar_system.attr("Pluto", tbss::Pluto());

  //  auto earth_system = py::module("earth_system");

  // martian_system
  // jovian_system
  // saturnian_system
  // neptunian_system
  // uranium_system
}

}// namespace tudatpy
