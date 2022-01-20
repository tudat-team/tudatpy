/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_interface/expose_spice.h"
//#include "interface/expose_json_interface.h"
//#include "interface/expose_sofa_interface.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy {
namespace interface {
using namespace spice;

void expose_interface(py::module &m) {

  auto spice = m.def_submodule("spice");
  expose_spice(spice);
  m.attr("spice_interface") = m.attr("spice");

  //  auto json_interface = m.def_submodule("json_interface");
  //  expose_json_interface(sofa_interface);

  //  auto sofa_interface = m.def_submodule("sofa_interface");
  //  expose_sofa_interface(sofa_interface);
}

}// namespace interface
}// namespace tudatpy
