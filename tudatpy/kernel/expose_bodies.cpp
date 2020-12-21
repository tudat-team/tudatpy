/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_bodies.h"
#include "expose_bodies/expose_default.h"
#include "expose_bodies/prototype/simpleBody.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace tudat::bodies;

namespace tudatpy {

void expose_bodies(py::module &m) {

  py::class_<
      SimpleBody,
      std::shared_ptr<SimpleBody>>(m, "SimpleBody")
      .def_property_readonly("parent", &SimpleBody::getParent)
      .def_property_readonly("children", &SimpleBody::getChildren)
      .def_property_readonly("gravitational_parameter", &SimpleBody::getGravitationalParameter)
      .def_property_readonly("name", &SimpleBody::getName)
      .def_property_readonly("symbol", &SimpleBody::getSymbol);

  py::class_<
      SimpleSystemOfBodies,
      std::shared_ptr<SimpleSystemOfBodies>>(m, "SimpleSystemOfBodies")
      .def("__getattr__", [=](std::shared_ptr<SimpleSystemOfBodies> &a, std::string &n) { return a->getBodyMap()[n]; })
      .def("body_map", &SimpleSystemOfBodies::getBodyMap)
      .def_property_readonly("gravitational_parameter", &SimpleSystemOfBodies::getGravitationalParameter);


  auto _default = m.def_submodule("default");
  expose_default(_default);
};

}// namespace tudatpy