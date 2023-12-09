/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_math/expose_interpolators.h"
#include "expose_math/expose_numerical_integrators.h"
#include "expose_math/expose_root_finders.h"
#include "expose_math/expose_geometry.h"
#include "expose_math/expose_statistics.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace tudatpy {
namespace math {

void expose_math(py::module &m) {

  auto interpolators = m.def_submodule("interpolators");
  interpolators::expose_interpolators(interpolators);

  auto numerical_integrators = m.def_submodule("numerical_integrators");
  expose_numerical_integrators(numerical_integrators);

  auto root_finders = m.def_submodule("root_finders");
  expose_root_finders(root_finders);

  auto geometry = m.def_submodule("geometry");
  expose_geometry(geometry);

  auto statistics = m.def_submodule("statistics");
  expose_statistics(statistics);

}
};

};// namespace tudatpy
