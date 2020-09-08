/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_interpolators.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace ti = tudat::interpolators;

namespace tudatpy {

void expose_interpolators(py::module &m) {

  py::class_<
      ti::LagrangeInterpolatorSettings,
      std::shared_ptr<ti::LagrangeInterpolatorSettings>>(m,
                                                         "LagrangeInterpolatorSettings",
                                                         "<no doc>");
};

}// namespace tudatpy
