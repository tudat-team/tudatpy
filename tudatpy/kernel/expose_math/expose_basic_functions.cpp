/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_basic_math.h"

#include <tudat/basics/basicTypedefs.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

namespace tbm = tudat::basic_mathematics;

namespace tudatpy {

void expose_basic_functions(py::module &m) {

    m.def("legendre_normalization_factor",
          &tbm::calculateLegendreGeodesyNormalizationFactor,
          py::arg("degree"),
          py::arg("order") );


};

}// namespace tudatpy
