/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_spice_interface.h"

#include "../docstrings.h"
#include "tudat/interface/spice.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tsi = tudat::spice_interface;

namespace tudatpy {

void expose_spice_interface(py::module &m) {

  m.def("load_standard_kernels",
        &tudat::spice_interface::loadStandardSpiceKernels,
        py::arg("alternative_kernels") = std::vector<std::string>(),// <pybind11/stl.h>
        tudatpy::load_standard_kernels_docstring().c_str());

  m.def("load_kernel",
        &tudat::spice_interface::loadSpiceKernelInTudat,
        "<no_doc>");

  m.def("clear_kernels",
        &tudat::spice_interface::clearSpiceKernels,
        tudatpy::clear_spice_kernels_docstring().c_str());

  m.def("get_average_radius",
        &tudat::spice_interface::getAverageRadius,
        "<no_doc>");

};

}// namespace tudatpy
