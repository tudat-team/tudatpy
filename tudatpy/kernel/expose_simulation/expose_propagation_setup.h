/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_PROPAGATION_SETUP_H
#define TUDATPY_EXPOSE_PROPAGATION_SETUP_H

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/simulation/propagation_setup.h"

#include "../docstrings.h"

namespace py = pybind11;

namespace tp = tudat::propagators;

namespace tudatpy {

void expose_propagation_setup(py::module &m);


}// namespace tudatpy

#endif// TUDATPY_EXPOSE_PROPAGATION_SETUP_H
