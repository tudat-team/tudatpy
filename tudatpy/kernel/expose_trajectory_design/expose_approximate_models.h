/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_APPROXIMATE_MODELS_H
#define TUDATPY_EXPOSE_APPROXIMATE_MODELS_H

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace tudatpy {
namespace trajectory_design {
namespace approximate_models {

    void expose_approximate_models(py::module &m);

} // namespace transfer_trajectory
} // namespace trajectory_design
} // namespace tudatpy

#endif// TUDATPY_EXPOSE_APPROXIMATE_MODELS_H
