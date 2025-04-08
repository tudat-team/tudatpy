/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_NUMERICAL_SIMULATION_H
#define TUDATPY_EXPOSE_NUMERICAL_SIMULATION_H

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "environment/environment.hpp"
#include "environment_setup/environment_setup.hpp"
#include "estimation/estimation.hpp"
#include "estimation_setup/estimation_setup.hpp"
#include "propagation/propagation.hpp"
#include "propagation_setup/propagation_setup.hpp"

namespace py = pybind11;

namespace tudatpy
{

namespace numerical_simulation
{

void expose_numerical_simulation( py::module &m );
void expose_numerical_simulation_simulator( py::module &m );
void expose_numerical_simulation_estimator( py::module &m );
void expose_numerical_simulation_variational( py::module &m );

}  // namespace numerical_simulation

}  // namespace tudatpy
#endif  // TUDATPY_EXPOSE_NUMERICAL_SIMULATION_H
