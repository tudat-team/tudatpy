/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_observable_models_setup.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"

namespace py = pybind11;
// namespace tss = tudat::simulation_setup;
// namespace tp = tudat::propagators;
// namespace tom = tudat::observation_models;
// namespace tep = tudat::estimatable_parameters;

namespace tudatpy
{
namespace estimation_refactoring
{
namespace observable_models_setup
{

void expose_observable_models_setup( py::module& m )
{

    // ************** Modules ***************
    auto biases = m.def_submodule( "biases" );
    biases::expose_biases( biases );

}

}  // namespace observable_models_setup
}  // namespace estimation_refactoring
}  // namespace tudatpy
