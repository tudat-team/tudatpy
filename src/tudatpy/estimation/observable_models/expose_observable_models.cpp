/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_observable_models.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "observables_simulation/expose_observables_simulation.h"

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace observable_models
{

void expose_observable_models( py::module& m )
{
    auto observables_simulation = m.def_submodule( "observables_simulation" );
    observables_simulation::expose_observables_simulation( observables_simulation );
}

}  // namespace observable_models
}  // namespace estimation
}  // namespace tudatpy
