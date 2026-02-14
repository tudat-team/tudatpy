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
#include "expose_propagation.h"

#include <tudat/astro/aerodynamics/aerodynamicGuidance.h>
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/propagators.h>

#include <tudat/io/serialization/base.h>

#include "scalarTypes.h"
#include "expose_propagation_bindings.h"

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace propagation
{

void expose_propagation( py::module& m )
{
    expose_propagation_state_utility_bindings( m );
    expose_propagation_results_bindings( m );
    expose_propagation_thrust_bindings( m );
}

}  // namespace propagation
}  // namespace dynamics
}  // namespace tudatpy
