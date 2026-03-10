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
#include "expose_simulator.h"

#include "expose_simulator_bindings.h"

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

void expose_simulator( py::module& m )
{
    expose_simulator_dynamics_bindings( m );
    expose_simulator_variational_bindings( m );
    expose_simulator_state_transition_bindings( m );
}

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy
