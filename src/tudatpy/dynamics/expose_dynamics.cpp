/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_dynamics.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"

namespace py = pybind11;

namespace tudatpy
{

namespace dynamics
{

void expose_dynamics( py::module &m )
{
    
    auto environment_setup_submodule = m.def_submodule( "environment_setup" );
    environment_setup::expose_environment_setup( environment_setup_submodule );

    auto environment_submodule = m.def_submodule( "environment" );
    environment::expose_environment( environment_submodule );

    auto propagation_setup_submodule = m.def_submodule( "propagation_setup" );
    propagation_setup::expose_propagation_setup( propagation_setup_submodule );

    auto propagation_submodule = m.def_submodule( "propagation" );
    propagation::expose_propagation( propagation_submodule );

    auto simulator_submodule = m.def_submodule( "simulator" );
    simulator::expose_simulator( simulator_submodule );

};

}  // namespace dynamics
}  // namespace tudatpy
