/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include "expose_dynamics.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#include "environment/expose_environment.h"
#include "environment_setup/expose_environment_setup.h"
#include "propagation/expose_propagation.h"
#include "propagation_setup/expose_propagation_setup.h"
#include "simulator/expose_simulator.h"
#include "parameters_setup/expose_parameters_setup.h"
#include "parameters/expose_parameters.h"
#include "scalarTypes.h"

namespace py = pybind11;

namespace tudatpy
{

namespace dynamics
{

void expose_dynamics_types( py::module& m )
{
    auto environment_submodule = m.def_submodule( "environment" );
    auto environment_setup_submodule = m.def_submodule( "environment_setup" );
    auto propagation_submodule = m.def_submodule( "propagation" );
    auto parameters_setup_submodule = m.def_submodule( "parameters_setup" );

    environment_setup::expose_environment_setup_types( environment_setup_submodule );
    propagation::expose_propagation_types( propagation_submodule );
    parameters_setup::expose_parameters_setup_types( parameters_setup_submodule );
    environment::expose_environment( environment_submodule );
}

void expose_dynamics( py::module& m )
{
    auto environment_setup_submodule = py::module_::import( "tudatpy.kernel.dynamics.environment_setup" );
    auto propagation_submodule = py::module_::import( "tudatpy.kernel.dynamics.propagation" );
    auto parameters_setup_submodule = py::module_::import( "tudatpy.kernel.dynamics.parameters_setup" );

    auto propagation_setup_submodule = m.def_submodule( "propagation_setup" );
    propagation_setup::expose_propagation_setup( propagation_setup_submodule );

    environment_setup::expose_environment_setup( environment_setup_submodule );
    propagation::expose_propagation( propagation_submodule );

    parameters_setup::expose_parameters_setup( parameters_setup_submodule );

    auto parameters_submodule = m.def_submodule( "parameters" );
    parameters::expose_parameters( parameters_submodule );

    auto simulator_submodule = m.def_submodule( "simulator" );
    simulator::expose_simulator( simulator_submodule );
};

}  // namespace dynamics
}  // namespace tudatpy
