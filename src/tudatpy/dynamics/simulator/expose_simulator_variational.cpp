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
#include "expose_simulator_bindings.h"

#include <pybind11/pybind11.h>

#include "expose_simulator_variational_helpers.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "kernelExternTemplates.h"

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

void expose_simulator_variational_bindings( py::module& m )
{
    py::class_< VariationalSolverBaseType, std::shared_ptr< VariationalSolverBaseType > >(
            m,
            "VariationalSimulator",
            R"doc(Base class for variational equations propagation.)doc" )
            .def_property_readonly(
                    "state_transition_interface",
                    &VariationalSolverBaseType::getStateTransitionMatrixInterface,
                    R"doc(State transition/sensitivity interface for covariance propagation.)doc" );

    SingleArcVariationalSimulatorPyClass singleArcClass(
            m,
            "SingleArcVariationalSimulator",
            R"doc(Class for single-arc variational equations propagation.)doc" );

    bind_single_arc_variational_integrate_methods( singleArcClass );
    bind_single_arc_variational_result_properties( singleArcClass );
    bind_create_variational_equations_solver_function( m );
}

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy
