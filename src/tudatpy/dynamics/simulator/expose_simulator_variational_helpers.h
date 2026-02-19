#ifndef TUDATPY_DYNAMICS_SIMULATOR_EXPOSE_SIMULATOR_VARIATIONAL_HELPERS_H
#define TUDATPY_DYNAMICS_SIMULATOR_EXPOSE_SIMULATOR_VARIATIONAL_HELPERS_H

#include <memory>

#include <pybind11/pybind11.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolverForwardDeclarations.h"

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

using VariationalSolverBaseType = tudat::propagators::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE, 0 >;
using SingleArcVariationalSolverType = tudat::propagators::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;

using SingleArcVariationalSimulatorPyClass = py::class_<
        SingleArcVariationalSolverType,
        std::shared_ptr< SingleArcVariationalSolverType >,
        VariationalSolverBaseType >;

void bind_single_arc_variational_integrate_methods( SingleArcVariationalSimulatorPyClass& singleArcClass );
void bind_single_arc_variational_result_properties( SingleArcVariationalSimulatorPyClass& singleArcClass );
void bind_create_variational_equations_solver_function( py::module& m );

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy

#endif
