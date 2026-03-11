#ifndef TUDATPY_DYNAMICS_PROPAGATION_EXPOSE_PROPAGATION_BINDINGS_H
#define TUDATPY_DYNAMICS_PROPAGATION_EXPOSE_PROPAGATION_BINDINGS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace propagation
{

void expose_propagation_state_utility_bindings( py::module& m );
void expose_propagation_results_bindings( py::module& m );
void expose_propagation_thrust_bindings( py::module& m );

}  // namespace propagation
}  // namespace dynamics
}  // namespace tudatpy

#endif
