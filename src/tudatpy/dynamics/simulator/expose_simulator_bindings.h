#ifndef TUDATPY_DYNAMICS_SIMULATOR_EXPOSE_SIMULATOR_BINDINGS_H
#define TUDATPY_DYNAMICS_SIMULATOR_EXPOSE_SIMULATOR_BINDINGS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

void expose_simulator_dynamics_bindings( py::module& m );
void expose_simulator_variational_bindings( py::module& m );
void expose_simulator_state_transition_bindings( py::module& m );

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy

#endif
