#ifndef TUDATPY_ESTIMATION_OBSERVATIONS_BINDINGS_H
#define TUDATPY_ESTIMATION_OBSERVATIONS_BINDINGS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace observations
{

void expose_observations_io_bindings( py::module& m );
void expose_observations_simulation_bindings( py::module& m );

}  // namespace observations
}  // namespace estimation
}  // namespace tudatpy

#endif
