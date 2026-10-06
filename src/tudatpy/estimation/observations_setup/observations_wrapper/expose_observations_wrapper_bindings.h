#ifndef TUDATPY_ESTIMATION_OBSERVATIONS_SETUP_OBSERVATIONS_WRAPPER_BINDINGS_H
#define TUDATPY_ESTIMATION_OBSERVATIONS_SETUP_OBSERVATIONS_WRAPPER_BINDINGS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{
namespace observations_wrapper
{

void expose_observations_wrapper_io_bindings( py::module& m );
void expose_observations_wrapper_simulation_bindings( py::module& m );
void expose_observations_wrapper_sum_lmk_bindings( py::module& m );

}  // namespace observations_wrapper
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy

#endif
