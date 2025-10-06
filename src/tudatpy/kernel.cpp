#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/pybind11.h>
#include <tudat/config.hpp>

namespace py = pybind11;

// Declarations from each binding .cpp file
void add_math_to_kernel(py::module_& m);
void add_astro_to_kernel(py::module_& m);
void add_trajectory_design_to_kernel(py::module_& m);
void add_constants_to_kernel(py::module_& m);
void add_interface_to_kernel(py::module_& m);
void add_data_to_kernel(py::module_& m);
void add_dynamics_to_kernel(py::module_& m);
void add_estimation_to_kernel(py::module_& m);
void add_exceptions_to_kernel(py::module_& m);

PYBIND11_MODULE(kernel, m)
{
    auto math = m.def_submodule("math");
    add_math_to_kernel(math);

    auto astro = m.def_submodule("astro");
    add_astro_to_kernel(astro);

    auto trajectory_design = m.def_submodule("trajectory_design");
    add_trajectory_design_to_kernel(trajectory_design);

    auto constants = m.def_submodule("constants");
    add_constants_to_kernel(constants);

    auto interface = m.def_submodule("interface");
    add_interface_to_kernel(interface);

    auto data = m.def_submodule("data");
    add_data_to_kernel(data);

    auto dynamics = m.def_submodule("dynamics");
    add_dynamics_to_kernel(dynamics);

    auto estimation = m.def_submodule("estimation");
    add_estimation_to_kernel(estimation);

    auto exceptions = m.def_submodule("exceptions");
    add_exceptions_to_kernel(exceptions);
}