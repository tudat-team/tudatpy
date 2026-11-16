#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include <pybind11/pybind11.h>
#include <tudat/config.hpp>

namespace py = pybind11;

// Declarations from each binding .cpp file
void add_math_to_kernel( py::module_& m );
void add_astro_types_to_kernel( py::module_& m );
void add_astro_to_kernel( py::module_& m );
void add_trajectory_design_to_kernel( py::module_& m );
void add_constants_to_kernel( py::module_& m );
void add_interface_to_kernel( py::module_& m );
void add_data_to_kernel( py::module_& m );
void add_data_input_to_kernel( py::module_& m );
void add_util_to_kernel( py::module_& m );
void add_dynamics_types_to_kernel( py::module_& m );
void add_dynamics_to_kernel( py::module_& m );
void add_estimation_types_to_kernel( py::module_& m );
void add_estimation_to_kernel( py::module_& m );
void add_exceptions_to_kernel( py::module_& m );

PYBIND11_MODULE( kernel, m )
{
    auto math = m.def_submodule( "math" );
    auto astro = m.def_submodule( "astro" );
    auto trajectory_design = m.def_submodule( "trajectory_design" );
    auto constants = m.def_submodule( "constants" );
    auto interface = m.def_submodule( "interface" );
    auto data = m.def_submodule( "data" );
    auto data_input = m.def_submodule( "data_input" );
    auto util = m.def_submodule( "util" );
    auto dynamics = m.def_submodule( "dynamics" );
    auto estimation = m.def_submodule( "estimation" );
    auto exceptions = m.def_submodule( "exceptions" );

    // Register types used across top-level module boundaries before binding
    // functions whose signatures contain them.
    add_astro_types_to_kernel( astro );
    add_estimation_types_to_kernel( estimation );

    add_math_to_kernel( math );
    add_dynamics_types_to_kernel( dynamics );
    add_astro_to_kernel( astro );
    add_trajectory_design_to_kernel( trajectory_design );
    add_constants_to_kernel( constants );
    add_interface_to_kernel( interface );
    add_data_to_kernel( data );
    add_data_input_to_kernel( data_input );
    add_util_to_kernel( util );
    add_dynamics_to_kernel( dynamics );
    add_estimation_to_kernel( estimation );
    add_exceptions_to_kernel( exceptions );
}
