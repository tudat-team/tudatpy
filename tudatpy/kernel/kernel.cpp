#include <pybind11/pybind11.h>

#include <tudat/config.hpp>

#include "expose_astro.h"
#include "expose_constants.h"
#include "expose_example.h"
#include "expose_interface.h"
#include "expose_io.h"
#include "expose_math.h"
#include "expose_numerical_simulation.h"
#include "expose_trajectory_design.h"
#include "expose_utils.h"

namespace py = pybind11;

PYBIND11_MODULE(kernel, m) {

    // Disable automatic function signatures in the docs.
    // NOTE: the 'options' object needs to stay alive
    // throughout the whole definition of the module.
    py::options options;
    options.enable_function_signatures();
    options.enable_user_defined_docstrings();

    // Export the tudat version.
    m.attr("_tudat_version") = TUDAT_VERSION;
    m.attr("_tudat_version_major") = TUDAT_VERSION_MAJOR;
    m.attr("_tudat_version_minor") = TUDAT_VERSION_MINOR;
    m.attr("_tudat_version_patch") = TUDAT_VERSION_PATCH;

    // math module
    auto utils = m.def_submodule("utils");
    tudatpy::utils::expose_utils(utils);

    // math module
    auto math = m.def_submodule("math");
    tudatpy::math::expose_math(math);

    // astro module
    auto astro = m.def_submodule("astro");
    tudatpy::astro::expose_astro(astro);

    // interface module
    auto interface = m.def_submodule("interface");
    tudatpy::interface::expose_interface(interface);

    // constants module
    auto constants = m.def_submodule("constants");
    tudatpy::constants::expose_constants(constants);

    // io module
    auto io = m.def_submodule("io");
    tudatpy::io::expose_io(io);

    // simulation module
    auto trajectory_design = m.def_submodule("trajectory_design");
    tudatpy::trajectory_design::expose_trajectory_design(trajectory_design);

    // simulation module
    auto numerical_simulation = m.def_submodule("numerical_simulation");
    tudatpy::numerical_simulation::expose_numerical_simulation(numerical_simulation);

//    // example module
//    auto example = m.def_submodule("example");
//    tudatpy::expose_example(example);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
