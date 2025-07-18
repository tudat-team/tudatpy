#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <pybind11/pybind11.h>

#include <tudat/config.hpp>

#include "astro/expose_astro.h"
#include "constants/expose_constants.h"
#include "data/expose_data.h"
#include "exceptions/expose_exceptions.h"
#include "interface/expose_interface.h"
#include "math/expose_math.h"
#include "trajectory_design/expose_trajectory_design.h"
#include "dynamics/expose_dynamics.h"
#include "estimation/expose_estimation.h"

namespace py = pybind11;

PYBIND11_MODULE( kernel, m )
{
    // Disable automatic function signatures in the docs.
    // NOTE: the 'options' object needs to stay alive
    // throughout the whole definition of the module.
    // py::options options;
    // options.enable_function_signatures( );
    // options.enable_user_defined_docstrings( );

    // // Export the tudat version.
    // m.attr( "_tudat_version" ) = TUDAT_VERSION;
    // m.attr( "_tudat_version_major" ) = TUDAT_VERSION_MAJOR;
    // m.attr( "_tudat_version_minor" ) = TUDAT_VERSION_MINOR;
    // m.attr( "_tudat_version_patch" ) = TUDAT_VERSION_PATCH;

    // math module
    auto math = m.def_submodule( "math" );
    tudatpy::math::expose_math( math );

    // astro module
    auto astro = m.def_submodule( "astro" );
    tudatpy::astro::expose_astro( astro );

    // simulation module
    auto trajectory_design = m.def_submodule( "trajectory_design" );
    tudatpy::trajectory_design::expose_trajectory_design( trajectory_design );

    // constants module
    auto constants = m.def_submodule( "constants" );
    tudatpy::constants::expose_constants( constants );

    // interface module
    auto interface = m.def_submodule( "interface" );
    tudatpy::interface::expose_interface( interface );

    // data module
    auto data = m.def_submodule( "data" );
    tudatpy::data::expose_data( data );

    // dynamics module
    auto dynamics = m.def_submodule( "dynamics" );
    tudatpy::dynamics::expose_dynamics( dynamics );

    // estimation module
    auto estimation = m.def_submodule( "estimation" );
    tudatpy::estimation::expose_estimation( estimation );

    // exceptions module
    auto exceptions = m.def_submodule( "exceptions" );
    tudatpy::exceptions::expose_exceptions( exceptions );

    //    // example module
    //    auto example = m.def_submodule("example");
    //    tudatpy::expose_example(example);
}
