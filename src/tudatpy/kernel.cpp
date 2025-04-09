#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <pybind11/pybind11.h>

#include <tudat/config.hpp>

#include "astro/expose_astro.h"
#include "constants/expose_constants.h"
#include "data/expose_data.h"
#include "interface/expose_interface.h"
#include "math/expose_math.h"
#include "numerical_simulation/expose_numerical_simulation.h"
#include "trajectory_design/expose_trajectory_design.h"

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

    // simulation module
    auto numerical_simulation = m.def_submodule( "numerical_simulation" );
    tudatpy::numerical_simulation::expose_numerical_simulation( numerical_simulation );
    tudatpy::numerical_simulation::expose_numerical_simulation_simulator( numerical_simulation );
    tudatpy::numerical_simulation::expose_numerical_simulation_variational( numerical_simulation );
    tudatpy::numerical_simulation::expose_numerical_simulation_estimator( numerical_simulation );

    //    // example module
    //    auto example = m.def_submodule("example");
    //    tudatpy::expose_example(example);
}
