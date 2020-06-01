#include <pybind11/pybind11.h>

#include "expose_constants.h"
#include "expose_root_finders.h"
#include "expose_interpolators.h"
#include "expose_spice_interface.h"
#include "expose_ephemerides.h"
#include "expose_reference_frames.h"
#include "expose_aerodynamics.h"
#include "expose_basic_astrodynamics.h"
#include "expose_gravitation.h"
#include "expose_numerical_integrators.h"
#include "expose_orbital_element_conversions.h"
#include "expose_propagators.h"
#include "expose_simulation_setup.h"
#include "expose_unit_tests.h"

#include "tudat/astro/aerodynamics/tests/testApolloCapsuleCoefficients.h"

namespace py = pybind11;

PYBIND11_MODULE (core, m) {

    // Disable automatic function signatures in the docs.
    // NOTE: the 'options' object needs to stay alive
    // throughout the whole definition of the module.
    py::options options;
    options.disable_function_signatures();
    options.enable_user_defined_docstrings();

    // Export the tudatBundle version.
//    m.attr("_tudatbundle_version_major") = TUDATBUNDLE_VERSION_MAJOR;
//    m.attr("_tudatbundle_version_minor") = TUDATBUNDLE_VERSION_MINOR;
//    m.attr("_tudatbundle_version_patch") = TUDATBUNDLE_VERSION_PATCH;

    // constants
    auto constants = m.def_submodule("_constants");
    tudatpy::expose_constants(constants);

    // root_finders module
    auto root_finders = m.def_submodule("_root_finders");
    tudatpy::expose_root_finders(root_finders);

    // interpolators
    auto interpolators = m.def_submodule("_interpolators");
    tudatpy::expose_interpolators(interpolators);

    // spice_interface module
    auto spice_interface = m.def_submodule("_spice_interface");
    tudatpy::expose_spice_interface(spice_interface);

    // ephemerides module
    auto ephemerides = m.def_submodule("_ephemerides");
    tudatpy::expose_ephemerides(ephemerides);

    // reference_frames module
    auto reference_frames = m.def_submodule("_reference_frames");
    tudatpy::expose_reference_frames(reference_frames);

    // aerodynamics module
    auto aerodynamics = m.def_submodule("_aerodynamics");
    tudatpy::expose_aerodynamics(aerodynamics);

    // basic_astrodynamics module
    auto basic_astrodynamics = m.def_submodule("_basic_astrodynamics");
    tudatpy::expose_basic_astrodynamics(basic_astrodynamics);

    // gravitation module
    auto gravitation = m.def_submodule("_gravitation");
    tudatpy::expose_gravitation(gravitation);

    // numerical_integrators module
    auto numerical_integrators = m.def_submodule("_numerical_integrators");
    tudatpy::expose_numerical_integrators(numerical_integrators);

    // orbital_element_conversions module
    auto orbital_element_conversions = m.def_submodule("_orbital_element_conversions");
    tudatpy::expose_orbital_element_conversions(orbital_element_conversions);

    // propagators module
    auto propagators = m.def_submodule("_propagators");
    tudatpy::expose_propagators(propagators);

    // simulation_setup module
    auto simulation_setup = m.def_submodule("_simulation_setup");
    tudatpy::expose_simulation_setup(simulation_setup);

    // unit_tests module
    auto unit_tests = m.def_submodule("_unit_tests");
    tudatpy::expose_unit_tests(unit_tests);

    unit_tests.def("get_apollo_coefficient_interface",
                   tudat::unit_tests::getApolloCoefficientInterface,
                   "<no doc>");


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif

}
