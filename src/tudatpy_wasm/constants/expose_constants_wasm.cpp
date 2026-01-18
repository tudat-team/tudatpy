/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Constants module bindings for WASM.
 *    Mirrors: src/tudatpy/constants/expose_constants.cpp
 *
 *    This module exposes physical and mathematical constants as
 *    module-level attributes accessible via tudat.constants.CONSTANT_NAME
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "expose_constants_wasm.h"
#include "../wasm_module.h"

// Tudat headers
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/constants.h"

namespace tbc = tudat::celestial_body_constants;
namespace tpc = tudat::physical_constants;
namespace tmc = tudat::mathematical_constants;
namespace tba = tudat::basic_astrodynamics;

// Register this module path
WASM_MODULE_PATH("constants")

EMSCRIPTEN_BINDINGS(tudatpy_constants) {
    using namespace emscripten;

    // ========================================================================
    // Physical Constants (from physicalConstants.h)
    // ========================================================================

    // Gravitational and geodetic constants
    constant("constants_SEA_LEVEL_GRAVITATIONAL_ACCELERATION",
             tpc::SEA_LEVEL_GRAVITATIONAL_ACCELERATION);
    constant("constants_GRAVITATIONAL_CONSTANT",
             tpc::GRAVITATIONAL_CONSTANT);

    // Time constants
    constant("constants_JULIAN_DAY", tpc::JULIAN_DAY);
    constant("constants_JULIAN_DAY_LONG", tpc::JULIAN_DAY_LONG);
    constant("constants_JULIAN_YEAR_IN_DAYS", tpc::JULIAN_YEAR_IN_DAYS);
    constant("constants_JULIAN_YEAR_IN_DAYS_LONG", tpc::JULIAN_YEAR_IN_DAYS_LONG);
    constant("constants_JULIAN_YEAR", tpc::JULIAN_YEAR);
    constant("constants_SIDEREAL_DAY", tpc::SIDEREAL_DAY);
    constant("constants_SIDEREAL_YEAR_IN_DAYS", tpc::SIDEREAL_YEAR_IN_DAYS);
    constant("constants_SIDEREAL_YEAR", tpc::SIDEREAL_YEAR);

    // Speed of light and derived
    constant("constants_SPEED_OF_LIGHT", tpc::SPEED_OF_LIGHT);
    constant("constants_SPEED_OF_LIGHT_LONG", tpc::SPEED_OF_LIGHT_LONG);
    constant("constants_INVERSE_SQUARE_SPEED_OF_LIGHT", tpc::INVERSE_SQUARE_SPEED_OF_LIGHT);
    constant("constants_INVERSE_CUBIC_SPEED_OF_LIGHT", tpc::INVERSE_CUBIC_SPEED_OF_LIGHT);
    constant("constants_INVERSE_QUARTIC_SPEED_OF_LIGHT", tpc::INVERSE_QUARTIC_SPEED_OF_LIGHT);
    constant("constants_INVERSE_QUINTIC_SPEED_OF_LIGHT", tpc::INVERSE_QUINTIC_SPEED_OF_LIGHT);

    // Distance and length constants
    constant("constants_ASTRONOMICAL_UNIT", tpc::ASTRONOMICAL_UNIT);

    // Gas constants
    constant("constants_SPECIFIC_GAS_CONSTANT_AIR", tpc::SPECIFIC_GAS_CONSTANT_AIR);
    constant("constants_MOLAR_GAS_CONSTANT", tpc::MOLAR_GAS_CONSTANT);

    // Quantum and thermodynamic constants
    constant("constants_PLANCK_CONSTANT", tpc::PLANCK_CONSTANT);
    constant("constants_BOLTZMANN_CONSTANT", tpc::BOLTZMANN_CONSTANT);
    constant("constants_STEFAN_BOLTZMANN_CONSTANT", tpc::STEFAN_BOLTZMANN_CONSTANT);

    // Electromagnetic constants
    constant("constants_VACUUM_PERMEABILITY", tpc::VACUUM_PERMEABILITY);
    constant("constants_VACUUM_PERMITTIVITY", tpc::VACUUM_PERMITTIVITY);

    // Time rate terms
    constant("constants_LG_TIME_RATE_TERM", tpc::LG_TIME_RATE_TERM);
    constant("constants_LG_TIME_RATE_TERM_LONG", tpc::LG_TIME_RATE_TERM_LONG);

    // ========================================================================
    // Time Epoch Constants (from timeConversions.h)
    // ========================================================================

    constant("constants_JULIAN_DAY_ON_J2000", tba::JULIAN_DAY_ON_J2000);
    constant("constants_JULIAN_DAY_AT_0_MJD", tba::JULIAN_DAY_AT_0_MJD);

    // ========================================================================
    // Mathematical Constants (from mathematicalConstants.h)
    // ========================================================================

    constant("constants_E", tmc::E);
    constant("constants_GOLDEN_RATIO", tmc::GOLDEN_RATIO);
    // Note: COMPLEX_I is a complex number, needs special handling
    // constant("constants_COMPLEX_I", tmc::COMPLEX_I);
    constant("constants_PI", tmc::PI);
    constant("constants_TUDAT_NAN", TUDAT_NAN);
}

#endif // __EMSCRIPTEN__
