//
// Created by ggarrett on 24-04-20.
//
#include <string>

namespace tudatpy {

// astro
std::string astro_docstring() {

  return 0;
}

// astro : two_body
std::string two_body_docstring() {

  return 0;
}



// constants
std::string constants_docstring() {
  return R"mydelimiter(
Summary
-------
Module containing all constants listed in the Tudat package.

Routine Listings
----------------
 - EARTH_EQUATORIAL_RADIUS
 - EARTH_FLATTENING_FACTOR
 - EARTH_GEODESY_NORMALIZED_J2
 - SUN_GRAVITATIONAL_PARAMETER
 - MERCURY_GRAVITATIONAL_PARAMETER
 - VENUS_GRAVITATIONAL_PARAMETER
 - EARTH_GRAVITATIONAL_PARAMETER
 - MOON_GRAVITATIONAL_PARAMETER
 - MARS_GRAVITATIONAL_PARAMETER
 - JUPITER_GRAVITATIONAL_PARAMETER
 - SATURN_GRAVITATIONAL_PARAMETER
 - URANUS_GRAVITATIONAL_PARAMETER
 - NEPTUNE_GRAVITATIONAL_PARAMETER
 - PLUTO_GRAVITATIONAL_PARAMETER
 - SEA_LEVEL_GRAVITATIONAL_ACCELERATION
 - JULIAN_DAY
 - JULIAN_DAY_LONG
 - JULIAN_YEAR_IN_DAYS
 - JULIAN_YEAR_IN_DAYS_LONG
 - JULIAN_YEAR
 - SIDEREAL_DAY
 - SIDEREAL_YEAR_IN_DAYS
 - SIDEREAL_YEAR
 - SPEED_OF_LIGHT
 - SPEED_OF_LIGHT_LONG
 - GRAVITATIONAL_CONSTANT
 - ASTRONOMICAL_UNIT
 - SPECIFIC_GAS_CONSTANT_AIR
 - MOLAR_GAS_CONSTANT
 - PLANCK_CONSTANT
 - BOLTZMANN_CONSTANT
 - STEFAN_BOLTZMANN_CONSTANT
 - INVERSE_SQUARE_SPEED_OF_LIGHT
 - INVERSE_CUBIC_SPEED_OF_LIGHT
 - INVERSE_QUARTIC_SPEED_OF_LIGHT
 - INVERSE_QUINTIC_SPEED_OF_LIGHT
 - VACUUM_PERMEABILITY
 - VACUUM_PERMITTIVITY
 - LG_TIME_RATE_TERM
 - LG_TIME_RATE_TERM_LONG
 - E
 - GOLDEN_RATIO
 - COMPLEX_I
 - PI
 - LONG_PI
 - TUDAT_NAN
        )mydelimiter";
};

// spice_interface
std::string load_standard_kernels_docstring() {
  return R"mydelimiter(
Definition at line 283 of file ``spiceInterface.cpp``.

References ``tudat::input_output::getSpiceKernelPath()``, and ``loadSpiceKernelInTudat()``.

Referenced by ``tudat::unit_tests::executeEarthOrbiterBiasEstimation()``, ``tudat::unit_tests::executeEarthOrbiterParameterEstimation()``, ``tudat::unit_tests::executePlanetaryParameterEstimation()``, and ``tudat::json_interface::loadSpiceKernels()``.
        )mydelimiter";
};

std::string clear_spice_kernels_docstring() {
  return R"mydelimiter(
This function removes all Spice kernels from the kernel pool. Wrapper for the kclear_c function.

Definition at line 281 of file ``spiceInterface.cpp``.

Referenced by ``tudat::json_interface::loadSpiceKernels()``.
        )mydelimiter";
};

// simulation_setup
std::string body_settings_docstring() {
  return R"mydelimiter(
Class in which the general properties of each environment model can be set (see above for the list of the available types of environment models). We note that for :class:`Body` objects that represent vehicles, the manual creation is typically used, as the vehicle conditions may depend on the celestial bodies, but not vice versa.

In many cases, default properties of (celestial) bodies may be used by calling the :literal:`getDefaultBodySettings` function, so that the user does not need to define all required properties line-by-line. At present, the following default settings are used (none if not in list):

Parameters
----------
- Ephemeris: Tabulated ephemeris created from Spice (valid in the interval that is specified by the input time-arguments to getDefaultBodySettings).
- Gravity field models: Point mass gravity field models, with gravitational parameter from Spice (if available). Exceptions are the Earth and Moon, for which the EGM96 and GGLP spherical harmonic gravity fields are loaded, respectively.
- Rotation model: For a given body (if available) the Spice rotation model, with ECLIPJ2000 as base frame, and for a body AAA frame IAU_AAA as target frame (the standard body-fixed frame for each body in Spice).
- Atmosphere model: 1976 US Standard Atmosphere for Earth (using pregenerated tables). For other bodies, no default shape model is given.
- Shape model: Spherical model with mean radius obtained from Spice (if avaiable).
        )mydelimiter";
};
}// namespace tudatpy