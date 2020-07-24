/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_constants.h"

#include "docstrings.h"
#include "tudat/constants.h"

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace tbc = tudat::celestial_body_constants;
namespace tpc = tudat::physical_constants;
namespace tmc = tudat::mathematical_constants;

namespace tudatpy {

void expose_constants(py::module &m) {

  m.attr("__doc__") = tudatpy::constants_docstring().c_str();

  // celestialBodyConstants.h
  m.attr("EARTH_EQUATORIAL_RADIUS") = tbc::EARTH_EQUATORIAL_RADIUS;
  m.attr("EARTH_FLATTENING_FACTOR") = tbc::EARTH_FLATTENING_FACTOR;
  m.attr("EARTH_GEODESY_NORMALIZED_J2") = tbc::EARTH_GEODESY_NORMALIZED_J2;
  m.attr("SUN_GRAVITATIONAL_PARAMETER") = tbc::SUN_GRAVITATIONAL_PARAMETER;
  m.attr("MERCURY_GRAVITATIONAL_PARAMETER") = tbc::MERCURY_GRAVITATIONAL_PARAMETER;
  m.attr("VENUS_GRAVITATIONAL_PARAMETER") = tbc::VENUS_GRAVITATIONAL_PARAMETER;
  m.attr("EARTH_GRAVITATIONAL_PARAMETER") = tbc::EARTH_GRAVITATIONAL_PARAMETER;
  m.attr("MOON_GRAVITATIONAL_PARAMETER") = tbc::MOON_GRAVITATIONAL_PARAMETER;
  m.attr("MARS_GRAVITATIONAL_PARAMETER") = tbc::MARS_GRAVITATIONAL_PARAMETER;
  m.attr("JUPITER_GRAVITATIONAL_PARAMETER") = tbc::JUPITER_GRAVITATIONAL_PARAMETER;
  m.attr("SATURN_GRAVITATIONAL_PARAMETER") = tbc::SATURN_GRAVITATIONAL_PARAMETER;
  m.attr("URANUS_GRAVITATIONAL_PARAMETER") = tbc::URANUS_GRAVITATIONAL_PARAMETER;
  m.attr("NEPTUNE_GRAVITATIONAL_PARAMETER") = tbc::NEPTUNE_GRAVITATIONAL_PARAMETER;
  m.attr("PLUTO_GRAVITATIONAL_PARAMETER") = tbc::PLUTO_GRAVITATIONAL_PARAMETER;

  // physicalConstants.h
  m.attr("SEA_LEVEL_GRAVITATIONAL_ACCELERATION") = tpc::SEA_LEVEL_GRAVITATIONAL_ACCELERATION;
  m.attr("JULIAN_DAY") = tpc::JULIAN_DAY;
  m.attr("JULIAN_DAY_LONG") = tpc::JULIAN_DAY_LONG;
  m.attr("JULIAN_YEAR_IN_DAYS") = tpc::JULIAN_YEAR_IN_DAYS;
  m.attr("JULIAN_YEAR_IN_DAYS_LONG") = tpc::JULIAN_YEAR_IN_DAYS_LONG;
  m.attr("JULIAN_YEAR") = tpc::JULIAN_YEAR;
  m.attr("SIDEREAL_DAY") = tpc::SIDEREAL_DAY;
  m.attr("SIDEREAL_YEAR_IN_DAYS") = tpc::SIDEREAL_YEAR_IN_DAYS;
  m.attr("SIDEREAL_YEAR") = tpc::SIDEREAL_YEAR;
  m.attr("SPEED_OF_LIGHT") = tpc::SPEED_OF_LIGHT;
  m.attr("SPEED_OF_LIGHT_LONG") = tpc::SPEED_OF_LIGHT_LONG;
  m.attr("GRAVITATIONAL_CONSTANT") = tpc::GRAVITATIONAL_CONSTANT;
  m.attr("ASTRONOMICAL_UNIT") = tpc::ASTRONOMICAL_UNIT;
  m.attr("SPECIFIC_GAS_CONSTANT_AIR") = tpc::SPECIFIC_GAS_CONSTANT_AIR;
  m.attr("MOLAR_GAS_CONSTANT") = tpc::MOLAR_GAS_CONSTANT;
  m.attr("PLANCK_CONSTANT") = tpc::PLANCK_CONSTANT;
  m.attr("BOLTZMANN_CONSTANT") = tpc::BOLTZMANN_CONSTANT;
  m.attr("STEFAN_BOLTZMANN_CONSTANT") = tpc::STEFAN_BOLTZMANN_CONSTANT;
  m.attr("INVERSE_SQUARE_SPEED_OF_LIGHT") = tpc::INVERSE_SQUARE_SPEED_OF_LIGHT;
  m.attr("INVERSE_CUBIC_SPEED_OF_LIGHT") = tpc::INVERSE_CUBIC_SPEED_OF_LIGHT;
  m.attr("INVERSE_QUARTIC_SPEED_OF_LIGHT") = tpc::INVERSE_QUARTIC_SPEED_OF_LIGHT;
  m.attr("INVERSE_QUINTIC_SPEED_OF_LIGHT") = tpc::INVERSE_QUINTIC_SPEED_OF_LIGHT;
  m.attr("VACUUM_PERMEABILITY") = tpc::VACUUM_PERMEABILITY;
  m.attr("VACUUM_PERMITTIVITY") = tpc::VACUUM_PERMITTIVITY;
  m.attr("LG_TIME_RATE_TERM") = tpc::LG_TIME_RATE_TERM;
  m.attr("LG_TIME_RATE_TERM_LONG") = tpc::LG_TIME_RATE_TERM_LONG;

  // mathematicalConstants.h
  m.attr("E") = tmc::E;
  m.attr("GOLDEN_RATIO") = tmc::GOLDEN_RATIO;
  m.attr("COMPLEX_I") = tmc::COMPLEX_I;
  m.attr("PI") = tmc::PI;
  m.attr("LONG_PI") = tmc::LONG_PI;
  m.attr("TUDAT_NAN") = TUDAT_NAN;
};

}// namespace tudatpy
