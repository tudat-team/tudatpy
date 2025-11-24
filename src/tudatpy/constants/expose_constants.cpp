/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_constants.h"

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/constants.h"

namespace py = pybind11;
namespace tbc = tudat::celestial_body_constants;
namespace tpc = tudat::physical_constants;
namespace tmc = tudat::mathematical_constants;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{
namespace constants
{
void expose_constants( py::module& m )
{
    m.attr( "__doc__" ) = R"doc(
Physical, mathematical, and astronomical constants.

This module provides fundamental constants used throughout tudatpy for astrodynamics 
calculations. Constants are organized into several categories:

- **Physical constants**: Fundamental physical quantities (speed of light, gravitational constant, etc.)
- **Time constants**: Time conversion factors (Julian day, sidereal time, etc.)
- **Mathematical constants**: Mathematical values (π, e, golden ratio, etc.)
- **Astronomical constants**: Astronomical unit, Julian day epochs

High-precision variants
-----------------------
Several constants have ``_LONG`` variants that provide extended precision using 
``long double`` type instead of ``double``. These are useful for high-accuracy 
calculations over long time spans:

- ``JULIAN_DAY_LONG``: High-precision seconds per Julian day
- ``JULIAN_YEAR_IN_DAYS_LONG``: High-precision days per Julian year  
- ``SPEED_OF_LIGHT_LONG``: High-precision speed of light [m/s]
- ``LG_TIME_RATE_TERM_LONG``: High-precision TCG-TT time rate term

Time Scales
-----------
The module includes constants for relativistic time scale transformations:

- ``LG_TIME_RATE_TERM``: Relative rate difference between Geocentric Coordinate Time (TCG) 
  and Terrestrial Time (TT). TCG differs from TT by approximately 0.7 parts per billion.
- ``LB_TIME_RATE_TERM``: Relative rate difference between Barycentric Coordinate Time (TCB) 
  and Barycentric Dynamical Time (TDB).

References
----------
Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards", 
in Highlights of Astronomy (I. Appenzeller, ed.), Table 1, Kluwer Academic Publishers, 
Dordrecht.

NASA SSD. Astrodynamic Parameters, http://ssd.jpl.nasa.gov/?constants#ref, 
6th September, 2011, last accessed: 21st February, 2012.

IAU 2012, Resolution B2. https://www.iau.org/static/resolutions/IAU2012_English.pdf

Anderson, J.D. Jr. Hypersonic and High-Temperature Gas Dynamics, Second Edition, 2006.

NIST. NIST reference on constants, units and uncertainty.
http://physics.nist.gov/cuu/Constants/index.html

)doc";

    // physicalConstants.h
    m.attr( "SEA_LEVEL_GRAVITATIONAL_ACCELERATION" ) = tpc::SEA_LEVEL_GRAVITATIONAL_ACCELERATION;
    m.attr( "JULIAN_DAY" ) = tpc::JULIAN_DAY;
    m.attr( "JULIAN_DAY_LONG" ) = tpc::JULIAN_DAY_LONG;
    m.attr( "JULIAN_YEAR_IN_DAYS" ) = tpc::JULIAN_YEAR_IN_DAYS;
    m.attr( "JULIAN_YEAR_IN_DAYS_LONG" ) = tpc::JULIAN_YEAR_IN_DAYS_LONG;
    m.attr( "JULIAN_YEAR" ) = tpc::JULIAN_YEAR;
    m.attr( "SIDEREAL_DAY" ) = tpc::SIDEREAL_DAY;
    m.attr( "SIDEREAL_YEAR_IN_DAYS" ) = tpc::SIDEREAL_YEAR_IN_DAYS;
    m.attr( "SIDEREAL_YEAR" ) = tpc::SIDEREAL_YEAR;
    m.attr( "SPEED_OF_LIGHT" ) = tpc::SPEED_OF_LIGHT;
    m.attr( "SPEED_OF_LIGHT_LONG" ) = tpc::SPEED_OF_LIGHT_LONG;
    m.attr( "GRAVITATIONAL_CONSTANT" ) = tpc::GRAVITATIONAL_CONSTANT;
    m.attr( "ASTRONOMICAL_UNIT" ) = tpc::ASTRONOMICAL_UNIT;
    m.attr( "SPECIFIC_GAS_CONSTANT_AIR" ) = tpc::SPECIFIC_GAS_CONSTANT_AIR;
    m.attr( "MOLAR_GAS_CONSTANT" ) = tpc::MOLAR_GAS_CONSTANT;
    m.attr( "PLANCK_CONSTANT" ) = tpc::PLANCK_CONSTANT;
    m.attr( "BOLTZMANN_CONSTANT" ) = tpc::BOLTZMANN_CONSTANT;
    m.attr( "STEFAN_BOLTZMANN_CONSTANT" ) = tpc::STEFAN_BOLTZMANN_CONSTANT;
    m.attr( "INVERSE_SQUARE_SPEED_OF_LIGHT" ) = tpc::INVERSE_SQUARE_SPEED_OF_LIGHT;
    m.attr( "INVERSE_CUBIC_SPEED_OF_LIGHT" ) = tpc::INVERSE_CUBIC_SPEED_OF_LIGHT;
    m.attr( "INVERSE_QUARTIC_SPEED_OF_LIGHT" ) = tpc::INVERSE_QUARTIC_SPEED_OF_LIGHT;
    m.attr( "INVERSE_QUINTIC_SPEED_OF_LIGHT" ) = tpc::INVERSE_QUINTIC_SPEED_OF_LIGHT;
    m.attr( "VACUUM_PERMEABILITY" ) = tpc::VACUUM_PERMEABILITY;
    m.attr( "VACUUM_PERMITTIVITY" ) = tpc::VACUUM_PERMITTIVITY;
    m.attr( "LG_TIME_RATE_TERM" ) = tpc::LG_TIME_RATE_TERM;
    m.attr( "LG_TIME_RATE_TERM_LONG" ) = tpc::LG_TIME_RATE_TERM_LONG;

    // time constants
    m.attr( "JULIAN_DAY_ON_J2000" ) = tba::JULIAN_DAY_ON_J2000;
    m.attr( "JULIAN_DAY_AT_0_MJD" ) = tba::JULIAN_DAY_AT_0_MJD;

    // mathematicalConstants.h
    m.attr( "E" ) = tmc::E;
    m.attr( "GOLDEN_RATIO" ) = tmc::GOLDEN_RATIO;
    m.attr( "COMPLEX_I" ) = tmc::COMPLEX_I;
    m.attr( "PI" ) = tmc::PI;
    m.attr( "TUDAT_NAN" ) = TUDAT_NAN;
};
}  // namespace constants
}  // namespace tudatpy
