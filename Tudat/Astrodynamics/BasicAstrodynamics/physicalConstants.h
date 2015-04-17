/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100906    J. Melman         First creation of code.
 *      110124    K. Kumar          Split into .h and .cpp files.
 *      110629    F.M. Engelen      Added specific gas constants.
 *      120127    D. Dirkx          Moved to Tudat core, removed variables related to
 *                                  obliquity of ecliptic.
 *      120203    K. Kumar          Added missing specific gas constant value; need unit test.
 *      121212    K. Kumar          Migrated namespace to directory-based protocol and added
 *                                  backwards compatibility.
 *      130111    D. Dirkx          Added Planck, Boltzmann, and Stefan-Boltzmann constants.
 *
 *    References
 *      Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards",
 *          in Highlights of Astronomy (I. Appenzeller, ed.), Table 1, Kluwer Academic Publishers,
 *          Dordrecht.
 *      Standish, E.M. (1998) "JPL Planetary and Lunar Ephemerides, DE405/LE405",
 *          JPL IOM 312. F-98-048.
 *      Anderson, J.D. Jr. Hypersonic and High-Temperature Gas Dynamics, Second Edition, p469,
 *          2006.
 *      NASA. Astrodynamics Constants, http://ssd.jpl.nasa.gov/?constants#ref, 6th September, 2011,
 *          last accessed: 21st February, 2012.
 *      NIST. NIST reference on constants, units and uncertainty.
 *          http://physics.nist.gov/cuu/Constants/index.html, last accessed: 11th January, 2013.
 *      Wolfram Research, http://scienceworld.wolfram.com/physics/Stefan-BoltzmannLaw.html,
 *          last accessed: 11th January 2013.
 *
 *    Notes
 *      Backwards compatibility of namespaces is implemented for Tudat Core 2 in this file. The
 *      code block marked "DEPRECATED!" at the end of the file should be removed in Tudat Core 3.
 *
 */

#ifndef TUDAT_PHYSICAL_CONSTANTS_H
#define TUDAT_PHYSICAL_CONSTANTS_H

#include <cmath>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace physical_constants
{

//! Julian day.
/*!
 * Julian day in seconds (NASA, 2012).
 */
const static double JULIAN_DAY = 86400.0;

//! Julian year in days.
/*!
 * Julian year in Julian days (NASA, 2012).
 */
const static double JULIAN_YEAR_IN_DAYS = 365.25;

//! Julian year.
/*!
 * Julian year in seconds. Result of JULIAN_YEAR_IN_DAYS * JULIAN_DAY.
 */
const static double JULIAN_YEAR = 3.15576e7;

//! Sidereal day.
/*!
 * Sidereal day in seconds (NASA, 2012).
 */
const static double SIDEREAL_DAY = 86164.09054;

//! Sidereal year in days.
/*!
 * Sidereal year in Julian days in quasar reference frame (NASA, 2012).
 */
const static double SIDEREAL_YEAR_IN_DAYS = 365.25636;

//! Sidereal year.
/*!
 * Sidereal year in seconds in quasar reference frame. Result of
 * SIDEREAL_YEAR_IN_DAYS * JULIAN_DAY.
 */
const static double SIDEREAL_YEAR = 3.1558149504e7;

//! Speed of light.
/*!
 * Speed of light in meters per second (Standish, 1995).
 */
const static double SPEED_OF_LIGHT = 299792458.0;

//! Gravitational constant.
/*!
 * Gravitational constant in meter^3 per kilogram per second^2, (Standish, 1995).
 */
const static double GRAVITATIONAL_CONSTANT = 6.67259e-11;

//! Astronomical Unit.
/*!
 * Astronomical Unit in meters (Standish, 1998).
 */
const static double ASTRONOMICAL_UNIT = 1.49597870691e11;

//! Specific gas constant of air.
/*!
 * The specific gas constant of air in J per kg Kelvin (J/(kg K)) (Anderson, 2006).
 */
const static double SPECIFIC_GAS_CONSTANT_AIR = 2.87e2;

//! Planck constant.
/*!
 * Planck's constant in m^{2} kg/s, (NIST, 2013).
 */
const static double PLANCK_CONSTANT = 6.62606957E-34;

//! Boltzmann constant.
/*!
 * The Boltzmann constant (gas constant per particle) in  m^{2} kg / ( s^{2} K ), (NIST, 2013).
 */
const static double BOLTZMANN_CONSTANT = 1.3806488E-23;

//! Stefan-Boltzmann constant.
/*!
 * Stefan-Boltzmann constant, used for calculating black body radiation intensity, (Wolfram, 2013)
 * in J / (s m^{2} K{4} )
 */
const static double STEFAN_BOLTZMANN_CONSTANT = 2.0 *
        std::pow( mathematical_constants::PI, 5.0 ) *
        std::pow( BOLTZMANN_CONSTANT, 4.0 ) /
        ( 15.0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT *
          PLANCK_CONSTANT * PLANCK_CONSTANT * PLANCK_CONSTANT );

} // namespace physical_constants

} // namespace tudat

#endif // TUDAT_PHYSICAL_CONSTANTS_H
