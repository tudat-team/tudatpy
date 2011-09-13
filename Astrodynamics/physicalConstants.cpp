/*! \file physicalConstants.cpp
 *    This file contains the implementation of a data structure with selected
 *    constants commonly used in astrodynamics.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@student.tudelft.nl
 *
 *    Date created      : 24 January, 2011
 *    Last modified     : 29 June, 2011
 *
 *    References
 *      Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical
 *          Standards", in Highlights of Astronomy (I. Appenzeller, ed.),
 *          Table 1, Kluwer Academic Publishers, Dordrecht.
 *      Standish, E.M. (1998) "JPL Planetary and Lunar Ephemerides,
 *          DE405/LE405", JPL IOM 312.F-98-048.
 *
 *    Notes
 *      The reference for the sidereal day and year should be updated.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110124    K. Kumar          First creation of code.
 *      110629    F.M. Engelen      Added specific gas constants.
 */

// Include directives.
#include "physicalConstants.h"

//! Julian day.
const double PhysicalConstants::JULIAN_DAY = 86400.0;

//! Julian year in days.
const double PhysicalConstants::JULIAN_YEAR_IN_DAYS = 365.25;

//! Julian year.
const double PhysicalConstants::JULIAN_YEAR = 3.15576e7;

//! Sidereal day.
const double PhysicalConstants::SIDEREAL_DAY = 86164.09054;

//! Sidereal year in days.
const double PhysicalConstants::SIDEREAL_YEAR_IN_DAYS = 365.25636;

//! Sidereal year.
const double PhysicalConstants::SIDEREAL_YEAR = 3.1558149504e7;

//! Speed of light.
const double PhysicalConstants::SPEED_OF_LIGHT = 299792458;

//! Gravitational constant.
const double PhysicalConstants::GRAVITATIONAL_CONSTANT = 6.67259e-11;

//! Uncertainty of the gravitational constant.
const double PhysicalConstants
    ::UNCERTAINTY_GRAVITATIONAL_CONSTANT = 0.00030e-11;

//! Obliquity of the ecliptic in arcseconds.
const double PhysicalConstants::OBLIQUITY_ECLIPTIC_IN_ARCSECONDS = 84381.412;

//! Uncertainty of the obliquity of the ecliptic in arcseconds.
const double PhysicalConstants
    ::UNCERTAINTY_OBLIQUITY_ECLIPTIC_IN_ARCSECONDS = 0.005;

//! Obliquity of the ecliptic in degrees.
const double PhysicalConstants::OBLIQUITY_ECLIPTIC_IN_DEGREES = 23.439281;

//! Obliquity of the ecliptic in radians.
const double PhysicalConstants::OBLIQUITY_ECLIPTIC = 0.40909263;

//! Astronomical Unit.
const double PhysicalConstants::ASTRONOMICAL_UNIT = 1.49597870691e11;

//! Uncertainty of the Astronomical Unit.
const double PhysicalConstants::UNCERTAINTY_ASTRONOMICAL_UNIT = 3.0;

//! Specific gas constant of air.
const double PhysicalConstants::SPECIFIC_GAS_CONSTANT_AIR = 2.87e2;

//! Ratio of specific heats of air/diatomic gases.
const double PhysicalConstants::RATIO_OF_SPECIFIC_HEATS_AIR = 1.4;

// End of file.
