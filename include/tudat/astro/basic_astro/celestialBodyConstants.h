/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      International Earth Rotation and Reference System Service, Conventions 2010 (IERS2010),
 *        http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
 *      JPL, NASA. Astrodynamic Constants, http://ssd.jpl.nasa.gov/?constants,
 *        last updated: 13 Dec, 2012, last accessed: 19th March, 2013.
 *      JPL, NASA. Introduction to SPICE, 
 *        http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/info/intrdctn.html, last updated:
 *        29th Feb, 2008, last accessed: 22nd March, 2012. 
 *      EGM 2008 gravity field model, http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/
 *        last accessed: 21 March, 2013.
 *      IAU, NOMINAL VALUES FOR SELECTED SOLAR AND PLANETARY QUANTITIES: IAU 2015 RESOLUTION B3, 2016.
 *        https://doi.org/10.3847/0004-6256/152/2/41
 *
 *    Notes
 *      All the gravitational constants given here are taken from NASA JPL's list of astrodynamic
 *      constants (JPL, 2012), except for the Earth (IERS, 2010).
 *
 *      The data in this file is not meant to be used for precise orbit applications, such as
 *      orbit determination or Earth-orbiting mission operations analysis. For such applications
 *      specific, up-to-date constants, should be retrieved from, for instance the IERS or IAU.
 *      The Spice libraries provide a range of additional constants too (JPL, 2008).
 *
 */

#ifndef TUDAT_CELESTIAL_BODY_CONSTANTS_H
#define TUDAT_CELESTIAL_BODY_CONSTANTS_H

#include <map>

namespace tudat
{
namespace celestial_body_constants
{

//! Earth equatorial radius.
/*!
 * Earth equatorial radius in meters [m] (IERS, 2010).
 */
const static double EARTH_EQUATORIAL_RADIUS = 6378136.6;

//! Earth flattening factor.
/*!
 * Earth flattening factor [-] (IERS, 2010).
 */
const static double EARTH_FLATTENING_FACTOR = 298.25642;

//! Earth spherical harmonic gravitational field J2
/*!
 * Earth spherical harmonic gravitational field J2 [-] (EGM2008, 2013)
 */
const static double EARTH_GEODESY_NORMALIZED_J2 = -0.484165143790815E-03;

//! Solar gravitational constant.
/*!
 * Solar gravitational constant in squared meters per second squared [m^3 s^-2] (JPL, 2012).
 */
const static double SUN_GRAVITATIONAL_PARAMETER = 1.32712440018e20;

//! Mercury gravitational constant.
/*!
 * Gravitational constant of Mercury, derived from Sun/Mercury mass ratio [m^3 s^-2] (JPL, 2012).
 */
const static double MERCURY_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / 6023600.0;

//! Venus gravitational constant.
/*!
 * Gravitational constant of Venus, derived from Sun/Venus mass ratio [m^3 s^-2] (JPL, 2012).
 */
const static double VENUS_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / 408523.71;

//! Earth gravitational constant.
/*!
 * Earth gravitational constant in squared meters per second squared [m^3 s^-2] (IERS, 2012).
 */
const static double EARTH_GRAVITATIONAL_PARAMETER = 3.986004418E14;

//! Moon gravitational constant.
/*!
 * Moon gravitational constant in squared meters per second squared [m^3 s^-2] (JPL, 2012).
 */
const static double MOON_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER /
        ( 328900.56 * ( 1.0 + 81.30059 ) );

//! Mars gravitational constant.
/*!
 * Gravitational constant of Mars, derived from Sun/Mars mass ratio [m^3 s^-2] (JPL, 2012).
 */
const static double MARS_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / 3098708.0;

//! Jupiter gravitational constant.
/*!
 * Gravitational constant of Jupiter, derived from Sun/Jupiter mass ratio [m^3 s^-2] (JPL, 2012).
 */
const static double JUPITER_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / 1047.3486;

//! Saturn gravitational constant.
/*!
 * Gravitational constant of Saturn, derived from Sun/Saturn mass ratio [m^3 s^-2] (JPL, 2012).
 */
const static double SATURN_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / 3497.898 ;

//! Uranus gravitational constant.
/*!
 * Gravitational constant of Uranus, derived from Sun/Uranus mass ratio [m^3 s^-2] (JPL, 2012).
 */
const static double URANUS_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / 22902.98;

//! Neptune gravitational constant.
/*!
 * Gravitational constant of Neptune, derived from Sun/Neptune mass ratio [m^3 s^-2] (JPL, 2012).
 */
const static double NEPTUNE_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / 19412.24;

//! Pluto gravitational constant.
/*!
 * Gravitational constant of Pluto, derived from Sun/Pluto mass ratio [m^3 s^-2] (JPL, 2012).
 */
const static double PLUTO_GRAVITATIONAL_PARAMETER = SUN_GRAVITATIONAL_PARAMETER / 1.35e8;

//! Solar luminosity.
/*!
 * Solar luminosity (or radiant power) [W] (IAU, 2016).
 */
const static double SUN_LUMINOSITY = 3.828e26;

//static std::map< int, std::string > planetNames =
//{
//    { 0, "Sun" },
//    { 1, "Mercury" },
//    { 2, "Venus" },
//    { 3, "Earth" },
//    { 4, "Mars" },
//    { 5, "Jupiter" },
//    { 6, "Saturn" },
//    { 7, "Uranus" },
//    { 8, "Neptune" },
//    { 9, "Pluto" }
//};

//static std::map< std::string, int > planetIdNumbers =
//{
//    { "Sun", 0 },
//    { "Mercury", 1 },
//    { "Venus", 2 },
//    { "Earth", 3 },
//    { "Mars", 4 },
//    { "Jupiter", 5 },
//    { "Saturn", 6 },
//    { "Uranus", 7 },
//    { "Neptune", 8 },
//    { "Pluto", 9 }
//};

} // namespace celestial_body_constants

} // namespace tudat

#endif // TUDAT_CELESTIAL_BODY_CONSTANTS_H
