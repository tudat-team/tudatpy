/*! \file physicalConstants.h
 *    This file contains a class with selected constants commonly used
 *    in astrodynamics.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 1
 *    Check status      : Unchecked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 3 September, 2010
 *    Last modified     : 24 January, 2011
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
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      YYMMDD    author        comment
 *      100906    J. Melman     First creation of code.
 *      110124    K. Kumar      Split into .h and .cpp files ( C++ standard ).
 */

#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

//! Physical constants structure.
/*!
 * Physical constants structure.
 */
struct PhysicalConstants
{
public:

    //! Julian day.
    /*!
     *  Julian day in seconds.
     */
    const static double JULIAN_DAY;

    //! Julian year in days.
    /*!
     *  Julian year in Julian days.
     */
    const static double JULIAN_YEAR_IN_DAYS;

    //! Julian year.
    /*!
     *  Julian year in seconds.
     *  Result of JULIAN_YEAR_IN_DAYS * JULIAN_DAY.
     */
    const static double JULIAN_YEAR;

    //! Sidereal day.
    /*!
     *  Sidereal day in seconds.
     *  Reference: http://ssd.jpl.nasa.gov/?constants#ref (temporary)
     */
    const static double SIDEREAL_DAY;

    //! Sidereal year in days.
    /*!
     *  Sidereal year in Julian days in quasar reference frame.
     *  Reference: http://ssd.jpl.nasa.gov/?constants#ref (temporary)
     */
    const static double SIDEREAL_YEAR_IN_DAYS;

    //! Sidereal year.
    /*!
     *  Sidereal year in seconds in quasar reference frame.
     *  Result of SIDEREAL_YEAR_IN_DAYS * JULIAN_DAY
     */
    const static double SIDEREAL_YEAR;

    //! Speed of light.
    /*!
     *  Speed of light in meters per second.
     *  Reference: Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards",
     *             in Highlights of Astronomy (I. Appenzeller, ed.), Table 1,
     *             Kluwer Academic Publishers, Dordrecht.
     *             http://iau-comm4.jpl.nasa.gov/iausgnsrpt.pdf
     */
    const static double SPEED_OF_LIGHT;

    //! Gravitational constant.
    /*!
     *  Gravitational constant in meter^3 per kilogram per second^2.
     *  Reference: Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards",
     *             in Highlights of Astronomy (I. Appenzeller, ed.), Table 1,
     *             Kluwer Academic Publishers, Dordrecht.
     *             http://iau-comm4.jpl.nasa.gov/iausgnsrpt.pdf
     */
    const static double GRAVITATIONAL_CONSTANT;

    //! Uncertainty of the gravitational constant.
    /*!
     *  Uncertainty of the gravitational constant in meter^3 per kilogram per second^2.
     *  Reference: Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards",
     *             in Highlights of Astronomy (I. Appenzeller, ed.), Table 1,
     *             Kluwer Academic Publishers, Dordrecht.
     *             http://iau-comm4.jpl.nasa.gov/iausgnsrpt.pdf
     */
    const static double UNCERTAINTY_GRAVITATIONAL_CONSTANT;

    //! Obliquity of the ecliptic in arcseconds.
    /*!
     *  Obliquity of the ecliptic in arcseconds at epoch J2000.
     *  Reference: Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards",
     *             in Highlights of Astronomy (I. Appenzeller, ed.), Table 1,
     *             Kluwer Academic Publishers, Dordrecht.
     *             http://iau-comm4.jpl.nasa.gov/iausgnsrpt.pdf
     */
    const static double OBLIQUITY_ECLIPTIC_IN_ARCSECONDS;

    //! Uncertainty of the obliquity of the ecliptic in arcseconds.
    /*!
     *  Uncertainty of the obliquity of the ecliptic in arcseconds at epoch J2000.
     *  Reference: Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards",
     *             in Highlights of Astronomy (I. Appenzeller, ed.), Table 1,
     *             Kluwer Academic Publishers, Dordrecht.
     *             http://iau-comm4.jpl.nasa.gov/iausgnsrpt.pdf
     */
    const static double UNCERTAINTY_OBLIQUITY_ECLIPTIC_IN_ARCSECONDS;

    //! Obliquity of the ecliptic in degrees.
    /*!
     *  Obliquity of the ecliptic in degrees at epoch J2000.
     *  Result of OBLIQUITY_ECLIPTIC_IN_ARCSECONDS * 60 * 60
     */
    const static double OBLIQUITY_ECLIPTIC_IN_DEGREES;

    //! Obliquity of the ecliptic in radians.
    /*!
     *  Obliquity of the ecliptic in radians at epoch J2000.
     *  Result of OBLIQUITY_ECLIPTIC_IN_ARCSECONDS * 60 * 60 * pi / 180.
     */
    const static double OBLIQUITY_ECLIPTIC;

    //! Astronomical Unit.
    /*!
     *  Astronomical Unit in meters.
     *  Reference: Standish, E.M. (1998) "JPL Planetary and Lunar Ephemerides, DE405/LE405",
     *             JPL IOM 312.F-98-048.
     *             http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
     */
    const static double ASTRONOMICAL_UNIT;

    //! Uncertainty of the Astronomical Unit.
    /*!
     *  Uncertainty of the Astronomical Unit in meters.
     *  Reference: Standish, E.M. (1998) "JPL Planetary and Lunar Ephemerides, DE405/LE405",
     *             JPL IOM 312.F-98-048.
     *             http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
     */
    const static double UNCERTAINTY_ASTRONOMICAL_UNIT;
};

#endif // PHYSICAL_CONSTANTS_H

// End of file.
