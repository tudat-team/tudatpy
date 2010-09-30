/*! \file predefinedPlanets.h
 *    This file contains instantiations of the CelestialBody class that
 *    assign default values to the characteristics of certain planets.
 *
 *    Path              : /Astrodynamics/Bodies/CelestialBodies/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 20 September, 2010
 *    Last modified     : 20 September, 2010
 *
 *    References
 *
 *    Notes
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
 *      100920    J. Melman     First creation of code.
 *      100929    K. Kumar      Minor comment changes.
 */

#ifndef PREDEFINED_PLANETS_H
#define PREDEFINED_PLANETS_H

// Include statements.
#include "celestialBody.h"
#include "predefinedGravityFieldParameters.h"

//! Predefined planets namespace.
/*!
 * Predefined planets namespace.
 */
namespace predefined_planets
{

// Definition of Earth using object of CelestialBody class.
const CelestialBody earth(
    predefined_gravity_field_parameters::defaultEarthGravityFieldParameters );

}

#endif // PREDEFINED_PLANETS_H

// End of file.
