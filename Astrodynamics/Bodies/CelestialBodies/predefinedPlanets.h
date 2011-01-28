/*! \file predefinedPlanets.h
 *    Header file that defines a namespace that contains predefined planets in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Bodies/CelestialBodies/
 *    Version           : 5
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
 *    Last modified     : 28 January, 2011
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
 *      YYMMDD    Author        Comment
 *      100920    J. Melman     First creation of code.
 *      100929    K. Kumar      Minor comment changes.
 *      110112    K. Kumar      Updated code to use enum and
 *                              createPredefinedPlanet() function.
 *      110114    J. Melman     Corrected Doxygen comments.
 *      110128    K. Kumar      Added Mars.
 */

#ifndef PREDEFINED_PLANETS_H
#define PREDEFINED_PLANETS_H

// Include statements.
#include "celestialBody.h"
#include "sphericalHarmonicsGravityField.h"
#include "predefinedGravityFieldModels.h"

// Using declarations.
using std::cerr;
using std::endl;

//! Predefined planets namespace.
/*!
 * Predefined planets namespace.
 */
namespace predefined_planets
{

//! Predefined planets.
/*!
 * Predefined planets.
 */
enum PredefinedPlanets
{
    earth,
    mars
};

//! Create predefined planet.
/*!
 * This function creates a predefined planet.
 * \param predefinedPlanet Name of predefined planet.
 * \return CelestialBody pointer.
 */
CelestialBody* createPredefinedPlanet(
        PredefinedPlanets predefinedPlanet );

}

#endif // PREDEFINED_PLANETS_H

// End of file.
