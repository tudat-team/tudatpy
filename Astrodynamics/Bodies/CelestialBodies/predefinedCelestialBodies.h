/*! \file predefinedCelestialBodies.h
 *    Header file that defines a namespace that contains predefined celestial
 *    bodies in Tudat.
 *
 *    Path              : /Astrodynamics/Bodies/CelestialBodies/
 *    Version           : 7
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
 *    Last modified     : 21 April, 2011
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
 *      YYMMDD    Author            Comment
 *      100920    J. Melman         First creation of code.
 *      100929    K. Kumar          Minor comment changes.
 *      110112    K. Kumar          Updated code to use enum and
 *                                  createPredefinedCelestialBody() function.
 *      110114    J. Melman         Corrected Doxygen comments.
 *      110128    K. Kumar          Added Mars.
 *      110310    K. Kumar          Added Sun; changed filename and name of
 *                                  namespace.
 *      110421    E. Iorfida        Added Jupiter and Venus.
 */

#ifndef PREDEFINED_CELESTIAL_BODIES_H
#define PREDEFINED_CELESTIAL_BODIES_H

// Include statements.
#include "approximatePlanetPositions.h"
#include "celestialBody.h"
#include "ephemeris.h"
#include "sphericalHarmonicsGravityField.h"
#include "predefinedGravityFieldModels.h"

// Using declarations.
using std::cerr;
using std::endl;

//! Predefined celestial bodies namespace.
/*!
 * Predefined celestial bodies namespace.
 */
namespace predefined_celestial_bodies
{

//! Predefined celestial bodies.
/*!
 * Predefined celestial bodies.
 */
enum PredefinedCelestialBodies
{
    earth,
    mars,
    sun,
    jupiter,
    venus
};

//! Create predefined celestial body.
/*!
 * This function creates a predefined celestial.
 * \param predefinedCelestialBody Name of predefined celestial.
 * \return CelestialBody pointer.
 */
CelestialBody* createPredefinedCelestialBody(
        PredefinedCelestialBodies predefinedCelestialBody );

}

#endif // PREDEFINED_CELESTIAL_BODIES_H

// End of file.
