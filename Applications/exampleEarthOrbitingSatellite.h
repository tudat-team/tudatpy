/*! \file exampleEarthOrbitingSatellite.h
 *    Header file that contains an example application of an Earth-orbiting
 *    satellite mission scenario using Tudat.
 *
 *    Path              : /Applications/
 *    Version           : 1
 *    Check status      : Unchecked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : B. Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : bart.romgens@gmail.com
 *
 *    Date created      : 16 February, 2011
 *    Last modified     : 16 February, 2011
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
 *      YYMMDD    Author              Comment
 *      110216    K. Kumar            File created.
 */

#ifndef EXAMPLEEARTHORBITINGSATELLITE_H
#define EXAMPLEEARTHORBITINGSATELLITE_H

// Include statements.
#include <cmath>
#include <iostream>
#include "cartesianElements.h"
#include "celestialBody.h"
#include "gravity.h"
#include "numericalPropagator.h"
#include "predefinedPlanets.h"
#include "rungeKutta4thOrderFixedStepsize.h"
#include "unitConversions.h"
#include "vehicle.h"
#include "writingOutputToFile.h"

//! Namespace for all applications.
/*!
 * Namespace containing all applications.
 */
namespace applications
{

//! Example of an Earth-orbiting satellite.
/*!
 * Example of an Earth-orbiting satellite.
 */
void exampleEarthOrbitingSatellite( );

}

#endif // EXAMPLEEARTHORBITINGSATELLITE_H

// End of file.
