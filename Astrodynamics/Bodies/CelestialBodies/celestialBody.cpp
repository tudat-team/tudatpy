/*! \file celestialBody.cpp
 *    This file contains a class that describes celestial bodies
 *    with all their characteristics.
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
 *    Date created      : 6 September, 2010
 *    Last modified     : 10 March, 2011
 *
 *    References
 *
 *    Notes
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
 *      100906    J. Melman         First creation of code.
 *      100910    J. Melman         No more switch statement and enum.
 *      100929    K. Kumar          Minor comment changes.
 *      110113    K. Kumar          Added getGravityFieldModel() function;
 *                                  updated path.
 *      110310    K. Kumar          Added ephemeris; added missing destructor.
 */

// Include statements.
#include "Astrodynamics/Bodies/CelestialBodies/celestialBody.h"

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          CelestialBody& celestialBody )
{
    // Using declarations.
    using std::endl;

    stream << "This is a CelestialBody object." << endl;
    stream << "The gravitational parameter is set to: "
           << celestialBody.getGravitationalParameter( ) << endl;

    // Return stream.
    return stream;
}

// End of file.
