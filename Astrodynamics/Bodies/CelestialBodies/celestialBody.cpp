/*! \file celestialBody.cpp
 *    This file contains a class that describes celestial bodies
 *    with all their characteristics.
 *
 *    Path              : /Astrodynamics/Environment/CelestialBodies/
 *    Version           : 3
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
 *    Last modified     : 29 September, 2010
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
 *      100906    J. Melman     First creation of code.
 *      100910    J. Melman     No more switch statement and enum.
 *      100929    K. Kumar      Minor comment changes.
 *
 */

// Include statements.
#include "celestialBody.h"

//! Default constructor.
CelestialBody::CelestialBody( )
{
}

//! Customized constructor.
CelestialBody::CelestialBody( const GravityFieldParameters&
                              gravityFieldParameters )
{
    // Set gravity field parameters.
    gravityFieldParameters_ = gravityFieldParameters;
}

//! Sets the gravity field parameters.
void CelestialBody::
        setGravityFieldParameters( const GravityFieldParameters&
                                   gravityFieldParameters )
{
    // Set gravity field parameters.
    gravityFieldParameters_ = gravityFieldParameters;
}

//! Gets the gravitational parameter.
const double CelestialBody::getGravitationalParameter( ) const
{
    // Return gravitational parameter.
    return gravityFieldParameters_.getGravitationalParameter( );
}

//! Gets the reference radius.
const double CelestialBody::getReferenceRadius( ) const
{
    // Return reference radius.
    return gravityFieldParameters_.getReferenceRadius( );
}

// End of file.
