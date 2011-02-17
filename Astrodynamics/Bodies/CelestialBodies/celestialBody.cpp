/*! \file celestialBody.cpp
 *    This file contains a class that describes celestial bodies
 *    with all their characteristics.
 *
 *    Path              : /Astrodynamics/Bodies/CelestialBodies/
 *    Version           : 4
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
 *    Last modified     : 13 January, 2011
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
 *      110113    K. Kumar      Added getGravityFieldModel() function;
 *                              updated path.
 */

// Include statements.
#include "celestialBody.h"

// Using declarations.
using std::endl;

//! Default constructor.
CelestialBody::CelestialBody( )
{
}

//! Set gravity field model.
void CelestialBody::setGravityFieldModel( GravityFieldModel*
                                          pointerToGravityFieldModel)

{
    // Set gravity field model.
    pointerToGravityFieldModel_ = pointerToGravityFieldModel;
}

//! Get gravitational parameter.
const double CelestialBody::getGravitationalParameter( ) const
{
    // Return gravitational parameter.
    return pointerToGravityFieldModel_->getGravitationalParameter( );
}

//! Get gravity field model.
GravityFieldModel* CelestialBody::getGravityFieldModel( )
 {
     return pointerToGravityFieldModel_;
 }

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          CelestialBody* pointerToCelestialBody )
{
    stream << "This is a CelestialBody object." << endl;
    stream << "The gravitational parameter is set to: "
           << pointerToCelestialBody->getGravitationalParameter( )
           << endl;

    // Return stream.
    return stream;
}

// End of file.
