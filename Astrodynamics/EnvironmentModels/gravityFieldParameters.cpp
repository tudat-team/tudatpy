/*! \file gravityFieldParameters.cpp
 *    This file contains a class that define a gravity field
 *    with all their characteristics.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
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
 *    Last modified     : 29 September, 2010
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modifiation is unlawful and
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

// Include statements.
#include "gravityFieldParameters.h"

//! Default constructor.
GravityFieldParameters::GravityFieldParameters( )
{
    // Set default value for gravitational parameter.
    gravitationalParameter_ = -0.0;

    // Set default value for reference radius.
    referenceRadius_        = -0.0;
}

//! Customized constructor.
GravityFieldParameters::GravityFieldParameters(
    const double& gravitationalParameter,
    const double& referenceRadius )
{
    // Set gravitational parameter.
    gravitationalParameter_ = gravitationalParameter;

    // Set reference radius.
    referenceRadius_ = referenceRadius;
}

//! Sets the gravitational parameter.
void GravityFieldParameters::
        setGravitationalParameter( const double& gravitationalParameter )
{
    // Set gravitational parameter.
    gravitationalParameter_ = gravitationalParameter;
}

//! Gets the gravitational parameter.
const double GravityFieldParameters::getGravitationalParameter( ) const
{
    // Return gravitational parameter.
    return gravitationalParameter_;
}

//! Sets the reference radius.
void GravityFieldParameters::
        setReferenceRadius( const double& referenceRadius )
{
    // Return reference radius.
    referenceRadius_ = referenceRadius;
}

//! Gets the reference radius.
const double GravityFieldParameters::getReferenceRadius( ) const
{
    // Return reference radius.
    return referenceRadius_;
}

// End of file.
