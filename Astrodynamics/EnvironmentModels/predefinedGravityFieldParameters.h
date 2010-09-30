/*! \file predefinedGravityFieldParameters.h
 *    This file contains a class that predefines several gravity fields.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
 *    Version           : 2
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

#ifndef PREDEFINED_GRAVITY_FIELD_PARAMETERS_H
#define PREDEFINED_GRAVITY_FIELD_PARAMETERS_H

// Include statements.
#include "gravityFieldParameters.h"

//! Predefined gravity field parameters namespace.
/*!
 * Predefined gravity field parameters namespace.
 */
namespace predefined_gravity_field_parameters
{

// Earth gravity field parameters using object of GravityFieldParameters class.
const GravityFieldParameters defaultEarthGravityFieldParameters(
        398600.4415e9, 6378136.3 );

}

#endif // PREDEFINED_GRAVITY_FIELD_PARAMETERS_H

// End of file.
