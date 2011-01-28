/*! \file orbitalElementConversions.h
 *    This header file contains a namespace with orbital element conversion
 *    functions.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 20 October, 2010
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
 *      101020    K. Kumar      First creation of code.
 *      101025    E. Iorfida    Fulfillment of the code with gravitational
 *                              parameter.
 *      101130    E. Iorfida    Gravitational parameter removed.
 *      101202    J. Melman     Replaced #endif statement and changed Doxygen
 *                              return statement.
 *      101203    E. Iorfida    Gravitational parameter added.
 *      101219    J. Melman     Doxygen comment on gravitational parameter
 *                              added.
 *      110128    K. Kumar      Changed references to pointers for functions.
 */

// Include statements.
#include "cartesianElements.h"
#include "keplerianElements.h"
#include "linearAlgebra.h"
#include "basicMathematicsFunctions.h"
#include "celestialBody.h"

#ifndef ORBITALELEMENTCONVERSIONS_H
#define ORBITALELEMENTCONVERSIONS_H

//! Orbital element conversions namespace.
/*!
 *  Orbital element conversions namespace.
 */
namespace orbital_element_conversions
{

//! Convert Keplerian to Cartesian orbital elements.
/*!
 * This function converts Keplerian to Cartesian orbital elements.
 * \param pointerToKeplerianElements Pointer to KeplerianElements object.
 * \param pointerToCelestialBody Pointer to CelestialBody object.
 * \return Pointer to CartesianElements object.
 */
CartesianElements* convertKeplerianToCartesianElements(
        KeplerianElements* pointerToKeplerianElements,
        CelestialBody* pointerToCelestialBody );

//! Convert Cartesian to Keplerian orbital elements.
/*!
 * This function converts Cartesian to Keplerian orbital elements.
 * \param pointerToCartesianElements Pointer to CartesianElements object.
 * \param pointerToCelestialBody Pointer to CelestialBody object.
 * \return Pointer to KeplerianElements object.
 */
KeplerianElements* convertCartesianToKeplerianElements(
        CartesianElements* pointerToCartesianElements,
        CelestialBody* pointerToCelestialBody );

}

#endif // ORBITALELEMENTCONVERSIONS_H

// End of file.
