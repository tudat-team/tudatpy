/*! \file predefinedGravityFieldModels.h
 *    Header file that defines a namespace containing pre-defined gravity field
 *    models.
 *
 *    Path              : /Astrodynamics/EnvironmentModel/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 16 November, 2010
 *    Last modified     : 15 January, 2011
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
 *      YYMMDD    author              comment
 *      101116    K. Kumar            File created.
 *      101208    K. Kumar            Added predefined Earth central body
 *                                    gravity field.
 *      101215    K. Kumar            Removed createPredefinedPlanet function;
 *                                    not functioning correctly. Added
 *                                    predefined Earth central body gravity
 *                                    field.
 *      110106    K. Kumar            Created new functions and split between
 *                                    .h and .cpp files.
 *      110107    J. Melman           Changed centralBodyGravityField into
 *                                    centralGravityField.
 *      110107    K. Kumar            Change enum and updated code.
 *      110114    J. Melman           Corrected Doxygen comments.
 */

#ifndef PREDEFINEDGRAVITYFIELDMODELS_H
#define PREDEFINEDGRAVITYFIELDMODELS_H

// Include statements.
#include "gravityFieldModel.h"
#include "sphericalHarmonicsGravityField.h"

//! Pre-defined gravity field models namespace.
/*!
 * Pre-defined gravity field models namespace.
 */
namespace predefined_gravity_field_models
{

//! Bodies with a predefined gravity field.
/*!
 * Bodies with a predefined gravity field.
*/
enum BodiesWithPredefinedGravityFields
{
    earth
};

//! Create predefined central gravity field.
/*!
 * This function creates a predefined central gravity field.
 * \param bodyWithPredefinedCentralGravityField Name of body with
 *          predefined central gravity field.
 * \return SphericalHarmonicsGravityField pointer.
 */
SphericalHarmonicsGravityField* createPredefinedCentralGravityField(
        BodiesWithPredefinedGravityFields
        bodyWithPredefinedCentralGravityField );

//! Create predefined spherical harmonics gravity field.
/*!
 * This function creates a predefined spherical harmonics gravity field.
 * \param bodyWithPredefinedSphericalHarmonicsGravityField Name of body with
 *          predefined spherical harmonics gravity field.
 * \return SphericalHarmonicsGravityField pointer.
 */
SphericalHarmonicsGravityField* createPredefinedSphericalHarmonicsGravityField(
        BodiesWithPredefinedGravityFields
        bodyWithPredefinedSphericalHarmonicsGravityField );

}

#endif // PREDEFINEDGRAVITYFIELDMODELS_H

// End of file.
