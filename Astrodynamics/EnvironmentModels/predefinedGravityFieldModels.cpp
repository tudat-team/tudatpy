/*! \file predefinedGravityFieldModels.cpp
 *    Source file that defines a namespace containing pre-defined gravity field
 *    models.
 *
 *    Path              : /Astrodynamics/EnvironmentModel/
 *    Version           : 2
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
 *    Date created      : 6 January, 2011
 *    Last modified     : 7 January, 2011
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
 *      110106    K. Kumar            File created.
 *      110107    K. Kumar            Updated code after comments by J. Melman.
 */

// Include statements.
#include "predefinedGravityFieldModels.h"

namespace predefined_gravity_field_models
{
    
//! Create predefined central gravity field.
SphericalHarmonicsGravityField* createPredefinedCentralGravityField(
        BodiesWithPredefinedGravityFields
        bodyWithPredefinedCentralGravityField )
{
    // Declare local variables.
    // Declare pointer to a spherical harmonics gravity field model.
   SphericalHarmonicsGravityField* pointerToSphericalHarmonicsGravityField;

   // Create spherical harmonics gravity field object.
   pointerToSphericalHarmonicsGravityField
           = new SphericalHarmonicsGravityField;

   // Select predefined central gravity field model based on input.
   switch( bodyWithPredefinedCentralGravityField )
   {
    case earth:

        // Set gravitational parameter.
        pointerToSphericalHarmonicsGravityField
                ->setGravitationalParameter( 398600.4415e9  );

        // Set reference radius.
        pointerToSphericalHarmonicsGravityField
                ->setReferenceRadius( 6378136.3 );

        // Set degree of expansion.
        pointerToSphericalHarmonicsGravityField->setDegreeOfExpansion( 0 );

        // Set order of expansion.
        pointerToSphericalHarmonicsGravityField->setOrderOfExpansion( 0 );
   }

   // Return pointer to a spherical harmonics gravity field model.
   return pointerToSphericalHarmonicsGravityField;
}

}
// End of file.
