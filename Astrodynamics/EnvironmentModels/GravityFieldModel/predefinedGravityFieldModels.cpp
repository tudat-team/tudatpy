/*! \file predefinedGravityFieldModels.cpp
 *    Source file that defines a namespace containing pre-defined gravity field
 *    models.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
 *    Version           : 5
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
 *    Last modified     : 21 April, 2011
 *
 *    References
 *       http://solarsystem.nasa.gov/planets/index.cfm
 *       http://en.wikipedia.org/wiki/Standard_gravitational_parameter
 *
 *    Notes
 *      The reference radii used in this code, are the "mean" radii of the planets
 *      on NASA website. They should be changed with better values of reference
 *      radii for the gravity fields.
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
 *      110106    K. Kumar          File created.
 *      110107    K. Kumar          Updated code after comments by J. Melman.
 *      110128    K. Kumar          Added Mars central gravity field.
 *      110310    K. Kumar          Added Sun central gravity field.
 *      110421    E. Iorfida        Added Jupiter and Venus central gravity fields.
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
                ->setGravitationalParameter( 398600.4418e9  );

        // Set reference radius.
        pointerToSphericalHarmonicsGravityField
                ->setReferenceRadius( 6371.0e3 );

        // Set degree of expansion.
        pointerToSphericalHarmonicsGravityField->setDegreeOfExpansion( 0 );

        // Set order of expansion.
        pointerToSphericalHarmonicsGravityField->setOrderOfExpansion( 0 );

        break;

   case mars:

       // Set gravitational parameter.
       pointerToSphericalHarmonicsGravityField
               ->setGravitationalParameter( 42828.0e9  );

       // Set reference radius.
       pointerToSphericalHarmonicsGravityField
               ->setReferenceRadius( 3389.5e3 );

       // Set degree of expansion.
       pointerToSphericalHarmonicsGravityField->setDegreeOfExpansion( 0 );

       // Set order of expansion.
       pointerToSphericalHarmonicsGravityField->setOrderOfExpansion( 0 );

       break;

   case sun:

       // Set gravitational parameter.
       // Reference: http://ssd.jpl.nasa.gov/?constants#ref
       pointerToSphericalHarmonicsGravityField
               ->setGravitationalParameter( 1.32712440018e20 );

       // Set reference radius.
       pointerToSphericalHarmonicsGravityField
               ->setReferenceRadius( 695508.0e3 );

       // Set degree of expansion.
       pointerToSphericalHarmonicsGravityField->setDegreeOfExpansion( 0 );

       // Set order of expansion.
       pointerToSphericalHarmonicsGravityField->setOrderOfExpansion( 0 );

       break;

   case jupiter:

       // Set gravitational parameter.
       pointerToSphericalHarmonicsGravityField
               ->setGravitationalParameter( 126686534.0e9  );

       // Set reference radius.
       pointerToSphericalHarmonicsGravityField
               ->setReferenceRadius( 69911.0e3 );

       // Set degree of expansion.
       pointerToSphericalHarmonicsGravityField->setDegreeOfExpansion( 0 );

       // Set order of expansion.
       pointerToSphericalHarmonicsGravityField->setOrderOfExpansion( 0 );

       break;

   case venus:

       // Set gravitational parameter.
       pointerToSphericalHarmonicsGravityField
               ->setGravitationalParameter( 324859.0e9  );

       // Set reference radius.
       pointerToSphericalHarmonicsGravityField
               ->setReferenceRadius( 6051.8e3 );

       // Set degree of expansion.
       pointerToSphericalHarmonicsGravityField->setDegreeOfExpansion( 0 );

       // Set order of expansion.
       pointerToSphericalHarmonicsGravityField->setOrderOfExpansion( 0 );

       break;

   default:

       // Print cerr statement.
       cerr << "Desired predefined central gravity field does not exist."
            << endl;
   }

   // Return pointer to a spherical harmonics gravity field model.
   return pointerToSphericalHarmonicsGravityField;
}

}
// End of file.
