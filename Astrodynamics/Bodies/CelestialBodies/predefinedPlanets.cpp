/*! \file predefinedPlanets.cpp
 *    Source file that defines a namespace that contains predefined planets in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Bodies/CelestialBodies/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 20 September, 2010
 *    Last modified     : 12 January, 2011
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
 *      110112    K. Kumar     First creation of code.
 */

// Include statements.
#include "predefinedPlanets.h"

//! Predefined planets namespace.
namespace predefined_planets
{

//! Create predefined planet.
CelestialBody* createPredefinedPlanet( PredefinedPlanets predefinedPlanet )
{
    // Declare local variables.
    // Declare pointer to a celestial body.
   CelestialBody* pointerToCelestialBody_;
   // Declare pointer to a gravity field model object.
   GravityFieldModel* pointerToGravityFieldModel_;

   // Create celestial body object.
   pointerToCelestialBody_ = new CelestialBody;

   // Select predefined planet on input.
   switch( predefinedPlanet )
   {
   case earth:

       // Set pointer to gravity field model to predefined Earth gravity field.
       pointerToGravityFieldModel_ =
               predefined_gravity_field_models
               ::createPredefinedCentralGravityField(
                       predefined_gravity_field_models::earth );

       // Set spherical harmonics gravity field.
       pointerToCelestialBody_->setGravityFieldModel( pointerToGravityFieldModel_ );
   }

   // Return pointer to a spherical harmonics gravity field model.
   return pointerToCelestialBody_;
}

}

// End of file.
