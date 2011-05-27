/*! \file predefinedCelestialBodies.cpp
 *    Source file that defines a namespace that contains predefined celestial
 *    bodies Tudat.
 *
 *    Path              : /Astrodynamics/Bodies/CelestialBodies/
 *    Version           : 4
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
 *    Last modified     : 21 April, 2011
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
 *      YYMMDD    Author            Comment
 *      110112    K. Kumar          First creation of code.
 *      110128    K. Kumar          Added Mars; added default case with cerr.
 *      110310    K. Kumar          Added Sun; changed filename and name of
 *                                  namespace.
 *      110421    E. Iorfida        Added Jupiter and Venus.
 */

// Include statements.
#include "predefinedCelestialBodies.h"

//! Predefined celestial bodies namespace.
namespace predefined_celestial_bodies
{

//! Create predefined celestial body.
CelestialBody* createPredefinedCelestialBody( PredefinedCelestialBodies
                                              predefinedCelestialBody )
{
    // Declare local variables.
    // Declare pointer to a celestial body and create object.
   CelestialBody* pointerToCelestialBody_ = new CelestialBody;

   // Declare pointer to a gravity field model object.
   GravityFieldModel* pointerToGravityFieldModel_;

   // Declare pointer to a JPL approximate planet postions object and create
   // object.
   ApproximatePlanetPositions* pointerToApproximatePlanetPositions_
           = new ApproximatePlanetPositions;

   // Select predefined planet on input.
   switch( predefinedCelestialBody )
   {
   case earth:

       // Set pointer to gravity field model to predefined Earth gravity field.
       pointerToGravityFieldModel_ =
               predefined_gravity_field_models
               ::createPredefinedCentralGravityField(
                       predefined_gravity_field_models::earth );

       // Set Earth as planet for ephemeris.
       pointerToApproximatePlanetPositions_->setPlanet(
               ApproximatePlanetPositions::earthMoonBarycenter );

       break;

   case mars:

       // Set pointer to gravity field model to predefined Mars gravity field.
       pointerToGravityFieldModel_ =
               predefined_gravity_field_models
               ::createPredefinedCentralGravityField(
                       predefined_gravity_field_models::mars );

       // Set Mars as planet for ephemeris.
       pointerToApproximatePlanetPositions_->setPlanet(
               ApproximatePlanetPositions::mars );

       break;

   case sun:

       // Set pointer to gravity field model to predefined Mars gravity field.
       pointerToGravityFieldModel_ =
               predefined_gravity_field_models
               ::createPredefinedCentralGravityField(
                       predefined_gravity_field_models::sun );

       break;

   case jupiter:

          // Set pointer to gravity field model to predefined Jupiter gravity field.
          pointerToGravityFieldModel_ =
                  predefined_gravity_field_models
                  ::createPredefinedCentralGravityField(
                          predefined_gravity_field_models::jupiter );

          // Set Jupiter as planet for ephemeris.
          pointerToApproximatePlanetPositions_->setPlanet(
                  ApproximatePlanetPositions::jupiter );
          break;

   case venus:

          // Set pointer to gravity field model to predefined Venus gravity field.
          pointerToGravityFieldModel_ =
                  predefined_gravity_field_models
                  ::createPredefinedCentralGravityField(
                          predefined_gravity_field_models::venus );

          // Set Venus as planet for ephemeris.
          pointerToApproximatePlanetPositions_->setPlanet(
                  ApproximatePlanetPositions::venus );

          break;

   default:

       // Print cerr statement.
       cerr << "Desired predefined planet does not exist." << endl;
   }

   // Set spherical harmonics gravity field.
   pointerToCelestialBody_
           ->setGravityFieldModel( pointerToGravityFieldModel_ );

   // Set JPL approximate planet positions ephemeris in celestial body.
   pointerToCelestialBody_->setEphemeris( pointerToApproximatePlanetPositions_ );

   // Return pointer to a celestial body.
   return pointerToCelestialBody_;
}

}

// End of file.
