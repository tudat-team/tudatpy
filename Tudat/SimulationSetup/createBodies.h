/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEBODIES_H
#define TUDAT_CREATEBODIES_H

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/SimulationSetup/createEphemeris.h"
#include "Tudat/SimulationSetup/createAtmosphereModel.h"
#include "Tudat/SimulationSetup/createBodyShapeModel.h"
#include "Tudat/SimulationSetup/createEphemeris.h"
#include "Tudat/SimulationSetup/createGravityField.h"
#include "Tudat/SimulationSetup/createRotationModel.h"
#include "Tudat/SimulationSetup/createRadiationPressureInterface.h"
#include "Tudat/SimulationSetup/createFlightConditions.h"

namespace tudat
{

namespace simulation_setup
{

//! Struct holding settings for a body to be created.
/*!
 *  Struct holding settings for a body to be created. From the settings, a CelestialBody object is
 *  created by the createBodies function. Default values can be generated from the function in
 *  defaultBodies.h.
 */
struct BodySettings
{
    //! Settings for the atmosphere model that the body is to contain.
    boost::shared_ptr< AtmosphereSettings > atmosphereSettings;

    //! Settings for the ephemeris model that the body is to contain.
    boost::shared_ptr< EphemerisSettings > ephemerisSettings;

    //! Settings for the gravity field model that the body is to contain.
    boost::shared_ptr< GravityFieldSettings > gravityFieldSettings;

    //! Settings for the rotation model that the body is to contain.
    boost::shared_ptr< RotationModelSettings > rotationModelSettings;

    //! Settings for the shape model that the body is to contain.
    boost::shared_ptr< BodyShapeSettings > shapeModelSettings;

    //! Settings for the radiations pressure interfaces that the body is to contain (source body as key).
    std::map< std::string,
              boost::shared_ptr< RadiationPressureInterfaceSettings > > radiationPressureSettings;

    //! Settings for the aerodynamic coefficients that the body is to contain.
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;

    std::vector< boost::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;

};

std::vector< std::pair< std::string, boost::shared_ptr< BodySettings > > > determineBodyCreationOrder(
        const std::map< std::string, boost::shared_ptr< BodySettings > >& bodySettings );

//! Function to create a map of bodies objects.
/*!
 *  Function to create a msap of body objects based on model-specific settings for the bodies,
 *  containing settings for each relevant environment model.
 *  \param bodySettings List of settings for the bodies that are to be created, defined as a map of
 *  pointers to an object of class BodySettings
 *  \return List of bodies created according to settings in bodySettings.
 */
NamedBodyMap createBodies(
        const std::map< std::string, boost::shared_ptr< BodySettings > >& bodySettings );


//! Function to define the global origin and orientation of the reference frame
/*!
 * Function to define the global origin and orientation of the reference frame that is to be used in
 * the simulations.  This function checks the origin and orientation of the Ephemeris and
 * RotationalEphemeris, and checks whether their origin/orientation is the same as that
 * globalFrameOrigin and globalFrameOrientation provided as input.  In particular, this function
 * sets the ephemerisFrameToBaseFrameFunction_ anf ephemerisFrameToBaseFrameLongFunction_ variables
 * of the Body objects, which provide a time-dependent translation of the global origin to the
 * body's ephemeris origin. In case of an inconsistency in the current and requried frames, this
 * function throws an error.
 * \param bodyMap List of body objects that constitute the environment.
 * \param globalFrameOrigin Global reference frame origin.
 * \param globalFrameOrientation Global referencef frame orientation.
 */
void setGlobalFrameBodyEphemerides( const NamedBodyMap& bodyMap,
                                    const std::string& globalFrameOrigin,
                                    const std::string& globalFrameOrientation );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEBODIES_H
