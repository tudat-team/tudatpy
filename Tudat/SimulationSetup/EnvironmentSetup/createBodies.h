/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodyShapeModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRadiationPressureInterface.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createFlightConditions.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"

namespace tudat
{

namespace simulation_setup
{

//! Struct holding settings for a body to be created.
/*!
 *  Struct holding settings for a body to be created. From the settings, a Body object is
 *  created by the createBodies function. Default values can be generated from the function in
 *  defaultBodies.h.
 */
struct BodySettings
{
    //! Constant mass.
    double constantMass = TUDAT_NAN;

    //! Settings for the atmosphere model that the body is to contain.
    std::shared_ptr< AtmosphereSettings > atmosphereSettings;

    //! Settings for the ephemeris model that the body is to contain.
    std::shared_ptr< EphemerisSettings > ephemerisSettings;

    //! Settings for the gravity field model that the body is to contain.
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings;

    //! Settings for the rotation model that the body is to contain.
    std::shared_ptr< RotationModelSettings > rotationModelSettings;

    //! Settings for the shape model that the body is to contain.
    std::shared_ptr< BodyShapeSettings > shapeModelSettings;

    //! Settings for the radiations pressure interfaces that the body is to contain (source body as key).
    std::map< std::string,
    std::shared_ptr< RadiationPressureInterfaceSettings > > radiationPressureSettings;

    //! Settings for the aerodynamic coefficients that the body is to contain.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;

    //! Settings for variations of the gravity field of the body.
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;

    std::vector< std::shared_ptr< GroundStationSettings > > groundStationSettings;

};

class BodyListSettings
{
public:

    BodyListSettings( const std::string frameOrigin = "SSB", const std::string frameOrientation = "ECLIPJ2000" ):
        _bodySettings_( std::map< std::string, std::shared_ptr< BodySettings > >( ) ),
        frameOrigin_( frameOrigin ), frameOrientation_( frameOrientation ){ }

    BodyListSettings( const std::map< std::string, std::shared_ptr< BodySettings > >& bodySettings,
                      const std::string frameOrigin = "SSB", const std::string frameOrientation = "ECLIPJ2000" ):
        _bodySettings_( bodySettings ), frameOrigin_( frameOrigin ), frameOrientation_( frameOrientation ){ }

    std::shared_ptr< BodySettings > at( const std::string& bodyName ) const
    {
        return _bodySettings_.at( bodyName );
    }

    void addSettings( std::shared_ptr< BodySettings > settingsToAdd, const std::string bodyName )
    {
        _bodySettings_[ bodyName ] = settingsToAdd;
    }

    void addSettings( const std::string bodyName )
    {
        _bodySettings_[ bodyName ] = std::make_shared< BodySettings >( );
    }

    std::string getFrameOrigin( ) const { return frameOrigin_; }

    std::string getFrameOrientation( ) const { return frameOrientation_; }

    std::map< std::string, std::shared_ptr< BodySettings > > get( ) const { return _bodySettings_; }


private:

    std::map< std::string, std::shared_ptr< BodySettings > > _bodySettings_;

    std::string frameOrigin_;

    std::string frameOrientation_;
};


//! Function that determines the order in which bodies are to be created
/*!
 * Function that determines the order in which bodies are to be created, to ensure that any dependency between body models are
 * correctly handled. Currently, no dependencies exist that force any particular creation order.
 * \param bodySettings Map of body settings (with map key the body name)
 * \return List of pairs: name and body settings of that body
 */
std::vector< std::pair< std::string, std::shared_ptr< BodySettings > > > determineBodyCreationOrder(
        const std::map< std::string, std::shared_ptr< BodySettings > >& bodySettings );

//! Function to create a map of bodies objects.
/*!
 *  Function to create a msap of body objects based on model-specific settings for the bodies,
 *  containing settings for each relevant environment model.
 *  \param bodySettings List of settings for the bodies that are to be created, defined as a map of
 *  pointers to an object of class BodySettings
 *  \return List of bodies created according to settings in bodySettings.
 */
NamedBodyMap createBodies(
        const BodyListSettings& bodySettings );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEBODIES_H
