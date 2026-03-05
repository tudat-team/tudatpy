/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEBODIESSETTINGS_H
#define TUDAT_CREATEBODIESSETTINGS_H

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/simulation/environment_setup/createAerodynamicCoefficientInterface.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/simulation/environment_setup/createBodyShapeModel.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createRadiationPressureInterface.h"
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/createRadiationSourceModel.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"

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
struct BodySettings {
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
    // RP-OLD
    std::map< std::string, std::shared_ptr< RadiationPressureInterfaceSettings > > radiationPressureSettings;

    //! Settings for the radiation source model that the body is to contain.
    std::shared_ptr< RadiationSourceModelSettings > radiationSourceModelSettings;

    //! Settings for the radiation pressure target model that the body is to contain.
    std::shared_ptr< RadiationPressureTargetModelSettings > radiationPressureTargetModelSettings;

    //! Settings for the aerodynamic coefficients that the body is to contain.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;

    std::shared_ptr< RigidBodyPropertiesSettings > rigidBodyPropertiesSettings;

    std::shared_ptr< FullPanelledBodySettings > bodyExteriorPanelSettings_;

    //! Settings for variations of the gravity field of the body.
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;

    std::vector< std::shared_ptr< BodyDeformationSettings > > bodyDeformationSettings;

    std::vector< std::shared_ptr< GroundStationSettings > > groundStationSettings;
};

class BodyListSettings
{
public:
    BodyListSettings( const std::string frameOrigin = "SSB", const std::string frameOrientation = "ECLIPJ2000" ):
        bodySettings_( std::map< std::string, std::shared_ptr< BodySettings > >( ) ), frameOrigin_( frameOrigin ),
        frameOrientation_( frameOrientation )
    { }

    BodyListSettings( const std::map< std::string, std::shared_ptr< BodySettings > >& bodySettings,
                      const std::string frameOrigin = "SSB",
                      const std::string frameOrientation = "ECLIPJ2000" ):
        bodySettings_( bodySettings ), frameOrigin_( frameOrigin ), frameOrientation_( frameOrientation )
    { }

    std::shared_ptr< BodySettings > at( const std::string& bodyName ) const
    {
        if( bodySettings_.count( bodyName ) == 0 )
        {
            throw std::runtime_error( "Error when retrieving body settings for " + bodyName + ", no such body exist in this object." );
        }

        return bodySettings_.at( bodyName );
    }

    std::shared_ptr< BodySettings > get( const std::string& bodyName ) const
    {
        return at( bodyName );
    }

    int count( const std::string& bodyName ) const
    {
        return bodySettings_.count( bodyName );
    }

    void addSettings( std::shared_ptr< BodySettings > settingsToAdd, const std::string bodyName )
    {
        bodySettings_[ bodyName ] = settingsToAdd;
    }

    void addSettings( const std::string bodyName )
    {
        bodySettings_[ bodyName ] = std::make_shared< BodySettings >( );
    }

    void clear( )
    {
        bodySettings_.clear( );
    }

    void resetFrames( const std::string frameOrigin = "SSB", const std::string frameOrientation = "ECLIPJ2000" )
    {
        frameOrigin_ = frameOrigin;
        frameOrientation_ = frameOrientation;
    }

    std::string getFrameOrigin( ) const
    {
        return frameOrigin_;
    }

    std::string getFrameOrientation( ) const
    {
        return frameOrientation_;
    }

    std::map< std::string, std::shared_ptr< BodySettings > > getMap( ) const
    {
        return bodySettings_;
    }

private:
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings_;

    std::string frameOrigin_;

    std::string frameOrientation_;
};

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATEBODIESSETTINGS_H
