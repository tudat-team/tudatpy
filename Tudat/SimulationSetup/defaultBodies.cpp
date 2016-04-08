/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/SimulationSetup/defaultBodies.h"

namespace tudat
{
namespace simulation_setup
{

//! Function to create default settings for a body's atmosphere model.
boost::shared_ptr< AtmosphereSettings > getDefaultAtmosphereModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime )
{

    boost::shared_ptr< AtmosphereSettings > atmosphereSettings;

    // A default atmosphere is only implemented for Earth.
    if( bodyName == "Earth" )
    {
        atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
                    input_output::getTudatRootPath( ) + "/External/AtmosphereTables/" +
                    "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );
    }


    return atmosphereSettings;
}

//! Function to create default settings for a body's ephemeris.
boost::shared_ptr< EphemerisSettings > getDefaultEphemerisSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime )
{
    // Create settings for an interpolated Spice ephemeris.
    return boost::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialTime, finalTime, 300.0, "SSB", "ECLIPJ2000" );
}

//! Function to create default settings for a body's gravity field model.
boost::shared_ptr< GravityFieldSettings > getDefaultGravityFieldSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime )
{
    // Create settings for a point mass gravity with data from Spice
    return boost::make_shared< GravityFieldSettings >( central_spice );
}

//! Function to create default settings from which to create a single body object.
boost::shared_ptr< RotationModelSettings > getDefaultRotationModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime )
{
    // Create settings for a rotation model taken directly from Spice.
    return boost::make_shared< RotationModelSettings >(
                spice_rotation_model, "ECLIPJ2000", "IAU_" + bodyName );
}

//! Function to create default settings for a body's shape model.
boost::shared_ptr< BodyShapeSettings > getDefaultBodyShapeSettings(
        const std::string& body,
        const double initialTime, const double finalTime )
{
    return boost::make_shared< SphericalBodyShapeSettings >(
                spice_interface::getAverageRadius( body ) );
}

//! Function to create default settings for a body's rotation model.
boost::shared_ptr< BodySettings > getDefaultSingleBodySettings(
        const std::string& body,
        const double initialTime,
        const double finalTime )
{
    boost::shared_ptr< BodySettings > singleBodySettings = boost::make_shared< BodySettings >( );

    // Get default settings for each of the environment models in the body.
    singleBodySettings->atmosphereSettings = getDefaultAtmosphereModelSettings(
                body, initialTime, finalTime );
    singleBodySettings->rotationModelSettings = getDefaultRotationModelSettings(
                body, initialTime, finalTime );
    singleBodySettings->ephemerisSettings = getDefaultEphemerisSettings(
                body, initialTime, finalTime );
    singleBodySettings->gravityFieldSettings = getDefaultGravityFieldSettings(
                body, initialTime, finalTime );
    singleBodySettings->shapeModelSettings = getDefaultBodyShapeSettings(
                body, initialTime, finalTime );

    return singleBodySettings;
}


//! Function to create default settings from which to create a set of body objects.
std::map< std::string, boost::shared_ptr< BodySettings > > getDefaultBodySettings(
        const std::vector< std::string >& bodies,
        const double initialTime,
        const double finalTime )
{
    std::map< std::string, boost::shared_ptr< BodySettings > > settingsMap;

    // Iterative over all bodies and get default settings.
    for( unsigned int i = 0; i < bodies.size( ); i++ )
    {
        settingsMap[ bodies.at( i ) ] = getDefaultSingleBodySettings(
                    bodies.at( i ), initialTime, finalTime );

    }
    return settingsMap;
}

}

}
