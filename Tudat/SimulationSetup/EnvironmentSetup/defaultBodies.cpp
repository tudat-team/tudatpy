/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#if USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#endif

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"

namespace tudat
{
namespace simulation_setup
{

//! Function to create default settings for a body's atmosphere model.
std::shared_ptr< AtmosphereSettings > getDefaultAtmosphereModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime )
{

    std::shared_ptr< AtmosphereSettings > atmosphereSettings;

    // A default atmosphere is only implemented for Earth.
    if( bodyName == "Earth" )
    {
        std::string atmosphereTableFile = input_output::getAtmosphereTablesPath( ) + "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
        atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >( atmosphereTableFile );
    }

    return atmosphereSettings;
}

//! Function to create default settings for a body's ephemeris.
std::shared_ptr< EphemerisSettings > getDefaultEphemerisSettings(
        const std::string& bodyName )
{
#if USE_CSPICE
    // Create settings for an interpolated Spice ephemeris.
    return std::make_shared< DirectSpiceEphemerisSettings >(
                "SSB", "ECLIPJ2000", false, false, false );
#else
    throw std::runtime_error( "Default ephemeris settings can only be used together with the SPICE library" );
#endif
}

//! Function to create default settings for a body's ephemeris.
std::shared_ptr< EphemerisSettings > getDefaultEphemerisSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime,
        const double timeStep )
{
#if USE_CSPICE
    // Create settings for an interpolated Spice ephemeris.
    return std::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialTime, finalTime, timeStep, "SSB", "ECLIPJ2000" );
#else
    throw std::runtime_error( "Default ephemeris settings can only be used together with the SPICE library" );
#endif
}

//! Function to create default settings for a body's gravity field model.
std::shared_ptr< GravityFieldSettings > getDefaultGravityFieldSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime )
{
    if( bodyName == "Earth" )
    {
        return std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( egm96 );
    }
    else if( bodyName == "Moon" )
    {
        return std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( lpe200 );
    }
    else if( bodyName == "Mars" )
    {
        return std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    }
    else
    {
#if USE_CSPICE
        // Create settings for a point mass gravity with data from Spice
        return std::make_shared< GravityFieldSettings >( central_spice );
#else
        throw std::runtime_error( "Default gravity field settings can only be used together with the SPICE library" );
#endif
    }
}


//! Function to create default settings from which to create a single body object.
std::shared_ptr< RotationModelSettings > getDefaultRotationModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime )
{
    TUDAT_UNUSED_PARAMETER( initialTime );
    TUDAT_UNUSED_PARAMETER( finalTime );

#if USE_CSPICE
    // Create settings for a rotation model taken directly from Spice.
    return std::make_shared< RotationModelSettings >(
                spice_rotation_model, "ECLIPJ2000", "IAU_" + bodyName );
#else
    throw std::runtime_error( "Default rotational model settings can only be used together with the SPICE library" );
#endif
}

//! Function to create default settings for a body's shape model.
std::shared_ptr< BodyShapeSettings > getDefaultBodyShapeSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime )
{
    TUDAT_UNUSED_PARAMETER( initialTime );
    TUDAT_UNUSED_PARAMETER( finalTime );

#if USE_CSPICE
    return std::make_shared< SphericalBodyShapeSettings >(
                spice_interface::getAverageRadius( bodyName ) );
#else
    throw std::runtime_error( "Default bodyName settings can only be used together with the SPICE library" );
#endif
}

//! Function to create default settings for a body's rotation model.
std::shared_ptr< BodySettings > getDefaultSingleBodySettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime,
        const double timeStep )
{
    std::shared_ptr< BodySettings > singleBodySettings = std::make_shared< BodySettings >( );

    // Get default settings for each of the environment models in the body.
    singleBodySettings->atmosphereSettings = getDefaultAtmosphereModelSettings(
                bodyName, initialTime, finalTime );
    singleBodySettings->rotationModelSettings = getDefaultRotationModelSettings(
                bodyName, initialTime, finalTime );

    if( ( !( initialTime == initialTime ) && ( finalTime == finalTime ) ) ||
            ( ( initialTime == initialTime ) && !( finalTime == finalTime ) ) )
    {
        throw std::runtime_error( "Error when getting default body settings, only one input time is NaN" );
    }
    else if( !( initialTime == initialTime ) )
    {
        singleBodySettings->ephemerisSettings = getDefaultEphemerisSettings(
                    bodyName );
    }
    else
    {
        singleBodySettings->ephemerisSettings = getDefaultEphemerisSettings(
                    bodyName, initialTime, finalTime, timeStep );
    }
    singleBodySettings->gravityFieldSettings = getDefaultGravityFieldSettings(
                bodyName, initialTime, finalTime );
    singleBodySettings->shapeModelSettings = getDefaultBodyShapeSettings(
                bodyName, initialTime, finalTime );

    return singleBodySettings;
}


//! Function to create default settings from which to create a set of body objects.
std::map< std::string, std::shared_ptr< BodySettings > > getDefaultBodySettings(
        const std::vector< std::string >& bodies,
        const double initialTime,
        const double finalTime,
        const double timeStep )
{
    std::map< std::string, std::shared_ptr< BodySettings > > settingsMap;

    // Iterative over all bodies and get default settings.
    for( unsigned int i = 0; i < bodies.size( ); i++ )
    {
        settingsMap[ bodies.at( i ) ] = getDefaultSingleBodySettings(
                    bodies.at( i ), initialTime, finalTime, timeStep );

    }
    return settingsMap;
}

//! Function to create default settings from which to create a set of body objects, without stringent limitations on
//! time-interval of validity of environment.
std::map< std::string, std::shared_ptr< BodySettings > > getDefaultBodySettings(
        const std::vector< std::string >& bodies )
{
    std::map< std::string, std::shared_ptr< BodySettings > > settingsMap;

    // Iterative over all bodies and get default settings.
    for( unsigned int i = 0; i < bodies.size( ); i++ )
    {
        settingsMap[ bodies.at( i ) ] = getDefaultSingleBodySettings(
                    bodies.at( i ), TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );

    }
    return settingsMap;
}

} // namespace simulation_setup

} // namespace tudat
