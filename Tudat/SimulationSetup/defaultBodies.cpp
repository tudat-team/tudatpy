/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      150501    D. Dirkx          Ported from personal code
 *
 *    References
 *
 *    Notes
 *
 */

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
