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

#include <iostream>

#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"
#include "Tudat/SimulationSetup/createRotationModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a rotation model.
boost::shared_ptr< ephemerides::RotationalEphemeris > createRotationModel(
        const boost::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body )
{
    using namespace tudat::ephemerides;

    // Declare return object.
    boost::shared_ptr< RotationalEphemeris > rotationalEphemeris;

    // Check which type of rotation model is to be created.
    switch( rotationModelSettings->getRotationType( ) )
    {
    case simple_rotation_model:
    {
        // Check whether settings for simple rotation model are consistent with its type.
        boost::shared_ptr< SimpleRotationModelSettings > simpleRotationSettings =
                boost::dynamic_pointer_cast< SimpleRotationModelSettings >( rotationModelSettings );
        if( simpleRotationSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected simple rotation model settings for " + body );
        }
        else
        {
            // Create and initialize simple rotation model.
            rotationalEphemeris = boost::make_shared< SimpleRotationalEphemeris >(
                        simpleRotationSettings->getInitialOrientation( ),
                        simpleRotationSettings->getRotationRate( ),
                        simpleRotationSettings->getInitialTime( ),
                        basic_astrodynamics::JULIAN_DAY_ON_J2000,
                        simpleRotationSettings->getOriginalFrame( ),
                        simpleRotationSettings->getTargetFrame( ) );
        }
        break;
    }
    case spice_rotation_model:
    {
        // Create rotational ephemeris directly from Spice.
        rotationalEphemeris = boost::make_shared< SpiceRotationalEphemeris >(
                    rotationModelSettings->getOriginalFrame( ),
                    rotationModelSettings->getTargetFrame( ) );
        break;
    }
    default:
        throw std::runtime_error(
                 "Error, did not recognize rotation model settings type " +
                  boost::lexical_cast< std::string >( rotationModelSettings->getRotationType( ) ) );
    }

    return rotationalEphemeris;
}

}

}
