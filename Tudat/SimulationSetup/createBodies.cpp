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

#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/SimulationSetup/createBodies.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;
using namespace gravitation;

//! Function to create a map of bodies objects.
NamedBodyMap createBodies(
        const std::map< std::string, boost::shared_ptr< BodySettings > >& bodySettings )
{
    // Declare map of bodies that is to be returned.
    NamedBodyMap bodyMap;

    // Create empty body objects.
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        bodyMap[ settingIterator->first ] = boost::make_shared< Body >( );
    }

    // Create ephemeris objects for each body (if required).
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        if( settingIterator->second->ephemerisSettings != NULL )
        {
            bodyMap[ settingIterator->first ]->setEphemeris(
                        createBodyEphemeris( settingIterator->second->ephemerisSettings,
                                                 settingIterator->first ) );
        }
        else
        {
            std::cerr<<"Warning, no ephemeris data found for body "<<
                       settingIterator->first<<std::endl;
        }
    }

    // Create atmosphere model objects for each body (if required).
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        if( settingIterator->second->atmosphereSettings != NULL )
        {
            bodyMap[ settingIterator->first ]->setAtmosphereModel(
                        createAtmosphereModel( settingIterator->second->atmosphereSettings,
                                               settingIterator->first ) );
        }
    }

    // Create rotation model objects for each body (if required).
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        if( settingIterator->second->rotationModelSettings != NULL )
        {
            bodyMap[ settingIterator->first ]->setRotationalEphemeris(
                        createRotationModel( settingIterator->second->rotationModelSettings,
                                             settingIterator->first ) );
        }
    }

    // Create gravity field model objects for each body (if required).
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        if( settingIterator->second->gravityFieldSettings != NULL )
        {
            bodyMap[ settingIterator->first ]->setGravityFieldModel(
                        createGravityFieldModel( settingIterator->second->gravityFieldSettings,
                                                 settingIterator->first ) );
        }
    }

    return bodyMap;

}

}

}
