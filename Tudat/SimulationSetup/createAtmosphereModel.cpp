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
#include <boost/lexical_cast.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"
#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"
#include "Tudat/SimulationSetup/createAtmosphereModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create an atmosphere model.
boost::shared_ptr< aerodynamics::AtmosphereModel > createAtmosphereModel(
        const boost::shared_ptr< AtmosphereSettings > atmosphereSettings,
        const std::string& body )
{
    using namespace tudat::aerodynamics;

    // Declare return object.
    boost::shared_ptr< AtmosphereModel > atmosphereModel;

    // Check which type of atmosphere model is to be created.
    switch( atmosphereSettings->getAtmosphereType( ) )
    {
    case exponential_atmosphere:
    {
        // Check whether settings for atmosphere are consistent with its type.
        boost::shared_ptr< ExponentialAtmosphereSettings > exponentialAtmosphereSettings =
                boost::dynamic_pointer_cast< ExponentialAtmosphereSettings >( atmosphereSettings );
        if( exponentialAtmosphereSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected exponential atmosphere settings for body " + body );
        }
        else
        {
            // Create and initialize exponential atmosphere model.
            boost::shared_ptr< ExponentialAtmosphere > exponentialAtmosphereModel =
                    boost::make_shared< ExponentialAtmosphere >(
                        exponentialAtmosphereSettings->getDensityScaleHeight( ) ,
                        exponentialAtmosphereSettings->getConstantTemperature( ),
                        exponentialAtmosphereSettings->getDensityAtZeroAltitude( ),
                        exponentialAtmosphereSettings->getSpecificGasConstant( ) );
            atmosphereModel = exponentialAtmosphereModel;
        }
        break;
    }
    case tabulated_atmosphere:
    {
        // Check whether settings for atmosphere are consistent with its type
        boost::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
                boost::dynamic_pointer_cast< TabulatedAtmosphereSettings >( atmosphereSettings );
        if( tabulatedAtmosphereSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected tabulated atmosphere settings for body " + body );
        }
        else
        {
            // Create and initialize tabulatedl atmosphere model.
            atmosphereModel = boost::make_shared< TabulatedAtmosphere >(
                        tabulatedAtmosphereSettings->getAtmosphereFile( ) );
        }
        break;
    }
    default:
        throw std::runtime_error(
                 "Error, did not recognize atmosphere model settings type " +
                  boost::lexical_cast< std::string >( atmosphereSettings->getAtmosphereType( ) ) );
    }
    return atmosphereModel;
}


}

}
