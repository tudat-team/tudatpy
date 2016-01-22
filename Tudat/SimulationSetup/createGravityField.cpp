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

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/SimulationSetup/createGravityField.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a gravity field model.
boost::shared_ptr< gravitation::GravityFieldModel > createGravityFieldModel(
        const boost::shared_ptr< GravityFieldSettings > gravityFieldSettings,
        const std::string& body )
{
    using namespace tudat::gravitation;

    // Declare return object.
    boost::shared_ptr< GravityFieldModel > gravityFieldModel;

    // Check which type of gravity field model is to be created.
    switch( gravityFieldSettings->getGravityFieldType( ) )
    {
    case central:
    {
        // Check whether settings for point mass gravity field model are consistent with its type.
        boost::shared_ptr< CentralGravityFieldSettings > centralFieldSettings =
                boost::dynamic_pointer_cast< CentralGravityFieldSettings >( gravityFieldSettings );
        if( centralFieldSettings == NULL )
        {
            throw std::runtime_error(
                "Error, expected central field settings when making gravity field model for body " +
                 body);
        }
        else
        {
            // Create and initialize point mass gravity field model.
            gravityFieldModel = boost::make_shared< GravityFieldModel >(
                        centralFieldSettings->getGravitationalParameter( ) );
        }
        break;
    }
    case central_spice:
    {
        // Create and initialize point mass gravity field model from Spice.
        gravityFieldModel = boost::make_shared< GravityFieldModel >(
                    spice_interface::getBodyGravitationalParameter( body ) );

        break;
    }
    case spherical_harmonic:
    {
        // Check whether settings for spherical harmonic gravity field model are consistent with
        // its type.
        boost::shared_ptr< SphericalHarmonicsGravityFieldSettings > sphericalHarmonicFieldSettings =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                    gravityFieldSettings );

        if( sphericalHarmonicFieldSettings == NULL )
        {
            throw std::runtime_error(
             "Error, expected spherical harmonic field settings when making gravity field model of "
                        + body );
        }
        else
        {
            // Check consistency of cosine and sine coefficients.
            if( ( sphericalHarmonicFieldSettings->getCosineCoefficients( ).rows( ) !=
                  sphericalHarmonicFieldSettings->getSineCoefficients( ).rows( ) ) ||
                    ( sphericalHarmonicFieldSettings->getCosineCoefficients( ).cols( ) !=
                      sphericalHarmonicFieldSettings->getSineCoefficients( ).cols( ) ) )
            {
                throw std::runtime_error(
                            std::string( "Error when making spherical harmonic field, sine and " ) +
                            std::string( "cosine matrix  sizes are not equal for body " ) + body );
            }
            else
            {
                // Create and initialize spherical harmonic gravity field model.
                gravityFieldModel = boost::make_shared< SphericalHarmonicsGravityField >(
                            sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                            sphericalHarmonicFieldSettings->getReferenceRadius( ),
                            sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                            sphericalHarmonicFieldSettings->getSineCoefficients( ) );
            }


        }
        break;
    }
    default:
        throw std::runtime_error(
                 "Error, did not recognize gravity field model settings type " +
                  boost::lexical_cast< std::string >(
                        gravityFieldSettings->getGravityFieldType( ) ) );
    }

    return gravityFieldModel;
}

}

}
