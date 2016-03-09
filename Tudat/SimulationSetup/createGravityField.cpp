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

#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/SimulationSetup/createGravityField.h"

namespace tudat
{

namespace simulation_setup
{


boost::shared_ptr< gravitation::GravityFieldVariationsSet > createGravityFieldModelVariationsSet(
        const std::string& body,
        const NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings )
{

    using namespace tudat::gravitation;

    // Declare lists for input to GravityFieldVariationsSet
    std::vector< boost::shared_ptr< GravityFieldVariations > > variationObjects;
    std::vector< BodyDeformationTypes > variationTypes;
    std::vector< std::string > variationIdentifiers;
    std::map< int, boost::shared_ptr< interpolators::InterpolatorSettings > > createInterpolators;
    std::map< int, double > initialTimes;
    std::map< int, double > finalTimes;
    std::map< int, double > timeSteps;

    // Iterate over all variations to create.
    for( unsigned int i = 0; i < gravityFieldVariationSettings.size( ); i++ )
    {
        // Get current type of deformation
        variationTypes.push_back( gravityFieldVariationSettings[ i ]->getBodyDeformationType( ) );

        // Set current variation object in list.
        variationObjects.push_back( createGravityFieldVariationsModel(
                                        gravityFieldVariationSettings[ i ], body, bodyMap  ) );

        variationIdentifiers.push_back( "" );

        // Check if current variation is interpolated, and set settings if necessary.
        if( gravityFieldVariationSettings[ i ]->getInterpolateVariation( ) == true )
        {
            createInterpolators[ i ] = gravityFieldVariationSettings[ i ]->getInterpolatorSettings( );
            initialTimes[ i ] = gravityFieldVariationSettings[ i ]->getInitialTime( );
            finalTimes[ i ] = gravityFieldVariationSettings[ i ]->getFinalTime( );
            timeSteps[ i ] = gravityFieldVariationSettings[ i ]->getTimeStep( );
        }
    }

    // Create object with settings for updating variations from new parameter values.
    boost::shared_ptr< GravityFieldVariationsSet > fieldVariationsSet =
            boost::make_shared< GravityFieldVariationsSet >(
                variationObjects, variationTypes, variationIdentifiers,
                createInterpolators, initialTimes, finalTimes, timeSteps );

    if( boost::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >(
                bodyMap.at( body )->getGravityFieldModel( ) ) == NULL )
    {
        std::cerr<<"Error when making gravity field variations of body "<<body<<", base type is not time dependent"<<std::endl;
    }

    boost::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >(
                bodyMap.at( body )->getGravityFieldModel( ) )->setFieldVariationSettings( fieldVariationsSet, 1 );


    return fieldVariationsSet;
}

//! Function to create a gravity field model.
boost::shared_ptr< gravitation::GravityFieldModel > createGravityFieldModel(
        const boost::shared_ptr< GravityFieldSettings > gravityFieldSettings,
        const std::string& body,
        const NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings,
        const bool immediatelySetVariations )
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
        else if( gravityFieldVariationSettings.size( ) != 0 )
        {
            std::cerr<<"Error, requested central gravity field, but field variations settings are not empty."<<std::endl;
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
        if( gravityFieldVariationSettings.size( ) != 0 )
        {
            std::cerr<<"Error, requested central gravity field, but field variations settings are not empty."<<std::endl;
        }
        else
        {
            // Create and initialize point mass gravity field model from Spice.
            gravityFieldModel = boost::make_shared< GravityFieldModel >(
                        spice_interface::getBodyGravitationalParameter( body ) );
        }

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
            std::cerr<<"Warning, spherical harmonic reference frame rotation not found"<<std::endl;
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

                if( gravityFieldVariationSettings.size( ) == 0 &&
                        sphericalHarmonicFieldSettings->getCreateTimeDependentField( ) == 0 )
                {
                    // Create and initialize spherical harmonic gravity field model.
                    gravityFieldModel = boost::make_shared< SphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ) );
                }
                else
                {
                    if( bodyMap.at( body )->getGravityFieldModel( ) != NULL )
                    {
                        std::cerr<<"Warning when making time-dependent gravity field model for body "<<body<<" existing gravity field "
                                <<" is not empty but overwritten in CelestialBody! "<<std::endl;
                    }

                    // Create preliminary TimeDependentSphericalHarmonicsGravityField, without actual variation settings.
                    gravityFieldModel = boost::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ) );
                    if( immediatelySetVariations && gravityFieldVariationSettings.size( ) > 0 )
                    {
                        bodyMap.at( body )->setGravityFieldModel( gravityFieldModel );
                        boost::shared_ptr< gravitation::GravityFieldVariationsSet > gravityFieldVariations =
                                createGravityFieldModelVariationsSet( body, bodyMap, gravityFieldVariationSettings );
                        bodyMap.at( body )->setGravityFieldVariationSet( gravityFieldVariations );
                    }
                }
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
