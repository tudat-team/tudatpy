/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

#if USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#endif

#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"

#include "Tudat/SimulationSetup/createGravityField.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a gravity field model.
boost::shared_ptr< gravitation::GravityFieldModel > createGravityFieldModel(
        const boost::shared_ptr< GravityFieldSettings > gravityFieldSettings,
        const std::string& body,
        const NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings )
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
    #if C_SPICE
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
    #endif
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

                if( gravityFieldVariationSettings.size( ) == 0 &&
                        sphericalHarmonicFieldSettings->getCreateTimeDependentField( ) == 0 )
                {
                    // Create and initialize spherical harmonic gravity field model.
                    gravityFieldModel = boost::make_shared< SphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ),
                                sphericalHarmonicFieldSettings->getAssociatedReferenceFrame( ) );
                }
                else
                {
                    if( bodyMap.at( body )->getGravityFieldModel( ) != NULL )
                    {
                        std::cerr<<"Warning when making time-dependent gravity field model for body "<<body<<" existing gravity field "
                                <<" is not empty but overwritten in Body! "<<std::endl;
                    }

                    // Create preliminary TimeDependentSphericalHarmonicsGravityField, without actual variation settings.
                    gravityFieldModel = boost::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ),
                                sphericalHarmonicFieldSettings->getAssociatedReferenceFrame( ) );
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

} // namespace simulation_setup

} // namespace tudat
