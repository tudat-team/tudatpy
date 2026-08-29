/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>
#include <functional>
#include <memory>

#include "tudat/simulation/propagation_setup/createGravityDeformationModels.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"

namespace tudat
{

namespace simulation_setup
{

using namespace aerodynamics;
using namespace gravitation;
using namespace basic_astrodynamics;
using namespace electromagnetism;
using namespace ephemerides;

std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > createMaxwellGravityFieldDeformationModel(
        const std::shared_ptr< simulation_setup::Body > deformingBody,
        const std::vector< std::shared_ptr< simulation_setup::Body > > perturbingBody,
        const std::string& nameOfDeformingBody,
        const std::vector< std::string >& nameOfPerturbingBody,
        const std::shared_ptr< GravityDeformationSettings > deformationSettings )
{
    // Declare pointer to return object
    std::shared_ptr< MaxwellGravityDeformationModel > deformationModel;

    // Dynamic cast deformation settings to required type and check consistency.
    std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings =
            std::dynamic_pointer_cast< MaxwellDeformationSettings >( deformationSettings );
    if( maxwellDeformationSettings == nullptr )
    {
        throw std::runtime_error( std::string( "Error, deformation settings inconsistent " ) +
                                  " when making maxwell gravity deformation of " + nameOfDeformingBody );
    }
    else
    {
        if( maxwellDeformationSettings->maximumDegree_ != 2 || maxwellDeformationSettings->maximumOrder_ != 2 )
        {
            throw std::runtime_error( "Error when creating Maxwell gravity deformation of " + nameOfDeformingBody +
                                      ": numerical gravity propagation currently supports degree and order 2 only." );
        }
        if( maxwellDeformationSettings->staticCoefficients_.size( ) != 5 )
        {
            throw std::runtime_error( "Error when creating Maxwell gravity deformation of " + nameOfDeformingBody +
                                      ": static degree-two coefficient vector must contain [C20, C21, C22, S21, S22]." );
        }

        // Get pointer to gravity field and rotational ephemeris of deforming body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( deformingBody->getGravityFieldModel( ) );

        std::shared_ptr< RotationalEphemeris > rotationalEphemeris = deformingBody->getRotationalEphemeris( );
        if( sphericalHarmonicsGravityField == nullptr )
        {
            throw std::runtime_error( std::string( "Error, spherical harmonic gravity field model not set when " ) +
                                      " creating Maxwell gravity deformation model of " + nameOfDeformingBody );
        }
        else
        {
            if( rotationalEphemeris == nullptr )
            {
                throw std::runtime_error( "Warning when creating Maxwell gravity deformation of " + nameOfDeformingBody +
                                          "no rotation model found for " + nameOfDeformingBody );
            }

            // Create gravity deformation object.
            std::vector< std::function< double( ) > > gravitationalParameterFunctionsPerturbingBodies;
            std::vector< basic_astrodynamics::MaxwellGravityDeformationModel::StateFunction > stateFunctionPerturbingBodies;
            for( unsigned int k = 0; k < perturbingBody.size( ); k++ )
            {
                const std::shared_ptr< gravitation::GravityFieldModel > perturbingGravityField =
                        perturbingBody.at( k )->getGravityFieldModel( );
                if( perturbingGravityField == nullptr )
                {
                    throw std::runtime_error( "Error when creating Maxwell gravity deformation of " + nameOfDeformingBody +
                                              ", no gravity field found for perturbing body " + nameOfPerturbingBody.at( k ) );
                }
                gravitationalParameterFunctionsPerturbingBodies.push_back(
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter, perturbingGravityField ) );
                stateFunctionPerturbingBodies.push_back(
                        std::bind( &Body::getStateByReference, perturbingBody.at( k ), std::placeholders::_1 ) );
            }

            std::function< double( ) > gravitationalParameterFunctionDeformingBody =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter, sphericalHarmonicsGravityField );

            deformationModel = std::make_shared< MaxwellGravityDeformationModel >(
                    std::bind( &Body::getStateByReference, deformingBody, std::placeholders::_1 ),
                    nameOfPerturbingBody,
                    maxwellDeformationSettings->maxwellRelaxationTime_,
                    maxwellDeformationSettings->globalRelaxationTime_,
                    gravitationalParameterFunctionDeformingBody,
                    gravitationalParameterFunctionsPerturbingBodies,
                    sphericalHarmonicsGravityField->getReferenceRadius( ),
                    std::bind( &Body::getCurrentAngularVelocityVectorInLocalFrame, deformingBody ),
                    std::bind( &Body::getCurrentAngularVelocityDerivativeVectorInLocalFrame, deformingBody ),
                    maxwellDeformationSettings->loveNumber_,
                    std::bind( &SphericalHarmonicsGravityField::getCosineCoefficientsBlock,
                               sphericalHarmonicsGravityField,
                               maxwellDeformationSettings->maximumDegree_,
                               maxwellDeformationSettings->maximumOrder_ ),
                    std::bind( &SphericalHarmonicsGravityField::getSineCoefficientsBlock,
                               sphericalHarmonicsGravityField,
                               maxwellDeformationSettings->maximumDegree_,
                               maxwellDeformationSettings->maximumOrder_ ),
                    stateFunctionPerturbingBodies,
                    std::bind( &Body::getCurrentRotationToGlobalFrame, deformingBody ),
                    std::bind( &Body::getCurrentRotationMatrixDerivativeToLocalFrame, deformingBody ),
                    maxwellDeformationSettings->staticCoefficients_,
                    maxwellDeformationSettings->includeOrder1_,
                    maxwellDeformationSettings->includeCentrifugalPotential_ );
        }
    }
    return deformationModel;
};

//! Function to create a list of gravity deformation models for a list of bodies.
basic_astrodynamics::GravityDeformationModelMap createGravityDeformationModelsMap(
        const SystemOfBodies& bodies,
        const SelectedGravityDeformationModelMap& gravityDeformationSettings )
{
    // Iterate over all bodies
    std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > > gravityDeformationModels;
    for( std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > >::const_iterator settingsIterator =
                 gravityDeformationSettings.begin( );
         settingsIterator != gravityDeformationSettings.end( );
         settingsIterator++ )
    {
        ensureIntegratedGravityFieldVariation( bodies.at( settingsIterator->first ), settingsIterator->first );

        // Iterate over all mass model settings for current body.
        for( unsigned int i = 0; i < settingsIterator->second.size( ); i++ )
        {
            switch( settingsIterator->second.at( i )->deformationType_ )
            {
                case maxwell_deformation: {
                    std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings =
                            std::dynamic_pointer_cast< MaxwellDeformationSettings >( settingsIterator->second.at( i ) );
                    std::string deformingBody = settingsIterator->first;

                    std::vector< std::string > perturbingBodyNames = maxwellDeformationSettings->perturbingBody_;
                    std::vector< std::shared_ptr< simulation_setup::Body > > perturbingBodies;
                    for( unsigned int k = 0; k < perturbingBodyNames.size( ); k++ )
                    {
                        perturbingBodies.push_back( bodies.at( perturbingBodyNames.at( k ) ) );
                    }

                    gravityDeformationModels[ settingsIterator->first ].push_back(
                            createMaxwellGravityFieldDeformationModel( bodies.at( deformingBody ),
                                                                       perturbingBodies,
                                                                       deformingBody,
                                                                       perturbingBodyNames,
                                                                       settingsIterator->second.at( i ) ) );
                    break;
                }
                default:
                    break;
            }
        }
    }
    return gravityDeformationModels;
}

}  // namespace simulation_setup

}  // namespace tudat
