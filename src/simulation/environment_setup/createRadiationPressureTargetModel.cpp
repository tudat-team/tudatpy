/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"

#include <map>
#include <memory>
#include <vector>

#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/body.h"


namespace tudat
{
namespace simulation_setup
{


std::shared_ptr<electromagnetism::RadiationPressureTargetModel> getRadiationPressureTargetModelOfType(
    const std::shared_ptr< Body > target,
    const RadiationPressureTargetModelType targetModelType,
    const std::string errorOutput )
{
    std::string targetName = target->getBodyName( );
    std::shared_ptr<electromagnetism::RadiationPressureTargetModel> targetModel;
    if( target->getRadiationPressureTargetModels( ).size( ) == 1 )
    {
        std::shared_ptr<electromagnetism::RadiationPressureTargetModel> availableTargetModel = target->getRadiationPressureTargetModels( ).at( 0 );
        if( getTargetModelType( availableTargetModel ) == targetModelType || targetModelType == undefined_target )
        {
            targetModel = availableTargetModel;
        }
        else
        {
            throw std::runtime_error( "Error " + errorOutput + ", body " +
                                      targetName +
                                      " has no radiation pressure target model of type ." + std::to_string( targetModelType ) );
        }
    }
    else if( targetModelType == undefined_target )
    {
        throw std::runtime_error( "Error " + errorOutput + ", body " +
                                  targetName +
                                  " has multiple radiation pressure target models" );
    }
    else
    {
        std::vector< std::shared_ptr<electromagnetism::RadiationPressureTargetModel> > availableTargetModels = target->getRadiationPressureTargetModels( );
        for( unsigned int i = 0; i < availableTargetModels.size( ); i++ )
        {
            if( getTargetModelType( availableTargetModels.at( i ) ) == targetModelType )
            {
                if( targetModel == nullptr )
                {
                    targetModel = availableTargetModels.at( i );
                }
                else
                {
                    throw std::runtime_error( "Error " + errorOutput + ", body " +
                                              targetName +
                                              " has multiple martching radiation pressure target models" );
                }
            }
        }

        if( targetModel == nullptr )
        {
            throw std::runtime_error( "Error " + errorOutput + ", body " +
                                      targetName +
                                      " has multiple radiation pressure target models, but none match required type" );
        }
    }
    return targetModel;
}

std::shared_ptr< electromagnetism::ReflectionLaw > createReflectionLaw(
    const std::shared_ptr< BodyPanelReflectionLawSettings > reflectionLawSettings )
{
    std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw;
    switch ( reflectionLawSettings->bodyPanelReflectionLawType_ )
    {
    case specular_diffuse_reflection_law:
    {
        std::shared_ptr< SpecularDiffuseBodyPanelReflectionLawSettings > specularDiffuseReflectionLawSettings =
            std::dynamic_pointer_cast< SpecularDiffuseBodyPanelReflectionLawSettings >( reflectionLawSettings );
        if( specularDiffuseReflectionLawSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating specular-diffuse reflection law; settings are incompatible" );
        }
        reflectionLaw = std::make_shared< electromagnetism::SpecularDiffuseMixReflectionLaw >(
            specularDiffuseReflectionLawSettings->absorptivity_,
            specularDiffuseReflectionLawSettings->specularReflectivity_,
            specularDiffuseReflectionLawSettings->diffuseReflectivity_,
            specularDiffuseReflectionLawSettings->withInstantaneousReradiation_ );
        break;
    }
    default:
        throw std::runtime_error( "Error when creating panel reflection law; type not recognzied " +
                                  std::to_string( static_cast< int>( reflectionLawSettings->bodyPanelReflectionLawType_ ) ) );
    }
    return reflectionLaw;
}

RadiationPressureTargetModelType getTargetModelType( const std::shared_ptr<electromagnetism::RadiationPressureTargetModel> targetModel )
{
    RadiationPressureTargetModelType targetModelType;
    if( std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >( targetModel ) != nullptr )
    {
        targetModelType = cannonball_target;
    }
    else if( std::dynamic_pointer_cast< electromagnetism::PaneledRadiationPressureTargetModel >( targetModel ) != nullptr )
    {
        targetModelType = paneled_target;
    }
    else
    {
        throw std::runtime_error( "Error when finding radiation pressure target model type, type not recognized." );
    }
    return targetModelType;
}


std::vector< std::shared_ptr<electromagnetism::RadiationPressureTargetModel> > createRadiationPressureTargetModel(
        const std::shared_ptr<RadiationPressureTargetModelSettings>& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies)
{
    using namespace tudat::electromagnetism;

    std::vector< std::shared_ptr<electromagnetism::RadiationPressureTargetModel> > radiationPressureTargetModels;

    // Validate occulting bodies map: can either have
    //  - multiple entries with source body names as keys
    //  - a single entry with the empty string as key (use same occulting bodies for all sources)
    auto sourceToTargetOccultingBodies = modelSettings->getSourceToTargetOccultingBodies();
    if (sourceToTargetOccultingBodies.count("") > 0 && sourceToTargetOccultingBodies.size() > 1)
    {
        throw std::runtime_error("Error, invalid occulting bodies map for " + body );
    }

    switch(modelSettings->getRadiationPressureTargetModelType())
    {
        case RadiationPressureTargetModelType::cannonball_target:
        {
            auto cannonballTargetModelSettings =
                    std::dynamic_pointer_cast< CannonballRadiationPressureTargetModelSettings >(modelSettings);

            if(cannonballTargetModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected cannonball radiation pressure target for body " + body );
            }

            radiationPressureTargetModels.push_back( std::make_shared<CannonballRadiationPressureTargetModel>(
                    cannonballTargetModelSettings->getArea(),
                    cannonballTargetModelSettings->getCoefficient(),
                    sourceToTargetOccultingBodies) );
            break;
        }
        case RadiationPressureTargetModelType::paneled_target:
        {
            if( bodies.at( body )->getVehicleSystems( ) == nullptr )
            {
                throw std::runtime_error( "Error, requested panelled radiation pressure model for " + body + ", but no system models found" );
            }

            std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > > sortedBodyPanelMap =
                bodies.at( body )->getVehicleSystems()->getVehicleExteriorPanels( );
            if( sortedBodyPanelMap.size( ) == 0 )
            {
                throw std::runtime_error( "Error, requested panelled radiation pressure model for " + body + ", no panels defined" );
            }

            std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > bodyFixedPanels;
            std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > > segmentFixedPanels;
            std::map< std::string, std::function< Eigen::Quaterniond( ) > > segmentFixedToBodyFixedRotations;

            for( auto it : sortedBodyPanelMap )
            {
                if( it.first != "" )
                {
                    segmentFixedPanels[ it.first ] = it.second;
                    segmentFixedToBodyFixedRotations[ it.first ] = std::bind(
                        &system_models::VehicleSystems::getPartRotationToBaseFrame, bodies.at( body )->getVehicleSystems(), it.first );
                }
                else
                {
                    bodyFixedPanels = it.second;
                }
            }


            radiationPressureTargetModels.push_back( std::make_shared<PaneledRadiationPressureTargetModel>(
                bodyFixedPanels, segmentFixedPanels, segmentFixedToBodyFixedRotations, sourceToTargetOccultingBodies) );
            break;
        }
        case RadiationPressureTargetModelType::multi_type_target:
        {
            std::shared_ptr< MultiRadiationPressureTargetModelSettings > multiTargetModelSettings =
                std::dynamic_pointer_cast< MultiRadiationPressureTargetModelSettings >(modelSettings);

            if(multiTargetModelSettings == nullptr)
            {
                throw std::runtime_error(
                    "Error, expected multi-type radiation pressure target for body " + body );
            }
            for( unsigned int i = 0; i < multiTargetModelSettings->radiationPressureTargetModelSettings_.size( ); i++ )
            {
                if( multiTargetModelSettings->radiationPressureTargetModelSettings_.at( i )->getRadiationPressureTargetModelType( ) == multi_type_target )
                {
                    throw std::runtime_error( "Error, cannot have multiple nested multi-type target model settings" );
                }
                radiationPressureTargetModels.push_back( createRadiationPressureTargetModel(
                    multiTargetModelSettings->radiationPressureTargetModelSettings_.at( i ), body, bodies ).at( 0 ) );
            }
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize radiation pressure target model settings for " + body );
    }

    return radiationPressureTargetModels;
}

} // tudat
} // electromagnetism
