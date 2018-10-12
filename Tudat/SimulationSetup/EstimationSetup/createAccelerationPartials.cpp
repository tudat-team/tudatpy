/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EstimationSetup/createAccelerationPartials.h"
#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a list of objects that can be used to compute partials of tidal gravity field variations
std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > createTidalLoveNumberInterfaces(
        const NamedBodyMap& bodyMap,
        const std::string& acceleratingBodyName )
{
    // Create return map.
    std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > loveNumberInterfaces;

    // Check if any gravity field variations are present
    if( bodyMap.at( acceleratingBodyName )->getGravityFieldVariationSet( ) != nullptr )
    {
        // Get list of tidal gravity field variations.
        std::vector< std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > >  variationObjectList =
                utilities::dynamicCastSVectorToTVector< gravitation::GravityFieldVariations,
                gravitation::BasicSolidBodyTideGravityFieldVariations >(
                    bodyMap.at( acceleratingBodyName )->getGravityFieldVariationSet( )->
                    getDirectTidalGravityFieldVariations( ) );

        // Create partial for each gravity field variation objet
        if( variationObjectList.size( ) > 0 )
        {
            // Get state/rotation functions for deformed body
            std::function< Eigen::Vector3d( ) > deformedBodyPositionFunction =
                    std::bind( &Body::getPosition, bodyMap.at( acceleratingBodyName ) );
            std::function< Eigen::Quaterniond( ) > rotationToDeformedBodyFrameFrameFunction =
                    std::bind( &Body::getCurrentRotationToLocalFrame, bodyMap.at( acceleratingBodyName ) );

            for( unsigned int i = 0; i < variationObjectList.size( ); i++ )
            {
                if( variationObjectList.at( i ) != nullptr )
                {
                    // Get state/rotation functions for deforming bodyies
                    std::vector< std::function< Eigen::Vector3d( ) > > deformingBodyStateFunctions;
                    std::vector< std::string > deformingBodies = variationObjectList.at( i )->getDeformingBodies( );
                    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
                    {
                        deformingBodyStateFunctions.push_back(
                                    std::bind( &Body::getPosition, bodyMap.at( deformingBodies.at( i ) ) ) );
                    }
                    // Get state/rotation functions for deformed body

                    // Create partial object
                    loveNumberInterfaces.push_back(
                                std::make_shared< orbit_determination::TidalLoveNumberPartialInterface >(
                                    variationObjectList.at( i ),
                                    deformedBodyPositionFunction,
                                    deformingBodyStateFunctions,
                                    rotationToDeformedBodyFrameFrameFunction,
                                    acceleratingBodyName ) );
                }
            }
        }
    }
    return loveNumberInterfaces;
}

template std::shared_ptr< acceleration_partials::AccelerationPartial > createAnalyticalAccelerationPartial< double >(
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratingBody,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate );

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template std::shared_ptr< acceleration_partials::AccelerationPartial > createAnalyticalAccelerationPartial< long double >(
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratingBody,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
        parametersToEstimate );
#endif

template orbit_determination::StateDerivativePartialsMap createAccelerationPartialsMap< double >(
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate );
        template orbit_determination::StateDerivativePartialsMap createAccelerationPartialsMap< long double >(
                const basic_astrodynamics::AccelerationMap& accelerationMap,
                const simulation_setup::NamedBodyMap& bodyMap,
                const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
                parametersToEstimate );

}

}
