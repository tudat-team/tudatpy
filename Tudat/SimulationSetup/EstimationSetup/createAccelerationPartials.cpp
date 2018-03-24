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
std::vector< boost::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > createTidalLoveNumberInterfaces(
        const NamedBodyMap& bodyMap,
        const std::string& acceleratingBodyName )
{
    // Create return map.
    std::vector< boost::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > loveNumberInterfaces;

    // Check if any gravity field variations are present
    if( bodyMap.at( acceleratingBodyName )->getGravityFieldVariationSet( ) != NULL )
    {
        // Get list of tidal gravity field variations.
        std::vector< boost::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > >  variationObjectList =
                utilities::dynamicCastSVectorToTVector< gravitation::GravityFieldVariations,
                gravitation::BasicSolidBodyTideGravityFieldVariations >(
                    bodyMap.at( acceleratingBodyName )->getGravityFieldVariationSet( )->
                    getDirectTidalGravityFieldVariations( ) );

        // Create partial for each gravity field variation objet
        if( variationObjectList.size( ) > 0 )
        {
            // Get state/rotation functions for deformed body
            boost::function< Eigen::Vector3d( ) > deformedBodyPositionFunction =
                    boost::bind( &Body::getPosition, bodyMap.at( acceleratingBodyName ) );
            boost::function< Eigen::Quaterniond( ) > rotationToDeformedBodyFrameFrameFunction =
                    boost::bind( &Body::getCurrentRotationToLocalFrame, bodyMap.at( acceleratingBodyName ) );

            for( unsigned int i = 0; i < variationObjectList.size( ); i++ )
            {
                if( variationObjectList.at( i ) != NULL )
                {
                    // Get state/rotation functions for deforming bodyies
                    std::vector< boost::function< Eigen::Vector3d( ) > > deformingBodyStateFunctions;
                    std::vector< std::string > deformingBodies = variationObjectList.at( i )->getDeformingBodies( );
                    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
                    {
                        deformingBodyStateFunctions.push_back(
                                    boost::bind( &Body::getPosition, bodyMap.at( deformingBodies.at( i ) ) ) );
                    }
                    // Get state/rotation functions for deformed body

                    // Create partial object
                    loveNumberInterfaces.push_back(
                                boost::make_shared< orbit_determination::TidalLoveNumberPartialInterface >(
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

}

}
