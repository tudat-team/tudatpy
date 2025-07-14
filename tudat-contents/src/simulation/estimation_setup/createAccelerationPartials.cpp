/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a list of objects that can be used to compute partials of tidal gravity field variations
std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > createTidalLoveNumberInterfaces(
        const SystemOfBodies& bodies,
        const std::string& acceleratingBodyName )
{
    // Create return map.
    std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > loveNumberInterfaces;

    // Check if any gravity field variations are present
    if( bodies.at( acceleratingBodyName )->getGravityFieldVariationSet( ) != nullptr )
    {
        // Get list of tidal gravity field variations.
        std::vector< std::shared_ptr< gravitation::SolidBodyTideGravityFieldVariations > > variationObjectList =
                utilities::dynamicCastSVectorToTVector< gravitation::GravityFieldVariations,
                                                        gravitation::SolidBodyTideGravityFieldVariations >(
                        bodies.at( acceleratingBodyName )->getGravityFieldVariationSet( )->getDirectTidalGravityFieldVariations( ) );

        // Create partial for each gravity field variation objet
        if( variationObjectList.size( ) > 0 )
        {
            // Get state/rotation functions for deformed body
            std::function< Eigen::Vector3d( ) > deformedBodyPositionFunction =
                    std::bind( &Body::getPosition, bodies.at( acceleratingBodyName ) );
            std::function< Eigen::Quaterniond( ) > rotationToDeformedBodyFrameFrameFunction =
                    std::bind( &Body::getCurrentRotationToLocalFrame, bodies.at( acceleratingBodyName ) );

            for( unsigned int i = 0; i < variationObjectList.size( ); i++ )
            {
                if( variationObjectList.at( i ) != nullptr )
                {
                    // Get state/rotation functions for deforming bodyies
                    std::vector< std::function< Eigen::Vector3d( ) > > deformingBodyStateFunctions;
                    std::vector< std::string > deformingBodies = variationObjectList.at( i )->getDeformingBodies( );
                    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
                    {
                        deformingBodyStateFunctions.push_back( std::bind( &Body::getPosition, bodies.at( deformingBodies.at( i ) ) ) );
                    }
                    // Create partial object
                    auto newLoveNumberInterface = std::make_shared< orbit_determination::TidalLoveNumberPartialInterface >(
                            variationObjectList.at( i ),
                            deformedBodyPositionFunction,
                            deformingBodyStateFunctions,
                            rotationToDeformedBodyFrameFrameFunction,
                            acceleratingBodyName );

                    bool addInterface = true;
                    for( unsigned int j = 0; j < loveNumberInterfaces.size( ); j++ )
                    {
                        if( acceleratingBodyName == loveNumberInterfaces.at( j )->getDeformedBody( ) )
                        {
                            if( utilities::compareStlVectors< std::string >( newLoveNumberInterface->getDeformingBodies( ),
                                                                             loveNumberInterfaces.at( j )->getDeformingBodies( ) ) )
                            {
                                addInterface = false;
                            }
                        }
                    }

                    if( addInterface )
                    {
                        loveNumberInterfaces.push_back( newLoveNumberInterface );
                    }
                }
            }
        }
    }
    return loveNumberInterfaces;
}

std::map< std::string, std::shared_ptr< acceleration_partials::AccelerationPartial > > createEihAccelerationPartials(
        const std::map< std::string, std::shared_ptr< relativity::EinsteinInfeldHoffmannAcceleration > > eihAccelerations )
{
    std::map< std::string, std::shared_ptr< acceleration_partials::AccelerationPartial > > partialsList;
    if( eihAccelerations.size( ) > 0 )
    {
        std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations = eihAccelerations.begin( )->second->getEihEquations( );

        for( auto it: eihAccelerations )
        {
            if( it.second->getEihEquations( ) != eihEquations )
            {
                throw std::runtime_error( "Error when making acceleration partials, different governing EIH equations found" );
            }
        }

        std::shared_ptr< acceleration_partials::EihEquationsPartials > eihPartials =
                std::make_shared< acceleration_partials::EihEquationsPartials >( eihEquations );

        for( auto it: eihAccelerations )
        {
            partialsList[ it.first ] = std::make_shared< acceleration_partials::EihAccelerationPartial >( eihPartials, it.first );
        }
    }
    return partialsList;
}

std::shared_ptr< estimatable_parameters::NumericalAccelerationPartialSettings > getDefaultPanelledSurfaceRadiationPressurePartialSettings(
        const std::string bodyUndergoingAcceleration,
        const std::string bodyExertingAcceleration )
{
    return std::make_shared< estimatable_parameters::NumericalAccelerationPartialSettings >(
            ( Eigen::VectorXd( 6 ) << 100.0, 100.0, 100.0, 1.0E-3, 1.0E-3, 1.0E-3 ).finished( ),
            bodyUndergoingAcceleration,
            bodyExertingAcceleration,
            basic_astrodynamics::radiation_pressure );
}

}  // namespace simulation_setup

}  // namespace tudat
