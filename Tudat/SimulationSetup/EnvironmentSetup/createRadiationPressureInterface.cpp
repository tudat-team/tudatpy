/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/bind.hpp>

#include "Tudat/SimulationSetup/EnvironmentSetup/createRadiationPressureInterface.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to obtain (by reference) the position functions and radii of occulting bodies
void getOccultingBodiesInformation(
        const NamedBodyMap& bodyMap, const std::vector< std::string >& occultingBodies,
        std::vector< boost::function< Eigen::Vector3d( ) > >& occultingBodyPositions,
        std::vector< double >& occultingBodyRadii )
{
    // Iterate over occulring bodies and retrieve radius and position function.
    for( unsigned int i = 0; i < occultingBodies.size( ); i++ )
    {
        if( bodyMap.count( occultingBodies[ i ] ) == 0 )
        {
           throw std::runtime_error( "Error, could not find body " + occultingBodies[ i ] +
                                     " in body map when making occulting body settings" );
        }
        else
        {
            boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel =
                    bodyMap.at( occultingBodies.at( i ) )->getShapeModel( );
            if( shapeModel == NULL )
            {
                throw std::runtime_error( "Error, no shape model for " + occultingBodies[ i ] +
                                          " when making occulting body settings" );
            }
            occultingBodyPositions.push_back(
                        boost::bind( &Body::getPosition, bodyMap.at( occultingBodies[ i ] ) ) );
            occultingBodyRadii.push_back( shapeModel->getAverageRadius( ) );
        }
    }
}

//! Function to create a radiation pressure interface.
boost::shared_ptr< electro_magnetism::RadiationPressureInterface > createRadiationPressureInterface(
        const boost::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureInterfaceSettings,
        const std::string& bodyName, const NamedBodyMap& bodyMap )
{
    boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface;

    // Check type of radiation pressure interface
    switch( radiationPressureInterfaceSettings->getRadiationPressureType( ) )
    {
    case cannon_ball:
    {
        // Check type consistency.
        boost::shared_ptr< CannonBallRadiationPressureInterfaceSettings > cannonBallSettings =
                boost::dynamic_pointer_cast< CannonBallRadiationPressureInterfaceSettings >(
                    radiationPressureInterfaceSettings );
        if( cannonBallSettings == NULL )
        {
            throw std::runtime_error( "Error when making cannon ball radiation interface, type does not match object" );
        }

        // Retrieve source body and check consistency.
        if( bodyMap.count( radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw std::runtime_error( "Error when making cannon ball radiation interface, source not found.");
        }

        boost::shared_ptr< Body > sourceBody =
                bodyMap.at( radiationPressureInterfaceSettings->getSourceBody( ) );

        // Get reqruied data for occulting bodies.
        std::vector< std::string > occultingBodies = cannonBallSettings->getOccultingBodies( );
        std::vector< boost::function< Eigen::Vector3d( ) > > occultingBodyPositions;
        std::vector< double > occultingBodyRadii;
        getOccultingBodiesInformation(
                    bodyMap, occultingBodies, occultingBodyPositions, occultingBodyRadii );

        // Retrive radius of source if occultations are used.
        double sourceRadius;
        if( occultingBodyPositions.size( ) > 0 )
        {
            boost::shared_ptr< basic_astrodynamics::BodyShapeModel > sourceShapeModel =
                    sourceBody->getShapeModel( );

            if( sourceShapeModel == NULL )
            {
                throw std::runtime_error( "Error when making occulted body, source body " +
                                          radiationPressureInterfaceSettings->getSourceBody( ) +
                                          " does not have a shape" );
            }
            else
            {
                sourceRadius = sourceShapeModel->getAverageRadius( );
            }
        }
        else
        {
            sourceRadius = 0.0;
        }

        // Create function returning radiated power.
        boost::function< double( ) > radiatedPowerFunction;
        if( defaultRadiatedPowerValues.count(
                    radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw std::runtime_error( "Error, no radiated power found for " +
                                      radiationPressureInterfaceSettings->getSourceBody( ) );
        }
        else
        {
            radiatedPowerFunction = boost::lambda::constant(
                        defaultRadiatedPowerValues.at(
                            radiationPressureInterfaceSettings->getSourceBody( ) ) );
        }

        // Create radiation pressure interface.
        radiationPressureInterface =
                boost::make_shared< electro_magnetism::RadiationPressureInterface >(
                    radiatedPowerFunction,
                    boost::bind( &Body::getPosition, sourceBody ),
                    boost::bind( &Body::getPosition, bodyMap.at( bodyName ) ),
                    cannonBallSettings->getRadiationPressureCoefficient( ),
                    cannonBallSettings->getArea( ), occultingBodyPositions, occultingBodyRadii,
                    sourceRadius );
        break;

    }
    default:
        throw std::runtime_error(
                    "Error, radiation pressure type" + std::to_string(
                        radiationPressureInterfaceSettings->getRadiationPressureType( ) ) +
                    "not recognized for body" + bodyName );
    }

    return radiationPressureInterface;
}

} // namespace simulation_setup

} // namespace tudat
