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
        std::vector< std::function< Eigen::Vector3d( ) > >& occultingBodyPositions,
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
            std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel =
                    bodyMap.at( occultingBodies.at( i ) )->getShapeModel( );
            if( shapeModel == nullptr )
            {
                throw std::runtime_error( "Error, no shape model for " + occultingBodies[ i ] +
                                          " when making occulting body settings" );
            }
            occultingBodyPositions.push_back(
                        std::bind( &Body::getPosition, bodyMap.at( occultingBodies[ i ] ) ) );
            occultingBodyRadii.push_back( shapeModel->getAverageRadius( ) );
        }
    }
}


//! Function to obtain (by reference) the position and velocity functions of the central body.
void getCentralBodyInformation(
    const NamedBodyMap& bodyMap, const std::string& centralBody,
    std::function< Eigen::Vector3d( ) >& centralBodyPosition,
    std::function< Eigen::Vector3d( ) >& centralBodyVelocity )
{

    // Check that the central body is defined.
    if( bodyMap.count( centralBody ) == 0 )
    {
        throw std::runtime_error( "Error, could not find body " + centralBody +
                                 " in body map when making central body body settings" );
    }
    // Retrieve position and velocity functions.
    else
    {
        centralBodyPosition = std::bind( &Body::getPosition, bodyMap.at( centralBody ) );
        centralBodyVelocity = std::bind( &Body::getVelocity, bodyMap.at( centralBody ) );
    }
}


//! Function to create a radiation pressure interface.
std::shared_ptr< electro_magnetism::RadiationPressureInterface > createRadiationPressureInterface(
        const std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureInterfaceSettings,
        const std::string& bodyName, const NamedBodyMap& bodyMap )
{
    std::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface;

    // Check type of radiation pressure interface
    switch( radiationPressureInterfaceSettings->getRadiationPressureType( ) )
    {
    case cannon_ball_radiation_pressure_interface:
    {
        // Check type consistency.
        std::shared_ptr< CannonBallRadiationPressureInterfaceSettings > cannonBallSettings =
                std::dynamic_pointer_cast< CannonBallRadiationPressureInterfaceSettings >(
                    radiationPressureInterfaceSettings );
        if( cannonBallSettings == nullptr )
        {
            throw std::runtime_error( "Error when making cannon ball radiation interface, type does not match object" );
        }

        // Retrieve source body and check consistency.
        if( bodyMap.count( radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw std::runtime_error( "Error when making cannon ball radiation interface, source not found.");
        }

        std::shared_ptr< Body > sourceBody =
                bodyMap.at( radiationPressureInterfaceSettings->getSourceBody( ) );

        // Get reqruied data for occulting bodies.
        std::vector< std::string > occultingBodies = cannonBallSettings->getOccultingBodies( );
        std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions;
        std::vector< double > occultingBodyRadii;
        getOccultingBodiesInformation(
                    bodyMap, occultingBodies, occultingBodyPositions, occultingBodyRadii );

        // Retrive radius of source if occultations are used.
        double sourceRadius;
        if( occultingBodyPositions.size( ) > 0 )
        {
            std::shared_ptr< basic_astrodynamics::BodyShapeModel > sourceShapeModel =
                    sourceBody->getShapeModel( );

            if( sourceShapeModel == nullptr )
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
        std::function< double( ) > radiatedPowerFunction;
        if( defaultRadiatedPowerValues.count(
                    radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw std::runtime_error( "Error, no radiated power found for " +
                                      radiationPressureInterfaceSettings->getSourceBody( ) );
        }
        else
        {
            radiatedPowerFunction = [ = ]( ){ return
                        defaultRadiatedPowerValues.at(
                            radiationPressureInterfaceSettings->getSourceBody( ) ); };
        }

        // Create radiation pressure interface.
        radiationPressureInterface =
                std::make_shared< electro_magnetism::RadiationPressureInterface >(
                    radiatedPowerFunction,
                    std::bind( &Body::getPosition, sourceBody ),
                    std::bind( &Body::getPosition, bodyMap.at( bodyName ) ),
                    cannonBallSettings->getRadiationPressureCoefficient( ),
                    cannonBallSettings->getArea( ), occultingBodyPositions, occultingBodyRadii,
                    sourceRadius );
        break;
    }
    case panelled_radiation_pressure_interface:
    {
        // Check type consistency.
        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > panelledSettings =
                std::dynamic_pointer_cast< PanelledRadiationPressureInterfaceSettings >(
                    radiationPressureInterfaceSettings );
        if( panelledSettings == NULL )
        {
            throw std::runtime_error( "Error when making panelled radiation interface, type does not match object" );
        }

        // Retrieve source body and check consistency.
        if( bodyMap.count( radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw std::runtime_error( "Error when making panelled radiation interface, source not found.");
        }

        std::shared_ptr< Body > sourceBody =
                bodyMap.at( radiationPressureInterfaceSettings->getSourceBody( ) );

        // Get required data for occulting bodies.
        std::vector< std::string > occultingBodies = panelledSettings->getOccultingBodies( );
        std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions;
        std::vector< double > occultingBodyRadii;
        getOccultingBodiesInformation( bodyMap, occultingBodies, occultingBodyPositions, occultingBodyRadii );

        // Retrieve radius of source if occultations are used.
        double sourceRadius;
        if( occultingBodyPositions.size( ) > 0 )
        {
            std::shared_ptr< basic_astrodynamics::BodyShapeModel > sourceShapeModel =
                    sourceBody->getShapeModel( );

            if( sourceShapeModel == NULL )
            {
                throw std::runtime_error(
                            "Error when making occulted body, source body " + radiationPressureInterfaceSettings->getSourceBody( ) +
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
        std::function< double( ) > radiatedPowerFunction;
        if( defaultRadiatedPowerValues.count(
                    radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw std::runtime_error( "Error, no radiated power found for " +
                                      radiationPressureInterfaceSettings->getSourceBody( ) );
        }
        else
        {
            radiatedPowerFunction = [ = ]( ){ return defaultRadiatedPowerValues.at(
                            radiationPressureInterfaceSettings->getSourceBody( ) ); };
        }

        std::vector< std::function< Eigen::Vector3d( const double ) > > localFrameSurfaceNormalFunctions = panelledSettings->getSurfaceNormalsInBodyFixedFrameFunctions();


        // Create radiation pressure interface.
        radiationPressureInterface =
                std::make_shared< electro_magnetism::PanelledRadiationPressureInterface >(
                    radiatedPowerFunction,
                    std::bind( &Body::getPosition, sourceBody ),
                    std::bind( &Body::getPosition, bodyMap.at( bodyName ) ),
                    localFrameSurfaceNormalFunctions,
                    panelledSettings->getEmissivities( ),
                    panelledSettings->getAreas( ),
                    panelledSettings->getDiffusionCoefficients( ),
                    std::bind( &Body::getCurrentRotationToGlobalFrame, bodyMap.at( bodyName ) ),
                    occultingBodyPositions, occultingBodyRadii,
                    sourceRadius );
        break;
    }
    case solar_sailing_radiation_pressure_interface:
    {
        // Check type consistency.
        std::shared_ptr< SolarSailRadiationInterfaceSettings > solarSailRadiationSettings =
                std::dynamic_pointer_cast< SolarSailRadiationInterfaceSettings >( radiationPressureInterfaceSettings );

        if( solarSailRadiationSettings == nullptr )
        {
            throw std::runtime_error( "Error when making solar sail radiation interface, type does not match object" );
        }

        // Retrieve source body and check consistency.
        if( bodyMap.count( radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw std::runtime_error( "Error when making solar sail radiation interface, source not found.");
        }
        std::shared_ptr< Body > sourceBody = bodyMap.at( radiationPressureInterfaceSettings->getSourceBody( ) );


        // Get required data for occulting bodies.
        std::vector< std::string > occultingBodies = solarSailRadiationSettings->getOccultingBodies( );
        std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions;
        std::vector< double > occultingBodyRadii;
        getOccultingBodiesInformation(
            bodyMap, occultingBodies, occultingBodyPositions, occultingBodyRadii );

        // Retrieve radius of source if occultations are used.
        double sourceRadius;
        if( occultingBodyPositions.size( ) > 0 )
        {
            std::shared_ptr< basic_astrodynamics::BodyShapeModel > sourceShapeModel =
                sourceBody->getShapeModel( );

            if( sourceShapeModel == nullptr )
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


        // Get required data for central bodies.
        std::string centralBody = solarSailRadiationSettings->getCentralBody();
        std::function< Eigen::Vector3d( ) > centralBodyPosition;
        std::function< Eigen::Vector3d( ) > centralBodyVelocity;
        getCentralBodyInformation( bodyMap, centralBody, centralBodyPosition, centralBodyVelocity );

        // Create function returning radiated power.
        std::function< double( ) > radiatedPowerFunction;
        if( defaultRadiatedPowerValues.count(
                radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw std::runtime_error( "Error, no radiated power found for " +
                                     radiationPressureInterfaceSettings->getSourceBody( ) );
        }
        else
        {
            radiatedPowerFunction = [ = ]( ){ return defaultRadiatedPowerValues.at( radiationPressureInterfaceSettings->getSourceBody( ) );};
        }

        // Create solar sailing radiation pressure interface.
        radiationPressureInterface =
            std::make_shared< electro_magnetism::SolarSailingRadiationPressureInterface >(
                radiatedPowerFunction,
                std::bind( &Body::getPosition, sourceBody ),
                std::bind( &Body::getPosition, bodyMap.at( bodyName ) ),
                std::bind( &Body::getVelocity, bodyMap.at( bodyName ) ),
                solarSailRadiationSettings->getArea( ),
                solarSailRadiationSettings->getConeAngle(),
                solarSailRadiationSettings->getClockAngle(),
                solarSailRadiationSettings->getFrontEmissivityCoefficient(),
                solarSailRadiationSettings->getBackEmissivityCoefficient(),
                solarSailRadiationSettings->getFrontLambertianCoefficient(),
                solarSailRadiationSettings->getBackLambertianCoefficient(),
                solarSailRadiationSettings->getReflectivityCoefficient(),
                solarSailRadiationSettings->getSpecularReflectionCoefficient(),
                occultingBodyPositions, centralBodyVelocity,
                occultingBodyRadii, sourceRadius);
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
