/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */





#include "tudat/simulation/environment_setup/createRadiationPressureInterface.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to obtain (by reference) the position functions and radii of occulting bodies
void getOccultingBodiesInformation(
        const SystemOfBodies& bodies, const std::vector< std::string >& occultingBodies,
        std::vector< std::function< Eigen::Vector3d( ) > >& occultingBodyPositions,
        std::vector< double >& occultingBodyRadii )
{
    // Iterate over occulring bodies and retrieve radius and position function.
    for( unsigned int i = 0; i < occultingBodies.size( ); i++ )
    {
        if( bodies.count( occultingBodies[ i ] ) == 0 )
        {
            throw std::runtime_error( "Error, could not find body " + occultingBodies[ i ] +
                                      " in system of bodies when making occulting body settings" );
        }
        else
        {
            std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel =
                    bodies.at( occultingBodies.at( i ) )->getShapeModel( );
            if( shapeModel == nullptr )
            {
                throw std::runtime_error( "Error, no shape model for " + occultingBodies[ i ] +
                                          " when making occulting body settings" );
            }
            occultingBodyPositions.push_back(
                        std::bind( &Body::getPosition, bodies.at( occultingBodies[ i ] ) ) );
            occultingBodyRadii.push_back( shapeModel->getAverageRadius( ) );
        }
    }
}


//! Function to obtain (by reference) the position and velocity functions of the central body.
void getCentralBodyInformation(
    const SystemOfBodies& bodies, const std::string& centralBody,
    std::function< Eigen::Vector3d( ) >& centralBodyPosition,
    std::function< Eigen::Vector3d( ) >& centralBodyVelocity )
{

    // Check that the central body is defined.
    if( bodies.count( centralBody ) == 0 )
    {
        throw std::runtime_error( "Error, could not find body " + centralBody +
                                 " in system of bodies when making central body body settings" );
    }
    // Retrieve position and velocity functions.
    else
    {
        centralBodyPosition = std::bind( &Body::getPosition, bodies.at( centralBody ) );
        centralBodyVelocity = std::bind( &Body::getVelocity, bodies.at( centralBody ) );
    }
}


//! Function to create a radiation pressure interface.
std::shared_ptr< electromagnetism::RadiationPressureInterface > createRadiationPressureInterface(
        const std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureInterfaceSettings,
        const std::string& bodyName, const SystemOfBodies& bodies )
{
    std::shared_ptr< electromagnetism::RadiationPressureInterface > radiationPressureInterface;

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
        if( bodies.count( radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw std::runtime_error( "Error when making cannon ball radiation interface, source not found.");
        }

        std::shared_ptr< Body > sourceBody =
                bodies.at( radiationPressureInterfaceSettings->getSourceBody( ) );

        // Get reqruied data for occulting bodies.
        std::vector< std::string > occultingBodies = cannonBallSettings->getOccultingBodies( );
        std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions;
        std::vector< double > occultingBodyRadii;
        getOccultingBodiesInformation(
                    bodies, occultingBodies, occultingBodyPositions, occultingBodyRadii );

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

        if( cannonBallSettings->getRadiationPressureCoefficient( ) !=
                cannonBallSettings->getRadiationPressureCoefficient( ) )
        {
        // Create radiation pressure interface.
        radiationPressureInterface =
                std::make_shared< electromagnetism::RadiationPressureInterface >(
                    radiatedPowerFunction,
                    std::bind( &Body::getPosition, sourceBody ),
                    std::bind( &Body::getPosition, bodies.at( bodyName ) ),
                    cannonBallSettings->getRadiationPressureCoefficientFunction( ),
                    cannonBallSettings->getArea( ), occultingBodyPositions, occultingBodyRadii,
                    sourceRadius,
                    occultingBodies );
        }
        else
        {
            // Create radiation pressure interface.
            radiationPressureInterface =
                    std::make_shared< electromagnetism::RadiationPressureInterface >(
                        radiatedPowerFunction,
                        std::bind( &Body::getPosition, sourceBody ),
                        std::bind( &Body::getPosition, bodies.at( bodyName ) ),
                        cannonBallSettings->getRadiationPressureCoefficient( ),
                        cannonBallSettings->getArea( ), occultingBodyPositions, occultingBodyRadii,
                        sourceRadius,
                        occultingBodies );
        }
        break;
    }
//    case panelled_radiation_pressure_interface:
//    {
//        // Check type consistency.
//        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > panelledSettings =
//                std::dynamic_pointer_cast< PanelledRadiationPressureInterfaceSettings >(
//                    radiationPressureInterfaceSettings );
//        if( panelledSettings == NULL )
//        {
//            throw std::runtime_error( "Error when making panelled radiation interface, type does not match object" );
//        }
//
//        // Retrieve source body and check consistency.
//        if( bodies.count( radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
//        {
//            throw std::runtime_error( "Error when making panelled radiation interface, source not found.");
//        }
//
//        std::shared_ptr< Body > sourceBody =
//                bodies.at( radiationPressureInterfaceSettings->getSourceBody( ) );
//
//        // Get required data for occulting bodies.
//        std::vector< std::string > occultingBodies = panelledSettings->getOccultingBodies( );
//        std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions;
//        std::vector< double > occultingBodyRadii;
//        getOccultingBodiesInformation( bodies, occultingBodies, occultingBodyPositions, occultingBodyRadii );
//
//        // Retrieve radius of source if occultations are used.
//        double sourceRadius;
//        if( occultingBodyPositions.size( ) > 0 )
//        {
//            std::shared_ptr< basic_astrodynamics::BodyShapeModel > sourceShapeModel =
//                    sourceBody->getShapeModel( );
//
//            if( sourceShapeModel == NULL )
//            {
//                throw std::runtime_error(
//                            "Error when making occulted body, source body " + radiationPressureInterfaceSettings->getSourceBody( ) +
//                            " does not have a shape" );
//            }
//            else
//            {
//                sourceRadius = sourceShapeModel->getAverageRadius( );
//            }
//        }
//        else
//        {
//            sourceRadius = 0.0;
//        }
//
//        // Create function returning radiated power.
//        std::function< double( ) > radiatedPowerFunction;
//        if( defaultRadiatedPowerValues.count(
//                    radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
//        {
//            throw std::runtime_error( "Error, no radiated power found for " +
//                                      radiationPressureInterfaceSettings->getSourceBody( ) );
//        }
//        else
//        {
//            radiatedPowerFunction = [ = ]( ){ return defaultRadiatedPowerValues.at(
//                            radiationPressureInterfaceSettings->getSourceBody( ) ); };
//        }
//
//        std::vector< std::function< Eigen::Vector3d( const double ) > > localFrameSurfaceNormalFunctions = panelledSettings->getSurfaceNormalsInBodyFixedFrameFunctions();
//
//
//        // Create radiation pressure interface.
//        radiationPressureInterface =
//                std::make_shared< electromagnetism::PanelledRadiationPressureInterface >(
//                    radiatedPowerFunction,
//                    std::bind( &Body::getPosition, sourceBody ),
//                    std::bind( &Body::getPosition, bodies.at( bodyName ) ),
//                    localFrameSurfaceNormalFunctions,
//                    panelledSettings->getEmissivities( ),
//                    panelledSettings->getAreas( ),
//                    panelledSettings->getDiffusionCoefficients( ),
//                    std::bind( &Body::getCurrentRotationToGlobalFrame, bodies.at( bodyName ) ),
//                    occultingBodyPositions, occultingBodyRadii,
//                    sourceRadius );
//        break;
//    }
//    case solar_sailing_radiation_pressure_interface:
//    {
//        // Check type consistency.
//        std::shared_ptr< SolarSailRadiationInterfaceSettings > solarSailRadiationSettings =
//                std::dynamic_pointer_cast< SolarSailRadiationInterfaceSettings >( radiationPressureInterfaceSettings );
//
//        if( solarSailRadiationSettings == nullptr )
//        {
//            throw std::runtime_error( "Error when making solar sail radiation interface, type does not match object" );
//        }
//
//        // Retrieve source body and check consistency.
//        if( bodies.count( radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
//        {
//            throw std::runtime_error( "Error when making solar sail radiation interface, source not found.");
//        }
//        std::shared_ptr< Body > sourceBody = bodies.at( radiationPressureInterfaceSettings->getSourceBody( ) );
//
//
//        // Get required data for occulting bodies.
//        std::vector< std::string > occultingBodies = solarSailRadiationSettings->getOccultingBodies( );
//        std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions;
//        std::vector< double > occultingBodyRadii;
//        getOccultingBodiesInformation(
//            bodies, occultingBodies, occultingBodyPositions, occultingBodyRadii );
//
//        // Retrieve radius of source if occultations are used.
//        double sourceRadius;
//        if( occultingBodyPositions.size( ) > 0 )
//        {
//            std::shared_ptr< basic_astrodynamics::BodyShapeModel > sourceShapeModel =
//                sourceBody->getShapeModel( );
//
//            if( sourceShapeModel == nullptr )
//            {
//                throw std::runtime_error( "Error when making occulted body, source body " +
//                                         radiationPressureInterfaceSettings->getSourceBody( ) +
//                                         " does not have a shape" );
//            }
//            else
//            {
//                sourceRadius = sourceShapeModel->getAverageRadius( );
//            }
//        }
//        else
//        {
//            sourceRadius = 0.0;
//        }
//
//
//        // Get required data for central bodies.
//        std::string centralBody = solarSailRadiationSettings->getCentralBody();
//        std::function< Eigen::Vector3d( ) > centralBodyPosition;
//        std::function< Eigen::Vector3d( ) > centralBodyVelocity;
//        getCentralBodyInformation( bodies, centralBody, centralBodyPosition, centralBodyVelocity );
//
//        // Create function returning radiated power.
//        std::function< double( ) > radiatedPowerFunction;
//        if( defaultRadiatedPowerValues.count(
//                radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
//        {
//            throw std::runtime_error( "Error, no radiated power found for " +
//                                     radiationPressureInterfaceSettings->getSourceBody( ) );
//        }
//        else
//        {
//            radiatedPowerFunction = [ = ]( ){ return defaultRadiatedPowerValues.at( radiationPressureInterfaceSettings->getSourceBody( ) );};
//        }
//
//        // Create solar sailing radiation pressure interface.
//        radiationPressureInterface =
//            std::make_shared< electromagnetism::SolarSailingRadiationPressureInterface >(
//                radiatedPowerFunction,
//                std::bind( &Body::getPosition, sourceBody ),
//                std::bind( &Body::getPosition, bodies.at( bodyName ) ),
//                std::bind( &Body::getVelocity, bodies.at( bodyName ) ),
//                solarSailRadiationSettings->getArea( ),
//                solarSailRadiationSettings->getConeAngle(),
//                solarSailRadiationSettings->getClockAngle(),
//                solarSailRadiationSettings->getFrontEmissivityCoefficient(),
//                solarSailRadiationSettings->getBackEmissivityCoefficient(),
//                solarSailRadiationSettings->getFrontLambertianCoefficient(),
//                solarSailRadiationSettings->getBackLambertianCoefficient(),
//                solarSailRadiationSettings->getReflectivityCoefficient(),
//                solarSailRadiationSettings->getSpecularReflectionCoefficient(),
//                occultingBodyPositions, centralBodyVelocity,
//                occultingBodyRadii, sourceRadius);
//        break;
//    }
    default:
        throw std::runtime_error(
                    "Error, radiation pressure type" + std::to_string(
                        radiationPressureInterfaceSettings->getRadiationPressureType( ) ) +
                    "not recognized for body" + bodyName );
    }

    return radiationPressureInterface;
}

std::function< double( const double ) > getOccultationFunction(
        const SystemOfBodies& bodyMap,
        const std::string& sourceBody,
        const std::string& occultingBody,
        const std::string& shadowedBody )
{
    std::function< Eigen::Vector3d( ) > sourceBodyPositionFunction =
            std::bind( &Body::getPosition, bodyMap.at( sourceBody ) );
    std::function< Eigen::Vector3d( ) > occultingBodyPositionFunction =
            std::bind( &Body::getPosition, bodyMap.at( occultingBody ) );
    std::function< Eigen::Vector3d( ) > shadowedBodyPositionFunction =
            std::bind( &Body::getPosition, bodyMap.at( shadowedBody ) );
    double sourceBodyRadius = bodyMap.at( sourceBody )->getShapeModel( )->getAverageRadius( );
    double occultingBodyRadius = bodyMap.at( occultingBody )->getShapeModel( )->getAverageRadius( );

    return [=]( const double ){
        return mission_geometry::computeShadowFunction(
                    sourceBodyPositionFunction( ), sourceBodyRadius,
                    occultingBodyPositionFunction( ), occultingBodyRadius,
                    shadowedBodyPositionFunction( ) ); };
}


} // namespace simulation_setup

} // namespace tudat
