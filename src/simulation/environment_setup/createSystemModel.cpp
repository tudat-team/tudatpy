/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/simulation/environment_setup/createThrustModelGuidance.h"
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"

namespace tudat
{

namespace simulation_setup
{



void addEngineModel(
        const std::string& bodyName,
        const std::string& engineName,
        const std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const Eigen::Vector3d bodyFixedThrustDirection )
{
    addVariableDirectionEngineModel( bodyName, engineName, thrustSettings, bodies,
                                     [=](const double){return bodyFixedThrustDirection; } );
}


void addVariableDirectionEngineModel(
        const std::string& bodyName,
        const std::string& engineName,
        const std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::function< Eigen::Vector3d( const double ) > bodyFixedThrustDirection )
{
    {
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > magnitudeUpdateSettings;
        std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper = createThrustMagnitudeWrapper(
                    thrustSettings,
                    bodies, bodyName, magnitudeUpdateSettings );


        std::shared_ptr< system_models::EngineModel > vehicleEngineModel =
                std::make_shared< system_models::EngineModel >(
                    thrustMagnitudeWrapper, engineName, bodyFixedThrustDirection);

        if( bodies.at( bodyName )->getVehicleSystems( ) == nullptr )
        {
            std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared<
                    system_models::VehicleSystems >(  );
            bodies.at( bodyName )->setVehicleSystems( vehicleSystems );
        }
        bodies.at( bodyName )->getVehicleSystems( )->setEngineModel( vehicleEngineModel );

    }
}

std::pair< std::shared_ptr< system_models::VehicleExteriorPanel >, std::string > createBodyExteriorPanel(
    const std::shared_ptr< BodyPanelSettings > panelSettings,
    const std::string& bodyName,
    const simulation_setup::SystemOfBodies& bodies )
{
    std::string panelFixedFrame = "";
    if( panelSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating body exterior panel settings, no panel settings provided" );
    }

    if( panelSettings->panelGeometry_ == nullptr )
    {
        throw std::runtime_error( "Error when creating body exterior panel settings, no panel geometry settings provided" );
    }

    double panelArea = TUDAT_NAN;
    double panelTemperature = TUDAT_NAN;
    std::function< Eigen::Vector3d( ) > localFrameSurfaceNormal = nullptr;
    std::function< Eigen::Vector3d( ) > localFramePositionVector = nullptr;
    std::string trackedBodyName = "";

    if( std::dynamic_pointer_cast< FrameFixedBodyPanelGeometrySettings >( panelSettings->panelGeometry_ ) != nullptr )
    {
        std::shared_ptr< FrameFixedBodyPanelGeometrySettings > fixedBodyPanelGeometrySettings =
            std::dynamic_pointer_cast< FrameFixedBodyPanelGeometrySettings >( panelSettings->panelGeometry_ );
        panelArea = fixedBodyPanelGeometrySettings->area_;
        panelTemperature = fixedBodyPanelGeometrySettings->panelTemperature_;

        localFrameSurfaceNormal = [=]( ){ return fixedBodyPanelGeometrySettings->surfaceNormal_; };
        localFramePositionVector = [=]( ){ return fixedBodyPanelGeometrySettings->positionVector_; };
    }
    else if( std::dynamic_pointer_cast< FrameVariableBodyPanelGeometrySettings >( panelSettings->panelGeometry_ ) != nullptr )
    {
        std::shared_ptr< FrameVariableBodyPanelGeometrySettings > variableBodyPanelGeometrySettings =
            std::dynamic_pointer_cast< FrameVariableBodyPanelGeometrySettings >( panelSettings->panelGeometry_ );
        panelArea = variableBodyPanelGeometrySettings->area_;
        panelTemperature = variableBodyPanelGeometrySettings->panelTemperature_;
        if( variableBodyPanelGeometrySettings->surfaceNormalFunction_ != nullptr )
        {
            if( variableBodyPanelGeometrySettings->bodyToTrack_.first != "" )
            {
                throw std::runtime_error( "Error, explicit panel orientation, and tracking target are defined" );
            }
            localFrameSurfaceNormal = variableBodyPanelGeometrySettings->surfaceNormalFunction_;
            localFramePositionVector = variableBodyPanelGeometrySettings->positionVectorFunction_; 
        }
        else if( variableBodyPanelGeometrySettings->bodyToTrack_.first != "" )
        {
            // Tracking a body means setting the surface normal towards the tracked body in the source local frame
            trackedBodyName = variableBodyPanelGeometrySettings->bodyToTrack_.first;
            const auto bodyToTrack = bodies.at(trackedBodyName);
            const auto targetBody = bodies.at(bodyName);
            const auto sign = variableBodyPanelGeometrySettings->bodyToTrack_.second ? +1.0 : -1.0;

            // Construct surface normal function always pointing towards/away from tracked body
            localFrameSurfaceNormal = [=] () {
                const Eigen::Quaterniond rotationFromPropagationToLocalFrame =
                    targetBody->getCurrentRotationToLocalFrame();
                const Eigen::Vector3d relativeSourcePositionInPropagationFrame =
                    bodyToTrack->getPosition() - targetBody->getPosition();
                const Eigen::Vector3d relativeSourcePositionInLocalFrame =
                    rotationFromPropagationToLocalFrame * relativeSourcePositionInPropagationFrame;
                Eigen::Vector3d surfaceNormal = sign * relativeSourcePositionInLocalFrame.normalized();
                return surfaceNormal;
            };
        }
        else
        {
            throw std::runtime_error( "Error, neither explicit panel orientation, nor tracking target are defined" );
        }
    }

    std::shared_ptr< system_models::VehicleExteriorPanel > exteriorPanel =
        std::make_shared< system_models::VehicleExteriorPanel >( 
            localFrameSurfaceNormal,
            localFramePositionVector,
            panelArea, 
            panelTemperature,
            trackedBodyName);

    if( panelSettings->reflectionLawSettings_ != nullptr )
    {
        std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = createReflectionLaw(
            panelSettings->reflectionLawSettings_ );
        exteriorPanel->setReflectionLaw( reflectionLaw );
        if( panelSettings->panelTypeId_ != "" )
        {
            exteriorPanel->setPanelTypeId( panelSettings->panelTypeId_ );
        }
    }

    if( bodies.at( bodyName )->getRotationalEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when creating body exterior panel model for body " + bodyName + ", no body rotation model provided" );
    }

    if( ( panelSettings->panelGeometry_->frameOrientation_ != "" ) &&
        ( panelSettings->panelGeometry_->frameOrientation_ != bodies.at( bodyName )->getRotationalEphemeris( )->getTargetFrameOrientation( ) ) )
    {
        panelFixedFrame = panelSettings->panelGeometry_->frameOrientation_;
    }

    return std::make_pair( exteriorPanel, panelFixedFrame );
}


void addBodyExteriorPanelledShape(
    const std::shared_ptr< FullPanelledBodySettings > panelSettings,
    const std::string& bodyName,
    const simulation_setup::SystemOfBodies& bodies )
{
    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when creating body exterior paneled shape for body " + bodyName + ", body does not exist" );
    }

    if( panelSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating body exterior paneled shape for body " + bodyName + ", provided settings are empty" );
    }

    if( bodies.at( bodyName )->getVehicleSystems( ) == nullptr )
    {
        bodies.at( bodyName )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
    }

    std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > > vehicleExteriorPanels;

    for( unsigned int i= 0; i < panelSettings->panelSettingsList_.size( ); i++ )
    {
        std::pair< std::shared_ptr< system_models::VehicleExteriorPanel >, std::string > currentPanel = createBodyExteriorPanel(
            panelSettings->panelSettingsList_.at( i ), bodyName, bodies );
        vehicleExteriorPanels[ currentPanel.second ].push_back( currentPanel.first );
    }

    std::map< std::string, std::shared_ptr< ephemerides::RotationalEphemeris > > vehiclePartOrientation;
    for( auto it : panelSettings->partRotationModelSettings_ )
    {
        if( it.second->getOriginalFrame( ) == "" )
        {
            it.second->resetOriginalFrame( bodies.at( bodyName )->getRotationalEphemeris( )->getTargetFrameOrientation( ) );
        }
        vehiclePartOrientation[ it.first ] = createRotationModel(
            it.second, bodyName + "_" + it.first, bodies );
    }
    bodies.at( bodyName )->getVehicleSystems( )->setVehicleExteriorPanels( vehicleExteriorPanels );
    bodies.at( bodyName )->getVehicleSystems( )->setVehiclePartOrientation( vehiclePartOrientation );

}

} // namespace simulation_setup

} // namespace tudat

