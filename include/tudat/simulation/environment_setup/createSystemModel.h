/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATESSYTEMMODEL_H
#define TUDAT_CREATESSYTEMMODEL_H

#include <memory>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/thrustSettings.h"

namespace tudat
{

namespace simulation_setup
{

enum BodyPanelPropertyType
{
    body_panel_reflection_law
};

class BodyPanelPropertySettings
{
public:

    virtual ~BodyPanelPropertySettings( ){ }

    virtual BodyPanelPropertyType getBodyPanelPropertyType( ) = 0;

};

class BodyPanelGeometry
{
public:

    BodyPanelGeometry( const std::string& frameOrientation ):
        frameOrientation_( frameOrientation ){ }

    virtual ~BodyPanelGeometry( ){ }

    std::string frameOrientation_;
};

class FrameFixedBodyPanelGeometry: public BodyPanelGeometry
{
public:

    FrameFixedBodyPanelGeometry(
        const Eigen::Vector3d& surfaceNormal,
        const double area,
        const std::string& frameOrientation ):BodyPanelGeometry( frameOrientation ),
        surfaceNormal_( surfaceNormal ), area_( area ){ }

    Eigen::Vector3d surfaceNormal_;

    double area_;
};

class FrameVariableBodyPanelGeometry: public BodyPanelGeometry
{
public:

    FrameVariableBodyPanelGeometry(
        const std::function< Eigen::Vector3d( ) >& surfaceNormalFunction,
        const double area,
        const std::string& frameOrientation ):
        BodyPanelGeometry( frameOrientation ), surfaceNormalFunction_( surfaceNormalFunction ), area_( area ),
        bodyToTrack_( std::make_pair( nullptr, 0 ) ){ }

    FrameVariableBodyPanelGeometry(
        const std::string& bodyToTrack,
        const bool towardsTrackedBody,
        const double area,
        const std::string& frameOrientation ):
        BodyPanelGeometry( frameOrientation ), surfaceNormalFunction_( nullptr ), area_( area ),
        bodyToTrack_( std::make_pair( bodyToTrack, towardsTrackedBody ) ){ }

    std::function< Eigen::Vector3d( ) > surfaceNormalFunction_;

    double area_;

    std::pair< std::string, bool > bodyToTrack_;

};


class BodyPanelSettings
{
public:
    BodyPanelSettings(
        const std::shared_ptr< BodyPanelGeometry > panelGeometry,
        const std::shared_ptr< BodyPanelPropertySettings > panelPropertySettings,
        const std::string panelTypeId = "" ):
        panelGeometry_( panelGeometry )
    {
        panelProperties_[ panelPropertySettings->getBodyPanelPropertyType( ) ] = panelPropertySettings;
        panelTypeId_ = panelTypeId;
    }

    std::shared_ptr< BodyPanelGeometry > panelGeometry_;

    std::map< BodyPanelPropertyType, std::shared_ptr< BodyPanelPropertySettings > > panelProperties_;

    std::string panelTypeId_;
};

void addEngineModel(
        const std::string& bodyName,
        const std::string& engineName,
        const std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) );

void addVariableDirectionEngineModel(
        const std::string& bodyName,
        const std::string& engineName,
        const std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::function< Eigen::Vector3d( const double ) > bodyFixedThrustDirection );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATESSYTEMMODEL_H
