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
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"

namespace tudat
{

namespace simulation_setup
{

class BodyPanelGeometrySettings
{
public:

    BodyPanelGeometrySettings( const std::string& frameOrientation = "" ):
        frameOrientation_( frameOrientation ){ }

    virtual ~BodyPanelGeometrySettings( ){ }

    std::string frameOrientation_;
};

class FrameFixedBodyPanelGeometrySettings: public BodyPanelGeometrySettings
{
public:

    FrameFixedBodyPanelGeometrySettings(
        const Eigen::Vector3d& surfaceNormal,
        const double area,
        const std::string& frameOrientation = "" ):BodyPanelGeometrySettings( frameOrientation ),
        surfaceNormal_( surfaceNormal.normalized( ) ), area_( area ){ }

    Eigen::Vector3d surfaceNormal_;

    double area_;
};

class FrameVariableBodyPanelGeometrySettings: public BodyPanelGeometrySettings
{
public:

    FrameVariableBodyPanelGeometrySettings(
        const std::function< Eigen::Vector3d( ) >& surfaceNormalFunction,
        const double area,
        const std::string& frameOrientation = "" ):
        BodyPanelGeometrySettings( frameOrientation ), surfaceNormalFunction_( surfaceNormalFunction ), area_( area ),
        bodyToTrack_( std::make_pair( "", 0 ) ){ }

    FrameVariableBodyPanelGeometrySettings(
        const std::string& bodyToTrack,
        const bool towardsTrackedBody,
        const double area,
        const std::string& frameOrientation = "" ):
        BodyPanelGeometrySettings( frameOrientation ), surfaceNormalFunction_( nullptr ), area_( area ),
        bodyToTrack_( std::make_pair( bodyToTrack, towardsTrackedBody ) ){ }

    std::function< Eigen::Vector3d( ) > surfaceNormalFunction_;

    double area_;

    std::pair< std::string, bool > bodyToTrack_;

};

inline std::shared_ptr< BodyPanelGeometrySettings > frameFixedPanelGeometry(
    const Eigen::Vector3d& surfaceNormal,
    const double area,
    const std::string& frameOrientation = "" )
{
    return std::make_shared< FrameFixedBodyPanelGeometrySettings >(
        surfaceNormal, area, frameOrientation );
}

inline std::shared_ptr< BodyPanelGeometrySettings > timeVaryingPanelGeometry(
    const std::function< Eigen::Vector3d( ) >& surfaceNormalFunction,
    const double area,
    const std::string& frameOrientation = "" )
{
    return std::make_shared< FrameVariableBodyPanelGeometrySettings >(
        surfaceNormalFunction, area, frameOrientation );
}

inline std::shared_ptr< BodyPanelGeometrySettings > bodyTrackingPanelGeometry(
    const std::string& bodyToTrack,
    const bool towardsTrackedBody,
    const double area,
    const std::string& frameOrientation = "" )
{
    return std::make_shared< FrameVariableBodyPanelGeometrySettings >(
        bodyToTrack, towardsTrackedBody, area, frameOrientation );
}


class BodyPanelSettings
{
public:
    BodyPanelSettings(
        const std::shared_ptr< BodyPanelGeometrySettings > panelGeometry,
        const std::shared_ptr< BodyPanelReflectionLawSettings > reflectionLawSettings,
        const std::string panelTypeId = "" ):
        panelGeometry_( panelGeometry ),
        reflectionLawSettings_( reflectionLawSettings ),
        panelTypeId_( panelTypeId ){ }

    std::shared_ptr< BodyPanelGeometrySettings > panelGeometry_;

    std::shared_ptr< BodyPanelReflectionLawSettings > reflectionLawSettings_;

    std::string panelTypeId_;
};


inline std::shared_ptr< BodyPanelSettings > bodyPanelSettings(
    const std::shared_ptr< BodyPanelGeometrySettings > panelGeometry,
    const std::shared_ptr< BodyPanelReflectionLawSettings > reflectionLawSettings = nullptr,
    const std::string panelTypeId = "" )
{
    return std::make_shared< BodyPanelSettings >( panelGeometry, reflectionLawSettings, panelTypeId );
}

class FullPanelledBodySettings
{
public:
    FullPanelledBodySettings(
        const std::vector< std::shared_ptr< BodyPanelSettings > >&  panelSettingsList,
        const std::map< std::string, std::shared_ptr< RotationModelSettings > >& partRotationModelSettings  =
            std::map< std::string, std::shared_ptr< RotationModelSettings > >( )):
        panelSettingsList_( panelSettingsList ), partRotationModelSettings_( partRotationModelSettings ){ }

    void addPanelRotationSettings(
        const std::shared_ptr< RotationModelSettings > rotationSettings,
        const std::string partName )
    {
        partRotationModelSettings_[ partName ] = rotationSettings;
    }

    std::vector< std::shared_ptr< BodyPanelSettings > >  panelSettingsList_;

    std::map< std::string, std::shared_ptr< RotationModelSettings > > partRotationModelSettings_;
};


inline std::shared_ptr< FullPanelledBodySettings > fullPanelledBodySettings(
    const std::vector< std::shared_ptr< BodyPanelSettings > >&  panelSettingsList,
    const std::map< std::string, std::shared_ptr< RotationModelSettings > >& partRotationModelSettings  =
    std::map< std::string, std::shared_ptr< RotationModelSettings > >( ) )
{
    return std::make_shared< FullPanelledBodySettings >( panelSettingsList, partRotationModelSettings );
}

inline std::shared_ptr< FullPanelledBodySettings > bodyWingPanelledGeometry(
    const double length,
    const double width,
    const double height,
    const double totalSolarArrayArea,
    const double busSpecularReflectivity,
    const double busDiffuseReflectivity,
    const double solarArraySpecularReflectivity,
    const double solarArrayDiffuseReflectivity,
    const bool busInstantaneousReradiation = true,
    const bool solarArrayInstantaneousReradiation = true )
{
    std::vector< std::shared_ptr< BodyPanelGeometrySettings > > panelGeometrySettingsList =
        std::vector< std::shared_ptr< BodyPanelGeometrySettings > >(
            { frameFixedPanelGeometry( Eigen::Vector3d::UnitX( ), width * height ),
              frameFixedPanelGeometry( -Eigen::Vector3d::UnitX( ), width * height ),
              frameFixedPanelGeometry( Eigen::Vector3d::UnitY( ), length * height ),
              frameFixedPanelGeometry( -Eigen::Vector3d::UnitY( ), length * height ),
              frameFixedPanelGeometry( Eigen::Vector3d::UnitZ( ), length * width  ),
              frameFixedPanelGeometry( -Eigen::Vector3d::UnitZ( ), length * width ),
              bodyTrackingPanelGeometry( "Sun", true, totalSolarArrayArea ),
              bodyTrackingPanelGeometry( "Sun", false, totalSolarArrayArea ) } );
    std::vector< std::shared_ptr< BodyPanelSettings > > panelSettings;
    for( unsigned int i = 0; i < 6; i++ )
    {
        panelSettings.push_back(
            bodyPanelSettings( panelGeometrySettingsList.at( i ), specularDiffuseBodyPanelReflectionLawSettings(
                busSpecularReflectivity, busDiffuseReflectivity, busInstantaneousReradiation ), "Bus" ) );
    }

    for( unsigned int i = 0; i < 2; i++ )
    {
        panelSettings.push_back(
            bodyPanelSettings( panelGeometrySettingsList.at( i + 6 ), specularDiffuseBodyPanelReflectionLawSettings(
                solarArraySpecularReflectivity, solarArrayDiffuseReflectivity, solarArrayInstantaneousReradiation ), "SolarPanel" ) );
    }

    return fullPanelledBodySettings( panelSettings );
}



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

std::pair< std::shared_ptr< system_models::VehicleExteriorPanel >, std::string > createBodyExteriorPanel(
    const std::shared_ptr< BodyPanelSettings > panelSettings,
    const std::string& bodyName,
    const simulation_setup::SystemOfBodies& bodies );

void addBodyExteriorPanelledShape(
    const std::shared_ptr< FullPanelledBodySettings > panelSettings,
    const std::string& bodyName,
    const simulation_setup::SystemOfBodies& bodies );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATESSYTEMMODEL_H
