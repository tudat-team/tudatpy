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
#include <iostream>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/thrustSettings.h"
#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"

namespace tudat
{

namespace simulation_setup
{

class MaterialProperties
{
public:
    MaterialProperties( const double specularReflectivity, const double diffuseReflectivity,
                        const double energyAccomodationCoefficient, const double normalAccomodationCoefficient, 
                        const double tangentialAccomodationCoefficient, const double normalVelocityAtWallRatio ):
                        specularReflectivity_( specularReflectivity ), diffuseReflectivity_( diffuseReflectivity ),
                        energyAccomodationCoefficient_( energyAccomodationCoefficient ), 
                        normalAccomodationCoefficient_( normalAccomodationCoefficient ),
                        tangentialAccomodationCoefficient_( tangentialAccomodationCoefficient ),
                        normalVelocityAtWallRatio_( normalVelocityAtWallRatio )
    { }

    double specularReflectivity_;    
    double diffuseReflectivity_;
    double energyAccomodationCoefficient_;
    double normalAccomodationCoefficient_;
    double tangentialAccomodationCoefficient_;
    double normalVelocityAtWallRatio_;
};

inline std::shared_ptr< MaterialProperties > materialProperties( const double specularReflectivity, const double diffuseReflectivity,
                    const double energyAccomodationCoefficient, const double normalAccomodationCoefficient, 
                    const double tangentialAccomodationCoefficient, const double normalVelocityAtWallRatio )
{
    return std::make_shared< MaterialProperties >( specularReflectivity, diffuseReflectivity, energyAccomodationCoefficient,
                                                   normalAccomodationCoefficient, tangentialAccomodationCoefficient, normalVelocityAtWallRatio );
}

class BodyPanelGeometrySettings
{
public:
    BodyPanelGeometrySettings( const std::string& frameOrientation = "" ): frameOrientation_( frameOrientation ) { }

    virtual ~BodyPanelGeometrySettings( ) { }

    std::string frameOrientation_;

};

class FrameFixedBodyPanelGeometrySettings : public BodyPanelGeometrySettings
{
public:
    FrameFixedBodyPanelGeometrySettings( const Eigen::Vector3d& surfaceNormal,
                                         const Eigen::Vector3d& positionVector,
                                         const double area,
                                         const double panelTemperature = 273.0,
                                         const std::string& frameOrientation = "" ):
        BodyPanelGeometrySettings( frameOrientation ), surfaceNormal_( surfaceNormal.normalized( ) ), positionVector_( positionVector ),
        area_( area ), panelTemperature_( panelTemperature ), geometry3dLoaded_( false )
    { }

    FrameFixedBodyPanelGeometrySettings( const Eigen::Vector3d& surfaceNormal,
                                         const double area,
                                         const std::string& frameOrientation = "" ):
        BodyPanelGeometrySettings( frameOrientation ), surfaceNormal_( surfaceNormal.normalized( ) ), area_( area ), geometry3dLoaded_( false )
    { }

    // added constructor for triangle
    FrameFixedBodyPanelGeometrySettings( const Eigen::Vector3d vertexA, const Eigen::Vector3d vertexB, const Eigen::Vector3d vertexC,
                                         const Eigen::Vector3d& surfaceNormal, const double area,
                                         const Eigen::Vector3d frameOrigin = Eigen::Vector3d::Zero( ), const std::string& frameOrientation = "" ):
        BodyPanelGeometrySettings( frameOrientation ), surfaceNormal_( surfaceNormal ), area_( area ), 
        vertexA_( vertexA ), vertexB_( vertexB ), vertexC_( vertexC ), 
        frameOrigin_( frameOrigin ), geometry3dLoaded_( true )
    {   
        Eigen::Vector3d positionVector;
        positionVector(0) = (vertexA_(0)+vertexB_(0)+vertexC_(0))/3;
        positionVector(1) = (vertexA_(1)+vertexB_(1)+vertexC_(1))/3;
        positionVector(2) = (vertexA_(2)+vertexB_(2)+vertexC_(2))/3;

        positionVector_ = positionVector;
    }

    Eigen::Vector3d surfaceNormal_;
    Eigen::Vector3d positionVector_;

    double area_;
    double panelTemperature_;

    Eigen::Vector3d vertexA_;
    Eigen::Vector3d vertexB_;
    Eigen::Vector3d vertexC_;
    Eigen::Vector3d frameOrigin_;
    bool geometry3dLoaded_;
};

class FrameVariableBodyPanelGeometrySettings : public BodyPanelGeometrySettings
{
public:
    FrameVariableBodyPanelGeometrySettings( const std::function< Eigen::Vector3d( ) >& surfaceNormalFunction,
                                            const std::function< Eigen::Vector3d( ) >& positionVectorFunction,
                                            const double area,
                                            const double panelTemperature = 273.0,
                                            const std::string& frameOrientation = "" ):
        BodyPanelGeometrySettings( frameOrientation ), surfaceNormalFunction_( surfaceNormalFunction ),
        positionVectorFunction_( positionVectorFunction ), area_( area ), panelTemperature_( panelTemperature ),
        bodyToTrack_( std::make_pair( "", 0 ) )
    { }

    FrameVariableBodyPanelGeometrySettings( const std::string& bodyToTrack,
                                            const bool towardsTrackedBody,
                                            const std::function< Eigen::Vector3d( ) >& positionVectorFunction,
                                            const double area,
                                            const double panelTemperature = 273.0,
                                            const std::string& frameOrientation = "" ):
        BodyPanelGeometrySettings( frameOrientation ), surfaceNormalFunction_( nullptr ), positionVectorFunction_( positionVectorFunction ),
        area_( area ), panelTemperature_( panelTemperature ), bodyToTrack_( std::make_pair( bodyToTrack, towardsTrackedBody ) )
    { }

    FrameVariableBodyPanelGeometrySettings( const std::function< Eigen::Vector3d( ) >& surfaceNormalFunction,
                                            const double area,
                                            const std::string& frameOrientation = "" ):
        BodyPanelGeometrySettings( frameOrientation ), surfaceNormalFunction_( surfaceNormalFunction ), area_( area ),
        bodyToTrack_( std::make_pair( "", 0 ) )
    { }

    FrameVariableBodyPanelGeometrySettings( const std::string& bodyToTrack,
                                            const bool towardsTrackedBody,
                                            const double area,
                                            const std::string& frameOrientation = "" ):
        BodyPanelGeometrySettings( frameOrientation ), surfaceNormalFunction_( nullptr ), area_( area ),
        bodyToTrack_( std::make_pair( bodyToTrack, towardsTrackedBody ) )
    { }

    std::function< Eigen::Vector3d( ) > surfaceNormalFunction_;
    std::function< Eigen::Vector3d( ) > positionVectorFunction_;



    double area_;
    double panelTemperature_;

    std::pair< std::string, bool > bodyToTrack_;


};

inline std::shared_ptr< BodyPanelGeometrySettings > frameFixedPanelGeometry( const Eigen::Vector3d& surfaceNormal,
                                                                             const Eigen::Vector3d& positionVector,
                                                                             const double area,
                                                                             const double panelTemperature = 273.0,
                                                                             const std::string& frameOrientation = "" )
{
    return std::make_shared< FrameFixedBodyPanelGeometrySettings >(
            surfaceNormal, positionVector, area, panelTemperature, frameOrientation );
}

inline std::shared_ptr< BodyPanelGeometrySettings > timeVaryingPanelGeometry(
        const std::function< Eigen::Vector3d( ) >& surfaceNormalFunction,
        const std::function< Eigen::Vector3d( ) >& positionVectorFunction,
        const double area,
        const double panelTemperature = 273.0,
        const std::string& frameOrientation = "" )
{
    return std::make_shared< FrameVariableBodyPanelGeometrySettings >(
            surfaceNormalFunction, positionVectorFunction, area, panelTemperature, frameOrientation );
}

inline std::shared_ptr< BodyPanelGeometrySettings > bodyTrackingPanelGeometry(
        const std::string& bodyToTrack,
        const bool towardsTrackedBody,
        const std::function< Eigen::Vector3d( ) >& positionVectorFunction,
        const double area,
        const double panelTemperature = 273.0,
        const std::string& frameOrientation = "" )
{
    return std::make_shared< FrameVariableBodyPanelGeometrySettings >(
            bodyToTrack, towardsTrackedBody, positionVectorFunction, area, panelTemperature, frameOrientation );
}

inline std::shared_ptr< BodyPanelGeometrySettings > frameFixedPanelGeometry( const Eigen::Vector3d& surfaceNormal,
                                                                             const double area,
                                                                             const std::string& frameOrientation = "" )
{
    return std::make_shared< FrameFixedBodyPanelGeometrySettings >( surfaceNormal, area, frameOrientation );
}

inline std::shared_ptr< BodyPanelGeometrySettings > timeVaryingPanelGeometry(
        const std::function< Eigen::Vector3d( ) >& surfaceNormalFunction,
        const double area,
        const std::string& frameOrientation = "" )
{
    return std::make_shared< FrameVariableBodyPanelGeometrySettings >( surfaceNormalFunction, area, frameOrientation );
}

inline std::shared_ptr< BodyPanelGeometrySettings > bodyTrackingPanelGeometry( const std::string& bodyToTrack,
                                                                               const bool towardsTrackedBody,
                                                                               const double area,
                                                                               const std::string& frameOrientation = "" )
{
    return std::make_shared< FrameVariableBodyPanelGeometrySettings >( bodyToTrack, towardsTrackedBody, area, frameOrientation );
}

class BodyPanelSettings
{
public:
    BodyPanelSettings( const std::shared_ptr< BodyPanelGeometrySettings > panelGeometry,
                       const std::shared_ptr< BodyPanelReflectionLawSettings > reflectionLawSettings,
                       const std::string panelTypeId = "",
                       const std::shared_ptr< MaterialProperties > materialProperties = nullptr ):
        panelGeometry_( panelGeometry ), reflectionLawSettings_( reflectionLawSettings ), panelTypeId_( panelTypeId ),
        materialProperties_( materialProperties )
    { }

    std::shared_ptr< BodyPanelGeometrySettings > panelGeometry_;

    std::shared_ptr< BodyPanelReflectionLawSettings > reflectionLawSettings_;

    std::string panelTypeId_;

    std::shared_ptr< MaterialProperties > materialProperties_;

};

inline std::shared_ptr< BodyPanelSettings > bodyPanelSettings(
        const std::shared_ptr< BodyPanelGeometrySettings > panelGeometry,
        const std::shared_ptr< BodyPanelReflectionLawSettings > reflectionLawSettings,
        const std::string panelTypeId = "",
        const std::shared_ptr< MaterialProperties > materialProperties = nullptr )
{
    return std::make_shared< BodyPanelSettings >( panelGeometry, reflectionLawSettings, panelTypeId, materialProperties );
}

class FullPanelledBodySettings
{
public:
    FullPanelledBodySettings( const std::vector< std::shared_ptr< BodyPanelSettings > >& panelSettingsList,
                              const std::map< std::string, std::shared_ptr< RotationModelSettings > >& partRotationModelSettings =
                                      std::map< std::string, std::shared_ptr< RotationModelSettings > >( ) ):
        panelSettingsList_( panelSettingsList ), partRotationModelSettings_( partRotationModelSettings )
    { }

    void addPanelRotationSettings( const std::shared_ptr< RotationModelSettings > rotationSettings, const std::string partName )
    {
        partRotationModelSettings_[ partName ] = rotationSettings;
    }

    std::vector< std::shared_ptr< BodyPanelSettings > > panelSettingsList_;

    std::map< std::string, std::shared_ptr< RotationModelSettings > > partRotationModelSettings_;
};

inline std::vector< std::shared_ptr< BodyPanelSettings > > bodyPanelSettingsListFromDae(
        const std::string filePath,
        const Eigen::Vector3d frameOrigin,
        std::map< std::string, std::shared_ptr< MaterialProperties > > materialPropertiesMap,
        std::map< std::string, bool > instantaneousReradiation,
        const std::string inputUnit = "m",
        const std::string frameOrientation = "" )
{
    std::vector< std::shared_ptr< BodyPanelSettings > > bodyPanelSettingsList;
    std::ifstream file(filePath);
    if (!file) 
    {
        throw std::runtime_error("Error, wrong path to DAE file!");
    }
    // collada parsing
    std::stringstream buffer;
    buffer << file.rdbuf(); 
    std::string collada = buffer.str();
    // library materials
    std::string tag = "<library_materials>";
    size_t start_materials = collada.find(tag) + tag.size();
    size_t end_materials = collada.find("</library_materials>", start_materials)-1;
    size_t pos = start_materials;
    std::vector<std::string> panelMaterialIdList;
    tag = "name=\"";
    while (pos<end_materials) {
        size_t start_mat = collada.find(tag, pos);
        if (start_mat == std::string::npos || start_mat >= end_materials) {
            break;
        }
        start_mat += tag.size();
        size_t end_mat = collada.find('"', start_mat);
        if (end_mat == std::string::npos || end_mat > end_materials) {
            break;
        }
        panelMaterialIdList.push_back(collada.substr(start_mat, end_mat-start_mat));
        pos = end_mat + 1;
    }
    // error handling
    for (const auto& materialId : panelMaterialIdList) 
    {
        if ( materialPropertiesMap.count( materialId ) == 0 )
        {
            throw std::runtime_error("Material ID " + materialId + " not found in material properties settings!");
        }
        if ( materialPropertiesMap.at( materialId )->specularReflectivity_ == -1 ||
             materialPropertiesMap.at( materialId )->diffuseReflectivity_ == -1 )
        {
             throw std::runtime_error("Material ID " + materialId + " has no specular or diffuse reflectivity coefficients!");
        }
        auto reradIt = instantaneousReradiation.find(materialId);
        if (reradIt == instantaneousReradiation.end())
        {
            throw std::runtime_error("Material ID " + materialId + " not found in re-radiation settings!");
        }
    }
    // library geometries
    tag = "<library_geometries>";
    size_t start_geometries = collada.find(tag) + tag.size();
    tag = "<source";
    size_t start_source = collada.find(tag, start_geometries) + tag.size();
    tag = "<float_array";
    size_t start_float = collada.find(tag, start_source) + tag.size();
    size_t start_array = collada.find(">", start_float) + 1;
    size_t end_array = collada.find("<", start_array);
    std::istringstream dummy1(collada.substr(start_array, end_array-start_array));
    std::vector<Eigen::Vector3d> panelVerticesList;

    float conversionFactor;

    if (inputUnit=="mm"){conversionFactor=1/1000.;}
    else if (inputUnit=="m"){conversionFactor=1;}
    else if (inputUnit=="in"){conversionFactor=0.0254;}
    else{ throw std::runtime_error("Input Unit " + inputUnit + " not recognized. Valid options are mm, m, in.");}


    double x1, y1, z1;
    while (dummy1 >> x1 >> y1 >> z1) {
        panelVerticesList.push_back(Eigen::Vector3d(x1, y1, z1)*conversionFactor);
    }
    pos = end_array;
    // sub library: create settings
    tag = "<triangles";
    for (size_t i = 0; i<panelMaterialIdList.size(); i++) {
        size_t start_triangle = collada.find(tag, pos) + tag.size();
        std::string tag_triangle = "<p>";
        size_t start_p = collada.find("<p>", start_triangle) + tag_triangle.size();
        size_t end_p = collada.find("<", start_p);
        std::istringstream dummy2(collada.substr(start_p, end_p-start_p));
        int a, b, c, n;
        while (dummy2 >> a >> n >> b >> n >> c >> n) {
            Eigen::Vector3d edgeAB = panelVerticesList[b] - panelVerticesList[a];
            Eigen::Vector3d edgeAC = panelVerticesList[c] - panelVerticesList[a];
            double area = 0.5*edgeAB.cross(edgeAC).norm();
            Eigen::Vector3d surfaceNormal = edgeAB.cross(edgeAC).normalized();
            std::shared_ptr< FrameFixedBodyPanelGeometrySettings > currentGeometrySettings = 
                std::make_shared< FrameFixedBodyPanelGeometrySettings >(panelVerticesList[a],
                     panelVerticesList[b], panelVerticesList[c], surfaceNormal, area, frameOrigin, frameOrientation);
            bodyPanelSettingsList.push_back(bodyPanelSettings(currentGeometrySettings, 
                specularDiffuseBodyPanelReflectionLawSettings(materialPropertiesMap[panelMaterialIdList[i]]->specularReflectivity_, 
                                                              materialPropertiesMap[panelMaterialIdList[i]]->diffuseReflectivity_,
                                                              instantaneousReradiation[panelMaterialIdList[i]]),
                                                              panelMaterialIdList[i],
                                                              materialPropertiesMap[ panelMaterialIdList[i] ]) );
                                                            
        }
        pos = end_p;
    }

    return bodyPanelSettingsList;

}

inline std::vector< std::shared_ptr< BodyPanelSettings > > mergeBodyPanelSettingsLists( std::vector< std::vector< std::shared_ptr< BodyPanelSettings > > > listOfLists )
{
    std::vector< std::shared_ptr<BodyPanelSettings > > mergedList;
    for ( auto it : listOfLists ) 
    {   
        mergedList.insert(mergedList.end(), it.begin(), it.end());
    }
    return mergedList;
    
}

inline std::shared_ptr< FullPanelledBodySettings > fullPanelledBodySettings(
        const std::vector< std::shared_ptr< BodyPanelSettings > >& panelSettingsList,
        const std::map< std::string, std::shared_ptr< RotationModelSettings > >& partRotationModelSettings =
                std::map< std::string, std::shared_ptr< RotationModelSettings > >( ) )
{
    return std::make_shared< FullPanelledBodySettings >( panelSettingsList, partRotationModelSettings );
}

inline std::shared_ptr< FullPanelledBodySettings > bodyWingPanelledGeometry( const double length,
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
                    { frameFixedPanelGeometry( Eigen::Vector3d::UnitX( ), Eigen::Vector3d::UnitX( ) * height / 2, width * height ),
                      frameFixedPanelGeometry( -Eigen::Vector3d::UnitX( ), -Eigen::Vector3d::UnitX( ) * height / 2, width * height ),
                      frameFixedPanelGeometry( Eigen::Vector3d::UnitY( ), Eigen::Vector3d::UnitY( ) * width / 2, length * height ),
                      frameFixedPanelGeometry( -Eigen::Vector3d::UnitY( ), -Eigen::Vector3d::UnitY( ) * width / 2, length * height ),
                      frameFixedPanelGeometry( Eigen::Vector3d::UnitZ( ), Eigen::Vector3d::UnitZ( ) * length / 2, length * width ),
                      frameFixedPanelGeometry( -Eigen::Vector3d::UnitZ( ), Eigen::Vector3d::UnitZ( ) * length / 2, length * width ),
                      bodyTrackingPanelGeometry(
                              "Sun", true, [ =, this ]( ) { return -Eigen::Vector3d::UnitY( ) * width / 2; }, totalSolarArrayArea ),
                      bodyTrackingPanelGeometry(
                              "Sun", false, [ =, this ]( ) { return Eigen::Vector3d::UnitY( ) * width / 2; }, totalSolarArrayArea ) } );
    std::vector< std::shared_ptr< BodyPanelSettings > > panelSettings;
    for( unsigned int i = 0; i < 6; i++ )
    {
        panelSettings.push_back( bodyPanelSettings( panelGeometrySettingsList.at( i ),
                                                    specularDiffuseBodyPanelReflectionLawSettings(
                                                            busSpecularReflectivity, busDiffuseReflectivity, busInstantaneousReradiation ),
                                                    "Bus" ) );
    }

    for( unsigned int i = 0; i < 2; i++ )
    {
        panelSettings.push_back( bodyPanelSettings(
                panelGeometrySettingsList.at( i + 6 ),
                specularDiffuseBodyPanelReflectionLawSettings(
                        solarArraySpecularReflectivity, solarArrayDiffuseReflectivity, solarArrayInstantaneousReradiation ),
                "SolarPanel" ) );
    }

    return fullPanelledBodySettings( panelSettings );
}

void addEngineModel( const std::string& bodyName,
                     const std::string& engineName,
                     const std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustSettings,
                     const simulation_setup::SystemOfBodies& bodies,
                     const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) );

void addVariableDirectionEngineModel( const std::string& bodyName,
                                      const std::string& engineName,
                                      const std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustSettings,
                                      const simulation_setup::SystemOfBodies& bodies,
                                      const std::function< Eigen::Vector3d( const double ) > bodyFixedThrustDirection );

std::pair< std::shared_ptr< system_models::VehicleExteriorPanel >, std::string > createBodyExteriorPanel(
        const std::shared_ptr< BodyPanelSettings > panelSettings,
        const std::string& bodyName,
        const simulation_setup::SystemOfBodies& bodies );

void addBodyExteriorPanelledShape( const std::shared_ptr< FullPanelledBodySettings > panelSettings,
                                   const std::string& bodyName,
                                   const simulation_setup::SystemOfBodies& bodies );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATESSYTEMMODEL_H
