/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H
#define TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H

#include <map>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/electromagnetism/reflectionLaw.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/body.h"
//#include "tudat/simulation/environment_setup/createSystemModel.h"

namespace tudat
{
namespace simulation_setup
{

/*!
 * Types of radiation pressure target models.
 */
enum RadiationPressureTargetModelType
{
    cannonball_target,
    paneled_target,
    multi_type_target,
    undefined_target
};

std::shared_ptr<electromagnetism::RadiationPressureTargetModel> getRadiationPressureTargetModelOfType(
    const std::shared_ptr< Body > target,
    const RadiationPressureTargetModelType targetModelType = undefined_target,
    const std::string errorOutput = "" );

/*!
 * Settings for a radiation pressure target model.
 *
 * @see RadiationPressureTargetModel
 */
class RadiationPressureTargetModelSettings
{
public:
    explicit RadiationPressureTargetModelSettings(
            RadiationPressureTargetModelType radiationPressureTargetModelType,
            const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies = {}) :
            radiationPressureTargetModelType_(radiationPressureTargetModelType),
            sourceToTargetOccultingBodies_(sourceToTargetOccultingBodies) {}

    virtual ~RadiationPressureTargetModelSettings() = default;

    RadiationPressureTargetModelType getRadiationPressureTargetModelType() const
    {
        return radiationPressureTargetModelType_;
    }

    std::map<std::string, std::vector<std::string>> getSourceToTargetOccultingBodies() const
    {
        return sourceToTargetOccultingBodies_;
    }

private:
    RadiationPressureTargetModelType radiationPressureTargetModelType_;
    // Map (source name -> list of occulting body names) of bodies to occult sources as seen from this target
    // If the same occulting bodies are to be used for all sources, there will be a single entry
    // with an emptry string as key
    std::map<std::string, std::vector<std::string>> sourceToTargetOccultingBodies_;
};

class MultiRadiationPressureTargetModelSettings: public RadiationPressureTargetModelSettings
{
public:
    explicit MultiRadiationPressureTargetModelSettings(
        std::vector< std::shared_ptr< RadiationPressureTargetModelSettings > > radiationPressureTargetModelSettings ) :
        RadiationPressureTargetModelSettings( multi_type_target ),
        radiationPressureTargetModelSettings_( radiationPressureTargetModelSettings ){ }

    std::vector< std::shared_ptr< RadiationPressureTargetModelSettings > > radiationPressureTargetModelSettings_;
};
/*!
 * Settings for a cannonball radiation pressure target model.
 *
 * @see CannonballRadiationPressureTargetModel
 */
class CannonballRadiationPressureTargetModelSettings : public RadiationPressureTargetModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param area Cross-sectional area of cannonball [m²]
     * @param coefficient Radiation pressure coefficient [-]
     * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
     *      to occult sources as seen from this target
     */
    explicit CannonballRadiationPressureTargetModelSettings(
            double area, double coefficient,
            const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies = {}) :
            RadiationPressureTargetModelSettings(
                    RadiationPressureTargetModelType::cannonball_target, sourceToTargetOccultingBodies),
            area_(area),
            coefficient_(coefficient) {}

    double getArea() const
    {
        return area_;
    }

    double getCoefficient() const
    {
        return coefficient_;
    }

private:
    double area_;

    double coefficient_;
};
//
///*!
// * Settings for a paneled radiation pressure target model.
// *
// * @see PaneledRadiationPressureTargetModel
// */
//class PaneledRadiationPressureTargetModelSettings : public RadiationPressureTargetModelSettings
//{
//public:
////    class Panel;
//
//    /*!
//     * Constructor.
//     *
//     * @param panels Vector of settings for panels comprising the paneled target
//     * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
//     *      to occult sources as seen from this target
//     */
//    explicit PaneledRadiationPressureTargetModelSettings(
//            const std::vector<Panel>& panels,
//            const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies) :
//            RadiationPressureTargetModelSettings(
//                    RadiationPressureTargetModelType::paneled_target, sourceToTargetOccultingBodies),
//            panels_(panels) {}
//
//    const std::vector<Panel>& getPanels() const
//    {
//        return panels_;
//    }
//
//private:
//    std::vector<Panel> panels_;
//};

enum BodyPanelReflectionLawType
{
    specular_diffuse_reflection_law
};

class BodyPanelReflectionLawSettings
{
public:
    BodyPanelReflectionLawSettings(
        const BodyPanelReflectionLawType bodyPanelReflectionLawType ):
        bodyPanelReflectionLawType_( bodyPanelReflectionLawType ){ }

    virtual ~BodyPanelReflectionLawSettings( ){ }

    BodyPanelReflectionLawType bodyPanelReflectionLawType_;
};

class SpecularDiffuseBodyPanelReflectionLawSettings: public BodyPanelReflectionLawSettings
{
public:

    SpecularDiffuseBodyPanelReflectionLawSettings(
        const double specularReflectivity,
        const double diffuseReflectivity,
        const bool withInstantaneousReradiation ):
            BodyPanelReflectionLawSettings( specular_diffuse_reflection_law ),
            specularReflectivity_(specularReflectivity),
            diffuseReflectivity_(diffuseReflectivity),
            withInstantaneousReradiation_(withInstantaneousReradiation)
    {
        absorptivity_ = 1.0 - ( specularReflectivity_ + diffuseReflectivity_ );
    }

    ~SpecularDiffuseBodyPanelReflectionLawSettings( ){ }

    double specularReflectivity_;

    double diffuseReflectivity_;

    double absorptivity_;

    bool withInstantaneousReradiation_;
};

inline std::shared_ptr< BodyPanelReflectionLawSettings > specularDiffuseBodyPanelReflectionLawSettings(
    const double specularReflectivity,
    const double diffuseReflectivity,
    const bool withInstantaneousReradiation = true )
{
    return std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings>(
        specularReflectivity, diffuseReflectivity, withInstantaneousReradiation );
}

inline std::shared_ptr< BodyPanelReflectionLawSettings > lambertainBodyPanelReflectionLawSettings(
    const double reflectivity )
{
    return std::make_shared< SpecularDiffuseBodyPanelReflectionLawSettings>(
        0.0, reflectivity, false );
}

std::shared_ptr< electromagnetism::ReflectionLaw > createReflectionLaw(
    const std::shared_ptr<BodyPanelReflectionLawSettings> modelSettings );

//
///*!
// * Settings for a single panel of a paneled radiation pressure target model.
// *
// * @see PaneledRadiationPressureTargetModel::Panel
// */
//class PaneledRadiationPressureTargetModelSettings::Panel
//{
//public:
//    /*!
//     * Constructor for a panel with a given surface normal vector.
//     *
//     * @param area Area of the panel [m²]
//     * @param specularReflectivity Specular reflectivity coefficient of the panel [-]
//     * @param diffuseReflectivity Diffuse reflectivity coefficient of the panel [-]
//     * @param withInstantaneousReradiation Whether to instantaneously reradiate absorbed radiation
//     * @param surfaceNormalFunction Function returning the surface normal vector of the panel [-]
//     */
//    explicit Panel(double area,
//                   double specularReflectivity,
//                   double diffuseReflectivity,
//                   bool withInstantaneousReradiation,
//                   const std::function<Eigen::Vector3d()>& surfaceNormalFunction) :
//            area_(area),
//            specularReflectivity_(specularReflectivity),
//            diffuseReflectivity_(diffuseReflectivity),
//            withInstantaneousReradiation_(withInstantaneousReradiation),
//            surfaceNormalFunction_(surfaceNormalFunction)  {}
//
//    /*!
//    * Constructor for a panel with a given surface normal vector.
//    *
//    * @param area Area of the panel [m²]
//    * @param specularReflectivity Specular reflectivity coefficient of the panel [-]
//    * @param diffuseReflectivity Diffuse reflectivity coefficient of the panel [-]
//    * @param withInstantaneousReradiation Whether to instantaneously reradiate absorbed radiation
//    * @param surfaceNormal Surface normal vector of the panel [-]
//    */
//    explicit Panel(double area,
//                   double specularReflectivity,
//                   double diffuseReflectivity,
//                   bool withInstantaneousReradiation,
//                   const Eigen::Vector3d& surfaceNormal) :
//            Panel(area,
//                  specularReflectivity,
//                  diffuseReflectivity,
//                  withInstantaneousReradiation,
//                  [=] () { return surfaceNormal; }) {}
//
//    /*!
//    * Constructor for a panel tracking a body (e.g., a solar panel tracking the Sun, or an antenna tracking Earth).
//     *
//    * @param area Area of the panel [m²]
//    * @param specularReflectivity Specular reflectivity coefficient of the panel [-]
//    * @param diffuseReflectivity Diffuse reflectivity coefficient of the panel [-]
//    * @param withInstantaneousReradiation Whether to instantaneously reradiate absorbed radiation
//    * @param bodyToTrack Name of the body to track
//    * @param towardsTrackedBody Whether the panel points towards or away from tracked body
//    */
//    explicit Panel(double area,
//                   double specularReflectivity,
//                   double diffuseReflectivity,
//                   bool withInstantaneousReradiation,
//                   const std::string& bodyToTrack,
//                   const bool towardsTrackedBody = true) :
//            area_(area),
//            specularReflectivity_(specularReflectivity),
//            diffuseReflectivity_(diffuseReflectivity),
//            withInstantaneousReradiation_(withInstantaneousReradiation),
//            bodyToTrack_(bodyToTrack),
//            towardsTrackedBody_(towardsTrackedBody) {}
//
//    double getArea() const
//    {
//        return area_;
//    }
//
//    const std::function<Eigen::Vector3d()>& getSurfaceNormalFunction() const
//    {
//        return surfaceNormalFunction_;
//    }
//
//    double getSpecularReflectivity() const
//    {
//        return specularReflectivity_;
//    }
//
//    double getDiffuseReflectivity() const
//    {
//        return diffuseReflectivity_;
//    }
//
//    bool isWithInstantaneousReradiation() const
//    {
//        return withInstantaneousReradiation_;
//    }
//
//    const std::string &getBodyToTrack() const
//    {
//        return bodyToTrack_;
//    }
//
//    bool isTowardsTrackedBody() const
//    {
//        return towardsTrackedBody_;
//    }
//
//private:
//    double area_;
//    double specularReflectivity_;
//    double diffuseReflectivity_;
//    bool withInstantaneousReradiation_;
//    std::function<Eigen::Vector3d()> surfaceNormalFunction_;
//    std::string bodyToTrack_;
//    bool towardsTrackedBody_{true};
//};

//typedef PaneledRadiationPressureTargetModelSettings::Panel TargetPanelSettings;

/*!
 * Create settings for a cannonball radiation pressure target model. Each source can have its own set
 * of occulting bodies.
 *
 * @param area Cross-sectional area of cannonball [m²]
 * @param coefficient Radiation pressure coefficient [-]
 * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
 *      to occult sources as seen from this target
 * @return Shared pointer to settings for a cannonball radiation pressure target model
 */
inline std::shared_ptr<RadiationPressureTargetModelSettings>
        cannonballRadiationPressureTargetModelSettingsWithOccultationMap(
            double area, double coefficient,
            const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies)
{
    return std::make_shared<CannonballRadiationPressureTargetModelSettings>(
        area, coefficient, sourceToTargetOccultingBodies);
}

/*!
 * Create settings for a cannonball radiation pressure target model. All sources are occulted by the
 * same set of bodies.
 *
 * @param area Cross-sectional area of cannonball [m²]
 * @param coefficient Radiation pressure coefficient [-]
 * @param sourceToTargetOccultingBodies Names of bodies to occult the source as seen from this target
 * @return Shared pointer to settings for a cannonball radiation pressure target model
 */
inline std::shared_ptr<RadiationPressureTargetModelSettings>
        cannonballRadiationPressureTargetModelSettings(
            double area, double coefficient,
            const std::vector<std::string>& sourceToTargetOccultingBodies = {})
{
    const std::map<std::string, std::vector<std::string>> occultingBodiesMap {{"", sourceToTargetOccultingBodies}};
    return cannonballRadiationPressureTargetModelSettingsWithOccultationMap(
        area, coefficient, occultingBodiesMap);
}

/*!
 * Create settings for a paneled radiation pressure target model. Each source can have its own set
 * of occulting bodies.
 *
 * @param panels List of settings for panels comprising the paneled target
 * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
 *      to occult sources as seen from this target
 * @return Shared pointer to settings for a paneled radiation pressure target model
 */
inline std::shared_ptr< RadiationPressureTargetModelSettings >
        paneledRadiationPressureTargetModelSettingsWithOccultationMap(
//            std::initializer_list<PaneledRadiationPressureTargetModelSettings::Panel> panels,
            const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies)
{
    return std::make_shared<RadiationPressureTargetModelSettings>(
            paneled_target,
            sourceToTargetOccultingBodies);
}

/*!
 * Create settings for a paneled radiation pressure target model. All sources are occulted by the
 * same set of bodies.
 *
 * @param panels List of settings for panels comprising the paneled target
 * @param sourceToTargetOccultingBodies Names of bodies to occult the source as seen from the target
 * @param originalSourceToSourceOccultingBodies Names of bodies to occult the source as seen from this target
 * @return Shared pointer to settings for a paneled radiation pressure target model
 */
inline std::shared_ptr<RadiationPressureTargetModelSettings>
        paneledRadiationPressureTargetModelSettings(
//            std::initializer_list<PaneledRadiationPressureTargetModelSettings::Panel> panels,
            const std::vector<std::string>& sourceToTargetOccultingBodies = {})
{
    const std::map<std::string, std::vector<std::string>> occultingBodiesMap {{"", sourceToTargetOccultingBodies}};
    return std::make_shared<RadiationPressureTargetModelSettings>(
        paneled_target,
        occultingBodiesMap);
}

RadiationPressureTargetModelType getTargetModelType( const std::shared_ptr<electromagnetism::RadiationPressureTargetModel> targetModel );

/*!
 * Create radiation pressure target model from its settings.
 *
 * @param modelSettings Settings of the radiation pressure target model
 * @param body Body to which the radiation pressure target model belongs
 * @param bodies System of bodies
 * @return Shared pointer to radiation pressure target model
 */
std::vector< std::shared_ptr<electromagnetism::RadiationPressureTargetModel> > createRadiationPressureTargetModel(
        const std::shared_ptr< RadiationPressureTargetModelSettings >& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies);

} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H
