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

#include <memory>

#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/body.h"

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
    paneled_target
};

/*!
 * Settings for a radiation pressure target model.
 *
 * @see RadiationPressureTargetModel
 */
class RadiationPressureTargetModelSettings
{
public:
    explicit RadiationPressureTargetModelSettings(
            RadiationPressureTargetModelType radiationPressureTargetModelType) :
            radiationPressureTargetModelType_(radiationPressureTargetModelType) {}

    virtual ~RadiationPressureTargetModelSettings() = default;

    RadiationPressureTargetModelType getRadiationPressureTargetModelType() const
    {
        return radiationPressureTargetModelType_;
    }

private:
    RadiationPressureTargetModelType radiationPressureTargetModelType_;
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
     */
    explicit CannonballRadiationPressureTargetModelSettings(
            double area,
            double coefficient) :
            RadiationPressureTargetModelSettings(cannonball_target),
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

/*!
 * Settings for a paneled radiation pressure target model.
 *
 * @see PaneledRadiationPressureTargetModel
 */
class PaneledRadiationPressureTargetModelSettings : public RadiationPressureTargetModelSettings
{
public:
    class Panel;

    /*!
     * Constructor.
     *
     * @param panels Vector of settings for panels comprising the paneled target
     */
    explicit PaneledRadiationPressureTargetModelSettings(
            const std::vector<Panel>& panels) :
            RadiationPressureTargetModelSettings(paneled_target),
            panels_(panels) {}

    const std::vector<Panel>& getPanels() const
    {
        return panels_;
    }

private:
    std::vector<Panel> panels_;
};

/*!
 * Settings for a single panel of a paneled radiation pressure target model.
 *
 * @see PaneledRadiationPressureTargetModel::Panel
 */
class PaneledRadiationPressureTargetModelSettings::Panel
{
public:
    /*!
     * Constructor for a panel with a given surface normal vector.
     *
     * @param area Area of the panel [m²]
     * @param specularReflectivity Specular reflectivity coefficient of the panel [-]
     * @param diffuseReflectivity Diffuse reflectivity coefficient of the panel [-]
     * @param surfaceNormalFunction Function returning the surface normal vector of the panel [-]
     */
    explicit Panel(double area,
                   double specularReflectivity,
                   double diffuseReflectivity,
                   const std::function<Eigen::Vector3d()>& surfaceNormalFunction) :
            area_(area),
            specularReflectivity_(specularReflectivity),
            diffuseReflectivity_(diffuseReflectivity),
            surfaceNormalFunction_(surfaceNormalFunction)  {}

    /*!
    * Constructor for a panel with a given surface normal vector.
     *
    * @param area Area of the panel [m²]
    * @param specularReflectivity Specular reflectivity coefficient of the panel [-]
    * @param diffuseReflectivity Diffuse reflectivity coefficient of the panel [-]
    * @param surfaceNormal Surface normal vector of the panel [-]
    */
    explicit Panel(double area,
                   double specularReflectivity,
                   double diffuseReflectivity,
                   const Eigen::Vector3d& surfaceNormal) :
            Panel(area,
                  specularReflectivity,
                  diffuseReflectivity,
                  [=] () { return surfaceNormal; }) {}

    /*!
    * Constructor for a panel tracking a body (e.g., a solar panel tracking the Sun, or an antenna tracking Earth).
     *
    * @param area Area of the panel [m²]
    * @param specularReflectivity Specular reflectivity coefficient of the panel [-]
    * @param diffuseReflectivity Diffuse reflectivity coefficient of the panel [-]
    * @param bodyToTrack Name of the body to track
    * @param towardsTrackedBody Whether the panel points towards or away from tracked body
    */
    explicit Panel(double area,
                   double specularReflectivity,
                   double diffuseReflectivity,
                   const std::string& bodyToTrack,
                   const bool towardsTrackedBody = true) :
            area_(area),
            specularReflectivity_(specularReflectivity),
            diffuseReflectivity_(diffuseReflectivity),
            bodyToTrack_(bodyToTrack),
            towardsTrackedBody_(towardsTrackedBody) {}

    double getArea() const
    {
        return area_;
    }

    const std::function<Eigen::Vector3d()>& getSurfaceNormalFunction() const
    {
        return surfaceNormalFunction_;
    }

    double getSpecularReflectivity() const
    {
        return specularReflectivity_;
    }

    double getDiffuseReflectivity() const
    {
        return diffuseReflectivity_;
    }

    const std::string &getBodyToTrack() const
    {
        return bodyToTrack_;
    }

    bool isTowardsTrackedBody() const
    {
        return towardsTrackedBody_;
    }

private:
    double area_;
    double specularReflectivity_;
    double diffuseReflectivity_;
    std::function<Eigen::Vector3d()> surfaceNormalFunction_;
    std::string bodyToTrack_;
    bool towardsTrackedBody_{true};
};

/*!
 * Create settings for a cannonball radiation pressure target model.
 *
 * @param area Cross-sectional area of cannonball [m²]
 * @param coefficient Radiation pressure coefficient [-]
 * @return Shared pointer to settings for a cannonball radiation pressure target model
 */
inline std::shared_ptr<CannonballRadiationPressureTargetModelSettings>
        cannonballRadiationPressureTargetModelSettings(double area, double coefficient)
{
    return std::make_shared< CannonballRadiationPressureTargetModelSettings >(area, coefficient);
}

/*!
 * Create settings for a paneled radiation pressure target model.
 *
 * @param panels List of settings for panels comprising the paneled target
 * @return Shared pointer to settings for a paneled radiation pressure target model
 */
inline std::shared_ptr<PaneledRadiationPressureTargetModelSettings>
        paneledRadiationPressureTargetModelSettings(std::initializer_list<PaneledRadiationPressureTargetModelSettings::Panel> panels)
{
    return std::make_shared< PaneledRadiationPressureTargetModelSettings >(
            std::vector<PaneledRadiationPressureTargetModelSettings::Panel>(panels));
}

/*!
 * Create radiation pressure target model from its settings.
 *
 * @param modelSettings Settings of the radiation pressure target model
 * @param body Body to which the radiation pressure target model belongs
 * @param bodies System of bodies
 * @return Shared pointer to radiation pressure target model
 */
std::shared_ptr<electromagnetism::RadiationPressureTargetModel> createRadiationPressureTargetModel(
        const std::shared_ptr< RadiationPressureTargetModelSettings >& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies);

typedef PaneledRadiationPressureTargetModelSettings::Panel TargetPanelSettings;

} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H
