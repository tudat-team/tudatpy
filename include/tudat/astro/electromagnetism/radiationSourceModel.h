/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      E. B. Saff, et al. "Distributing many points on a sphere".
 *          The Mathematical Intelligencer 19. 1(1997): 5–11.
 *      Frank G. Lemoine, et al. "High‒degree gravity models from GRAIL primary mission data".
 *          Journal of Geophysical Research: Planets 118. 8(2013): 1676–1698.
 *      Knocke, Philip et al. "Earth radiation pressure effects on satellites."
 *          Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988.
 */

#ifndef TUDAT_RADIATIONSOURCEMODEL_H
#define TUDAT_RADIATIONSOURCEMODEL_H

#include <vector>
#include <map>
#include <utility>
#include <functional>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/electromagnetism/luminosityModel.h"
#include "tudat/astro/electromagnetism/reflectionLaw.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"

namespace tudat
{
namespace electromagnetism
{

typedef std::vector<std::pair<double, Eigen::Vector3d>> IrradianceWithSourceList;


class RadiationSourceModel
{
public:
    explicit RadiationSourceModel() = default;

    virtual ~RadiationSourceModel() = default;

    void updateMembers(double currentTime );

    /*!
     * Evaluate the irradiance [W/m²] at a certain position due to this source.
     * @param targetPosition Position where to evaluate the irradiance in local (i.e. source-fixed) coordinates
     * @param originalSourceIrradiance Irradiance from the original source (if applicable)
     * @param originalSourceToSourceDirection Direction of incoming radiation in local (i.e. source-fixed) coordinates
     * @return A list of irradiances and their source-fixed origin. Every element can be thought of as a ray.
     * Single element for point sources, multiple elements for paneled sources.
     */
    virtual IrradianceWithSourceList evaluateIrradianceAtPosition(
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const = 0;

protected:
    virtual void updateMembers_(const double currentTime) {};

    double currentTime_{TUDAT_NAN};
};

//*********************************************************************************************
//   Isotropic point radiation source
//*********************************************************************************************

class IsotropicPointRadiationSourceModel : public RadiationSourceModel
{
public:
    explicit IsotropicPointRadiationSourceModel(
            const std::shared_ptr<LuminosityModel>& luminosityModel) :
        luminosityModel_(luminosityModel) {}

    IrradianceWithSourceList evaluateIrradianceAtPosition(
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    double evaluateIrradianceAtPosition(const Eigen::Vector3d& targetPosition) const;

    const std::shared_ptr<LuminosityModel>& getLuminosityModel() const
    {
        return luminosityModel_;
    }

private:
    void updateMembers_(double currentTime) override;

    std::shared_ptr<LuminosityModel> luminosityModel_;
};

//*********************************************************************************************
//   Paneled radiation source
//*********************************************************************************************

/*!
 * A source that is discretized into panels, each with its own radiosity model. A radiosity model describes the
 * radiation emitted (thermal) or reflected (albedo) by each panel.
 */
class PaneledRadiationSourceModel : public RadiationSourceModel
{
public:
    class Panel;
    class PanelRadiosityModel; // cannot make this nested class of Panel since forward declaration is impossible

    explicit PaneledRadiationSourceModel(
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& sourceBodyShapeModel) :
        sourceBodyShapeModel_(sourceBodyShapeModel) {}

    explicit PaneledRadiationSourceModel() :
            PaneledRadiationSourceModel(nullptr) {}

    IrradianceWithSourceList evaluateIrradianceAtPosition(
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    virtual const std::vector<Panel>& getPanels() const = 0;

protected:
    std::shared_ptr<basic_astrodynamics::BodyShapeModel> sourceBodyShapeModel_;
};

/*!
 * A paneled source with paneling of the whole body that is constant and generated only once. Panel properties such as
 * albedo and emissivity are constant in time as well.
 */
class StaticallyPaneledRadiationSourceModel : public PaneledRadiationSourceModel
{
public:
    explicit StaticallyPaneledRadiationSourceModel(
            const std::vector<Panel>& panels) :
            PaneledRadiationSourceModel(),
            numberOfPanels(panels.size()),
            panels_(panels) {}

    explicit StaticallyPaneledRadiationSourceModel(
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& sourceBodyShapeModel,
            const std::vector<std::function<std::shared_ptr<PanelRadiosityModel>(double, double)>>& radiosityModelFunctions,
            int numberOfPanels) :
            PaneledRadiationSourceModel(sourceBodyShapeModel),
            numberOfPanels(numberOfPanels),
            radiosityModelFunctions_(radiosityModelFunctions) {}

    const std::vector<Panel>& getPanels() const override
    {
        return panels_;
    }

private:
    void updateMembers_(double currentTime) override;

    unsigned int numberOfPanels;
    std::vector<std::function<std::shared_ptr<PanelRadiosityModel>(double, double)>> radiosityModelFunctions_;
    std::vector<Panel> panels_;
};

//class DynamicallyPaneledRadiationSourceModel : public PaneledRadiationSourceModel
//{
//private:
//    // keep target-specific panel lists, otherwise may have to regenerate a lot
//    std::map<std::string, std::vector<Panel>> panels_;
//};

class PaneledRadiationSourceModel::Panel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    explicit Panel(
            double area,
            const Eigen::Vector3d& relativeCenter,
            const Eigen::Vector3d& surfaceNormal,
            const std::vector<std::shared_ptr<PanelRadiosityModel>>& radiosityModels) :
        relativeCenter_(relativeCenter),
        surfaceNormal_(surfaceNormal),
        radiosityModels_(radiosityModels),
        area_(area) {}

    double getArea() const
    {
        return area_;
    }

    const Eigen::Vector3d& getRelativeCenter() const
    {
        return relativeCenter_;
    }

    const Eigen::Vector3d& getSurfaceNormal() const
    {
        return surfaceNormal_;
    }

    const std::vector<std::shared_ptr<PanelRadiosityModel>>& getRadiosityModels() const
    {
        return radiosityModels_;
    }

private:
    Eigen::Vector3d relativeCenter_;
    Eigen::Vector3d surfaceNormal_;
    std::vector<std::shared_ptr<PanelRadiosityModel>> radiosityModels_;
    double area_;
};

class PaneledRadiationSourceModel::PanelRadiosityModel
{
public:
    virtual ~PanelRadiosityModel() = default;

    /*!
     * Evaluate the irradiance [W/m²] at a certain position due to this panel.
     * @param panel The panel this radiosity model belongs to
     * @param targetPosition Position where to evaluate the irradiance in local coordinates (source rotation, centered in panel)
     * @return The irradiance
     */
    virtual double evaluateIrradianceAtPosition(
            const Panel& panel,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const = 0;
};

/*!
 * Panel radiosity model for albedo radiation with an arbitrary reflection law. This model was introduced in
 * Knocke (1988) for Earth thermal radiation, but assuming Lambertian reflectance.
 *
 * For most cases, albedo radiation with a diffuse-only Lambertian reflectance law is sufficient.
 */
class AlbedoPanelRadiosityModel : public PaneledRadiationSourceModel::PanelRadiosityModel
{
public:
    explicit AlbedoPanelRadiosityModel(
            const std::shared_ptr<ReflectionLaw>& reflectionLaw) :
            reflectionLaw_(reflectionLaw) {}

    double evaluateIrradianceAtPosition(
            const PaneledRadiationSourceModel::Panel& panel,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    const std::shared_ptr<ReflectionLaw>& getReflectionLaw() const
    {
        return reflectionLaw_;
    }

private:
    std::shared_ptr<ReflectionLaw> reflectionLaw_;
};

/*!
 * Panel radiosity model for thermal emissions, based on delayed, isotropic and constant flux. This model was introduced
 * in Knocke (1988) for Earth thermal radiation.
 */
class DelayedThermalPanelRadiosityModel : public PaneledRadiationSourceModel::PanelRadiosityModel
{
public:
    explicit DelayedThermalPanelRadiosityModel(
            double emissivity) :
            emissivity_(emissivity) {}

    double evaluateIrradianceAtPosition(
            const PaneledRadiationSourceModel::Panel& panel,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    double getEmissivity() const
    {
        return emissivity_;
    }

private:
    double emissivity_;
};

/*!
 * Panel radiosity model for thermal emissions, based on angle to subsolar point. This model was introduced in
 * Lemoine (2013) for lunar thermal radiation. At the subsolar point, the thermal radiation is due to
 * emissivity-corrected black-body radiation at maxTemperature. At positions away from the sub-solar point, the
 * temperature drops, until it reaches minTemperature on the backside.
 */
class AngleBasedThermalPanelRadiosityModel : public PaneledRadiationSourceModel::PanelRadiosityModel
{
public:
    explicit AngleBasedThermalPanelRadiosityModel(
            double minTemperature,
            double maxTemperature,
            double emissivity) :
            minTemperature_(minTemperature),
            maxTemperature_(maxTemperature),
            emissivity_(emissivity) {}

    double evaluateIrradianceAtPosition(
            const PaneledRadiationSourceModel::Panel& panel,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    double getMinTemperature() const
    {
        return minTemperature_;
    }

    double getMaxTemperature() const
    {
        return maxTemperature_;
    }

    double getEmissivity() const
    {
        return emissivity_;
    }

private:
    double minTemperature_;
    double maxTemperature_;
    double emissivity_;
};

/*!
 * Generate evenly spaced points on a sphere using algorithm from Saff (1997).
 * @param n number of points to generate
 * @return a pair of vectors, first vector are polar angles, second vector are azimuth angles
 */
std::pair<std::vector<double>, std::vector<double>> generateEvenlySpacedPoints(unsigned int n);


} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONSOURCEMODEL_H
