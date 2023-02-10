/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SOURCESourcePanelRadiosityModel_H
#define TUDAT_SOURCESourcePanelRadiosityModel_H

#include <memory>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/electromagnetism/reflectionLaw.h"
#include "tudat/astro/electromagnetism/surfacePropertyDistribution.h"


namespace tudat
{
namespace electromagnetism
{

/*!
 * Class modeling the radiosity emitted by a panel of a paneled radiation source.
 *
 * @see PaneledRadiationSourceModel::Panel
 */
class SourcePanelRadiosityModel
{
public:
    virtual ~SourcePanelRadiosityModel() = default;

    /*!
     * Evaluate the irradiance [W/m²] at a certain position due to this panel.
     *
     * @param targetPosition Position where to evaluate the irradiance in panel-local coordinates (source rotation,
     *        centered in panel)
     * @return Irradiance due to this radiosity model for single panel
     */
    virtual double evaluateIrradianceAtPosition(
            double panelArea,
            const Eigen::Vector3d& panelSurfaceNormal,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const = 0;

    /*!
     * Update class members.
     *
     * @param panelLatitude Latitude of the panel this radiosity model belongs to
     * @param panelLongitude Longitude of the panel this radiosity model belongs to
     * @param currentTime Current simulation time
     */
    void updateMembers(
            const double panelLatitude,
            const double panelLongitude,
            const double currentTime);

    /*!
     * Clone this object. Surface property distributions, e.g., for albedo and
     * emissivity, should be shared by the clone.
     *
     * @return A clone of this object
     */
    virtual std::unique_ptr<SourcePanelRadiosityModel> clone() const = 0;

protected:
    virtual void updateMembers_(
            const double panelLatitude,
            const double panelLongitude,
            const double currentTime) {};

    /*!
     * Whether the radiosity model is invariant with time. If yes, its members will not be updated even if the time
     * changed, it will only be updated when the latitude/longitude change. The time invariance is usually determined
     * by the albedo/emissivity distribution.
     *
     * @return Whether the radiosity model is invariant with time
     */
    virtual bool isTimeInvariant() = 0;

    double currentTime_{TUDAT_NAN};
    double panelLatitude_{TUDAT_NAN};
    double panelLongitude_{TUDAT_NAN};
};

/*!
 * Panel radiosity model with constant Lambertian radiosity.
 */
class ConstantSourcePanelRadiosityModel : public SourcePanelRadiosityModel
{
public:
    /*!
     * Constructor
     *
     * @param constantRadiosity Constant radiosity
     */
    explicit ConstantSourcePanelRadiosityModel(const double constantRadiosity) :
        constantRadiosity_(constantRadiosity) {}

    double evaluateIrradianceAtPosition(
            double panelArea,
            const Eigen::Vector3d& panelSurfaceNormal,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    std::unique_ptr<SourcePanelRadiosityModel> clone() const override
    {
        return std::make_unique<ConstantSourcePanelRadiosityModel>(*this);
    }

    double getConstantRadiosity() const
    {
        return constantRadiosity_;
    }

private:
    bool isTimeInvariant() override
    {
        return true;
    }

    const double constantRadiosity_;
};

/*!
 * Panel radiosity model for albedo radiation. This model was introduced in Knocke (1988) for Earth thermal radiation,
 * assuming Lambertian reflectance.
 *
 * For most cases, albedo radiation with a diffuse-only Lambertian reflection law is sufficient. Only a small fraction
 * of Earth's albedo is actually specular, since only calm bodies of water are truly specular and specular reflection
 * occurs mostly at low solar zenith angles where the reflectance is low anyway (Knocke, 1988).
 *
 * More sophisticated reflection laws (e.g., SpecularDiffuseMixReflectionLaw, or any BRDF via a custom ReflectionLaw)
 * can be implemented to take into account surface properties like different vegetation and ground types on Earth.
 * A separate panel radiosity model has to be created for these laws.
 */
class AlbedoSourcePanelRadiosityModel : public SourcePanelRadiosityModel
{
public:
    /*!
     * Constructor.
     *
     * @param albedoDistribution Albedo distribution
     */
    explicit AlbedoSourcePanelRadiosityModel(
            const std::shared_ptr<SurfacePropertyDistribution>& albedoDistribution) :
            albedoDistribution_(albedoDistribution),
            reflectionLaw_(std::make_shared<LambertianReflectionLaw>(TUDAT_NAN)) {}

    double evaluateIrradianceAtPosition(
            double panelArea,
            const Eigen::Vector3d& panelSurfaceNormal,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    std::unique_ptr<SourcePanelRadiosityModel> clone() const override
    {
        return std::make_unique<AlbedoSourcePanelRadiosityModel>(*this);
    }

    const std::shared_ptr<LambertianReflectionLaw>& getReflectionLaw() const
    {
        return reflectionLaw_;
    }

private:
    void updateMembers_(
            double panelLatitude,
            double panelLongitude,
            double currentTime) override;

    bool isTimeInvariant() override
    {
        return albedoDistribution_->isTimeInvariant();
    }

    std::shared_ptr<SurfacePropertyDistribution> albedoDistribution_;

    // Reflection law governing reflection of original source radiation
    std::shared_ptr<LambertianReflectionLaw> reflectionLaw_;
};

/*!
 * Panel radiosity model for thermal emissions based on delayed, isotropic and constant flux. This model was introduced
 * in Knocke (1988) for Earth thermal radiation.
 *
 * As opposed to instantaneous reradiation, the body is assumed to act as heat buffer, absorbing incoming radiation and
 * reradiating it in a delayed fashion as longwave infrared radiation. For most bodies, and especially Earth, this is a
 * good assumption (Knocke, 1988).
 */
class DelayedThermalSourcePanelRadiosityModel : public SourcePanelRadiosityModel
{
public:
    /*!
     * Constructor.
     *
     * @param emissivityDistribution Emissivity distribution
     */
    explicit DelayedThermalSourcePanelRadiosityModel(
            const std::shared_ptr<SurfacePropertyDistribution>& emissivityDistribution) :
            emissivityDistribution_(emissivityDistribution) {}

    double evaluateIrradianceAtPosition(
            double panelArea,
            const Eigen::Vector3d& panelSurfaceNormal,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    std::unique_ptr<SourcePanelRadiosityModel> clone() const override
    {
        return std::make_unique<DelayedThermalSourcePanelRadiosityModel>(*this);
    }

    double getEmissivity() const
    {
        return emissivity;
    }

private:
    void updateMembers_(
            double panelLatitude,
            double panelLongitude,
            double currentTime) override;

    bool isTimeInvariant() override
    {
        return emissivityDistribution_->isTimeInvariant();
    }

    std::shared_ptr<SurfacePropertyDistribution> emissivityDistribution_;
    double emissivity{TUDAT_NAN};
};

/*!
 * Panel radiosity model for thermal emissions based on the angle to sub-solar point and emissivity-corrected
 * black-body radiation. This model was introduced in Lemoine (2013) for lunar thermal radiation with
 * minTemperature = 100 K and maxTemperature = 375 K.
 *
 * The surface temperature is approximated by interpolation between minTemperature and maxTemperature based on the angle
 * to the subsolar point. From the surface temperature, the black-body radiation is then calculated and corrected for
 * emissivity. The radiation is maximum if the target is above the sub-solar point. At positions away from the sub-solar
 * point, the temperature drops, until it reaches minTemperature on the backside.
 */
class AngleBasedThermalSourcePanelRadiosityModel : public SourcePanelRadiosityModel
{
public:
    /*!
     * Constructor.
     *
     * @param minTemperature Minimum surface temperature (in shade) [K]
     * @param maxTemperature Maximum surface temperature (at subsolar point) [K]
     * @param emissivityDistribution Emissivity distribution
     */
    explicit AngleBasedThermalSourcePanelRadiosityModel(
            double minTemperature,
            double maxTemperature,
            const std::shared_ptr<SurfacePropertyDistribution>& emissivityDistribution) :
            minTemperature_(minTemperature),
            maxTemperature_(maxTemperature),
            emissivityDistribution_(emissivityDistribution) {}

    double evaluateIrradianceAtPosition(
            double panelArea,
            const Eigen::Vector3d& panelSurfaceNormal,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;


    std::unique_ptr<SourcePanelRadiosityModel> clone() const override
    {
        return std::make_unique<AngleBasedThermalSourcePanelRadiosityModel>(*this);
    }

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
        return emissivity;
    }

private:
    void updateMembers_(
            double panelLatitude,
            double panelLongitude,
            double currentTime) override;

    bool isTimeInvariant() override
    {
        return emissivityDistribution_->isTimeInvariant();
    }

    double minTemperature_;
    double maxTemperature_;
    std::shared_ptr<SurfacePropertyDistribution> emissivityDistribution_;
    double emissivity{TUDAT_NAN};
};

} // tudat
} // electromagnetism

#endif //TUDAT_SOURCESourcePanelRadiosityModel_H
