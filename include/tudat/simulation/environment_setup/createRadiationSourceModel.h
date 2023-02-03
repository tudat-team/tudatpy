/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATERADIATIONSOURCEMODEL_H
#define TUDAT_CREATERADIATIONSOURCEMODEL_H

#include <memory>
#include <vector>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/astro/electromagnetism/luminosityModel.h"
#include "tudat/astro/electromagnetism/surfacePropertyDistribution.h"


namespace tudat
{
namespace simulation_setup
{

/*!
 * Types of radiation source models.
 */
enum RadiationSourceModelType
{
    isotropic_point_source,
    statically_paneled_source
};

/*!
 * Settings for a radiation source model.
 */
class RadiationSourceModelSettings
{
public:
    explicit RadiationSourceModelSettings(
            const RadiationSourceModelType radiationSourceModelType) :
            radiationSourceModelType_(radiationSourceModelType) {}

    virtual ~RadiationSourceModelSettings() = default;

    RadiationSourceModelType getRadiationSourceModelType() const
    {
        return radiationSourceModelType_;
    }

private:
    RadiationSourceModelType radiationSourceModelType_;
};

//*********************************************************************************************
//   Isotropic point radiation source
//*********************************************************************************************

/*!
 * Types of luminosity models.
 */
enum LuminosityModelType
{
    constant_radiant_power,
    irradiance_based_radiant_power
};

/*!
 * Settings for a luminosity model.
 */
class LuminosityModelSettings
{
public:
    explicit LuminosityModelSettings(
            const LuminosityModelType luminosityModelType) :
            luminosityModelType_(luminosityModelType) {}

    virtual ~LuminosityModelSettings() = default;

    LuminosityModelType getLuminosityModelType() const
    {
        return luminosityModelType_;
    }

private:
    LuminosityModelType luminosityModelType_;
};

/*!
 * Settings for a constant luminosity model.
 *
 * @see ConstantLuminosityModel
 */
class ConstantLuminosityModelSettings : public LuminosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param luminosity Constant luminosity of the source [W]
     */
    explicit ConstantLuminosityModelSettings(
            const double luminosity) :
            LuminosityModelSettings(constant_radiant_power),
            luminosity_(luminosity) {}

    double getLuminosity() const
    {
        return luminosity_;
    }

private:
    double luminosity_;
};

/*!
 * Settings for an irradiance-based luminosity model.
 *
 * @see IrradianceBasedLuminosityModel
 */
class IrradianceBasedLuminosityModelSettings : public LuminosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param irradianceAtDistanceFunction Function returning the irradiance at a given time [W/m²]
     * @param distance Distance from the source at which the irradiance was evaluated/measured
     */
    explicit IrradianceBasedLuminosityModelSettings(
            const std::function<double(double)>& irradianceAtDistanceFunction,
            const double distance) :
            LuminosityModelSettings(irradiance_based_radiant_power),
            irradianceAtDistanceFunction_(irradianceAtDistanceFunction),
            distance_(distance) {}

    const std::function<double(double)>& getIrradianceAtDistanceFunction() const
    {
        return irradianceAtDistanceFunction_;
    }

    double getDistance() const
    {
        return distance_;
    }

private:
    std::function<double(double)> irradianceAtDistanceFunction_;
    double distance_;
};

/*!
 * Settings for an isotropic point radiation source model.
 */
class IsotropicPointRadiationSourceModelSettings : public RadiationSourceModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param luminosityModelSettings Luminosity of this source
     */
    explicit IsotropicPointRadiationSourceModelSettings(
            const std::shared_ptr<LuminosityModelSettings>& luminosityModelSettings) :
            RadiationSourceModelSettings(isotropic_point_source),
            luminosityModelSettings_(luminosityModelSettings) {}

    const std::shared_ptr<LuminosityModelSettings>& getLuminosityModelSettings() const
    {
        return luminosityModelSettings_;
    }

private:
    std::shared_ptr<LuminosityModelSettings> luminosityModelSettings_;
};

/*!
 * Create settings for a constant luminosity model.
 *
 * @param luminosity Luminosity of this source
 * @return Shared pointer to settings for a constant luminosity model.
 */
inline std::shared_ptr<ConstantLuminosityModelSettings>
        constantLuminosityModelSettings(double luminosity)
{
    return std::make_shared< ConstantLuminosityModelSettings >(luminosity);
}

/*!
 * Create settings for an irradiance-based luminosity model.
 *
 * @param irradianceAtDistance Irradiance [W/m²]
 * @param distance Distance from the source at which the irradiance was evaluated/measured
 * @return Shared pointer to settings for an irradiance-based luminosity model.
 */
inline std::shared_ptr<IrradianceBasedLuminosityModelSettings>
        irradianceBasedLuminosityModelSettings(double irradianceAtDistance, double distance)
{
    return std::make_shared< IrradianceBasedLuminosityModelSettings >(
            [=] (double) { return irradianceAtDistance; },
            distance);
}

/*!
 * Create settings for an isotropic point radiation source model.
 *
 * @param luminosityModel Luminosity of this source
 * @return Shared pointer to settings for an isotropic point radiation source model.
 */
inline std::shared_ptr<IsotropicPointRadiationSourceModelSettings>
        isotropicPointRadiationSourceModelSettings(const std::shared_ptr<LuminosityModelSettings>& luminosityModel)
{
    return std::make_shared< IsotropicPointRadiationSourceModelSettings >(luminosityModel);
}

//*********************************************************************************************
//   Paneled radiation source
//*********************************************************************************************

/*!
 * Types of surface property distributions.
 */
enum SurfacePropertyDistributionType
{
    constant,
    spherical_harmonics
};

/*!
 * Settings for distribution of a property on the surface of a sphere such as albedo or emissivity.
 *
 * @see SurfacePropertyDistribution
 */
class SurfacePropertyDistributionSettings
{
public:
    explicit SurfacePropertyDistributionSettings(
            const SurfacePropertyDistributionType surfacePropertyDistributionType) :
            surfacePropertyDistributionType_(surfacePropertyDistributionType) {}

    virtual ~SurfacePropertyDistributionSettings() = default;

    SurfacePropertyDistributionType getSurfacePropertyDistributionType() const
    {
        return surfacePropertyDistributionType_;
    }

private:
    SurfacePropertyDistributionType surfacePropertyDistributionType_;
};

/*!
 * Settings for a constant surface property distribution.
 *
 * @see ConstantSurfacePropertyDistribution
 */
class ConstantSurfacePropertyDistributionSettings : public SurfacePropertyDistributionSettings
{
public:
    /*!
     * Constructor.
     *
     * @param constantValue Constant value
     */
    explicit ConstantSurfacePropertyDistributionSettings(
            const double constantValue) :
            SurfacePropertyDistributionSettings(constant),
            constantValue_(constantValue) {}

    double getConstantValue() const
    {
        return constantValue_;
    }

private:
    double constantValue_;
};

enum SphericalHarmonicsSurfacePropertyDistributionModel
{
    custom,
    albedo_dlam1 /**< DLAM-1 lunar albedo model: Floberghagen, R. et al. "Lunar Albedo Force Modeling and its Effect on Low Lunar Orbit and Gravity Field Determination". ASR 23. 4(1999): 733-738. */
};

/*!
 * Settings for a surface property distribution described by a spherical harmonics expansion. The reference frame of the
 * body and spherical harmonics must be identical.
 *
 * @see SphericalHarmonicsSurfacePropertyDistribution
 */
class SphericalHarmonicsSurfacePropertyDistributionSettings : public SurfacePropertyDistributionSettings
{
public:
    /*!
     * Constructor with custom model.
     *
     * @param cosineCoefficients Cosine spherical harmonic coefficients (not normalized)
     * @param sineCoefficients Sine spherical harmonic coefficients (not normalized)
     */
    explicit SphericalHarmonicsSurfacePropertyDistributionSettings(
            const Eigen::MatrixXd& cosineCoefficients,
            const Eigen::MatrixXd& sineCoefficients) :
            SurfacePropertyDistributionSettings(spherical_harmonics),
            model_(custom),
            cosineCoefficients_(cosineCoefficients),
            sineCoefficients_(sineCoefficients) {}

    /*!
    * Constructor with model included in Tudat.
    *
    * @param model Spherical harmonics model to be used
    */
    explicit SphericalHarmonicsSurfacePropertyDistributionSettings(
            SphericalHarmonicsSurfacePropertyDistributionModel model);

    const Eigen::MatrixXd& getCosineCoefficients() const
    {
        return cosineCoefficients_;
    }

    const Eigen::MatrixXd& getSineCoefficients() const
    {
        return sineCoefficients_;
    }

private:
    SphericalHarmonicsSurfacePropertyDistributionModel model_;

    // Cosine spherical harmonic coefficients (not normalized)
    Eigen::MatrixXd cosineCoefficients_;

    // Sine spherical harmonic coefficients (not normalized)
    Eigen::MatrixXd sineCoefficients_;
};

/*!
 * Types of panel radiosity models.
 */
enum PanelRadiosityModelType
{
    albedo,
    thermal_delayed,
    thermal_angle_based
};

/*!
 * Settings for a radiosity model of a paneled source panel.
 *
 * @see PaneledRadiationSourceModel::PanelRadiosityModel
 */
class PanelRadiosityModelSettings
{
public:
    explicit PanelRadiosityModelSettings(
            const PanelRadiosityModelType panelRadiosityModelType) :
            panelRadiosityModelType_(panelRadiosityModelType) {}

    virtual ~PanelRadiosityModelSettings() = default;

    PanelRadiosityModelType getPanelRadiosityModelType() const
    {
        return panelRadiosityModelType_;
    }

private:
    PanelRadiosityModelType panelRadiosityModelType_;
};

/*!
 * Settings for an albedo radiosity model of a paneled source panel.
 *
 * The panel will diffusely reflect incoming radiation with a Lambertian reflectance law, using the albedo as
 * diffuse reflectivity coefficient.
 *
 * @see AlbedoPanelRadiosityModel
 * @see LambertianReflectanceLaw
 */
class AlbedoPanelRadiosityModelSettings : public PanelRadiosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param albedoDistribution Albedo distribution
     * @param withInstantaneousReradiation Whether to instantaneously reradiate absorbed radiation
     */
    explicit AlbedoPanelRadiosityModelSettings(
            const std::shared_ptr<SurfacePropertyDistributionSettings>& albedoDistribution,
            bool withInstantaneousReradiation) :
            PanelRadiosityModelSettings(albedo),
            albedoDistribution_(albedoDistribution),
            withInstantaneousReradiation_(withInstantaneousReradiation) {}

    const std::shared_ptr<SurfacePropertyDistributionSettings>& getAlbedoDistribution() const
    {
        return albedoDistribution_;
    }

    bool getWithInstantaneousReradiation() const
    {
        return withInstantaneousReradiation_;
    }

private:
    std::shared_ptr<SurfacePropertyDistributionSettings> albedoDistribution_;
    bool withInstantaneousReradiation_;
};

/*!
 * Settings for a delayed thermal radiosity model of a paneled source panel.
 *
 * @see DelayedThermalPanelRadiosityModel
 */
class DelayedThermalPanelRadiosityModelSettings : public PanelRadiosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param emissivityDistribution Emissivity distribution
     */
    explicit DelayedThermalPanelRadiosityModelSettings(
            const std::shared_ptr<SurfacePropertyDistributionSettings>& emissivityDistribution) :
            PanelRadiosityModelSettings(thermal_delayed),
            emissivityDistribution_(emissivityDistribution) {}

    const std::shared_ptr<SurfacePropertyDistributionSettings>& getEmissivityDistribution() const
    {
        return emissivityDistribution_;
    }

private:
    std::shared_ptr<SurfacePropertyDistributionSettings> emissivityDistribution_;
};

/*!
 * Settings for an angle-based thermal radiosity model of a paneled source panel.
 *
 * @see AngleBasedThermalPanelRadiosityModel
 */
class AngleBasedThermalPanelRadiosityModelSettings : public PanelRadiosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param minTemperature Minimum surface temperature (in shade) [K]
     * @param maxTemperature Maximum surface temperature (at subsolar point) [K]
     * @param emissivityDistribution Emissivity distribution
     */
    explicit AngleBasedThermalPanelRadiosityModelSettings(
            double minTemperature,
            double maxTemperature,
            const std::shared_ptr<SurfacePropertyDistributionSettings>& emissivityDistribution) :
            PanelRadiosityModelSettings(thermal_angle_based),
            minTemperature_(minTemperature),
            maxTemperature_(maxTemperature),
            emissivityDistribution_(emissivityDistribution) {}

    double getMinTemperature() const
    {
        return minTemperature_;
    }

    double getMaxTemperature() const
    {
        return maxTemperature_;
    }

    const std::shared_ptr<SurfacePropertyDistributionSettings>& getEmissivityDistribution() const
    {
        return emissivityDistribution_;
    }

private:
    double minTemperature_;
    double maxTemperature_;
    std::shared_ptr<SurfacePropertyDistributionSettings> emissivityDistribution_;
};

/*!
 * Settings for a statically paneled radiation source model.
 */
class StaticallyPaneledRadiationSourceModelSettings : public RadiationSourceModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param originalSourceName Name of the original source body
     * @param panelRadiosityModelSettings Vector of settings for radiosity model of all panels
     * @param numberOfPanels Number of panels for automatic source discretization
     * @param originalSourceToSourceOccultingBodies Names of bodies to occult the original source as seen from this source
     */
    explicit StaticallyPaneledRadiationSourceModelSettings(
            const std::string& originalSourceName,
            const std::vector<std::shared_ptr<PanelRadiosityModelSettings>>& panelRadiosityModelSettings,
            unsigned int numberOfPanels,
            const std::vector<std::string>& originalSourceToSourceOccultingBodies = {}) :
            RadiationSourceModelSettings(statically_paneled_source),
            originalSourceName_(originalSourceName),
            panelRadiosityModelSettings_(panelRadiosityModelSettings),
            numberOfPanels_(numberOfPanels),
            originalSourceToSourceOccultingBodies_(originalSourceToSourceOccultingBodies) {}

    std::string getOriginalSourceName() const
    {
        return originalSourceName_;
    }

    unsigned int getNumberOfPanels() const
    {
        return numberOfPanels_;
    }

    const std::vector<std::shared_ptr<PanelRadiosityModelSettings>>& getPanelRadiosityModelSettings() const
    {
        return panelRadiosityModelSettings_;
    }

    std::vector<std::string> getOriginalSourceToSourceOccultingBodies() const
    {
        return originalSourceToSourceOccultingBodies_;
    }

private:
    std::string originalSourceName_;
    std::vector<std::shared_ptr<PanelRadiosityModelSettings>> panelRadiosityModelSettings_;
    unsigned int numberOfPanels_;
    std::vector<std::string> originalSourceToSourceOccultingBodies_;
};

/*!
 * Create settings for constant surface property distribution.
 *
 * @param constantValue Constant value
 * @return Shared pointer to settings for a constant surface property distribution.
 */
inline std::shared_ptr<ConstantSurfacePropertyDistributionSettings>
        constantSurfacePropertyDistributionSettings(double constantValue)
{
    return std::make_shared< ConstantSurfacePropertyDistributionSettings >(constantValue);
}

/*!
 * Create settings for constant surface property distribution from coefficients.
 *
 * @param cosineCoefficients Cosine spherical harmonic coefficients (not normalized)
 * @param sineCoefficients Sine spherical harmonic coefficients (not normalized)
 * @return Shared pointer to settings for a spherical harmonics surface property distribution.
 */
inline std::shared_ptr<SphericalHarmonicsSurfacePropertyDistributionSettings>
        sphericalHarmonicsSurfacePropertyDistributionSettings(
                const Eigen::MatrixXd& cosineCoefficients,
                const Eigen::MatrixXd& sineCoefficients)
{
    return std::make_shared< SphericalHarmonicsSurfacePropertyDistributionSettings >(
            cosineCoefficients, sineCoefficients);
}

/*!
 * Create settings for constant surface property distribution from model included in Tudat.
 *
 * @param model Spherical harmonics model to be used
 * @return Shared pointer to settings for a spherical harmonics surface property distribution.
 */
inline std::shared_ptr<SphericalHarmonicsSurfacePropertyDistributionSettings>
        sphericalHarmonicsSurfacePropertyDistributionSettings(
                SphericalHarmonicsSurfacePropertyDistributionModel model)
{
    return std::make_shared< SphericalHarmonicsSurfacePropertyDistributionSettings >(model);
}

/*!
 * Create settings for an albedo panel radiosity model with same albedo at any point on surface.
 *
 * @param albedo Constant albedo
 * @param withInstantaneousReradiation Whether to instantaneously reradiate absorbed radiation
 * @return Shared pointer to settings for an albedo panel radiosity model
 */
inline std::shared_ptr<AlbedoPanelRadiosityModelSettings>
        albedoPanelRadiosityModelSettings(double albedo, bool withInstantaneousReradiation = false)
{
    return std::make_shared< AlbedoPanelRadiosityModelSettings >(
        constantSurfacePropertyDistributionSettings(albedo), withInstantaneousReradiation);

}

/*!
 * Create settings for an albedo panel radiosity model with spherical harmonics albedo.
 *
 * @param albedo Constant albedo
 * @param withInstantaneousReradiation Whether to instantaneously reradiate absorbed radiation
 * @return Shared pointer to settings for an albedo panel radiosity model
 */
inline std::shared_ptr<AlbedoPanelRadiosityModelSettings>
        albedoPanelRadiosityModelSettings(
                SphericalHarmonicsSurfacePropertyDistributionModel albedoModel,
                bool withInstantaneousReradiation = false)
{
    return std::make_shared< AlbedoPanelRadiosityModelSettings >(
            sphericalHarmonicsSurfacePropertyDistributionSettings(albedoModel), withInstantaneousReradiation);

}

/*!
 * Create settings for a delayed thermal panel radiosity model with same emissivity at any point on surface.
 *
 * @param emissivity Constant emissivity
 * @return Shared pointer to settings for a delayed thermal panel radiosity model
 */
inline std::shared_ptr<DelayedThermalPanelRadiosityModelSettings>
        delayedThermalPanelRadiosityModelSettings(double emissivity)
{
    return std::make_shared< DelayedThermalPanelRadiosityModelSettings >(
            constantSurfacePropertyDistributionSettings(emissivity));
}

/*!
 * Create settings for an angle-based thermal panel radiosity model with same emissivity at any point on surface.
 *
 * @param minTemperature Minimum surface temperature (in shade) [K]
 * @param maxTemperature Maximum surface temperature (at subsolar point) [K]
 * @param emissivity Constant emissivity
 * @return Shared pointer to settings for an angle-based panel radiosity model
 */
inline std::shared_ptr<AngleBasedThermalPanelRadiosityModelSettings>
        angleBasedThermalPanelRadiosityModelSettings(double minTemperature,
                                             double maxTemperature,
                                             double emissivity)
{
    return std::make_shared< AngleBasedThermalPanelRadiosityModelSettings >(
            minTemperature,
            maxTemperature,
            constantSurfacePropertyDistributionSettings(emissivity));
}

/*!
 * Create settings for a statically paneled radiation source model.
 *
 * @param originalSourceName Name of the original source body
 * @param panelRadiosityModels List of settings for radiosity models of all panels
 * @param numberOfPanels Number of panels for automatic source discretization
 * @param originalSourceToSourceOccultingBodies Names of bodies to occult the original source as seen from this source
 * @return Shared pointer to settings for a statically paneled radiation source model
 */
inline std::shared_ptr<StaticallyPaneledRadiationSourceModelSettings>
        staticallyPaneledRadiationSourceModelSettings(
                const std::string& originalSourceName,
                std::initializer_list<std::shared_ptr<PanelRadiosityModelSettings>> panelRadiosityModels,
                unsigned int numberOfPanels,
                const std::vector<std::string>& originalSourceToSourceOccultingBodies = {})
{
    return std::make_shared< StaticallyPaneledRadiationSourceModelSettings >(
            originalSourceName,
            std::vector<std::shared_ptr<PanelRadiosityModelSettings>>(panelRadiosityModels),
            numberOfPanels, originalSourceToSourceOccultingBodies);
}

/*!
 * Create luminosity model from its settings.
 *
 * @param modelSettings Settings of the luminosity model
 * @param body Body to which the luminosity model belongs
 * @return Shared pointer to luminosity model
 */
std::shared_ptr<electromagnetism::LuminosityModel> createLuminosityModel(
        const std::shared_ptr<LuminosityModelSettings>& modelSettings,
        const std::string& body);

/*!
 * Create surface property distribution from its settings.
 *
 * @param distributionSettings Settings of the surface property distribution
 * @param body Body to which the surface property distribution belongs
 * @return Shared pointer to surface property distribution
 */
std::shared_ptr<electromagnetism::SurfacePropertyDistribution> createSurfacePropertyDistribution(
        const std::shared_ptr<SurfacePropertyDistributionSettings>& distributionSettings,
        const std::string& body);

/*!
 * Create function returning panel radiosity model at a given polar and azimuth angle from its settings.
 *
 * @param modelSettings Settings of the panel radiosity model
 * @param body Body to which the panel radiosity model belongs
 * @return Panel radiosity model
 */
std::unique_ptr<electromagnetism::PaneledRadiationSourceModel::PanelRadiosityModel>
        createPanelRadiosityModel(
        const std::shared_ptr<PanelRadiosityModelSettings>& modelSettings,
        const std::string& body);

/*!
 * Create radiation source model from its settings.
 *
 * @param modelSettings Settings of the radiation source model
 * @param body Body to which the radiation source model belongs
 * @param bodies System of bodies
 * @return Shared pointer to radiation source model
 */
std::shared_ptr<electromagnetism::RadiationSourceModel> createRadiationSourceModel(
        const std::shared_ptr<RadiationSourceModelSettings >& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies);


} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONSOURCEMODEL_H
