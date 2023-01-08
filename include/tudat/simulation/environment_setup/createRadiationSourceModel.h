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

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"

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
irradianceBasedLuminosityModelSettings(double irradianceAtDistance,
                                       double distance)
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
     * @param albedoDistribution Function returning the albedo at a given polar and azimuth angle on the body
     */
    explicit AlbedoPanelRadiosityModelSettings(
            const std::function<double(double,double)>& albedoDistribution) :
            PanelRadiosityModelSettings(albedo),
            albedoDistribution_(albedoDistribution) {}

    const std::function<double(double,double)>& getAlbedoDistribution() const
    {
        return albedoDistribution_;
    }

private:
    std::function<double(double, double)> albedoDistribution_;
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
     * @param emissivityDistribution Function returning the emissivity at a given polar and azimuth angle on the body
     */
    explicit DelayedThermalPanelRadiosityModelSettings(
            const std::function<double(double,double)>& emissivityDistribution) :
            PanelRadiosityModelSettings(thermal_delayed),
            emissivityDistribution_(emissivityDistribution) {}

    const std::function<double(double,double)>& getEmissivityDistribution() const
    {
        return emissivityDistribution_;
    }

private:
    std::function<double(double, double)> emissivityDistribution_;
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
     * @param emissivityDistribution Function returning the emissivity at a given polar and azimuth angle on the body
     */
    explicit AngleBasedThermalPanelRadiosityModelSettings(
            double minTemperature,
            double maxTemperature,
            const std::function<double(double,double)>& emissivityDistribution) :
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

    const std::function<double(double,double)>& getEmissivityDistribution() const
    {
        return emissivityDistribution_;
    }

private:
    double minTemperature_;
    double maxTemperature_;
    std::function<double(double, double)> emissivityDistribution_;
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
     * @param panelRadiosityModelSettings Vector of settings for radiosity model of all panels
     * @param numberOfPanels Number of panels for automatic source discretization
     */
    explicit StaticallyPaneledRadiationSourceModelSettings(
            const std::vector<std::shared_ptr<PanelRadiosityModelSettings>>& panelRadiosityModelSettings,
            unsigned int numberOfPanels) :
            RadiationSourceModelSettings(statically_paneled_source),
            panelRadiosityModelSettings_(panelRadiosityModelSettings),
            numberOfPanels_(numberOfPanels) {}

    unsigned int getNumberOfPanels() const
    {
        return numberOfPanels_;
    }

    const std::vector<std::shared_ptr<PanelRadiosityModelSettings>>& getPanelRadiosityModelSettings() const
    {
        return panelRadiosityModelSettings_;
    }

private:
    std::vector<std::shared_ptr<PanelRadiosityModelSettings>> panelRadiosityModelSettings_;
    unsigned int numberOfPanels_;
};

// TODO-DOMINIK provide convenience functions for SH albedo/emissivity distributions

/*!
 * Create settings for an albedo panel radiosity model with same albedo at any point on surface.
 *
 * @param albedo Constant albedo
 * @return Shared pointer to settings for an albedo panel radiosity model
 */
inline std::shared_ptr<AlbedoPanelRadiosityModelSettings>
albedoPanelRadiosityModelSettings(double albedo)
{
    return std::make_shared< AlbedoPanelRadiosityModelSettings >([=] (double, double) { return albedo; });

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
    return std::make_shared< DelayedThermalPanelRadiosityModelSettings >([=] (double, double) { return emissivity; });
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
            [=] (double, double) { return emissivity; } );
}

/*!
 * Create settings for a statically paneled radiation source model.
 *
 * @param panelRadiosityModels List of settings for radiosity models of all panels
 * @param numberOfPanels Number of panels for automatic source discretization
 * @return Shared pointer to settings for a statically paneled radiation source model
 */
inline std::shared_ptr<StaticallyPaneledRadiationSourceModelSettings>
staticallyPaneledRadiationSourceModelSettings(
        std::initializer_list<std::shared_ptr<PanelRadiosityModelSettings>> panelRadiosityModels,
        int numberOfPanels)
{
    return std::make_shared< StaticallyPaneledRadiationSourceModelSettings >(
            std::vector<std::shared_ptr<PanelRadiosityModelSettings>>(panelRadiosityModels),
            numberOfPanels);
}

/*!
 * Create luminosity model from its settings.
 *
 * @param modelSettings Settings of the luminosity model
 * @param body Body to which the luminosity model belongs
 * @return Shared pointer to luminosity model
 */
std::shared_ptr<electromagnetism::LuminosityModel> createLuminosityModel(
        const std::shared_ptr< LuminosityModelSettings >& modelSettings,
        const std::string& body);

/*!
 * Create function returning panel radiosity model at a given polar and azimuth angle from its settings.
 *
 * @param modelSettings Settings of the panel radiosity model
 * @param body Body to which the panel radiosity model belongs
 * @return Function returning panel radiosity model at a given polar and azimuth angle
 */
std::function<std::shared_ptr<electromagnetism::PaneledRadiationSourceModel::PanelRadiosityModel>(double, double)>
        createPanelRadiosityModel(
        const std::shared_ptr< PanelRadiosityModelSettings >& modelSettings,
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
        const std::shared_ptr< RadiationSourceModelSettings >& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies);


} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONSOURCEMODEL_H
