#ifndef TUDAT_CREATERADIATIONSOURCEMODEL_H
#define TUDAT_CREATERADIATIONSOURCEMODEL_H

#include <memory>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"

namespace tudat
{
namespace simulation_setup
{

//! List of radiation source models available in simulations
/*!
 *  List of radiation source models available in simulations. Radiation source models
 *  not defined by this given enum cannot be used for automatic model setup.
 */
enum RadiationSourceModelType
{
    isotropic_point_source,
    statically_paneled_source
};

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

enum LuminosityModelType
{
    constant_radiant_power,
    irradiance_based_radiant_power
};

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

class ConstantLuminosityModelSettings : public LuminosityModelSettings
{
public:
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

class IrradianceBasedLuminosityModelSettings : public LuminosityModelSettings
{
public:
    explicit IrradianceBasedLuminosityModelSettings(
            const std::function<double()>& irradianceAtDistanceFunction,
            const double distance) :
            LuminosityModelSettings(irradiance_based_radiant_power),
            irradianceAtDistanceFunction_(irradianceAtDistanceFunction),
            distance_(distance) {}

    const std::function<double()>& getIrradianceAtDistanceFunction() const
    {
        return irradianceAtDistanceFunction_;
    }

    double getDistance() const
    {
        return distance_;
    }

private:
    std::function<double()> irradianceAtDistanceFunction_;
    double distance_;
};

class IsotropicPointRadiationSourceModelSettings : public RadiationSourceModelSettings
{
public:
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

inline std::shared_ptr<ConstantLuminosityModelSettings>
constantLuminosityModelSettings(double luminosity)
{
    return std::make_shared< ConstantLuminosityModelSettings >(luminosity);
}

inline std::shared_ptr<IrradianceBasedLuminosityModelSettings>
irradianceBasedLuminosityModelSettings(double irradianceAtDistance,
                                       double distance)
{
    return std::make_shared< IrradianceBasedLuminosityModelSettings >(
            [=] () { return irradianceAtDistance; },
            distance);
}

inline std::shared_ptr<IsotropicPointRadiationSourceModelSettings>
isotropicPointRadiationSourceModelSettings(const std::shared_ptr<LuminosityModelSettings>& luminosityModel)
{
    return std::make_shared< IsotropicPointRadiationSourceModelSettings >(luminosityModel);
}

//*********************************************************************************************
//   Paneled radiation source
//*********************************************************************************************

enum PanelRadiosityModelType
{
    albedo,
    thermal_delayed,
    thermal_angle_based
};

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

class AlbedoPanelRadiosityModelSettings : public PanelRadiosityModelSettings
{
public:
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

class DelayedThermalPanelRadiosityModelSettings : public PanelRadiosityModelSettings
{
public:
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

class AngleBasedThermalPanelRadiosityModelSettings : public PanelRadiosityModelSettings
{
public:
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

class StaticallyPaneledRadiationSourceModelSettings : public RadiationSourceModelSettings
{
public:
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

inline std::shared_ptr<AlbedoPanelRadiosityModelSettings>
albedoPanelRadiosityModelSettings(double albedo)
{
    return std::make_shared< AlbedoPanelRadiosityModelSettings >([=] (double, double) { return albedo; });

}

inline std::shared_ptr<DelayedThermalPanelRadiosityModelSettings>
delayedThermalPanelRadiosityModelSettings(double emissivity)
{
    return std::make_shared< DelayedThermalPanelRadiosityModelSettings >([=] (double, double) { return emissivity; });
}

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

inline std::shared_ptr<StaticallyPaneledRadiationSourceModelSettings>
staticallyPaneledRadiationSourceModelSettings(
        std::initializer_list<std::shared_ptr<PanelRadiosityModelSettings>> panelRadiosityModels,
        int numberOfPanels)
{
    return std::make_shared< StaticallyPaneledRadiationSourceModelSettings >(
            std::vector<std::shared_ptr<PanelRadiosityModelSettings>>(panelRadiosityModels),
            numberOfPanels);
}


std::shared_ptr<electromagnetism::LuminosityModel> createLuminosityModel(
        const std::shared_ptr< LuminosityModelSettings >& modelSettings,
        const std::string& body);

std::function<std::shared_ptr<electromagnetism::PaneledRadiationSourceModel::PanelRadiosityModel>(double, double)>
        createPanelRadiosityModel(
        const std::shared_ptr< PanelRadiosityModelSettings >& modelSettings,
        const std::string& body);

std::shared_ptr<electromagnetism::RadiationSourceModel> createRadiationSourceModel(
        const std::shared_ptr< RadiationSourceModelSettings >& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies);


} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONSOURCEMODEL_H
