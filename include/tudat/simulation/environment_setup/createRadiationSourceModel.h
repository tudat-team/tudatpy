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
#include <utility>
#include <vector>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/astro/electromagnetism/luminosityModel.h"
#include "tudat/simulation/environment_setup/createSurfacePropertyDistribution.h"


namespace tudat
{
namespace simulation_setup
{

/*!
 * Types of radiation source models.
 */
enum class RadiationSourceModelType
{
    isotropic_point_source,
    extended_source
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
enum class LuminosityModelType
{
    constant_radiant_power,
    time_variable_isotropic_radiant_power
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
            LuminosityModelSettings(LuminosityModelType::constant_radiant_power),
            luminosity_(luminosity) {}

    double getLuminosity() const
    {
        return luminosity_;
    }

private:
    double luminosity_;
};

class TimeVariableLuminosityModelSettings : public LuminosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param luminosity Constant luminosity of the source [W]
     */
    explicit TimeVariableLuminosityModelSettings(
        const std::function< double( const double ) > luminosityFunction) :
        LuminosityModelSettings(LuminosityModelType::time_variable_isotropic_radiant_power),
        luminosityFunction_(luminosityFunction) {}

    std::function< double( const double ) > getLuminosityFuntion() const
    {
        return luminosityFunction_;
    }

private:
    std::function< double( const double ) > luminosityFunction_;
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
            RadiationSourceModelSettings(RadiationSourceModelType::isotropic_point_source),
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
inline std::shared_ptr<LuminosityModelSettings>
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
inline std::shared_ptr<LuminosityModelSettings>
        irradianceBasedLuminosityModelSettings(double irradianceAtDistance, double distance)
{
    return std::make_shared< ConstantLuminosityModelSettings >(
        electromagnetism::computeLuminosityFromIrradiance( irradianceAtDistance, distance ) );
}

inline std::shared_ptr<LuminosityModelSettings>
timeVariableLuminosityModelSettings(const std::function< double( const double ) > luminosityFunction )
{
    return std::make_shared< TimeVariableLuminosityModelSettings >(luminosityFunction);
}


inline std::shared_ptr<LuminosityModelSettings>
timeVariableIrradianceBasedLuminosityModelSettings(const std::function< double( const double ) >  irradianceAtDistanceFunction, double distance)
{
    return std::make_shared< TimeVariableLuminosityModelSettings >(
        [=](const double time){ return electromagnetism::computeLuminosityFromIrradiance(
            irradianceAtDistanceFunction( time ), distance ); }
        );
}

/*!
 * Create settings for an isotropic point radiation source model.
 *
 * @param luminosityModel Luminosity of this source
 * @return Shared pointer to settings for an isotropic point radiation source model.
 */
inline std::shared_ptr<RadiationSourceModelSettings>
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
enum class PanelRadiosityModelType
{
    constant,
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
 * Settings for an inherent radiosity model of a paneled source panel.
 *
 * The radiosity is independent of an original source. This is the case for observed fluxes or internal (e.g., tidal)
 * heating.
 *
 * @see InherentSourcePanelRadiosityModel
 */
class InherentPanelRadiosityModelSettings : public PanelRadiosityModelSettings
{
public:
    explicit InherentPanelRadiosityModelSettings(
            const PanelRadiosityModelType panelRadiosityModelType) :
            PanelRadiosityModelSettings(panelRadiosityModelType) {}
};

/*!
 * Settings for an original-source-dependent radiosity model of a paneled source panel.
 *
 * The radiosity depends on an original source.  This is the case for albedo or thermal radiation, which require the
 * irradiance and direction of the incident solar radiation.
 *
 * @see OriginalSourceDependentSourcePanelRadiosityModel
 */
class OriginalSourceDependentPanelRadiosityModelSettings : public PanelRadiosityModelSettings
{
public:
    explicit OriginalSourceDependentPanelRadiosityModelSettings(
            const PanelRadiosityModelType panelRadiosityModelType,
            const std::string& originalSourceName) :
            PanelRadiosityModelSettings(panelRadiosityModelType),
            originalSourceName_(originalSourceName) {}

    const std::string& getOriginalSourceName() const
    {
        return originalSourceName_;
    }

private:
    std::string originalSourceName_;
};

/*!
 * Settings for a constant radiosity model of a paneled source panel.
 *
 * The panel will emit radiation with a Lambertian law.
 *
 * @see ConstantPanelRadiosityModel
 */
class ConstantPanelRadiosityModelSettings : public InherentPanelRadiosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param constantRadiosity Constant radiosity
     */
    explicit ConstantPanelRadiosityModelSettings(const double constantRadiosity) :
            InherentPanelRadiosityModelSettings(PanelRadiosityModelType::constant),
            constantRadiosity_(constantRadiosity) {}

    double getConstantRadiosity() const
    {
        return constantRadiosity_;
    }

private:
    const double constantRadiosity_;
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
class AlbedoPanelRadiosityModelSettings : public OriginalSourceDependentPanelRadiosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param albedoDistribution Albedo distribution
     * @param originalSourceName Name of the original source
     */
    explicit AlbedoPanelRadiosityModelSettings(
            const std::shared_ptr<SurfacePropertyDistributionSettings>& albedoDistribution,
            const std::string& originalSourceName) :
            OriginalSourceDependentPanelRadiosityModelSettings(
                    PanelRadiosityModelType::albedo, originalSourceName),
            albedoDistribution_(albedoDistribution) {}

    const std::shared_ptr<SurfacePropertyDistributionSettings>& getAlbedoDistribution() const
    {
        return albedoDistribution_;
    }

private:
    std::shared_ptr<SurfacePropertyDistributionSettings> albedoDistribution_;
};

/*!
 * Settings for a delayed thermal radiosity model of a paneled source panel.
 *
 * @see DelayedThermalPanelRadiosityModel
 */
class DelayedThermalPanelRadiosityModelSettings : public OriginalSourceDependentPanelRadiosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param emissivityDistribution Emissivity distribution
     * @param originalSourceName Name of the original source
     */
    explicit DelayedThermalPanelRadiosityModelSettings(
            const std::shared_ptr<SurfacePropertyDistributionSettings>& emissivityDistribution,
            const std::string& originalSourceName) :
            OriginalSourceDependentPanelRadiosityModelSettings(
                    PanelRadiosityModelType::thermal_delayed, originalSourceName),
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
class AngleBasedThermalPanelRadiosityModelSettings : public OriginalSourceDependentPanelRadiosityModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param minTemperature Minimum surface temperature (in shade) [K]
     * @param maxTemperature Maximum surface temperature (at subsolar point) [K]
     * @param emissivityDistribution Emissivity distribution
     * @param originalSourceName Name of the original source
     */
    explicit AngleBasedThermalPanelRadiosityModelSettings(
            double minTemperature,
            double maxTemperature,
            const std::shared_ptr<SurfacePropertyDistributionSettings>& emissivityDistribution,
            const std::string& originalSourceName) :
            OriginalSourceDependentPanelRadiosityModelSettings(
                    PanelRadiosityModelType::thermal_angle_based, originalSourceName),
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
 * Settings for an extended (dynamically paneled) radiation source model.
 *
 * @see DynamicallyPaneledRadiationSourceModel
 */
class ExtendedRadiationSourceModelSettings : public RadiationSourceModelSettings
{
public:
    /*!
     * Constructor.
     *
     * @param panelRadiosityModelSettings Vector of settings for radiosity model of all panels
     * @param numberOfPanelsPerRing Number of panels for each ring, excluding the central cap
     * @param originalSourceToSourceOccultingBodies Names of bodies to occult original sources as seen from this source
     */
    explicit ExtendedRadiationSourceModelSettings(
            const std::vector<std::shared_ptr<PanelRadiosityModelSettings>>& panelRadiosityModelSettings,
            const std::vector<int>& numberOfPanelsPerRing,
            const std::map<std::string, std::vector<std::string>>& originalSourceToSourceOccultingBodies) :
            RadiationSourceModelSettings(RadiationSourceModelType::extended_source),
            panelRadiosityModelSettings_(panelRadiosityModelSettings),
            numberOfPanelsPerRing_(numberOfPanelsPerRing),
            originalSourceToSourceOccultingBodies_(originalSourceToSourceOccultingBodies) {}

    const std::vector<int>& getNumberOfPanelsPerRing() const
    {
        return numberOfPanelsPerRing_;
    }

    const std::vector<std::shared_ptr<PanelRadiosityModelSettings>>& getPanelRadiosityModelSettings() const
    {
        return panelRadiosityModelSettings_;
    }

    const std::map<std::string, std::vector<std::string>>& getOriginalSourceToSourceOccultingBodies() const
    {
        return originalSourceToSourceOccultingBodies_;
    }

private:
    std::vector<std::shared_ptr<PanelRadiosityModelSettings>> panelRadiosityModelSettings_;
    const std::vector<int> numberOfPanelsPerRing_;
    // Map (original source name -> list of occulting body names) of bodies to occult original sources as seen from this source
    // If the same occulting bodies are to be used for all original sources, there will be a single entry
    // with an emptry string as key
    std::map<std::string, std::vector<std::string>> originalSourceToSourceOccultingBodies_;
};

/*!
 * Create settings for a constant panel radiosity model.
 *
 * @param constantRadiosity Constant radiosity
 * @return Shared pointer to settings for a constant panel radiosity model
 */
inline std::shared_ptr<PanelRadiosityModelSettings>
        constantPanelRadiosityModelSettings(double constantRadiosity)
{
    return std::make_shared< ConstantPanelRadiosityModelSettings >(constantRadiosity);
}

/*!
 * Create settings for an albedo panel radiosity model with same albedo at any point on surface.
 *
 * @param albedo Constant albedo
 * @param originalSourceName Name of the original source
 * @return Shared pointer to settings for an albedo panel radiosity model
 */
inline std::shared_ptr<PanelRadiosityModelSettings> albedoPanelRadiosityModelSettings(
        double albedo,
        const std::string& originalSourceName)
{
    return std::make_shared< AlbedoPanelRadiosityModelSettings >(
            constantSurfacePropertyDistributionSettings(albedo),
            originalSourceName);
}

/*!
 * Create settings for an albedo panel radiosity model with spherical harmonics albedo from model included in Tudat.
 *
 * @param albedoModel Spherical harmonics model to be used
 * @param originalSourceName Name of the original source
 * @return Shared pointer to settings for an albedo panel radiosity model
 */
inline std::shared_ptr<PanelRadiosityModelSettings>
        albedoPanelRadiosityModelSettings(
                SphericalHarmonicsSurfacePropertyDistributionModel albedoModel,
                const std::string& originalSourceName)
{
    return std::make_shared< AlbedoPanelRadiosityModelSettings >(
            sphericalHarmonicsSurfacePropertyDistributionSettings(albedoModel),
            originalSourceName);
}

/*!
 * Create settings for an albedo panel radiosity model with second-degree zonal periodic spherical harmonics albedo
 * from model included in Tudat.
 *
 * @param albedoModel Model to be used
 * @param originalSourceName Name of the original source
 * @return Shared pointer to settings for an albedo panel radiosity model
 */
inline std::shared_ptr<PanelRadiosityModelSettings>
        albedoPanelRadiosityModelSettings(
                KnockeTypeSurfacePropertyDistributionModel albedoModel,
                const std::string& originalSourceName)
{
    return std::make_shared< AlbedoPanelRadiosityModelSettings >(
            secondDegreeZonalPeriodicSurfacePropertyDistributionSettings(albedoModel),
            originalSourceName);
}

inline std::shared_ptr<PanelRadiosityModelSettings>
albedoPanelRadiosityModelSettingsGeneric(
    const std::shared_ptr<SurfacePropertyDistributionSettings>& albedoDistribution,
    const std::string& originalSourceName)
{

    return std::make_shared< AlbedoPanelRadiosityModelSettings >(
        albedoDistribution,
        originalSourceName);
}

/*!
 * Create settings for a delayed thermal panel radiosity model with same emissivity at any point on surface.
 *
 * @param emissivity Constant emissivity
 * @param originalSourceName Name of the original source
 * @return Shared pointer to settings for a delayed thermal panel radiosity model
 */
inline std::shared_ptr<PanelRadiosityModelSettings>
        delayedThermalPanelRadiosityModelSettings(
                double emissivity,
                const std::string& originalSourceName)
{
    return std::make_shared< DelayedThermalPanelRadiosityModelSettings >(
            constantSurfacePropertyDistributionSettings(emissivity),
            originalSourceName);
}

/*!
 * Create settings for a delayed thermal panel radiosity model with second-degree zonal periodic spherical harmonics
 * emissivity from model included in Tudat.
 *
 * @param emissivityModel Model to be used
 * @param originalSourceName Name of the original source
 * @return Shared pointer to settings for a delayed thermal panel radiosity model
 */
inline std::shared_ptr<PanelRadiosityModelSettings>
        delayedThermalPanelRadiosityModelSettings(
                KnockeTypeSurfacePropertyDistributionModel emissivityModel,
                const std::string& originalSourceName)
{
    return std::make_shared< DelayedThermalPanelRadiosityModelSettings >(
            secondDegreeZonalPeriodicSurfacePropertyDistributionSettings(emissivityModel),
            originalSourceName);
}

inline std::shared_ptr<PanelRadiosityModelSettings>
delayedThermalPanelRadiosityModelSettingsGeneric(
    const std::shared_ptr<SurfacePropertyDistributionSettings>& emissivityDistribution,
    const std::string& originalSourceName)
{
    return std::make_shared< DelayedThermalPanelRadiosityModelSettings >(
        emissivityDistribution,
        originalSourceName);
}


/*!
 * Create settings for an angle-based thermal panel radiosity model with same emissivity at any point on surface.
 *
 * @param minTemperature Minimum surface temperature (in shade) [K]
 * @param maxTemperature Maximum surface temperature (at subsolar point) [K]
 * @param emissivity Constant emissivity
 * @param originalSourceName Name of the original source
 * @return Shared pointer to settings for an angle-based panel radiosity model
 */
inline std::shared_ptr<PanelRadiosityModelSettings>
        angleBasedThermalPanelRadiosityModelSettings(
                double minTemperature,
                double maxTemperature,
                double emissivity,
                const std::string& originalSourceName)
{
    return std::make_shared< AngleBasedThermalPanelRadiosityModelSettings >(
            minTemperature,
            maxTemperature,
            constantSurfacePropertyDistributionSettings(emissivity),
            originalSourceName);
}

/*!
 * Create settings for an extended (dynamically paneled) radiation source model.
 *
 * @param panelRadiosityModels List of settings for radiosity models of all panels
 * @param numberOfPanelsPerRing Number of panels for each ring, excluding the central cap
 * @param originalSourceToSourceOccultingBodies Map (original source name -> list of occulting body names) of bodies
 *      to occult original sources as seen from this source
 * @return Shared pointer to settings for an extended radiation source model
 */
inline std::shared_ptr<RadiationSourceModelSettings>
        extendedRadiationSourceModelSettingsWithOccultationMap(
                std::vector<std::shared_ptr<PanelRadiosityModelSettings>> panelRadiosityModels,
                const std::vector<int>& numberOfPanelsPerRing,
                const std::map<std::string, std::vector<std::string>>& originalSourceToSourceOccultingBodies)
{
    return std::make_shared< ExtendedRadiationSourceModelSettings >(
            panelRadiosityModels, numberOfPanelsPerRing, originalSourceToSourceOccultingBodies);
}

/*!
 * Create settings for an extended (dynamically paneled) radiation source model.
 *
 * @param panelRadiosityModels List of settings for radiosity models of all panels
 * @param numberOfPanelsPerRing Number of panels for each ring, excluding the central cap
 * @param originalSourceToSourceOccultingBodies Names of bodies to occult original sources as seen from this source
 * @return Shared pointer to settings for an extended radiation source model
 */
inline std::shared_ptr<RadiationSourceModelSettings>
        extendedRadiationSourceModelSettings(
                std::vector<std::shared_ptr<PanelRadiosityModelSettings>> panelRadiosityModels,
                const std::vector<int>& numberOfPanelsPerRing,
                const std::vector<std::string>& originalSourceToSourceOccultingBodies = {})
{
    const std::map<std::string, std::vector<std::string>> occultingBodiesMap {{"", originalSourceToSourceOccultingBodies}};
    return extendedRadiationSourceModelSettingsWithOccultationMap(
            std::move(panelRadiosityModels), numberOfPanelsPerRing, occultingBodiesMap);
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
 * Create panel radiosity model from its settings.
 *
 * @param modelSettings Settings of the panel radiosity model
 * @param sourceBodyName Body to which the panel radiosity model belongs
 * @return Unique pointer to panel radiosity model
 */
std::unique_ptr<electromagnetism::SourcePanelRadiosityModel> createPanelRadiosityModel(
        const std::shared_ptr<PanelRadiosityModelSettings>& modelSettings,
        const std::string& sourceBodyName);

/*!
 * Create panel radiosity model updater from a list of panel radiosity model settings.
 *
 * @param modelSettings Settings of the panel radiosity models
 * @param originalSourceToSourceOccultingBodiesMap Occultation map
 * @param sourceBodyName Body to which the panel radiosity models belong
 * @param bodies System of bodies
 * @return Shared pointer to panel radiosity model updater
 */
std::unique_ptr<electromagnetism::SourcePanelRadiosityModelUpdater> createSourcePanelRadiosityModelUpdater(
        const std::vector<std::shared_ptr<PanelRadiosityModelSettings>>& modelSettings,
        const std::map<std::string, std::vector<std::string>>& originalSourceToSourceOccultingBodiesMap,
        const std::string& sourceBodyName,
        const SystemOfBodies& bodies);

/*!
 * Create radiation source model from its settings.
 *
 * @param modelSettings Settings of the radiation source model
 * @param sourceBodyName Body to which the radiation source model belongs
 * @param bodies System of bodies
 * @return Shared pointer to radiation source model
 */
std::shared_ptr<electromagnetism::RadiationSourceModel> createRadiationSourceModel(
        const std::shared_ptr<RadiationSourceModelSettings >& modelSettings,
        const std::string& sourceBodyName,
        const SystemOfBodies& bodies);


} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONSOURCEMODEL_H
