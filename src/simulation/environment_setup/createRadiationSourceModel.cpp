/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createRadiationSourceModel.h"

namespace tudat
{
namespace simulation_setup
{

std::shared_ptr<electromagnetism::LuminosityModel> createLuminosityModel(
        const std::shared_ptr<LuminosityModelSettings>& modelSettings,
        const std::string &body)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<LuminosityModel> luminosityModel;

    switch(modelSettings->getLuminosityModelType())
    {
        case constant_radiant_power:
        {
            auto constantLuminosityModelSettings =
                    std::dynamic_pointer_cast< ConstantLuminosityModelSettings >(modelSettings);
            if(constantLuminosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected constant luminosity model for body " + body );
            }

            luminosityModel = std::make_shared<ConstantLuminosityModel>(
                constantLuminosityModelSettings->getLuminosity()
            );
            break;
        }
        case irradiance_based_radiant_power:
        {
            auto irradianceBasedLuminosityModelSettings =
                    std::dynamic_pointer_cast< IrradianceBasedLuminosityModelSettings >(modelSettings);
            if(irradianceBasedLuminosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected irradiance-based luminosity model for body " + body );
            }

            luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
                irradianceBasedLuminosityModelSettings->getIrradianceAtDistanceFunction(),
                irradianceBasedLuminosityModelSettings->getDistance()
            );
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize luminosity model settings for " + body );
    }

    return luminosityModel;
}

std::function<std::shared_ptr<electromagnetism::PaneledRadiationSourceModel::PanelRadiosityModel>(double, double)> createPanelRadiosityModelFunction(
        const std::shared_ptr<PanelRadiosityModelSettings>& modelSettings,
        const std::string& body)
{
    using namespace tudat::electromagnetism;

    std::function<std::shared_ptr<PaneledRadiationSourceModel::PanelRadiosityModel>(double, double)> panelRadiosityModelFunction;

    switch(modelSettings->getPanelRadiosityModelType())
    {
        case albedo:
        {
            auto albedoPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< AlbedoPanelRadiosityModelSettings >(modelSettings);
            if(albedoPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected albedo panel radiosity model for body " + body );
            }

            // Create function returning albedo radiosity model for panel at given polar and azimuth angle
            // Radiosity model uses Lambertian reflectance with given albedo as diffuse reflectivity coefficient
            panelRadiosityModelFunction = [=] (double polarAngle, double azimuthAngle) {
                const auto albedoDistribution = albedoPanelRadiosityModelSettings->getAlbedoDistribution();
                const auto albedo = albedoDistribution(polarAngle, azimuthAngle);
                return std::make_shared<AlbedoPanelRadiosityModel>(
                        std::make_shared<LambertianReflectionLaw>(
                            albedo,
                            albedoPanelRadiosityModelSettings->getWithInstantaneousReradiation()));
            };
            break;
        }
        case thermal_delayed:
        {
            auto delayedThermalPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< DelayedThermalPanelRadiosityModelSettings >(modelSettings);
            if(delayedThermalPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected delayed thermal panel radiosity model for body " + body );
            }

            // Create function returning delayed thermal radiosity model for panel at given polar and azimuth angle
            panelRadiosityModelFunction = [=] (double polarAngle, double azimuthAngle) {
                const auto emissivityDistribution = delayedThermalPanelRadiosityModelSettings->getEmissivityDistribution();
                const auto emissivity = emissivityDistribution(polarAngle, azimuthAngle);
                return std::make_shared<DelayedThermalPanelRadiosityModel>(emissivity);
            };
            break;
        }
        case thermal_angle_based:
        {
            auto angleBasedThermalPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< AngleBasedThermalPanelRadiosityModelSettings >(modelSettings);
            if(angleBasedThermalPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected angle-based thermal panel radiosity model for body " + body );
            }

            // Create function returning angle-based radiosity model for panel at given polar and azimuth angle
            panelRadiosityModelFunction = [=] (double polarAngle, double azimuthAngle) {
                const auto emissivityDistribution = angleBasedThermalPanelRadiosityModelSettings->getEmissivityDistribution();
                const auto emissivity = emissivityDistribution(polarAngle, azimuthAngle);
                return std::make_shared<AngleBasedThermalPanelRadiosityModel>(
                        angleBasedThermalPanelRadiosityModelSettings->getMinTemperature(),
                        angleBasedThermalPanelRadiosityModelSettings->getMaxTemperature(),
                        emissivity);
            };
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize panel radiosity model settings for " + body );
    }

    return panelRadiosityModelFunction;
}

std::shared_ptr<electromagnetism::RadiationSourceModel> createRadiationSourceModel(
        const std::shared_ptr<RadiationSourceModelSettings>& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<RadiationSourceModel> radiationSourceModel;

    switch(modelSettings->getRadiationSourceModelType())
    {
    case isotropic_point_source:
    {
        auto isotropicPointModelSettings =
                std::dynamic_pointer_cast< IsotropicPointRadiationSourceModelSettings >(modelSettings);

        if(isotropicPointModelSettings == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected isotropic point radiation source for body " + body );
        }
        if(isotropicPointModelSettings->getLuminosityModelSettings() == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected isotropic point radiation source to have a luminosity model for body " + body);
        }

        auto luminosityModel = createLuminosityModel(
                isotropicPointModelSettings->getLuminosityModelSettings(), body);

        radiationSourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
        break;
    }
    case statically_paneled_source:
    {
        auto paneledModelSettings =
                std::dynamic_pointer_cast< StaticallyPaneledRadiationSourceModelSettings >(modelSettings);

        if(paneledModelSettings == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected statically paneled radiation source for body " + body );
        }
        if( paneledModelSettings->getOriginalSourceName().empty() )
        {
            throw std::runtime_error(
                    "Error, expected statically paneled radiation source to have an original source for body " + body);
        }
        if(paneledModelSettings->getNumberOfPanels() == 0)
        {
            throw std::runtime_error(
                    "Error, expected statically paneled radiation source to have at least one panel for body " + body);
        }
        if(paneledModelSettings->getPanelRadiosityModelSettings().empty())
        {
            throw std::runtime_error(
                    "Error, expected statically paneled radiation source to have at least one panel radiosity model for body " + body);
        }

        auto sourceBody = bodies.getBody(body);

        // Create list of functions, together returning all panel radiosity models at a given point on body
        std::vector<std::function<std::shared_ptr<PaneledRadiationSourceModel::PanelRadiosityModel>(double, double)>> radiosityModelFunctions;
        for (auto& radiosityModelSetting : paneledModelSettings->getPanelRadiosityModelSettings())
        {
            radiosityModelFunctions.push_back(createPanelRadiosityModelFunction(radiosityModelSetting, body));
        }

        radiationSourceModel = std::make_shared<StaticallyPaneledRadiationSourceModel>(
                paneledModelSettings->getOriginalSourceName(),
                sourceBody->getShapeModel(),
                radiosityModelFunctions,
                paneledModelSettings->getNumberOfPanels(),
                paneledModelSettings->getOriginalSourceToSourceOccultingBodies());
        break;
    }
    default:
        throw std::runtime_error( "Error, do not recognize radiation source model settings for " + body );
    }

    return radiationSourceModel;
}

} // tudat
} // electromagnetism
