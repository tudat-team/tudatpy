/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Ashwin R. Vasavada et al.. "Lunar equatorial surface temperatures and regolith properties from the Diviner Lunar
 *          Radiometer Experiment".Journal of Geophysical Research: Planets 117, no.E12 (2012).
 */

#include "tudat/simulation/environment_setup/createRadiationSourceModel.h"

#include <memory>
#include <set>
#include <string>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/simulation/environment_setup/createOccultationModel.h"
#include "tudat/simulation/environment_setup/createSurfacePropertyDistribution.h"
#include "tudat/astro/electromagnetism/occultationModel.h"


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
        case LuminosityModelType::constant_radiant_power:
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
        case LuminosityModelType::time_variable_isotropic_radiant_power :
        {
            auto variableLuminosityModelSettings =
                std::dynamic_pointer_cast< TimeVariableLuminosityModelSettings >(modelSettings);
            if(variableLuminosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                    "Error, expected time-variable luminosity model for body " + body );
            }

            luminosityModel = std::make_shared<VariableLuminosityModel>(
                variableLuminosityModelSettings->getLuminosityFuntion( )
            );
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize luminosity model settings for " + body );
    }

    return luminosityModel;
}

std::unique_ptr<electromagnetism::SourcePanelRadiosityModel> createPanelRadiosityModel(
        const std::shared_ptr<PanelRadiosityModelSettings>& modelSettings,
        const std::string& sourceBodyName)
{
    using namespace tudat::electromagnetism;

    std::unique_ptr<SourcePanelRadiosityModel> panelRadiosityModel;

    switch(modelSettings->getPanelRadiosityModelType())
    {
        case PanelRadiosityModelType::constant:
        {
            auto constantPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< ConstantPanelRadiosityModelSettings >(modelSettings);
            if(constantPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected constant panel radiosity model for sourceBodyName " + sourceBodyName );
            }

            panelRadiosityModel = std::make_unique<ConstantSourcePanelRadiosityModel>(
                    constantPanelRadiosityModelSettings->getConstantRadiosity());
            break;
        }
        case PanelRadiosityModelType::albedo:
        {
            auto albedoPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< AlbedoPanelRadiosityModelSettings >(modelSettings);
            if(albedoPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected albedo panel radiosity model for sourceBodyName " + sourceBodyName );
            }

            panelRadiosityModel = std::make_unique<AlbedoSourcePanelRadiosityModel>(
                    albedoPanelRadiosityModelSettings->getOriginalSourceName(),
                    createSurfacePropertyDistribution(
                            albedoPanelRadiosityModelSettings->getAlbedoDistribution(), sourceBodyName));
            break;
        }
        case PanelRadiosityModelType::thermal_delayed:
        {
            auto delayedThermalPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< DelayedThermalPanelRadiosityModelSettings >(modelSettings);
            if(delayedThermalPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected delayed thermal panel radiosity model for sourceBodyName " + sourceBodyName );
            }

            panelRadiosityModel = std::make_unique<DelayedThermalSourcePanelRadiosityModel>(
                    delayedThermalPanelRadiosityModelSettings->getOriginalSourceName(),
                    createSurfacePropertyDistribution(
                            delayedThermalPanelRadiosityModelSettings->getEmissivityDistribution(), sourceBodyName));
            break;
        }
        case PanelRadiosityModelType::thermal_angle_based:
        {
            auto angleBasedThermalPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< AngleBasedThermalPanelRadiosityModelSettings >(modelSettings);
            if(angleBasedThermalPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected angle-based thermal panel radiosity model for sourceBodyName " + sourceBodyName );
            }

            panelRadiosityModel = std::make_unique<AngleBasedThermalSourcePanelRadiosityModel>(
                    angleBasedThermalPanelRadiosityModelSettings->getOriginalSourceName(),
                    angleBasedThermalPanelRadiosityModelSettings->getMinTemperature(),
                    angleBasedThermalPanelRadiosityModelSettings->getMaxTemperature(),
                    createSurfacePropertyDistribution(
                            angleBasedThermalPanelRadiosityModelSettings->getEmissivityDistribution(), sourceBodyName));
            break;
        }
        default:
            throw std::runtime_error("Error, do not recognize panel radiosity model settings for " + sourceBodyName );
    }

    return panelRadiosityModel;
}

std::unique_ptr<electromagnetism::SourcePanelRadiosityModelUpdater> createSourcePanelRadiosityModelUpdater(
        const std::vector<std::shared_ptr<PanelRadiosityModelSettings>>& modelSettings,
        const std::map<std::string, std::vector<std::string>>& originalSourceToSourceOccultingBodiesMap,
        const std::string& sourceBodyName,
        const SystemOfBodies& bodies)
{

    auto sourceBody = bodies.at(sourceBodyName);

    // Collect original source body names
    std::set<std::string> originalSourceBodyNames;
    for (auto& radiosityModelSetting : modelSettings)
    {
        auto originalSourceDependentRadiosityModelSetting =
            std::dynamic_pointer_cast<OriginalSourceDependentPanelRadiosityModelSettings>(radiosityModelSetting);
        if (originalSourceDependentRadiosityModelSetting != nullptr)
        {
            originalSourceBodyNames.insert(originalSourceDependentRadiosityModelSetting->getOriginalSourceName());
        }
    }

    // Collect original source properties
    std::map<std::string, std::shared_ptr<electromagnetism::IsotropicPointRadiationSourceModel>> originalSourceModels;
    std::map<std::string, std::shared_ptr<basic_astrodynamics::BodyShapeModel>> originalSourceBodyShapeModels;
    std::map<std::string, std::function<Eigen::Vector3d()>> originalSourcePositionFunctions;
    std::map<std::string, std::shared_ptr<electromagnetism::OccultationModel>> originalSourceToSourceOccultationModels;
    for (auto& originalSourceBodyName : originalSourceBodyNames)
    {
        auto originalSourceBody = bodies.at(originalSourceBodyName);

        if( originalSourceBodyName.empty())
        {
            throw std::runtime_error(
                    "Error, expected panel radiosity model to have an original source for sourceBody " + sourceBodyName);
        }
        if( bodies.count( originalSourceBodyName ) == 0 )
        {
            throw std::runtime_error( "Error, body " + originalSourceBodyName + " (original source) was not found. " +
                                      "This may happen when a default radiation source model is used.");
        }

        auto originalIsotropicPointRadiationSourceModel =
                std::dynamic_pointer_cast<electromagnetism::IsotropicPointRadiationSourceModel>(
                        originalSourceBody->getRadiationSourceModel());

        if( originalIsotropicPointRadiationSourceModel == nullptr )
        {
            throw std::runtime_error( "Error, body " + originalSourceBodyName +
                                      " (original source) has no isotropic point radiation source model." );
        }

        std::vector<std::string> originalSourceToSourceOccultingBodies = {};
        if (originalSourceToSourceOccultingBodiesMap.count(originalSourceBodyName) > 0)
        {
            originalSourceToSourceOccultingBodies = originalSourceToSourceOccultingBodiesMap.at(originalSourceBodyName);
        }
        else if (originalSourceToSourceOccultingBodiesMap.count("") > 0)
        {
            // Use the same occulting bodies for all original sources
            originalSourceToSourceOccultingBodies = originalSourceToSourceOccultingBodiesMap.at("");
        }

        // Check if occulting bodies are not original source or source
        for (auto& occultingBodyName : originalSourceToSourceOccultingBodies)
        {
            if (occultingBodyName == originalSourceBodyName)
            {
                throw std::runtime_error( "Error, original source body cannot act as occulting body.");
            }
            if (occultingBodyName == sourceBodyName)
            {
                throw std::runtime_error( "Error, source body cannot act as occulting body.");
            }
        }
        auto originalSourceToSourceOccultationModel = createOccultationModel(originalSourceToSourceOccultingBodies, bodies);

        originalSourceModels[originalSourceBodyName] = originalIsotropicPointRadiationSourceModel;
        originalSourceBodyShapeModels[originalSourceBodyName] = originalSourceBody->getShapeModel();
        originalSourcePositionFunctions[originalSourceBodyName] = std::bind(&Body::getPosition, originalSourceBody);
        originalSourceToSourceOccultationModels[originalSourceBodyName] = originalSourceToSourceOccultationModel;
    }

    return std::make_unique<electromagnetism::SourcePanelRadiosityModelUpdater>(
            [sourceBody] { return sourceBody->getPosition(); },
            [sourceBody] { return sourceBody->getCurrentRotationToGlobalFrame(); },
            originalSourceModels,
            originalSourceBodyShapeModels,
            originalSourcePositionFunctions,
            originalSourceToSourceOccultationModels);
}

std::shared_ptr<electromagnetism::RadiationSourceModel> createRadiationSourceModel(
        const std::shared_ptr<RadiationSourceModelSettings>& modelSettings,
        const std::string& sourceBodyName,
        const SystemOfBodies& bodies)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<RadiationSourceModel> radiationSourceModel;

    switch(modelSettings->getRadiationSourceModelType())
    {
    case RadiationSourceModelType::isotropic_point_source:
    {
        auto isotropicPointModelSettings =
                std::dynamic_pointer_cast< IsotropicPointRadiationSourceModelSettings >(modelSettings);

        if(isotropicPointModelSettings == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected isotropic point radiation source for sourceBody " + sourceBodyName );
        }
        if(isotropicPointModelSettings->getLuminosityModelSettings() == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected isotropic point radiation source to have a luminosity model for sourceBody " + sourceBodyName);
        }

        auto luminosityModel = createLuminosityModel(
                isotropicPointModelSettings->getLuminosityModelSettings(), sourceBodyName);

        radiationSourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel, sourceBodyName);
        break;
    }
    case RadiationSourceModelType::extended_source:
    {
        auto paneledModelSettings =
                std::dynamic_pointer_cast< ExtendedRadiationSourceModelSettings >(modelSettings);

        if(paneledModelSettings == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected extended radiation source for sourceBody " + sourceBodyName );
        }
        if(paneledModelSettings->getPanelRadiosityModelSettings().empty())
        {
            throw std::runtime_error(
                    "Error, expected extended radiation source to have at least one panel radiosity model for sourceBody " + sourceBodyName);
        }

        auto sourceBody = bodies.getBody(sourceBodyName);

        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels;
        for (auto& radiosityModelSetting : paneledModelSettings->getPanelRadiosityModelSettings())
        {
            radiosityModels.push_back(createPanelRadiosityModel(radiosityModelSetting, sourceBodyName));
        }

        auto sourcePanelRadiosityModelUpdater = createSourcePanelRadiosityModelUpdater(
                paneledModelSettings->getPanelRadiosityModelSettings(),
                paneledModelSettings->getOriginalSourceToSourceOccultingBodies(),
                sourceBodyName,
                bodies);

        radiationSourceModel = std::make_shared<DynamicallyPaneledRadiationSourceModel>(
                sourceBody->getShapeModel(),
                std::move(sourcePanelRadiosityModelUpdater),
                radiosityModels,
                paneledModelSettings->getNumberOfPanelsPerRing(),
                sourceBodyName);
        break;
    }
    default:
        throw std::runtime_error("Error, do not recognize radiation source model settings for " + sourceBodyName );
    }

    return radiationSourceModel;
}

} // tudat
} // electromagnetism
