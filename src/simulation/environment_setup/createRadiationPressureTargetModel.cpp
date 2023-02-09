/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"

#include <map>
#include <memory>
#include <vector>

#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/body.h"


namespace tudat
{
namespace simulation_setup
{

std::shared_ptr<electromagnetism::RadiationPressureTargetModel> createRadiationPressureTargetModel(
        const std::shared_ptr<RadiationPressureTargetModelSettings>& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<electromagnetism::RadiationPressureTargetModel> radiationPressureTargetModel;

    // Validate occulting bodies map: can either have
    //  - multiple entries with source body names as keys
    //  - a single entry with the empty string as key (use same occulting bodies for all sources)
    auto sourceToTargetOccultingBodies = modelSettings->getSourceToTargetOccultingBodies();
    if (sourceToTargetOccultingBodies.count("") > 0 && sourceToTargetOccultingBodies.size() > 1)
    {
        throw std::runtime_error("Error, invalid occulting bodies map for " + body );
    }

    switch(modelSettings->getRadiationPressureTargetModelType())
    {
        case RadiationPressureTargetModelType::cannonball_target:
        {
            auto cannonballTargetModelSettings =
                    std::dynamic_pointer_cast< CannonballRadiationPressureTargetModelSettings >(modelSettings);

            if(cannonballTargetModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected cannonball radiation pressure target for body " + body );
            }

            radiationPressureTargetModel = std::make_shared<CannonballRadiationPressureTargetModel>(
                    cannonballTargetModelSettings->getArea(),
                    cannonballTargetModelSettings->getCoefficient(),
                    sourceToTargetOccultingBodies);
            break;
        }
        case RadiationPressureTargetModelType::paneled_target:
        {
            auto paneledTargetModelSettings =
                    std::dynamic_pointer_cast< PaneledRadiationPressureTargetModelSettings >(modelSettings);

            if(paneledTargetModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected paneled radiation pressure target for body " + body );
            }

            // Create panels from panel settings
            std::vector<PaneledRadiationPressureTargetModel::Panel> panels;
            for (auto& panel : paneledTargetModelSettings->getPanels())
            {
                std::function<Eigen::Vector3d()> surfaceNormalFunction;
                if((panel.getSurfaceNormalFunction() && !panel.getBodyToTrack().empty()) || // both given
                        (!panel.getSurfaceNormalFunction() && panel.getBodyToTrack().empty())) // none given
                {
                    throw std::runtime_error(
                            "Error, must specify either surface normal or body to track for all "
                            "paneled radiation pressure target panels for body " + body );
                }
                else if (panel.getSurfaceNormalFunction())
                {
                    surfaceNormalFunction = [=] () { return panel.getSurfaceNormalFunction()().normalized(); };
                }
                else {
                    // Tracking a body means setting the surface normal towards the tracked body in the source local frame
                    const auto bodyToTrack = bodies.at(panel.getBodyToTrack());
                    const auto targetBody = bodies.at(body);
                    const auto sign = panel.isTowardsTrackedBody() ? +1 : -1;

                    // Construct surface normal function always pointing towards/away from tracked body
                    surfaceNormalFunction = [=] () {
                        const Eigen::Quaterniond rotationFromPropagationToLocalFrame =
                                targetBody->getCurrentRotationToLocalFrame();
                        const Eigen::Vector3d relativeSourcePositionInPropagationFrame =
                                bodyToTrack->getPosition() - targetBody->getPosition();
                        const Eigen::Vector3d relativeSourcePositionInLocalFrame =
                                rotationFromPropagationToLocalFrame * relativeSourcePositionInPropagationFrame;
                        Eigen::Vector3d surfaceNormal = sign * relativeSourcePositionInLocalFrame.normalized();
                        return surfaceNormal;
                    };
                }

                // Create panel with specular/diffuse-mix reflection law
                panels.emplace_back(
                        panel.getArea(),
                        surfaceNormalFunction,
                        reflectionLawFromSpecularAndDiffuseReflectivity(
                                panel.getSpecularReflectivity(),
                                panel.getDiffuseReflectivity(),
                                panel.isWithInstantaneousReradiation()));
            }

            radiationPressureTargetModel = std::make_shared<PaneledRadiationPressureTargetModel>(
                panels, sourceToTargetOccultingBodies);
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize radiation pressure target model settings for " + body );
    }

    return radiationPressureTargetModel;
}

} // tudat
} // electromagnetism
