#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"


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

    switch(modelSettings->getRadiationPressureTargetModelType())
    {
        case cannonball_target:
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
                    cannonballTargetModelSettings->getCoefficient()
            );
            break;
        }
        case paneled_target:
        {
            auto paneledTargetModelSettings =
                    std::dynamic_pointer_cast< PaneledRadiationPressureTargetModelSettings >(modelSettings);

            if(paneledTargetModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected paneled radiation pressure target for body " + body );
            }

            std::vector<PaneledRadiationPressureTargetModel::Panel> panels;
            for (auto& panel : paneledTargetModelSettings->getPanels())
            {
                std::function<Eigen::Vector3d()> surfaceNormalFunction;
                if((panel.getSurfaceNormalFunction() && !panel.getBodyToTrack().empty()) ||
                        (!panel.getSurfaceNormalFunction() && panel.getBodyToTrack().empty()))
                {
                    throw std::runtime_error(
                            "Error, must specify either surface normal or body to track for all"
                            "paneled radiation pressure target panel for body " + body );
                }
                else if (panel.getSurfaceNormalFunction())
                {
                    surfaceNormalFunction = [=] () { return panel.getSurfaceNormalFunction()().normalized(); };
                }
                else {
                    // Tracking a body means setting the surface normal towards the tracked body in the local frame
                    const auto bodyToTrack = bodies.at(panel.getBodyToTrack());
                    const auto targetBody = bodies.at(body);
                    const auto sign = panel.isTowardsTrackedBody() ? +1 : -1;
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

                panels.emplace_back(
                        panel.getArea(),
                        surfaceNormalFunction,
                        SpecularDiffuseMixReflectionLaw::fromSpecularAndDiffuseReflectivity(
                                panel.getSpecularReflectivity(),
                                panel.getDiffuseReflectivity()));
            }

            radiationPressureTargetModel = std::make_shared<PaneledRadiationPressureTargetModel>(panels);
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize radiation pressure target model settings for " + body );
    }

    return radiationPressureTargetModel;
}
} // tudat
} // electromagnetism
