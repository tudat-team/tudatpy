#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"

namespace tudat
{
namespace electromagnetism
{

Eigen::Vector3d CannonballRadiationPressureTargetModel::evaluateRadiationPressureForce(double sourceIrradiance,
                                                                                       Eigen::Vector3d sourceToTargetDirection) const
{
    const auto radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
    const auto forceMagnitude = coefficient_ * area_ * radiationPressure;
    const auto force = forceMagnitude * sourceToTargetDirection;
    return force;
}

Eigen::Vector3d PaneledRadiationPressureTargetModel::evaluateRadiationPressureForce(
        double sourceIrradiance,
        Eigen::Vector3d sourceToTargetDirection) const
{
    auto force = Eigen::Vector3d::Zero().eval();
    for (auto& panel : panels_)
    {
        const auto surfaceNormal = panel.getSurfaceNormalFunction()();
        // No need to check if panel is illuminated because reaction vector is zero if not
        const auto effectiveArea = panel.getArea() * -sourceToTargetDirection.dot(surfaceNormal);
        const auto radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
        const auto reactionVector =
                panel.getReflectionLaw()->evaluateReactionVector(surfaceNormal, sourceToTargetDirection);
        force += radiationPressure * effectiveArea * reactionVector;
    }
    return force;
}
} // tudat
} // electromagnetism
