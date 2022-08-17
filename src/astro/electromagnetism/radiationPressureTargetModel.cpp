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
    Eigen::Vector3d force = forceMagnitude * sourceToTargetDirection;
    return force;
}

Eigen::Vector3d PaneledRadiationPressureTargetModel::evaluateRadiationPressureForce(
        double sourceIrradiance,
        Eigen::Vector3d sourceToTargetDirection) const
{
    Eigen::Vector3d force = Eigen::Vector3d::Zero();
    for (auto& panel : panels_)
    {
        const auto surfaceNormal = panel.getSurfaceNormal();
        const double cosBetweenNormalAndIncoming = (-sourceToTargetDirection).dot(surfaceNormal);
        if (cosBetweenNormalAndIncoming >= 0)
        {
            const double effectiveArea = panel.getArea() * cosBetweenNormalAndIncoming;
            const auto radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
            const auto reactionVector =
                    panel.getReflectionLaw()->evaluateReactionVector(surfaceNormal, sourceToTargetDirection);
            force += radiationPressure * effectiveArea * reactionVector;
        }
    }
    return force;
}

void RadiationPressureTargetModel::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

void PaneledRadiationPressureTargetModel::updateMembers_(double currentTime)
{
    for (auto& panel : panels_)
    {
        panel.updateMembers();
    }
}

void PaneledRadiationPressureTargetModel::Panel::updateMembers()
{
    // Evaluate only once per timestep since surface normal function could be expensive to evaluate
    surfaceNormal_ = surfaceNormalFunction_();
}
} // tudat
} // electromagnetism
