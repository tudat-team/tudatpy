#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"

namespace tudat
{
namespace electromagnetism
{

void RadiationPressureAcceleration::updateMembers(const double currentTime)
{
    if(this->currentTime_ != currentTime)
    {
        this->currentAcceleration_ = calculateAcceleration();
        this->currentTime_ = currentTime;
    }
}

Eigen::Vector3d RadiationPressureAcceleration::calculateAcceleration()
{
    // Radiation pressure force is evaluated in local frame, then rotated to propagation frame
    auto targetPosition = targetPositionFunction_();
    auto targetRotationFromLocalToPropagationFrame = targetRotationFromLocalToPropagationFrameFunction_();
    auto targetRotationFromPropagationToLocalFrame = targetRotationFromLocalToPropagationFrame.inverse();

    auto totalForce = Eigen::Vector3d::Zero().eval();
    auto irradiancesFromSource = sourceModel_->evaluateIrradianceAtPosition(targetPosition);
    for (auto sourceIrradianceAndPosition : irradiancesFromSource) {
        auto sourceIrradiance = std::get<0>(sourceIrradianceAndPosition);
        auto sourcePosition = std::get<1>(sourceIrradianceAndPosition);

        auto sourceToTargetDirection = (targetPosition - sourcePosition).normalized();
        auto sourceToTargetDirectionInLocalFrame = targetRotationFromPropagationToLocalFrame * sourceToTargetDirection;
        totalForce += targetModel_->evaluateRadiationPressureForce(sourceIrradiance,
                                                                   sourceToTargetDirectionInLocalFrame);
    }

    auto acceleration = targetRotationFromLocalToPropagationFrame * totalForce / targetMassFunction_();
    return acceleration;
}


} // tudat
} // electromagnetism
