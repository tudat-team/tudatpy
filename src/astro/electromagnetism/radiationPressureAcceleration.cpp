#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"

namespace tudat
{
namespace electromagnetism
{

void RadiationPressureAcceleration::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        currentAcceleration_ = calculateAcceleration();
    }
}

Eigen::Vector3d RadiationPressureAcceleration::calculateAcceleration()
{
    auto sourceCenterPositionInGlobalFrame = sourcePositionFunction_(); // position of center of source (e.g. planet)
    auto sourceRotationFromLocalToGlobalFrame = sourceRotationFromLocalToGlobalFrameFunction_();
    auto sourceRotationFromGlobalToLocalFrame = sourceRotationFromLocalToGlobalFrame.inverse();
    
    auto targetCenterPositionInGlobalFrame = targetPositionFunction_();
    auto targetRotationFromLocalToGlobalFrame = targetRotationFromLocalToGlobalFrameFunction_();
    auto targetRotationFromGlobalToLocalFrame = targetRotationFromLocalToGlobalFrame.inverse();

    // Evaluate irradiances at target position in source frame
    auto targetCenterPositionInSourceFrame =
            sourceRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame);
    auto irradiancesFromSource = sourceModel_->evaluateIrradianceAtPosition(targetCenterPositionInSourceFrame);

    // Calculate radiation pressure force due to all sub-sources
    auto totalForceInTargetFrame = Eigen::Vector3d::Zero().eval();
    for (auto sourceIrradianceAndPosition : irradiancesFromSource) {
        auto sourceIrradiance = std::get<0>(sourceIrradianceAndPosition);
        auto sourcePositionInSourceFrame = std::get<1>(sourceIrradianceAndPosition); // position of sub-source (e.g. panel)
        auto sourcePositionInGlobalFrame = sourceCenterPositionInGlobalFrame + sourceRotationFromLocalToGlobalFrame * sourcePositionInSourceFrame;

        auto sourceToTargetDirectionInGlobalFrame = (targetCenterPositionInGlobalFrame - sourcePositionInGlobalFrame).normalized();
        auto sourceToTargetDirectionInTargetFrame = targetRotationFromGlobalToLocalFrame * sourceToTargetDirectionInGlobalFrame;
        totalForceInTargetFrame += targetModel_->evaluateRadiationPressureForce(sourceIrradiance,
                                                                                sourceToTargetDirectionInTargetFrame);
    }

    auto acceleration = targetRotationFromLocalToGlobalFrame * totalForceInTargetFrame / targetMassFunction_();
    return acceleration;
}

} // tudat
} // electromagnetism
