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
        // TODO for dynamic paneling, set target position here
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

    double irradianceFromOriginalSource;
    Eigen::Vector3d originalSourceToSourceDirectionInSourceFrame;

    if (originalSourceModel_ != nullptr)
    {
        // Evaluate irradiances from original source at source position in original source frame (for albedo-reflected radiation)
        auto originalSourceCenterPositionInGlobalFrame = originalSourcePositionFunction_();
        auto originalSourceRotationFromLocalToGlobalFrame = originalSourceRotationFromLocalToGlobalFrameFunction_();
        auto originalSourceRotationFromGlobalToLocalFrame = originalSourceRotationFromLocalToGlobalFrame.inverse();

        auto sourceCenterPositionInOriginalSourceFrame =
                originalSourceRotationFromGlobalToLocalFrame * (sourceCenterPositionInGlobalFrame - originalSourceCenterPositionInGlobalFrame);
        irradianceFromOriginalSource = originalSourceModel_->evaluateIrradianceAtPosition(sourceCenterPositionInOriginalSourceFrame);

        originalSourceToSourceDirectionInSourceFrame =
                sourceRotationFromGlobalToLocalFrame * (sourceCenterPositionInGlobalFrame - originalSourceCenterPositionInGlobalFrame).normalized();
    }
    else
    {
        irradianceFromOriginalSource = 0;
        originalSourceToSourceDirectionInSourceFrame = Eigen::Vector3d::Zero();
    }

    // Evaluate irradiances at target position in source frame
    auto targetCenterPositionInSourceFrame =
            sourceRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame);
    auto irradiancesFromSource = sourceModel_->evaluateIrradianceAtPosition(
            targetCenterPositionInSourceFrame,
            irradianceFromOriginalSource,
            originalSourceToSourceDirectionInSourceFrame);

    // Calculate radiation pressure force due to all sub-sources
    auto totalForceInTargetFrame = Eigen::Vector3d::Zero().eval();
    for (auto sourceIrradianceAndPosition : irradiancesFromSource) {
        auto sourceIrradiance = std::get<0>(sourceIrradianceAndPosition);
        auto sourcePositionInSourceFrame = std::get<1>(sourceIrradianceAndPosition); // position of sub-source (e.g. panel)
        auto sourcePositionInGlobalFrame = sourceCenterPositionInGlobalFrame + sourceRotationFromLocalToGlobalFrame * sourcePositionInSourceFrame;

        auto sourceToTargetDirectionInTargetFrame =
                targetRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourcePositionInGlobalFrame).normalized();
        totalForceInTargetFrame += targetModel_->evaluateRadiationPressureForce(sourceIrradiance,
                                                                                sourceToTargetDirectionInTargetFrame);
    }

    auto acceleration = targetRotationFromLocalToGlobalFrame * totalForceInTargetFrame / targetMassFunction_();
    return acceleration;
}

} // tudat
} // electromagnetism
