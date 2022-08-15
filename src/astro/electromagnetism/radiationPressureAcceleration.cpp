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
    Eigen::Vector3d sourceCenterPositionInGlobalFrame = sourcePositionFunction_(); // position of center of source (e.g. planet)
    Eigen::Quaterniond sourceRotationFromLocalToGlobalFrame = sourceRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond sourceRotationFromGlobalToLocalFrame = sourceRotationFromLocalToGlobalFrame.inverse();

    Eigen::Vector3d targetCenterPositionInGlobalFrame = targetPositionFunction_();
    Eigen::Quaterniond targetRotationFromLocalToGlobalFrame = targetRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond targetRotationFromGlobalToLocalFrame = targetRotationFromLocalToGlobalFrame.inverse();

    double irradianceFromOriginalSource;
    Eigen::Vector3d originalSourceToSourceDirectionInSourceFrame;

    if (originalSourceModel_ != nullptr)
    {
        // Evaluate irradiances from original source at source position in original source frame (for albedo-reflected radiation)
        Eigen::Vector3d originalSourceCenterPositionInGlobalFrame = originalSourcePositionFunction_();
        Eigen::Quaterniond originalSourceRotationFromLocalToGlobalFrame = originalSourceRotationFromLocalToGlobalFrameFunction_();
        Eigen::Quaterniond originalSourceRotationFromGlobalToLocalFrame = originalSourceRotationFromLocalToGlobalFrame.inverse();

        Eigen::Vector3d sourceCenterPositionInOriginalSourceFrame =
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
    Eigen::Vector3d targetCenterPositionInSourceFrame =
            sourceRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame);
    auto irradiancesFromSource = sourceModel_->evaluateIrradianceAtPosition(
            targetCenterPositionInSourceFrame,
            irradianceFromOriginalSource,
            originalSourceToSourceDirectionInSourceFrame);

    // Calculate radiation pressure force due to all sub-sources
    Eigen::Vector3d totalForceInTargetFrame = Eigen::Vector3d::Zero();
    for (auto sourceIrradianceAndPosition : irradiancesFromSource) {
        auto sourceIrradiance = std::get<0>(sourceIrradianceAndPosition);
        Eigen::Vector3d sourcePositionInSourceFrame = std::get<1>(sourceIrradianceAndPosition); // position of sub-source (e.g. panel)
        Eigen::Vector3d sourcePositionInGlobalFrame = sourceCenterPositionInGlobalFrame + sourceRotationFromLocalToGlobalFrame * sourcePositionInSourceFrame;

        Eigen::Vector3d sourceToTargetDirectionInTargetFrame =
                targetRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourcePositionInGlobalFrame).normalized();
        totalForceInTargetFrame += targetModel_->evaluateRadiationPressureForce(sourceIrradiance,
                                                                                sourceToTargetDirectionInTargetFrame);
    }

    Eigen::Vector3d acceleration = targetRotationFromLocalToGlobalFrame * totalForceInTargetFrame / targetMassFunction_();
    return acceleration;
}

} // tudat
} // electromagnetism
