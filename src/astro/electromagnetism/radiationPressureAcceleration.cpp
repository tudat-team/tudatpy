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

Eigen::Vector3d IsotropicPointSourceRadiationPressureAcceleration::calculateAcceleration()
{
    Eigen::Vector3d sourceCenterPositionInGlobalFrame = sourcePositionFunction_(); // position of center of source (e.g. planet)
    Eigen::Quaterniond sourceRotationFromLocalToGlobalFrame = sourceRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond sourceRotationFromGlobalToLocalFrame = sourceRotationFromLocalToGlobalFrame.inverse();

    Eigen::Vector3d targetCenterPositionInGlobalFrame = targetPositionFunction_();
    Eigen::Quaterniond targetRotationFromLocalToGlobalFrame = targetRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond targetRotationFromGlobalToLocalFrame = targetRotationFromLocalToGlobalFrame.inverse();

    // Evaluate irradiances at target position in source frame
    Eigen::Vector3d targetCenterPositionInSourceFrame =
            sourceRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame);
    auto occultationFactor = occultationModel_->evaluateReceivedFraction(
            sourceCenterPositionInGlobalFrame, sourceBodyShapeModel_, targetCenterPositionInGlobalFrame);
    auto sourceIrradiance = isotropicPointSourceModel_->evaluateIrradianceAtPosition(targetCenterPositionInSourceFrame);
    auto occultedSourceIrradiance = sourceIrradiance * occultationFactor;

    if (occultedSourceIrradiance > 0)
    {
        // No body is occluding source seen from target
        Eigen::Vector3d sourceToTargetDirectionInTargetFrame =
                targetRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame).normalized();
        // Calculate radiation pressure force due to source
        Eigen::Vector3d totalForceInTargetFrame =
                targetModel_->evaluateRadiationPressureForce(sourceIrradiance, sourceToTargetDirectionInTargetFrame);
        // Calculate acceleration due to radiation pressure in global frame
        Eigen::Vector3d acceleration = targetRotationFromLocalToGlobalFrame * totalForceInTargetFrame / targetMassFunction_();
        return acceleration;
    }
    else
    {
        return Eigen::Vector3d::Zero();
    }
}

Eigen::Vector3d PaneledSourceRadiationPressureAcceleration::calculateAcceleration()
{
    Eigen::Vector3d sourceCenterPositionInGlobalFrame = sourcePositionFunction_(); // position of center of source (e.g. planet)
    Eigen::Quaterniond sourceRotationFromLocalToGlobalFrame = sourceRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond sourceRotationFromGlobalToLocalFrame = sourceRotationFromLocalToGlobalFrame.inverse();

    Eigen::Vector3d targetCenterPositionInGlobalFrame = targetPositionFunction_();
    Eigen::Quaterniond targetRotationFromLocalToGlobalFrame = targetRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond targetRotationFromGlobalToLocalFrame = targetRotationFromLocalToGlobalFrame.inverse();

    Eigen::Vector3d originalSourceCenterPositionInGlobalFrame = originalSourcePositionFunction_();
    Eigen::Quaterniond originalSourceRotationFromLocalToGlobalFrame = originalSourceRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond originalSourceRotationFromGlobalToLocalFrame = originalSourceRotationFromLocalToGlobalFrame.inverse();

    // Evaluate irradiances from original source at source position in original source frame (for albedo-reflected radiation)
    Eigen::Vector3d sourceCenterPositionInOriginalSourceFrame =
            originalSourceRotationFromGlobalToLocalFrame * (sourceCenterPositionInGlobalFrame - originalSourceCenterPositionInGlobalFrame);
    auto originalSourceIrradiance = originalSourceModel_->evaluateIrradianceAtPosition(sourceCenterPositionInOriginalSourceFrame);
    Eigen::Vector3d originalSourceToSourceDirectionInSourceFrame =
            sourceRotationFromGlobalToLocalFrame * (sourceCenterPositionInGlobalFrame - originalSourceCenterPositionInGlobalFrame).normalized();

    // Evaluate irradiances at target position in source frame
    Eigen::Vector3d targetCenterPositionInSourceFrame =
            sourceRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame);
    auto sourceIrradiancesAndPositions = sourceModel_->evaluateIrradianceAtPosition(
            targetCenterPositionInSourceFrame,
            originalSourceIrradiance,
            originalSourceToSourceDirectionInSourceFrame);

    // Calculate radiation pressure force due to all sub-sources in target frame
    Eigen::Vector3d totalForceInTargetFrame = Eigen::Vector3d::Zero();
    for (auto sourceIrradianceAndPosition : sourceIrradiancesAndPositions) {
        auto sourceIrradiance = std::get<0>(sourceIrradianceAndPosition);
        Eigen::Vector3d sourcePositionInSourceFrame =
                std::get<1>(sourceIrradianceAndPosition); // position of sub-source (e.g. panel)
        Eigen::Vector3d sourcePositionInGlobalFrame =
                sourceCenterPositionInGlobalFrame + sourceRotationFromLocalToGlobalFrame * sourcePositionInSourceFrame;

        auto originalSourceToSourceOccultationFactor = occultationModel_->evaluateReceivedFraction(
                originalSourceCenterPositionInGlobalFrame, originalSourceBodyShapeModel_, targetCenterPositionInGlobalFrame);
        auto sourceToTargetOccultationFactor = static_cast<double>(occultationModel_->evaluateVisibility(
                sourcePositionInSourceFrame, targetCenterPositionInGlobalFrame));
        auto occultedSourceIrradiance =
                sourceIrradiance * originalSourceToSourceOccultationFactor * sourceToTargetOccultationFactor;

        if (occultedSourceIrradiance > 0)
        {
            // No body is occluding source seen from target
            Eigen::Vector3d sourceToTargetDirectionInTargetFrame =
                    targetRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourcePositionInGlobalFrame).normalized();
            totalForceInTargetFrame += targetModel_->evaluateRadiationPressureForce(occultedSourceIrradiance,
                                                                                    sourceToTargetDirectionInTargetFrame);
        }
    }

    // Calculate acceleration due to radiation pressure in global frame
    Eigen::Vector3d acceleration = targetRotationFromLocalToGlobalFrame * totalForceInTargetFrame / targetMassFunction_();
    return acceleration;
}
} // tudat
} // electromagnetism
