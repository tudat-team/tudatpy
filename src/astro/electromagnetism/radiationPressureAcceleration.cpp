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

    Eigen::Vector3d originalSourceCenterPositionInGlobalFrame;
    if (originalSourceModel_ != nullptr)
    {
        originalSourceCenterPositionInGlobalFrame = originalSourcePositionFunction_();
    }
    else
    {
        originalSourceCenterPositionInGlobalFrame = Eigen::Vector3d::Zero();
    }

    // Evaluate original source irradiance in source frame (if applicable, e.g. for albedo)
    const auto originalSourceIrradianceAndPosition = calculateOriginalSourceIrradiance(
            sourceCenterPositionInGlobalFrame, sourceRotationFromGlobalToLocalFrame);
    double originalSourceIrradiance = std::get<0>(originalSourceIrradianceAndPosition);
    Eigen::Vector3d originalSourceToSourceDirectionInSourceFrame = std::get<1>(originalSourceIrradianceAndPosition);

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

        auto occultationFactor = calculateTotalOccultationFactor(
                originalSourceCenterPositionInGlobalFrame, sourcePositionInGlobalFrame, targetCenterPositionInGlobalFrame);
        auto occultedSourceIrradiance = sourceIrradiance * occultationFactor;

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

std::pair<double, Eigen::Vector3d> RadiationPressureAcceleration::calculateOriginalSourceIrradiance(
        const Eigen::Vector3d& sourceCenterPositionInGlobalFrame,
        const Eigen::Quaterniond& sourceRotationFromGlobalToLocalFrame) const
{
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

    return std::make_pair(irradianceFromOriginalSource, originalSourceToSourceDirectionInSourceFrame);
}

double RadiationPressureAcceleration::calculateTotalOccultationFactor(
        const Eigen::Vector3d& originalSourcePosition,
        const Eigen::Vector3d& sourcePosition,
        const Eigen::Vector3d& targetPosition) const
{
    double totalOccultationFactor = 1;

    if (originalSourceModel_ != nullptr)
    {
        totalOccultationFactor *= originalSourceOccultationFactorFunction_(originalSourcePosition, sourcePosition);
    }

    totalOccultationFactor *= sourceOccultationFactorFunction_(sourcePosition, targetPosition);

    return totalOccultationFactor;
}

} // tudat
} // electromagnetism
