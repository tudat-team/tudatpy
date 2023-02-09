/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"

#include <functional>

#include <Eigen/Core>
#include <Eigen/Geometry>


namespace tudat
{
namespace electromagnetism
{

void RadiationPressureAcceleration::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;

        sourceToTargetOccultationModel_->updateMembers(currentTime);
        updateMembers_(currentTime);

        currentAcceleration_ = calculateAcceleration();
    }
}

Eigen::Vector3d IsotropicPointSourceRadiationPressureAcceleration::calculateAcceleration()
{
    Eigen::Vector3d sourceCenterPositionInGlobalFrame = sourcePositionFunction_(); // position of center of source (e.g. planet)

    Eigen::Vector3d targetCenterPositionInGlobalFrame = targetPositionFunction_();
    Eigen::Quaterniond targetRotationFromLocalToGlobalFrame = targetRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond targetRotationFromGlobalToLocalFrame = targetRotationFromLocalToGlobalFrame.inverse();

    // Evaluate irradiances at target position in source frame
    // No rotation to source frame is necessary because isotropic sources are rotation-invariant
    Eigen::Vector3d targetCenterPositionInSourceFrame = targetCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame;
    sourceToTargetReceivedFraction = sourceToTargetOccultationModel_->evaluateReceivedFractionFromExtendedSource(
            sourceCenterPositionInGlobalFrame, sourceBodyShapeModel_, targetCenterPositionInGlobalFrame);
    auto sourceIrradiance = sourceModel_->evaluateIrradianceAtPosition(targetCenterPositionInSourceFrame);
    auto occultedSourceIrradiance = sourceIrradiance * sourceToTargetReceivedFraction;

    // Update dependent variable
    receivedIrradiance = occultedSourceIrradiance;

    if (occultedSourceIrradiance <= 0)
    {
        // Some body is occluding source seen from target
        return Eigen::Vector3d::Zero();
    }

    // Acceleration points from source to target
    Eigen::Vector3d sourceToTargetDirectionInTargetFrame =
            targetRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame).normalized();
    // Calculate radiation pressure force due to source
    Eigen::Vector3d totalForceInTargetFrame =
            targetModel_->evaluateRadiationPressureForce(sourceIrradiance, sourceToTargetDirectionInTargetFrame);
    // Calculate acceleration due to radiation pressure in global frame
    Eigen::Vector3d acceleration = targetRotationFromLocalToGlobalFrame * totalForceInTargetFrame / targetMassFunction_();
    return acceleration;
}

Eigen::Vector3d PaneledSourceRadiationPressureAcceleration::calculateAcceleration()
{
    // Could use class member to avoid allocation every call, but profiling shows allocation is by far
    // dominated by algebraic operations
    Eigen::Vector3d sourceCenterPositionInGlobalFrame = sourcePositionFunction_(); // position of center of source (e.g. planet)
    Eigen::Quaterniond sourceRotationFromLocalToGlobalFrame = sourceRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond sourceRotationFromGlobalToLocalFrame = sourceRotationFromLocalToGlobalFrame.inverse();

    Eigen::Vector3d targetCenterPositionInGlobalFrame = targetPositionFunction_();
    Eigen::Quaterniond targetRotationFromLocalToGlobalFrame = targetRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond targetRotationFromGlobalToLocalFrame = targetRotationFromLocalToGlobalFrame.inverse();

    Eigen::Vector3d originalSourceCenterPositionInGlobalFrame = originalSourcePositionFunction_();

    // Evaluate irradiances from original source at source position in original source frame (for albedo-reflected radiation)
    // If other types than isotropic point sources are supported as original source, rotate to original source frame here
    Eigen::Vector3d sourceCenterPositionInOriginalSourceFrame = sourceCenterPositionInGlobalFrame - originalSourceCenterPositionInGlobalFrame;
    auto originalSourceIrradiance = originalSourceModel_->evaluateIrradianceAtPosition(sourceCenterPositionInOriginalSourceFrame);
    Eigen::Vector3d originalSourceToSourceDirectionInSourceFrame =
            sourceRotationFromGlobalToLocalFrame * (sourceCenterPositionInGlobalFrame - originalSourceCenterPositionInGlobalFrame).normalized();

    // Evaluate irradiances from all sub-sources at target position in source frame
    Eigen::Vector3d targetCenterPositionInSourceFrame =
            sourceRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame);
    auto sourceIrradiancesAndPositions = sourceModel_->evaluateIrradianceAtPosition(
            targetCenterPositionInSourceFrame,
            originalSourceIrradiance,
            originalSourceToSourceDirectionInSourceFrame);

    // For dependent variables
    double totalReceivedIrradiance = 0;
    unsigned int visibleSourcePanelCounter = 0;
    unsigned int illuminatedSourcePanelCounter = 0;
    unsigned int visibleAndIlluminatedSourcePanelCounter = 0;

    // Calculate radiation pressure force due to all sub-sources in target frame
    Eigen::Vector3d totalForceInTargetFrame = Eigen::Vector3d::Zero();
    for (auto sourceIrradianceAndPosition : sourceIrradiancesAndPositions) {
        auto sourceIrradiance = std::get<0>(sourceIrradianceAndPosition);
        Eigen::Vector3d sourcePositionInSourceFrame =
                std::get<1>(sourceIrradianceAndPosition); // position of sub-source (e.g. panel)
        Eigen::Vector3d sourcePositionInGlobalFrame =
                sourceCenterPositionInGlobalFrame + sourceRotationFromLocalToGlobalFrame * sourcePositionInSourceFrame;

        auto originalSourceToSourceReceivedFraction =
                originalSourceToSourceOccultationModel_->evaluateReceivedFractionFromExtendedSource(
                originalSourceCenterPositionInGlobalFrame, originalSourceBodyShapeModel_, sourcePositionInGlobalFrame);
        auto sourceToTargetReceivedFraction =
                sourceToTargetOccultationModel_->evaluateReceivedFractionFromPointSource(sourcePositionInGlobalFrame,
                                                                                         targetCenterPositionInGlobalFrame);
        auto occultedSourceIrradiance =
                sourceIrradiance * originalSourceToSourceReceivedFraction * sourceToTargetReceivedFraction;

        if (occultedSourceIrradiance > 0)
        {
            // No body is occluding source seen from target
            Eigen::Vector3d sourceToTargetDirectionInTargetFrame =
                    targetRotationFromGlobalToLocalFrame * (targetCenterPositionInGlobalFrame - sourcePositionInGlobalFrame).normalized();
            totalForceInTargetFrame +=
                    targetModel_->evaluateRadiationPressureForce(occultedSourceIrradiance, sourceToTargetDirectionInTargetFrame);
        }

        totalReceivedIrradiance += occultedSourceIrradiance;
        if (sourceToTargetReceivedFraction > 0)
        {
            visibleSourcePanelCounter += 1;
        }
        if (originalSourceToSourceReceivedFraction > 0)
        {
            illuminatedSourcePanelCounter += 1;
        }
        if (sourceToTargetReceivedFraction > 0 && originalSourceToSourceReceivedFraction > 0)
        {
            visibleAndIlluminatedSourcePanelCounter += 1;
        }
    }

    // Update dependent variables
    receivedIrradiance = totalReceivedIrradiance;
    visibleSourcePanelCount = visibleSourcePanelCounter;
    illuminatedSourcePanelCount = illuminatedSourcePanelCounter;
    visibleAndIlluminatedSourcePanelCount = visibleAndIlluminatedSourcePanelCounter;

    // Calculate acceleration due to radiation pressure in global frame
    Eigen::Vector3d acceleration = targetRotationFromLocalToGlobalFrame * totalForceInTargetFrame / targetMassFunction_();
    return acceleration;
}

void PaneledSourceRadiationPressureAcceleration::updateMembers_(double currentTime)
{
    originalSourceToSourceOccultationModel_->updateMembers(currentTime);
}

} // tudat
} // electromagnetism
