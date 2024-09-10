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
        computeAcceleration( );
    }
}

void IsotropicPointSourceRadiationPressureAcceleration::computeAcceleration()
{
    sourceCenterPositionInGlobalFrame_ = sourcePositionFunction_();
    targetCenterPositionInGlobalFrame_ = targetPositionFunction_();
    targetCenterPositionInSourceFrame_ = targetCenterPositionInGlobalFrame_ - sourceCenterPositionInGlobalFrame_;
    currentTargetMass_ = targetMassFunction_();
    // Evaluate irradiances at target position in source frame
    // No rotation to source frame is necessary because isotropic sources are rotation-invariant
    sourceToTargetReceivedFraction = sourceToTargetOccultationModel_->evaluateReceivedFractionFromExtendedSource(
            sourceCenterPositionInGlobalFrame_, sourceBodyShapeModel_, targetCenterPositionInGlobalFrame_ );

    receivedIrradiance =
        sourceModel_->evaluateIrradianceAtPosition( targetCenterPositionInSourceFrame_).front().first * sourceToTargetReceivedFraction;

    if (receivedIrradiance <= 0)
    {
        // Some body is occluding source as seen from target
        currentUnscaledAcceleration_ = Eigen::Vector3d::Zero();
        targetModel_->resetComputations( sourceName_ );
    }
    else
    {

        if ( targetModel_->forceFunctionRequiresLocalFrameInputs( ) )
        {
            targetRotationFromLocalToGlobalFrame_ = targetRotationFromLocalToGlobalFrameFunction_( );
            targetRotationFromGlobalToLocalFrame_ = targetRotationFromLocalToGlobalFrame_.inverse( );

            // Calculate acceleration due to radiation pressure in global frame
            targetModel_->updateRadiationPressureForcing(
                receivedIrradiance, targetRotationFromGlobalToLocalFrame_ *
                                    targetCenterPositionInSourceFrame_.normalized( ), true, sourceName_ );
            targetModel_->saveLocalComputations( sourceName_, true );
            currentUnscaledAcceleration_ = targetRotationFromLocalToGlobalFrame_ *
                                   targetModel_->getCurrentRadiationPressureForce( sourceName_ ) /
                                   currentTargetMass_;
        }
        else
        {
            targetModel_->updateRadiationPressureForcing(
                receivedIrradiance, targetCenterPositionInSourceFrame_.normalized( ), true, sourceName_ );
            currentUnscaledAcceleration_ = targetModel_->getCurrentRadiationPressureForce( sourceName_ ) / currentTargetMass_;
        }
    }
    scaleRadiationPressureAcceleration( );
}

void PaneledSourceRadiationPressureAcceleration::computeAcceleration()
{
    // Could use class member to avoid allocation every call, but profiling shows allocation is by far
    // dominated by algebraic operations
    sourceCenterPositionInGlobalFrame_ = sourcePositionFunction_();
    targetCenterPositionInGlobalFrame_ = targetPositionFunction_();
    targetCenterPositionInSourceFrame_ = targetCenterPositionInGlobalFrame_ - sourceCenterPositionInGlobalFrame_;

    Eigen::Quaterniond sourceRotationFromLocalToGlobalFrame = sourceRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond sourceRotationFromGlobalToLocalFrame = sourceRotationFromLocalToGlobalFrame.inverse();

    Eigen::Quaterniond targetRotationFromLocalToGlobalFrame = targetRotationFromLocalToGlobalFrameFunction_();
    Eigen::Quaterniond targetRotationFromGlobalToLocalFrame = targetRotationFromLocalToGlobalFrame.inverse();

    // Evaluate irradiances from all sub-sources at target position in source frame
    Eigen::Vector3d targetCenterPositionInSourceFrame =
            sourceRotationFromGlobalToLocalFrame * targetCenterPositionInSourceFrame_;
    auto sourceIrradiancesAndPositions = sourceModel_->evaluateIrradianceAtPosition(targetCenterPositionInSourceFrame);

    // For dependent variables
    double totalReceivedIrradiance = 0;
    unsigned int visibleAndEmittingSourcePanelCounter = 0;

    // Calculate radiation pressure force due to all sub-sources in target frame
    Eigen::Vector3d totalForceInTargetFrame = Eigen::Vector3d::Zero();
    targetModel_->resetComputations( sourceName_ );
    int counter = 0;
    for (auto sourceIrradianceAndPosition : sourceIrradiancesAndPositions)
    {
        double sourceIrradiance = std::get<0>(sourceIrradianceAndPosition);
        Eigen::Vector3d sourcePositionInSourceFrame =
                std::get<1>(sourceIrradianceAndPosition); // position of sub-source (e.g. panel)
        Eigen::Vector3d sourcePositionInGlobalFrame =
                sourceCenterPositionInGlobalFrame_ + sourceRotationFromLocalToGlobalFrame * sourcePositionInSourceFrame;

        auto sourceToTargetReceivedFraction =
                sourceToTargetOccultationModel_->evaluateReceivedFractionFromPointSource(sourcePositionInGlobalFrame,
                                                                                         targetCenterPositionInGlobalFrame_);
        auto occultedSourceIrradiance =
                sourceIrradiance * sourceToTargetReceivedFraction;

        if (occultedSourceIrradiance > 0)
        {
            // No body is occluding source as seen from target
            Eigen::Vector3d sourceToTargetDirectionInTargetFrame =
                targetRotationFromGlobalToLocalFrame * ( targetCenterPositionInGlobalFrame_ - sourcePositionInGlobalFrame ).normalized();

            targetModel_->updateRadiationPressureForcing( occultedSourceIrradiance, sourceToTargetDirectionInTargetFrame, false, sourceName_ );

            totalReceivedIrradiance += occultedSourceIrradiance;
            visibleAndEmittingSourcePanelCounter += 1;
        }
        if( savePanellingIrradiance_ )
        {
            savedPanelIrradiances_[ counter ] = sourceIrradiance;
        }
        counter++;
    }
    targetModel_->saveLocalComputations( sourceName_, false );
    if( savePanellingGeometry_ )
    {
        savedPanelGeometries_ = sourceModel_->getCurrentPanelGeomtry( );
    }
    totalForceInTargetFrame = targetModel_->getCurrentRadiationPressureForce( sourceName_ );


    // Update dependent variables
    receivedIrradiance = totalReceivedIrradiance;
    visibleAndEmittingSourcePanelCount = visibleAndEmittingSourcePanelCounter;

    // Calculate acceleration due to radiation pressure in global frame
    currentUnscaledAcceleration_ = targetRotationFromLocalToGlobalFrame * totalForceInTargetFrame / targetMassFunction_();
    scaleRadiationPressureAcceleration( );
}

} // tudat
} // electromagnetism
