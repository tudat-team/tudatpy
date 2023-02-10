/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/electromagnetism/sourcePanelRadiosityModel.h"

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/physicalConstants.h"


namespace tudat
{
namespace electromagnetism
{

using mathematical_constants::PI;

void SourcePanelRadiosityModel::updateMembers(
        const double panelLatitude,
        const double panelLongitude,
        const double currentTime)
{
    if (
        // New location
            (panelLatitude_ != panelLatitude || panelLongitude_ != panelLongitude) ||
            // New time and not time-invariant
            (!isTimeInvariant() && currentTime_ != currentTime))
    {
        panelLatitude_ = panelLatitude;
        panelLongitude_ = panelLongitude;
        currentTime_ = currentTime;

        updateMembers_(panelLatitude, panelLongitude, currentTime);
    }
}

double ConstantSourcePanelRadiosityModel::evaluateIrradianceAtPosition(
        double panelArea,
        const Eigen::Vector3d& panelSurfaceNormal,
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    const double cosBetweenNormalAndTarget = panelSurfaceNormal.dot(targetPosition.normalized());
    if (cosBetweenNormalAndTarget <= 0)
    {
        // Target is on backside of panel
        return 0;
    }

    // Determine irradiance from exitance based on source-to-target distance and emitting panel area
    const double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    const auto effectiveEmittingArea = cosBetweenNormalAndTarget * panelArea;
    const auto irradiance = constantRadiosity_ * effectiveEmittingArea / (PI * distanceSourceToTargetSquared);
    return irradiance;
}

double AlbedoSourcePanelRadiosityModel::evaluateIrradianceAtPosition(
        double panelArea,
        const Eigen::Vector3d& panelSurfaceNormal,
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    const Eigen::Vector3d targetDirection = targetPosition.normalized();
    const double cosBetweenNormalAndOriginalSource = panelSurfaceNormal.dot(-originalSourceToSourceDirection);
    const double cosBetweenNormalAndTarget = panelSurfaceNormal.dot(targetDirection);
    if (cosBetweenNormalAndOriginalSource <= 0 || cosBetweenNormalAndTarget <= 0)
    {
        // Target or original source are on backside of panel
        // This is checked by evaluateReflectedFraction as well, but check here to avoid unnecessary calculations/calls
        return 0;
    }

    const auto receivedIrradiance = cosBetweenNormalAndOriginalSource * originalSourceIrradiance;

    // Reflected fraction is given in [1/sr] (i.e. per solid angle)
    const auto reflectedFraction =
            reflectionLaw_->evaluateReflectedFraction(panelSurfaceNormal, originalSourceToSourceDirection, targetDirection);

    // Determine irradiance from reflected radiosity based on source-to-target distance and emitting panel area
    const double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    const auto effectiveEmittingArea = cosBetweenNormalAndTarget * panelArea;
    const auto albedoIrradiance =
            receivedIrradiance * reflectedFraction * effectiveEmittingArea / distanceSourceToTargetSquared;
    return albedoIrradiance;
}

void AlbedoSourcePanelRadiosityModel::updateMembers_(
        double panelLatitude,
        double panelLongitude,
        double currentTime)
{
    albedoDistribution_->updateMembers(currentTime);
    const auto albedo = albedoDistribution_->getValue(panelLatitude, panelLongitude);
    reflectionLaw_->setDiffuseReflectivity(albedo);
}

double DelayedThermalSourcePanelRadiosityModel::evaluateIrradianceAtPosition(
        double panelArea,
        const Eigen::Vector3d& panelSurfaceNormal,
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    const double cosBetweenNormalAndTarget = panelSurfaceNormal.dot(targetPosition.normalized());
    if (cosBetweenNormalAndTarget <= 0)
    {
        // Target is on backside of panel
        return 0;
    }

    // Exitance is only one quarter of original irradiance due to ratio of absorbing to emitting area (1:4)
    const auto emittedExitance = emissivity * originalSourceIrradiance / 4;

    // Determine irradiance from exitance based on source-to-target distance and emitting panel area
    const double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    const auto effectiveEmittingArea = cosBetweenNormalAndTarget * panelArea;
    const auto thermalIrradiance = emittedExitance * effectiveEmittingArea / (PI * distanceSourceToTargetSquared);
    return thermalIrradiance;
}

void DelayedThermalSourcePanelRadiosityModel::updateMembers_(
        double panelLatitude,
        double panelLongitude,
        double currentTime)
{
    emissivityDistribution_->updateMembers(currentTime);
    emissivity = emissivityDistribution_->getValue(panelLatitude, panelLongitude);
}

double AngleBasedThermalSourcePanelRadiosityModel::evaluateIrradianceAtPosition(
        double panelArea,
        const Eigen::Vector3d& panelSurfaceNormal,
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    const double cosBetweenNormalAndTarget = panelSurfaceNormal.dot(targetPosition.normalized());
    if (cosBetweenNormalAndTarget <= 0)
    {
        // Target is on backside of panel
        return 0;
    }

    const double cosBetweenNormalAndOriginalSource = panelSurfaceNormal.dot(-originalSourceToSourceDirection);
    const double positiveCosBetweenNormalAndOriginalSource = std::max(cosBetweenNormalAndOriginalSource, 0.);

    // Interpolate temperature using Lemoine (2013) Eq. 3
    const auto temperature = std::max(
            maxTemperature_ * pow(positiveCosBetweenNormalAndOriginalSource, 1./4),
            minTemperature_);
    // Calculate emissivity-corrected black-body radiation using Stefan-Boltzmann law
    const auto emittedExitance = emissivity * physical_constants::STEFAN_BOLTZMANN_CONSTANT * pow(temperature, 4);

    // Determine irradiance from exitance based on source-to-target distance and emitting panel area
    const double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    const auto effectiveEmittingArea = cosBetweenNormalAndTarget * panelArea;
    const auto thermalIrradiance = emittedExitance * effectiveEmittingArea / (PI * distanceSourceToTargetSquared);
    return thermalIrradiance;
}

void AngleBasedThermalSourcePanelRadiosityModel::updateMembers_(
        double panelLatitude,
        double panelLongitude,
        double currentTime)
{
    emissivityDistribution_->updateMembers(currentTime);
    emissivity = emissivityDistribution_->getValue(panelLatitude, panelLongitude);
}

} // tudat
} // electromagnetism
