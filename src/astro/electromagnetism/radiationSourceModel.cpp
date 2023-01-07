/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{
namespace electromagnetism
{

using mathematical_constants::PI;

void RadiationSourceModel::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

double RadiationSourceModel::evaluateTotalIrradianceAtPosition(
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    // Calculate irradiances due to all sub-sources
    auto irradiances = evaluateIrradianceAtPosition(
            targetPosition, originalSourceIrradiance, originalSourceToSourceDirection);

    // Sum contributions of all sub-sources
    double totalIrradiance = 0;
    for (auto& e: irradiances)
    {
        totalIrradiance += std::get<0>(e);
    }

    return totalIrradiance;
}

//*********************************************************************************************
//   Isotropic point radiation source
//*********************************************************************************************

IrradianceWithSourceList IsotropicPointRadiationSourceModel::evaluateIrradianceAtPosition(
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    // The radiation of an isotropic point source originates from the source center
    return IrradianceWithSourceList {
        std::make_pair(evaluateIrradianceAtPosition(targetPosition), Eigen::Vector3d::Zero()) };

}
double IsotropicPointRadiationSourceModel::evaluateIrradianceAtPosition(const Eigen::Vector3d& targetPosition) const
{
    double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    auto luminosity = luminosityModel_->getLuminosity();

    auto sphereArea = 4 * PI * distanceSourceToTargetSquared;

    // Since the source is isotropic, the radiation is uniformly distributed in all directions
    auto irradiance = luminosity / sphereArea;
    return irradiance;
}

void IsotropicPointRadiationSourceModel::updateMembers_(double currentTime)
{
    luminosityModel_->updateMembers(currentTime);
}

//*********************************************************************************************
//   Paneled radiation source
//*********************************************************************************************

void StaticallyPaneledRadiationSourceModel::updateMembers_(double currentTime)
{
    // Check if panels have been properly initialized, generate them once it not
    if (panels_.size() != numberOfPanels)
    {
        generatePanels();
    }
}

void StaticallyPaneledRadiationSourceModel::generatePanels()
{
    panels_.clear();

    // Assume that all panels have the same area since they are evenly spaced on the sphere
    const auto bodyAverageRadius = sourceBodyShapeModel_->getAverageRadius();
    const auto totalBodySurfaceArea = 4 * PI * bodyAverageRadius * bodyAverageRadius;
    const auto panelArea = totalBodySurfaceArea / numberOfPanels;

    // Generate center points of panels in spherical coordinates
    const auto pairOfAngleVectors = generateEvenlySpacedPoints(numberOfPanels);
    const auto polarAngles = std::get<0>(pairOfAngleVectors);
    const auto azimuthAngles = std::get<1>(pairOfAngleVectors);

    for (unsigned int i = 0; i < numberOfPanels; ++i)
    {
        const auto polarAngle = polarAngles[i];
        const auto azimuthAngle = azimuthAngles[i];
        //TODO-DOMINIK if oblate spheroid body, use actual position on surface
        const auto distanceFromSourceCenter = sourceBodyShapeModel_->getAverageRadius();

        // Calculate panel center relative to source center and surface normal in Cartesian coordinates
        // Surface normal is just vector from source center to panel center
        const Eigen::Vector3d relativeCenter = coordinate_conversions::convertSphericalToCartesian(
                Eigen::Vector3d(distanceFromSourceCenter, polarAngle, azimuthAngle));
        const Eigen::Vector3d surfaceNormal = relativeCenter.normalized();

        // Create radiosity models for panel
        std::vector<std::shared_ptr<PanelRadiosityModel>> radiosityModels;
        for (auto& radiosityModelFunction: radiosityModelFunctions_)
        {
            radiosityModels.push_back(radiosityModelFunction(polarAngle, azimuthAngle));
        }

        panels_.emplace_back(
                panelArea,
                relativeCenter, surfaceNormal,
                radiosityModels);
    }
}

IrradianceWithSourceList PaneledRadiationSourceModel::evaluateIrradianceAtPosition(
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    IrradianceWithSourceList irradiances{};

    for (auto& panel : getPanels())
    {
        const Eigen::Vector3d targetPositionRelativeToPanel = targetPosition - panel.getRelativeCenter();
        if (targetPositionRelativeToPanel.dot(panel.getSurfaceNormal()) <= 0)
        {
            // Avoids unnecessary panel radiosity model evaluations
            // No need to normalize target position here
            // TODO-DOMINIK find better way to eliminate unnecessary evaluations.
            //      maybe use data structure for fast lookup on sphere (spherical quad-tree, HEALPix)
            continue;
        }

        // The irradiance from a panel is the sum of the irradiances from all of its radiosity models
        double irradiance = 0;
        for (auto& radiosityModel : panel.getRadiosityModels())
        {
            irradiance += radiosityModel->evaluateIrradianceAtPosition(
                    panel,
                    targetPositionRelativeToPanel,
                    originalSourceIrradiance,
                    originalSourceToSourceDirection);
        }

        if (irradiance != 0)
        {
            // Do not add panels to list if they do not contribute to irradiance at target location
            // This prevents unnecessary evaluations in the radiation pressure acceleration evaluation
            irradiances.emplace_back(irradiance, panel.getRelativeCenter());
        }
    }

    return irradiances;
}

double AlbedoPanelRadiosityModel::evaluateIrradianceAtPosition(
        const PaneledRadiationSourceModel::Panel& panel,
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    const auto& surfaceNormal = panel.getSurfaceNormal();
    const Eigen::Vector3d targetDirection = targetPosition.normalized();
    const double cosBetweenNormalAndOriginalSource = surfaceNormal.dot(-originalSourceToSourceDirection);
    const double cosBetweenNormalAndTarget = surfaceNormal.dot(targetDirection);
    if (cosBetweenNormalAndOriginalSource <= 0 || cosBetweenNormalAndTarget <= 0)
    {
        // Target or original source are on backside of panel
        // This is checked by evaluateReflectedFraction as well, but check here to avoid unnecessary calculations/calls
        return 0;
    }

    const auto receivedIrradiance = cosBetweenNormalAndOriginalSource * originalSourceIrradiance;

    // Reflected fraction is given in [1/sr] (i.e. per solid angle)
    const auto reflectedFraction =
            reflectionLaw_->evaluateReflectedFraction(surfaceNormal, originalSourceToSourceDirection, targetDirection);

    // Determine irradiance from reflected radiosity based on source-to-target distance and emitting panel area
    const double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    const auto effectiveEmittingArea = cosBetweenNormalAndTarget * panel.getArea();
    const auto albedoIrradiance =
            receivedIrradiance * reflectedFraction * effectiveEmittingArea / distanceSourceToTargetSquared;
    return albedoIrradiance;
}

double DelayedThermalPanelRadiosityModel::evaluateIrradianceAtPosition(
        const PaneledRadiationSourceModel::Panel& panel,
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    const auto& surfaceNormal = panel.getSurfaceNormal();
    const double cosBetweenNormalAndTarget = surfaceNormal.dot(targetPosition.normalized());
    if (cosBetweenNormalAndTarget <= 0)
    {
        // Target is on backside of panel
        return 0;
    }

    // Exitance is only one quarter of original irradiance due to ratio of absorbing to emitting area (1:4)
    const auto emittedExitance = emissivity_ * originalSourceIrradiance / 4;

    // Determine irradiance from exitance based on source-to-target distance and emitting panel area
    const double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    const auto effectiveEmittingArea = cosBetweenNormalAndTarget * panel.getArea();
    const auto thermalIrradiance = emittedExitance * effectiveEmittingArea / (PI * distanceSourceToTargetSquared);
    return thermalIrradiance;
}

double AngleBasedThermalPanelRadiosityModel::evaluateIrradianceAtPosition(
        const PaneledRadiationSourceModel::Panel& panel,
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    const auto& surfaceNormal = panel.getSurfaceNormal();
    const double cosBetweenNormalAndTarget = surfaceNormal.dot(targetPosition.normalized());
    if (cosBetweenNormalAndTarget <= 0)
    {
        // Target is on backside of panel
        return 0;
    }

    const double cosBetweenNormalAndOriginalSource = surfaceNormal.dot(-originalSourceToSourceDirection);
    const double positiveCosBetweenNormalAndOriginalSource = std::max(cosBetweenNormalAndOriginalSource, 0.);

    // Interpolate temperature using Lemoine (2013) Eq. 3
    const auto temperature = std::max(
            maxTemperature_ * pow(positiveCosBetweenNormalAndOriginalSource, 1./4),
            minTemperature_);
    // Calculate emissivity-corrected black-body radiation using Stefan-Boltzmann law
    const auto emittedExitance = emissivity_ * physical_constants::STEFAN_BOLTZMANN_CONSTANT * pow(temperature, 4);

    // Determine irradiance from exitance based on source-to-target distance and emitting panel area
    const double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    const auto effectiveEmittingArea = cosBetweenNormalAndTarget * panel.getArea();
    const auto thermalIrradiance = emittedExitance * effectiveEmittingArea / (PI * distanceSourceToTargetSquared);
    return thermalIrradiance;
}

std::pair<std::vector<double>, std::vector<double>> generateEvenlySpacedPoints(unsigned int n)
{
    std::vector<double> polarAngles;
    std::vector<double> azimuthAngles;

    // Saff (1997) Eq. 8
    double previousAzimuthAngle{};
    for (unsigned int k = 1; k <= n; ++k)
    {
        double h = -1 + 2. * (k-1) / (n-1);
        double polarAngle = acos(h);
        double azimuthAngle;
        if (k == 1 || k == n)
        {
            azimuthAngle = 0.;
        }
        else
        {
            azimuthAngle = fmod(previousAzimuthAngle + 3.6 / sqrt(n * (1 - h * h)), 2 * PI);
        }

        polarAngles.push_back(polarAngle);
        azimuthAngles.push_back(azimuthAngle);

        previousAzimuthAngle = azimuthAngle;
    }

    return std::make_pair(polarAngles, azimuthAngles);
}

} // tudat
} // electromagnetism
