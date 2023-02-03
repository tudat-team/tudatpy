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

#include <vector>
#include <utility>
#include <memory>

#include <Eigen/Core>

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
    // No need to check if model has been updated this timestep, already done by RadiationSourceModel
    for (auto& panel : panels_)
    {
        panel.updateMembers(currentTime);
    }
}

void StaticallyPaneledRadiationSourceModel::generatePanels(
        const std::vector<std::unique_ptr<SourcePanelRadiosityModel>>& baseRadiosityModels)
{
    panels_.clear();

    // Assume that all panels have the same area since they are evenly spaced on the sphere
    const auto bodyAverageRadius = sourceBodyShapeModel_->getAverageRadius();
    const auto totalBodySurfaceArea = 4 * PI * bodyAverageRadius * bodyAverageRadius;
    const auto panelArea = totalBodySurfaceArea / numberOfPanels;

    // Generate center points of panels in spherical coordinates
    const auto pairOfAngleVectors = generateEvenlySpacedPoints_Staggered(numberOfPanels);
    const auto polarAngles = std::get<0>(pairOfAngleVectors);
    const auto azimuthAngles = std::get<1>(pairOfAngleVectors);

    for (unsigned int i = 0; i < numberOfPanels; ++i)
    {
        const auto polarAngle = polarAngles[i];
        const auto azimuthAngle = azimuthAngles[i];
        //TODO-DOMINIK if oblate spheroid body, use actual position on surface
        const auto distanceFromSourceCenter = sourceBodyShapeModel_->getAverageRadius();

        // Calculate panel center relative to source center and surface normal in Cartesian coordinates
        // Surface normal is just vector from source center to panel center for sphere
        const Eigen::Vector3d relativeCenter = coordinate_conversions::convertSphericalToCartesian(
                Eigen::Vector3d(distanceFromSourceCenter, polarAngle, azimuthAngle));
        const Eigen::Vector3d surfaceNormal = relativeCenter.normalized();

        // Create radiosity models for panel
        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels;
        for (auto& baseRadiosityModel: baseRadiosityModels)
        {
            // Radiosity models are cloned such that they share their albedo/emissivity distribution
            radiosityModels.push_back(baseRadiosityModel->clone());
        }

        panels_.emplace_back(
                panelArea,
                relativeCenter, surfaceNormal,
                std::move(radiosityModels));
    }
}

IrradianceWithSourceList PaneledRadiationSourceModel::evaluateIrradianceAtPosition(
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    IrradianceWithSourceList irradiances{};

    for (const auto& panel : getPanels())
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
                    panel.getArea(),
                    panel.getSurfaceNormal(),
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

void PaneledRadiationSourceModel::Panel::updateMembers(double currentTime)
{
    // No need to check if model has been updated this timestep, already done by RadiationSourceModel
    for (const auto& radiosityModel : radiosityModels_)
    {
        radiosityModel->updateMembers(latitude_, longitude_, currentTime);
    }
}

std::pair<std::vector<double>, std::vector<double>> generateEvenlySpacedPoints_Spiraling(unsigned int n)
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

std::pair<std::vector<double>, std::vector<double>> generateEvenlySpacedPoints_Staggered(unsigned int n)
{
    std::vector<double> polarAngles;
    std::vector<double> azimuthAngles;

    double azimuthStep = PI * (3 - sqrt(5));

    // Wetterer (2014) Eqs. 34 + 35
    double previousAzimuthAngle{};
    double z = -1 + 1./n;
    for (unsigned int j = 1; j <= n; ++j)
    {
        // Numerator and denominator of the atan argument seem to be switched in the given equation, which produced
        // points that are skewed towards the poles
        double polarAngle = M_PI_2 - atan(z / sqrt(1 - z*z));
        double azimuthAngle;
        if (j == 1)
        {
            azimuthAngle = 0.;
        }
        else
        {
            azimuthAngle = fmod(previousAzimuthAngle + azimuthStep, 2 * PI);
        }

        polarAngles.push_back(polarAngle);
        azimuthAngles.push_back(azimuthAngle);

        previousAzimuthAngle = azimuthAngle;
        // We use increasing z instead of decreasing z as in the given equation to produce points from southmost to
        // northmost, which is consistent with Saff's algorithm
        z += 2./n;
    }

    return std::make_pair(polarAngles, azimuthAngles);
}

} // tudat
} // electromagnetism
