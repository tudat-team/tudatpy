#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/basics/utilityMacros.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/linearAlgebra.h"
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

//*********************************************************************************************
//   Isotropic point radiation source
//*********************************************************************************************

IrradianceWithSourceList IsotropicPointRadiationSourceModel::evaluateIrradianceAtPosition(
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    return IrradianceWithSourceList {
        std::make_pair(evaluateIrradianceAtPosition(targetPosition), Eigen::Vector3d::Zero()) };

}
double IsotropicPointRadiationSourceModel::evaluateIrradianceAtPosition(const Eigen::Vector3d& targetPosition) const
{
    double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    auto luminosity = luminosityModel_->getLuminosity();

    auto sphereArea = 4 * PI * distanceSourceToTargetSquared;
    auto irradiance = luminosity / sphereArea;
    return irradiance;
}

void IsotropicPointRadiationSourceModel::updateMembers_(double currentTime)
{
    luminosityModel_->updateMembers();
}

//*********************************************************************************************
//   Paneled radiation source
//*********************************************************************************************

IrradianceWithSourceList PaneledRadiationSourceModel::evaluateIrradianceAtPosition(
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    IrradianceWithSourceList irradiances{};

    for (auto& panel : getPanels())
    {
        const Eigen::Vector3d targetPositionRelativeToPanel = targetPosition - panel.getRelativeCenter();

        double irradiance = 0;
        for (auto& radiosityModel : panel.getRadiosityModels())
        {
            irradiance += radiosityModel->evaluateIrradianceAtPosition(
                    panel,
                    targetPositionRelativeToPanel,
                    originalSourceIrradiance,
                    originalSourceToSourceDirection);
        }

        irradiances.emplace_back(irradiance, panel.getRelativeCenter());
    }

    return irradiances;
}

void StaticallyPaneledRadiationSourceModel::updateMembers_(double currentTime)
{
    // Check if panels have been initialized, generate them once it not
    if (panels_.size() != numberOfPanels)
    {
        panels_.clear();

        // Assume that all panels have the same area since they are evenly spaced on the sphere
        const auto bodyAverageRadius = sourceBodyShapeModel_->getAverageRadius();
        const auto totalBodyArea = 4 * PI * bodyAverageRadius * bodyAverageRadius;
        const auto panelArea = totalBodyArea / numberOfPanels;

        const auto pairOfAngleVectors = generateEvenlySpacedPoints(numberOfPanels);
        const auto polarAngles = std::get<0>(pairOfAngleVectors);
        const auto azimuthAngles = std::get<1>(pairOfAngleVectors);

        for (unsigned int i = 0; i < numberOfPanels; ++i)
        {
            const auto polarAngle = polarAngles[i];
            const auto azimuthAngle = azimuthAngles[i];
            //TODO-DOMINIK if oblate spheroid body, use actual position on surface
            const Eigen::Vector3d relativeCenter = coordinate_conversions::convertSphericalToCartesian(
                    Eigen::Vector3d(sourceBodyShapeModel_->getAverageRadius(), polarAngle, azimuthAngle));
            const Eigen::Vector3d surfaceNormal = relativeCenter.normalized();

            std::vector<std::shared_ptr<PanelRadiosityModel>> radiosityModels;
            for (auto& radiosityModelFunction : radiosityModelFunctions_)
            {
                radiosityModels.push_back(radiosityModelFunction(polarAngle, azimuthAngle));
            }

            panels_.emplace_back(
                    panelArea,
                    relativeCenter, surfaceNormal,
                    radiosityModels);
        }
    }
}

std::pair<std::vector<double>, std::vector<double>> generateEvenlySpacedPoints(unsigned int n)
{
    std::vector<double> polarAngles;
    std::vector<double> azimuthAngles;

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

double AlbedoPanelRadiosityModel::evaluateIrradianceAtPosition(
        const PaneledRadiationSourceModel::Panel& panel,
        const Eigen::Vector3d& targetPosition,
        double originalSourceIrradiance,
        const Eigen::Vector3d& originalSourceToSourceDirection) const
{
    const auto& surfaceNormal = panel.getSurfaceNormal();
    const auto cosBetweenNormalAndOriginalSource =
            linear_algebra::computeCosineOfAngleBetweenVectors(surfaceNormal, -originalSourceToSourceDirection);
    const auto cosBetweenNormalAndTarget =
            linear_algebra::computeCosineOfAngleBetweenVectors(surfaceNormal, targetPosition);

    const auto receivedIrradiance = cosBetweenNormalAndOriginalSource * originalSourceIrradiance;

    const Eigen::Vector3d targetDirection = targetPosition.normalized();
    const auto reflectedFraction =  // [1/sr]
            reflectionLaw_->evaluateReflectedFraction(surfaceNormal, originalSourceToSourceDirection, targetDirection);

    const double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    const auto effectiveEmittingArea = cosBetweenNormalAndTarget * panel.getArea();
    // No need to check if panel is illuminated or target would receive radiation because reaction vector is zero if not
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
    const auto cosBetweenNormalAndTarget =
            linear_algebra::computeCosineOfAngleBetweenVectors(surfaceNormal, targetPosition);
    if (cosBetweenNormalAndTarget <= 0)
    {
        // Target is on backside of panel
        return 0;
    }

    const auto emittedExitance = emissivity_ * originalSourceIrradiance / 4;

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
    const auto cosBetweenNormalAndOriginalSource =
            linear_algebra::computeCosineOfAngleBetweenVectors(surfaceNormal, -originalSourceToSourceDirection);
    const auto cosBetweenNormalAndTarget =
            linear_algebra::computeCosineOfAngleBetweenVectors(surfaceNormal, targetPosition);
    if (cosBetweenNormalAndTarget <= 0)
    {
        // Target is on backside of panel
        return 0;
    }

    const auto temperature = std::max(maxTemperature_ * pow(cosBetweenNormalAndOriginalSource, 1./4), minTemperature_);
    const auto emittedExitance = emissivity_ * physical_constants::STEFAN_BOLTZMANN_CONSTANT * pow(temperature, 4);

    const double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    const auto effectiveEmittingArea = cosBetweenNormalAndTarget * panel.getArea();
    const auto thermalIrradiance = emittedExitance * effectiveEmittingArea / (PI * distanceSourceToTargetSquared);
    return thermalIrradiance;
}

} // tudat
} // electromagnetism
