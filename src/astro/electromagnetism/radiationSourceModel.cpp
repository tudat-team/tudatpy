#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/basics/utilityMacros.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{
namespace electromagnetism
{

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

    auto sphereArea = 4 * mathematical_constants::PI * distanceSourceToTargetSquared;
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

void PaneledRadiationSourceModel::updateMembers_(double currentTime)
{
}

void StaticallyPaneledRadiationSourceModel::updateMembers_(double currentTime)
{
    PaneledRadiationSourceModel::updateMembers_(currentTime);

    // Check if panels have been initialized
    if (panels_.size() != n_)
    {
        panels_.clear();

        // Assume that all panels have the same area since they are evenly spaced on the sphere
        const auto bodyAverageRadius = sourceBodyShapeModel_->getAverageRadius();
        const auto totalBodyArea = 4 * mathematical_constants::PI * bodyAverageRadius * bodyAverageRadius;
        const auto panelArea = totalBodyArea / n_;

        const auto pairOfAngleVectors = generateEvenlySpacedPoints(n_);
        const auto polarAngles = std::get<0>(pairOfAngleVectors);
        const auto azimuthAngles = std::get<1>(pairOfAngleVectors);

        for (unsigned int i = 0; i < n_; ++i)
        {
            const auto polarAngle = polarAngles[i];
            const auto azimuthAngle = azimuthAngles[i];
            const Eigen::Vector3d relativeCenter = coordinate_conversions::convertSphericalToCartesian(
                    Eigen::Vector3d(sourceBodyShapeModel_->getAverageRadius(), polarAngle, azimuthAngle));
            const Eigen::Vector3d surfaceNormal = relativeCenter.normalized();

            panels_.emplace_back(
                    panelArea,
                    relativeCenter, surfaceNormal,
                    radiosityModelFunction_(polarAngle, azimuthAngle));
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
            azimuthAngle = fmod(previousAzimuthAngle + 3.6 / sqrt(n * (1 - h * h)), 2 * mathematical_constants::PI);
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

} // tudat
} // electromagnetism
