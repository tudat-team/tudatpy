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
#include <Eigen/Geometry>

#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"


namespace tudat
{
namespace electromagnetism
{

using mathematical_constants::PI;
using basic_mathematics::computeModulo;

void RadiationSourceModel::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

double RadiationSourceModel::evaluateTotalIrradianceAtPosition(
        const Eigen::Vector3d& targetPosition)
{
    // Calculate irradiances due to all sub-sources
    auto irradiances = evaluateIrradianceAtPosition(targetPosition);

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
        const Eigen::Vector3d& targetPosition)
{
    double distanceSourceToTargetSquared = targetPosition.squaredNorm();
    auto luminosity = luminosityModel_->getLuminosity();

    auto sphereArea = 4 * PI * distanceSourceToTargetSquared;

    // Since the source is isotropic, the radiation is uniformly distributed in all directions
    auto irradiance = luminosity / sphereArea;
    // The radiation of an isotropic point source originates from the source center
    return IrradianceWithSourceList { std::make_pair(irradiance, Eigen::Vector3d::Zero()) };
}

void IsotropicPointRadiationSourceModel::updateMembers_(double currentTime)
{
    luminosityModel_->updateMembers(currentTime);
}

//*********************************************************************************************
//   Paneled radiation source
//*********************************************************************************************

IrradianceWithSourceList PaneledRadiationSourceModel::evaluateIrradianceAtPosition(
        const Eigen::Vector3d& targetPosition)
{
    IrradianceWithSourceList irradiances{};

    visibleArea = 0;
    for (const auto& panel : getPanels())
    {
        const Eigen::Vector3d targetPositionRelativeToPanel = targetPosition - panel.getRelativeCenter();
        if (targetPositionRelativeToPanel.dot(panel.getSurfaceNormal()) <= 0)
        {
            // Avoids unnecessary panel radiosity model evaluations
            // No need to normalize target position here
            continue;
        }

        visibleArea += panel.getArea();

        // The irradiance from a panel is the sum of the irradiances from all of its radiosity models
        double irradiance = 0;
        for (auto& radiosityModel : panel.getRadiosityModels())
        {
            irradiance += radiosityModel->evaluateIrradianceAtPosition(
                    panel.getArea(),
                    panel.getSurfaceNormal(),
                    targetPositionRelativeToPanel);
        }

        if (irradiance > 0)
        {
            // Do not add panels to list if they do not contribute to irradiance at target location
            // This prevents unnecessary evaluations in the radiation pressure acceleration
            irradiances.emplace_back(irradiance, panel.getRelativeCenter());
        }
    }

    return irradiances;
}

void StaticallyPaneledRadiationSourceModel::updateMembers_(double currentTime)
{
    sourcePanelRadiosityModelUpdater_->updateMembers(currentTime);
    for (auto& panel : panels_)
    {
        panel.updateMembers(currentTime);
        sourcePanelRadiosityModelUpdater_->updatePanel(panel);
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

DynamicallyPaneledRadiationSourceModel::DynamicallyPaneledRadiationSourceModel(
        const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& sourceBodyShapeModel,
        std::unique_ptr<SourcePanelRadiosityModelUpdater> sourcePanelRadiosityModelUpdater,
        const std::vector<std::unique_ptr<SourcePanelRadiosityModel>>& baseRadiosityModels,
        const std::vector<int>& numberOfPanelsPerRing) :
        PaneledRadiationSourceModel(sourceBodyShapeModel, std::move(sourcePanelRadiosityModelUpdater)),
        numberOfPanelsPerRing_(numberOfPanelsPerRing)
{
    if( sourceBodyShapeModel == nullptr )
    {
        throw std::runtime_error( "Error when creating dynamically panelled radiation source model; no shape model defined" );
    }
    numberOfPanels = 1;
    for (const auto& numberOfPanelsInCurrentRing : numberOfPanelsPerRing)
    {
        numberOfPanels += numberOfPanelsInCurrentRing;
    }

    for (unsigned int i = 0; i < numberOfPanels; ++i)
    {
        // Create radiosity models for panel
        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels;
        for (auto& baseRadiosityModel: baseRadiosityModels)
        {
            // Radiosity models are cloned such that they share their albedo/emissivity distribution
            radiosityModels.push_back(baseRadiosityModel->clone());
        }

        // Initialize with empty panels
        panels_.emplace_back(
                TUDAT_NAN,
                Eigen::Vector3d(TUDAT_NAN, TUDAT_NAN, TUDAT_NAN),
                Eigen::Vector3d(TUDAT_NAN, TUDAT_NAN, TUDAT_NAN),
                std::move(radiosityModels));
    }
}

IrradianceWithSourceList DynamicallyPaneledRadiationSourceModel::evaluateIrradianceAtPosition(
        const Eigen::Vector3d& targetPosition)
{
    // Generate center points of panels in spherical coordinates
    const auto panelProperties = generatePaneledSphericalCap_EqualProjectedAttenuatedArea(
            targetPosition, numberOfPanelsPerRing_, sourceBodyShapeModel_->getAverageRadius());
    const auto panelCenters = std::get<0>(panelProperties);
    const auto polarAngles = std::get<1>(panelProperties);
    const auto azimuthAngles = std::get<2>(panelProperties);
    const auto areas = std::get<3>(panelProperties);

    for (unsigned int i = 0; i < numberOfPanels; ++i)
    {
        const auto& relativeCenter = panelCenters[i];
        const Eigen::Vector3d surfaceNormal = relativeCenter.normalized();
        const auto polarAngle = polarAngles[i];
        const auto azimuthAngle = azimuthAngles[i];
        const auto area = areas[i];

        panels_[i].setRelativeCenter(relativeCenter, polarAngle, azimuthAngle);
        panels_[i].setSurfaceNormal(surfaceNormal);
        panels_[i].setArea(area);

        // Always update (independently of current time) because evaluation may come from different targets each call
        panels_[i].updateMembers(currentTime_);
        sourcePanelRadiosityModelUpdater_->updatePanel(panels_[i]);
    }

    return PaneledRadiationSourceModel::evaluateIrradianceAtPosition(targetPosition);
}

void DynamicallyPaneledRadiationSourceModel::updateMembers_(double currentTime)
{
    sourcePanelRadiosityModelUpdater_->updateMembers(currentTime);
}

void RadiationSourcePanel::updateMembers(double currentTime)
{
    for (const auto& radiosityModel : radiosityModels_)
    {
        radiosityModel->updateMembers(latitude_, longitude_, currentTime);
    }
}

void SourcePanelRadiosityModelUpdater::updateMembers(const double currentTime)
{
    // Update all properties that only depend on the source center, not on the panels
    // This includes the direction and irradiance from the original source, since it is assumed far enough away from
    // the source so that the panel position does not make a difference

    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;

        Eigen::Vector3d sourceCenterPositionInGlobalFrame = sourcePositionFunction_(); // position of center of source (e.g. planet)
        Eigen::Quaterniond sourceRotationFromLocalToGlobalFrame = sourceRotationFromLocalToGlobalFrameFunction_();
        Eigen::Quaterniond sourceRotationFromGlobalToLocalFrame = sourceRotationFromLocalToGlobalFrame.inverse();

        for ( const auto &kv : originalSourcePositionFunctions_ ) {
            auto originalSourceName = kv.first;
            auto originalSourcePositionFunction = kv.second;
            Eigen::Vector3d originalSourceCenterPositionInGlobalFrame = originalSourcePositionFunction();

            // Evaluate irradiances from original source at source center in original source frame (for albedo-reflected radiation)
            // If other types than isotropic point sources are supported as original source, rotate to original source frame here
            Eigen::Vector3d sourceCenterPositionInOriginalSourceFrame = sourceCenterPositionInGlobalFrame - originalSourceCenterPositionInGlobalFrame;
            originalSourceUnoccultedIrradiances_[originalSourceName] =
                    originalSourceModels_[originalSourceName]->evaluateIrradianceAtPosition(sourceCenterPositionInOriginalSourceFrame).front().first;

            Eigen::Vector3d originalSourceToSourceDirectionInSourceFrame =
                    sourceRotationFromGlobalToLocalFrame
                    * -(originalSourceCenterPositionInGlobalFrame - sourceCenterPositionInGlobalFrame).normalized();
            originalSourceToSourceCenterDirections_[originalSourceName] = originalSourceToSourceDirectionInSourceFrame;

            originalSourceToSourceOccultationModels_[originalSourceName]->updateMembers(currentTime);
        }
    }
}

void SourcePanelRadiosityModelUpdater::updatePanel(
        RadiationSourcePanel& panel)
{
    // Update all properties that depend on the panels
    Eigen::Vector3d sourceCenterPositionInGlobalFrame = sourcePositionFunction_(); // position of center of source (e.g. planet)
    Eigen::Quaterniond sourceRotationFromLocalToGlobalFrame = sourceRotationFromLocalToGlobalFrameFunction_();

    for(auto& radiosityModel : panel.getRadiosityModels())
    {
        if (!radiosityModel->dependsOnOriginalSource())
        {
            continue;
        }

        Eigen::Vector3d sourcePositionInGlobalFrame =
                sourceCenterPositionInGlobalFrame + sourceRotationFromLocalToGlobalFrame * panel.getRelativeCenter();

        auto* originalSourceDependentRadiosityModel =
                static_cast<OriginalSourceDependentSourcePanelRadiosityModel*>(radiosityModel.get());
        auto originalSourceName = originalSourceDependentRadiosityModel->getOriginalSourceName();
        Eigen::Vector3d originalSourceCenterPositionInGlobalFrame = originalSourcePositionFunctions_[originalSourceName]();

        auto originalSourceToSourceReceivedFraction =
                originalSourceToSourceOccultationModels_[originalSourceName]->evaluateReceivedFractionFromExtendedSource(
                originalSourceCenterPositionInGlobalFrame,
                originalSourceBodyShapeModels_[originalSourceName],
                sourcePositionInGlobalFrame);
        auto originalSourceOccultedIrradiance =
                originalSourceUnoccultedIrradiances_[originalSourceName] * originalSourceToSourceReceivedFraction;
        originalSourceDependentRadiosityModel->updateOriginalSourceProperties(
                originalSourceUnoccultedIrradiances_[originalSourceName],
                originalSourceOccultedIrradiance,
                originalSourceToSourceCenterDirections_[originalSourceName]);
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
            azimuthAngle = computeModulo(previousAzimuthAngle + 3.6 / sqrt(n * (1 - h * h)), 2 * PI);
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
        double polarAngle = mathematical_constants::PI / 2.0 - atan(z / sqrt(1 - z*z));
        double azimuthAngle;
        if (j == 1)
        {
            azimuthAngle = 0.;
        }
        else
        {
            azimuthAngle = computeModulo(previousAzimuthAngle + azimuthStep, 2 * PI);
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

std::tuple<std::vector<Eigen::Vector3d>, std::vector<double>, std::vector<double>, std::vector<double>>
generatePaneledSphericalCap_EqualAngularResolution(
        const Eigen::Vector3d& targetPosition,
        const std::vector<int>& numberOfPanelsPerRing,
        double bodyRadius)
{
    // The panels are first generated as if the target were above the north pole ("pole-aligned frame"),
    // then rotated to the actual position ("target-aligned frame"). This works because a spherical body
    // is assumed so that panel areas do not change upon rotation.

    std::vector<Eigen::Vector3d> panelCenters;
    std::vector<double> polarAngles;
    std::vector<double> azimuthAngles;
    std::vector<double> areas;

    const Eigen::Quaterniond rotationFromPoleAlignedToTargetAlignedFrame =
            Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), targetPosition);

    const auto numberOfRings = numberOfPanelsPerRing.size();
    const auto sphericalCapAngle = acos(bodyRadius / targetPosition.norm());
    // Angular distance between rings
    const auto angularResolutionPolar = sphericalCapAngle / (numberOfRings + 1);

    // Create central cap
    const Eigen::Vector3d centralCapCenterInTargetAlignedFrameCartesian = targetPosition.normalized() * bodyRadius;
    const Eigen::Vector3d centralCapCenterInTargetAlignedFrameSpherical =
            coordinate_conversions::convertCartesianToSpherical(centralCapCenterInTargetAlignedFrameCartesian);
    const auto centralCapCenterPolarAngleInTargetAlignedFrame = centralCapCenterInTargetAlignedFrameSpherical[1];
    const auto centralCapCenterAzimuthAngleInTargetAlignedFrame = computeModulo(centralCapCenterInTargetAlignedFrameSpherical[2], 2 * PI);
    const auto centralCapArea = 2 * PI * bodyRadius * bodyRadius * (1 - cos(angularResolutionPolar));

    panelCenters.push_back(centralCapCenterInTargetAlignedFrameCartesian);
    polarAngles.push_back(centralCapCenterPolarAngleInTargetAlignedFrame);
    azimuthAngles.push_back(centralCapCenterAzimuthAngleInTargetAlignedFrame);
    areas.push_back(centralCapArea);

    // Create rings
    for (unsigned int currentRingNumber = 0; currentRingNumber < numberOfRings; currentRingNumber++)
    {
        int numberOfPanelsInCurrentRing = numberOfPanelsPerRing[currentRingNumber];
        // First ring stretches from 1*angularResolutionPolar to 2*angularResolutionPolar, so its
        // center is at 1.5*angularResolutionPolar
        double panelCenterPolarAngleInPoleAlignedFrame = (1.5 + currentRingNumber) * angularResolutionPolar;
        // Angular distance between panels within ring
        double angularResolutionAzimuth = 2 * PI / numberOfPanelsInCurrentRing;

        // Area of a sphere sector bounded by constant-polar/constant-azimuth angle lines
        // Panels within the same ring have the same area
        double panelArea =
                2 * bodyRadius * bodyRadius * angularResolutionAzimuth
                * sin(angularResolutionPolar / 2) * sin(panelCenterPolarAngleInPoleAlignedFrame);

        // Create panels of ring
        for (int currentPanelNumber = 0; currentPanelNumber < numberOfPanelsInCurrentRing; currentPanelNumber++)
        {
            double panelCenterAzimuthAngleInPoleAlignedFrame = currentPanelNumber * angularResolutionAzimuth;

            // Rotate from pole-aligned to target-aligned frame in Cartesian coordinates
            Eigen::Vector3d panelCenterInPoleAlignedFrameSpherical(
                    bodyRadius, panelCenterPolarAngleInPoleAlignedFrame, panelCenterAzimuthAngleInPoleAlignedFrame);
            Eigen::Vector3d panelCenterInTargetAlignedFrameCartesian =
                    rotationFromPoleAlignedToTargetAlignedFrame * coordinate_conversions::convertSphericalToCartesian(panelCenterInPoleAlignedFrameSpherical);
            Eigen::Vector3d panelCenterInTargetAlignedFrameSpherical =
                    coordinate_conversions::convertCartesianToSpherical(panelCenterInTargetAlignedFrameCartesian);

            double panelCenterPolarAngleInTargetAlignedFrame = panelCenterInTargetAlignedFrameSpherical[1];
            double panelCenterAzimuthAngleInTargetAlignedFrame = computeModulo(panelCenterInTargetAlignedFrameSpherical[2], 2 * PI);

            panelCenters.push_back(panelCenterInTargetAlignedFrameCartesian);
            polarAngles.push_back(panelCenterPolarAngleInTargetAlignedFrame);
            azimuthAngles.push_back(panelCenterAzimuthAngleInTargetAlignedFrame);
            areas.push_back(panelArea);
        }
    }

    return std::make_tuple(panelCenters, polarAngles, azimuthAngles, areas);
}

std::tuple<std::vector<Eigen::Vector3d>, std::vector<double>, std::vector<double>, std::vector<double>>
generatePaneledSphericalCap_EqualProjectedAttenuatedArea(
        const Eigen::Vector3d& targetPosition,
        const std::vector<int>& numberOfPanelsPerRing,
        double R_e)
{
    // The panels are first generated as if the target were above the north pole ("polar-centric frame"),
    // then rotated to the actual position ("target-centric frame"). This works because a spherical body
    // is assumed so that panel areas do not change upon rotation.

    // Algorithm adapted from Knocke (1989), Appendix A
    // Nomenclature from original algorithm:
    //  - N: total number of panels
    //  - N_s: number of panels in current ring
    //  - A_prime: projected, attenuated area of each panel
    //  - r_s: distance between target and body center
    //  - R_e: body radius
    //  - k: number of panels in current and more inward rings
    //  - zeta, beta, beta_star: defined in Fig. 2.4 and Fig. A.1
    //  - gamma: identical to the viewing angle alpha

    std::vector<Eigen::Vector3d> panelCenters;
    std::vector<double> polarAngles;
    std::vector<double> azimuthAngles;
    std::vector<double> areas;

    std::vector<double> betas;

    const Eigen::Quaterniond rotationFromPoleAlignedToTargetAlignedFrame =
            Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), targetPosition);

    const auto numberOfRings = numberOfPanelsPerRing.size();
    int N = 1;
    for (const auto& N_s : numberOfPanelsPerRing) {
        N += N_s;
    }

    const auto r_s = targetPosition.norm();

    // Calculate ring boundaries
    const auto zeta_m = asin(R_e / r_s);
    const auto zeta_1 = acos((N - 1 + cos(zeta_m)) / N);
    const auto gamma_1 = asin(std::min(1.0, r_s * sin(zeta_1) / R_e));
    betas.push_back(gamma_1 - zeta_1);

    int k = 1;
    for (const auto& N_s : numberOfPanelsPerRing) {
        k += N_s;
        auto zeta_i = acos(k * cos(zeta_1) - k + 1);
        // min is necessary because argument may slightly exceed 1.0 due to floating point errors
        auto gamma_i = asin(std::min(1.0, r_s * sin(zeta_i) / R_e));
        betas.push_back(gamma_i - zeta_i);
    }

    // Create central cap
    const Eigen::Vector3d centralCapCenterInTargetAlignedFrameCartesian = targetPosition.normalized() * R_e;
    const Eigen::Vector3d centralCapCenterInTargetAlignedFrameSpherical =
            coordinate_conversions::convertCartesianToSpherical(centralCapCenterInTargetAlignedFrameCartesian);
    const auto centralCapCenterPolarAngleInTargetAlignedFrame = centralCapCenterInTargetAlignedFrameSpherical[1];
    const auto centralCapCenterAzimuthAngleInTargetAlignedFrame = computeModulo(centralCapCenterInTargetAlignedFrameSpherical[2], 2 * PI);
    const auto centralCapArea = 2 * PI * R_e * R_e * (1 - cos(betas.front()));

    panelCenters.push_back(centralCapCenterInTargetAlignedFrameCartesian);
    polarAngles.push_back(centralCapCenterPolarAngleInTargetAlignedFrame);
    azimuthAngles.push_back(centralCapCenterAzimuthAngleInTargetAlignedFrame);
    areas.push_back(centralCapArea);

    // Create rings
    for (unsigned int currentRingNumber = 0; currentRingNumber < numberOfRings; currentRingNumber++)
    {
        int N_s = numberOfPanelsPerRing[currentRingNumber];

        // Angular distance between panels within ring
        double angularResolutionAzimuth = 2 * PI / N_s;

        // Ring center is polar-angle-wise halfway between both boundaries
        double beta_star = (betas[currentRingNumber] + betas[currentRingNumber + 1]) / 2;

        // The panel area could also be calculated from the constant A'. This has been implemented here:
        // https://github.com/DominikStiller/tudat/blob/d58c9840af0bac16026e313bb95461cbda290c3e/src/astro/electromagnetism/radiationSourceModel.cpp#L495
        // However, this leads to extremely large areas for the outer panels, since a constant viewing angle (alpha) is
        // assumed. Calculating the panel area from the sphere geometry, as done here, gives a realistic panel area.
        // Experiments for LAGEOS-1 showed that the resulting RP accelerations for both area calculation approaches
        // agree within 2%. Both converge for a large number of rings, since the outer panels are smaller then.
        double panelArea =  2 * PI * R_e * R_e * (cos(betas[currentRingNumber]) - cos(betas[currentRingNumber + 1])) / N_s;

        // Create panels of ring
        for (int currentPanelNumber = 0; currentPanelNumber < N_s; currentPanelNumber++)
        {
            double panelCenterAzimuthAngleInPoleAlignedFrame = currentPanelNumber * angularResolutionAzimuth;

            // Rotate from pole-aligned to target-aligned frame in Cartesian coordinates
            Eigen::Vector3d panelCenterInPoleAlignedFrameSpherical(
                    R_e, beta_star, panelCenterAzimuthAngleInPoleAlignedFrame);
            Eigen::Vector3d panelCenterInTargetAlignedFrameCartesian =
                    rotationFromPoleAlignedToTargetAlignedFrame * coordinate_conversions::convertSphericalToCartesian(panelCenterInPoleAlignedFrameSpherical);
            Eigen::Vector3d panelCenterInTargetAlignedFrameSpherical =
                    coordinate_conversions::convertCartesianToSpherical(panelCenterInTargetAlignedFrameCartesian);

            double panelCenterPolarAngleInTargetAlignedFrame = panelCenterInTargetAlignedFrameSpherical[1];
            double panelCenterAzimuthAngleInTargetAlignedFrame = computeModulo(panelCenterInTargetAlignedFrameSpherical[2], 2 * PI);

            panelCenters.push_back(panelCenterInTargetAlignedFrameCartesian);
            polarAngles.push_back(panelCenterPolarAngleInTargetAlignedFrame);
            azimuthAngles.push_back(panelCenterAzimuthAngleInTargetAlignedFrame);
            areas.push_back(panelArea);
        }
    }

    return std::make_tuple(panelCenters, polarAngles, azimuthAngles, areas);
}

} // tudat
} // electromagnetism
