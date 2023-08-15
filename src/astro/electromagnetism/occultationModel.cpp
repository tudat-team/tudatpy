/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/electromagnetism/occultationModel.h"

#include <iostream>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/missionGeometry.h"
#include "tudat/math/basic/linearAlgebra.h"


namespace tudat
{
namespace electromagnetism
{

void OccultationModel::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

double SingleOccultingBodyOccultationModel::evaluateReceivedFractionFromExtendedSource(
        const Eigen::Vector3d& occultedSourcePosition,
        const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedSourceShapeModel,
        const Eigen::Vector3d& targetPosition) const
{
    const auto shadowFunction = mission_geometry::computeShadowFunction(
            occultedSourcePosition,
            occultedSourceShapeModel->getAverageRadius(),
            occultingBodyPosition,
            occultingBodyShapeModel_->getAverageRadius(),
            targetPosition);
    return shadowFunction;
}

double SingleOccultingBodyOccultationModel::evaluateReceivedFractionFromPointSource(
        const Eigen::Vector3d& occultedSourcePosition,
        const Eigen::Vector3d& targetPosition) const
{
    auto isVisible = evaluatePointToPointVisibilityWithOccultation(
            occultedSourcePosition,
            occultingBodyPosition,
            occultingBodyShapeModel_->getAverageRadius(),
            targetPosition);
    return static_cast<double>(isVisible);
}

void SingleOccultingBodyOccultationModel::updateMembers_(double currentTime)
{
    occultingBodyPosition = occultingBodyPositionFunction_();
}

double SimpleMultipleOccultingBodyOccultationModel::evaluateReceivedFractionFromExtendedSource(
        const Eigen::Vector3d& occultedSourcePosition,
        const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedSourceShapeModel,
        const Eigen::Vector3d& targetPosition) const
{
    double totalShadowFunction = 1.0;
    double shadowFunctionOfBody;
    unsigned int numberOfCurrentlyOccultingBodies = 0;
    for (unsigned int i = 0; i < getNumberOfOccultingBodies(); i++)
    {
        shadowFunctionOfBody = mission_geometry::computeShadowFunction(
                occultedSourcePosition,
                occultedSourceShapeModel->getAverageRadius(),
                occultingBodyPositions[i],
                occultingBodyShapeModels_[i]->getAverageRadius(),
                targetPosition);
        totalShadowFunction *= shadowFunctionOfBody;

        if (shadowFunctionOfBody < 1.0)
        {
            numberOfCurrentlyOccultingBodies++;
        }
    }

    if (numberOfCurrentlyOccultingBodies > 1)
    {
        std::cerr << "Warning, multiple occultation occurred, radiation pressure may be slightly underestimated" << std::endl;
    }

    return totalShadowFunction;
}

double SimpleMultipleOccultingBodyOccultationModel::evaluateReceivedFractionFromPointSource(
        const Eigen::Vector3d& occultedSourcePosition,
        const Eigen::Vector3d& targetPosition) const
{
    bool isVisible = true;
    for (unsigned int i = 0; i < getNumberOfOccultingBodies(); i++)
    {
        // Visibility for multiple occulting bodies is just the logical conjunction (AND-ing) of the individual
        // visibilities
        isVisible = isVisible && evaluatePointToPointVisibilityWithOccultation(
                occultedSourcePosition,
                occultingBodyPositions[i],
                occultingBodyShapeModels_[i]->getAverageRadius(),
                targetPosition);
    }
    return static_cast<double>(isVisible);
}

void SimpleMultipleOccultingBodyOccultationModel::updateMembers_(double currentTime)
{
    for (unsigned int i = 0; i < getNumberOfOccultingBodies(); i++)
    {
        occultingBodyPositions[i] = occultingBodyPositionFunctions_[i]();
    }
}

bool evaluatePointToPointVisibilityWithOccultation(
        const Eigen::Vector3d& occultedSourcePosition,
        const Eigen::Vector3d& occultingBodyPosition,
        double occultingBodyRadius,
        const Eigen::Vector3d& targetPosition)
{
    // Vallado (2013), Sec. 5.3.3
    const Eigen::Vector3d sourceToOccultingBodyVector = occultedSourcePosition - occultingBodyPosition;
    const Eigen::Vector3d targetToOccultingBodyVector = targetPosition - occultingBodyPosition;
    const double theta =
            acos(linear_algebra::computeCosineOfAngleBetweenVectors(sourceToOccultingBodyVector, targetToOccultingBodyVector));
    const double theta1 = acos(occultingBodyRadius / sourceToOccultingBodyVector.norm());
    const double theta2 = acos(occultingBodyRadius / targetToOccultingBodyVector.norm());

    const bool isSourceVisibleFromTarget = theta1 + theta2 > theta;
    return isSourceVisibleFromTarget;
}

} // electromagnetism
} // tudat
