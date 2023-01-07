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

double SingleOccultingBodyOccultationModel::evaluateReceivedFraction(
        const Eigen::Vector3d& occultedSourcePosition,
        const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedSourceShapeModel,
        const Eigen::Vector3d& targetPosition)
{
    const auto shadowFunction = mission_geometry::computeShadowFunction(
            occultedSourcePosition,
            occultedSourceShapeModel->getAverageRadius(),
            occultingBodyPosition,
            occultingBodyShapeModel_->getAverageRadius(),
            targetPosition);
    return shadowFunction;
}

bool SingleOccultingBodyOccultationModel::evaluateVisibility(
        const Eigen::Vector3d& sourcePosition,
        const Eigen::Vector3d& targetPosition)
{
    return evaluateVisibilityWithOccultation(
            sourcePosition,
            occultingBodyPosition,
            occultingBodyShapeModel_->getAverageRadius(),
            targetPosition);
}

void SingleOccultingBodyOccultationModel::updateMembers_(double currentTime)
{
    occultingBodyPosition = occultingBodyPositionFunction_();
}

bool evaluateVisibilityWithOccultation(
        const Eigen::Vector3d& sourcePosition,
        const Eigen::Vector3d& occultingBodyPosition,
        double occultingBodyRadius,
        const Eigen::Vector3d& targetPosition)
{
    // Geometry between source and occulting body as seen from target
    const Eigen::Vector3d occultingBodyToTargetVector = occultingBodyPosition - targetPosition;
    const Eigen::Vector3d sourceToTargetVector = sourcePosition - targetPosition;
    const double cosAngleBetweenOccultingBodyAndSource =
            linear_algebra::computeCosineOfAngleBetweenVectors(sourceToTargetVector, occultingBodyToTargetVector);
    const double angleBetweenOccultingBodyAndSource = acos(cosAngleBetweenOccultingBodyAndSource);
    const double apparentSeparation = occultingBodyToTargetVector.norm() * sin(angleBetweenOccultingBodyAndSource);

    // Geometry between source and target as seen from occulting body
    const Eigen::Vector3d sourceToOccultingBodyVector = sourcePosition - occultingBodyPosition;
    const Eigen::Vector3d targetToOccultingBodyVector = targetPosition - occultingBodyPosition;
    const double cosAngleBetweenSourceAndTarget =
            linear_algebra::computeCosineOfAngleBetweenVectors(sourceToOccultingBodyVector, targetToOccultingBodyVector);

    // To be visible, the source must be in front or to the side of the occulting body as seen from target
    const bool isSourceOutsideBodyRadiusSeenFromTarget = apparentSeparation > occultingBodyRadius;  // to the side?
    const bool isTargetInFrontOfOccultingBodyCenterSeenFromSource = cosAngleBetweenSourceAndTarget > 0;  // in front?
    const bool isSourceVisibleFromTarget = isSourceOutsideBodyRadiusSeenFromTarget || isTargetInFrontOfOccultingBodyCenterSeenFromSource;
    return isSourceVisibleFromTarget;
}

} // electromagnetism
} // tudat
