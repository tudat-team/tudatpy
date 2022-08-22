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
        const Eigen::Vector3d& occultedBodyPosition,
        const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedBodyShapeModel,
        const Eigen::Vector3d& targetPosition)
{
    const auto shadowFunction = mission_geometry::computeShadowFunction(
            occultedBodyPosition,
            occultedBodyShapeModel->getAverageRadius(),
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
    const Eigen::Vector3d occultingBodyToTargetVector = occultingBodyPosition - targetPosition;
    const Eigen::Vector3d sourceToTargetVector = sourcePosition - targetPosition;
    const double cosAngleBetweenOccultingBodyAndSource =
            linear_algebra::computeCosineOfAngleBetweenVectors(sourceToTargetVector, occultingBodyToTargetVector);
    const double angleBetweenOccultingBodyAndSource = acos(cosAngleBetweenOccultingBodyAndSource);
    const double apparentSeparation = occultingBodyToTargetVector.norm() * sin(angleBetweenOccultingBodyAndSource);

    const Eigen::Vector3d sourceToOccultingBodyVector = sourcePosition - occultingBodyPosition;
    const Eigen::Vector3d targetToOccultingBodyVector = targetPosition - occultingBodyPosition;
    const double cosAngleBetweenSourceAndTarget =
            linear_algebra::computeCosineOfAngleBetweenVectors(sourceToOccultingBodyVector, targetToOccultingBodyVector);

    const bool isSourceOutsideBodyRadiusSeenFromTarget = apparentSeparation > occultingBodyRadius;
    const bool isTargetInFrontOfOccultingBodyCenterSeenFromSource = cosAngleBetweenSourceAndTarget > 0;
    const bool isSourceVisibleFromTarget = isSourceOutsideBodyRadiusSeenFromTarget || isTargetInFrontOfOccultingBodyCenterSeenFromSource;
    return isSourceVisibleFromTarget;
}

} // tudat
} // electromagnetism
