/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Oliver Montenbruck, et al. Satellite Orbits. Springer Berlin Heidelberg, 2000.
 */

#ifndef TUDAT_OCCULTATIONMODEL_H
#define TUDAT_OCCULTATIONMODEL_H

#include <Eigen/Core>
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"

namespace tudat
{
namespace electromagnetism
{

/*!
 * Only aware of occulting bodies, not occulted or radius
 */
class OccultationModel
{
public:
    explicit OccultationModel() = default;

    virtual ~OccultationModel() = default;

    void updateMembers(double currentTime);

    virtual double evaluateReceivedFraction(
            const Eigen::Vector3d& occultedBodyPosition,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedBodyShapeModel,
            const Eigen::Vector3d& targetPosition) = 0;

    virtual bool evaluateVisibility(
            const Eigen::Vector3d& sourcePosition,
            const Eigen::Vector3d& targetPosition) = 0;

private:
    virtual void updateMembers_(double currentTime) {}

    double currentTime_{TUDAT_NAN};
};

class NoOccultingBodyOccultationModel : public OccultationModel
{
public:
    double evaluateReceivedFraction(
            const Eigen::Vector3d& occultedBodyPosition,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedBodyShapeModel,
            const Eigen::Vector3d& targetPosition) override
    {
        return 1.0;
    }

    bool evaluateVisibility(
            const Eigen::Vector3d& sourcePosition,
            const Eigen::Vector3d& targetPosition) override
    {
        return true;
    }
};

class SingleOccultingBodyOccultationModel : public OccultationModel
{
public:
    SingleOccultingBodyOccultationModel(
            const std::function<Eigen::Vector3d()>& occultingBodyPositionFunction,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultingBodyShapeModel) :
            occultingBodyPositionFunction_(occultingBodyPositionFunction),
            occultingBodyShapeModel_(occultingBodyShapeModel) {}

    double evaluateReceivedFraction(
            const Eigen::Vector3d& occultedBodyPosition,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedBodyShapeModel,
            const Eigen::Vector3d& targetPosition) override;

    bool evaluateVisibility(
            const Eigen::Vector3d& sourcePosition,
            const Eigen::Vector3d& targetPosition) override;

private:
    void updateMembers_(double currentTime) override;

    std::function<Eigen::Vector3d()> occultingBodyPositionFunction_;
    std::shared_ptr<basic_astrodynamics::BodyShapeModel> occultingBodyShapeModel_;
    Eigen::Vector3d occultingBodyPosition;
};

/*!
 * Evaluate whether two points have a line of sight with an occulting spherical body in between. There is a line of
 * sight if the apparent separation (distance between source and occulting body centers as seen from target) is larger
 * than the occulting body radius, or if the target is in front of the occulting body center (i.e. in the source side
 * of the plane through the occulting body center perpendicular to the source-occulting body vector).
 *
 * Since point-to-point visibility is commutative, sourcePosition and targetPosition are exchangeable.
 *
 * @param sourcePosition
 * @param occultingBodyPosition
 * @param occultingBodyRadius
 * @param targetPosition
 * @return
 */
bool evaluateVisibilityWithOccultation(
        const Eigen::Vector3d& sourcePosition,
        const Eigen::Vector3d& occultingBodyPosition,
        double occultingBodyRadius,
        const Eigen::Vector3d& targetPosition);

} // tudat
} // electromagnetism

#endif //TUDAT_OCCULTATIONMODEL_H
