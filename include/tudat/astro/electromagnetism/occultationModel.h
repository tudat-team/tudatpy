/*    Copyright (c) 2010-2022, Delft University of Technology
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

#include <vector>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"


namespace tudat
{
namespace electromagnetism
{

/*!
 * Class modeling the occultation of an occulted body due to occulting bodies as seen from a target position. This class
 * is only aware of the occulting bodies, not the occulted body or target.
 *
 * evaluateReceivedFractionFromExtendedSource() is used for occultation of an extended source (e.g., a sphere like the
 * sun), while evaluateReceivedFractionFromPointSource() is used for occultation of a point source (e.g., the center of
 * a panel on a paneled source). Therefore, both are used as occultation factor for received irradiance calculations in
 * the radiation pressure acceleration.
 */
class OccultationModel
{
public:
    /*!
    * Constructor.
    *
    * @param occultingBodyNames Names of the bodies that are occulting
    */
    explicit OccultationModel(const std::vector<std::string>& occultingBodyNames) :
        numberOfOccultingBodies(occultingBodyNames.size()),
        occultingBodyNames_(occultingBodyNames) {}

    virtual ~OccultationModel() = default;

    void updateMembers(double currentTime);

    /*!
     *  Evaluate how much of an occulted extended source is visible (i.e. the shadow function).
     *
     * @param occultedSourcePosition Position of the occulted source in global coordinates
     * @param occultedSourceShapeModel Shape model of the occulted source
     * @param targetPosition Position of the target from which occultation is observed in global coordinates
     * @return Visible fraction of the occulted source (between 0 and 1)
     */
    virtual double evaluateReceivedFractionFromExtendedSource(
            const Eigen::Vector3d& occultedSourcePosition,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedSourceShapeModel,
            const Eigen::Vector3d& targetPosition) const = 0;

    /*!
     * Evaluate how much of an occulted point source is visible (i.e. point-to-point visibility). This function is
     * commutative w.r.t. sourcePosition and targetPosition.
     *
     * @param occultedSourcePosition Position of the occulted source in global coordinates
     * @param targetPosition Position of the target from which occultation is observed in global coordinates
     * @return Visible fraction of the occulted source (either 0 and 1)
     */
    virtual double evaluateReceivedFractionFromPointSource(
            const Eigen::Vector3d& occultedSourcePosition,
            const Eigen::Vector3d& targetPosition) const = 0;

    std::vector<std::string> getOccultingBodyNames() const
    {
        return occultingBodyNames_;
    }

    unsigned int getNumberOfOccultingBodies() const
    {
        return numberOfOccultingBodies;
    }

private:
    virtual void updateMembers_(double currentTime) {}

    double currentTime_{TUDAT_NAN};
    unsigned int numberOfOccultingBodies;
    std::vector<std::string> occultingBodyNames_;
};

/*!
 * Occultation model without occultation, i.e. the received fraction is 1.0 and the visibility is true.
 */
class NoOccultingBodyOccultationModel : public OccultationModel
{
public:
    explicit NoOccultingBodyOccultationModel() :
            OccultationModel({}) {}

    double evaluateReceivedFractionFromExtendedSource(
            const Eigen::Vector3d& occultedSourcePosition,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedSourceShapeModel,
            const Eigen::Vector3d& targetPosition) const override
    {
        return 1.0;
    }

    double evaluateReceivedFractionFromPointSource(
            const Eigen::Vector3d& occultedSourcePosition,
            const Eigen::Vector3d& targetPosition) const override
    {
        return 1.0;
    }
};

/*!
 * Occultation model with a single occulting body. This evaluates the standard shadow function. The occulting body
 * shape is approximated as sphere of its actual shape model's average radius.
 */
class SingleOccultingBodyOccultationModel : public OccultationModel
{
public:
    SingleOccultingBodyOccultationModel(
            const std::string& occultingBodyName,
            const std::function<Eigen::Vector3d()>& occultingBodyPositionFunction,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultingBodyShapeModel) :
            OccultationModel({occultingBodyName}),
            occultingBodyPositionFunction_(occultingBodyPositionFunction),
            occultingBodyShapeModel_(occultingBodyShapeModel) {}

    double evaluateReceivedFractionFromExtendedSource(
            const Eigen::Vector3d& occultedSourcePosition,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedSourceShapeModel,
            const Eigen::Vector3d& targetPosition) const override;

    double evaluateReceivedFractionFromPointSource(
            const Eigen::Vector3d& occultedSourcePosition,
            const Eigen::Vector3d& targetPosition) const override;

private:
    void updateMembers_(double currentTime) override;

    std::function<Eigen::Vector3d()> occultingBodyPositionFunction_;
    std::shared_ptr<basic_astrodynamics::BodyShapeModel> occultingBodyShapeModel_;
    Eigen::Vector3d occultingBodyPosition;
};

/*!
 * Occultation model with multiple occulting bodies. This evaluates the product of the standard shadow function for each
 * occulting body. For non-concurrent occultations, this class behaves as expected. Multiple concurrent occultations of
 * extended sources may result in overestimated occultation (i.e. underestimated received irradiance and thus radiation
 * pressure). Occultation of multiple point sources is always handled correctly. The occulting body shapes are
 * approximated as spheres of their actual shape models' average radius.
 */
class SimpleMultipleOccultingBodyOccultationModel : public OccultationModel
{
public:
    SimpleMultipleOccultingBodyOccultationModel(
            const std::vector<std::string>& occultingBodyNames,
            const std::vector<std::function<Eigen::Vector3d()>>& occultingBodyPositionFunctions,
            const std::vector<std::shared_ptr<basic_astrodynamics::BodyShapeModel>>& occultingBodyShapeModels) :
            OccultationModel(occultingBodyNames),
            occultingBodyPositionFunctions_(occultingBodyPositionFunctions),
            occultingBodyShapeModels_(occultingBodyShapeModels),
            occultingBodyPositions(occultingBodyNames.size()) {}

    double evaluateReceivedFractionFromExtendedSource(
            const Eigen::Vector3d& occultedSourcePosition,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& occultedSourceShapeModel,
            const Eigen::Vector3d& targetPosition) const override;

    double evaluateReceivedFractionFromPointSource(
            const Eigen::Vector3d& occultedSourcePosition,
            const Eigen::Vector3d& targetPosition) const override;

private:
    void updateMembers_(double currentTime) override;

    std::vector<std::function<Eigen::Vector3d()>> occultingBodyPositionFunctions_;
    std::vector<std::shared_ptr<basic_astrodynamics::BodyShapeModel>> occultingBodyShapeModels_;
    std::vector<Eigen::Vector3d> occultingBodyPositions;
};

// TODO Realistic two-body occultation (DOI: 10.1016/j.asr.2018.02.002)

/*!
 * Evaluate whether two points have a line of sight with an occulting spherical body in between. There is a line of
 * sight if the apparent separation (distance between source and occulting body centers as seen from target) is larger
 * than the occulting body radius, or if the target is in front of the occulting body center (i.e. in the source side
 * of the plane through the occulting body center perpendicular to the source-occulting body vector).
 *
 * Since point-to-point visibility is commutative, sourcePosition and targetPosition are exchangeable.
 *
 * @param occultedSourcePosition Position of the occulted source in global coordinates
 * @param occultingBodyPosition Position of the occulting body in global coordinates
 * @param occultingBodyRadius Radius of the occulting body
 * @param targetPosition Position of the target from which occultation is observed in global coordinates
 * @return Whether the target can see the source
 */
bool evaluatePointToPointVisibilityWithOccultation(
        const Eigen::Vector3d& occultedSourcePosition,
        const Eigen::Vector3d& occultingBodyPosition,
        double occultingBodyRadius,
        const Eigen::Vector3d& targetPosition);

} // electromagnetism
} // tudat

#endif //TUDAT_OCCULTATIONMODEL_H
