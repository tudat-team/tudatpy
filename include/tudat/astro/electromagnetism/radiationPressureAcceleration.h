/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIATIONPRESSUREACCELERATION_H
#define TUDAT_RADIATIONPRESSUREACCELERATION_H

#include <functional>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/astro/electromagnetism/occultationModel.h"


namespace tudat
{
namespace electromagnetism
{

/*!
 * Class modeling radiation pressure acceleration. Radiation pressure accelerates a target due to electromagnetic
 * radiation from a source.
 */
class RadiationPressureAcceleration: public basic_astrodynamics::AccelerationModel3d
{
public:
    /*!
     * Update class members.
     *
     * @param currentTime Current simulation time
     */
    void updateMembers(double currentTime) override;

    virtual std::shared_ptr<RadiationSourceModel> getSourceModel() const = 0;

    std::shared_ptr<RadiationPressureTargetModel> getTargetModel() const
    {
        return targetModel_;
    }

    double getReceivedIrradiance() const
    {
        return receivedIrradiance;
    }

protected:
    RadiationPressureAcceleration(const std::function<Eigen::Vector3d()>& sourcePositionFunction,
                                  const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
                                  const std::function<Eigen::Vector3d()>& targetPositionFunction,
                                  const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
                                  const std::function<double()>& targetMassFunction,
                                  const std::shared_ptr<OccultationModel>& sourceToTargetOccultationModel) :
            sourcePositionFunction_(sourcePositionFunction),
            targetModel_(targetModel),
            targetPositionFunction_(targetPositionFunction),
            targetRotationFromLocalToGlobalFrameFunction_(targetRotationFromLocalToGlobalFrameFunction),
            targetMassFunction_(targetMassFunction),
            sourceToTargetOccultationModel_(sourceToTargetOccultationModel),
            receivedIrradiance(TUDAT_NAN) {}

    virtual Eigen::Vector3d calculateAcceleration() = 0;

    virtual void updateMembers_(const double currentTime) {};

    std::function<Eigen::Vector3d()> sourcePositionFunction_;
    std::shared_ptr<RadiationPressureTargetModel> targetModel_;
    std::function<Eigen::Vector3d()> targetPositionFunction_;
    std::function<Eigen::Quaterniond()> targetRotationFromLocalToGlobalFrameFunction_;
    std::function<double()> targetMassFunction_;
    std::shared_ptr<OccultationModel> sourceToTargetOccultationModel_;

    // For dependent variable
    double receivedIrradiance;
};

/*!
 * Class modeling radiation pressure acceleration from an isotropic point source.
 *
 * Assumptions:
 *  - Radiation from the source is evaluated at the target center (i.e. the target extent is neglected,
 *    which is acceptable, e.g., when evaluating solar radiation at Earth since rays are
 *    essentially parallel and the target extent vanishes compared to the 1 AU distance).
 */
class IsotropicPointSourceRadiationPressureAcceleration : public RadiationPressureAcceleration
{
public:
    /*!
    * Creates a radiation pressure acceleration model with an isotropic point source.
    *
    * @param sourceModel Radiation source model of the source body
    * @param sourceBodyShapeModel Body shape model of the source body
    * @param sourcePositionFunction Position function of the source body
    * @param targetModel Radiation pressure target model of the target body
    * @param targetPositionFunction Position function of the target body
    * @param targetRotationFromLocalToGlobalFrameFunction Local-to-global rotation function of the target body
    * @param targetMassFunction Mass function of the target body
    * @param sourceToTargetOccultationModel Occultation model between source and target
    */
    IsotropicPointSourceRadiationPressureAcceleration(
            const std::shared_ptr<IsotropicPointRadiationSourceModel>& sourceModel,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& sourceBodyShapeModel,
            const std::function<Eigen::Vector3d()>& sourcePositionFunction,
            const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
            const std::function<Eigen::Vector3d()>& targetPositionFunction,
            const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
            const std::function<double()>& targetMassFunction,
            const std::shared_ptr<OccultationModel>& sourceToTargetOccultationModel) :
            RadiationPressureAcceleration(sourcePositionFunction, targetModel,
                                          targetPositionFunction, targetRotationFromLocalToGlobalFrameFunction,
                                          targetMassFunction, sourceToTargetOccultationModel),
                                          sourceModel_(sourceModel),
                                          sourceBodyShapeModel_(sourceBodyShapeModel),
            sourceToTargetReceivedFraction(TUDAT_NAN) {}

    std::shared_ptr<RadiationSourceModel> getSourceModel() const override
    {
        return sourceModel_;
    }

    double getSourceToTargetReceivedFraction() const
    {
        return sourceToTargetReceivedFraction;
    }

private:
    Eigen::Vector3d calculateAcceleration() override;

    std::shared_ptr<IsotropicPointRadiationSourceModel> sourceModel_;
    std::shared_ptr<basic_astrodynamics::BodyShapeModel> sourceBodyShapeModel_;

    // For dependent variable
    double sourceToTargetReceivedFraction;
};

/*!
 * Class modeling radiation pressure acceleration from a paneled source. An original source illuminates the
 * source, which emits the radiation that accelerates the target. Only isotropic point sources are supported as
 * original sources for now. For occultation, panel-center-to-target-center visibility is used for each panel instead
 * of a shadow function.
 *
 * Assumptions:
 *  - Radiation from the original source is evaluated at the source center (i.e. the source extent is
 *    neglected, which is acceptable, e.g., when evaluating solar radiation at Earth for albedo). Occultation using
 *    the shadow function is done for each panel at the panel center, however.
 *  - Radiation from each panel/subsource of the paneled source is evaluated at the target center (i.e. the target
 *    extent is neglected). This also applies for occultation calculations.
 */
 // If other original source types are to be supported in the future, add the rotational state update to
 // createEnvironmentUpdater.cpp under case radiation_pressure
class PaneledSourceRadiationPressureAcceleration : public RadiationPressureAcceleration
{
public:
    /*!
    * Creates a radiation pressure acceleration model with a paneled source.
    *
    * @param sourceModel Radiation source model of the source body
    * @param sourcePositionFunction Position function of the source body
    * @param sourceRotationFromLocalToGlobalFrameFunction Local-to-global rotation function of the source body
    * @param targetModel Radiation pressure target model of the target body
    * @param targetPositionFunction Position function of the target body
    * @param targetRotationFromLocalToGlobalFrameFunction Local-to-global rotation function of the target body
    * @param targetMassFunction Mass function of the target body
    * @param originalSourceModel Radiation source model of the original source body
    * @param originalSourceBodyShapeModel Body shape model of the original source body
    * @param originalSourcePositionFunction Position function of the original source body
    * @param sourceToTargetOccultationModel Occultation model between source and target
    * @param originalSourceToSourceOccultationModel Occultation model between original source and source
    */
    PaneledSourceRadiationPressureAcceleration(
            const std::shared_ptr<PaneledRadiationSourceModel>& sourceModel,
            const std::function<Eigen::Vector3d()>& sourcePositionFunction,
            const std::function<Eigen::Quaterniond()>& sourceRotationFromLocalToGlobalFrameFunction,
            const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
            const std::function<Eigen::Vector3d()>& targetPositionFunction,
            const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
            const std::function<double()>& targetMassFunction,
            const std::shared_ptr<IsotropicPointRadiationSourceModel>& originalSourceModel,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& originalSourceBodyShapeModel,
            const std::function<Eigen::Vector3d()>& originalSourcePositionFunction,
            const std::shared_ptr<OccultationModel>& sourceToTargetOccultationModel,
            const std::shared_ptr<OccultationModel>& originalSourceToSourceOccultationModel) :
            RadiationPressureAcceleration(
                    sourcePositionFunction,
                    targetModel,
                    targetPositionFunction,
                    targetRotationFromLocalToGlobalFrameFunction,
                    targetMassFunction,
                    sourceToTargetOccultationModel),
            sourceModel_(sourceModel),
            sourceRotationFromLocalToGlobalFrameFunction_(sourceRotationFromLocalToGlobalFrameFunction),
            originalSourceModel_(originalSourceModel),
            originalSourceBodyShapeModel_(originalSourceBodyShapeModel),
            originalSourcePositionFunction_(originalSourcePositionFunction),
            originalSourceToSourceOccultationModel_(originalSourceToSourceOccultationModel),
            visibleSourcePanelCount(-1),
            illuminatedSourcePanelCount(-1),
            visibleAndIlluminatedSourcePanelCount(-1) {}

    std::shared_ptr<RadiationSourceModel> getSourceModel() const override
    {
        return sourceModel_;
    }

    const std::shared_ptr<IsotropicPointRadiationSourceModel>& getOriginalSourceModel() const
    {
        return originalSourceModel_;
    }

    unsigned int getVisibleSourcePanelCount() const
    {
        return visibleSourcePanelCount;
    }

    unsigned int getIlluminatedSourcePanelCount() const
    {
        return illuminatedSourcePanelCount;
    }

    unsigned int getVisibleAndIlluminatedSourcePanelCount() const
    {
        return visibleAndIlluminatedSourcePanelCount;
    }

    double getVisibleSourceArea() const
    {
        return sourceModel_->getVisibleArea();
    }

private:
    Eigen::Vector3d calculateAcceleration() override;

    void updateMembers_(double currentTime) override;

    std::shared_ptr<PaneledRadiationSourceModel> sourceModel_;
    std::function<Eigen::Quaterniond()> sourceRotationFromLocalToGlobalFrameFunction_;
    std::shared_ptr<IsotropicPointRadiationSourceModel> originalSourceModel_;
    std::shared_ptr<basic_astrodynamics::BodyShapeModel> originalSourceBodyShapeModel_;
    std::function<Eigen::Vector3d()> originalSourcePositionFunction_;
    std::shared_ptr<OccultationModel> originalSourceToSourceOccultationModel_;

    // For dependent variables
    unsigned int visibleSourcePanelCount;
    unsigned int illuminatedSourcePanelCount;
    unsigned int visibleAndIlluminatedSourcePanelCount;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSUREACCELERATION_H
