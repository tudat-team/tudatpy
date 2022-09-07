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
#include "radiationSourceModel.h"
#include "radiationPressureTargetModel.h"
#include "occultationModel.h"

namespace tudat
{
namespace electromagnetism
{

class RadiationPressureAcceleration: public basic_astrodynamics::AccelerationModel3d
{
public:
    void updateMembers(double currentTime) override;

    std::shared_ptr<RadiationSourceModel> getSourceModel() const
    {
        return sourceModel_;
    }

    std::shared_ptr<RadiationPressureTargetModel> getTargetModel() const
    {
        return targetModel_;
    }

protected:
    RadiationPressureAcceleration(const std::shared_ptr<RadiationSourceModel>& sourceModel,
                                  const std::function<Eigen::Vector3d()>& sourcePositionFunction,
                                  const std::function<Eigen::Quaterniond()>& sourceRotationFromLocalToGlobalFrameFunction,
                                  const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
                                  const std::function<Eigen::Vector3d()>& targetPositionFunction,
                                  const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
                                  const std::function<double()>& targetMassFunction,
                                  const std::shared_ptr<OccultationModel>& occultationModel) :
            sourceModel_(sourceModel),
            sourcePositionFunction_(sourcePositionFunction),
            sourceRotationFromLocalToGlobalFrameFunction_(sourceRotationFromLocalToGlobalFrameFunction),
            targetModel_(targetModel),
            targetPositionFunction_(targetPositionFunction),
            targetRotationFromLocalToGlobalFrameFunction_(targetRotationFromLocalToGlobalFrameFunction),
            targetMassFunction_(targetMassFunction),
            occultationModel_(occultationModel) {}

    virtual Eigen::Vector3d calculateAcceleration() = 0;

    std::shared_ptr<RadiationSourceModel> sourceModel_;

    std::function<Eigen::Vector3d()> sourcePositionFunction_;

    std::function< Eigen::Quaterniond( ) > sourceRotationFromLocalToGlobalFrameFunction_;

    std::shared_ptr<RadiationPressureTargetModel> targetModel_;

    std::function<Eigen::Vector3d()> targetPositionFunction_;

    std::function< Eigen::Quaterniond( ) > targetRotationFromLocalToGlobalFrameFunction_;

    std::function<double()> targetMassFunction_;

    std::shared_ptr<OccultationModel> occultationModel_;
};

class IsotropicPointSourceRadiationPressureAcceleration : public RadiationPressureAcceleration
{
public:
    /*!
    * Creates a radiation pressure acceleration model with an isotropic point source.
    *
    * @param sourceModel the radiation source model of the source body
    * @param sourceBodyShapeModel the body shape model of the source body
    * @param sourcePositionFunction the position function of the source body
    * @param sourceRotationFromLocalToGlobalFrameFunction the local-to-global rotation function of the source body
    * @param targetModel the radiation pressure target model of the target body
    * @param targetPositionFunction the position function of the target body
    * @param targetRotationFromLocalToGlobalFrameFunction the local-to-global rotation function of the target body
    * @param targetMassFunction the mass function of the target body
    * @param occultationModel the occultation model between original source, source and target
    */
    IsotropicPointSourceRadiationPressureAcceleration(
            const std::shared_ptr<IsotropicPointRadiationSourceModel>& sourceModel,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& sourceBodyShapeModel,
            const std::function<Eigen::Vector3d()>& sourcePositionFunction,
            const std::function<Eigen::Quaterniond()>& sourceRotationFromLocalToGlobalFrameFunction,
            const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
            const std::function<Eigen::Vector3d()>& targetPositionFunction,
            const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
            const std::function<double()>& targetMassFunction,
            const std::shared_ptr<OccultationModel>& occultationModel)
            : RadiationPressureAcceleration(sourceModel, sourcePositionFunction,
                                            sourceRotationFromLocalToGlobalFrameFunction, targetModel,
                                            targetPositionFunction, targetRotationFromLocalToGlobalFrameFunction,
                                            targetMassFunction, occultationModel),
                                            isotropicPointSourceModel_(sourceModel),
                                            sourceBodyShapeModel_(sourceBodyShapeModel) {}

private:
    Eigen::Vector3d calculateAcceleration() override;

    std::shared_ptr<IsotropicPointRadiationSourceModel> isotropicPointSourceModel_;

    std::shared_ptr<basic_astrodynamics::BodyShapeModel> sourceBodyShapeModel_;
};

class PaneledSourceRadiationPressureAcceleration : public RadiationPressureAcceleration
{
public:
    /*!
    * Creates a radiation pressure acceleration model with a paneled source.
    *
    * @param sourceModel the radiation source model of the source body
    * @param sourcePositionFunction the position function of the source body
    * @param sourceRotationFromLocalToGlobalFrameFunction the local-to-global rotation function of the source body
    * @param targetModel the radiation pressure target model of the target body
    * @param targetPositionFunction the position function of the target body
    * @param targetRotationFromLocalToGlobalFrameFunction the local-to-global rotation function of the target body
    * @param targetMassFunction the mass function of the target body
    * @param originalSourceName the name of the original source body
    * @param originalSourceModel the radiation source model of the original source body
    * @param originalSourceBodyShapeModel the body shape model of the original source body
    * @param originalSourcePositionFunction the position function of the original source body
    * @param originalSourceRotationFromLocalToGlobalFrameFunction the local-to-global rotation function of the original source body
    * @param occultationModel the occultation model between original source, source and target
    */
    PaneledSourceRadiationPressureAcceleration(
            const std::shared_ptr<PaneledRadiationSourceModel>& sourceModel,
            const std::function<Eigen::Vector3d()>& sourcePositionFunction,
            const std::function<Eigen::Quaterniond()>& sourceRotationFromLocalToGlobalFrameFunction,
            const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
            const std::function<Eigen::Vector3d()>& targetPositionFunction,
            const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
            const std::function<double()>& targetMassFunction,
            const std::string& originalSourceName,
            const std::shared_ptr<IsotropicPointRadiationSourceModel>& originalSourceModel,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& originalSourceBodyShapeModel,
            const std::function<Eigen::Vector3d()>& originalSourcePositionFunction,
            const std::function<Eigen::Quaterniond()>& originalSourceRotationFromLocalToGlobalFrameFunction,
            const std::shared_ptr<OccultationModel>& occultationModel) :
            RadiationPressureAcceleration(
                    sourceModel,
                    sourcePositionFunction,
                    sourceRotationFromLocalToGlobalFrameFunction,
                    targetModel,
                    targetPositionFunction,
                    targetRotationFromLocalToGlobalFrameFunction,
                    targetMassFunction,
                    occultationModel),
            originalSourceName_(originalSourceName),
            originalSourceModel_(originalSourceModel),
            originalSourceBodyShapeModel_(originalSourceBodyShapeModel),
            originalSourcePositionFunction_(originalSourcePositionFunction),
            originalSourceRotationFromLocalToGlobalFrameFunction_(originalSourceRotationFromLocalToGlobalFrameFunction) {}

    const std::shared_ptr<IsotropicPointRadiationSourceModel>& getOriginalSourceModel() const
    {
        return originalSourceModel_;
    }

    const std::string& getOriginalSourceName() const
    {
        return originalSourceName_;
    }
private:
    Eigen::Vector3d calculateAcceleration() override;

    std::string originalSourceName_; // needed for environment updater

    std::shared_ptr<IsotropicPointRadiationSourceModel> originalSourceModel_;

    std::shared_ptr<basic_astrodynamics::BodyShapeModel> originalSourceBodyShapeModel_;

    std::function<Eigen::Vector3d()> originalSourcePositionFunction_;

    std::function< Eigen::Quaterniond( ) > originalSourceRotationFromLocalToGlobalFrameFunction_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSUREACCELERATION_H
