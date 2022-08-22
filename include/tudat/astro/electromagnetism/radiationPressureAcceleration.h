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

typedef std::function<double(
        const Eigen::Vector3d& originalSourcePosition,
        const Eigen::Vector3d& sourcePosition)> OriginalSourceOccultationFactorFunctionType;

typedef std::function<double(
        const Eigen::Vector3d& sourcePosition,
        const Eigen::Vector3d& targetPosition)> SourceOccultationFactorFunctionType;

class RadiationPressureAcceleration: public basic_astrodynamics::AccelerationModel3d
{
public:
    RadiationPressureAcceleration(const std::shared_ptr<RadiationSourceModel>& sourceModel,
                                  const std::function<Eigen::Vector3d()>& sourcePositionFunction,
                                  const std::function<Eigen::Quaterniond()>& sourceRotationFromLocalToGlobalFrameFunction,
                                  const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
                                  const std::function<Eigen::Vector3d()>& targetPositionFunction,
                                  const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
                                  const std::function<double()>& targetMassFunction,
                                  const SourceOccultationFactorFunctionType& sourceOccultationFactorFunction) :
            RadiationPressureAcceleration(sourceModel, sourcePositionFunction, sourceRotationFromLocalToGlobalFrameFunction,
                                          targetModel, targetPositionFunction, targetRotationFromLocalToGlobalFrameFunction,
                                          targetMassFunction,
                                          "", nullptr, nullptr, nullptr,
                                          sourceOccultationFactorFunction, nullptr) {}

    /*!
    * Creates a radiation pressure acceleration model. If the source model is paneled, original source members
     * must be populated such that albedo radiation can be reflected.
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
    * @param originalSourcePositionFunction the position function of the original source body
    * @param originalSourceRotationFromLocalToGlobalFrameFunction the local-to-global rotation function of the original source body
    * @param sourceOccultationFunction a function applying occultation effects to source irradiance
    * @param originalSourceVisibilityFunction a function checking the visibility between original source and source
    */
    RadiationPressureAcceleration(const std::shared_ptr<RadiationSourceModel>& sourceModel,
                                  const std::function<Eigen::Vector3d()>& sourcePositionFunction,
                                  const std::function<Eigen::Quaterniond()>& sourceRotationFromLocalToGlobalFrameFunction,
                                  const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
                                  const std::function<Eigen::Vector3d()>& targetPositionFunction,
                                  const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
                                  const std::function<double()>& targetMassFunction,
                                  const std::string& originalSourceName,
                                  const std::shared_ptr<IsotropicPointRadiationSourceModel>& originalSourceModel,
                                  const std::function<Eigen::Vector3d()>& originalSourcePositionFunction,
                                  const std::function<Eigen::Quaterniond()>& originalSourceRotationFromLocalToGlobalFrameFunction,
                                  const SourceOccultationFactorFunctionType& sourceOccultationFactorFunction,
                                  const OriginalSourceOccultationFactorFunctionType& originalSourceOccultationFactorFunction) :
            sourceModel_(sourceModel),
            sourcePositionFunction_(sourcePositionFunction),
            sourceRotationFromLocalToGlobalFrameFunction_(sourceRotationFromLocalToGlobalFrameFunction),
            targetModel_(targetModel),
            targetPositionFunction_(targetPositionFunction),
            targetRotationFromLocalToGlobalFrameFunction_(targetRotationFromLocalToGlobalFrameFunction),
            targetMassFunction_(targetMassFunction),
            originalSourceName_(originalSourceName),
            originalSourceModel_(originalSourceModel),
            originalSourcePositionFunction_(originalSourcePositionFunction),
            originalSourceRotationFromLocalToGlobalFrameFunction_(originalSourceRotationFromLocalToGlobalFrameFunction),
            sourceOccultationFactorFunction_(sourceOccultationFactorFunction),
            originalSourceOccultationFactorFunction_(originalSourceOccultationFactorFunction) {}

    void updateMembers(double currentTime) override;

    std::shared_ptr<RadiationSourceModel> getSourceModel() const
    {
        return sourceModel_;
    }

    std::shared_ptr<RadiationPressureTargetModel> getTargetModel() const
    {
        return targetModel_;
    }

    const std::shared_ptr<IsotropicPointRadiationSourceModel>& getOriginalSourceModel() const
    {
        return originalSourceModel_;
    }

    const std::string& getOriginalSourceName() const
    {
        return originalSourceName_;
    }

private:
    Eigen::Vector3d calculateAcceleration();

    std::pair<double, Eigen::Vector3d> calculateOriginalSourceIrradiance(
            const Eigen::Vector3d& sourceCenterPositionInGlobalFrame,
            const Eigen::Quaterniond& sourceRotationFromGlobalToLocalFrame) const;

    double calculateTotalOccultationFactor(
            const Eigen::Vector3d& originalSourcePosition,
            const Eigen::Vector3d& sourcePosition,
            const Eigen::Vector3d& targetPosition) const;


    std::shared_ptr<RadiationSourceModel> sourceModel_;

    std::function<Eigen::Vector3d()> sourcePositionFunction_;

    std::function< Eigen::Quaterniond( ) > sourceRotationFromLocalToGlobalFrameFunction_;

    std::shared_ptr<RadiationPressureTargetModel> targetModel_;

    std::function<Eigen::Vector3d()> targetPositionFunction_;

    std::function< Eigen::Quaterniond( ) > targetRotationFromLocalToGlobalFrameFunction_;

    std::function<double()> targetMassFunction_;

    // Original source members are only populated if source is paneled
    std::string originalSourceName_; // needed for environment updater

    std::shared_ptr<IsotropicPointRadiationSourceModel> originalSourceModel_;

    std::function<Eigen::Vector3d()> originalSourcePositionFunction_;

    std::function< Eigen::Quaterniond( ) > originalSourceRotationFromLocalToGlobalFrameFunction_;

    SourceOccultationFactorFunctionType sourceOccultationFactorFunction_;

    OriginalSourceOccultationFactorFunctionType originalSourceOccultationFactorFunction_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSUREACCELERATION_H
