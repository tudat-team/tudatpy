#ifndef TUDAT_RADIATIONPRESSUREACCELERATION_H
#define TUDAT_RADIATIONPRESSUREACCELERATION_H

#include <functional>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "radiationSourceModel.h"
#include "radiationPressureTargetModel.h"

namespace tudat
{
namespace electromagnetism
{

class RadiationPressureAcceleration: public basic_astrodynamics::AccelerationModel3d
{
public:
    RadiationPressureAcceleration(const std::shared_ptr<RadiationSourceModel> sourceModel,
                                  const std::shared_ptr<RadiationPressureTargetModel> targetModel,
                                  const std::function<Eigen::Vector3d()> targetPositionFunction,
                                  const std::function<double()> targetMassFunction,
                                  const std::function<Eigen::Quaterniond()> targetRotationFromLocalToPropagationFrameFunction)
            :
            sourceModel_(sourceModel),
            targetModel_(targetModel),
            targetPositionFunction_(targetPositionFunction),
            targetMassFunction_(targetMassFunction),
            targetRotationFromLocalToPropagationFrameFunction_(targetRotationFromLocalToPropagationFrameFunction)
    {}

    void updateMembers(const double currentTime) override;

    const std::shared_ptr<RadiationSourceModel> getSourceModel() const
    {
        return sourceModel_;
    }

    const std::shared_ptr<RadiationPressureTargetModel> getTargetModel() const
    {
        return targetModel_;
    }

private:
    Eigen::Vector3d calculateAcceleration();


    std::shared_ptr<RadiationSourceModel> sourceModel_;

    std::shared_ptr<RadiationPressureTargetModel> targetModel_;

    std::function<Eigen::Vector3d()> targetPositionFunction_;

    std::function<double()> targetMassFunction_;

    std::function< Eigen::Quaterniond( ) > targetRotationFromLocalToPropagationFrameFunction_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSUREACCELERATION_H
