#ifndef TUDAT_RADIATIONPRESSUREACCELERATION_H
#define TUDAT_RADIATIONPRESSUREACCELERATION_H

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/simulation/environment_setup/body.h"
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
