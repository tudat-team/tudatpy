#ifndef TUDATBUNDLE_RADIATIONPRESSUREACCELERATION_H
#define TUDATBUNDLE_RADIATIONPRESSUREACCELERATION_H

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/simulation/environment_setup/body.h"
#include "radiationSourceInterface.h"
#include "radiationPressureTargetInterface.h"

namespace tudat
{
namespace electromagnetism
{

class RadiationPressureAcceleration: public basic_astrodynamics::AccelerationModel3d
{
public:
    RadiationPressureAcceleration(const std::shared_ptr<RadiationSourceInterface> sourceInterface,
                                  const std::shared_ptr<RadiationPressureTargetInterface> targetInterface,
                                  const std::function<Eigen::Vector3d()> targetPositionFunction,
                                  const std::function<double()> targetMassFunction,
                                  const std::function<Eigen::Quaterniond()> targetRotationFromLocalToPropagationFrameFunction)
            :
            sourceInterface_(sourceInterface),
            targetInterface_(targetInterface),
            targetPositionFunction_(targetPositionFunction),
            targetMassFunction_(targetMassFunction),
            targetRotationFromLocalToPropagationFrameFunction_(targetRotationFromLocalToPropagationFrameFunction)
    {}

    void updateMembers(const double currentTime) override;

private:
    Eigen::Vector3d calculateAcceleration();


    std::shared_ptr<RadiationSourceInterface> sourceInterface_;

    std::shared_ptr<RadiationPressureTargetInterface> targetInterface_;

    std::function<Eigen::Vector3d()> targetPositionFunction_;

    std::function<double()> targetMassFunction_;

    std::function< Eigen::Quaterniond( ) > targetRotationFromLocalToPropagationFrameFunction_;
};

} // tudat
} // electromagnetism

#endif //TUDATBUNDLE_RADIATIONPRESSUREACCELERATION_H
