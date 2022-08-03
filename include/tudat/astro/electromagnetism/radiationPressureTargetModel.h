#ifndef TUDATBUNDLE_RADIATIONPRESSURETARGETMODEL_H
#define TUDATBUNDLE_RADIATIONPRESSURETARGETMODEL_H

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{
namespace electromagnetism
{

// All calculations in body-fixed frame
class RadiationPressureTargetModel
{
public:
    virtual ~RadiationPressureTargetModel() = default;

    virtual Eigen::Vector3d evaluateRadiationPressureForce(
            double sourceIrradiance, Eigen::Vector3d sourceToTargetDirection) const = 0;
};


class CannonballRadiationPressureTargetModel : public RadiationPressureTargetModel
{
public:
    CannonballRadiationPressureTargetModel(double area, double coefficient)
            : area_(area), coefficient_(coefficient) {}

    Eigen::Vector3d evaluateRadiationPressureForce(
            double sourceIrradiance,
            Eigen::Vector3d sourceToTargetDirection) const override
    {
        auto radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
        auto forceMagnitude = coefficient_ * area_ * radiationPressure;
        auto force = forceMagnitude * sourceToTargetDirection;
        return force;
    }

private:
    double area_{};
    double coefficient_{};
};


} // tudat
} // electromagnetism

#endif //TUDATBUNDLE_RADIATIONPRESSURETARGETMODEL_H
