#ifndef TUDATBUNDLE_RADIANTPOWERMODEL_H
#define TUDATBUNDLE_RADIANTPOWERMODEL_H

#include "tudat/math/basic/mathematicalConstants.h"


namespace tudat
{
namespace electromagnetism
{


class RadiantPowerModel
{
public:
    RadiantPowerModel() = default;

    virtual ~RadiantPowerModel() = default;

    virtual double getRadiantPower() const = 0;
};


class ConstantRadiantPowerModel : public RadiantPowerModel
{
public:
    explicit ConstantRadiantPowerModel(double radiantPower) : radiantPower(radiantPower) {}

    double getRadiantPower() const override
    {
        return radiantPower;
    }

private:
    double radiantPower;
};


class IrradianceBasedRadiantPowerModel : public RadiantPowerModel
{
public:
    IrradianceBasedRadiantPowerModel(const std::function<double()> irradianceAtDistanceFunction, double distance)
            : irradianceAtDistanceFunction(irradianceAtDistanceFunction), distance(distance) {}

    double getRadiantPower() const override
    {
        auto sphereArea = 4 * mathematical_constants::PI * distance * distance;
        auto radiantPower = irradianceAtDistanceFunction() * sphereArea;
        return radiantPower;
    }

private:
    std::function<double()> irradianceAtDistanceFunction;
    double distance;
};


}
}


#endif //TUDATBUNDLE_RADIANTPOWERMODEL_H
