#ifndef TUDATBUNDLE_RADIATIONSOURCEINTERFACE_H
#define TUDATBUNDLE_RADIATIONSOURCEINTERFACE_H

#include <iostream>
#include <vector>
#include <tuple>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/electromagnetism/radiantPowerModel.h"


namespace tudat
{
namespace electromagnetism
{

typedef std::vector<std::tuple<double, Eigen::Vector3d>> IrradianceWithSourceList;


class RadiationSourceInterface
{
public:
    explicit RadiationSourceInterface(const std::function< Eigen::Vector3d( ) > sourcePositionFunction):
            sourcePositionFunction_(sourcePositionFunction) {}

    virtual ~RadiationSourceInterface() = default;

    virtual IrradianceWithSourceList evaluateIrradianceAtPosition(Eigen::Vector3d) = 0;

protected:
    std::function< Eigen::Vector3d( ) > sourcePositionFunction_;
};


class IsotropicPointRadiationSourceInterface : public RadiationSourceInterface
{
public:
    explicit IsotropicPointRadiationSourceInterface(
            const std::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const std::shared_ptr<RadiantPowerModel> radiantPowerModel):
        RadiationSourceInterface(sourcePositionFunction),
        radiantPowerModel_(radiantPowerModel) {}

    IrradianceWithSourceList evaluateIrradianceAtPosition(Eigen::Vector3d targetPosition) override
    {
        auto sourcePosition = sourcePositionFunction_();
        auto distanceSourceToTarget = (targetPosition - sourcePosition).norm();
        auto radiantPower = radiantPowerModel_->getRadiantPower();

        auto sphereArea = 4 * mathematical_constants::PI * distanceSourceToTarget * distanceSourceToTarget;
        auto irradiance = radiantPower / sphereArea;

        return IrradianceWithSourceList { std::make_tuple(irradiance, sourcePosition) };
    }

private:
    std::shared_ptr<RadiantPowerModel> radiantPowerModel_;
};

} // tudat
} // electromagnetism

#endif //TUDATBUNDLE_RADIATIONSOURCEINTERFACE_H
