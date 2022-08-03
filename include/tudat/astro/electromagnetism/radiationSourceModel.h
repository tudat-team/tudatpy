#ifndef TUDAT_RADIATIONSOURCEMODEL_H
#define TUDAT_RADIATIONSOURCEMODEL_H

#include <vector>
#include <tuple>
#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/electromagnetism/radiantPowerModel.h"

namespace tudat
{
namespace electromagnetism
{

typedef std::vector<std::tuple<double, Eigen::Vector3d>> IrradianceWithSourceList;


class RadiationSourceModel
{
public:
    explicit RadiationSourceModel(const std::function< Eigen::Vector3d( ) > sourcePositionFunction):
            sourcePositionFunction_(sourcePositionFunction) {}

    virtual ~RadiationSourceModel() = default;

    virtual IrradianceWithSourceList evaluateIrradianceAtPosition(Eigen::Vector3d) const = 0;

protected:
    std::function< Eigen::Vector3d( ) > sourcePositionFunction_;
};


class IsotropicPointRadiationSourceModel : public RadiationSourceModel
{
public:
    explicit IsotropicPointRadiationSourceModel(
            const std::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const std::shared_ptr<RadiantPowerModel> radiantPowerModel):
            RadiationSourceModel(sourcePositionFunction),
            radiantPowerModel_(radiantPowerModel) {}

    IrradianceWithSourceList evaluateIrradianceAtPosition(Eigen::Vector3d targetPosition) const override
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

#endif //TUDAT_RADIATIONSOURCEMODEL_H
