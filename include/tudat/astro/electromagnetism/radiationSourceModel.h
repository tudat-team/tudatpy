#ifndef TUDAT_RADIATIONSOURCEMODEL_H
#define TUDAT_RADIATIONSOURCEMODEL_H

#include <vector>
#include <tuple>
#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/electromagnetism/luminosityModel.h"

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
            const std::shared_ptr<LuminosityModel> luminosityModel):
            RadiationSourceModel(sourcePositionFunction),
            luminosityModel_(luminosityModel) {}

    IrradianceWithSourceList evaluateIrradianceAtPosition(Eigen::Vector3d targetPosition) const override
    {
        auto sourcePosition = sourcePositionFunction_();
        auto distanceSourceToTarget = (targetPosition - sourcePosition).norm();
        auto luminosity = luminosityModel_->getLuminosity();

        auto sphereArea = 4 * mathematical_constants::PI * distanceSourceToTarget * distanceSourceToTarget;
        auto irradiance = luminosity / sphereArea;

        return IrradianceWithSourceList { std::make_tuple(irradiance, sourcePosition) };
    }

    std::shared_ptr<LuminosityModel> getLuminosityModel() const
    {
        return luminosityModel_;
    }

private:
    std::shared_ptr<LuminosityModel> luminosityModel_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONSOURCEMODEL_H
