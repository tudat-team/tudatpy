#include "tudat/astro/electromagnetism/radiationSourceModel.h"

namespace tudat
{
namespace electromagnetism
{


IrradianceWithSourceList IsotropicPointRadiationSourceModel::evaluateIrradianceAtPosition(Eigen::Vector3d targetPosition) const
{
    auto sourcePosition = sourcePositionFunction_();
    auto distanceSourceToTarget = (targetPosition - sourcePosition).norm();
    auto luminosity = luminosityModel_->getLuminosity();

    auto sphereArea = 4 * mathematical_constants::PI * distanceSourceToTarget * distanceSourceToTarget;
    auto irradiance = luminosity / sphereArea;

    return IrradianceWithSourceList { std::make_tuple(irradiance, sourcePosition) };
}

} // tudat
} // electromagnetism
