#include "tudat/astro/electromagnetism/radiationSourceModel.h"

namespace tudat
{
namespace electromagnetism
{

void RadiationSourceModel::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

//*********************************************************************************************
//   Isotropic point radiation source
//*********************************************************************************************

IrradianceWithSourceList IsotropicPointRadiationSourceModel::evaluateIrradianceAtPosition(Eigen::Vector3d targetPosition) const
{
    auto distanceSourceToTarget = targetPosition.norm();
    auto luminosity = luminosityModel_->getLuminosity();

    auto sphereArea = 4 * mathematical_constants::PI * distanceSourceToTarget * distanceSourceToTarget;
    auto irradiance = luminosity / sphereArea;

    return IrradianceWithSourceList { std::make_pair(irradiance, Eigen::Vector3d::Zero()) };
}

void RadiationSourceModel::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        sourcePosition_ = sourcePositionFunction_();
        updateMembers_(currentTime);
    }
}

void IsotropicPointRadiationSourceModel::updateMembers_(double currentTime)
{
    luminosityModel_->updateMembers();
}
} // tudat
} // electromagnetism
