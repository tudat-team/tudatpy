#include "tudat/astro/electromagnetism/luminosityModel.h"

namespace tudat
{
namespace electromagnetism
{

double IrradianceBasedLuminosityModel::getLuminosity() const
{
    auto sphereArea = 4 * mathematical_constants::PI * distance_ * distance_;
    auto luminosity = irradianceAtDistance_ * sphereArea;
    return luminosity;
}

void IrradianceBasedLuminosityModel::updateMembers()
{
    // Evaluate only once per timestep since irradiance function could be expensive to evaluate
    irradianceAtDistance_ = irradianceAtDistanceFunction_();
}
} // tudat
} // electromagnetism
