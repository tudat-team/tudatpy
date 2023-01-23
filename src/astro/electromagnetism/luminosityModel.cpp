/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Oliver Montenbruck, et al. Satellite Orbits. Springer Berlin Heidelberg, 2000.
 */

#include "tudat/astro/electromagnetism/luminosityModel.h"

#include "tudat/math/basic/mathematicalConstants.h"


namespace tudat
{
namespace electromagnetism
{

void LuminosityModel::updateMembers(double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

double IrradianceBasedLuminosityModel::getLuminosity() const
{
    // The source is assumed to be isotropic
    auto sphereArea = 4 * mathematical_constants::PI * distance_ * distance_;
    auto luminosity = irradianceAtDistance_ * sphereArea;
    return luminosity;
}

void IrradianceBasedLuminosityModel::updateMembers_(const double currentTime)
{
    // Evaluate only once per timestep since irradiance function could be expensive to evaluate
    irradianceAtDistance_ = irradianceAtDistanceFunction_(currentTime);
}
} // tudat
} // electromagnetism
