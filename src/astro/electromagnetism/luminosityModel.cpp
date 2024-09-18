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

double computeLuminosityFromIrradiance( const double irradiance, const double distance )
{
    return 4.0 * mathematical_constants::PI * distance * distance * irradiance;
}

void LuminosityModel::updateMembers(double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

} // tudat
} // electromagnetism
