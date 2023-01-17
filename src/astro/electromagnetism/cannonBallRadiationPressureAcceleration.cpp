/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */
//RP-OLD

#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/astro/electromagnetism/cannonBallRadiationPressureAcceleration.h"
#include "tudat/astro/electromagnetism/cannonBallRadiationPressureForce.h"

namespace tudat
{
namespace electromagnetism
{

//! Compute radiation pressure acceleration using a cannon-ball model.
Eigen::Vector3d computeCannonBallRadiationPressureAcceleration(
        const double radiationPressure,
        const Eigen::Vector3d& vectorToSource,
        const double area,
        const double radiationPressureCoefficient,
        const double mass )
{
    return computeCannonBallRadiationPressureForce(
                radiationPressure, vectorToSource, area, radiationPressureCoefficient ) / mass;
}

} // namespace electromagnetism
} // namespace tudat
