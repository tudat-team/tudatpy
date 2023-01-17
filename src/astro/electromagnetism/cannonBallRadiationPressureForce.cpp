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

#include "tudat/astro/electromagnetism/cannonBallRadiationPressureForce.h"

namespace tudat
{
namespace electromagnetism
{

//! Compute radiation pressure force using a cannon-ball model.
Eigen::Vector3d computeCannonBallRadiationPressureForce(
        const double radiationPressure,
        const Eigen::Vector3d& vectorToSource,
        const double area,
        const double radiationPressureCoefficient )
{
    return -radiationPressure * radiationPressureCoefficient * area * vectorToSource;
}

} // namespace electromagnetism
} // namespace tudat
